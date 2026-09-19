// landscape.cpp -- fast simulator for "Economic Landscapes".
//
// Profit is the quadratic form over the integer budget simplex
// (x_i >= 0, sum_i x_i = B):
//
//   Pi(x) = sum_i sq_i x_i^2 + 0.5 * sum_{i,j} cross_ij x_i x_j
//
// identical to the R reference engine in R/simulations.R with the linear
// terms dropped (Pi = x' Q x, Q = diag(sq) + 0.5*cross, cross symmetric,
// zero diagonal).  Search uses the cyclic one-dollar-transfer operator:
// with C connections, input i may give a dollar to (i+1),...,(i+C) mod N.
// Hill-climbing strategies: SA (steepest), LA (least), MA (median improving).
//
// The hot path uses an incremental Delta-profit: moving one dollar from a
// to b changes profit by
//   sq_a(1-2 x_a) + sq_b(1+2 x_b) + (h_b - h_a) - cross_ab,   h = cross * x
// which is O(1) per candidate move (O(N) to refresh h after a committed
// move), instead of recomputing the whole quadratic.
//
// Usage:
//   landscape check <case.txt>
//   landscape sweep <N> <B> <regime> <C1,C2,...> <nland> <nstart> <seed>
//     regime in {uniform, allpos, signed, leontief}; prints CSV to stdout.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <set>
#include <numeric>
#include <functional>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
using namespace std;

static mt19937_64 RNG;   // declared early: build_moves_topo (random topology) uses it

struct Landscape {
  int N;
  vector<double> sq;     // length N
  vector<double> cross;  // N*N, symmetric, zero diagonal
  inline double cij(int i, int j) const { return cross[(size_t)i * N + j]; }
};

static double profit_full(const vector<int>& x, const Landscape& L) {
  int N = L.N; double s = 0.0, c = 0.0;
  for (int i = 0; i < N; i++) s += L.sq[i] * (double)x[i] * x[i];
  for (int i = 0; i < N; i++) {
    double xi = x[i]; if (xi == 0) continue;
    const double* row = &L.cross[(size_t)i * N];
    for (int j = 0; j < N; j++) c += row[j] * xi * x[j];
  }
  return s + 0.5 * c;
}

struct Moves { vector<int> from, to; };
static Moves build_moves(int N, int C) {
  int k = min(C, N - 1); Moves m;
  for (int i = 0; i < N; i++)
    for (int off = 1; off <= k; off++) { m.from.push_back(i); m.to.push_back((i + off) % N); }
  return m;
}

// Neighbourhood topologies, for the operator-robustness test:
//   0 = directed ring  (default): i -> i+1,...,i+C  (mod N), one-way
//   1 = symmetric ring         : i <-> i±1,...,i±C, every connection reversible
//   2 = random partners        : C distinct random targets per input (directed)
static Moves build_moves_topo(int N, int C, int topo) {
  int k = min(C, N - 1); Moves m;
  if (topo == 1) {
    for (int i = 0; i < N; i++)
      for (int off = 1; off <= k; off++) {
        m.from.push_back(i); m.to.push_back((i + off) % N);
        m.from.push_back(i); m.to.push_back((i - off + N) % N);
      }
  } else if (topo == 2) {
    vector<int> pool(N - 1);
    for (int i = 0; i < N; i++) {
      int t = 0; for (int j = 0; j < N; j++) if (j != i) pool[t++] = j;   // t = N-1
      for (int s = 0; s < k; s++) {                                       // partial shuffle
        int r = s + (int)(uniform_real_distribution<double>(0, 1)(RNG) * (t - s));
        swap(pool[s], pool[r]); m.from.push_back(i); m.to.push_back(pool[s]);
      }
    }
  } else {
    for (int i = 0; i < N; i++)
      for (int off = 1; off <= k; off++) { m.from.push_back(i); m.to.push_back((i + off) % N); }
  }
  return m;
}

// strat: 0=SA (max delta), 1=LA (min positive delta), 2=MA (median improving)
// Returns the local-optimum profit; if walk_len != nullptr, writes the number
// of accepted moves (walk length) there; if xout != nullptr, writes the final
// allocation there (for optima-location statistics).
static double hillclimb(vector<int> x, const Landscape& L, const Moves& mv, int strat,
                        long* walk_len = nullptr, vector<int>* xout = nullptr) {
  int N = L.N;
  vector<double> h(N, 0.0);
  int M = (int)mv.from.size();
  vector<double> dz; vector<int> idx; dz.reserve(M); idx.reserve(M);
  long steps = 0;
  while (true) {
    // Recompute h = cross * x exactly each step. This costs O(N * active) but
    // avoids the drift that, under a *reversible* operator, can make both a->b
    // and b->a look improving and trap the walk in an infinite ping-pong.
    fill(h.begin(), h.end(), 0.0);
    for (int i = 0; i < N; i++) {
      double xi = x[i]; if (xi == 0) continue;
      const double* row = &L.cross[(size_t)i * N];
      for (int j = 0; j < N; j++) h[j] += row[j] * xi;
    }
    dz.clear(); idx.clear();
    for (int m = 0; m < M; m++) {
      int a = mv.from[m]; if (x[a] < 1) continue; int b = mv.to[m];
      double d = L.sq[a] * (1 - 2.0 * x[a]) + L.sq[b] * (1 + 2.0 * x[b])
               + (h[b] - h[a]) - L.cij(a, b);
      if (d > 0) { dz.push_back(d); idx.push_back(m); }
    }
    if (dz.empty()) break;
    int pick;
    if (strat == 0) { int am = 0; for (int t = 1; t < (int)dz.size(); t++) if (dz[t] > dz[am]) am = t; pick = am; }
    else if (strat == 1) { int am = 0; for (int t = 1; t < (int)dz.size(); t++) if (dz[t] < dz[am]) am = t; pick = am; }
    else { // median improving: ascending order, R index ceil(n/2) (1-based)
      vector<int> ord(dz.size()); iota(ord.begin(), ord.end(), 0);
      sort(ord.begin(), ord.end(), [&](int p, int q) { return dz[p] < dz[q]; });
      int mid = (int)ceil(dz.size() / 2.0) - 1; pick = ord[mid];
    }
    int m = idx[pick], a = mv.from[m], b = mv.to[m];
    x[a]--; x[b]++; steps++;
    if (steps > 2000000) break;     // safety net; never reached with exact h
  }
  if (walk_len) *walk_len = steps;
  if (xout) *xout = x;
  return profit_full(x, L);
}

static inline double U(double a, double b) {
  return a + (b - a) * uniform_real_distribution<double>(0.0, 1.0)(RNG);
}

static Landscape draw(const string& reg, int N) {
  Landscape L; L.N = N; L.sq.assign(N, 0.0); L.cross.assign((size_t)N * N, 0.0);
  auto sc = [&](int i, int j, double v) { L.cross[(size_t)i * N + j] = v; L.cross[(size_t)j * N + i] = v; };
  if (reg == "uniform") {
    for (int i = 0; i < N; i++) L.sq[i] = U(-1, 1);
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) sc(i, j, U(-1, 1));
  } else if (reg == "allpos") {
    for (int i = 0; i < N; i++) L.sq[i] = U(0, 1);
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) sc(i, j, U(0, 1));
  } else if (reg == "signed") {
    for (int i = 0; i < N; i++) L.sq[i] = U(0, 1);
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) sc(i, j, U(-1, 0));
  } else if (reg == "leontief") {
    for (int i = 0; i < N; i++) { double s = 0; for (int t = 0; t < N - 1; t++) s += U(-1, 0); L.sq[i] = s; }
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) sc(i, j, U(0, 1) + U(0, 1));
  } else { cerr << "unknown regime: " << reg << "\n"; exit(1); }
  return L;
}

// Tunable regime: the diagonal of Q (returns to scale) is centred at r,
// the off-diagonal (complementarities) ~ U[-1,1]. r = 0 reproduces the
// uniform regime. Sweeping r moves Q from negative definite (concave,
// diminishing returns) through indefinite to positive definite (convex,
// increasing returns).
static Landscape draw_tuned(int N, double r) {
  Landscape L; L.N = N; L.sq.assign(N, 0.0); L.cross.assign((size_t)N * N, 0.0);
  for (int i = 0; i < N; i++) L.sq[i] = r + U(-1, 1);
  for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) {
    double v = U(-1, 1); L.cross[(size_t)i * N + j] = v; L.cross[(size_t)j * N + i] = v;
  }
  return L;
}

static vector<int> random_plan(int N, int B) {
  vector<int> x(N, 0); uniform_int_distribution<int> box(0, N - 1);
  for (int b = 0; b < B; b++) x[box(RNG)]++; return x;
}

int main(int argc, char** argv) {
  if (argc < 2) { cerr << "usage: landscape check <file> | sweep <N> <B> <regime> <Clist> <nland> <nstart> <seed>\n"; return 1; }
  string cmd = argv[1];

  if (cmd == "check") {
    ifstream f(argv[2]); int N, B, C; f >> N >> B >> C;
    Landscape L; L.N = N; L.sq.resize(N); L.cross.assign((size_t)N * N, 0.0);
    for (int i = 0; i < N; i++) f >> L.sq[i];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { double v; f >> v; L.cross[(size_t)i * N + j] = v; }
    vector<int> st(N); for (int i = 0; i < N; i++) f >> st[i];
    Moves mv = build_moves(N, C);
    printf("SA %.10f LA %.10f MA %.10f\n",
           hillclimb(st, L, mv, 0), hillclimb(st, L, mv, 1), hillclimb(st, L, mv, 2));
    return 0;
  }

  if (cmd == "sweep") {
    if (argc < 9) { cerr << "sweep needs 7 args\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]); string reg = argv[4];
    vector<int> Cs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[6]), nstart = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    double denom = (double)B * B;
    printf("regime,N,C,C_over_N,SA,SA_se,LA,LA_se,gap,gap_se\n");
    for (int C : Cs) {
      Moves mv = build_moves(N, C);
      vector<double> SAv(nland), LAv(nland), GAPv(nland);
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        double msa = 0, mla = 0;
        for (int s = 0; s < nstart; s++) {
          vector<int> st = random_plan(N, B);
          msa += 100.0 * hillclimb(st, L, mv, 0) / denom;
          mla += 100.0 * hillclimb(st, L, mv, 1) / denom;
        }
        msa /= nstart; mla /= nstart;
        SAv[l] = msa; LAv[l] = mla; GAPv[l] = mla - msa;
      }
      auto mean = [](vector<double>& v) { double m = 0; for (double z : v) m += z; return m / v.size(); };
      auto se = [](vector<double>& v, double m) { double s = 0; for (double z : v) s += (z - m) * (z - m); return sqrt(s / (v.size() - 1) / v.size()); };
      double mSA = mean(SAv), mLA = mean(LAv), mG = mean(GAPv);
      printf("%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
             reg.c_str(), N, C, (double)C / N, mSA, se(SAv, mSA), mLA, se(LAv, mLA), mG, se(GAPv, mG));
    }
    return 0;
  }

  // grid: per (N,C) cell, report mean profit per strategy, the within-strategy
  // SD of outcomes (variability of results), and the mean cross-strategy SD
  // (how much the choice of strategy matters). Used to map the (N,C) plane.
  if (cmd == "grid") {
    if (argc < 9) { cerr << "grid needs 7 args\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]); string reg = argv[4];
    vector<int> Cs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[6]), nstart = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    double denom = (double)B * B;
    printf("regime,N,C,C_over_N,mean_SA,sd_SA_pooled,sd_SA_within,mean_LA,mean_MA,cross_sd\n");
    for (int C : Cs) {
      Moves mv = build_moves(N, C);
      vector<double> SA, LA, MA;                 // pooled over landscapes x starts
      SA.reserve((size_t)nland * nstart); LA.reserve((size_t)nland * nstart); MA.reserve((size_t)nland * nstart);
      double cross_sum = 0; long ncross = 0;     // cross-strategy SD per observation
      double within_sum = 0; long nwithin = 0;   // SA SD across starts, per landscape
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        vector<double> saL(nstart);
        for (int s = 0; s < nstart; s++) {
          vector<int> st = random_plan(N, B);
          double a = 100.0 * hillclimb(st, L, mv, 0) / denom;
          double b = 100.0 * hillclimb(st, L, mv, 1) / denom;
          double m = 100.0 * hillclimb(st, L, mv, 2) / denom;
          SA.push_back(a); LA.push_back(b); MA.push_back(m); saL[s] = a;
          double mu = (a + b + m) / 3.0;
          cross_sum += sqrt(((a-mu)*(a-mu) + (b-mu)*(b-mu) + (m-mu)*(m-mu)) / 2.0);
          ncross++;
        }
        if (nstart > 1) {                        // within-landscape SD of SA across starts
          double m = 0; for (double z : saL) m += z; m /= nstart;
          double v = 0; for (double z : saL) v += (z - m) * (z - m);
          within_sum += sqrt(v / (nstart - 1)); nwithin++;
        }
      }
      auto mean = [](vector<double>& v) { double s = 0; for (double z : v) s += z; return s / v.size(); };
      auto sd = [](vector<double>& v, double m) { double s = 0; for (double z : v) s += (z - m) * (z - m); return sqrt(s / (v.size() - 1)); };
      double mSA = mean(SA), mLA = mean(LA), mMA = mean(MA);
      printf("%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
             reg.c_str(), N, C, (double)C / N, mSA, sd(SA, mSA),
             nwithin ? within_sum / nwithin : 0.0, mLA, mMA, cross_sum / ncross);
    }
    return 0;
  }
  // optest: like grid, but for a chosen neighbourhood topology (0 directed
  // ring, 1 symmetric ring, 2 random partners). Used to check whether the
  // C/N transition is a property of the directed ring or survives reversible
  // and non-lattice operators. Reports mean profit per strategy, the mean
  // out-degree per input, and cross-strategy SD.
  if (cmd == "optest") {
    if (argc < 10) { cerr << "optest needs: <topo> <N> <B> <regime> <Clist> <nland> <nstart> <seed>\n"; return 1; }
    int topo = atoi(argv[2]), N = atoi(argv[3]), B = atoi(argv[4]); string reg = argv[5];
    vector<int> Cs; { stringstream ss(argv[6]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[7]), nstart = atoi(argv[8]);
    RNG.seed(strtoull(argv[9], nullptr, 10));
    double denom = (double)B * B;
    printf("topo,N,C,C_over_N,outdeg,mean_SA,mean_LA,mean_MA,cross_sd\n");
    for (int C : Cs) {
      Moves mv; if (topo != 2) mv = build_moves_topo(N, C, topo);
      vector<double> SA, LA, MA; double cross_sum = 0; long ncross = 0;
      double outdeg_sum = 0; long nmv = 0;
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        if (topo == 2) mv = build_moves_topo(N, C, topo);   // fresh random graph per landscape
        outdeg_sum += (double)mv.from.size(); nmv++;
        for (int s = 0; s < nstart; s++) {
          vector<int> st = random_plan(N, B);
          double a = 100.0 * hillclimb(st, L, mv, 0) / denom;
          double b = 100.0 * hillclimb(st, L, mv, 1) / denom;
          double m = 100.0 * hillclimb(st, L, mv, 2) / denom;
          SA.push_back(a); LA.push_back(b); MA.push_back(m);
          double mu = (a + b + m) / 3.0;
          cross_sum += sqrt(((a-mu)*(a-mu) + (b-mu)*(b-mu) + (m-mu)*(m-mu)) / 2.0); ncross++;
        }
      }
      auto mean = [](vector<double>& v){ double s=0; for(double z:v) s+=z; return s/v.size(); };
      printf("%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
             topo, N, C, (double)C / N, outdeg_sum / nmv / N,
             mean(SA), mean(LA), mean(MA), cross_sum / ncross);
    }
    return 0;
  }
  // diag: per (regime, C), report for each strategy (SA/LA/MA) the mean profit,
  // mean walk length (accepted moves), mean number of active inputs at the
  // optimum, and mean concentration (max budget share on one input). Feeds the
  // walk-length and optima-location results.
  if (cmd == "diag") {
    if (argc < 9) { cerr << "diag needs: <N> <B> <regime> <Clist> <nland> <nstart> <seed>\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]); string reg = argv[4];
    vector<int> Cs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[6]), nstart = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    double denom = (double)B * B;
    const char* sname[3] = {"SA", "LA", "MA"};
    printf("regime,N,C,C_over_N,strat,mean_profit,mean_walk,mean_active,mean_maxshare\n");
    for (int C : Cs) {
      Moves mv = build_moves(N, C);
      double sp[3] = {0,0,0}, sw[3] = {0,0,0}, sa[3] = {0,0,0}, sm[3] = {0,0,0}; long cnt = 0;
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        for (int s = 0; s < nstart; s++) {
          vector<int> st = random_plan(N, B);
          for (int strat = 0; strat < 3; strat++) {
            long wl = 0; vector<int> xf;
            double p = 100.0 * hillclimb(st, L, mv, strat, &wl, &xf) / denom;
            int active = 0, mx = 0; for (int v : xf) { if (v > 0) active++; if (v > mx) mx = v; }
            sp[strat] += p; sw[strat] += wl; sa[strat] += active; sm[strat] += (double)mx / B;
          }
          cnt++;
        }
      }
      for (int strat = 0; strat < 3; strat++)
        printf("%s,%d,%d,%.4f,%s,%.4f,%.2f,%.3f,%.4f\n", reg.c_str(), N, C, (double)C / N,
               sname[strat], sp[strat]/cnt, sw[strat]/cnt, sa[strat]/cnt, sm[strat]/cnt);
    }
    return 0;
  }
  // autocorr: a standard ruggedness metric. For each C, average over landscapes
  // the lag-1 autocorrelation of normalized profit along a random walk of
  // length wlen (Weinberger). Lower autocorrelation = more rugged.
  if (cmd == "autocorr") {
    if (argc < 9) { cerr << "autocorr needs: <N> <B> <regime> <Clist> <nland> <wlen> <seed>\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]); string reg = argv[4];
    vector<int> Cs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[6]), wlen = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    double denom = (double)B * B;
    printf("regime,N,C,C_over_N,autocorr\n");
    for (int C : Cs) {
      Moves mv = build_moves(N, C);
      double ac_sum = 0; long ac_n = 0;
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        vector<int> x = random_plan(N, B);
        vector<double> p; p.reserve(wlen + 1);
        p.push_back(100.0 * profit_full(x, L) / denom);
        for (int t = 0; t < wlen; t++) {
          vector<int> ok;
          for (int m = 0; m < (int)mv.from.size(); m++) if (x[mv.from[m]] >= 1) ok.push_back(m);
          if (ok.empty()) break;
          int m = ok[(int)(uniform_real_distribution<double>(0,1)(RNG) * ok.size())];
          x[mv.from[m]]--; x[mv.to[m]]++;
          p.push_back(100.0 * profit_full(x, L) / denom);
        }
        int n = (int)p.size(); if (n < 3) continue;
        double mean = 0; for (double z : p) mean += z; mean /= n;
        double num = 0, den = 0;
        for (int t = 0; t < n; t++) { double d = p[t] - mean; den += d * d; if (t < n-1) num += d * (p[t+1] - mean); }
        if (den > 0) { ac_sum += num / den; ac_n++; }
      }
      printf("%s,%d,%d,%.4f,%.4f\n", reg.c_str(), N, C, (double)C / N, ac_n ? ac_sum/ac_n : 0.0);
    }
    return 0;
  }
  // lon: a local-optima-network size proxy. For each C, the mean number of
  // *distinct* local optima reached by Steepest Ascent from `nstart` random
  // starts on a landscape (distinct final allocations), averaged over
  // landscapes. Falls toward 1 as the landscape smooths.
  if (cmd == "lon") {
    if (argc < 9) { cerr << "lon needs: <N> <B> <regime> <Clist> <nland> <nstart> <seed>\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]); string reg = argv[4];
    vector<int> Cs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[6]), nstart = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    printf("regime,N,C,C_over_N,mean_distinct_optima\n");
    for (int C : Cs) {
      Moves mv = build_moves(N, C);
      double sum = 0; long cnt = 0;
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        set<vector<int>> opt;
        for (int s = 0; s < nstart; s++) {
          vector<int> st = random_plan(N, B), xf;
          hillclimb(st, L, mv, 0, nullptr, &xf);
          opt.insert(xf);
        }
        sum += (double)opt.size(); cnt++;
      }
      printf("%s,%d,%d,%.4f,%.4f\n", reg.c_str(), N, C, (double)C / N, sum / cnt);
    }
    return 0;
  }
  // qsweep: fix (N, C) and sweep the returns-to-scale level r of the tuned
  // regime, reporting mean profit per strategy and cross-strategy SD.
  if (cmd == "qsweep") {
    if (argc < 9) { cerr << "qsweep needs 7 args\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]), C = atoi(argv[4]);
    vector<double> Rs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Rs.push_back(atof(t.c_str())); }
    int nland = atoi(argv[6]), nstart = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    double denom = (double)B * B;
    Moves mv = build_moves(N, C);
    printf("r,N,C,C_over_N,mean_SA,sd_SA_within,mean_LA,mean_MA,cross_sd\n");
    for (double r : Rs) {
      vector<double> SA, LA, MA;
      double cross_sum = 0; long ncross = 0, nwithin = 0; double within_sum = 0;
      for (int l = 0; l < nland; l++) {
        Landscape L = draw_tuned(N, r);
        vector<double> saL(nstart);
        for (int s = 0; s < nstart; s++) {
          vector<int> st = random_plan(N, B);
          double a = 100.0 * hillclimb(st, L, mv, 0) / denom;
          double b = 100.0 * hillclimb(st, L, mv, 1) / denom;
          double m = 100.0 * hillclimb(st, L, mv, 2) / denom;
          SA.push_back(a); LA.push_back(b); MA.push_back(m); saL[s] = a;
          double mu = (a + b + m) / 3.0;
          cross_sum += sqrt(((a-mu)*(a-mu) + (b-mu)*(b-mu) + (m-mu)*(m-mu)) / 2.0); ncross++;
        }
        if (nstart > 1) { double m = 0; for (double z : saL) m += z; m /= nstart;
          double v = 0; for (double z : saL) v += (z - m) * (z - m); within_sum += sqrt(v / (nstart - 1)); nwithin++; }
      }
      auto mean = [](vector<double>& v) { double s = 0; for (double z : v) s += z; return s / v.size(); };
      printf("%.3f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
             r, N, C, (double)C / N, mean(SA), nwithin ? within_sum / nwithin : 0.0,
             mean(LA), mean(MA), cross_sum / ncross);
    }
    return 0;
  }
  // enumerate: exhaustive ground truth on small instances. For each C, over
  // nland landscapes, enumerate every composition of B into N to find the exact
  // global optimum and the exact count of local optima (same directed-ring
  // operator as the main sweeps), then hill-climb from nstart random starts and
  // report attained profit as a fraction of the worst-to-best range. Feeds the
  // small-instance validation (exact optima counts, fraction-of-range, crossover
  // where ground truth is known). Output ratios are in [0,1]; cross_sd is the
  // dispersion across the three rules of that ratio.
  if (cmd == "enumerate") {
    if (argc < 9) { cerr << "enumerate needs: <N> <B> <regime> <Clist> <nland> <nstart> <seed>\n"; return 1; }
    int N = atoi(argv[2]), B = atoi(argv[3]); string reg = argv[4];
    vector<int> Cs; { stringstream ss(argv[5]); string t; while (getline(ss, t, ',')) Cs.push_back(atoi(t.c_str())); }
    int nland = atoi(argv[6]), nstart = atoi(argv[7]);
    RNG.seed(strtoull(argv[8], nullptr, 10));
    printf("regime,N,B,C,C_over_N,npoints,mean_nopt,frac_opt,ratio_SA,ratio_LA,ratio_MA,cross_sd\n");
    for (int C : Cs) {
      Moves mv = build_moves(N, C);
      int M = (int)mv.from.size();
      double sum_nopt = 0; long npoints = 0;
      double sR[3] = {0,0,0}; double sumcross = 0; long ncross = 0;
      vector<int> x(N, 0); vector<double> h(N, 0.0);
      for (int l = 0; l < nland; l++) {
        Landscape L = draw(reg, N);
        double gmax = -1e300, gmin = 1e300; long nopt = 0, npts = 0;
        // recursive enumeration of every composition of B into N
        function<void(int,int)> rec = [&](int idx, int rem) {
          if (idx == N - 1) {
            x[N-1] = rem;
            for (int i = 0; i < N; i++) h[i] = 0.0;            // h = cross * x
            for (int i = 0; i < N; i++) { double xi = x[i]; if (xi == 0) continue;
              const double* row = &L.cross[(size_t)i * N];
              for (int j = 0; j < N; j++) h[j] += row[j] * xi; }
            double p = 0.0;                                    // profit = x'Qx
            for (int i = 0; i < N; i++) p += L.sq[i] * (double)x[i] * x[i];
            for (int i = 0; i < N; i++) p += 0.5 * x[i] * h[i];
            if (p > gmax) gmax = p; if (p < gmin) gmin = p;
            bool improving = false;                            // local-optimum test
            for (int m = 0; m < M; m++) {
              int a = mv.from[m]; if (x[a] < 1) continue; int b = mv.to[m];
              double d = L.sq[a] * (1 - 2.0 * x[a]) + L.sq[b] * (1 + 2.0 * x[b])
                       + (h[b] - h[a]) - L.cij(a, b);
              if (d > 0) { improving = true; break; }
            }
            if (!improving) nopt++;
            npts++;
            return;
          }
          for (int v = 0; v <= rem; v++) { x[idx] = v; rec(idx + 1, rem - v); }
        };
        rec(0, B);
        npoints = npts; sum_nopt += nopt;
        if (gmax > gmin) {
          for (int s = 0; s < nstart; s++) {
            vector<int> st = random_plan(N, B);
            double r[3];
            for (int strat = 0; strat < 3; strat++) {
              double pf = hillclimb(st, L, mv, strat);
              r[strat] = (pf - gmin) / (gmax - gmin); sR[strat] += r[strat];
            }
            double mu = (r[0] + r[1] + r[2]) / 3.0;
            sumcross += sqrt(((r[0]-mu)*(r[0]-mu) + (r[1]-mu)*(r[1]-mu) + (r[2]-mu)*(r[2]-mu)) / 2.0);
            ncross++;
          }
        }
      }
      double mean_nopt = sum_nopt / nland;
      printf("%s,%d,%d,%d,%.4f,%ld,%.4f,%.6g,%.4f,%.4f,%.4f,%.4f\n",
             reg.c_str(), N, B, C, (double)C / N, npoints, mean_nopt, (double)mean_nopt / npoints,
             ncross ? sR[0]/ncross : 0.0, ncross ? sR[1]/ncross : 0.0, ncross ? sR[2]/ncross : 0.0,
             ncross ? sumcross/ncross : 0.0);
    }
    return 0;
  }
  cerr << "unknown command: " << cmd << "\n"; return 1;
}
