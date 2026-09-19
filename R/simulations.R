# ==================================================================
# simulations.R
# ------------------------------------------------------------------
# Simulation engine for "Economic Landscapes."
#
# A firm has a fixed integer budget B that it allocates among N inputs.
# A "plan" is an integer vector x = (x_1, ..., x_N) with x_i >= 0 and
# sum(x_i) = B.  Profit is given by a quadratic ("landscape") function
#
#   Pi(x) = sum_i c_i x_i                    (linear terms)
#         + sum_i d_i x_i^2                  (squared terms)
#         + sum_{i<j} e_{ij} x_i x_j         (cross-product terms)
#
# whose N + N + N(N-1)/2 = N(N+3)/2 coefficients are drawn at random.
# A different draw of coefficients is a different "landscape."
#
# NEIGHBOURHOOD / OPERATOR.  A move takes one dollar from input i and
# gives it to a partner input.  With `connections` = k partners per
# input, input i may transfer a dollar to inputs i+1, ..., i+k (indices
# taken cyclically, mod N).  The neighbourhood of a plan is every plan
# reachable by one such transfer (the source input must have >= 1
# dollar to give).  This is the cyclic one-dollar-redistribution
# operator described in the paper, generalised to k connections.
#
# HILL-CLIMBING STRATEGIES.  Starting from a random plan, a firm
# surveys its neighbourhood and steps to a higher-profit neighbour,
# repeating until it reaches a local maximum (no neighbour is higher).
#   - Steepest Ascent (SA): step to the single highest neighbour.
#   - Median  Ascent (MA): of the neighbours above the current point,
#                          step to the median one.
#   - Least   Ascent (LA): step to the lowest neighbour that is still
#                          an improvement (the most gradual climb).
#
# NORMALISATION.  Profit is a quadratic form, so it scales as B^2; we divide
# by B^2 so that normalised profit is comparable across budgets (with order-one
# coefficients and the budget on a single input, x'Qx is of order B^2).
# Reported figures multiply the normalised profit by 100 (a harmless rescaling).
#
# Where sensible the code is vectorised over the neighbourhood (all
# neighbours of a plan are scored in one matrix multiply); the climb
# itself and the loop over landscapes are explicit.
# ==================================================================

suppressPackageStartupMessages(library(ggplot2))

# ------------------------------------------------------------------
# Coefficient draws ("regimes") used by the three tables.
# Each returns a list with components used by profit():
#   lin  : length-N vector of linear coefficients
#   sq   : length-N vector of squared-term coefficients
#   cross: N x N symmetric matrix of cross-product coefficients
#          (only i<j entries are independent; diagonal is 0)
# ------------------------------------------------------------------

# pairs() : the i<j index pairs, in column-major order.
.pair_index <- function(N) {
  ij <- which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
  ij[order(ij[, "col"], ij[, "row"]), , drop = FALSE]
}

# Regime 1 (Table 2): every coefficient ~ Uniform[-1, 1].
draw_coef_uniform <- function(N) {
  np <- N * (N - 1) / 2
  cross <- matrix(0, N, N)
  ij <- .pair_index(N)
  cross[ij] <- runif(np, -1, 1)
  cross <- cross + t(cross)
  list(lin = runif(N, -1, 1), sq = runif(N, -1, 1), cross = cross)
}

# Regime 2 (Table 3): linear and squared ~ U[0,1]; cross ~ U[-1,0].
draw_coef_signed <- function(N) {
  np <- N * (N - 1) / 2
  cross <- matrix(0, N, N)
  ij <- .pair_index(N)
  cross[ij] <- runif(np, -1, 0)
  cross <- cross + t(cross)
  list(lin = runif(N, 0, 1), sq = runif(N, 0, 1), cross = cross)
}

# Regime 3 (Table 4): a Leontief-style penalty.  Linear terms are 0;
# each squared coefficient is the sum of (N-1) draws from U[-1,0]
# (so it lies in [-(N-1), 0]); each cross coefficient is the sum of
# two draws from U[0,1] (so it lies in [0, 2]).
draw_coef_leontief <- function(N) {
  np <- N * (N - 1) / 2
  cross <- matrix(0, N, N)
  ij <- .pair_index(N)
  cross[ij] <- runif(np, 0, 1) + runif(np, 0, 1)
  cross <- cross + t(cross)
  sq <- rowSums(matrix(runif(N * (N - 1), -1, 0), N, N - 1))
  list(lin = rep(0, N), sq = sq, cross = cross)
}

# ------------------------------------------------------------------
# profit(): profit of one plan, or of a matrix of plans (one per row).
# coef is a list as returned by the draw_coef_* functions above.
# Vectorised: scores all rows of `x` in one pass.
# ------------------------------------------------------------------
profit <- function(x, coef) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  lin_term <- as.vector(x %*% coef$lin)
  sq_term  <- as.vector((x * x) %*% coef$sq)
  # sum_{i<j} e_ij x_i x_j = 0.5 * ( x' E x ) since E is symmetric, 0-diag
  cross_term <- 0.5 * rowSums((x %*% coef$cross) * x)
  lin_term + sq_term + cross_term
}

# ------------------------------------------------------------------
# Neighbourhood / operator.
# build_moves(N, connections): the list of (from, to) input pairs that
# define a one-dollar transfer.  With k connections, input i may give a
# dollar to i+1, ..., i+k (cyclic, 1-based).  Returns an M x 2 matrix.
# ------------------------------------------------------------------
build_moves <- function(N, connections) {
  k <- min(connections, N - 1L)
  from <- rep(seq_len(N), each = k)
  off  <- rep(seq_len(k), times = N)
  to   <- ((from - 1L + off) %% N) + 1L
  cbind(from = from, to = to)
}

# neighbours(): given a plan x and a move table, return the matrix of
# neighbour plans (one per row).  A move is legal only if the source
# input currently holds at least one dollar; illegal moves are dropped.
neighbours <- function(x, moves) {
  legal <- x[moves[, "from"]] >= 1
  mv <- moves[legal, , drop = FALSE]
  if (nrow(mv) == 0L) return(NULL)
  nb <- matrix(x, nrow = nrow(mv), ncol = length(x), byrow = TRUE)
  idx <- seq_len(nrow(mv))
  nb[cbind(idx, mv[, "from"])] <- nb[cbind(idx, mv[, "from"])] - 1L
  nb[cbind(idx, mv[, "to"])]   <- nb[cbind(idx, mv[, "to"])]   + 1L
  nb
}

# ------------------------------------------------------------------
# random_plan(): a random integer allocation of B dollars over N inputs
# (a uniform composition: drop B balls into N boxes).
# ------------------------------------------------------------------
random_plan <- function(N, B) {
  tabulate(sample.int(N, B, replace = TRUE), nbins = N)
}

# ------------------------------------------------------------------
# hillclimb(): climb from a starting plan to a local maximum under one
# of the three strategies.  Returns the profit of the local maximum.
#   strategy %in% c("SA", "MA", "LA")
# ------------------------------------------------------------------
hillclimb <- function(start, coef, moves, strategy = c("SA", "MA", "LA")) {
  strategy <- match.arg(strategy)
  x  <- start
  px <- profit(x, coef)
  repeat {
    nb <- neighbours(x, moves)
    if (is.null(nb)) break
    pn <- profit(nb, coef)
    up <- which(pn > px)
    if (length(up) == 0L) break          # local maximum reached
    o <- up[order(pn[up])]               # improving neighbours, ascending
    pick <- switch(strategy,
                   SA = o[length(o)],            # highest
                   LA = o[1L],                   # smallest improvement
                   MA = o[ceiling(length(o) / 2)]) # median improvement
    x  <- nb[pick, ]
    px <- pn[pick]
  }
  px
}

# ------------------------------------------------------------------
# climb_normalised(): one landscape, one strategy -> normalised profit
# (x100), the quantity tabulated in Tables 2-4.
# ------------------------------------------------------------------
climb_normalised <- function(coef, N, B, connections, strategy) {
  moves <- build_moves(N, connections)
  start <- random_plan(N, B)
  100 * hillclimb(start, coef, moves, strategy) / (B^2)
}

# ==================================================================
# TABLE 2-4: mean (and s.e.) normalised profit over many landscapes,
# for each strategy x connections cell, under a given coefficient regime.
# ==================================================================
strategy_table <- function(draw_coef, N = 20L, B = 50L,
                           connections = 1:5,
                           strategies = c("SA", "MA", "LA"),
                           n_land = 1000L) {
  out <- list()
  for (k in connections) {
    # one fresh landscape per replication; all strategies share it so the
    # comparison is paired (as in the paper, where the same landscapes
    # are climbed by every strategy).
    mat <- matrix(NA_real_, n_land, length(strategies),
                  dimnames = list(NULL, strategies))
    for (r in seq_len(n_land)) {
      coef <- draw_coef(N)
      for (s in strategies)
        mat[r, s] <- climb_normalised(coef, N, B, k, s)
    }
    for (s in strategies) {
      out[[length(out) + 1L]] <- data.frame(
        strategy    = s,
        connections = k,
        mean        = mean(mat[, s]),
        se          = sd(mat[, s]) / sqrt(n_land)
      )
    }
  }
  do.call(rbind, out)
}

# format_table(): turn a strategy_table() result into a wide
# "mean (se)" data frame with strategies as rows, connections as columns.
format_table <- function(tab, strategies = c("SA", "MA", "LA"),
                         labels = c(SA = "Steepest Ascent",
                                    MA = "Median Ascent",
                                    LA = "Least Ascent"),
                         digits = 1) {
  ks <- sort(unique(tab$connections))
  cell <- function(s, k) {
    row <- tab[tab$strategy == s & tab$connections == k, ]
    sprintf("%.*f (%.1f)", digits, row$mean, row$se)
  }
  df <- data.frame(Strategy = labels[strategies], stringsAsFactors = FALSE)
  for (k in ks) df[[as.character(k)]] <- vapply(strategies, cell, "", k = k)
  rownames(df) <- NULL
  df
}

# ==================================================================
# LOCAL OPTIMA: of `points_per` random plans on each of `n_land`
# landscapes, what fraction are local maxima?  (1,000 points x
# 2,000 landscapes -> 8 local maxima out of 2,000,000.)
# ==================================================================
count_local_optima <- function(N = 20L, B = 50L, connections = 1L,
                               n_land = 2000L, points_per = 1000L) {
  moves <- build_moves(N, connections)
  is_max <- 0L
  tested <- 0L
  for (l in seq_len(n_land)) {
    coef <- draw_coef_uniform(N)
    for (p in seq_len(points_per)) {
      x  <- random_plan(N, B)
      nb <- neighbours(x, moves)
      if (is.null(nb)) next
      if (all(profit(nb, coef) < profit(x, coef))) is_max <- is_max + 1L
      tested <- tested + 1L
    }
  }
  list(tested = tested, local_maxima = is_max,
       fraction = is_max / tested)
}

# ==================================================================
# COMPARATIVE STATICS sweeps (Figures 15-16 and the budget result).
# Each returns mean normalised profit (x100) of the Steepest-Ascent
# local maximum, with a standard error, as one parameter varies.
# ==================================================================

# Helper: mean +/- se of SA-climbed normalised profit over n_land draws.
.sa_profit <- function(draw_coef, N, B, connections, n_land) {
  v <- numeric(n_land)
  for (r in seq_len(n_land))
    v[r] <- climb_normalised(draw_coef(N), N, B, connections, "SA")
  c(mean = mean(v), se = sd(v) / sqrt(n_land))
}

# Profit vs number of INPUTS (Figure 15).  Budget and connections fixed.
sweep_inputs <- function(input_grid, B = 50L, connections = 1L,
                         n_land = 1000L, draw_coef = draw_coef_uniform) {
  do.call(rbind, lapply(input_grid, function(N) {
    s <- .sa_profit(draw_coef, N, B, connections, n_land)
    data.frame(inputs = N, mean = s["mean"], se = s["se"],
               row.names = NULL)
  }))
}

# Profit vs CONNECTIONS per input (Figure 16).  Inputs and budget fixed.
sweep_connections <- function(conn_grid, N = 20L, B = 50L,
                              n_land = 1000L, draw_coef = draw_coef_uniform) {
  do.call(rbind, lapply(conn_grid, function(k) {
    s <- .sa_profit(draw_coef, N, B, k, n_land)
    data.frame(connections = k, mean = s["mean"], se = s["se"],
               row.names = NULL)
  }))
}

# Profit vs BUDGET (the "budget" comparative-static: normalised profit is
# roughly flat as B grows, inputs and connections fixed).
sweep_budget <- function(budget_grid, N = 20L, connections = 1L,
                         n_land = 100L, draw_coef = draw_coef_uniform) {
  do.call(rbind, lapply(budget_grid, function(B) {
    s <- .sa_profit(draw_coef, N, B, connections, n_land)
    data.frame(budget = B, mean = s["mean"], se = s["se"],
               row.names = NULL)
  }))
}

# ==================================================================
# Combinatorial term: number of plans satisfying the budget constraint,
# C(N + B - 1, B).  Returned as a double (it overflows integers).
# ==================================================================
n_plans <- function(N, B) choose(N + B - 1, B)

# ------------------------------------------------------------------
# Shared plotting theme (matches the companion clinical-trials paper).
# ------------------------------------------------------------------
theme_trials <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          legend.position = "right")
}

# ==================================================================
# Stylized profit landscapes (Section 2.2).  Three conceptual 3-D
# terrain surfaces over the N = 3, B = 50 budget simplex, drawn with
# base graphics.  No profit function is assumed; the surfaces are
# illustrative -- a smooth single peak, a bimodal surface, and a
# rugged surface with 1,000+ local maxima but one visible global max.
# ==================================================================
.land_sm <- function(M) {                       # light separable blur
  f <- function(v) (c(v[1], head(v, -1)) + v + c(tail(v, -1), v[length(v)])) / 3
  M <- t(apply(M, 1, f)); apply(M, 2, f)
}
.land_norm <- function(M)
  (M - min(M, na.rm = TRUE)) / (max(M, na.rm = TRUE) - min(M, na.rm = TRUE))

.land_facets <- function(z, pal) {              # per-facet terrain colours
  nr <- nrow(z); nc <- ncol(z)
  zf <- (z[-nr, -nc] + z[-nr, -1] + z[-1, -nc] + z[-1, -1]) / 4
  matrix(pal[floor(.land_norm(zf) * (length(pal) - 1)) + 1], nrow = nr - 1)
}

.land_persp <- function(x, y, Z, main)
  persp(x, y, Z, theta = -35, phi = 26, expand = 0.5,
        col = .land_facets(Z, hcl.colors(140, "Terrain")), border = NA,
        shade = 0.45, ltheta = -30, lphi = 45, box = FALSE,
        main = main, cex.main = 1.0, r = 2.2)

.land_flag <- function(pmat, gx, gy, gz) {      # mark the global maximum
  top <- gz + 0.16
  lines(trans3d(c(gx, gx), c(gy, gy), c(gz, top), pmat), col = "grey25", lwd = 1.3)
  points(trans3d(gx, gy, top, pmat), pch = 21, bg = "firebrick",
         col = "white", cex = 1.3, lwd = 1)
  text(trans3d(gx, gy, top, pmat), labels = "global", pos = 3,
       cex = 0.8, col = "grey20")
}

fig_stylized_landscapes <- function(B = 50) {
  op <- par(mfrow = c(1, 3), mar = c(0, 0, 2, 0)); on.exit(par(op))
  n <- 150; x <- seq(0, B, length.out = n); y <- x; off <- outer(x, y, `+`) > B
  # (a) concave: single central peak
  Za <- outer(x, y, function(a, b) exp(-((a - B/3)^2 + (b - B/3)^2) / (2 * 11^2)))
  Za[off] <- NA; Za <- .land_norm(Za)
  ga <- which(Za == max(Za, na.rm = TRUE), arr.ind = TRUE)[1, ]
  pm <- .land_persp(x, y, Za, "(a) Concave")
  .land_flag(pm, x[ga[1]], y[ga[2]], Za[ga[1], ga[2]])
  # (b) bimodal: taller peak is the global
  Zb <- outer(x, y, function(a, b)
    1.00 * exp(-((a - 12)^2 + (b - 30)^2) / (2 * 6.5^2)) +
    0.82 * exp(-((a - 31)^2 + (b - 10)^2) / (2 * 6.5^2)))
  Zb[off] <- NA; Zb <- .land_norm(Zb)
  gb <- which(Zb == max(Zb, na.rm = TRUE), arr.ind = TRUE)[1, ]
  pm <- .land_persp(x, y, Zb, "(b) Bimodal")
  .land_flag(pm, x[gb[1]], y[gb[2]], Zb[gb[1], gb[2]])
  # (c) rugged: 1,000+ local maxima, one visible global
  set.seed(1998)
  n3 <- 300; x3 <- seq(0, B, length.out = n3); y3 <- x3
  noise <- .land_sm(.land_sm(matrix(rnorm(n3 * n3), n3)))
  dome  <- outer(x3, y3, function(a, b) exp(-((a - 19)^2 + (b - 15)^2) / (2 * 15^2)))
  Zc <- .land_norm(.land_norm(noise) * 0.75 + .land_norm(dome) * 0.95)
  Zc[outer(x3, y3, `+`) > B] <- NA; Zc <- .land_norm(Zc)
  gc <- which(Zc == max(Zc, na.rm = TRUE), arr.ind = TRUE)[1, ]
  pm <- .land_persp(x3, y3, Zc, "(c) Rugged")
  .land_flag(pm, x3[gc[1]], y3[gc[2]], Zc[gc[1], gc[2]])
}

# ------------------------------------------------------------------
# A single random realization of the quadratic profit function of
# Section 3.1, for N = 3 and B = 50, plotted over the budget simplex.
# The linear, squared, and cross-product coefficients are drawn
# independently from U(-1, 1); the global maximum is flagged.
# ------------------------------------------------------------------
fig_random_landscape <- function(B = 50, seed = 12) {
  set.seed(seed)
  cc <- runif(6, -1, 1)            # c1..c3 squared (diagonal of Q), c4..c6 cross
  x <- 0:B; y <- 0:B
  Z <- matrix(NA_real_, B + 1, B + 1)
  for (i in seq_along(x)) for (j in seq_along(y)) {
    x1 <- x[i]; x2 <- y[j]; x3 <- B - x1 - x2
    if (x3 < 0) next
    Z[i, j] <- cc[1]*x1^2 + cc[2]*x2^2 + cc[3]*x3^2 +
               cc[4]*x1*x2 + cc[5]*x1*x3 + cc[6]*x2*x3
  }
  Z <- .land_norm(Z)
  op <- par(mar = c(0, 0, 1, 0)); on.exit(par(op))
  pm <- .land_persp(x, y, Z, "")
  g <- which(Z == max(Z, na.rm = TRUE), arr.ind = TRUE)[1, ]
  .land_flag(pm, x[g[1]], y[g[2]], Z[g[1], g[2]])
}
