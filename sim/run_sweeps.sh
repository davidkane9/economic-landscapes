#!/bin/bash
# Reproduce every simulation CSV for "Economic Landscapes".
# Requires a C++17 compiler (set CXX to override; default c++). Writes
# sim/{transition,regime,budget,plane,operator,diag,lon,autocorr,enumerate,qplane}.csv.
# Profit is normalized by B^2. All jobs run in parallel; the seeds are fixed below.
set -e
cd "$(dirname "$0")/.."
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
${CXX:-c++} -O3 -std=c++17 -o "$T/landscape" src/landscape.cpp
BIN="$T/landscape"; NL=1000; NS=8
$BIN sweep 10 50 uniform $(seq -s, 1 9)                     $NL $NS 101 > $T/tr_N10.csv &
$BIN sweep 20 50 uniform $(seq -s, 1 19)                    $NL $NS 102 > $T/tr_N20.csv &
$BIN sweep 30 50 uniform 1,2,3,4,5,6,8,10,12,15,18,21,24,29 $NL $NS 103 > $T/tr_N30.csv &
$BIN sweep 40 50 uniform 1,2,3,4,6,8,10,13,16,20,26,32,39   $NL $NS 104 > $T/tr_N40.csv &
RG=1,2,3,4,6,8,10,12,14,16,18,19
$BIN sweep 20 50 allpos   $RG $NL $NS 201 > $T/rg_allpos.csv &
$BIN sweep 20 50 signed   $RG $NL $NS 202 > $T/rg_signed.csv &
$BIN sweep 20 50 leontief $RG $NL $NS 203 > $T/rg_leontief.csv &
# Budget-robustness: the complexity band (cross-strategy SD) across the (C/N)
# axis at several budgets, to show the transition is governed by C/N, not B.
for B in 30 50 80 100; do
  $BIN grid 20 $B uniform $(seq -s, 1 19) 800 $NS $((300 + B)) > $T/bud_$B.csv &
done
# The (N, C) plane (uniform regime, B = 50): mean profit and cross-strategy SD
# for every N from 4 to 40, C from 1 to N-1. Feeds fig-complexity-plane/-band.
for N in $(seq 4 40); do
  $BIN grid $N 50 uniform $(seq -s, 1 $((N - 1))) 800 $NS $((500 + N)) > $T/pl_$N.csv &
done
wait
head -1 $T/tr_N10.csv > sim/transition.csv
for f in N10 N20 N30 N40; do tail -n +2 $T/tr_$f.csv >> sim/transition.csv; done
head -1 $T/tr_N20.csv > sim/regime.csv
tail -n +2 $T/tr_N20.csv >> sim/regime.csv
for f in allpos signed leontief; do tail -n +2 $T/rg_$f.csv >> sim/regime.csv; done
head -1 $T/bud_50.csv | sed 's/^regime,N,/regime,N,B,/' > sim/budget.csv
for B in 30 50 80 100; do tail -n +2 $T/bud_$B.csv | sed "s/^\(uniform,[0-9]*\),/\1,$B,/"; done >> sim/budget.csv
head -1 $T/pl_4.csv > sim/plane.csv
for N in $(seq 4 40); do tail -n +2 $T/pl_$N.csv >> sim/plane.csv; done

# Referee-revision sweeps: operator robustness, walk-length/optima-location
# diagnostics, the local-optima-network ruggedness metric, and the (null)
# autocorrelation metric. These use the exact-h engine and so are slower.
for topo in 0 1 2; do $BIN optest $topo 20 50 uniform $(seq -s, 1 19) 1000 $NS $((740 + topo)) > $T/op_$topo.csv & done
for reg in leontief uniform allpos signed; do $BIN diag 20 50 $reg 1,2,3,5 500 $NS 860 > $T/dg_$reg.csv & done
$BIN lon      20 50 uniform 1,2,3,5,7,9,11,13,15,17,19 400 16  870 > $T/lon.csv &
$BIN autocorr 20 50 uniform 1,3,5,9,13,19             300 300 880 > $T/ac.csv &
wait
head -1 $T/op_0.csv > sim/operator.csv
for topo in 0 1 2; do tail -n +2 $T/op_$topo.csv >> sim/operator.csv; done
head -1 $T/dg_uniform.csv > sim/diag.csv
for reg in leontief uniform allpos signed; do tail -n +2 $T/dg_$reg.csv >> sim/diag.csv; done
cp $T/lon.csv sim/lon.csv; cp $T/ac.csv sim/autocorr.csv

# Small-instance exact enumeration: where the simplex is small enough to walk in
# full (N=8,B=20 -> 888,030 allocations; N=6,B=18 -> 33,649), compute ground
# truth -- exact global optimum, exact local-optima counts, and hill-climber
# performance as a fraction of the worst-to-best range -- across regimes and C.
$BIN enumerate 8 20 uniform  1,2,3,4,5,6,7 200 8 901 > $T/en_u8.csv &
$BIN enumerate 8 20 signed   1,2,3,4,5,6,7 200 8 904 > $T/en_s8.csv &
$BIN enumerate 8 20 leontief 1,2,3,4,5,6,7 200 8 903 > $T/en_l8.csv &
$BIN enumerate 6 18 uniform  1,2,3,4,5     400 8 902 > $T/en_u6.csv &
wait
head -1 $T/en_u8.csv > sim/enumerate.csv
for f in u8 s8 l8 u6; do tail -n +2 $T/en_$f.csv >> sim/enumerate.csv; done

# The two-dial map (own-input curvature r x connection density C, at N = 24): how
# the complex band migrates through the plane as the economics turn from
# diminishing (r<0) to increasing (r>0) returns. Feeds fig-complexity-qplane.
for C in $(seq 1 23); do
  $BIN qsweep 24 50 $C -2,-1.5,-1,-0.5,0,0.5,1,1.5 400 $NS $((600 + C)) > $T/qp_$C.csv &
done
wait
head -1 $T/qp_1.csv > sim/qplane.csv
for C in $(seq 1 23); do tail -n +2 $T/qp_$C.csv >> sim/qplane.csv; done
echo "wrote sim/{transition,regime,budget,plane,operator,diag,lon,autocorr,enumerate,qplane}.csv"
