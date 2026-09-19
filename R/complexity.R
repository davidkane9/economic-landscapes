# Figures for "The Geometry of Complexity".
#
# These read CSVs precomputed by the C++ engine (src/landscape.cpp, driven by
# sim/run_sweeps.sh) rather than running the simulation at render time, so the
# document builds quickly and the heavy computation is reproducible and
# separately depositable. Each row of a CSV is a (regime, N, C) cell with the
# mean and standard error of Steepest- and Least-Ascent normalised profit and
# of their paired difference (the "Least Ascent advantage", x100).

suppressPackageStartupMessages(library(ggplot2))

# Least-Ascent advantage vs C/N, one line per coefficient regime (N = 20),
# ordered by the convexity of the profit matrix Q.
fig_complexity_regime <- function(path = "sim/regime.csv") {
  rg <- read.csv(path)
  lab <- c(leontief = "Diminishing returns (concave Q)",
           uniform  = "Uniform U[-1,1] (indefinite Q)",
           allpos   = "All-positive U[0,1]",
           signed   = "Increasing returns (convex Q)")
  rg$regime <- factor(rg$regime, levels = names(lab), labels = lab)
  ggplot(rg, aes(C_over_N, gap, colour = regime, fill = regime,
                 shape = regime, linetype = regime)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_ribbon(aes(ymin = gap - 2 * gap_se, ymax = gap + 2 * gap_se),
                alpha = 0.13, colour = NA, linetype = 0) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
    labs(x = "Connection density  C / N",
         y = "Least Ascent advantage  (LA - SA, x100)",
         colour = NULL, fill = NULL, shape = NULL, linetype = NULL) +
    theme_trials()
}

# Least-Ascent advantage vs C/N for several N (uniform regime): the crossover
# collapses onto C/N.
fig_complexity_collapse <- function(path = "sim/transition.csv") {
  tr <- read.csv(path)
  tr$N <- factor(tr$N, levels = sort(unique(tr$N)))
  ggplot(tr, aes(C_over_N, gap, colour = N, fill = N)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_ribbon(aes(ymin = gap - 2 * gap_se, ymax = gap + 2 * gap_se),
                alpha = 0.13, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
    labs(x = "Connection density  C / N",
         y = "Least Ascent advantage  (LA - SA, x100)",
         colour = "Inputs N", fill = "Inputs N") +
    theme_trials()
}

# Average profit across the (N, C) plane (uniform regime): the easy/hard
# gradient, with the complexity ridge C = 0.45 N marked.
fig_complexity_plane <- function(path = "sim/plane.csv") {
  d <- read.csv(path)
  ggplot(d, aes(N, C, fill = mean_SA)) +
    geom_tile() +
    geom_abline(slope = 0.45, intercept = 0, linetype = "dashed",
                colour = "white", linewidth = 0.5) +
    scale_fill_viridis_c(option = "C", name = "Mean\nprofit") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Number of inputs  N", y = "Connections per input  C") +
    theme_trials() + theme(legend.position = "right")
}

# Cross-strategy dispersion vs C/N, one line per N: the complexity band,
# peaking near C/N = 0.45 for every N.
fig_complexity_band <- function(path = "sim/plane.csv") {
  d <- read.csv(path); d$Nf <- factor(d$N)
  ggplot(d, aes(C_over_N, cross_sd, colour = N, group = Nf)) +
    geom_vline(xintercept = 0.45, linetype = "dashed", colour = "grey55") +
    geom_line(alpha = 0.8, linewidth = 0.4) +
    scale_colour_viridis_c(option = "D", name = "N") +
    labs(x = "Connection density  C / N",
         y = "Cross-strategy dispersion  (x100)") +
    theme_trials() + theme(legend.position = "right")
}

# Cross-strategy dispersion for three neighbourhood operators (directed ring,
# symmetric ring, random partners), plotted against out-degree per input / N --
# the natural common scale, since the symmetric ring has twice the out-degree of
# the others at a given C/N. If degree alone set the peak the three curves would
# align; they do not (the reversible ring peaks much later), so wiring matters
# beyond raw degree. (N = 20, uniform regime, B = 50.)
fig_complexity_operator <- function(path = "sim/operator.csv") {
  d <- read.csv(path)
  lab <- c("0" = "Directed ring (baseline)", "1" = "Symmetric ring",
           "2" = "Random partners")
  d$op <- factor(lab[as.character(d$topo)], levels = lab)
  ggplot(d, aes(outdeg / N, cross_sd, colour = op, shape = op)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    labs(x = "Out-degree per input  /  N",
         y = "Cross-strategy dispersion  (x100)", colour = NULL, shape = NULL) +
    theme_trials() + theme(legend.position = "right")
}

# Small-instance validation table (N = 8, B = 20, uniform regime): exact
# local-optima counts and hill-climber performance as a fraction of the
# worst-to-best profit range, where the 888,030-point simplex is enumerated in
# full. Reads sim/enumerate.csv.
tbl_enumerate <- function(path = "sim/enumerate.csv") {
  d <- read.csv(path); d <- d[d$regime == "uniform" & d$N == 8, ]
  d <- d[order(d$C_over_N), ]
  data.frame(
    `C/N`              = sprintf("%.2f", d$C_over_N),
    `Local optima`     = formatC(round(d$mean_nopt), big.mark = ",", format = "d"),
    `Steepest (%)`     = sprintf("%.1f", 100 * d$ratio_SA),
    `Least (%)`        = sprintf("%.1f", 100 * d$ratio_LA),
    check.names = FALSE)
}

# The two-dial complexity map: cross-strategy dispersion over the plane of
# connection density C/N (structural) and own-input curvature r (economic, the
# diagonal level of Q). Complexity grows monotonically toward increasing
# returns; the white line traces the C/N at which it peaks for each r, the
# "border" between the too-hard and too-easy regimes, which drifts upward as
# returns increase. (N = 24, B = 50; precomputed by src/landscape.cpp.)
fig_complexity_qplane <- function(path = "sim/qplane.csv") {
  d <- read.csv(path)
  ridge <- do.call(rbind, by(d, d$r, function(s) s[which.max(s$cross_sd),
                                                    c("r", "C_over_N")]))
  ggplot(d, aes(C_over_N, r, fill = cross_sd)) +
    geom_tile() +
    geom_path(data = ridge, aes(C_over_N, r), inherit.aes = FALSE,
              colour = "white", linewidth = 0.8) +
    geom_point(data = ridge, aes(C_over_N, r), inherit.aes = FALSE,
               colour = "white", size = 1.3) +
    scale_fill_viridis_c(option = "C", trans = "sqrt",
                         name = "Cross-\nstrategy SD") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Connection density  C / N",
         y = "Own-input curvature  r   (concave < 0 < convex)") +
    theme_trials() + theme(legend.position = "right")
}
