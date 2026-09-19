# Fraction of negative eigenvalues of the profit form's curvature, by regime.
#
# The firm allocates a fixed budget, so it moves only within the constraint
# hyperplane sum_i x_i = B. The curvature that governs whether profit is locally
# concave (one peak, trivial search) or convex (rugged) is therefore that of the
# quadratic form Pi(x) = x' Q x restricted to the constraint's tangent space
# {v : sum_i v_i = 0}, i.e. the *reduced* (projected) Hessian Z' Q Z for an
# orthonormal basis Z of that subspace -- not the raw spectrum of Q. This script
# reports both, so the numbers quoted in the text and in the caption of
# fig-complexity-regime can be reproduced. (Q = diag(sq) + 0.5*cross, matching
# src/landscape.cpp and R/simulations.R.)
#
#   Rscript R/eigen_fractions.R

set.seed(2024)
N <- 20L
ndraw <- 8000L
U <- function(n, a, b) runif(n, a, b)

draw_Q <- function(reg, N) {
  sq <- numeric(N); cr <- matrix(0, N, N)
  off <- function(f) for (i in 1:(N - 1)) for (j in (i + 1):N) {
    v <- f(); cr[i, j] <<- v; cr[j, i] <<- v
  }
  if (reg == "uniform")       { sq <- U(N, -1, 1); off(function() U(1, -1, 1)) }
  else if (reg == "allpos")   { sq <- U(N,  0, 1); off(function() U(1,  0, 1)) }
  else if (reg == "signed")   { sq <- U(N,  0, 1); off(function() U(1, -1, 0)) }
  else if (reg == "leontief") { for (i in 1:N) sq[i] <- sum(U(N - 1, -1, 0))
                                off(function() U(1, 0, 1) + U(1, 0, 1)) }
  diag(sq) + 0.5 * cr
}

# Orthonormal basis Z (N x (N-1)) for the tangent space {v : sum v = 0}.
P <- diag(N) - matrix(1 / N, N, N)
ez <- eigen(P, symmetric = TRUE)
Z <- ez$vectors[, ez$values > 0.5]

cat(sprintf("%-10s  full-Q   projected\n", "regime"))
for (reg in c("leontief", "uniform", "allpos", "signed")) {
  ff <- fp <- numeric(ndraw)
  for (d in seq_len(ndraw)) {
    Q <- draw_Q(reg, N)
    ff[d] <- mean(eigen(Q, symmetric = TRUE, only.values = TRUE)$values < 0)
    fp[d] <- mean(eigen(t(Z) %*% Q %*% Z, symmetric = TRUE, only.values = TRUE)$values < 0)
  }
  cat(sprintf("%-10s  %.3f    %.3f\n", reg, mean(ff), mean(fp)))
}
