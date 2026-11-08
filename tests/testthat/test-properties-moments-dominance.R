# Property tests for NNS.moments and the stochastic-dominance verdicts.
#
# Each test_that names the Lean theorem it mirrors in nns-proofs/NNSProofs:
#   Moments.lean       NNS.moments as partial moments about the mean.
#   Integration.lean   LPM_{d+1} = (d+1) * integral of LPM_d, the lifting
#                      FSD => SSD => TSD, knot sufficiency (LPM_one_le_of_le_at_knots,
#                      LPM_two_le_of_le_at_knots) and the tsdCex / tsdVertexCex
#                      counterexamples.
#   PrefixBackend.lean the >32-target prefix backend equals the direct kernel
#                      (prefixLPM_eq_LPM, prefixUPM_eq_UPM, usesPrefix_iff).
#   DominanceAlgo.lean the C++ verdict kernel sd_dom_pair is exact: sdDomPair_fsd_iff,
#                      sdDomPair_ssd_iff, sdDomPair_tsd_iff, sdDomPair_omega_iff
#                      (continuous FSD), sdDomPair_eq_zero_or_one, centring and
#                      translation (sdDomPairCentred_eq, sdDomPairCentred_omega_eq,
#                      LPMDominates_shift), tsdKnotOnly_vertexCex.
# Tests without a Lean theorem (the bidirectional string labels) cite SPEC.md.
#
# Data regimes:
#   exact  : small integer / dyadic data with n a power of 2, so sums and
#            divisions are exact in double; expect_identical.
#   random : fixed seed, rnorm / rexp / rt(df = 3), rounded to force ties;
#            absolute tolerance (n + k) * eps * scale.
#   real   : datasets shipped with R, targets at mean, median, quantiles and
#            exact data points (ties).
#
# Tests that fail on the baseline build are kept as the specification and
# carry a "# fails on 8898017: see BUGS.md #<id>" comment.

local_edition(3)

# ---- helpers ---------------------------------------------------------------

eps <- .Machine$double.eps

# Absolute scale of a degree-k moment about t: mean(|x - t|^k).
abs_mom <- function(x, k, t = mean(x)) mean(abs(x - t)^k)

# Naive summation of n terms of degree k: error below (n + k) * eps * scale.
mom_tol <- function(x, k, t = mean(x)) 4 * (length(x) + k + 2) * eps * abs_mom(x, k, t)

expect_within <- function(lhs, rhs, tol, label = NULL) {
  expect_lte(max(abs(lhs - rhs) - tol), 0, label = label)
}

# Partial-moment assembly of NNS.moments (momentsPM in Moments.lean).
moments_pm <- function(x, population = TRUE) {
  n <- length(x)
  mu <- mean(x)
  V <- LPM(2, mu, x) + UPM(2, mu, x)
  S <- UPM(3, mu, x) - LPM(3, mu, x)
  K <- UPM(4, mu, x) + LPM(4, mu, x)
  if (population) {
    return(list(mean = mu, variance = V, skewness = S / V^1.5, kurtosis = K / V^2 - 3))
  }
  s2 <- V * n / (n - 1)
  list(mean = mu, variance = s2,
       skewness = n / ((n - 1) * (n - 2)) * (n * S) / s2^1.5,
       kurtosis = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * (n * K) / s2^2 -
         3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}

# Absolute tolerances for the standardized moments (both sides carry the
# summation error of the degree-k sum, relative to V^(k/2)).
moments_tol <- function(x) {
  V <- abs_mom(x, 2)
  c(mean = 4 * eps * max(abs(x)),
    variance = mom_tol(x, 2),
    skewness = 4 * mom_tol(x, 3) / V^1.5 + 8 * length(x) * eps * abs_mom(x, 3) / V^1.5,
    kurtosis = 4 * mom_tol(x, 4) / V^2 + 8 * length(x) * eps * abs_mom(x, 4) / V^2)
}

expect_moments_equal <- function(got, want, tol, label) {
  for (nm in c("mean", "variance", "skewness", "kurtosis")) {
    expect_within(got[[nm]], want[[nm]], tol[[nm]], label = paste(label, nm))
  }
}

targets_of <- function(x) {
  c(mean(x), median(x), quantile(x, c(0.05, 0.25, 0.75, 0.95), names = FALSE),
    x[1], min(x), max(x))
}

# Exact quadrature of a piecewise polynomial on consecutive knots k[1] < k[2] < ...
# (the LPM of degree d is a polynomial of degree d between observations).
# Returns cumulative integrals from k[1] to each knot.
cum_step <- function(k, f) c(0, cumsum(head(f(k), -1) * diff(k)))            # d = 0, right-continuous
cum_trap <- function(k, f) { v <- f(k); c(0, cumsum((head(v, -1) + v[-1]) / 2 * diff(k))) }
cum_simpson_half <- function(k, f) {  # returns 3 * integral, kept dyadic: h * (fa + 4 fm + fb) / 2
  v <- f(k); m <- f((head(k, -1) + k[-1]) / 2)
  c(0, cumsum(diff(k) * (head(v, -1) + 4 * m + v[-1]) / 2))
}

# Knot grid for the integration identities: one point below every observation,
# the sorted distinct observations, one point above. refine = TRUE adds the
# midpoints (still a partition on whose cells LPM_d is one polynomial), which
# also pushes small real samples past 32 targets onto the prefix backend.
knot_grid <- function(x, refine = FALSE) {
  u <- sort(unique(x))
  k <- c(u[1] - 1, u, u[length(u)] + 1)
  if (refine) k <- sort(c(k, (head(k, -1) + k[-1]) / 2))
  k
}

# The R verdict for SSD on the merged knots (SPEC "sd_dom_pair, SSD"), as a
# function of the knot differences d = LPM1_x - LPM1_y. NA = inside the band.
ssd_expected <- function(d, tol) {
  if (max(d) > tol) return(0L)
  if (all(d <= 0) && any(d < -tol)) return(1L)
  NA_integer_
}

# Implementation-free LPM(1, t, x) = mean(max(t - x, 0)), used as the SSD oracle
# (identical() to LPM(1, t, x) on the integer pairs, but shares no code with NNS).
lpm1_ref <- function(ts, x) vapply(ts, function(t) mean(pmax(t - x, 0)), numeric(1))

# Exact TSD criterion of sdDomPair_tsd_iff for integer samples. Sufficiency is
# LPM_two_le_of_le_at_knots (knots, interior vertices s* = b / a, mean gate);
# strictness is checked at the knots only, which LPM_two_eq_of_eq_at_knots
# makes exact. All quantities are integers scaled by mx * my, so the test is
# exact. tsd_parts returns the three components separately so tests can count
# which one decided.
tsd_parts <- function(x, y) {
  mx <- length(x); my <- length(y)
  k <- sort(unique(c(x, y)))
  A <- vapply(k, function(t) sum(x <= t) * my - sum(y <= t) * mx, numeric(1))
  B <- vapply(k, function(t) sum(x[x <= t]) * my - sum(y[y <= t]) * mx, numeric(1))
  C <- vapply(k, function(t) sum(x[x <= t]^2) * my - sum(y[y <= t]^2) * mx, numeric(1))
  q_knot <- A * k^2 - 2 * B * k + C
  lo <- head(k, -1); hi <- k[-1]; a <- head(A, -1); b <- head(B, -1); cc <- head(C, -1)
  inside <- a < 0 & lo * a > b & b > hi * a
  list(knots_ok = all(q_knot <= 0),
       vertex_bad = any(inside & (cc * a - b^2 < 0)),
       mean_ok = sum(x) * my >= sum(y) * mx,
       differs = any(q_knot != 0))
}
tsd_exact <- function(x, y) {
  p <- tsd_parts(x, y)
  as.integer(p$knots_ok && !p$vertex_bad && p$mean_ok && p$differs)
}

# Exact continuous-FSD (Omega) criterion of sdDomPair_omega_iff for integer
# samples: R_x(t) <= R_y(t) for all t, with R = LPM.ratio(1, t, .) =
# L1 / (L1 + U1) (0/0 -> k/m, SPEC "sd_dom_pair, continuous FSD"), and the
# curves differ somewhere. Independent of NNS: counts and prefix sums only.
# On a knot cell both denominators are positive, so R_x <= R_y iff
# g(t) = L1x * U1y - L1y * U1x <= 0, a quadratic al t^2 + be t + ga; it is
# checked at the knots and at an interior vertex. A quadratic that vanishes
# at both knots and at the midpoint is zero, so "differs" is knots + midpoints.
omega_parts <- function(x, y) {
  mx <- length(x); my <- length(y); sx <- sum(x); sy <- sum(y)
  k <- sort(unique(c(x, y)))
  ax <- vapply(k, function(t) sum(x <= t), numeric(1))
  px <- vapply(k, function(t) sum(x[x <= t]), numeric(1))
  ay <- vapply(k, function(t) sum(y <= t), numeric(1))
  py <- vapply(k, function(t) sum(y[y <= t]), numeric(1))
  lx <- ax * k - px; dx <- sx - 2 * px + (2 * ax - mx) * k   # m * L1, m * (L1 + U1)
  ly <- ay * k - py; dy <- sy - 2 * py + (2 * ay - my) * k
  nx <- ifelse(dx == 0, ax, lx); dx <- ifelse(dx == 0, mx, dx)
  ny <- ifelse(dy == 0, ay, ly); dy <- ifelse(dy == 0, my, dy)
  g_knot <- nx * dy - ny * dx                                 # sign(R_x - R_y)
  i <- seq_len(length(k) - 1)
  a1 <- ax[i]; p1 <- px[i]; a2 <- ay[i]; p2 <- py[i]; lo <- k[i]; hi <- k[i + 1]
  al <- -a1 * (my - a2) + a2 * (mx - a1)
  be <- a1 * (sy - p2) + p1 * (my - a2) - a2 * (sx - p1) - p2 * (mx - a1)
  ga <- -p1 * (sy - p2) + p2 * (sx - p1)
  inside <- al < 0 & -be > 2 * al * hi & -be < 2 * al * lo   # -be / (2 al) in (lo, hi)
  mid <- lo + hi
  list(knots_ok = all(g_knot <= 0),
       vertex_bad = any(inside & (4 * al * ga - be^2 < 0)),
       knot_differs = any(g_knot != 0),
       mid_differs = any(al * mid^2 + 2 * be * mid + 4 * ga != 0))
}
omega_exact <- function(x, y) {
  p <- omega_parts(x, y)
  as.integer(p$knots_ok && !p$vertex_bad && (p$knot_differs || p$mid_differs))
}

lpm_each <- function(d, ts, x) vapply(ts, function(t) LPM(d, t, x), numeric(1))
upm_each <- function(d, ts, x) vapply(ts, function(t) UPM(d, t, x), numeric(1))

# ---- data ------------------------------------------------------------------

# exact regime: n = 8, integers with ties at 1 and 4; mean 3
x_int <- c(0, 1, 1, 3, 4, 4, 4, 7)
# dyadic sample with mean exactly 0 and a tie at 0.25
x_dy <- c(-2, -1, -0.5, 0, 0.5, 0.5, 1.5, 3) - 0.25
t_dy <- c(-3, -1.25, -0.25, 0, 0.25, 1, 2.75, 4)

set.seed(20261108)
random_samples <- list(
  normal = round(rnorm(400), 1),
  expo   = round(rexp(300, 2), 2),
  t3     = round(rt(500, df = 3), 1)
)

real_samples <- list(
  dax       = as.numeric(diff(log(EuStockMarkets[, "DAX"]))),
  eruptions = faithful$eruptions,
  ozone     = as.numeric(na.omit(airquality$Ozone)),
  mag       = quakes$mag,
  mpg       = mtcars$mpg
)

all_samples <- c(random_samples, real_samples)

# Random small integer pairs for the exact dominance tests (values 0..6, sizes 1..6).
set.seed(7)
int_pairs <- replicate(400, list(x = sample(0:6, sample(1:6, 1), replace = TRUE),
                                 y = sample(0:6, sample(1:6, 1), replace = TRUE)),
                       simplify = FALSE)

# Pairs that are FSD by construction: y_i <= x_i, some strictly, ties kept.
set.seed(11)
fsd_pairs <- replicate(100, {
  x <- round(rt(sample(5:40, 1), df = 3), 1)
  list(x = x, y = round(x - rexp(length(x)) * rbinom(length(x), 1, 0.5), 1))
}, simplify = FALSE)

verdicts <- function(p) c(
  fsd = NNS.FSD.uni(p$x, p$y, type = "discrete"),
  ssd = NNS.SSD.uni(p$x, p$y),
  tsd = NNS.TSD.uni(p$x, p$y)
)
int_verdicts <- t(vapply(int_pairs, verdicts, integer(3)))
int_verdicts_rev <- t(vapply(int_pairs, function(p) verdicts(list(x = p$y, y = p$x)), integer(3)))
fsd_verdicts <- t(vapply(fsd_pairs, verdicts, integer(3)))
int_omega <- vapply(int_pairs, function(p) NNS.FSD.uni(p$x, p$y, type = "continuous"), integer(1))
int_omega_rev <- vapply(int_pairs, function(p) NNS.FSD.uni(p$y, p$x, type = "continuous"), integer(1))

# Pairs on which the interior TSD vertex check fires (some knot cell has its
# vertex s* = b / a inside with a positive LPM2 difference). On int_pairs the
# vertex never decides a verdict, so reject-sample a wider regime (values 0..5,
# sizes 2..15) until 25 pairs pass the knots, the min gate and the mean gate and
# are rejected by the vertex alone; keep up to 200 vertex-firing pairs overall.
# The 1e5 iteration cap is a safety net only: seed 23 reaches 25 in about 2.5e4.
set.seed(23)
vertex_pairs <- list(); n_vertex_decided <- 0L; iter <- 0L
while (n_vertex_decided < 25L && iter < 100000L) {
  iter <- iter + 1L
  x <- sample(0:5, sample(2:15, 1), replace = TRUE)
  y <- sample(0:5, sample(2:15, 1), replace = TRUE)
  pt <- tsd_parts(x, y)
  if (!pt$vertex_bad) next
  decided <- pt$knots_ok && pt$mean_ok && pt$differs && min(x) >= min(y)
  if (decided || length(vertex_pairs) < 200L) {
    vertex_pairs[[length(vertex_pairs) + 1L]] <- list(x = x, y = y, decided = decided)
    n_vertex_decided <- n_vertex_decided + decided
  }
}
rm(x, y, pt, iter)

# ============================================================================
# Moments.lean
# ============================================================================

test_that("centralMoment_two_eq / three_eq / four_eq: mean((x-t)^k) = UPM(k) +- LPM(k), exact on dyadic data", {
  for (t in t_dy) {
    expect_identical(LPM(2, t, x_dy) + UPM(2, t, x_dy), mean((x_dy - t)^2))
    expect_identical(UPM(3, t, x_dy) - LPM(3, t, x_dy), mean((x_dy - t)^3))
    expect_identical(UPM(4, t, x_dy) + LPM(4, t, x_dy), mean((x_dy - t)^4))
  }
})

test_that("centralMoment_two_eq / three_eq / four_eq: identity at mean, median, quantiles and ties", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    for (t in targets_of(x)) {
      expect_within(LPM(2, t, x) + UPM(2, t, x), mean((x - t)^2), mom_tol(x, 2, t), label = nm)
      expect_within(UPM(3, t, x) - LPM(3, t, x), mean((x - t)^3), mom_tol(x, 3, t), label = nm)
      expect_within(UPM(4, t, x) + LPM(4, t, x), mean((x - t)^4), mom_tol(x, 4, t), label = nm)
    }
  }
})

test_that("UPM_one_mean_eq_LPM_one_mean: the mean balances degree-1 partial moments", {
  expect_identical(UPM(1, mean(x_dy), x_dy), LPM(1, mean(x_dy), x_dy))
  expect_identical(UPM(1, mean(x_int), x_int), LPM(1, mean(x_int), x_int))
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    expect_within(UPM(1, mean(x), x), LPM(1, mean(x), x), mom_tol(x, 1), label = nm)
  }
})

test_that("popVariance_eq_PM: population variance is LPM(2, mu) + UPM(2, mu), exact on dyadic data", {
  m <- NNS.moments(x_dy)
  expect_identical(m$mean, 0)
  expect_identical(m$variance, LPM(2, 0, x_dy) + UPM(2, 0, x_dy))
  m <- NNS.moments(x_int)
  expect_identical(m$variance, LPM(2, 3, x_int) + UPM(2, 3, x_int))
})

test_that("moments_eq_momentsPM (population): NNS.moments equals the partial-moment assembly", {
  expect_moments_equal(NNS.moments(x_dy), moments_pm(x_dy), moments_tol(x_dy), "dyadic")
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    expect_moments_equal(NNS.moments(x), moments_pm(x), moments_tol(x), nm)
  }
})

test_that("moments_eq_momentsPM (sample): NNS.moments(population = FALSE) equals the PM assembly", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    expect_moments_equal(NNS.moments(x, population = FALSE), moments_pm(x, FALSE),
                         2 * moments_tol(x), nm)
  }
})

test_that("moments_probe: NNS.moments(c(2,4,4,7,9)) matches the exact rational values", {
  m <- NNS.moments(c(2, 4, 4, 7, 9))
  expect_identical(m$mean, 5.2)
  expect_equal(m$variance, 154 / 25, tolerance = 4 * eps)
  expect_equal(m$skewness, 306 / (77 * sqrt(154)), tolerance = 8 * eps)
  expect_equal(m$kurtosis, -15073 / 11858, tolerance = 8 * eps)
})

test_that("moments_probe_sample: NNS.moments(c(2,4,4,7,9), FALSE) matches the exact values", {
  m <- NNS.moments(c(2, 4, 4, 7, 9), population = FALSE)
  expect_identical(m$mean, 5.2)
  expect_equal(m$variance, 77 / 10, tolerance = 4 * eps)
  expect_equal(m$skewness, sqrt(20) / 3 * 306 / (77 * sqrt(154)), tolerance = 8 * eps)
  expect_equal(m$kurtosis, -6430 / 5929, tolerance = 8 * eps)
})

test_that("sampleVariance_eq: sample variance is sum((x - mu)^2) / (n - 1)", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    expect_within(NNS.moments(x, FALSE)$variance, sum((x - mean(x))^2) / (length(x) - 1),
                  mom_tol(x, 2), label = nm)
  }
})

test_that("sampleSkewness_eq_G1 / sampleKurtosis_eq_G2: Joanes-Gill G1 and G2 (n > 3)", {
  g1g2 <- function(x) {
    n <- length(x); z <- x - mean(x)
    m2 <- mean(z^2); m3 <- mean(z^3); m4 <- mean(z^4)
    g1 <- m3 / m2^1.5; g2 <- m4 / m2^2 - 3
    c(skewness = sqrt(n * (n - 1)) / (n - 2) * g1,
      kurtosis = (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6))
  }
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    m <- NNS.moments(x, FALSE)
    tol <- 2 * moments_tol(x)
    expect_within(m$skewness, g1g2(x)[["skewness"]], tol[["skewness"]], label = nm)
    expect_within(m$kurtosis, g1g2(x)[["kurtosis"]], tol[["kurtosis"]], label = nm)
  }
})

test_that("popVariance_eq_zero_iff: a constant sample has variance exactly 0 (skew, kurt NaN by R convention)", {
  for (v in c(3, -0.1, 1e6)) {
    m <- NNS.moments(rep(v, 7))
    expect_identical(m$variance, 0)
    expect_true(is.nan(m$skewness))
    expect_true(is.nan(m$kurtosis))
  }
  expect_gt(NNS.moments(c(3, 3, 3, 3.5))$variance, 0)
})

test_that("popVariance_eq_zero_iff (n = 1): population variance 0, sample variance NaN", {
  m <- NNS.moments(5)
  expect_identical(m$mean, 5)
  expect_identical(m$variance, 0)
  expect_true(is.nan(NNS.moments(5, FALSE)$variance))
})

test_that("centralMoment_three_mean_of_card_two: n = 2 has m3 = 0, so sample skewness is NaN", {
  expect_identical(UPM(3, 1.5, c(1, 2)) - LPM(3, 1.5, c(1, 2)), 0)
  expect_identical(NNS.moments(c(1, 2))$skewness, 0)
  expect_true(is.nan(NNS.moments(c(1, 2), FALSE)$skewness))
})

test_that("sampleKurtosis_eq_G2 guard (n > 3): n = 3 sample kurtosis is NaN", {
  m <- NNS.moments(c(1, 2, 3), FALSE)
  expect_identical(m$skewness, 0)
  expect_true(is.nan(m$kurtosis))
})

test_that("popSkewness_pos_iff / popSkewness_neg_iff: sign(skewness) = sign(UPM(3, mu) - LPM(3, mu))", {
  samples <- c(all_samples, list(right = c(0, 0, 0, 0, 0, 0, 0, 8), left = -c(0, 0, 0, 0, 0, 0, 0, 8)))
  for (nm in names(samples)) {
    x <- samples[[nm]]
    expect_identical(sign(NNS.moments(x)$skewness), sign(UPM(3, mean(x), x) - LPM(3, mean(x), x)),
                     label = nm)
  }
})

test_that("popSkewness_sq_sub_two_le_popKurtosis: Pearson's inequality kurtosis >= skewness^2 - 2", {
  samples <- c(all_samples, list(two_point = c(0, 0, 0, 1), bernoulli = rep(c(0, 1), 8)))
  for (nm in names(samples)) {
    m <- NNS.moments(samples[[nm]])
    expect_gte(m$kurtosis - (m$skewness^2 - 2), -64 * eps * (1 + m$skewness^2), label = nm)
  }
  # equality case: two-point distributions attain the bound
  m <- NNS.moments(c(0, 0, 0, 1))
  expect_equal(m$kurtosis, m$skewness^2 - 2, tolerance = 16 * eps)
})

test_that("popKurtosis_le / popSkewness_sq_le: kurtosis <= n - 3 and skewness^2 <= n - 1", {
  samples <- c(all_samples, list(spike = c(rep(0, 15), 1)))
  for (nm in names(samples)) {
    x <- samples[[nm]]; n <- length(x)
    m <- NNS.moments(x)
    expect_lte(m$kurtosis, (n - 3) * (1 + 64 * eps), label = nm)
    expect_lte(m$skewness^2, (n - 1) * (1 + 64 * eps), label = nm)
  }
  # a single spike of n = 16: skewness^2 = 196/15, kurtosis = 166/15 (exact values)
  m <- NNS.moments(c(rep(0, 15), 1))
  expect_equal(m$skewness^2, 196 / 15, tolerance = 16 * eps)
  expect_equal(m$kurtosis, 166 / 15, tolerance = 16 * eps)
})

test_that("popVariance_affine / popSkewness_affine / popKurtosis_affine: exact dyadic map x -> -2 x + 8", {
  m <- NNS.moments(x_dy)
  ma <- NNS.moments(-2 * x_dy + 8)
  expect_identical(ma$mean, -2 * m$mean + 8)
  expect_identical(ma$variance, 4 * m$variance)
  expect_identical(ma$skewness, -m$skewness)
  expect_identical(ma$kurtosis, m$kurtosis)
})

test_that("popVariance_affine / popSkewness_affine / popKurtosis_affine: random and real data", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    for (a in c(-3, 0.5, 7)) {
      y <- a * x + 2.5
      m <- NNS.moments(x); ma <- NNS.moments(y)
      tol <- moments_tol(y)
      expect_within(ma$variance, a^2 * m$variance, tol[["variance"]] + a^2 * moments_tol(x)[["variance"]], label = nm)
      expect_within(ma$skewness, sign(a) * m$skewness, 2 * tol[["skewness"]], label = nm)
      expect_within(ma$kurtosis, m$kurtosis, 2 * tol[["kurtosis"]], label = nm)
    }
  }
})

test_that("sampleSkewness_affine / sampleKurtosis_affine: sample branch under x -> a x + c", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    m <- NNS.moments(x, FALSE); ma <- NNS.moments(-4 * x + 1, FALSE)
    tol <- 2 * moments_tol(-4 * x + 1)
    expect_within(ma$skewness, -m$skewness, 2 * tol[["skewness"]], label = nm)
    expect_within(ma$kurtosis, m$kurtosis, 2 * tol[["kurtosis"]], label = nm)
  }
})

# fails on 8898017: see BUGS.md #12 (same underflow family; the NNS.moments
# case is the JOURNAL.md "scale 1e-200" finding, no dedicated BUGS.md id):
# mean(z^3), mean(z^4) underflow to 0 and skew/kurt become NaN.
test_that("popSkewness_affine / popKurtosis_affine: invariance under scale a = 1e-200", {
  for (nm in c("normal", "eruptions", "mpg")) {
    x <- all_samples[[nm]]
    m <- NNS.moments(x); ms <- NNS.moments(1e-200 * x)
    expect_equal(ms$skewness, m$skewness, tolerance = 1e-12, label = nm)
    expect_equal(ms$kurtosis, m$kurtosis, tolerance = 1e-12, label = nm)
  }
})

# ============================================================================
# Integration.lean
# ============================================================================

test_that("LPM_one_eq_integral_LPM_zero: LPM(1, t) = integral of the ECDF LPM(0, .), exact with ties", {
  k <- knot_grid(x_int)
  expect_identical(LPM(1, k, x_int), cum_step(k, function(s) LPM(0, s, x_int)))
})

test_that("LPM_succ_eq_integral d = 1: LPM(2, t) = 2 * integral of LPM(1, .), exact", {
  k <- knot_grid(x_int)
  expect_identical(LPM(2, k, x_int), 2 * cum_trap(k, function(s) LPM(1, s, x_int)))
})

test_that("LPM_succ_eq_integral d = 2: LPM(3, t) = 3 * integral of LPM(2, .), exact (Simpson)", {
  k <- knot_grid(x_int)
  expect_identical(LPM(3, k, x_int), cum_simpson_half(k, function(s) LPM(2, s, x_int)))
})

test_that("LPM_two_eq_iterated_integral: LPM(2, t) = 2 * double integral of the ECDF (factor 2)", {
  k <- knot_grid(x_int)
  inner <- cum_step(k, function(s) LPM(0, s, x_int))  # piecewise linear, exact at knots
  outer <- c(0, cumsum((head(inner, -1) + inner[-1]) / 2 * diff(k)))
  expect_identical(LPM(2, k, x_int), 2 * outer)
})

test_that("UPM_one_eq_integral_UPM_zero / UPM_succ_eq_integral: upper analogues, exact", {
  k <- knot_grid(x_int)
  rev_cum <- function(v) rev(cumsum(rev(v)))
  h <- diff(k)
  u0 <- UPM(0, head(k, -1), x_int)          # strict upper side: constant on [k_i, k_i+1)
  expect_identical(UPM(1, k, x_int), c(rev_cum(u0 * h), 0))
  u1 <- UPM(1, k, x_int)
  expect_identical(UPM(2, k, x_int), 2 * c(rev_cum((head(u1, -1) + u1[-1]) / 2 * h), 0))
})

test_that("LPM_succ_eq_integral d = 0, 1, 2: random and real data via the >32-target prefix path", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    k <- knot_grid(x, refine = TRUE)
    expect_gt(length(k), 32)
    scale <- max(abs(x - mean(x)))
    tol <- function(d) 8 * (length(x) + length(k)) * eps * (scale + diff(range(k)))^(d + 1)
    expect_within(LPM(1, k, x), cum_step(k, function(s) LPM(0, s, x)), tol(0), label = nm)
    expect_within(LPM(2, k, x), 2 * cum_trap(k, function(s) LPM(1, s, x)), tol(1), label = nm)
    expect_within(LPM(3, k, x), cum_simpson_half(k, function(s) LPM(2, s, x)), tol(2), label = nm)
  }
})

test_that("prefixLPM_eq_LPM / prefixUPM_eq_UPM (workerLPM_eq_rLPMVec / workerUPM_eq_rUPMVec): >32-target prefix backend agrees with scalar targets", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- sort(c(targets_of(x), seq(min(x), max(x), length.out = 40)))
    for (d in 0:3) {
      tol <- 8 * (length(x) + d) * eps * max(abs_mom(x, d, min(x)), abs_mom(x, d, max(x)))
      expect_within(LPM(d, ts, x), lpm_each(d, ts, x), tol, label = sprintf("%s LPM d=%d", nm, d))
      expect_within(UPM(d, ts, x), upm_each(d, ts, x), tol, label = sprintf("%s UPM d=%d", nm, d))
    }
  }
})

test_that("usesPrefix_iff: 32 targets stay on the direct kernel (identical to scalar calls), 33 switch to the prefix backend", {
  # NNS_DIRECT_PATH_MAX_TARGETS = 32: a 32-target call runs the same scalar
  # kernel per target, so it is bit-identical; at 33 targets the prefix backend
  # takes over (a probe shows ~1e-16 differences at d >= 1), held to the
  # backend's rounding bound.
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    for (nt in c(32L, 33L)) {
      ts <- seq(min(x) - 0.5, max(x) + 0.5, length.out = nt)
      ts[c(2, nt - 1)] <- c(x[1], median(x))   # include ties with observations
      for (d in 0:3) {
        lab <- sprintf("%s nt=%d d=%d", nm, nt, d)
        if (nt == 32L) {
          expect_identical(LPM(d, ts, x), lpm_each(d, ts, x), label = paste(lab, "LPM"))
          expect_identical(UPM(d, ts, x), upm_each(d, ts, x), label = paste(lab, "UPM"))
        } else {
          tol <- 8 * (length(x) + d) * eps * max(abs_mom(x, d, min(ts)), abs_mom(x, d, max(ts)))
          expect_within(LPM(d, ts, x), lpm_each(d, ts, x), tol, label = paste(lab, "LPM"))
          expect_within(UPM(d, ts, x), upm_each(d, ts, x), tol, label = paste(lab, "UPM"))
        }
      }
    }
  }
})

test_that("LPM_monotone_target / UPM_antitone_target: monotone in the target, every degree", {
  expect_true(all(diff(LPM(0, t_dy, x_dy)) >= 0))
  expect_true(all(diff(UPM(0, t_dy, x_dy)) <= 0))
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- sort(c(x[1:5], seq(min(x) - 1, max(x) + 1, length.out = 100)))
    for (d in 0:3) {
      tol <- 8 * length(x) * eps * max(abs_mom(x, d, min(x) - 1), abs_mom(x, d, max(x) + 1))
      expect_gte(min(diff(LPM(d, ts, x))), -tol, label = sprintf("%s LPM d=%d", nm, d))
      expect_lte(max(diff(UPM(d, ts, x))), tol, label = sprintf("%s UPM d=%d", nm, d))
    }
  }
})

test_that("LPM_eq_zero_of_lt / UPM_eq_zero_of_le: tails vanish; degree 0 tie at max is UPM = 0", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    for (d in 0:3) {
      expect_identical(LPM(d, min(x) - 1, x), 0, label = nm)
      expect_identical(UPM(d, max(x), x), 0, label = nm)
    }
    expect_gt(LPM(0, min(x), x), 0)  # lower side counts the tie at min
  }
})

test_that("LPM_pos_of_lt: strictly above some observation the LPM is positive, every degree", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    t <- min(x) + diff(range(x)) / 100
    for (d in 0:3) expect_gt(LPM(d, t, x), 0, label = sprintf("%s d=%d", nm, d))
  }
})

test_that("LPM_one_of_forall_le: at or above every observation LPM(1, t) = t - mean", {
  expect_identical(LPM(1, c(7, 8, 100), x_int), c(7, 8, 100) - 3)
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- max(x) + c(0, 0.5, 10)
    expect_within(LPM(1, ts, x), ts - mean(x), mom_tol(x, 1, max(x)) + 4 * eps * max(abs(ts)),
                  label = nm)
  }
})

test_that("mean_sub_eq_integral: mean - t = integral_t^b (1 - F) - integral_a^t F, exact", {
  k <- knot_grid(x_int)
  F_int <- cum_step(k, function(s) LPM(0, s, x_int))
  S_int <- rev(c(0, cumsum(rev((1 - LPM(0, head(k, -1), x_int)) * diff(k)))))
  expect_identical(mean(x_int) - k, S_int - F_int)
})

# ---- dominance: lifting (FSD => SSD => TSD) --------------------------------

test_that("LPM_one_le_of_LPM_zero_le: FSD.uni == 1 implies SSD.uni == 1 (constructed FSD pairs)", {
  expect_gte(sum(fsd_verdicts[, "fsd"]), 90)
  expect_identical(fsd_verdicts[fsd_verdicts[, "fsd"] == 1L, "ssd"],
                   rep(1L, sum(fsd_verdicts[, "fsd"])))
})

test_that("LPM_two_le_of_LPM_zero_le: FSD.uni == 1 implies TSD.uni == 1 (constructed FSD pairs)", {
  expect_identical(fsd_verdicts[fsd_verdicts[, "fsd"] == 1L, "tsd"],
                   rep(1L, sum(fsd_verdicts[, "fsd"])))
})

test_that("LPM_one_le_of_LPM_zero_le / LPM_two_le_of_LPM_one_le: FSD => SSD => TSD on random integer pairs", {
  v <- int_verdicts
  expect_gt(sum(v[, "fsd"]), 20)
  expect_true(all(v[v[, "fsd"] == 1L, "ssd"] == 1L))
  expect_true(all(v[v[, "ssd"] == 1L, "tsd"] == 1L))
})

test_that("LPM_add_le_of_forall_le: FSD => SSD => TSD on real data (mtcars, faithful)", {
  pairs <- list(
    mpg = list(x = mtcars$mpg[mtcars$cyl == 4], y = mtcars$mpg[mtcars$cyl == 8]),
    eruptions = list(x = faithful$eruptions[faithful$waiting > 70],
                     y = faithful$eruptions[faithful$waiting <= 70]),
    dax_shift = list(x = real_samples$dax, y = real_samples$dax - 0.001)
  )
  for (nm in names(pairs)) {
    expect_identical(verdicts(pairs[[nm]]), c(fsd = 1L, ssd = 1L, tsd = 1L), label = nm)
  }
})

# ---- dominance: R's knot walk is exact (DominanceAlgo.lean) -----------------

test_that("sdDomPair_ssd_iff: SSD.uni == 1 iff LPM1_x <= LPM1_y at all merged knots, strict at one (base-R oracle)", {
  # LPM_one_le_of_le_at_knots gives the knot sufficiency, LPM_one_eq_of_eq_at_knots
  # makes knot-only strictness exact; the oracle is mean(pmax(t - x, 0)).
  ssd_knot <- vapply(int_pairs, function(p) {
    k <- sort(unique(c(p$x, p$y)))
    d <- lpm1_ref(k, p$x) - lpm1_ref(k, p$y)
    as.integer(all(d <= 0) && any(d < 0))
  }, integer(1))
  expect_gt(sum(ssd_knot), 50)
  expect_identical(int_verdicts[, "ssd"], ssd_knot)
})

test_that("sdDomPair_ssd_iff: SSD knot criterion on random rounded and real data (base-R oracle)", {
  set.seed(3)
  pairs <- list(
    list(x = random_samples$normal[1:200] + 0.3, y = random_samples$normal[201:400]),
    list(x = random_samples$t3[1:250], y = random_samples$t3[251:500]),
    list(x = round(rexp(150), 2), y = round(rexp(120), 2) - 0.2),
    list(x = real_samples$dax, y = as.numeric(diff(log(EuStockMarkets[, "FTSE"])))),
    list(x = as.numeric(diff(log(EuStockMarkets[, "SMI"]))), y = as.numeric(diff(log(EuStockMarkets[, "CAC"])))),
    list(x = real_samples$mag + 0.05, y = real_samples$mag),
    list(x = mtcars$mpg[mtcars$am == 1], y = mtcars$mpg[mtcars$am == 0])
  )
  for (i in seq_along(pairs)) {
    p <- pairs[[i]]
    k <- sort(unique(c(p$x, p$y)))
    d <- lpm1_ref(k, p$x) - lpm1_ref(k, p$y)
    S <- max(abs(c(p$x, p$y)))
    want <- ssd_expected(d, 64 * eps * S + 4 * (length(p$x) + length(p$y)) * eps * S)
    expect_false(is.na(want), label = paste("pair", i, "outside the tolerance band"))
    expect_identical(NNS.SSD.uni(p$x, p$y), want, label = paste("pair", i))
  }
})

test_that("sdDomPair_tsd_iff: TSD.uni equals the exact knot + vertex + mean-gate criterion on random integer pairs", {
  tsd_ref <- vapply(int_pairs, function(p) tsd_exact(p$x, p$y), integer(1))
  expect_gt(sum(tsd_ref), 50)
  # on this regime only the knots and the mean gate decide
  mean_decided <- vapply(int_pairs, function(p) {
    q <- tsd_parts(p$x, p$y); q$knots_ok && q$differs && !q$mean_ok
  }, logical(1))
  expect_gt(sum(mean_decided), 0)
  expect_identical(int_verdicts[, "tsd"], tsd_ref)
})

test_that("sdDomPair_tsd_iff (vertex branch): TSD.uni agrees with the exact criterion where the interior vertex fires", {
  # LPM_two_le_of_le_at_knots needs the vertex hypothesis (src/SD.cpp:264-277);
  # here the knots, min and mean gates all pass and the vertex alone says 0.
  expect_identical(n_vertex_decided, 25L)
  decided <- vapply(vertex_pairs, `[[`, logical(1), "decided")
  expect_gt(sum(!decided), 0)
  for (i in seq_along(vertex_pairs)) {
    p <- vertex_pairs[[i]]
    expect_identical(tsd_exact(p$x, p$y), 0L, label = paste("vertex pair", i))
    expect_identical(NNS.TSD.uni(p$x, p$y), 0L, label = paste("vertex pair", i))
  }
  # the same pairs reversed exercise the vertex walk on the other orientation
  rev_ref <- vapply(vertex_pairs, function(p) tsd_exact(p$y, p$x), integer(1))
  rev_got <- vapply(vertex_pairs, function(p) NNS.TSD.uni(p$y, p$x), integer(1))
  expect_identical(rev_got, rev_ref)
})

test_that("tsdKnotOnly_vertexCex / sdDomPair_tsd_vertexCex: knot-only TSD check is insufficient; x=(0,4,5,1,4), y=(4,2,4,0) is no-TSD", {
  x <- c(0, 4, 5, 1, 4); y <- c(4, 2, 4, 0)
  k <- sort(unique(c(x, y)))
  expect_true(all(LPM(2, k, x) - LPM(2, k, y) <= 0))     # tsdVertexCex_le_at_knots
  expect_gte(mean(x), mean(y))                            # tsdVertexCex_mean_le
  expect_equal(LPM(2, 3, x) - LPM(2, 3, y), 0.1, tolerance = 8 * eps)  # tsdVertexCex_not_le_at_three
  q <- tsd_parts(x, y)
  expect_true(q$knots_ok && q$mean_ok && q$differs)       # tsdKnotOnly_vertexCex: knot-only says 1
  expect_true(q$vertex_bad)                               # the vertex branch decides
  expect_identical(tsd_exact(x, y), 0L)
  expect_identical(NNS.TSD.uni(x, y), 0L)                 # sdDomPair_tsd_vertexCex
  expect_identical(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
})

test_that("tsdCex (Integration.lean): the TSD mean gate is load-bearing; x = 0, y = (-1, 2) is no-TSD", {
  x <- 0; y <- c(-1, 2)
  s <- seq(-2, 2, by = 0.25)
  expect_true(all(LPM(2, s, x) <= LPM(2, s, y)))  # tsdCex_le_on_Iic
  expect_lt(mean(x), mean(y))                      # tsdCex_mean_lt
  expect_lt(LPM(2, 3, y), LPM(2, 3, x))            # tsdCex_not_le_at_three
  expect_identical(tsd_exact(x, y), 0L)            # the mean gate decides
  expect_identical(NNS.TSD.uni(x, y), 0L)
  expect_identical(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
})

test_that("exists_le_of_forall_LPM_le / mean_le_of_LPM_one_le (SSD) / mean_le_of_LPM_two_le (TSD): dominance forces min and mean gates", {
  for (i in seq_along(int_pairs)) {
    p <- int_pairs[[i]]
    if (any(int_verdicts[i, ] == 1L) || int_omega[i] == 1L) expect_gte(min(p$x), min(p$y))
    if (int_verdicts[i, "ssd"] == 1L) {                       # mean_le_of_LPM_one_le
      expect_gte(sum(p$x) * length(p$y), sum(p$y) * length(p$x))
    }
    if (int_verdicts[i, "tsd"] == 1L) {                       # mean_le_of_LPM_two_le
      expect_gte(sum(p$x) * length(p$y), sum(p$y) * length(p$x))
    }
  }
})

# ---- continuous FSD (Omega / area ratio) ------------------------------------

test_that("sdDomPair_omega_iff: FSD.uni(continuous) equals the exact LPM.ratio(1) criterion on random integer pairs", {
  parts <- lapply(int_pairs, function(p) omega_parts(p$x, p$y))
  ref <- vapply(int_pairs, function(p) omega_exact(p$x, p$y), integer(1))
  expect_gt(sum(ref), 50)
  # both interior branches of the kernel are exercised: the stationary point
  # decides (knots pass, interior crossing) and the midpoint alone shows the
  # curves differ (equal at every knot)
  expect_gt(sum(vapply(parts, function(q) q$knots_ok && q$vertex_bad, logical(1))), 0)
  expect_gt(sum(vapply(parts, function(q) q$knots_ok && !q$vertex_bad && !q$knot_differs &&
                         q$mid_differs, logical(1))), 0)
  expect_identical(int_omega, ref)
  expect_identical(int_omega_rev, vapply(int_pairs, function(p) omega_exact(p$y, p$x), integer(1)))
})

test_that("sdDomPair_omega_iff: book pair x = (0,1,4,5), y = (0,2,2.5,3) is Omega but not FSD", {
  x <- c(0, 1, 4, 5); y <- c(0, 2, 2.5, 3)
  expect_identical(omega_exact(x, y), 1L)
  expect_identical(NNS.FSD.uni(x, y, type = "continuous"), 1L)
  expect_identical(NNS.FSD.uni(x, y, type = "discrete"), 0L)
  expect_identical(NNS.FSD(x, y, type = "continuous", plot = FALSE), "X OMEGA Y")
  expect_identical(NNS.FSD(x, y, type = "discrete", plot = FALSE), "NO FSD EXISTS")
})

test_that("sdDomPair_omega_iff: crossing ratio curves c(1,1) vs c(0,2) give 0 in both directions", {
  expect_identical(omega_exact(c(1, 1), c(0, 2)), 0L)
  expect_identical(omega_exact(c(0, 2), c(1, 1)), 0L)
  expect_identical(NNS.FSD.uni(c(1, 1), c(0, 2), type = "continuous"), 0L)
  expect_identical(NNS.FSD.uni(c(0, 2), c(1, 1), type = "continuous"), 0L)
  expect_identical(NNS.FSD(c(1, 1), c(0, 2), type = "continuous", plot = FALSE),
                   "NO OMEGA DOMINANCE EXISTS")
})

test_that("sdDomPair_omega_iff (0/0 -> k/m fallback): constant samples", {
  # all observations at t: L1 + U1 = 0 and the ratio is k/m (LPM.ratio gives 1)
  expect_identical(LPM.ratio(1, 2, c(2, 2)), 1)
  expect_identical(NNS.FSD.uni(c(2, 2), 2, type = "continuous"), 0L)  # same distribution
  expect_identical(NNS.FSD.uni(2, c(2, 2), type = "continuous"), 0L)
  expect_identical(NNS.FSD.uni(c(3, 3), c(2, 2), type = "continuous"), 1L)
  expect_identical(NNS.FSD.uni(c(2, 2), c(3, 3), type = "continuous"), 0L)
  expect_identical(omega_exact(c(3, 3), c(2, 2)), 1L)
})

test_that("LPM_one_le_of_LPM_zero_le + UPM_succ_le_of_forall_le (d = 0) => sdDomPair_omega_iff: discrete FSD == 1 implies Omega == 1", {
  # F_x <= F_y gives L1_x <= L1_y (LPM_one_le_of_LPM_zero_le) and U1_x >= U1_y
  # (UPM_succ_le_of_forall_le at d = 0, since UPM0 = 1 - F), hence R_x <= R_y
  expect_gt(sum(int_verdicts[, "fsd"]), 20)
  expect_true(all(int_omega[int_verdicts[, "fsd"] == 1L] == 1L))
  fsd_omega <- vapply(fsd_pairs, function(p) NNS.FSD.uni(p$x, p$y, type = "continuous"), integer(1))
  expect_true(all(fsd_omega[fsd_verdicts[, "fsd"] == 1L] == 1L))
})

# ---- verdict shape ----------------------------------------------------------

test_that("sdDomPair_eq_zero_or_one: every .uni verdict is the integer 0 or 1", {
  for (v in list(int_verdicts, int_verdicts_rev, fsd_verdicts, int_omega, int_omega_rev)) {
    expect_type(v, "integer")
    expect_true(all(v %in% c(0L, 1L)))
  }
})

test_that("sdDomPair_{fsd,ssd,tsd,omega}_iff: antisymmetry, never both X dominates Y and Y dominates X", {
  # exact dominance is antisymmetric (curves <= both ways means equal curves)
  expect_false(any(int_verdicts == 1L & int_verdicts_rev == 1L))
  expect_false(any(int_omega == 1L & int_omega_rev == 1L))
})

test_that("LPM_eq_of_identical / sdDomPair_*_iff: a sample never dominates itself or a permutation of itself", {
  for (x in c(list(x_int, x_dy, 3), random_samples[1], real_samples["mpg"])) {
    expect_identical(verdicts(list(x = x, y = rev(x))), c(fsd = 0L, ssd = 0L, tsd = 0L))
    expect_identical(NNS.FSD.uni(x, rev(x), type = "continuous"), 0L)
  }
})

test_that("sdDomPair_*_iff (curves never differ): equal distributions of different sizes give 0", {
  x <- c(1, 2); y <- c(1, 1, 2, 2)
  expect_identical(verdicts(list(x = x, y = y)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(verdicts(list(x = y, y = x)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(NNS.FSD.uni(x, y, type = "continuous"), 0L)
  expect_identical(NNS.FSD.uni(y, x, type = "continuous"), 0L)
})

test_that("sdDomPair_*_iff: n = 1 samples, 2 dominates 1 in every order, 1 vs 1 in none", {
  expect_identical(verdicts(list(x = 2, y = 1)), c(fsd = 1L, ssd = 1L, tsd = 1L))
  expect_identical(verdicts(list(x = 1, y = 2)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(verdicts(list(x = 1, y = 1)), c(fsd = 0L, ssd = 0L, tsd = 0L))
})

test_that("sdDomPair_*_iff: constant samples, equal constants give 0, a higher constant dominates", {
  expect_identical(verdicts(list(x = c(2, 2), y = 2)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(verdicts(list(x = rep(3, 4), y = c(1, 3))), c(fsd = 1L, ssd = 1L, tsd = 1L))
  # mean-preserving spread: the constant SSD/TSD-dominates, FSD fails
  expect_identical(verdicts(list(x = c(1, 1), y = c(0, 2))), c(fsd = 0L, ssd = 1L, tsd = 1L))
})

test_that("LPMDominates_shift / sdDomPairCentred_eq / sdDomPairCentred_omega_eq: shifting both samples by 1024 keeps every verdict", {
  shifted <- t(vapply(int_pairs, function(p) verdicts(list(x = p$x + 1024, y = p$y + 1024)),
                      integer(3)))
  expect_identical(shifted, int_verdicts)
  shifted_omega <- vapply(int_pairs, function(p) NNS.FSD.uni(p$x + 1024, p$y + 1024, type = "continuous"),
                          integer(1))
  expect_identical(shifted_omega, int_omega)
})

test_that("bidirectional wrappers agree with the .uni verdicts and swap labels with the arguments (SPEC bidirectional wrappers)", {
  for (i in seq_len(60)) {
    p <- int_pairs[[i]]
    for (ord in c("FSD", "SSD", "TSD")) {
      fn <- switch(ord, FSD = NNS.FSD, SSD = NNS.SSD, TSD = NNS.TSD)
      col <- tolower(ord)
      want <- if (int_verdicts[i, col] == 1L) paste("X", ord, "Y")
              else if (int_verdicts_rev[i, col] == 1L) paste("Y", ord, "X")
              else paste("NO", ord, "EXISTS")
      expect_identical(fn(p$x, p$y, plot = FALSE), want)
    }
    want <- if (int_omega[i] == 1L) "X OMEGA Y"
            else if (int_omega_rev[i] == 1L) "Y OMEGA X"
            else "NO OMEGA DOMINANCE EXISTS"
    expect_identical(NNS.FSD(p$x, p$y, type = "continuous", plot = FALSE), want)
  }
})
