# Property tests for NNS.moments and the stochastic-dominance verdicts.
#
# Each test_that names the Lean theorem it mirrors in
# nns-proofs/NNSProofs/Moments.lean (NNS.moments as partial moments about the
# mean) or nns-proofs/NNSProofs/Integration.lean (LPM_{d+1} = (d+1) * integral
# of LPM_d, dominance lifting FSD => SSD => TSD, exactness of R's knot walk).
# Verdict-shape tests without a Lean theorem cite SPEC.md instead.
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

# Exact TSD criterion of LPM_two_le_of_le_at_knots for integer samples:
# knots, interior vertices s* = b / a, mean gate, strictness at a knot.
# All quantities are integers scaled by mx * my, so the test is exact.
tsd_exact <- function(x, y) {
  mx <- length(x); my <- length(y)
  k <- sort(unique(c(x, y)))
  A <- vapply(k, function(t) sum(x <= t) * my - sum(y <= t) * mx, numeric(1))
  B <- vapply(k, function(t) sum(x[x <= t]) * my - sum(y[y <= t]) * mx, numeric(1))
  C <- vapply(k, function(t) sum(x[x <= t]^2) * my - sum(y[y <= t]^2) * mx, numeric(1))
  q_knot <- A * k^2 - 2 * B * k + C
  lo <- head(k, -1); hi <- k[-1]; a <- head(A, -1); b <- head(B, -1); cc <- head(C, -1)
  inside <- a < 0 & lo * a > b & b > hi * a
  vertex_bad <- any(inside & (cc * a - b^2 < 0))
  mean_ok <- sum(x) * my >= sum(y) * mx
  as.integer(all(q_knot <= 0) && !vertex_bad && mean_ok && any(q_knot != 0))
}

lpm_each <- function(d, ts, x) vapply(ts, function(t) LPM(d, t, x), numeric(1))

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

test_that("LPM_succ_eq_integral: prefix backend (>32 targets) agrees with scalar targets", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- sort(c(targets_of(x), seq(min(x), max(x), length.out = 40)))
    for (d in 0:3) {
      expect_within(LPM(d, ts, x), lpm_each(d, ts, x),
                    8 * (length(x) + d) * eps * max(abs_mom(x, d, min(x)), abs_mom(x, d, max(x))),
                    label = sprintf("%s d=%d", nm, d))
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

# ---- dominance: R's knot walk is exact --------------------------------------

test_that("LPM_one_le_of_le_at_knots: SSD.uni == 1 iff LPM1_x <= LPM1_y at all merged knots (strict somewhere)", {
  ssd_knot <- vapply(int_pairs, function(p) {
    k <- sort(unique(c(p$x, p$y)))
    d <- LPM(1, k, p$x) - LPM(1, k, p$y)
    as.integer(all(d <= 0) && any(d < 0))
  }, integer(1))
  expect_identical(int_verdicts[, "ssd"], ssd_knot)
})

test_that("LPM_one_le_of_le_at_knots: SSD knot criterion on random rounded and real data", {
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
    d <- LPM(1, k, p$x) - LPM(1, k, p$y)
    S <- max(abs(c(p$x, p$y)))
    want <- ssd_expected(d, 64 * eps * S + 4 * (length(p$x) + length(p$y)) * eps * S)
    expect_false(is.na(want), label = paste("pair", i, "outside the tolerance band"))
    expect_identical(NNS.SSD.uni(p$x, p$y), want, label = paste("pair", i))
  }
})

test_that("LPM_two_le_of_le_at_knots: TSD.uni equals the exact knot + vertex + mean-gate criterion", {
  tsd_ref <- vapply(int_pairs, function(p) tsd_exact(p$x, p$y), integer(1))
  expect_gt(sum(tsd_ref), 50)
  expect_identical(int_verdicts[, "tsd"], tsd_ref)
})

test_that("tsdVertexCex: knot-only TSD check is insufficient; x=(0,4,5,1,4), y=(4,2,4,0) is no-TSD", {
  x <- c(0, 4, 5, 1, 4); y <- c(4, 2, 4, 0)
  k <- sort(unique(c(x, y)))
  expect_true(all(LPM(2, k, x) - LPM(2, k, y) <= 0))     # tsdVertexCex_le_at_knots
  expect_gte(mean(x), mean(y))                            # tsdVertexCex_mean_le
  expect_equal(LPM(2, 3, x) - LPM(2, 3, y), 0.1, tolerance = 8 * eps)  # tsdVertexCex_not_le_at_three
  expect_identical(tsd_exact(x, y), 0L)                  # the vertex branch decides
  expect_identical(NNS.TSD.uni(x, y), 0L)
  expect_identical(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
})

test_that("tsdCex: the TSD mean gate is load-bearing; x = 0, y = (-1, 2) is no-TSD", {
  x <- 0; y <- c(-1, 2)
  s <- seq(-2, 2, by = 0.25)
  expect_true(all(LPM(2, s, x) <= LPM(2, s, y)))  # tsdCex_le_on_Iic
  expect_lt(mean(x), mean(y))                      # tsdCex_mean_lt
  expect_lt(LPM(2, 3, y), LPM(2, 3, x))            # tsdCex_not_le_at_three
  expect_identical(tsd_exact(x, y), 0L)            # the mean gate decides
  expect_identical(NNS.TSD.uni(x, y), 0L)
  expect_identical(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
})

test_that("exists_le_of_forall_LPM_le / mean_le_of_LPM_one_le: dominance forces min and mean gates", {
  for (i in seq_along(int_pairs)) {
    p <- int_pairs[[i]]
    if (any(int_verdicts[i, ] == 1L)) expect_gte(min(p$x), min(p$y))
    if (int_verdicts[i, "ssd"] == 1L || int_verdicts[i, "tsd"] == 1L) {
      expect_gte(sum(p$x) * length(p$y), sum(p$y) * length(p$x))
    }
  }
})

# ---- verdict shape (SPEC "sd_dom_pair (all orders): verdict semantics") -----

test_that("verdict antisymmetry: never both X dominates Y and Y dominates X (SPEC verdict semantics)", {
  expect_false(any(int_verdicts == 1L & int_verdicts_rev == 1L))
})

test_that("verdict irreflexivity: a sample never dominates itself or a permutation of itself", {
  for (x in c(list(x_int, x_dy, 3), random_samples[1], real_samples["mpg"])) {
    expect_identical(verdicts(list(x = x, y = rev(x))), c(fsd = 0L, ssd = 0L, tsd = 0L))
  }
})

test_that("verdict on equal distributions of different sizes is 0 (curves identical)", {
  x <- c(1, 2); y <- c(1, 1, 2, 2)
  expect_identical(verdicts(list(x = x, y = y)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(verdicts(list(x = y, y = x)), c(fsd = 0L, ssd = 0L, tsd = 0L))
})

test_that("n = 1 samples: 2 dominates 1 in every order, 1 vs 1 in none", {
  expect_identical(verdicts(list(x = 2, y = 1)), c(fsd = 1L, ssd = 1L, tsd = 1L))
  expect_identical(verdicts(list(x = 1, y = 2)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(verdicts(list(x = 1, y = 1)), c(fsd = 0L, ssd = 0L, tsd = 0L))
})

test_that("constant samples: equal constants of different sizes give 0, a higher constant dominates", {
  expect_identical(verdicts(list(x = c(2, 2), y = 2)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  expect_identical(verdicts(list(x = rep(3, 4), y = c(1, 3))), c(fsd = 1L, ssd = 1L, tsd = 1L))
  # mean-preserving spread: the constant SSD/TSD-dominates, FSD fails
  expect_identical(verdicts(list(x = c(1, 1), y = c(0, 2))), c(fsd = 0L, ssd = 1L, tsd = 1L))
})

test_that("verdict translation invariance: shifting both samples by 1024 keeps every verdict", {
  shifted <- t(vapply(int_pairs, function(p) verdicts(list(x = p$x + 1024, y = p$y + 1024)),
                      integer(3)))
  expect_identical(shifted, int_verdicts)
})

test_that("bidirectional wrappers agree with the .uni verdicts and swap labels with the arguments", {
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
  }
})
