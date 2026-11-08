# Property tests mirroring the Lean theorems of nns-proofs:
#   NNSProofs/Ratios.lean  (LPM.ratio / UPM.ratio)
#   NNSProofs/CDF.lean     (univariate NNS.CDF: CDF, survival, cumulative hazard)
#   NNSProofs/VaRRoot.lean (LPM.VaR / UPM.VaR: d > 0 roots of the ratio via
#                           existsUnique_LPMratio_eq / LPMratio_strictMonoOn_Icc;
#                           d = 0 is the type-7 quantile, quantile7_eq).
# Each test_that names the theorem it mirrors; tests with no Lean counterpart say so
# ("BUGS #N regression", "test-derived bound"). Three data regimes:
#   exact  : small dyadic data, n a power of two, arithmetic exact in double;
#   random : fixed seed, normal / exponential / t(3), ties forced by rounding;
#   real   : datasets shipped with R, targets at mean, median, quantiles, data points.

eps <- .Machine$double.eps

# Absolute tolerance for a ratio in [0, 1] computed from n terms.
ratio_tol <- function(x) 8 * length(x) * eps

expect_close <- function(object, expected, tol, info = NULL) {
  expect_equal(length(object), length(expected), info = info)
  # NaN-safe: a NaN difference is a failure, not an error.
  dev <- max(abs(object - expected))
  expect_true(is.finite(dev) && dev <= tol,
              info = paste("max |diff| =", format(dev, digits = 3), "> tol =", format(tol, digits = 3), info))
}

# Exact regime: n = 8, dyadic values, ties at 0 and 1, dyadic mean 0.34375.
x_exact <- c(-2, -0.5, 0, 0, 0.25, 1, 1, 3)

set.seed(20261108)
samples <- list(
  exact   = x_exact,
  norm    = round(rnorm(200), 2),
  exp     = round(rexp(150), 1),
  t3      = round(rt(300, df = 3), 1),
  dax     = as.numeric(diff(log(EuStockMarkets))[, "DAX"]),
  waiting = faithful$waiting,
  ozone   = as.numeric(na.omit(airquality$Ozone)),
  mag     = quakes$mag,
  mpg     = mtcars$mpg
)

targets_for <- function(x) {
  unique(c(mean(x), median(x), quantile(x, c(.05, .25, .75, .95), names = FALSE),
           x[1], x[length(x) %/% 2], min(x), max(x)))
}

ecdf_count <- function(t, x) vapply(t, function(s) sum(x <= s), numeric(1)) / length(x)
surv_count <- function(t, x) vapply(t, function(s) sum(x > s), numeric(1)) / length(x)

# More than 32 targets: LPM_CPv switches to the prefix backend. The margin is relative
# to the range so that most grid points are strictly inside (min x, max x) for every
# sample (dax has range ~0.15: an absolute margin of 1 left it 4 interior points).
grid_for <- function(x, k = 64) {
  m <- 0.1 * diff(range(x))
  sort(c(seq(min(x) - m, max(x) + m, length.out = k), x[1:min(8, length(x))]))
}

# ---------------------------------------------------------------------------
# Ratios.lean: range and complementarity
# ---------------------------------------------------------------------------

test_that("LPMratio_mem_Icc / UPMratio_mem_Icc: both ratios lie in [0, 1]", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    t <- c(targets_for(x), grid_for(x))
    for (r in list(LPM.ratio(d, t, x), UPM.ratio(d, t, x))) {
      expect_true(all(r >= 0 & r <= 1), info = paste(nm, d))
    }
  }
})

test_that("BUGS #12 regression (no Lean theorem; under/overflow is outside the exact-real model of Ratios.lean): high degree on symmetric data gives 1/2, not NaN", {
  # Symmetric two-point sample: normalised deviations are +/-1, so 0.5 is exact.
  # fails on 8898017: see BUGS.md #12 (L and U underflow / overflow to 0/0 or Inf/Inf -> NaN)
  expect_identical(LPM.ratio(400, 0, c(-1e-2, 1e-2)), 0.5)
  expect_identical(LPM.ratio(400, 0, c(-1e3, 1e3)), 0.5)
})

test_that("LPMratio_add_UPMratio: degree 0 sums to 1 exactly on dyadic data", {
  t <- c(-3, x_exact, 0.1, 2, 4)
  expect_identical(LPM.ratio(0, t, x_exact) + UPM.ratio(0, t, x_exact), rep(1, length(t)))
})

test_that("UPMratio_eq_one_sub_LPMratio: UPM.ratio = 1 - LPM.ratio at every degree", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    t <- targets_for(x)
    expect_close(UPM.ratio(d, t, x), 1 - LPM.ratio(d, t, x), 4 * eps, info = paste(nm, d))
  }
})

test_that("UPMratio_eq_one_sub_LPMratio: also on the > 32 target prefix path", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    t <- grid_for(x)
    expect_close(UPM.ratio(d, t, x), 1 - LPM.ratio(d, t, x), 4 * eps, info = paste(nm, d))
  }
})

test_that("LPMratio_eq_sum_div: d > 0 ratio is sum (t-x)_+^d / sum |x-t|^d", {
  for (nm in names(samples)) for (d in 1:3) {
    x <- samples[[nm]]
    t <- targets_for(x)
    t <- t[t > min(x)]  # a target off the sample (or with mass on both sides)
    manual <- vapply(t, function(s) sum(pmax(s - x, 0)^d) / sum(abs(x - s)^d), numeric(1))
    expect_close(LPM.ratio(d, t, x), manual, ratio_tol(x), info = paste(nm, d))
  }
})

# ---------------------------------------------------------------------------
# Ratios.lean: degree 0 is the empirical CDF / survival
# ---------------------------------------------------------------------------

test_that("LPMratio_zero_eq_ecdf: degree-0 LPM.ratio is #{x <= t}/n (ties low)", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    t <- c(targets_for(x), x)  # exact data points: ties at the target
    expect_identical(LPM.ratio(0, t, x), ecdf_count(t, x), info = nm)
  }
})

test_that("LPMratio_zero_eq_ecdf: prefix path (> 32 targets) is the ECDF too", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    t <- grid_for(x)
    expect_identical(LPM.ratio(0, t, x), ecdf_count(t, x), info = nm)
  }
})

test_that("UPMratio_zero_eq_survival: degree-0 UPM.ratio is #{x > t}/n (strict)", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    t <- c(targets_for(x), x)
    expect_identical(UPM.ratio(0, t, x), surv_count(t, x), info = nm)
  }
})

test_that("LPMratio (elementwise target): > 32 targets agree with scalar calls", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    t <- grid_for(x)
    scalar <- vapply(t, function(s) LPM.ratio(d, s, x), numeric(1))
    expect_close(LPM.ratio(d, t, x), scalar, ratio_tol(x), info = paste(nm, d))
  }
})

# ---------------------------------------------------------------------------
# Ratios.lean: degenerate 0/0 convention, n = 1, NULL target
# ---------------------------------------------------------------------------

test_that("LPMratio_of_forall_eq / UPMratio_of_forall_eq: constant sample at its value gives 1 / 0", {
  for (d in 1:4) for (x in list(rep(5, 3), rep(-0.25, 8), 7)) {
    expect_identical(LPM.ratio(d, x[1], x), 1, info = paste(d, length(x)))
    expect_identical(UPM.ratio(d, x[1], x), 0, info = paste(d, length(x)))
  }
})

test_that("LPMratio_of_forall_eq: constant sample is the step 1[c <= t] at every degree", {
  for (d in 0:3) {
    expect_identical(LPM.ratio(d, c(4, 5, 6), rep(5, 3)), c(0, 1, 1), info = d)
    expect_identical(UPM.ratio(d, c(4, 5, 6), rep(5, 3)), c(1, 0, 0), info = d)
  }
})

test_that("LPMratio_of_forall_eq: n = 1 is constant, step at the single value", {
  for (d in 0:3) {
    expect_identical(LPM.ratio(d, c(0, 1, 2), 1), c(0, 1, 1), info = d)
    expect_identical(UPM.ratio(d, c(0, 1, 2), 1), c(1, 0, 0), info = d)
  }
})

test_that("LPMratioDefault_of_const: NULL target on a constant sample gives 1, not 1/2", {
  for (d in 0:3) {
    expect_identical(LPM.ratio(d, NULL, rep(2.5, 4)), 1, info = d)
    expect_identical(UPM.ratio(d, NULL, rep(2.5, 4)), 0, info = d)
  }
})

test_that("LPMratioDefault: NULL target is the sample mean", {
  for (nm in names(samples)) for (d in 0:3) {
    x <- samples[[nm]]
    expect_identical(LPM.ratio(d, NULL, x), LPM.ratio(d, mean(x), x), info = paste(nm, d))
  }
})

# ---------------------------------------------------------------------------
# Ratios.lean: boundary values
# ---------------------------------------------------------------------------

test_that("LPMratio_eq_one_of_forall_le: at or above the maximum LPM.ratio = 1 exactly", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    expect_identical(LPM.ratio(d, c(max(x), max(x) + 1), x), c(1, 1), info = paste(nm, d))
  }
})

test_that("LPMratio_eq_zero_of_forall_lt: strictly below the minimum LPM.ratio = 0 exactly", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    expect_identical(LPM.ratio(d, min(x) - 1, x), 0, info = paste(nm, d))
  }
})

test_that("LPMratio_eq_zero_of_forall_ge: d > 0 at the minimum gives 0, degree 0 gives the atom", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    for (d in 1:4) expect_identical(LPM.ratio(d, min(x), x), 0, info = paste(nm, d))
    expect_identical(LPM.ratio(0, min(x), x), sum(x == min(x)) / length(x), info = nm)
  }
})

test_that("UPMratio_eq_one_of_forall_ge: d > 0 at the minimum UPM.ratio = 1, degree 0 drops the atom", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    for (d in 1:4) expect_identical(UPM.ratio(d, min(x), x), 1, info = paste(nm, d))
    expect_identical(UPM.ratio(0, min(x), x), sum(x > min(x)) / length(x), info = nm)
  }
})

test_that("LPMratio_mem_Ioo: d > 0 and a target strictly inside the range give 0 < ratio < 1", {
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    t <- targets_for(x)
    t <- t[t > min(x) & t < max(x)]
    r <- LPM.ratio(d, t, x)
    expect_true(all(r > 0 & r < 1), info = paste(nm, d))
  }
})

# ---------------------------------------------------------------------------
# Ratios.lean: monotonicity in the target
# ---------------------------------------------------------------------------

test_that("LPMratio_monotone: exact dyadic data, nondecreasing with zero tolerance", {
  t <- seq(-3, 4, by = 0.125)
  for (d in 0:3) expect_true(all(diff(LPM.ratio(d, t, x_exact)) >= 0), info = d)
})

test_that("LPMratio_monotone / UPMratio_antitone: random and real data, both target paths", {
  for (nm in names(samples)) for (d in 0:4) {
    x <- samples[[nm]]
    t <- grid_for(x)
    tol <- ratio_tol(x)
    expect_gte(min(diff(LPM.ratio(d, t, x))), -tol, label = paste("LPM", nm, d))
    expect_lte(max(diff(UPM.ratio(d, t, x))), tol, label = paste("UPM", nm, d))
    scalar <- vapply(t, function(s) LPM.ratio(d, s, x), numeric(1))
    expect_gte(min(diff(scalar)), -tol, label = paste("scalar", nm, d))
  }
})

test_that("LPMratio_strictMonoOn: d > 0 strictly increasing on [min x, max x]", {
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    t <- sort(unique(c(x, seq(min(x), max(x), length.out = 40))))
    t <- t[c(TRUE, diff(t) > 1e-9 * diff(range(x)))]  # drop grid points within rounding of a datum
    expect_true(all(diff(LPM.ratio(d, t, x)) > 0), info = paste(nm, d))
  }
})

test_that("LPMratio_strictMonoOn fails at degree 0: the ECDF is flat between observations", {
  xs <- sort(unique(x_exact))
  mid <- (xs[-1] + xs[-length(xs)]) / 2
  expect_identical(LPM.ratio(0, mid, x_exact), LPM.ratio(0, xs[-length(xs)], x_exact))
})

# ---------------------------------------------------------------------------
# Ratios.lean: degree 1 and the mean
# ---------------------------------------------------------------------------

test_that("LPMratio_one_mean: exact dyadic data, LPM.ratio(1, mean) is exactly 1/2", {
  expect_identical(mean(x_exact), 0.34375)
  expect_identical(LPM.ratio(1, mean(x_exact), x_exact), 0.5)
  expect_identical(UPM.ratio(1, mean(x_exact), x_exact), 0.5)
})

test_that("LPMratio_one_mean / UPMratio_one_mean: 1/2 at the mean on random and real data", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    expect_close(LPM.ratio(1, mean(x), x), 0.5, ratio_tol(x), info = nm)
    expect_close(UPM.ratio(1, mean(x), x), 0.5, ratio_tol(x), info = nm)
  }
})

test_that("LPMratioDefault_one: NULL target at degree 1 gives 1/2 on non-constant data", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    expect_close(LPM.ratio(1, NULL, x), 0.5, ratio_tol(x), info = nm)
  }
})

test_that("LPMratio_one_eq_half_iff: away from the mean the degree-1 ratio is not 1/2 (sign of t - mean)", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    t <- targets_for(x)
    t <- t[abs(t - mean(x)) > 1e-6 * diff(range(x)) & t > min(x) & t < max(x)]
    expect_identical(sign(LPM.ratio(1, t, x) - 0.5), sign(t - mean(x)), info = nm)
  }
})

# ---------------------------------------------------------------------------
# Ratios.lean: invariances and reflection
# ---------------------------------------------------------------------------

test_that("LPMratio_affine: exact dyadic data, power-of-two scale and dyadic shift", {
  t <- c(-1, 0, 0.34375, 1, 2)
  for (d in 0:3) {
    expect_identical(LPM.ratio(d, 4 * t + 0.5, 4 * x_exact + 0.5), LPM.ratio(d, t, x_exact), info = d)
    expect_identical(UPM.ratio(d, 4 * t + 0.5, 4 * x_exact + 0.5), UPM.ratio(d, t, x_exact), info = d)
  }
})

test_that("LPMratio_affine: random and real data under a = 3.7, b = -1.3", {
  for (nm in names(samples)) for (d in 1:3) {
    x <- samples[[nm]]
    t <- targets_for(x)
    t <- t[t != min(x)]  # keep the target off the affine-rounded minimum
    expect_close(LPM.ratio(d, 3.7 * t - 1.3, 3.7 * x - 1.3), LPM.ratio(d, t, x),
                 1e-12, info = paste(nm, d))
  }
})

test_that("BUGS #12 regression (no Lean theorem; LPMratio_affine in exact reals): degree 2 is scale-free under a tiny power-of-two scale", {
  # fails on 8898017: see BUGS.md #12 (deviations^2 underflow at scale 2^-600 -> NaN)
  # A fix that rescales (e.g. by max |x - t|) rounds differently at the non-dyadic
  # targets (quantiles -1.475, 2.3): ~1e-16 apart, so compare to 4 eps, not bitwise.
  a <- 2^-600
  t <- targets_for(x_exact)
  expect_close(LPM.ratio(2, a * t, a * x_exact), LPM.ratio(2, t, x_exact), 4 * eps)
})

test_that("LPMratio_eq_UPMratio_neg: d > 0 reflection on exact dyadic data", {
  t <- c(-1, 0, 0.34375, 1, 2)  # 0 and 1 are tied observations
  for (d in 1:3) {
    expect_identical(LPM.ratio(d, t, x_exact), UPM.ratio(d, -t, -x_exact), info = d)
  }
})

test_that("LPMratio_eq_UPMratio_neg: d > 0 reflection on random and real data, ties at target", {
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    t <- targets_for(x)
    expect_close(LPM.ratio(d, t, x), UPM.ratio(d, -t, -x), 4 * eps, info = paste(nm, d))
  }
})

test_that("LPMratio_ne_UPMratio_neg_of_forall_eq: reflection breaks in the degenerate case", {
  for (d in 1:3) {
    expect_identical(LPM.ratio(d, 1, c(1, 1, 1)), 1, info = d)
    expect_identical(UPM.ratio(d, -1, -c(1, 1, 1)), 0, info = d)
  }
})

test_that("LPMratio_zero_eq_UPMratio_neg_add: degree-0 reflection picks up the atom (exact)", {
  t <- c(-3, x_exact, 0.1, 4)
  atom <- vapply(t, function(s) sum(x_exact == s), numeric(1)) / length(x_exact)
  expect_identical(LPM.ratio(0, t, x_exact), UPM.ratio(0, -t, -x_exact) + atom)
})

test_that("LPMratio_zero_eq_UPMratio_neg_add: degree-0 reflection with atom, random and real data", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    t <- c(targets_for(x), x[1:5])
    atom <- vapply(t, function(s) sum(x == s), numeric(1)) / length(x)
    expect_close(LPM.ratio(0, t, x), UPM.ratio(0, -t, -x) + atom, 2 * eps, info = nm)
  }
})

# ---------------------------------------------------------------------------
# CDF.lean: degree 0 is the right-continuous ECDF
# ---------------------------------------------------------------------------

cdf_fun <- function(x, d = 0, type = "CDF", target = NULL) {
  NNS.CDF(x, d, target = target, type = type, plot = FALSE)
}

test_that("cdfCurve_probe_zero: NNS.CDF(c(3,1,2,2,10)) is 0.2 0.6 0.6 0.8 1", {
  f <- cdf_fun(c(3, 1, 2, 2, 10))$Function
  expect_identical(f$x, c(1, 2, 2, 3, 10))
  expect_identical(f$CDF, c(1, 3, 3, 4, 5) / 5)
})

test_that("cdfCurve_zero_eq_ecdf: degree-0 NNS.CDF equals stats::ecdf at the sorted sample", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    f <- cdf_fun(x)$Function
    expect_identical(f$x, sort(x), info = nm)
    expect_identical(f$CDF, stats::ecdf(x)(sort(x)), info = nm)
  }
})

test_that("cdfTable_zero_ge: degree-0 row k has CDF >= k/n (ties push it up)", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    n <- length(x)
    f <- cdf_fun(x)$Function
    expect_true(all(f$CDF >= seq_len(n) / n), info = nm)
  }
})

test_that("cdfTable_zero_of_injective: tie-free sample, degree-0 row k has CDF = k/n exactly", {
  tie_free <- list(
    exact = c(3, -2, 0.25, -0.5, 1, 0),             # dyadic, n = 6 (k/6 not dyadic)
    dyad8 = c(-2, -0.5, 0, 0.125, 0.25, 1, 2, 3),    # n = 8: k/8 exact
    dax   = unique(samples$dax),                     # large n, real data
    norm  = unique(samples$norm)
  )
  for (nm in names(tie_free)) {
    x <- tie_free[[nm]]
    n <- length(x)
    expect_false(as.logical(anyDuplicated(x)), info = nm)  # the hypothesis really holds
    expect_gt(n, 1)
    f <- cdf_fun(x)$Function
    expect_identical(f$x, sort(x), info = nm)
    expect_identical(f$CDF, seq_len(n) / n, info = nm)
  }
})

test_that("cdfCurve_zero_eventually_const_right: the degree-0 CDF is right-continuous", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    xs <- sort(unique(x))
    delta <- min(diff(xs)) / 4
    expect_identical(LPM.ratio(0, xs + delta, x), LPM.ratio(0, xs, x), info = nm)
  }
})

test_that("cdfCurve_zero_tendsto_left / cdfCurve_zero_eq_lt_add_eq: left limit is #{x < t}/n, jump is the atom", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    xs <- sort(unique(x))
    delta <- min(diff(xs)) / 4
    left <- vapply(xs, function(s) sum(x < s), numeric(1)) / length(x)
    atom <- vapply(xs, function(s) sum(x == s), numeric(1)) / length(x)
    expect_identical(LPM.ratio(0, xs - delta, x), left, info = nm)
    expect_close(LPM.ratio(0, xs, x) - LPM.ratio(0, xs - delta, x), atom, 2 * eps, info = nm)
  }
})

# ---------------------------------------------------------------------------
# CDF.lean: degree > 0 table shape
# ---------------------------------------------------------------------------

test_that("cdfCurve_probe_one / cdfCurve_probe_tied_min / cdfTarget_probe", {
  expect_close(LPM.ratio(1, 2.5, c(1, 2, 2, 3, 10)), 5 / 21, 2 * eps)
  expect_close(cdf_fun(c(1, 1, 2, 5), 1)$Function$CDF, c(0, 0, 2 / 5, 1), 2 * eps)
  expect_close(cdf_fun(c(1, 1, 2, 5), 1, target = 1.5)$target.value, 1 / 5, 2 * eps)
})

test_that("cdfTable_snd_monotone / cdfTable_snd_last: CDF column nondecreasing, last row 1", {
  for (nm in names(samples)) for (d in 0:3) {
    f <- cdf_fun(samples[[nm]], d)$Function
    expect_gte(min(diff(f$CDF)), -ratio_tol(f$x), label = paste(nm, d))
    expect_identical(f$CDF[nrow(f)], 1, info = paste(nm, d))
  }
})

test_that("cdfTable_snd_first / cdfTable_snd_eq_zero_of_tied_first: d > 0 rows at the minimum are 0", {
  for (nm in names(samples)) for (d in 1:3) {
    x <- samples[[nm]]
    f <- cdf_fun(x, d)$Function
    expect_identical(f$CDF[f$x == min(x)], rep(0, sum(x == min(x))), info = paste(nm, d))
  }
})

test_that("cdfCurve_lt_cdfCurve: d > 0 strictly increasing across distinct rows, equal on tied rows", {
  for (nm in names(samples)) for (d in 1:3) {
    f <- cdf_fun(samples[[nm]], d)$Function
    step <- diff(f$CDF)
    tied <- diff(f$x) == 0
    expect_true(all(step[!tied] > 0), info = paste(nm, d))
    expect_true(all(step[tied] == 0), info = paste(nm, d))
  }
})

# CDF.lean:369 cdfCurve_continuous (and VaRRoot.lean:208 LPMratio_continuous) prove only
# `Continuous` for d > 0 on a non-constant sample. The two Lipschitz constants below are
# derived here and are NOT formalised. With L, U the lower / upper partial moments at t,
# F0 = #{x <= t}/n, F = L / (L + U), F' = (L' U - L U') / (L + U)^2:
#   d = 1: L' = F0, U' = -(1 - F0), numerator F0 U + (1 - F0) L <= L + U, so
#          F' <= 1 / (L + U) = 1 / mean|x - t| <= 1 / mean|x - median| (median minimises).
#   d = 2: L' = 2 l1, U' = -2 u1 with l1 = mean (t - x)_+, u1 = mean (x - t)_+;
#          Cauchy-Schwarz l1 <= sqrt(F0 L), u1 <= sqrt((1 - F0) U) gives
#          l1 U + u1 L <= sqrt(L U) sqrt(L + U) <= (L + U)^(3/2) / 2, so
#          F' <= 1 / sqrt(L + U) = 1 / sqrt(mean (x - t)^2) <= 1 / sd_pop.
#          The test uses 2 / sd_pop (2x slack).
test_that("cdfCurve_continuous (quantitative, test-derived bound): d = 1 CDF is Lipschitz with constant 1 / mean|x - median|", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    xs <- sort(unique(x))
    delta <- 1e-6 * diff(range(x))
    lip <- 1 / mean(abs(x - median(x)))
    jump <- LPM.ratio(1, xs + delta, x) - LPM.ratio(1, xs - delta, x)
    expect_lte(max(jump), 2 * delta * lip * (1 + 1e-6) + ratio_tol(x), label = nm)
  }
})

test_that("cdfCurve_continuous (quantitative, test-derived bound): d = 2 CDF is Lipschitz with constant 2 / sd_pop", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    xs <- sort(unique(x))
    delta <- 1e-6 * diff(range(x))
    lip <- 2 / sqrt(mean((x - mean(x))^2))
    jump <- LPM.ratio(2, xs + delta, x) - LPM.ratio(2, xs - delta, x)
    expect_lte(max(jump), 2 * delta * lip * (1 + 1e-6) + ratio_tol(x), label = nm)
  }
})

test_that("cdfCurve_continuous fails at degree 0: jump at each observation is its atom >= 1/n", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    xs <- sort(unique(x))
    delta <- min(diff(xs)) / 4
    jump <- LPM.ratio(0, xs, x) - LPM.ratio(0, xs - delta, x)
    expect_true(all(jump >= 1 / length(x) - 2 * eps), info = nm)
  }
})

test_that("cdfCurve_of_const / not_continuousAt_cdfCurve_of_const: constant sample CDF is a step", {
  for (d in 0:3) {
    expect_identical(cdf_fun(rep(5, 3), d)$Function$CDF, c(1, 1, 1), info = d)
    expect_identical(cdf_fun(5, d)$Function$CDF, 1, info = d)
    expect_identical(LPM.ratio(d, 5 - 1e-9, rep(5, 3)), 0, info = d)
  }
})

# ---------------------------------------------------------------------------
# CDF.lean: survival, cumulative hazard, target.value
# ---------------------------------------------------------------------------

test_that("survivalCurve: NNS.CDF survival column is 1 - CDF column", {
  for (nm in names(samples)) for (d in 0:3) {
    x <- samples[[nm]]
    s <- cdf_fun(x, d, "survival")$Function[["S(x)"]]
    expect_identical(s, 1 - cdf_fun(x, d)$Function$CDF, info = paste(nm, d))
  }
})

# survivalCurve_eq_UPMratio' / cdfCurve_eq_LPMratio (the CDF <-> Ratios bridge) appear only
# in the CDF.lean doc comment and are not compiled; the compiled statements are
# survivalCurve_eq_UPMratio (d > 0, `if all x = t then 0 else U/(L+U)`, which is exactly
# UPM.ratio's definition) and survivalCurve_zero_eq_UPM / survivalCurve_zero_eq_count (d = 0).
test_that("survivalCurve_eq_UPMratio / survivalCurve_zero_eq_UPM / survivalCurve_zero_eq_count: survival = UPM.ratio, = #{x > t}/n at d = 0", {
  for (nm in names(samples)) for (d in 0:3) {
    x <- samples[[nm]]
    s <- cdf_fun(x, d, "survival")$Function[["S(x)"]]
    expect_close(s, UPM.ratio(d, sort(x), x), 4 * eps, info = paste(nm, d))
    if (d == 0) expect_close(s, surv_count(sort(x), x), 2 * eps, info = nm)
  }
})

test_that("cumHazardCurve_eq: H = -log(max(S, 1e-10))", {
  for (nm in names(samples)) for (d in 0:3) {
    x <- samples[[nm]]
    h <- cdf_fun(x, d, "cumulative hazard")$Function[["H(x)"]]
    s <- 1 - LPM.ratio(d, sort(x), x)
    expect_close(h, -log(pmax(s, 1e-10)), 1e-12, info = paste(nm, d))
  }
})

test_that("cumHazardCurve_eq with survivalCurve_eq_UPMratio: H = -log(max(UPM.ratio, 1e-10))", {
  for (nm in names(samples)) for (d in 0:3) {
    x <- samples[[nm]]
    h <- cdf_fun(x, d, "cumulative hazard")$Function[["H(x)"]]
    s <- UPM.ratio(d, sort(x), x)
    # S is known to O(eps); -log amplifies an absolute error in S by 1 / S.
    expect_true(all(abs(h + log(pmax(s, 1e-10))) <= 8 * eps / pmax(s, 1e-10)), info = paste(nm, d))
  }
})

test_that("cumHazardCurve_nonneg / cumHazardCurve_le / cumHazardCurve_monotone: 0 <= H <= 10 log 10, nondecreasing", {
  for (nm in names(samples)) for (d in 0:3) {
    h <- cdf_fun(samples[[nm]], d, "cumulative hazard")$Function[["H(x)"]]
    expect_true(all(h >= 0 & h <= 10 * log(10) + 1e-12), info = paste(nm, d))
    expect_gte(min(diff(h)), -1e-12, label = paste(nm, d))
  }
})

test_that("cumHazardCurve_of_forall_le: at the sample maximum H is the cap 10 log 10", {
  for (nm in names(samples)) for (d in 0:3) {
    h <- cdf_fun(samples[[nm]], d, "cumulative hazard")$Function[["H(x)"]]
    expect_close(h[length(h)], 10 * log(10), 1e-12, info = paste(nm, d))
  }
})

test_that("cumHazardCurve_zero_eq_neg_log: degree 0, H = -log S wherever S > 0", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    h <- cdf_fun(x, 0, "cumulative hazard")$Function[["H(x)"]]
    s <- surv_count(sort(x), x)
    expect_close(h[s > 0], -log(s[s > 0]), 1e-12, info = nm)
  }
})

test_that("cumHazardCurve_ne_neg_log_example: d = 1, x = (0, 1), t = 1 - 1e-11 hits the floor", {
  t <- 1 - 1e-11
  s <- UPM.ratio(1, t, c(0, 1))
  expect_gt(s, 0)
  expect_lt(s, 1e-10)
  h <- cdf_fun(c(0, 1), 1, "cumulative hazard", target = t)$target.value
  expect_close(h, 10 * log(10), 1e-12)
})

test_that("cdfTarget_eq / survivalTarget_eq / cumHazardTarget_eq: target.value is the curve at target", {
  for (nm in names(samples)) for (d in 0:1) {
    x <- samples[[nm]]
    for (t in c(median(x), x[1])) {
      l <- LPM.ratio(d, t, x)
      expect_identical(cdf_fun(x, d, target = t)$target.value, l, info = paste(nm, d))
      expect_identical(cdf_fun(x, d, "survival", target = t)$target.value, 1 - l, info = paste(nm, d))
      expect_close(cdf_fun(x, d, "cumulative hazard", target = t)$target.value,
                   -log(max(1 - l, 1e-10)), 1e-12, info = paste(nm, d))
    }
  }
})

test_that("hazardFromDensity_mem_Icc: hazard output is clamped to [0, 1e6]", {
  for (x in list(c(1, 2, 2, 3, 4, 6, 7, 9, 10, 15), samples$mpg)) {
    h <- cdf_fun(x, 0, "hazard")$Function[["h(x)"]]
    expect_true(all(h >= 0 & h <= 1e6))
  }
})

# ---------------------------------------------------------------------------
# LPM.VaR / UPM.VaR: roots of the ratio
# ---------------------------------------------------------------------------

var_p <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

test_that("LPMratio_strictMonoOn (VaR root): LPM.ratio(d, LPM.VaR(p, d, x), x) = p, d = 1..4", {
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    v <- LPM.VaR(var_p, d, x)
    expect_close(LPM.ratio(d, v, x), var_p, ratio_tol(x), info = paste(nm, d))
  }
})

test_that("UPMratio_strictAntiOn (VaR root): UPM.ratio(d, UPM.VaR(p, d, x), x) = p, d = 1..4", {
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    v <- UPM.VaR(var_p, d, x)
    expect_close(UPM.ratio(d, v, x), var_p, ratio_tol(x), info = paste(nm, d))
  }
})

test_that("LPMratio_strictMonoOn (VaR root): fallback path d = 2.5, 5 within the uniroot tolerance", {
  # Fallback solves to sqrt(eps) * range in t (BUGS.md #15, design choice).
  for (nm in names(samples)) for (d in c(2.5, 5)) {
    x <- samples[[nm]]
    v <- LPM.VaR(var_p, d, x)
    expect_close(LPM.ratio(d, v, x), var_p, 1e-6, info = paste(nm, d))
  }
})

test_that("LPMratio_strictMonoOn (VaR root): exact dyadic probe with ties, z = (1,1,1,2,5,5)", {
  z <- c(1, 1, 1, 2, 5, 5)
  v <- LPM.VaR(c(0.1, 0.5, 0.9), 1, z)
  expect_close(v, c(1.3, 2.5, 4.318181818181818), 4 * eps * 5)
  expect_close(LPM.ratio(1, v, z), c(0.1, 0.5, 0.9), 4 * eps)
})

test_that("LPMratio_strictMonoOn (VaR uniqueness): the root at an order statistic is that order statistic", {
  # fails on 8898017: see BUGS.md #8 (median-anchored power sums cancel at d = 3, 4)
  set.seed(11)
  x <- c(rnorm(40, -1e5, 1), rnorm(61, 0, 1))
  target <- sort(x)[20]
  for (d in 3:4) {
    p <- LPM.ratio(d, target, x)
    expect_close(LPM.VaR(p, d, x), target, 1e-6, info = d)
  }
})

test_that("LPMratio_monotone (VaR): LPM.VaR is nondecreasing in p, UPM.VaR nonincreasing", {
  p <- seq(0, 1, by = 0.05)
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    expect_true(all(diff(LPM.VaR(p, d, x)) >= 0), info = paste(nm, d))
    expect_true(all(diff(UPM.VaR(p, d, x)) <= 0), info = paste(nm, d))
  }
})

test_that("LPMratio_eq_zero_of_forall_ge / LPMratio_eq_one_of_forall_le (VaR ends): p = 0 -> min, p = 1 -> max", {
  for (nm in names(samples)) for (d in 1:4) {
    x <- samples[[nm]]
    expect_identical(LPM.VaR(c(0, 1), d, x), c(min(x), max(x)), info = paste(nm, d))
    expect_identical(UPM.VaR(c(1, 0), d, x), c(min(x), max(x)), info = paste(nm, d))
  }
})

test_that("LPMratio_of_forall_eq (VaR): constant sample and n = 1 return the common value", {
  for (d in c(0, 1, 2, 2.5)) {
    expect_identical(unname(LPM.VaR(c(0.3, 0.6), d, rep(4, 3))), c(4, 4), info = d)
    expect_identical(unname(LPM.VaR(0.3, d, 7)), 7, info = d)
  }
})

test_that("quantile7_eq / UPMVaR0_eq_LPMVaR0 (VaRRoot.lean): degree-0 VaR is the type-7 sample quantile", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    expect_identical(LPM.VaR(var_p, 0, x), quantile(x, var_p, type = 7), info = nm)
    expect_identical(unname(UPM.VaR(var_p, 0, x)), unname(quantile(x, 1 - var_p, type = 7)), info = nm)
  }
})

# At p = k/(n-1) the type-7 index (n-1)p = k is an integer: floor_add_one_le_card_le_quantile7
# gives ECDF(VaR) >= (k+1)/n >= p and card_lt_quantile7_le_ceil gives #{x < VaR} <= k, so the
# ECDF just below VaR is <= k/n < p: VaR is the minimal t with ECDF(t) >= p there.
test_that("floor_add_one_le_card_le_quantile7 / card_lt_quantile7_le_ceil (d = 0): at p = k/(n-1) VaR is the minimal t with ECDF(t) >= p", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    n <- length(x)
    p <- unique(round(seq(1, n - 2, length.out = 9))) / (n - 1)
    p <- p[(n - 1) * p == round((n - 1) * p)]  # type-7 index exactly an integer in double
    v <- unname(LPM.VaR(p, 0, x))
    below <- vapply(v, function(s) max(c(min(x) - 1, x[x < s])), numeric(1))
    expect_true(all(LPM.ratio(0, v, x) >= p), info = nm)
    expect_true(all(LPM.ratio(0, below, x) < p), info = nm)
  }
})

test_that("LPMratio_zero_quantile7_gt / floor_add_one_le_card_le_quantile7 / card_lt_quantile7_le_ceil (d = 0): type-7 VaR count bounds", {
  for (nm in names(samples)) {
    x <- samples[[nm]]
    n <- length(x)
    v <- unname(LPM.VaR(var_p, 0, x))
    h <- (n - 1) * var_p  # the type-7 index, computed in double exactly as quantile() does
    n_le <- vapply(v, function(s) sum(x <= s), numeric(1))
    n_lt <- vapply(v, function(s) sum(x < s), numeric(1))
    expect_true(all(n_le >= floor(h) + 1), info = nm)                  # floor_add_one_le_card_le_quantile7
    expect_true(all(LPM.ratio(0, v, x) > var_p - var_p / n), info = nm) # LPMratio_zero_quantile7_gt
    expect_true(all(n_lt <= ceiling(h)), info = nm)                      # card_lt_quantile7_le_ceil
  }
})

test_that("LPMVaR0_counterexample / LPMVaR0_not_minimal: type 7 is not the ECDF inverse on (0, 1, 2)", {
  x <- c(0, 1, 2)
  v <- unname(LPM.VaR(c(0.1, 0.9), 0, x))
  expect_close(v, c(0.2, 1.8), 4 * eps)
  expect_identical(LPM.ratio(0, v[2], x), 2 / 3)   # 2/3 < 0.9: ECDF(VaR) >= p fails
  expect_identical(LPM.ratio(0, 0, x), 1 / 3)      # 0 < VaR(0.1) yet ECDF(0) >= 0.1
})
