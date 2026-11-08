# LPM.ratio, UPM.ratio and univariate NNS.CDF. Lean: nns-proofs Ratios.lean, CDF.lean.

# Lean: LPMratio_mem_Icc, UPMratio_mem_Icc, LPMratio_add_UPMratio, UPMratio_eq_one_sub_LPMratio, LPMratio_mem_Ioo | BUGS: none
test_that("R01: LPM.ratio and UPM.ratio are complementary probabilities in [0, 1] on both target paths", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_mem_Icc / UPMratio_mem_Icc: both ratios lie in [0, 1]", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      t <- c(aud_targets_for(x), aud_grid_for(x))
      for (r in list(LPM.ratio(d, t, x), UPM.ratio(d, t, x))) {
        expect_true(all(r >= 0 & r <= 1), info = paste(nm, d))
      }
    }
  })
  aud_case("LPMratio_add_UPMratio: degree 0 sums to 1 exactly on dyadic data", {
    t <- c(-3, fr$x_exact, 0.1, 2, 4)
    expect_identical(LPM.ratio(0, t, fr$x_exact) + UPM.ratio(0, t, fr$x_exact), rep(1, length(t)))
  })
  aud_case("UPMratio_eq_one_sub_LPMratio: UPM.ratio = 1 - LPM.ratio at every degree", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      t <- aud_targets_for(x)
      aud_close(UPM.ratio(d, t, x), 1 - LPM.ratio(d, t, x), 4 * aud_eps, info = paste(nm, d))
    }
  })
  aud_case("UPMratio_eq_one_sub_LPMratio: also on the > 32 target prefix path", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      t <- aud_grid_for(x)
      aud_close(UPM.ratio(d, t, x), 1 - LPM.ratio(d, t, x), 4 * aud_eps, info = paste(nm, d))
    }
  })
  aud_case("LPMratio_mem_Ioo: d > 0 and a target strictly inside the range give 0 < ratio < 1", {
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      t <- aud_targets_for(x)
      t <- t[t > min(x) & t < max(x)]
      r <- LPM.ratio(d, t, x)
      expect_true(all(r > 0 & r < 1), info = paste(nm, d))
    }
  })
})

# Lean: LPMratio_zero_eq_ecdf, UPMratio_zero_eq_survival | BUGS: none
test_that("R02: degree-0 ratios are the ECDF (ties low) and the strict survival share on both target paths", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_zero_eq_ecdf: degree-0 LPM.ratio is #{x <= t}/n (ties low)", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      t <- c(aud_targets_for(x), x)  # exact data points: ties at the target
      expect_identical(LPM.ratio(0, t, x), aud_ecdf_count(t, x), info = nm)
    }
  })
  aud_case("LPMratio_zero_eq_ecdf: prefix path (> 32 targets) is the ECDF too", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      t <- aud_grid_for(x)
      expect_identical(LPM.ratio(0, t, x), aud_ecdf_count(t, x), info = nm)
    }
  })
  aud_case("UPMratio_zero_eq_survival: degree-0 UPM.ratio is #{x > t}/n (strict)", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      t <- c(aud_targets_for(x), x)
      expect_identical(UPM.ratio(0, t, x), aud_surv_count(t, x), info = nm)
    }
  })
  aud_case("LPMratio (elementwise target): > 32 targets agree with scalar calls", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      t <- aud_grid_for(x)
      scalar <- vapply(t, function(s) LPM.ratio(d, s, x), numeric(1))
      aud_close(LPM.ratio(d, t, x), scalar, aud_ratio_tol(x), info = paste(nm, d))
    }
  })
})

# Lean: LPMratio_eq_sum_div, LPMratio_one_mean, UPMratio_one_mean, LPMratioDefault_one, LPMratio_one_eq_half_iff | BUGS: #12
# These cases also pass before the fix; the #12 failures there are in N04, C12, C15, R08 and V04.
test_that("R03: degree > 0 LPM.ratio is L / (L + U), exactly 1/2 at the mean, and a NULL target is the mean", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPM.ratio on ordinary data is still the plain L / (L + U)", {
    withr::local_seed(2)
    r <- rnorm(1000, 0, 0.02)
    expect_identical(LPM.ratio(2, 0, r), LPM(2, 0, r) / (LPM(2, 0, r) + UPM(2, 0, r)))
  })
  aud_case("LPMratio_eq_sum_div: d > 0 ratio is sum (t-x)_+^d / sum |x-t|^d", {
    for (nm in names(fr$kept)) for (d in 1:3) {
      x <- fr$kept[[nm]]
      t <- aud_targets_for(x)
      t <- t[t > min(x)]  # a target off the sample (or with mass on both sides)
      manual <- vapply(t, function(s) sum(pmax(s - x, 0)^d) / sum(abs(x - s)^d), numeric(1))
      aud_close(LPM.ratio(d, t, x), manual, aud_ratio_tol(x), info = paste(nm, d))
    }
  })
  aud_case("LPMratioDefault: NULL target is the sample mean", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      x <- fr$kept[[nm]]
      expect_identical(LPM.ratio(d, NULL, x), LPM.ratio(d, mean(x), x), info = paste(nm, d))
    }
  })
  aud_case("LPMratio_one_mean: exact dyadic data, LPM.ratio(1, mean) is exactly 1/2", {
    expect_identical(mean(fr$x_exact), 0.34375)
    expect_identical(LPM.ratio(1, mean(fr$x_exact), fr$x_exact), 0.5)
    expect_identical(UPM.ratio(1, mean(fr$x_exact), fr$x_exact), 0.5)
  })
  aud_case("LPMratio_one_mean / UPMratio_one_mean: 1/2 at the mean on random and real data", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      aud_close(LPM.ratio(1, mean(x), x), 0.5, aud_ratio_tol(x), info = nm)
      aud_close(UPM.ratio(1, mean(x), x), 0.5, aud_ratio_tol(x), info = nm)
    }
  })
  aud_case("LPMratioDefault_one: NULL target at degree 1 gives 1/2 on non-constant data", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      aud_close(LPM.ratio(1, NULL, x), 0.5, aud_ratio_tol(x), info = nm)
    }
  })
  aud_case("LPMratio_one_eq_half_iff: away from the mean the degree-1 ratio is not 1/2 (sign of t - mean)", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      t <- aud_targets_for(x)
      t <- t[abs(t - mean(x)) > 1e-6 * diff(range(x)) & t > min(x) & t < max(x)]
      expect_identical(sign(LPM.ratio(1, t, x) - 0.5), sign(t - mean(x)), info = nm)
    }
  })
})

# Lean: LPMratio_of_forall_eq, UPMratio_of_forall_eq, LPMratioDefault_of_const | BUGS: none
test_that("R04: constant and one-point samples follow the 0/0 step convention", {
  withr::local_seed(20261108)
  aud_case("LPMratio_of_forall_eq / UPMratio_of_forall_eq: constant sample at its value gives 1 / 0", {
    for (d in 1:4) for (x in list(rep(5, 3), rep(-0.25, 8), 7)) {
      expect_identical(LPM.ratio(d, x[1], x), 1, info = paste(d, length(x)))
      expect_identical(UPM.ratio(d, x[1], x), 0, info = paste(d, length(x)))
    }
  })
  aud_case("LPMratio_of_forall_eq: constant sample is the step 1[c <= t] at every degree", {
    for (d in 0:3) {
      expect_identical(LPM.ratio(d, c(4, 5, 6), rep(5, 3)), c(0, 1, 1), info = d)
      expect_identical(UPM.ratio(d, c(4, 5, 6), rep(5, 3)), c(1, 0, 0), info = d)
    }
  })
  aud_case("LPMratio_of_forall_eq: n = 1 is constant, step at the single value", {
    for (d in 0:3) {
      expect_identical(LPM.ratio(d, c(0, 1, 2), 1), c(0, 1, 1), info = d)
      expect_identical(UPM.ratio(d, c(0, 1, 2), 1), c(1, 0, 0), info = d)
    }
  })
  aud_case("LPMratioDefault_of_const: NULL target on a constant sample gives 1, not 1/2", {
    for (d in 0:3) {
      expect_identical(LPM.ratio(d, NULL, rep(2.5, 4)), 1, info = d)
      expect_identical(UPM.ratio(d, NULL, rep(2.5, 4)), 0, info = d)
    }
  })
  aud_case("the 0/0 fill at degree > 0, scalar and vector targets", {
    expect_equal(LPM.ratio(1, 3, rep(3, 4)), 1)
    expect_equal(UPM.ratio(2, 3, rep(3, 4)), 0)
    expect_equal(LPM.ratio(1, rep(3, 40), 3), rep(1, 40))
    expect_equal(NNS.CDF(rep(2, 4), 1, plot = FALSE)$Function$CDF, rep(1, 4))
    na_or_error <- function(expr) tryCatch(all(is.na(expr)), error = function(e) TRUE)
    expect_true(na_or_error(LPM.ratio(NA_real_, 2, 1:3)))
    expect_true(na_or_error(UPM.ratio(NA_real_, 2, 1:3)))
  })
})

# Lean: LPMratio_eq_one_of_forall_le, LPMratio_eq_zero_of_forall_lt, LPMratio_eq_zero_of_forall_ge, UPMratio_eq_one_of_forall_ge | BUGS: none
test_that("R05: ratios at and beyond the sample range are exactly 0 or 1, with the degree-0 atom", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_eq_one_of_forall_le: at or above the maximum LPM.ratio = 1 exactly", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      expect_identical(LPM.ratio(d, c(max(x), max(x) + 1), x), c(1, 1), info = paste(nm, d))
    }
  })
  aud_case("LPMratio_eq_zero_of_forall_lt: strictly below the minimum LPM.ratio = 0 exactly", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      expect_identical(LPM.ratio(d, min(x) - 1, x), 0, info = paste(nm, d))
    }
  })
  aud_case("LPMratio_eq_zero_of_forall_ge: d > 0 at the minimum gives 0, degree 0 gives the atom", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      for (d in 1:4) expect_identical(LPM.ratio(d, min(x), x), 0, info = paste(nm, d))
      expect_identical(LPM.ratio(0, min(x), x), sum(x == min(x)) / length(x), info = nm)
    }
  })
  aud_case("UPMratio_eq_one_of_forall_ge: d > 0 at the minimum UPM.ratio = 1, degree 0 drops the atom", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      for (d in 1:4) expect_identical(UPM.ratio(d, min(x), x), 1, info = paste(nm, d))
      expect_identical(UPM.ratio(0, min(x), x), sum(x > min(x)) / length(x), info = nm)
    }
  })
})

# Lean: LPMratio_monotone, UPMratio_antitone, LPMratio_strictMonoOn | BUGS: none
test_that("R06: LPM.ratio is monotone in the target, strictly so only for d > 0", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_monotone: exact dyadic data, nondecreasing with zero tolerance", {
    t <- seq(-3, 4, by = 0.125)
    for (d in 0:3) expect_true(all(diff(LPM.ratio(d, t, fr$x_exact)) >= 0), info = d)
  })
  aud_case("LPMratio_monotone / UPMratio_antitone: random and real data, both target paths", {
    for (nm in names(fr$kept)) for (d in 0:4) {
      x <- fr$kept[[nm]]
      t <- aud_grid_for(x)
      tol <- aud_ratio_tol(x)
      expect_gte(min(diff(LPM.ratio(d, t, x))), -tol, label = paste("LPM", nm, d))
      expect_lte(max(diff(UPM.ratio(d, t, x))), tol, label = paste("UPM", nm, d))
      scalar <- vapply(t, function(s) LPM.ratio(d, s, x), numeric(1))
      expect_gte(min(diff(scalar)), -tol, label = paste("scalar", nm, d))
    }
  })
  aud_case("LPMratio_strictMonoOn: d > 0 strictly increasing on [min x, max x]", {
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      t <- sort(unique(c(x, seq(min(x), max(x), length.out = 40))))
      t <- t[c(TRUE, diff(t) > 1e-9 * diff(range(x)))]  # drop grid points within rounding of a datum
      expect_true(all(diff(LPM.ratio(d, t, x)) > 0), info = paste(nm, d))
    }
  })
  aud_case("LPMratio_strictMonoOn fails at degree 0: the ECDF is flat between observations", {
    xs <- sort(unique(fr$x_exact))
    mid <- (xs[-1] + xs[-length(xs)]) / 2
    expect_identical(LPM.ratio(0, mid, fr$x_exact), LPM.ratio(0, xs[-length(xs)], fr$x_exact))
  })
})

# Lean: LPMratio_affine, LPMratio_eq_UPMratio_neg, LPMratio_ne_UPMratio_neg_of_forall_eq, LPMratio_zero_eq_UPMratio_neg_add | BUGS: none
test_that("R07: ratios are invariant under affine maps and reflect with the degree-0 atom", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_affine: exact dyadic data, power-of-two scale and dyadic shift", {
    t <- c(-1, 0, 0.34375, 1, 2)
    for (d in 0:3) {
      expect_identical(LPM.ratio(d, 4 * t + 0.5, 4 * fr$x_exact + 0.5), LPM.ratio(d, t, fr$x_exact), info = d)
      expect_identical(UPM.ratio(d, 4 * t + 0.5, 4 * fr$x_exact + 0.5), UPM.ratio(d, t, fr$x_exact), info = d)
    }
  })
  aud_case("LPMratio_affine: random and real data under a = 3.7, b = -1.3", {
    for (nm in names(fr$kept)) for (d in 1:3) {
      x <- fr$kept[[nm]]
      t <- aud_targets_for(x)
      t <- t[t != min(x)]  # keep the target off the affine-rounded minimum
      aud_close(LPM.ratio(d, 3.7 * t - 1.3, 3.7 * x - 1.3), LPM.ratio(d, t, x),
                1e-12, info = paste(nm, d))
    }
  })
  aud_case("LPMratio_eq_UPMratio_neg: d > 0 reflection on exact dyadic data", {
    t <- c(-1, 0, 0.34375, 1, 2)  # 0 and 1 are tied observations
    for (d in 1:3) {
      expect_identical(LPM.ratio(d, t, fr$x_exact), UPM.ratio(d, -t, -fr$x_exact), info = d)
    }
  })
  aud_case("LPMratio_eq_UPMratio_neg: d > 0 reflection on random and real data, ties at target", {
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      t <- aud_targets_for(x)
      aud_close(LPM.ratio(d, t, x), UPM.ratio(d, -t, -x), 4 * aud_eps, info = paste(nm, d))
    }
  })
  aud_case("LPMratio_ne_UPMratio_neg_of_forall_eq: reflection breaks in the degenerate case", {
    for (d in 1:3) {
      expect_identical(LPM.ratio(d, 1, c(1, 1, 1)), 1, info = d)
      expect_identical(UPM.ratio(d, -1, -c(1, 1, 1)), 0, info = d)
    }
  })
  aud_case("LPMratio_zero_eq_UPMratio_neg_add: degree-0 reflection picks up the atom (exact)", {
    t <- c(-3, fr$x_exact, 0.1, 4)
    atom <- vapply(t, function(s) sum(fr$x_exact == s), numeric(1)) / length(fr$x_exact)
    expect_identical(LPM.ratio(0, t, fr$x_exact), UPM.ratio(0, -t, -fr$x_exact) + atom)
  })
  aud_case("LPMratio_zero_eq_UPMratio_neg_add: degree-0 reflection with atom, random and real data", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      t <- c(aud_targets_for(x), x[1:5])
      atom <- vapply(t, function(s) sum(x == s), numeric(1)) / length(x)
      aud_close(LPM.ratio(0, t, x), UPM.ratio(0, -t, -x) + atom, 2 * aud_eps, info = nm)
    }
  })
})

# Lean: LPMratio_affine | BUGS: #5, #12
test_that("R08: BUGS #12: ratios are scale-free and NaN-free at extreme degree and scale", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPM.ratio on the prefix path is 0, not NaN, when only the upper sum overflows", {
    y <- c(0, 1e200)
    expect_identical(LPM.ratio(4, c(-1e-170, aud_pad), y)[1], LPM.ratio(4, -1e-170, y))
    expect_identical(LPM.ratio(4, c(-1e-170, aud_pad), y)[1], 0)
  })
  aud_case("UPM.ratio on the prefix path is 1, not NaN (Inf/Inf), when the upper sum overflows", {
    # Every observation lies above the target, so U/(L+U) = 1 exactly; the
    # overflowed raw sums are recomputed on deviations scaled by max|x - t|.
    y <- c(0, 1e200)
    expect_identical(UPM.ratio(4, c(-1e-170, aud_pad), y)[1], UPM.ratio(4, -1e-170, y))
    expect_identical(UPM.ratio(4, c(-1e-170, aud_pad), y)[1], 1)
  })
  aud_case("univariate NNS.CDF over 32 observations has no NaN when a power sum overflows", {
    d <- c(0, 1e200, seq(1, 2, length.out = 40))
    f <- NNS.CDF(d, degree = 4, plot = FALSE)$Function
    expect_identical(f$CDF[1:3], vapply(f$x[1:3], function(t) LPM.ratio(4, t, d), numeric(1)))
    expect_identical(f$CDF[1:3], c(0, 0, 0))
  })
  aud_case("UPM.ratio is 0.5, not NaN, when every deviation overflows", {
    expect_identical(UPM.ratio(400, 0, c(-1e3, 1e3)), 0.5)
  })
  aud_case("LPM.ratio does not depend on the units of the data at high degree", {
    x <- c(-2, -1, 0.5, 3)
    expect_equal(LPM.ratio(150, 0, 1e-3 * x), LPM.ratio(150, 0, x), tolerance = 1e-12)
    expect_equal(LPM.ratio(150, 0, 1e3 * x), LPM.ratio(150, 0, x), tolerance = 1e-12)
  })
  aud_case("LPM.ratio and UPM.ratio on more than 32 targets have no NaN at degree 300", {
    withr::local_seed(1)
    r <- rnorm(250, 0, 0.01)
    expect_false(anyNA(LPM.ratio(300, r, r)))
    expect_equal(LPM.ratio(300, r, r) + UPM.ratio(300, r, r), rep(1, 250))
  })
  aud_case("univariate NNS.CDF at degree 400 has no NaN on daily-return data", {
    withr::local_seed(1)
    r <- rnorm(250, 0, 0.01)
    expect_false(anyNA(NNS.CDF(r, degree = 400, plot = FALSE)$Function$CDF))
  })
  aud_case("BUGS #12 regression (no Lean theorem; under/overflow is outside the exact-real model of Ratios.lean): high degree on symmetric data gives 1/2, not NaN", {
    # Symmetric two-point sample: normalised deviations are +/-1, so 0.5 is exact.
    # fails on 8898017: see BUGS.md #12 (L and U underflow / overflow to 0/0 or Inf/Inf -> NaN)
    expect_identical(LPM.ratio(400, 0, c(-1e-2, 1e-2)), 0.5)
    expect_identical(LPM.ratio(400, 0, c(-1e3, 1e3)), 0.5)
  })
  aud_case("BUGS #12 regression (no Lean theorem; LPMratio_affine in exact reals): degree 2 is scale-free under a tiny power-of-two scale", {
    # fails on 8898017: see BUGS.md #12 (deviations^2 underflow at scale 2^-600 -> NaN)
    # A fix that rescales (e.g. by max |x - t|) rounds differently at the non-dyadic
    # targets (quantiles -1.475, 2.3): ~1e-16 apart, so compare to 4 eps, not bitwise.
    a <- 2^-600
    t <- aud_targets_for(fr$x_exact)
    aud_close(LPM.ratio(2, a * t, a * fr$x_exact), LPM.ratio(2, t, fr$x_exact), 4 * aud_eps)
  })
})

# Lean: cdfCurve_probe_zero, cdfCurve_zero_eq_ecdf, cdfTable_zero_ge, cdfTable_zero_of_injective, cdfCurve_zero_eventually_const_right, cdfCurve_zero_eq_lt_add_eq, cdfCurve_zero_tendsto_left, cdfCurve_continuous | BUGS: none
test_that("R09: univariate NNS.CDF at degree 0 is the right-continuous ECDF with atom jumps", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("cdfCurve_probe_zero: NNS.CDF(c(3,1,2,2,10)) is 0.2 0.6 0.6 0.8 1", {
    f <- aud_cdf_fun(c(3, 1, 2, 2, 10))$Function
    expect_identical(f$x, c(1, 2, 2, 3, 10))
    expect_identical(f$CDF, c(1, 3, 3, 4, 5) / 5)
  })
  aud_case("cdfCurve_zero_eq_ecdf: degree-0 NNS.CDF equals stats::ecdf at the sorted sample", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      f <- aud_cdf_fun(x)$Function
      expect_identical(f$x, sort(x), info = nm)
      expect_identical(f$CDF, stats::ecdf(x)(sort(x)), info = nm)
    }
  })
  aud_case("cdfTable_zero_ge: degree-0 row k has CDF >= k/n (ties push it up)", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      n <- length(x)
      f <- aud_cdf_fun(x)$Function
      expect_true(all(f$CDF >= seq_len(n) / n), info = nm)
    }
  })
  aud_case("cdfTable_zero_of_injective: tie-free sample, degree-0 row k has CDF = k/n exactly", {
    tie_free <- list(
      exact = c(3, -2, 0.25, -0.5, 1, 0),             # dyadic, n = 6 (k/6 not dyadic)
      dyad8 = c(-2, -0.5, 0, 0.125, 0.25, 1, 2, 3),    # n = 8: k/8 exact
      dax   = unique(fr$samples$dax),                     # large n, real data
      norm  = unique(fr$samples$norm)
    )
    for (nm in names(tie_free)) {
      x <- tie_free[[nm]]
      n <- length(x)
      expect_false(as.logical(anyDuplicated(x)), info = nm)  # the hypothesis really holds
      expect_gt(n, 1)
      f <- aud_cdf_fun(x)$Function
      expect_identical(f$x, sort(x), info = nm)
      expect_identical(f$CDF, seq_len(n) / n, info = nm)
    }
  })
  aud_case("cdfCurve_zero_eventually_const_right: the degree-0 CDF is right-continuous", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      xs <- sort(unique(x))
      delta <- min(diff(xs)) / 4
      expect_identical(LPM.ratio(0, xs + delta, x), LPM.ratio(0, xs, x), info = nm)
    }
  })
  aud_case("cdfCurve_zero_tendsto_left / cdfCurve_zero_eq_lt_add_eq: left limit is #{x < t}/n, jump is the atom", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      xs <- sort(unique(x))
      delta <- min(diff(xs)) / 4
      left <- vapply(xs, function(s) sum(x < s), numeric(1)) / length(x)
      atom <- vapply(xs, function(s) sum(x == s), numeric(1)) / length(x)
      expect_identical(LPM.ratio(0, xs - delta, x), left, info = nm)
      aud_close(LPM.ratio(0, xs, x) - LPM.ratio(0, xs - delta, x), atom, 2 * aud_eps, info = nm)
    }
  })
  aud_case("cdfCurve_continuous fails at degree 0: jump at each observation is its atom >= 1/n", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      xs <- sort(unique(x))
      delta <- min(diff(xs)) / 4
      jump <- LPM.ratio(0, xs, x) - LPM.ratio(0, xs - delta, x)
      expect_true(all(jump >= 1 / length(x) - 2 * aud_eps), info = nm)
    }
  })
})

# Lean: cdfCurve_probe_one, cdfCurve_probe_tied_min, cdfTarget_probe, cdfTable_snd_monotone, cdfTable_snd_last, cdfTable_snd_first, cdfTable_snd_eq_zero_of_tied_first, cdfCurve_lt_cdfCurve, cdfCurve_continuous, cdfCurve_of_const, not_continuousAt_cdfCurve_of_const | BUGS: none
test_that("R10: univariate NNS.CDF at degree > 0 rises continuously and strictly from 0 to 1", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("cdfCurve_probe_one / cdfCurve_probe_tied_min / cdfTarget_probe", {
    aud_close(LPM.ratio(1, 2.5, c(1, 2, 2, 3, 10)), 5 / 21, 2 * aud_eps)
    aud_close(aud_cdf_fun(c(1, 1, 2, 5), 1)$Function$CDF, c(0, 0, 2 / 5, 1), 2 * aud_eps)
    aud_close(aud_cdf_fun(c(1, 1, 2, 5), 1, target = 1.5)$target.value, 1 / 5, 2 * aud_eps)
  })
  aud_case("cdfTable_snd_monotone / cdfTable_snd_last: CDF column nondecreasing, last row 1", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      f <- aud_cdf_fun(fr$kept[[nm]], d)$Function
      expect_gte(min(diff(f$CDF)), -aud_ratio_tol(f$x), label = paste(nm, d))
      expect_identical(f$CDF[nrow(f)], 1, info = paste(nm, d))
    }
  })
  aud_case("cdfTable_snd_first / cdfTable_snd_eq_zero_of_tied_first: d > 0 rows at the minimum are 0", {
    for (nm in names(fr$kept)) for (d in 1:3) {
      x <- fr$kept[[nm]]
      f <- aud_cdf_fun(x, d)$Function
      expect_identical(f$CDF[f$x == min(x)], rep(0, sum(x == min(x))), info = paste(nm, d))
    }
  })
  aud_case("cdfCurve_lt_cdfCurve: d > 0 strictly increasing across distinct rows, equal on tied rows", {
    for (nm in names(fr$kept)) for (d in 1:3) {
      f <- aud_cdf_fun(fr$kept[[nm]], d)$Function
      step <- diff(f$CDF)
      tied <- diff(f$x) == 0
      expect_true(all(step[!tied] > 0), info = paste(nm, d))
      expect_true(all(step[tied] == 0), info = paste(nm, d))
    }
  })
  aud_case("cdfCurve_continuous (quantitative, test-derived bound): d = 1 CDF is Lipschitz with constant 1 / mean|x - median|", {
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
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      xs <- sort(unique(x))
      delta <- 1e-6 * diff(range(x))
      lip <- 1 / mean(abs(x - median(x)))
      jump <- LPM.ratio(1, xs + delta, x) - LPM.ratio(1, xs - delta, x)
      expect_lte(max(jump), 2 * delta * lip * (1 + 1e-6) + aud_ratio_tol(x), label = nm)
    }
  })
  aud_case("cdfCurve_continuous (quantitative, test-derived bound): d = 2 CDF is Lipschitz with constant 2 / sd_pop", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      xs <- sort(unique(x))
      delta <- 1e-6 * diff(range(x))
      lip <- 2 / sqrt(mean((x - mean(x))^2))
      jump <- LPM.ratio(2, xs + delta, x) - LPM.ratio(2, xs - delta, x)
      expect_lte(max(jump), 2 * delta * lip * (1 + 1e-6) + aud_ratio_tol(x), label = nm)
    }
  })
  aud_case("cdfCurve_of_const / not_continuousAt_cdfCurve_of_const: constant sample CDF is a step", {
    for (d in 0:3) {
      expect_identical(aud_cdf_fun(rep(5, 3), d)$Function$CDF, c(1, 1, 1), info = d)
      expect_identical(aud_cdf_fun(5, d)$Function$CDF, 1, info = d)
      expect_identical(LPM.ratio(d, 5 - 1e-9, rep(5, 3)), 0, info = d)
    }
  })
})

# Lean: cdfCurve_eq_LPMratio, survivalCurve_eq_UPMratio', survivalCurve_zero_eq_UPM, survivalCurve_zero_eq_count, survivalCurve_eq_UPMratio, cumHazardCurve_eq, cumHazardCurve_nonneg, cumHazardCurve_le, cumHazardCurve_monotone, cumHazardCurve_of_forall_le, cumHazardCurve_zero_eq_neg_log, cumHazardCurve_ne_neg_log_example, cdfTarget_eq, survivalTarget_eq, cumHazardTarget_eq, hazardFromDensity_mem_Icc | BUGS: none
test_that("R11: NNS.CDF survival, cumulative hazard, hazard and target.value agree with the CDF curve", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("survivalCurve: NNS.CDF survival column is 1 - CDF column", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      x <- fr$kept[[nm]]
      s <- aud_cdf_fun(x, d, "survival")$Function[["S(x)"]]
      expect_identical(s, 1 - aud_cdf_fun(x, d)$Function$CDF, info = paste(nm, d))
    }
  })
  aud_case("survivalCurve_eq_UPMratio / survivalCurve_zero_eq_UPM / survivalCurve_zero_eq_count: survival = UPM.ratio, = #{x > t}/n at d = 0", {
    # survivalCurve_eq_UPMratio' / cdfCurve_eq_LPMratio (the CDF <-> Ratios bridge) appear only
    # in the CDF.lean doc comment and are not compiled; the compiled statements are
    # survivalCurve_eq_UPMratio (d > 0, `if all x = t then 0 else U/(L+U)`, which is exactly
    # UPM.ratio's definition) and survivalCurve_zero_eq_UPM / survivalCurve_zero_eq_count (d = 0).
    for (nm in names(fr$kept)) for (d in 0:3) {
      x <- fr$kept[[nm]]
      s <- aud_cdf_fun(x, d, "survival")$Function[["S(x)"]]
      aud_close(s, UPM.ratio(d, sort(x), x), 4 * aud_eps, info = paste(nm, d))
      if (d == 0) aud_close(s, aud_surv_count(sort(x), x), 2 * aud_eps, info = nm)
    }
  })
  aud_case("cumHazardCurve_eq: H = -log(max(S, 1e-10))", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      x <- fr$kept[[nm]]
      h <- aud_cdf_fun(x, d, "cumulative hazard")$Function[["H(x)"]]
      s <- 1 - LPM.ratio(d, sort(x), x)
      aud_close(h, -log(pmax(s, 1e-10)), 1e-12, info = paste(nm, d))
    }
  })
  aud_case("cumHazardCurve_eq with survivalCurve_eq_UPMratio: H = -log(max(UPM.ratio, 1e-10))", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      x <- fr$kept[[nm]]
      h <- aud_cdf_fun(x, d, "cumulative hazard")$Function[["H(x)"]]
      s <- UPM.ratio(d, sort(x), x)
      # S is known to O(eps); -log amplifies an absolute error in S by 1 / S.
      expect_true(all(abs(h + log(pmax(s, 1e-10))) <= 8 * aud_eps / pmax(s, 1e-10)), info = paste(nm, d))
    }
  })
  aud_case("cumHazardCurve_nonneg / cumHazardCurve_le / cumHazardCurve_monotone: 0 <= H <= 10 log 10, nondecreasing", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      h <- aud_cdf_fun(fr$kept[[nm]], d, "cumulative hazard")$Function[["H(x)"]]
      expect_true(all(h >= 0 & h <= 10 * log(10) + 1e-12), info = paste(nm, d))
      expect_gte(min(diff(h)), -1e-12, label = paste(nm, d))
    }
  })
  aud_case("cumHazardCurve_of_forall_le: at the sample maximum H is the cap 10 log 10", {
    for (nm in names(fr$kept)) for (d in 0:3) {
      h <- aud_cdf_fun(fr$kept[[nm]], d, "cumulative hazard")$Function[["H(x)"]]
      aud_close(h[length(h)], 10 * log(10), 1e-12, info = paste(nm, d))
    }
  })
  aud_case("cumHazardCurve_zero_eq_neg_log: degree 0, H = -log S wherever S > 0", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      h <- aud_cdf_fun(x, 0, "cumulative hazard")$Function[["H(x)"]]
      s <- aud_surv_count(sort(x), x)
      aud_close(h[s > 0], -log(s[s > 0]), 1e-12, info = nm)
    }
  })
  aud_case("cumHazardCurve_ne_neg_log_example: d = 1, x = (0, 1), t = 1 - 1e-11 hits the floor", {
    t <- 1 - 1e-11
    s <- UPM.ratio(1, t, c(0, 1))
    expect_gt(s, 0)
    expect_lt(s, 1e-10)
    h <- aud_cdf_fun(c(0, 1), 1, "cumulative hazard", target = t)$target.value
    aud_close(h, 10 * log(10), 1e-12)
  })
  aud_case("cdfTarget_eq / survivalTarget_eq / cumHazardTarget_eq: target.value is the curve at target", {
    for (nm in names(fr$kept)) for (d in 0:1) {
      x <- fr$kept[[nm]]
      for (t in c(median(x), x[1])) {
        l <- LPM.ratio(d, t, x)
        expect_identical(aud_cdf_fun(x, d, target = t)$target.value, l, info = paste(nm, d))
        expect_identical(aud_cdf_fun(x, d, "survival", target = t)$target.value, 1 - l, info = paste(nm, d))
        aud_close(aud_cdf_fun(x, d, "cumulative hazard", target = t)$target.value,
                  -log(max(1 - l, 1e-10)), 1e-12, info = paste(nm, d))
      }
    }
  })
  aud_case("hazardFromDensity_mem_Icc: hazard output is clamped to [0, 1e6]", {
    for (x in list(c(1, 2, 2, 3, 4, 6, 7, 9, 10, 15), fr$samples$mpg)) {
      h <- aud_cdf_fun(x, 0, "hazard")$Function[["h(x)"]]
      expect_true(all(h >= 0 & h <= 1e6))
    }
  })
})

# pin: NNS-Python reference values
test_that("R12: NNS-Python reference values: LPM.ratio, UPM.ratio and the survival curve", {
  withr::local_seed(20261108)
  x <- nns_py$x
  x_df <- as.data.frame(x)
  ref <- data.frame(fn = rep(c("LPM.ratio", "UPM.ratio"), each = 3), degree = rep(c(0, 1, 2), 2),
                    expected = c(0.49, 0.5000000000000002, 0.49720627, 0.51, 0.4999999999999999, 0.5027937984146681))
  ref <- rbind(cbind(ref, input = "vector"), cbind(ref, input = "frame"))
  ref$label <- sprintf("%s(degree = %g) at the mean, %s input", ref$fn, ref$degree, ref$input)
  aud_cases(ref, function(r) {
    got <- if (r$input == "vector") {
      match.fun(r$fn)(degree = r$degree, target = mean(x), variable = x)
    } else {
      match.fun(r$fn)(degree = r$degree, target = colMeans(x_df), variable = x_df)
    }
    expect_equal(got, r$expected, tolerance = 1e-5)
  })
  aud_case("NNS.CDF survival curve of a sample with ties", {
    A <- c(1,1,2,2,3,3,4,4,5,5,2.5)
    T1 <- data.frame(matrix(
      c(
        1.0, 1.0, 2.0, 2.0, 2.5, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0,
        0.8181818, 0.8181818, 0.6363636, 0.6363636, 0.5454545, 0.3636364, 0.3636364, 0.1818182, 0.1818182,0.0000000,0.0000000
      ),
      ncol=2
    ))
    colnames(T1) <- c("x", "S(x)")
    B <- NNS.CDF(A, type="survival")
    expect_equal(as.data.frame(B$Function), T1, tolerance=1e-5)
    expect_equal(B$target.value, numeric(0), tolerance=1e-5)
  })
})

test_that("R13: multivariate NNS.CDF: degree-0 joint survival, degree > 0 survival, cumulative hazard and unsupported hazards", {
  withr::local_seed(20261108)
  cdf_moments_var <- fx_cdf_moments_var()
  aud_case("degree-0 survival for d = 3 matches strict joint counts with ties", env = cdf_moments_var, {
    set.seed(7)
    A <- matrix(round(rnorm(90), 1), ncol = 3)
    S <- NNS.CDF(A, 0, type = "survival", plot = FALSE)$Function[[4]]
    expect_equal(S, apply(A, 1, strict_joint_survival, A = A), tolerance = 1e-12)
    t <- colMeans(A)
    expect_equal(NNS.CDF(A, 0, target = t, type = "survival", plot = FALSE)$target.value,
                 strict_joint_survival(A, t), tolerance = 1e-12)
  })
  aud_case("unsupported multivariate hazard never returns a CDF divided by survival", {
    set.seed(7)
    A <- matrix(rnorm(80), ncol = 2)
    expect_error(NNS.CDF(A, 0, target = A[5, ], type = "hazard", plot = FALSE),
                 "joint density or coordinate-hazard definition")
  })
  aud_case("bivariate cumulative hazard target.value equals the row value at that target", {
    set.seed(7)
    A <- matrix(rnorm(80), ncol = 2)
    out <- NNS.CDF(A, 0, target = A[5, ], type = "cumulative hazard", plot = FALSE)
    expect_equal(out$target.value, out$Function[[3]][5], tolerance = 1e-12)
  })
  aud_case("degree > 0 multivariate survival is the normalised joint UPM", {
    set.seed(7)
    A <- matrix(rnorm(90), ncol = 3)
    S <- NNS.CDF(A, 1, type = "survival", plot = FALSE)$Function[[4]]
    expect_equal(S, apply(A, 1, function(t) Co.UPM_nD(A, t, degree = 1)), tolerance = 1e-12)
  })
  aud_case("degree > 0 survival is 0 with a coordinate at its maximum", {
    B <- cbind(c(0, 1, 2, 4), c(1, 0, 3, 2))
    expect_equal(NNS.CDF(B, 1, target = c(4, 1), type = "survival", plot = FALSE)$target.value, 0)
  })
  aud_case("degree-0 joint survival with many columns is the joint upper share", {
    set.seed(9)
    x <- matrix(sample(0:4, 400, TRUE), 20, 20)
    s <- sapply(1:20, function(i) mean(colSums(t(x) > x[i, ]) == 20))
    h <- NNS.CDF(x, degree = 0, type = "cumulative hazard", plot = FALSE)$Function$CDF
    expect_equal(h, pmax(-log(pmax(s, 1e-10)), 0))
  })
  aud_case("multivariate cumulative hazard needs degree 0 (#100)", {
    X <- cbind(c(0, 1, 2, 3), c(0, 2, 1, 3))
    expect_error(NNS.CDF(X, degree = 1, type = "hazard", plot = FALSE), "not implemented")
    expect_error(NNS.CDF(X, degree = 1, type = "cumulative hazard", plot = FALSE), "needs degree = 0")
    expect_true(is.list(NNS.CDF(X, degree = 0, type = "cumulative hazard", plot = FALSE)))
    # One dimension at degree 1 is a proper distribution function and keeps its hazards.
    set.seed(1)
    expect_true(all(is.finite(NNS.CDF(rnorm(50), degree = 1, type = "hazard", plot = FALSE)$Function[[2]])))
    # The degree-0 joint CDF gives every rectangle nonnegative mass.
    A <- matrix(rnorm(200), 100, 2)
    f <- function(t) NNS.CDF(A, degree = 0, target = t, plot = FALSE)$target.value
    expect_gte(f(c(.5, .5)) - f(c(-.5, .5)) - f(c(.5, -.5)) + f(c(-.5, -.5)), 0)
  })
})

test_that("R14: univariate NNS.CDF hazards: tie runs and constant data; NA dropped; large constant sample at degree > 0", {
  withr::local_seed(20261108)
  aud_case("hazard is finite with tie runs at the minimum, interior and maximum", {
    x <- c(rep(0, 12), 1, 2, rep(3, 12), 4, 5, rep(6, 12))
    expect_true(all(is.finite(NNS.CDF(x, type = "hazard", plot = FALSE)$Function[[2]])))
  })
  aud_case("hazard of a constant variable stops with a clear message", {
    expect_error(NNS.CDF(rep(2, 5), type = "hazard", plot = FALSE), "at least two distinct values")
  })
  aud_case("NNS.CDF drops missing observations", {
    y <- c(1, 2, NA, 4)
    expect_equal(NNS.CDF(y, plot = FALSE)$Function$CDF, c(1, 2, 3) / 3)
    expect_equal(NNS.CDF(y, target = 2, plot = FALSE)$target.value, 2 / 3)
  })
  aud_case("degree > 0 CDF of a large constant sample is 1 and not quadratic", {
    elapsed <- system.time(z <- NNS.CDF(rep(2, 1e5), 1, plot = FALSE))[["elapsed"]]
    expect_true(all(z$Function$CDF == 1))
    expect_lt(elapsed, 1)
  })
})

test_that("R15: NNS.CDF target.value equals the curve at data points, lies between neighbours and clamps the hazard", {
  withr::local_seed(20261108)
  cdf_target_value <- fx_cdf_target_value()
  aud_case("univariate target.value equals the curve at each data point", env = cdf_target_value, {
    for (x in cdf_target_data()) for (degree in 0:1) for (type in cdf_types) {
      curve <- NNS.CDF(x, degree, type = type, plot = FALSE)$Function[[2]]
      at <- vapply(sort(x), function(t) cdf_at(x, degree, type, t), numeric(1))
      expect_equal(at, curve, tolerance = 1e-12, info = paste(type, degree))
    }
  })
  aud_case("univariate target.value lies between the neighbouring curve values", env = cdf_target_value, {
    for (x in cdf_target_data()) for (degree in 0:1) for (type in cdf_types[-3]) {
      xs <- sort(unique(x))
      curve <- NNS.CDF(x, degree, type = type, plot = FALSE)$Function[[2]]
      curve <- curve[!duplicated(sort(x))]
      mid <- (xs[-1] + xs[-length(xs)]) / 2
      at <- vapply(mid, function(t) cdf_at(x, degree, type, t), numeric(1))
      lo <- pmin(curve[-1], curve[-length(curve)])
      hi <- pmax(curve[-1], curve[-length(curve)])
      expect_true(all(at >= lo - 1e-12 & at <= hi + 1e-12), info = paste(type, degree))
    }
  })
  aud_case("hazard target.value between data points is finite and clamped like the curve", env = cdf_target_value, {
    for (x in cdf_target_data()) for (degree in 0:1) {
      xs <- sort(unique(x))
      mid <- (xs[-1] + xs[-length(xs)]) / 2
      at <- vapply(mid, function(t) cdf_at(x, degree, "hazard", t), numeric(1))
      expect_true(all(is.finite(at) & at >= 0 & at <= 1e6))
    }
  })
  aud_case("cumulative hazard at a target is -log of the survival there", env = cdf_target_value, {
    expect_equal(cdf_at(1:20, 0, "cumulative hazard", 19), -log(.05))
    expect_equal(cdf_at(1:20, 0, "cumulative hazard", 10), -log(.5))
    h <- NNS.CDF(1:20, 0, target = 20, type = "cumulative hazard", plot = FALSE)
    expect_true(is.finite(h$target.value))
    expect_equal(h$target.value, h$Function[[2]][20])
  })
})
