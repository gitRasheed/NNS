# NNS.moments as partial moments about the mean. Lean: nns-proofs Moments.lean.

# Lean: UPM_one_mean_eq_LPM_one_mean, popVariance_eq_PM, moments_eq_momentsPM, moments_probe | BUGS: none
test_that("N01: NNS.moments (population) is the partial-moment assembly about the mean", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("UPM_one_mean_eq_LPM_one_mean: the mean balances degree-1 partial moments", {
    expect_identical(UPM(1, mean(fm$x_dy), fm$x_dy), LPM(1, mean(fm$x_dy), fm$x_dy))
    expect_identical(UPM(1, mean(fm$x_int), fm$x_int), LPM(1, mean(fm$x_int), fm$x_int))
  })
  aud_case("popVariance_eq_PM: population variance is LPM(2, mu) + UPM(2, mu), exact on dyadic data", {
    m <- NNS.moments(fm$x_dy)
    expect_identical(m$mean, 0)
    expect_identical(m$variance, LPM(2, 0, fm$x_dy) + UPM(2, 0, fm$x_dy))
    m <- NNS.moments(fm$x_int)
    expect_identical(m$variance, LPM(2, 3, fm$x_int) + UPM(2, 3, fm$x_int))
  })
  aud_case("moments_eq_momentsPM (population): NNS.moments equals the partial-moment assembly", {
    aud_moments_equal(NNS.moments(fm$x_dy), aud_moments_pm(fm$x_dy), aud_moments_tol(fm$x_dy), "dyadic")
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      aud_moments_equal(NNS.moments(x), aud_moments_pm(x), aud_moments_tol(x), nm)
    }
  })
  aud_case("moments_probe: NNS.moments(c(2,4,4,7,9)) matches the exact rational values", {
    m <- NNS.moments(c(2, 4, 4, 7, 9))
    expect_identical(m$mean, 5.2)
    expect_equal(m$variance, 154 / 25, tolerance = 4 * aud_eps)
    expect_equal(m$skewness, 306 / (77 * sqrt(154)), tolerance = 8 * aud_eps)
    expect_equal(m$kurtosis, -15073 / 11858, tolerance = 8 * aud_eps)
  })
})

# Lean: moments_eq_momentsPM, moments_probe_sample, sampleVariance_eq, sampleSkewness_eq_G1, sampleKurtosis_eq_G2, popVariance_eq_zero_iff, centralMoment_three_mean_of_card_two | BUGS: none
test_that("N02: NNS.moments (sample) uses the Joanes-Gill corrections and is NaN where they are undefined", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("moments_eq_momentsPM (sample): NNS.moments(population = FALSE) equals the PM assembly", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      aud_moments_equal(NNS.moments(x, population = FALSE), aud_moments_pm(x, FALSE),
                        2 * aud_moments_tol(x), nm)
    }
  })
  aud_case("moments_probe_sample: NNS.moments(c(2,4,4,7,9), FALSE) matches the exact values", {
    m <- NNS.moments(c(2, 4, 4, 7, 9), population = FALSE)
    expect_identical(m$mean, 5.2)
    expect_equal(m$variance, 77 / 10, tolerance = 4 * aud_eps)
    expect_equal(m$skewness, sqrt(20) / 3 * 306 / (77 * sqrt(154)), tolerance = 8 * aud_eps)
    expect_equal(m$kurtosis, -6430 / 5929, tolerance = 8 * aud_eps)
  })
  aud_case("sampleVariance_eq: sample variance is sum((x - mu)^2) / (n - 1)", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      aud_within_max(NNS.moments(x, FALSE)$variance, sum((x - mean(x))^2) / (length(x) - 1),
                     aud_mom_tol(x, 2), label = nm)
    }
  })
  aud_case("sampleSkewness_eq_G1 / sampleKurtosis_eq_G2: Joanes-Gill G1 and G2 (n > 3)", {
    g1g2 <- function(x) {
      n <- length(x); z <- x - mean(x)
      m2 <- mean(z^2); m3 <- mean(z^3); m4 <- mean(z^4)
      g1 <- m3 / m2^1.5; g2 <- m4 / m2^2 - 3
      c(skewness = sqrt(n * (n - 1)) / (n - 2) * g1,
        kurtosis = (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6))
    }
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      m <- NNS.moments(x, FALSE)
      tol <- 2 * aud_moments_tol(x)
      aud_within_max(m$skewness, g1g2(x)[["skewness"]], tol[["skewness"]], label = nm)
      aud_within_max(m$kurtosis, g1g2(x)[["kurtosis"]], tol[["kurtosis"]], label = nm)
    }
  })
  aud_case("popVariance_eq_zero_iff: a constant sample has variance exactly 0 (skew, kurt NaN by R convention)", {
    for (v in c(3, -0.1, 1e6)) {
      m <- NNS.moments(rep(v, 7))
      expect_identical(m$variance, 0)
      expect_true(is.nan(m$skewness))
      expect_true(is.nan(m$kurtosis))
    }
    expect_gt(NNS.moments(c(3, 3, 3, 3.5))$variance, 0)
  })
  aud_case("popVariance_eq_zero_iff (n = 1): population variance 0, sample variance NaN", {
    m <- NNS.moments(5)
    expect_identical(m$mean, 5)
    expect_identical(m$variance, 0)
    expect_true(is.nan(NNS.moments(5, FALSE)$variance))
  })
  aud_case("centralMoment_three_mean_of_card_two: n = 2 has m3 = 0, so sample skewness is NaN", {
    expect_identical(UPM(3, 1.5, c(1, 2)) - LPM(3, 1.5, c(1, 2)), 0)
    expect_identical(NNS.moments(c(1, 2))$skewness, 0)
    expect_true(is.nan(NNS.moments(c(1, 2), FALSE)$skewness))
  })
  aud_case("sampleKurtosis_eq_G2 guard (n > 3): n = 3 sample kurtosis is NaN", {
    m <- NNS.moments(c(1, 2, 3), FALSE)
    expect_identical(m$skewness, 0)
    expect_true(is.nan(m$kurtosis))
  })
  aud_case("NNS.moments(c(1, 2, 4, 8, 16), population = FALSE) skewness 1.325315", {
    expect_equal(NNS.moments(c(1, 2, 4, 8, 16), population = FALSE)$skewness, 1.325315, tolerance = 1e-6)
  })
})

# Lean: popSkewness_pos_iff, popSkewness_neg_iff, popSkewness_sq_sub_two_le_popKurtosis, popKurtosis_le, popSkewness_sq_le, popVariance_affine, popSkewness_affine, popKurtosis_affine, sampleSkewness_affine, sampleKurtosis_affine | BUGS: none
test_that("N03: standardized moments: sign of UPM3 - LPM3, Pearson's inequality, bounds and affine invariance", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("popSkewness_pos_iff / popSkewness_neg_iff: sign(skewness) = sign(UPM(3, mu) - LPM(3, mu))", {
    samples <- c(fm$all_samples, list(right = c(0, 0, 0, 0, 0, 0, 0, 8), left = -c(0, 0, 0, 0, 0, 0, 0, 8)))
    for (nm in names(samples)) {
      x <- samples[[nm]]
      expect_identical(sign(NNS.moments(x)$skewness), sign(UPM(3, mean(x), x) - LPM(3, mean(x), x)),
                       label = nm)
    }
  })
  aud_case("popSkewness_sq_sub_two_le_popKurtosis: Pearson's inequality kurtosis >= skewness^2 - 2", {
    samples <- c(fm$all_samples, list(two_point = c(0, 0, 0, 1), bernoulli = rep(c(0, 1), 8)))
    for (nm in names(samples)) {
      m <- NNS.moments(samples[[nm]])
      expect_gte(m$kurtosis - (m$skewness^2 - 2), -64 * aud_eps * (1 + m$skewness^2), label = nm)
    }
    # equality case: two-point distributions attain the bound
    m <- NNS.moments(c(0, 0, 0, 1))
    expect_equal(m$kurtosis, m$skewness^2 - 2, tolerance = 16 * aud_eps)
  })
  aud_case("popKurtosis_le / popSkewness_sq_le: kurtosis <= n - 3 and skewness^2 <= n - 1", {
    samples <- c(fm$all_samples, list(spike = c(rep(0, 15), 1)))
    for (nm in names(samples)) {
      x <- samples[[nm]]; n <- length(x)
      m <- NNS.moments(x)
      expect_lte(m$kurtosis, (n - 3) * (1 + 64 * aud_eps), label = nm)
      expect_lte(m$skewness^2, (n - 1) * (1 + 64 * aud_eps), label = nm)
    }
    # a single spike of n = 16: skewness^2 = 196/15, kurtosis = 166/15 (exact values)
    m <- NNS.moments(c(rep(0, 15), 1))
    expect_equal(m$skewness^2, 196 / 15, tolerance = 16 * aud_eps)
    expect_equal(m$kurtosis, 166 / 15, tolerance = 16 * aud_eps)
  })
  aud_case("popVariance_affine / popSkewness_affine / popKurtosis_affine: exact dyadic map x -> -2 x + 8", {
    m <- NNS.moments(fm$x_dy)
    ma <- NNS.moments(-2 * fm$x_dy + 8)
    expect_identical(ma$mean, -2 * m$mean + 8)
    expect_identical(ma$variance, 4 * m$variance)
    expect_identical(ma$skewness, -m$skewness)
    expect_identical(ma$kurtosis, m$kurtosis)
  })
  aud_case("popVariance_affine / popSkewness_affine / popKurtosis_affine: random and real data", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      for (a in c(-3, 0.5, 7)) {
        y <- a * x + 2.5
        m <- NNS.moments(x); ma <- NNS.moments(y)
        tol <- aud_moments_tol(y)
        aud_within_max(ma$variance, a^2 * m$variance, tol[["variance"]] + a^2 * aud_moments_tol(x)[["variance"]], label = nm)
        aud_within_max(ma$skewness, sign(a) * m$skewness, 2 * tol[["skewness"]], label = nm)
        aud_within_max(ma$kurtosis, m$kurtosis, 2 * tol[["kurtosis"]], label = nm)
      }
    }
  })
  aud_case("sampleSkewness_affine / sampleKurtosis_affine: sample branch under x -> a x + c", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      m <- NNS.moments(x, FALSE); ma <- NNS.moments(-4 * x + 1, FALSE)
      tol <- 2 * aud_moments_tol(-4 * x + 1)
      aud_within_max(ma$skewness, -m$skewness, 2 * tol[["skewness"]], label = nm)
      aud_within_max(ma$kurtosis, m$kurtosis, 2 * tol[["kurtosis"]], label = nm)
    }
  })
})

# Lean: popSkewness_affine, popKurtosis_affine | BUGS: #12
test_that("N04: NNS.moments stays finite at data scales 1e-200 and 1e200 (BUGS #12 family)", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("NNS.moments sample skewness and kurtosis are finite on data of scale 1e200", {
    x <- c(-3, 1, 2, 2, 5)
    big <- NNS.moments(1e200 * x, population = FALSE)
    ref <- NNS.moments(x, population = FALSE)
    expect_equal(big$skewness, ref$skewness, tolerance = 1e-14)
    expect_equal(big$kurtosis, ref$kurtosis, tolerance = 1e-14)
  })
  aud_case("NNS.moments keeps the underflowed variance of data of scale 1e-200", {
    expect_identical(NNS.moments(c(-3e-200, 1e-200, 2e-200, 2e-200))$variance, 0)
  })
  aud_case("popSkewness_affine / popKurtosis_affine: invariance under scale a = 1e-200", {
    # fails on 8898017: see BUGS.md #12 (same underflow family; the NNS.moments
    # case is the JOURNAL.md "scale 1e-200" finding, no dedicated BUGS.md id):
    # mean(z^3), mean(z^4) underflow to 0 and skew/kurt become NaN.
    for (nm in c("normal", "eruptions", "mpg")) {
      x <- fm$all_samples[[nm]]
      m <- NNS.moments(x); ms <- NNS.moments(1e-200 * x)
      expect_equal(ms$skewness, m$skewness, tolerance = 1e-12, label = nm)
      expect_equal(ms$kurtosis, m$kurtosis, tolerance = 1e-12, label = nm)
    }
  })
})
