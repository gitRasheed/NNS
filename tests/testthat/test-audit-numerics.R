# Numerical robustness of the partial-moment ratios, the VaR solvers and
# NNS.moments at extreme degrees and scales, and the zero-total convention of
# the normalised nD partial moments.

# --- LPM.ratio / UPM.ratio: scale-free at extreme degree --------------------

test_that("LPM.ratio is 0.5, not NaN, when every deviation underflows", {
  expect_identical(LPM.ratio(400, 0, c(-1e-2, 1e-2)), 0.5)
})

test_that("UPM.ratio is 0.5, not NaN, when every deviation overflows", {
  expect_identical(UPM.ratio(400, 0, c(-1e3, 1e3)), 0.5)
})

test_that("LPM.ratio does not depend on the units of the data at high degree", {
  x <- c(-2, -1, 0.5, 3)
  expect_equal(LPM.ratio(150, 0, 1e-3 * x), LPM.ratio(150, 0, x), tolerance = 1e-12)
  expect_equal(LPM.ratio(150, 0, 1e3 * x), LPM.ratio(150, 0, x), tolerance = 1e-12)
})

test_that("LPM.ratio and UPM.ratio on more than 32 targets have no NaN at degree 300", {
  set.seed(1)
  r <- rnorm(250, 0, 0.01)
  expect_false(anyNA(LPM.ratio(300, r, r)))
  expect_equal(LPM.ratio(300, r, r) + UPM.ratio(300, r, r), rep(1, 250))
})

test_that("univariate NNS.CDF at degree 400 has no NaN on daily-return data", {
  set.seed(1)
  r <- rnorm(250, 0, 0.01)
  expect_false(anyNA(NNS.CDF(r, degree = 400, plot = FALSE)$Function$CDF))
})

test_that("LPM.ratio on ordinary data is still the plain L / (L + U)", {
  set.seed(2)
  r <- rnorm(1000, 0, 0.02)
  expect_identical(LPM.ratio(2, 0, r), LPM(2, 0, r) / (LPM(2, 0, r) + UPM(2, 0, r)))
})

# --- VaR solvers --------------------------------------------------------------

test_that("the non-integer-degree VaR solver finds the symmetric median at degree 3e9 without warnings", {
  expect_no_warning(v <- LPM.VaR(0.5, 3e9, c(0, 1)))
  expect_equal(v, 0.5, tolerance = 1e-6)
})

test_that("UPM.VaR at degree 1e4 finds the symmetric median", {
  expect_equal(UPM.VaR(0.5, 1e4, c(0, 1)), 0.5, tolerance = 1e-6)
})

test_that("the integer-degree VaR solver is scale invariant on data of scale 1e-200", {
  x <- c(-3, 1, 2, 2)
  p <- c(0.1, 0.5, 0.9)
  expect_equal(LPM.VaR(p, 2, 1e-200 * x) * 1e200, LPM.VaR(p, 2, x), tolerance = 1e-12)
  expect_equal(UPM.VaR(p, 4, 1e200 * x) * 1e-200, UPM.VaR(p, 4, x), tolerance = 1e-12)
})

test_that("the non-integer-degree VaR solver is scale invariant on data of scale 1e-200", {
  x <- c(-3, 1, 2, 2)
  p <- c(0.1, 0.5, 0.9)
  expect_equal(LPM.VaR(p, 2.5, 1e-200 * x) * 1e200, LPM.VaR(p, 2.5, x), tolerance = 1e-6)
})

test_that("integer-degree VaR meets its absolute ratio tolerance at tiny percentiles", {
  # The documented contract is |ratio(VaR) - p| <= about eps^0.75, not a
  # relative error: at p near 1e-21 the returned threshold can be far from the
  # order statistic that has that ratio (see LPM.VaR's value section).
  set.seed(11)
  x <- c(rnorm(40, -1e5, 1), rnorm(61, 0, 1))
  p <- LPM.ratio(4, sort(x)[c(5, 20, 35)], x)
  expect_true(all(abs(LPM.ratio(4, LPM.VaR(p, 4, x), x) - p) <= .Machine$double.eps^0.75))
})

# --- NNS.moments --------------------------------------------------------------

test_that("NNS.moments skewness and kurtosis are finite on data of scale 1e-200", {
  x <- c(-3, 1, 2, 2)
  tiny <- NNS.moments(1e-200 * x)
  ref <- NNS.moments(x)
  expect_equal(tiny$skewness, ref$skewness, tolerance = 1e-14)
  expect_equal(tiny$kurtosis, ref$kurtosis, tolerance = 1e-14)
})

test_that("NNS.moments sample skewness and kurtosis are finite on data of scale 1e200", {
  x <- c(-3, 1, 2, 2, 5)
  big <- NNS.moments(1e200 * x, population = FALSE)
  ref <- NNS.moments(x, population = FALSE)
  expect_equal(big$skewness, ref$skewness, tolerance = 1e-14)
  expect_equal(big$kurtosis, ref$kurtosis, tolerance = 1e-14)
})

test_that("NNS.moments keeps the underflowed variance of data of scale 1e-200", {
  expect_identical(NNS.moments(c(-3e-200, 1e-200, 2e-200, 2e-200))$variance, 0)
})

# --- normalised nD partial moments: zero total ----------------------------------

test_that("normalised nD partial moments with every row tied take the degree-0 values", {
  X <- matrix(1, 3, 2)
  expect_identical(c(Co.LPM_nD(X, c(1, 1), 1), Co.UPM_nD(X, c(1, 1), 1), DPM_nD(X, c(1, 1), 1)),
                   c(1, 0, 0))
})

test_that("normalised nD partial moments with a tied column sum to 1", {
  X <- cbind(1:5, rep(1, 5))
  cells <- c(Co.LPM_nD(X, c(3, 1), 2), Co.UPM_nD(X, c(3, 1), 2), DPM_nD(X, c(3, 1), 2))
  expect_identical(cells, c(Co.LPM_nD(X, c(3, 1), 0), Co.UPM_nD(X, c(3, 1), 0), DPM_nD(X, c(3, 1), 0)))
  expect_equal(sum(cells), 1)
})

test_that("one-column normalised Co.LPM_nD on constant data agrees with LPM.ratio", {
  expect_identical(Co.LPM_nD(matrix(1, 3, 1), 1, 1), LPM.ratio(1, 1, rep(1, 3)))
  expect_identical(Co.UPM_nD(matrix(1, 3, 1), 1, 1), UPM.ratio(1, 1, rep(1, 3)))
})

test_that("the batched nD kernel takes the degree-0 value on a zero total", {
  set.seed(1)
  Z <- cbind(rnorm(6), 0)
  expect_identical(NNS:::Co.LPM_nD.batch(Z, Z, 1, TRUE), NNS:::Co.LPM_nD.batch(Z, Z, 0, TRUE))
})

test_that("multivariate NNS.CDF at degree 1 with a constant column equals the degree-0 CDF", {
  X <- cbind(1:5, rep(1, 5))
  expect_identical(NNS.CDF(X, 1, plot = FALSE)$Function$CDF, NNS.CDF(X, 0, plot = FALSE)$Function$CDF)
})

test_that("multivariate NNS.CDF target.value agrees with the curve on a constant column", {
  X <- cbind(1:5, rep(1, 5))
  expect_identical(NNS.CDF(X, 1, target = c(3, 1), plot = FALSE)$target.value,
                   NNS.CDF(X, 1, plot = FALSE)$Function$CDF[3])
})

test_that("multivariate NNS.CDF survival at degree 1 is Co.UPM_nD on a constant column", {
  X <- cbind(1:5, rep(1, 5))
  S <- NNS.CDF(X, 1, type = "survival", plot = FALSE)$Function$CDF
  expect_identical(S, apply(X, 1, function(t) Co.UPM_nD(X, t, 1)))
  expect_identical(S, rep(0, 5))
})

test_that("the batched upper kernel is the reflected lower kernel on ordinary data", {
  set.seed(3)
  Y <- matrix(rnorm(60), 20)
  expect_identical(NNS:::Co.LPM_nD.batch(Y, Y, 2, TRUE, upper = TRUE),
                   NNS:::Co.LPM_nD.batch(-Y, -Y, 2, TRUE))
  expect_equal(NNS:::Co.LPM_nD.batch(Y, Y, 2, TRUE, upper = TRUE),
               apply(Y, 1, function(t) Co.UPM_nD(Y, t, 2)))
})

# --- normalised nD partial moments: under/overflow ----------------------------

test_that("normalised Co.LPM_nD is not 0 when every product underflows", {
  y <- matrix(c(-1, 1, -2, 2), 2)
  expect_equal(Co.LPM_nD(1e-100 * y, c(0, 0), 4), Co.LPM_nD(y, c(0, 0), 4))
  expect_equal(Co.LPM_nD(1e-100 * y, c(0, 0), 4), 0.5)
})

test_that("normalised nD partial moments are invariant to per-column scale at degree 200", {
  set.seed(1)
  M <- cbind(rnorm(250, 0, 0.01), rnorm(250, 0, 0.01))
  cells <- function(D) c(Co.LPM_nD(D, c(0, 0), 200), Co.UPM_nD(D, c(0, 0), 200), DPM_nD(D, c(0, 0), 200))
  expect_equal(cells(M %*% diag(c(1e3, 1e-3))), cells(M), tolerance = 1e-10)
  expect_equal(sum(cells(M)), 1)
})

test_that("an underflowed nD total with one tied row uses the untied rows, not the degree-0 counts", {
  Z <- rbind(c(1, 0), c(-1e-200, -1e-200), c(2e-200, 1e-200))
  expect_equal(c(Co.LPM_nD(Z, c(0, 0), 2), Co.UPM_nD(Z, c(0, 0), 2), DPM_nD(Z, c(0, 0), 2)),
               c(0.2, 0.8, 0))
})
