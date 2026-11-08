# Degree handling in the partial-moment kernels.
# LPM()/UPM() with more than 32 targets use the sorted prefix backend; with 32
# or fewer they scan directly. `pad` adds 32 targets to force the prefix path.
pad <- seq(-1, 1, length.out = 32)

test_that("a degree just above 0 excludes ties on the prefix path as on the direct path", {
  x <- c(1, 2, 2, 3)
  expect_identical(LPM(1e-13, 2, x), 0.25)
  expect_identical(LPM(1e-13, c(2, pad), x)[1], LPM(1e-13, 2, x))
  expect_identical(LPM(0.1 * 3 - 0.3, c(2, pad), x)[1], LPM(0.1 * 3 - 0.3, 2, x))
})

test_that("a degree just above 0 gives the same ratios on both dispatch paths", {
  x <- c(1, 2, 2, 3)
  expect_identical(LPM.ratio(1e-13, c(2, pad), x)[1], LPM.ratio(1e-13, 2, x))
  expect_identical(UPM.ratio(1e-13, c(2, pad), x)[1], UPM.ratio(1e-13, 2, x))
})

test_that("near-integer degrees are evaluated with pow on every path, not rounded", {
  set.seed(1)
  r <- round(rnorm(250, 0, 0.01), 4)
  t <- sort(r)
  expect_identical(LPM.ratio(1e-13, t, r),
                   vapply(t, function(tt) LPM.ratio(1e-13, tt, r), numeric(1)))
  expect_identical(LPM(2 + 1e-13, t, r),
                   vapply(t, function(tt) LPM(2 + 1e-13, tt, r), numeric(1)))
})

test_that("NNS.CDF at a tiny positive degree is not the degree-0 ECDF on tied data", {
  set.seed(1)
  r <- round(rnorm(250, 0, 0.01), 4)
  expect_false(identical(NNS.CDF(r, degree = 1e-13, plot = FALSE)$Function$CDF,
                         NNS.CDF(r, degree = 0, plot = FALSE)$Function$CDF))
})

test_that("infinite degree is rejected by the univariate partial moments", {
  x <- c(0, 0.5, 2)
  expect_error(LPM(Inf, 1, x), "`degree` must be finite.", fixed = TRUE)
  expect_error(UPM(Inf, 1, x), "`degree` must be finite.", fixed = TRUE)
  expect_error(LPM(Inf, c(1, pad), x), "`degree` must be finite.", fixed = TRUE)
  expect_error(LPM(Inf, 1, x, excess_ret = TRUE), "`degree` must be finite.", fixed = TRUE)
  expect_error(LPM.ratio(Inf, 1, x), "`degree` must be finite.", fixed = TRUE)
  expect_error(UPM.ratio(Inf, 1, x), "`degree` must be finite.", fixed = TRUE)
  expect_error(NNS.CDF(x, degree = Inf, plot = FALSE), "`degree` must be finite.", fixed = TRUE)
})

test_that("infinite degree is rejected by the co- and divergent partial moments", {
  x <- c(0, 0.5, 2)
  expect_error(Co.LPM(Inf, x, x, 1, 1), "`degree_lpm` must be finite.", fixed = TRUE)
  expect_error(Co.LPM(1, x, x, 1, 1, degree_y = Inf), "`degree_y` must be finite.", fixed = TRUE)
  expect_error(Co.UPM(Inf, 3, 0, 1, 1), "`degree_upm` must be finite.", fixed = TRUE)
  expect_error(D.LPM(1, Inf, x, x, 1, 1), "`degree_upm` must be finite.", fixed = TRUE)
  expect_error(D.UPM(Inf, 1, x, x, 1, 1), "`degree_lpm` must be finite.", fixed = TRUE)
  expect_error(PM.matrix(Inf, 1, target = "mean", variable = cbind(1:3, c(3, 1, 2)), pop_adj = TRUE),
               "`LPM_degree` must be finite.", fixed = TRUE)
})

test_that("infinite degree is rejected by the n-dimensional partial moments", {
  m <- rbind(c(3, 1.5), c(0, 1), c(2, 2))
  expect_error(Co.LPM_nD(m, c(1, 1), Inf), "`degree` must be finite.", fixed = TRUE)
  expect_error(Co.UPM_nD(m, c(1, 1), Inf), "`degree` must be finite.", fixed = TRUE)
  expect_error(DPM_nD(m, c(1, 1), Inf), "`degree` must be finite.", fixed = TRUE)
  expect_error(NNS.CDF(m, degree = Inf, plot = FALSE), "`degree` must be finite.", fixed = TRUE)
})

test_that("negative and missing degrees keep the non-negative error", {
  x <- c(0, 0.5, 2)
  expect_error(LPM(-Inf, 1, x), "'degree' must be a non-negative number", fixed = TRUE)
  expect_error(LPM(NaN, 1, x), "'degree' must be a non-negative number", fixed = TRUE)
  expect_error(LPM(-1, 1, x), "'degree' must be a non-negative number", fixed = TRUE)
})

test_that("finite degrees above the int range are accepted and use pow", {
  x <- c(1, 1.5, 2)
  y <- c(1, 2, 2, 3)
  expect_identical(LPM(3e9, 3, y), Inf)
  expect_identical(LPM(2^31, 2, x), 1 / 3)
  expect_identical(LPM(2^31, 1.2, x), 0)
  expect_identical(UPM(3e9, 0, x), Inf)
  expect_identical(LPM(3e9, 3, y, excess_ret = TRUE), Inf)
  expect_identical(Co.LPM(3e9, x, x, 2, 2), 1 / 3)
})

test_that("integer degrees up to 32 keep the left-fold product bit for bit", {
  # Deliberate compatibility pin: earlier releases evaluated integer degrees as
  # the product 1 * v * v * ..., and results for degrees <= 32 must not move.
  expect_identical(LPM(32, 1, -0.3), Reduce(`*`, rep(1 - (-0.3), 32), 1))
  expect_identical(UPM(5, 0, 1.1), Reduce(`*`, rep(1.1, 5), 1))
  expect_identical(LPM(3, 1, -0.3, excess_ret = TRUE), Reduce(`*`, rep(1.3, 3), 1))
})

test_that("integer degrees above 32 agree with the power to a few ulps", {
  tol <- 8 * .Machine$double.eps
  expect_equal(LPM(33, 1, -0.3), 1.3^33, tolerance = tol)
  expect_equal(UPM(40, 0, 1.1), 1.1^40, tolerance = tol)
  expect_equal(Co.LPM(40, -0.3, -0.3, 1, 1, degree_y = 1), 1.3^40 * 1.3, tolerance = tol)
})

test_that("the cost of a huge integer degree does not grow with the degree", {
  # The left fold took about 11 s for degree 1e7 and n = 1000 and would take
  # hours for 2^31 - 1; pow takes microseconds. The bound only checks that the
  # calls do not hang, so it is loose enough for slow or instrumented builds.
  el <- system.time(r1 <- LPM(2^31 - 1, 2, rep(1, 1000)))[["elapsed"]]
  el <- el + system.time(r2 <- UPM(1e7, 0, rep(1, 1000)))[["elapsed"]]
  el <- el + system.time(r3 <- Co.LPM(2e8, c(0, 0.5, 2), c(0, 0.5, 2), 1, 1))[["elapsed"]]
  expect_lt(el, 60)
  expect_identical(c(r1, r2, r3), c(1, 1, 1 / 3))
})

test_that("the prefix path returns Inf, not NaN, when an anchored sum overflows", {
  x <- c(-1e308, 1e308)
  expect_identical(LPM(2, c(1e308, pad), x)[1], LPM(2, 1e308, x))
  expect_identical(LPM(2, c(1e308, pad), x)[1], Inf)
})

test_that("the prefix recurrence returns Inf, not NaN, when a gap power underflows", {
  w <- c(-1e200, 0, 1e-170, 1, 2)
  expect_identical(LPM(4, c(2, 1.5, pad), w)[1:2], c(LPM(4, 2, w), LPM(4, 1.5, w)))
  expect_identical(UPM(4, c(-1e-170, pad), c(0, 1e200))[1], UPM(4, -1e-170, c(0, 1e200)))
})

test_that("the prefix path has no NaN at a high degree on widely spread data", {
  set.seed(2)
  z <- runif(40) * 1e10
  out <- LPM(32, z, z)
  expect_false(any(is.nan(out)))
  expect_identical(is.infinite(out), is.infinite(vapply(z, function(t) LPM(32, t, z), numeric(1))))
})

test_that("a near-integer UPM degree is not rounded on the prefix path", {
  x <- c(0, 0.5, 1, 2)
  expect_identical(UPM(2 + 1e-13, c(0.5, pad), x)[1], UPM(2 + 1e-13, 0.5, x))
  expect_false(identical(UPM(2 + 1e-13, c(0.5, pad), x)[1], UPM(2, 0.5, x)))
})

test_that("LPM.ratio on the prefix path is 0, not NaN, when only the upper sum overflows", {
  y <- c(0, 1e200)
  expect_identical(LPM.ratio(4, c(-1e-170, pad), y)[1], LPM.ratio(4, -1e-170, y))
  expect_identical(LPM.ratio(4, c(-1e-170, pad), y)[1], 0)
})

test_that("UPM.ratio stays NaN (Inf/Inf) on the prefix path when the upper sum overflows", {
  y <- c(0, 1e200)
  expect_identical(UPM.ratio(4, c(-1e-170, pad), y)[1], UPM.ratio(4, -1e-170, y))
  expect_identical(UPM.ratio(4, c(-1e-170, pad), y)[1], NaN)
})

test_that("univariate NNS.CDF over 32 observations has no NaN when a power sum overflows", {
  d <- c(0, 1e200, seq(1, 2, length.out = 40))
  f <- NNS.CDF(d, degree = 4, plot = FALSE)$Function
  expect_identical(f$CDF[1:3], vapply(f$x[1:3], function(t) LPM.ratio(4, t, d), numeric(1)))
  expect_identical(f$CDF[1:3], c(0, 0, 0))
})
