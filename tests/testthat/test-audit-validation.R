# Input validation and coercion across the partial-moment family.
x <- c(1, 2, 2, 3)
x5 <- c(1, 2, 2, 3, 10)

# --- excess_ret flag (#3) ----------------------------------------------------

test_that("excess_ret = NA is an error for a scalar target", {
  expect_error(LPM(1, 2, x, excess_ret = NA), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
  expect_error(UPM(1, 2, x, excess_ret = NA), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
})

test_that("excess_ret = NA gives the same error for empty and vector targets", {
  expect_error(LPM(1, numeric(0), x, excess_ret = NA), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
  expect_error(LPM(1, c(1, 2), x, excess_ret = NA), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
  expect_error(UPM(1, c(1, 2), x, excess_ret = NA), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
})

test_that("excess_ret rejects vectors and non-logical types", {
  expect_error(LPM(1, 2, x, excess_ret = c(TRUE, FALSE)), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
  expect_error(UPM(1, 2, x, excess_ret = "yes"), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
  expect_error(LPM(1, 2, x, excess_ret = NULL), "`excess_ret` must be TRUE or FALSE.", fixed = TRUE)
})

test_that("excess_ret = 1 still means TRUE", {
  expect_identical(LPM(1, 2, x, excess_ret = 1), LPM(1, 2, x, excess_ret = TRUE))
  expect_identical(UPM(1, 2, x, excess_ret = 0), UPM(1, 2, x))
})

# --- NULL / empty target in LPM/UPM (#4) -------------------------------------

test_that("LPM/UPM with a NULL target use mean(variable), like LPM.ratio", {
  expect_identical(LPM(1, NULL, x), LPM(1, mean(x), x))
  expect_identical(UPM(2, NULL, x5), UPM(2, mean(x5), x5))
  expect_identical(LPM(1, NULL, x, excess_ret = TRUE), LPM(1, mean(x), x, excess_ret = TRUE))
})

test_that("LPM with a NULL target and empty variable reports a non-finite target", {
  expect_error(LPM(1, NULL, numeric(0)), "`target` must be finite.", fixed = TRUE)
})

test_that("a zero-length target still returns numeric(0)", {
  expect_identical(LPM(1, numeric(0), x), numeric(0))
  expect_identical(UPM(1, numeric(0), x), numeric(0))
})

# --- factor input (#6, #26) --------------------------------------------------

test_that("LPM/UPM reject a factor variable instead of using its level codes", {
  f <- factor(c("10", "20", "30"))
  expect_error(LPM(1, 15, f), "`variable` must be numeric, not a factor.", fixed = TRUE)
  expect_error(UPM(1, 15, f), "`variable` must be numeric, not a factor.", fixed = TRUE)
})

test_that("LPM/UPM reject a factor target", {
  expect_error(LPM(1, factor("2"), c(1, 2, 3)), "`target` must be numeric, not a factor.", fixed = TRUE)
  expect_error(UPM(1, factor("2"), c(1, 2, 3)), "`target` must be numeric, not a factor.", fixed = TRUE)
})

test_that("ratios and co-moments reject factor data", {
  f <- factor(c("10", "20", "30"))
  expect_error(LPM.ratio(1, 15, f), "`variable` must be numeric, not a factor.", fixed = TRUE)
  expect_error(UPM.ratio(1, 15, f), "`variable` must be numeric, not a factor.", fixed = TRUE)
  expect_error(Co.LPM(1, f, c(1, 2, 3), 1, 1), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(D.UPM(1, 1, c(1, 2, 3), f, 1, 1), "`y` must be numeric, not a factor.", fixed = TRUE)
})

test_that("a data.frame holding a factor column is rejected like a bare factor", {
  expect_error(LPM.ratio(1, 15, data.frame(a = factor(c(10, 20, 30)))),
               "`variable` must be numeric, not a factor.", fixed = TRUE)
})

test_that("co-moments reject a factor target", {
  expect_error(Co.LPM(1, x5, x5, factor("z"), 1), "`target_x` must be numeric, not a factor.", fixed = TRUE)
})

test_that("NNS.CDF rejects a factor variable", {
  expect_error(NNS.CDF(factor(c(10, 20, 30)), plot = FALSE),
               "`variable` must be numeric, not a factor.", fixed = TRUE)
})

test_that("NNS.moments rejects a factor instead of using its level codes", {
  expect_error(NNS.moments(factor(c(10, 20, 40))), "`x` must be numeric, not a factor.", fixed = TRUE)
})

test_that("stochastic dominance functions reject factor input", {
  f <- factor(c(10, 20))
  expect_error(NNS.FSD.uni(f, c(1.5, 1.5)), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.SSD.uni(c(1.5, 1.5), f), "`y` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.TSD.uni(f, c(1.5, 1.5)), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(NNS.FSD(f, c(1.5, 1.5), plot = FALSE), "`x` must be numeric, not a factor.", fixed = TRUE)
})

test_that("stochastic dominance rejects a factor column inside a data.frame", {
  expect_error(NNS.SSD(data.frame(a = factor(c(10, 20))), c(1.5, 1.5), plot = FALSE),
               "`x` must be numeric, not a factor.", fixed = TRUE)
})

test_that("stochastic dominance still flattens single-column numeric containers", {
  expect_identical(NNS.FSD.uni(data.frame(a = c(10, 20)), c(1.5, 1.5)),
                   NNS.FSD.uni(c(10, 20), c(1.5, 1.5)))
  expect_identical(NNS.SSD.uni(matrix(c(10, 20), ncol = 1), c(1.5, 1.5)), NNS.SSD.uni(c(10, 20), c(1.5, 1.5)))
})

test_that("LPM keeps accepting logical and integer data", {
  expect_identical(LPM(1, 2, c(TRUE, FALSE)), 1.5)
  expect_identical(LPM(1, 2, 1:3), LPM(1, 2, c(1, 2, 3)))
})

# --- non-numeric targets (#9, #16) -------------------------------------------

test_that("a logical target is coerced like LPM(), not replaced by the mean", {
  expect_identical(LPM.ratio(0, TRUE, x5), LPM.ratio(0, 1, x5))
  expect_identical(UPM.ratio(0, c(TRUE, FALSE), x5), UPM.ratio(0, c(1, 0), x5))
  expect_identical(Co.LPM(0, x5, x5, TRUE, TRUE), Co.LPM(0, x5, x5, 1, 1))
  expect_identical(LPM.ratio(0, TRUE, x5), LPM(0, TRUE, x5))
})

test_that("a logical NA target is non-finite, as NA_real_ is", {
  expect_error(LPM.ratio(0, NA, x5), "`target` must be finite.", fixed = TRUE)
  expect_error(UPM.ratio(0, NA, x5), "`target` must be finite.", fixed = TRUE)
  expect_error(Co.LPM(1, x5, x5, NA, 1), "`target_x` must be finite.", fixed = TRUE)
  expect_error(D.LPM(1, 1, x5, x5, 1, NA), "`target_y` must be finite.", fixed = TRUE)
})

test_that("character and list targets are errors in the ratio functions", {
  expect_error(LPM.ratio(0, "a", x5), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_error(UPM.ratio(0, "1", x5), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_error(LPM.ratio(0, list(1), x5), "`target` must be numeric or NULL.", fixed = TRUE)
})

test_that("character and list targets are errors in the co- and divergent moments", {
  expect_error(Co.LPM(1, x5, x5, "3", 4), "`target_x` must be numeric or NULL.", fixed = TRUE)
  expect_error(Co.UPM(1, x5, x5, 1, "3"), "`target_y` must be numeric or NULL.", fixed = TRUE)
  expect_error(D.LPM(1, 1, x5, x5, list(0), 1), "`target_x` must be numeric or NULL.", fixed = TRUE)
  expect_error(D.UPM(1, 1, x5, x5, 1, list(0)), "`target_y` must be numeric or NULL.", fixed = TRUE)
})

test_that("LPM/UPM reject character and list targets like the ratio functions", {
  expect_error(LPM(1, "1", x5), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_error(UPM(1, list(1), x5), "`target` must be numeric or NULL.", fixed = TRUE)
})

test_that("NULL and integer targets are unchanged", {
  expect_identical(LPM.ratio(0, NULL, x5), LPM.ratio(0, mean(x5), x5))
  expect_identical(LPM.ratio(0, 2L, x5), LPM.ratio(0, 2, x5))
  expect_identical(Co.LPM(1, x5, rev(x5), NULL, NULL), Co.LPM(1, x5, rev(x5), mean(x5), mean(rev(x5))))
  expect_identical(Co.LPM(1, x5, rev(x5), 2L, 2L), Co.LPM(1, x5, rev(x5), 2, 2))
})

test_that("PM.matrix rejects non-numeric targets", {
  A <- cbind(a = x5, b = rev(x5))
  expect_error(PM.matrix(1, 1, c("0", "0"), A, TRUE), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_error(PM.matrix(1, 1, data.frame(0, 0), A, TRUE), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_error(PM.matrix(1, 1, "median", A, TRUE), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_error(PM.matrix(1, 1, c(NA, NA), A, TRUE), "`target` must be finite.", fixed = TRUE)
})

test_that("PM.matrix keeps target = \"mean\"/\"means\" as an alias for NULL", {
  A <- cbind(a = x5, b = rev(x5))
  expect_identical(PM.matrix(1, 1, "mean", A, TRUE), PM.matrix(1, 1, NULL, A, TRUE))
  expect_identical(PM.matrix(1, 1, "means", A, TRUE), PM.matrix(1, 1, NULL, A, TRUE))
})

test_that("PM.matrix integer targets are unchanged", {
  A <- cbind(a = x5, b = rev(x5))
  expect_identical(PM.matrix(1, 1, c(2L, 3L), A, FALSE), PM.matrix(1, 1, c(2, 3), A, FALSE))
})

test_that("nD functions reject a factor target", {
  m <- cbind(x5, rev(x5))
  expect_error(Co.LPM_nD(m, factor(c(1, 2)), 1), "`target` must be numeric, not a factor.", fixed = TRUE)
  expect_error(DPM_nD(m, c("1", "2"), 1), "`target` must be numeric.", fixed = TRUE)
})

test_that("NNS.CDF validates a non-numeric or missing target", {
  expect_error(NNS.CDF(x5, target = NA, plot = FALSE), "`target` must be finite.", fixed = TRUE)
  expect_error(NNS.CDF(x5, target = "2", plot = FALSE), "`target` must be numeric or NULL.", fixed = TRUE)
  expect_identical(NNS.CDF(x5, target = TRUE, plot = FALSE)$target.value, LPM.ratio(0, 1, x5))
})

# --- VaR with infinite data (#10) --------------------------------------------

test_that("LPM.VaR rejects +Inf the same way at every degree", {
  y <- c(1, 2, Inf)
  expect_error(LPM.VaR(0.5, 0, y), "`x` must be finite.", fixed = TRUE)
  expect_error(LPM.VaR(0.5, 1, y), "`x` must be finite.", fixed = TRUE)
  expect_error(LPM.VaR(0.5, 2, y), "`x` must be finite.", fixed = TRUE)
  expect_error(LPM.VaR(0.5, 1.5, y), "`x` must be finite.", fixed = TRUE)
  expect_error(LPM.VaR(0.5, 5, y), "`x` must be finite.", fixed = TRUE)
})

test_that("UPM.VaR rejects -Inf the same way at every degree", {
  y <- c(-Inf, 1, 2)
  expect_error(UPM.VaR(0.5, 0, y), "`x` must be finite.", fixed = TRUE)
  expect_error(UPM.VaR(0.5, 1, y), "`x` must be finite.", fixed = TRUE)
  expect_error(UPM.VaR(0.5, 2.5, y), "`x` must be finite.", fixed = TRUE)
})

test_that("VaR still drops NA and NaN", {
  set.seed(7)
  r <- rnorm(100, 0, 0.01)
  expect_identical(LPM.VaR(0.05, 1, c(r, NA, NaN)), LPM.VaR(0.05, 1, r))
  expect_identical(UPM.VaR(0.05, 2, c(NA, r)), UPM.VaR(0.05, 2, r))
})

test_that("VaR rejects factor and character x", {
  expect_error(LPM.VaR(0.5, 1, factor(c(1, 2, 3))), "`x` must be numeric, not a factor.", fixed = TRUE)
  expect_error(UPM.VaR(0.5, 1, c("1", "2")), "`x` must be numeric.", fixed = TRUE)
})

# --- NNS.CDF input forms (#13) -----------------------------------------------

test_that("NNS.CDF accepts a one-column data.frame like a vector", {
  expect_identical(NNS.CDF(data.frame(a = x5), plot = FALSE), NNS.CDF(x5, plot = FALSE))
  expect_identical(NNS.CDF(data.frame(a = c(x5, NA)), target = 2, plot = FALSE),
                   NNS.CDF(x5, target = 2, plot = FALSE))
})

test_that("NNS.CDF evaluates a univariate vector target element-wise", {
  expect_identical(NNS.CDF(x5, target = c(2, 3), plot = FALSE)$target.value,
                   c(NNS.CDF(x5, target = 2, plot = FALSE)$target.value,
                     NNS.CDF(x5, target = 3, plot = FALSE)$target.value))
  expect_identical(NNS.CDF(x5, target = c(2, 3), type = "survival", plot = FALSE)$target.value,
                   1 - LPM.ratio(0, c(2, 3), x5))
})

test_that("NNS.CDF hazard with a vector target matches the scalar calls", {
  expect_equal(NNS.CDF(x5, target = c(2, 3), type = "hazard", plot = FALSE)$target.value,
               c(NNS.CDF(x5, target = 2, type = "hazard", plot = FALSE)$target.value,
                 NNS.CDF(x5, target = 3, type = "hazard", plot = FALSE)$target.value))
})

test_that("NNS.CDF still rejects a vector target with an out-of-support value", {
  expect_error(NNS.CDF(x5, target = c(2, 30), plot = FALSE), "target out of bounds", fixed = TRUE)
})

# --- zero-row data (#20, #21) ------------------------------------------------

test_that("PM.matrix rejects a zero-row variable whatever the target", {
  A0 <- matrix(numeric(0), 0, 3)
  expect_error(PM.matrix(1, 1, c(0, 0, 0), A0, TRUE), "`variable` must have at least one row", fixed = TRUE)
  expect_error(PM.matrix(1, 1, NULL, A0, TRUE), "`variable` must have at least one row", fixed = TRUE)
  expect_error(PM.matrix(0, 0, c(0, 0, 0), A0, FALSE, norm = TRUE), "`variable` must have at least one row", fixed = TRUE)
})

test_that("PM.matrix with one row still returns the documented list", {
  expect_named(PM.matrix(1, 1, c(0, 0, 0), matrix(c(1, -1, 2), 1), TRUE),
               c("cupm", "dupm", "dlpm", "clpm", "cov.matrix"))
})

test_that("NNS.copula on zero-row data still short-circuits to 0", {
  expect_identical(NNS.copula(matrix(numeric(0), 0, 3)), 0)
})

test_that("nD co-partial moments reject zero-row data like the batch backend", {
  m <- matrix(numeric(0), 0, 2)
  expect_error(Co.LPM_nD(m, c(0, 0), 1, TRUE), "`data` must have at least one row", fixed = TRUE)
  expect_error(Co.LPM_nD(m, c(0, 0), 0), "`data` must have at least one row", fixed = TRUE)
  expect_error(Co.UPM_nD(m, c(0, 0), 1, FALSE), "`data` must have at least one row", fixed = TRUE)
  expect_error(DPM_nD(m, c(0, 0), 1, TRUE), "`data` must have at least one row", fixed = TRUE)
  expect_error(NNS:::Co.LPM_nD.batch(m, matrix(0, 1, 2), 1), "`data` must have at least one row", fixed = TRUE)
})

test_that("one-row nD data agrees with the batch backend", {
  r1 <- matrix(c(-1, -1), 1, 2)
  expect_identical(Co.LPM_nD(r1, c(0, 0), 1, TRUE), NNS:::Co.LPM_nD.batch(r1, matrix(0, 1, 2), 1, TRUE))
})

# --- long vectors (#7) -------------------------------------------------------

test_that("LPM/UPM excess_ret handle vectors longer than .Machine$integer.max", {
  skip_on_cran()
  skip_if_not(isTRUE(as.numeric(Sys.getenv("NNS_BIGMEM_GB", "0")) >= 64),
              "needs about 48 GiB; set NNS_BIGMEM_GB to run")
  n <- 2^32 + 5 # an int length would wrap to 5
  y <- numeric(n)
  y[n] <- -1
  expect_equal(LPM(1, 0, y, excess_ret = TRUE), 1 / n)
  expect_equal(UPM(0, -0.5, y, excess_ret = TRUE), (n - 1) / n)
})

test_that("VaR rejects a factor column inside a data.table", {
  skip_if_not_installed("data.table")
  dt <- data.table::data.table(a = factor(c(10, 20, 30)))
  expect_error(LPM.VaR(0.5, 1, dt), "`x` must be numeric, not a factor.", fixed = TRUE)
})
