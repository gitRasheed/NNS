# Input contract of the partial-moment, VaR, CDF and dominance functions (BUGS.md fixes).

# Lean: none | BUGS: #1, #17
test_that("I01: partial moments reject infinite, negative and missing degrees by name (BUGS #1, #17)", {
  withr::local_seed(20261108)
  aud_case("infinite degree is rejected by the univariate partial moments", {
    x <- c(0, 0.5, 2)
    aud_errors("`degree` must be finite.",
      LPM(Inf, 1, x),
      UPM(Inf, 1, x),
      LPM(Inf, c(1, aud_pad), x),
      LPM(Inf, 1, x, excess_ret = TRUE),
      LPM.ratio(Inf, 1, x),
      UPM.ratio(Inf, 1, x),
      NNS.CDF(x, degree = Inf, plot = FALSE))
  })
  aud_case("infinite degree is rejected by the co- and divergent partial moments", {
    x <- c(0, 0.5, 2)
    expect_error(Co.LPM(Inf, x, x, 1, 1), "`degree_lpm` must be finite.", fixed = TRUE)
    expect_error(Co.LPM(1, x, x, 1, 1, degree_y = Inf), "`degree_y` must be finite.", fixed = TRUE)
    aud_errors("`degree_upm` must be finite.",
      Co.UPM(Inf, 3, 0, 1, 1),
      D.LPM(1, Inf, x, x, 1, 1))
    expect_error(D.UPM(Inf, 1, x, x, 1, 1), "`degree_lpm` must be finite.", fixed = TRUE)
    expect_error(PM.matrix(Inf, 1, target = "mean", variable = cbind(1:3, c(3, 1, 2)), pop_adj = TRUE),
                 "`LPM_degree` must be finite.", fixed = TRUE)
  })
  aud_case("infinite degree is rejected by the n-dimensional partial moments", {
    m <- rbind(c(3, 1.5), c(0, 1), c(2, 2))
    aud_errors("`degree` must be finite.",
      Co.LPM_nD(m, c(1, 1), Inf),
      Co.UPM_nD(m, c(1, 1), Inf),
      DPM_nD(m, c(1, 1), Inf),
      NNS.CDF(m, degree = Inf, plot = FALSE))
  })
  aud_case("negative and missing degrees keep the non-negative error", {
    x <- c(0, 0.5, 2)
    aud_errors("'degree' must be a non-negative number",
      LPM(-Inf, 1, x),
      LPM(NaN, 1, x))
  })
})

# Lean: none | BUGS: #3, #7, b724e92
test_that("I02: excess_ret must be TRUE or FALSE, and a zero-length target with it is a length error (BUGS #3, #7)", {
  withr::local_seed(20261108)
  aud_case("excess_ret = NA is an error for a scalar target", {
    aud_errors("`excess_ret` must be TRUE or FALSE.",
      LPM(1, 2, aud_x4, excess_ret = NA),
      UPM(1, 2, aud_x4, excess_ret = NA))
  })
  aud_case("excess_ret = NA gives the same error for empty and vector targets", {
    aud_errors("`excess_ret` must be TRUE or FALSE.",
      LPM(1, numeric(0), aud_x4, excess_ret = NA),
      LPM(1, c(1, 2), aud_x4, excess_ret = NA),
      UPM(1, c(1, 2), aud_x4, excess_ret = NA))
  })
  aud_case("excess_ret rejects vectors and non-logical types", {
    aud_errors("`excess_ret` must be TRUE or FALSE.",
      LPM(1, 2, aud_x4, excess_ret = c(TRUE, FALSE)),
      UPM(1, 2, aud_x4, excess_ret = "yes"),
      LPM(1, 2, aud_x4, excess_ret = NULL))
  })
  aud_case("excess_ret = 1 still means TRUE", {
    expect_identical(LPM(1, 2, aud_x4, excess_ret = 1), LPM(1, 2, aud_x4, excess_ret = TRUE))
    expect_identical(UPM(1, 2, aud_x4, excess_ret = 0), UPM(1, 2, aud_x4))
  })
  aud_case("a zero-length target with excess_ret = TRUE is a length error", {
    expect_error(LPM(1, numeric(0), c(1, 2, 3), excess_ret = TRUE), "target must be length 1 or same length as variable", fixed = TRUE)
  })
})

# Lean: none | BUGS: #7
# Skipped unless NNS_BIGMEM_GB >= 64, so the #7 failures before the fix that always run are in I02 and I06.
test_that("I03: LPM/UPM excess_ret handle vectors longer than .Machine$integer.max (BUGS #7)", {
  withr::local_seed(20261108)
  skip_on_cran()
  skip_if_not(isTRUE(as.numeric(Sys.getenv("NNS_BIGMEM_GB", "0")) >= 64),
              "needs about 48 GiB; set NNS_BIGMEM_GB to run")
  n <- 2^32 + 5 # an int length would wrap to 5 and miss the last element
  y <- numeric(n)
  # The only non-zero deviation is the last element, of size n, so the exact
  # answer is n / n = 1 and the truncated one is 0.
  y[n] <- -n
  expect_identical(LPM(1, 0, y, excess_ret = TRUE), 1)
  y[n] <- n
  expect_identical(UPM(1, 0, y, excess_ret = TRUE), 1)
})

# Lean: none | BUGS: #4, #6, #9, #16, #26
test_that("I04: targets: NULL is the mean, empty gives numeric(0), logical is coerced, anything else is rejected by name (BUGS #4, #9, #16)", {
  withr::local_seed(20261108)
  aud_case("LPM/UPM with a NULL target use mean(variable), like LPM.ratio", {
    expect_identical(LPM(1, NULL, aud_x4), LPM(1, mean(aud_x4), aud_x4))
    expect_identical(UPM(2, NULL, aud_x5), UPM(2, mean(aud_x5), aud_x5))
    expect_identical(LPM(1, NULL, aud_x4, excess_ret = TRUE), LPM(1, mean(aud_x4), aud_x4, excess_ret = TRUE))
  })
  aud_case("LPM with a NULL target and empty variable reports a non-finite target", {
    expect_error(LPM(1, NULL, numeric(0)), "`target` must be finite.", fixed = TRUE)
  })
  aud_case("a zero-length target still returns numeric(0)", {
    expect_identical(LPM(1, numeric(0), aud_x4), numeric(0))
    expect_identical(UPM(1, numeric(0), aud_x4), numeric(0))
  })
  aud_case("LPM/UPM reject a factor target", {
    aud_errors("`target` must be numeric, not a factor.",
      LPM(1, factor("2"), c(1, 2, 3)),
      UPM(1, factor("2"), c(1, 2, 3)))
  })
  aud_case("co-moments reject a factor target", {
    expect_error(Co.LPM(1, aud_x5, aud_x5, factor("z"), 1), "`target_x` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("a logical target is coerced like LPM(), not replaced by the mean", {
    expect_identical(LPM.ratio(0, TRUE, aud_x5), LPM.ratio(0, 1, aud_x5))
    expect_identical(UPM.ratio(0, c(TRUE, FALSE), aud_x5), UPM.ratio(0, c(1, 0), aud_x5))
    expect_identical(Co.LPM(0, aud_x5, aud_x5, TRUE, TRUE), Co.LPM(0, aud_x5, aud_x5, 1, 1))
    expect_identical(LPM.ratio(0, TRUE, aud_x5), LPM(0, TRUE, aud_x5))
  })
  aud_case("a logical NA target is non-finite, as NA_real_ is", {
    aud_errors("`target` must be finite.",
      LPM.ratio(0, NA, aud_x5),
      UPM.ratio(0, NA, aud_x5))
    expect_error(Co.LPM(1, aud_x5, aud_x5, NA, 1), "`target_x` must be finite.", fixed = TRUE)
    expect_error(D.LPM(1, 1, aud_x5, aud_x5, 1, NA), "`target_y` must be finite.", fixed = TRUE)
  })
  aud_case("character and list targets are errors in the ratio functions", {
    aud_errors("`target` must be numeric or NULL.",
      LPM.ratio(0, "a", aud_x5),
      UPM.ratio(0, "1", aud_x5),
      LPM.ratio(0, list(1), aud_x5))
  })
  aud_case("character and list targets are errors in the co- and divergent moments", {
    expect_error(Co.LPM(1, aud_x5, aud_x5, "3", 4), "`target_x` must be numeric or NULL.", fixed = TRUE)
    expect_error(Co.UPM(1, aud_x5, aud_x5, 1, "3"), "`target_y` must be numeric or NULL.", fixed = TRUE)
    expect_error(D.LPM(1, 1, aud_x5, aud_x5, list(0), 1), "`target_x` must be numeric or NULL.", fixed = TRUE)
    expect_error(D.UPM(1, 1, aud_x5, aud_x5, 1, list(0)), "`target_y` must be numeric or NULL.", fixed = TRUE)
  })
  aud_case("LPM/UPM reject character and list targets like the ratio functions", {
    aud_errors("`target` must be numeric or NULL.",
      LPM(1, "1", aud_x5),
      UPM(1, list(1), aud_x5))
  })
  aud_case("NULL and integer targets are unchanged", {
    expect_identical(LPM.ratio(0, NULL, aud_x5), LPM.ratio(0, mean(aud_x5), aud_x5))
    expect_identical(LPM.ratio(0, 2L, aud_x5), LPM.ratio(0, 2, aud_x5))
    expect_identical(Co.LPM(1, aud_x5, rev(aud_x5), NULL, NULL), Co.LPM(1, aud_x5, rev(aud_x5), mean(aud_x5), mean(rev(aud_x5))))
    expect_identical(Co.LPM(1, aud_x5, rev(aud_x5), 2L, 2L), Co.LPM(1, aud_x5, rev(aud_x5), 2, 2))
  })
  aud_case("PM.matrix rejects non-numeric targets", {
    A <- cbind(a = aud_x5, b = rev(aud_x5))
    aud_errors("`target` must be numeric or NULL.",
      PM.matrix(1, 1, c("0", "0"), A, TRUE),
      PM.matrix(1, 1, data.frame(0, 0), A, TRUE),
      PM.matrix(1, 1, "median", A, TRUE))
    aud_errors("`target` must be finite.",
      PM.matrix(1, 1, c(NA, NA), A, TRUE),
      PM.matrix(1, 1, NA, matrix(rnorm(9), 3), TRUE),
      PM.matrix(1, 1, Inf, A, TRUE))
    expect_error(PM.matrix(1, 1, c(0, 0, 0), A, TRUE), "variable matrix cols != target vector length", fixed = TRUE)
  })
  aud_case("PM.matrix keeps target = \"mean\"/\"means\" as an alias for NULL", {
    A <- cbind(a = aud_x5, b = rev(aud_x5))
    expect_identical(PM.matrix(1, 1, "mean", A, TRUE), PM.matrix(1, 1, NULL, A, TRUE))
    expect_identical(PM.matrix(1, 1, "means", A, TRUE), PM.matrix(1, 1, NULL, A, TRUE))
  })
  aud_case("PM.matrix integer targets are unchanged", {
    A <- cbind(a = aud_x5, b = rev(aud_x5))
    expect_identical(PM.matrix(1, 1, c(2L, 3L), A, FALSE), PM.matrix(1, 1, c(2, 3), A, FALSE))
  })
  aud_case("nD functions reject a factor target", {
    m <- cbind(aud_x5, rev(aud_x5))
    expect_error(Co.LPM_nD(m, factor(c(1, 2)), 1), "`target` must be numeric, not a factor.", fixed = TRUE)
    expect_error(DPM_nD(m, c("1", "2"), 1), "`target` must be numeric.", fixed = TRUE)
  })
  aud_case("NNS.CDF validates a non-numeric or missing target", {
    expect_error(NNS.CDF(aud_x5, target = NA, plot = FALSE), "`target` must be finite.", fixed = TRUE)
    expect_error(NNS.CDF(aud_x5, target = "2", plot = FALSE), "`target` must be numeric or NULL.", fixed = TRUE)
    expect_identical(NNS.CDF(aud_x5, target = TRUE, plot = FALSE)$target.value, LPM.ratio(0, 1, aud_x5))
  })
})

# Lean: none | BUGS: #4, #6, #26, 8e6e9b0, b724e92, daf6261
test_that("I05: data: factors and non-numeric data are rejected by name; logical, integer and numeric frames are accepted (BUGS #6, #26)", {
  withr::local_seed(20261108)
  aud_case("LPM/UPM reject a factor variable instead of using its level codes", {
    f <- factor(c("10", "20", "30"))
    aud_errors("`variable` must be numeric, not a factor.",
      LPM(1, 15, f),
      UPM(1, 15, f))
  })
  aud_case("ratios and co-moments reject factor data", {
    f <- factor(c("10", "20", "30"))
    aud_errors("`variable` must be numeric, not a factor.",
      LPM.ratio(1, 15, f),
      UPM.ratio(1, 15, f))
    expect_error(Co.LPM(1, f, c(1, 2, 3), 1, 1), "`x` must be numeric, not a factor.", fixed = TRUE)
    expect_error(D.UPM(1, 1, c(1, 2, 3), f, 1, 1), "`y` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("a data.frame holding a factor column is rejected like a bare factor", {
    expect_error(LPM.ratio(1, 15, data.frame(a = factor(c(10, 20, 30)))),
                 "`variable` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("a factor column of a multi-column data.frame is rejected, not flattened to codes", {
    mixed <- data.frame(a = 1:2, b = factor(c("x", "y")))
    first <- data.frame(b = factor(c("x", "y")), a = 1:2)
    for (df in list(mixed, first)) {
      aud_errors("`variable` must be numeric, not a factor.",
        LPM.ratio(0, 2, df),
        UPM.ratio(0, 2, df))
      expect_error(Co.LPM(1, df, df, 0, 0), "`x` must be numeric, not a factor.", fixed = TRUE)
      expect_error(Co.UPM(1, 1:4, df, 0, 0), "`y` must be numeric, not a factor.", fixed = TRUE)
      expect_error(D.LPM(1, 1, df, 1:4, 0, 0), "`x` must be numeric, not a factor.", fixed = TRUE)
      expect_error(D.UPM(1, 1, 1:4, df, 0, 0), "`y` must be numeric, not a factor.", fixed = TRUE)
      expect_error(PM.matrix(1, 1, NULL, df, TRUE), "`variable` must be numeric, not a factor.", fixed = TRUE)
    }
    numeric_df <- data.frame(a = 1:2, b = 3:4)
    expect_identical(LPM.ratio(0, 2, numeric_df), LPM.ratio(0, 2, c(1, 2, 3, 4)))
    expect_identical(Co.LPM(1, numeric_df, numeric_df, 2, 2), Co.LPM(1, 1:4, 1:4, 2, 2))
  })
  aud_case("NNS.CDF rejects a factor variable", {
    expect_error(NNS.CDF(factor(c(10, 20, 30)), plot = FALSE),
                 "`variable` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("NNS.moments rejects a factor instead of using its level codes", {
    expect_error(NNS.moments(factor(c(10, 20, 40))), "`x` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("stochastic dominance functions reject factor input", {
    f <- factor(c(10, 20))
    expect_error(NNS.FSD.uni(f, c(1.5, 1.5)), "`x` must be numeric, not a factor.", fixed = TRUE)
    expect_error(NNS.SSD.uni(c(1.5, 1.5), f), "`y` must be numeric, not a factor.", fixed = TRUE)
    aud_errors("`x` must be numeric, not a factor.",
      NNS.TSD.uni(f, c(1.5, 1.5)),
      NNS.FSD(f, c(1.5, 1.5), plot = FALSE))
  })
  aud_case("stochastic dominance rejects a factor column inside a data.frame", {
    expect_error(NNS.SSD(data.frame(a = factor(c(10, 20))), c(1.5, 1.5), plot = FALSE),
                 "`x` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("stochastic dominance still flattens single-column numeric containers", {
    expect_identical(NNS.FSD.uni(data.frame(a = c(10, 20)), c(1.5, 1.5)),
                     NNS.FSD.uni(c(10, 20), c(1.5, 1.5)))
    expect_identical(NNS.SSD.uni(matrix(c(10, 20), ncol = 1), c(1.5, 1.5)), NNS.SSD.uni(c(10, 20), c(1.5, 1.5)))
  })
  aud_case("LPM keeps accepting logical and integer data", {
    expect_identical(LPM(1, 2, c(TRUE, FALSE)), 1.5)
    expect_identical(LPM(1, 2, 1:3), LPM(1, 2, c(1, 2, 3)))
  })
  aud_case("LPM/UPM reject a NULL variable, like LPM.ratio", {
    aud_errors("`variable` must be numeric.",
      LPM(1, 2, NULL),
      UPM(1, 2, NULL))
  })
  aud_case("the compiled ratio and co-moment kernels reject non-numeric data in the R wording", {
    # Same message as .nns_as_numeric() on the R side (LPM(), LPM.VaR(), ...).
    aud_errors("`variable` must be numeric.",
      LPM.ratio(1, 1, "a"),
      UPM.ratio(1, 1, list(1, 2)),
      LPM.ratio(1, 2, NULL))
    expect_error(Co.LPM(1, "a", aud_x5, 1, 1), "`x` must be numeric.", fixed = TRUE)
    expect_error(Co.UPM(1, aud_x5, c("1", "2", "3", "4", "5"), 1, 1), "`y` must be numeric.", fixed = TRUE)
    expect_error(D.LPM(1, 1, "a", aud_x5, 1, 1), "`x` must be numeric.", fixed = TRUE)
  })
  aud_case("LPM/UPM keep NaN for an empty numeric variable", {
    expect_identical(LPM(1, 2, numeric(0)), NaN)
    expect_identical(UPM(1, 2, numeric(0)), NaN)
  })
  aud_case("nD functions reject a factor column in data", {
    df <- data.frame(a = factor(c(1, 2, 3)), b = c(1, 2, 3))
    aud_errors("`data` must be numeric, not a factor.",
      Co.LPM_nD(df, c(1, 1), 1),
      Co.UPM_nD(df, c(1, 1), 1),
      DPM_nD(df, c(1, 1), 1))
  })
  aud_case("nD functions reject bare factor data", {
    expect_error(Co.LPM_nD(factor(c(1, 2, 3)), 1, 1), "`data` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("nD functions still accept a numeric data frame", {
    df <- data.frame(a = c(-1, 0, 2), b = c(1, -2, 3))
    expect_identical(Co.LPM_nD(df, c(0, 0), 1, FALSE), Co.LPM_nD(as.matrix(df), c(0, 0), 1, FALSE))
  })
  aud_case("multivariate NNS.CDF rejects a factor column", {
    df <- data.frame(a = factor(c(1, 2, 3)), b = c(1, 2, 3))
    expect_error(NNS.CDF(df, plot = FALSE), "`variable` must be numeric, not a factor.", fixed = TRUE)
  })
})

# Lean: none | BUGS: #6, #7, #10, 4217b25, a08d04d, 61bc773, 86dd679, 8e6e9b0
test_that("I06: LPM.VaR and UPM.VaR input contract: data, percentile and degree (BUGS #10)", {
  withr::local_seed(20261108)
  aud_case("LPM.VaR rejects +Inf the same way at every degree", {
    y <- c(1, 2, Inf)
    aud_errors("`x` must be finite.",
      LPM.VaR(0.5, 0, y),
      LPM.VaR(0.5, 1, y),
      LPM.VaR(0.5, 1.5, y))
  })
  aud_case("UPM.VaR rejects -Inf the same way at every degree", {
    y <- c(-Inf, 1, 2)
    aud_errors("`x` must be finite.",
      UPM.VaR(0.5, 0, y),
      UPM.VaR(0.5, 1, y),
      UPM.VaR(0.5, 2.5, y))
  })
  aud_case("VaR still drops NA and NaN", {
    withr::local_seed(7)
    r <- rnorm(100, 0, 0.01)
    expect_identical(LPM.VaR(0.05, 1, c(r, NA, NaN)), LPM.VaR(0.05, 1, r))
    expect_identical(UPM.VaR(0.05, 2, c(NA, r)), UPM.VaR(0.05, 2, r))
  })
  aud_case("VaR rejects factor and character x", {
    expect_error(LPM.VaR(0.5, 1, factor(c(1, 2, 3))), "`x` must be numeric, not a factor.", fixed = TRUE)
    expect_error(UPM.VaR(0.5, 1, c("1", "2")), "`x` must be numeric.", fixed = TRUE)
  })
  aud_case("VaR rejects a missing or non-numeric percentile at every degree", {
    # A missing percentile used to give NA at degree 0, a silent 0 at degree 1 to
    # 4 when other percentiles were present, and "missing value where
    # TRUE/FALSE needed" otherwise.
    x <- c(1, 2, 3, 5)
    for (d in c(0, 1, 2, 4, 2.5)) {
      aud_errors("`percentile` must not be NA.",
        LPM.VaR(NA, d, x),
        UPM.VaR(c(0.1, NaN), d, x),
        LPM.VaR(c(0.1, NA), d, x))
      expect_error(UPM.VaR("0.1", d, x), "`percentile` must be numeric.", fixed = TRUE)
    }
    # Out-of-range percentiles are still clamped to [0, 1].
    expect_identical(LPM.VaR(c(-1, 2), 2, x), c(1, 5))
  })
  aud_case("VaR rejects a missing, non-scalar or non-numeric degree by name", {
    # NA and NaN used to fail with "missing value where TRUE/FALSE needed" from
    # the degree == 0 branch.
    x <- c(1, 2, 3, 5)
    msg <- "'degree' must be a non-negative number"
    for (f in list(LPM.VaR, UPM.VaR)) {
      aud_errors(msg,
        f(0.1, NA, x),
        f(0.1, NA_real_, x),
        f(0.1, NaN, x),
        f(0.1, -1, x),
        f(0.1, -Inf, x),
        f(0.1, c(1, 2), x),
        f(0.1, numeric(0), x),
        f(0.1, NULL, x),
        f(0.1, "1", x))
      expect_error(f(0.1, Inf, x), "`degree` must be finite.", fixed = TRUE)
    }
    # The degree is checked before the data, as the compiled kernels do.
    expect_error(LPM.VaR(0.1, NA, c(1, Inf)), msg, fixed = TRUE)
    # Valid degrees, including a logical, are unchanged.
    expect_identical(LPM.VaR(0.1, TRUE, x), LPM.VaR(0.1, 1, x))
    expect_identical(UPM.VaR(0.1, 0, x), stats::quantile(x, 0.9))
  })
  aud_case("VaR rejects a factor column inside a data.table", {
    # Class-tagged data.frame: the VaR input path dispatches on is.data.frame,
    # so this covers data.table input without declaring data.table in Suggests.
    dt <- structure(data.frame(a = factor(c(10, 20, 30))), class = c("data.table", "data.frame"))
    expect_error(LPM.VaR(0.5, 1, dt), "`x` must be numeric, not a factor.", fixed = TRUE)
  })
  aud_case("the VaR anchored sums reject a degree beyond the prefix limit", {
    expect_error(.Call("_NNS_VaR_anchor_sums_RCPP", c(0, 1), c(1, 1), 33L),
                 "invalid anchored-sum input", fixed = TRUE)
    expect_error(.Call("_NNS_VaR_anchor_sums_RCPP", c(0, 1), c(1, 1), .Machine$integer.max),
                 "invalid anchored-sum input", fixed = TRUE)
    expect_identical(dim(.Call("_NNS_VaR_anchor_sums_RCPP", c(0, 1), c(1, 1), 32L)$lower), c(1L, 33L))
  })
  aud_case("VaR treats a one-column data.frame like a numeric vector", {
    expect_identical(LPM.VaR(0.5, 1, data.frame(a = c(1, 2, 3))), LPM.VaR(0.5, 1, c(1, 2, 3)))
    expect_identical(UPM.VaR(0.25, 2, data.frame(a = c(1, 5, 3))), UPM.VaR(0.25, 2, c(1, 5, 3)))
  })
  aud_case("VaR flattens a multi-column data.frame like a multi-column tibble", {
    expect_identical(LPM.VaR(0.5, 1, data.frame(a = c(1, 2), b = c(3, 4))), LPM.VaR(0.5, 1, c(1, 2, 3, 4)))
  })
  aud_case("VaR rejects a factor column of a multi-column data.frame", {
    df <- data.frame(a = factor(c(10, 20)), b = c(1, 2))
    expect_error(LPM.VaR(0.5, 1, df), "`x` must be numeric, not a factor.", fixed = TRUE)
  })
})

# Lean: none | BUGS: #20, #21
test_that("I07: zero-row and one-row data (BUGS #20, #21)", {
  withr::local_seed(20261108)
  aud_case("PM.matrix rejects a zero-row variable whatever the target", {
    A0 <- matrix(numeric(0), 0, 3)
    aud_errors("`variable` must have at least one row",
      PM.matrix(1, 1, c(0, 0, 0), A0, TRUE),
      PM.matrix(1, 1, NULL, A0, TRUE),
      PM.matrix(0, 0, c(0, 0, 0), A0, FALSE, norm = TRUE))
  })
  aud_case("PM.matrix with one row still returns the documented list", {
    expect_named(PM.matrix(1, 1, c(0, 0, 0), matrix(c(1, -1, 2), 1), TRUE),
                 c("cupm", "dupm", "dlpm", "clpm", "cov.matrix"))
  })
  aud_case("NNS.copula on zero-row data still short-circuits to 0", {
    expect_identical(NNS.copula(matrix(numeric(0), 0, 3)), 0)
  })
  aud_case("nD co-partial moments reject zero-row data like the batch backend", {
    m <- matrix(numeric(0), 0, 2)
    aud_errors("`data` must have at least one row",
      Co.LPM_nD(m, c(0, 0), 1, TRUE),
      Co.LPM_nD(m, c(0, 0), 0),
      Co.UPM_nD(m, c(0, 0), 1, FALSE),
      DPM_nD(m, c(0, 0), 1, TRUE),
      NNS:::Co.LPM_nD.batch(m, matrix(0, 1, 2), 1))
  })
  aud_case("one-row nD data agrees with the batch backend", {
    r1 <- matrix(c(-1, -1), 1, 2)
    expect_identical(Co.LPM_nD(r1, c(0, 0), 1, TRUE), NNS:::Co.LPM_nD.batch(r1, matrix(0, 1, 2), 1, TRUE))
  })
})

# Lean: none | BUGS: #13
test_that("I08: NNS.CDF accepts one-column frames and vector targets (BUGS #13)", {
  withr::local_seed(20261108)
  aud_case("NNS.CDF accepts a one-column data.frame like a vector", {
    expect_identical(NNS.CDF(data.frame(a = aud_x5), plot = FALSE), NNS.CDF(aud_x5, plot = FALSE))
    expect_identical(NNS.CDF(data.frame(a = c(aud_x5, NA)), target = 2, plot = FALSE),
                     NNS.CDF(aud_x5, target = 2, plot = FALSE))
  })
  aud_case("NNS.CDF evaluates a univariate vector target element-wise", {
    expect_identical(NNS.CDF(aud_x5, target = c(2, 3), plot = FALSE)$target.value,
                     c(NNS.CDF(aud_x5, target = 2, plot = FALSE)$target.value,
                       NNS.CDF(aud_x5, target = 3, plot = FALSE)$target.value))
    expect_identical(NNS.CDF(aud_x5, target = c(2, 3), type = "survival", plot = FALSE)$target.value,
                     1 - LPM.ratio(0, c(2, 3), aud_x5))
  })
  aud_case("NNS.CDF hazard with a vector target matches the scalar calls", {
    expect_equal(NNS.CDF(aud_x5, target = c(2, 3), type = "hazard", plot = FALSE)$target.value,
                 c(NNS.CDF(aud_x5, target = 2, type = "hazard", plot = FALSE)$target.value,
                   NNS.CDF(aud_x5, target = 3, type = "hazard", plot = FALSE)$target.value))
  })
  aud_case("NNS.CDF still rejects a vector target with an out-of-support value", {
    expect_error(NNS.CDF(aud_x5, target = c(2, 30), plot = FALSE), "target out of bounds", fixed = TRUE)
  })
})

# Lean: none | BUGS: #23
test_that("I09: stochastic dominance: type must be discrete or continuous, any case (BUGS #23)", {
  withr::local_seed(20261108)
  aud_case("NNS.FSD rejects an invalid type instead of running the Omega order", {
    expect_error(NNS.FSD(aud_x_om, aud_y_om, type = "foo", plot = FALSE), aud_type_msg, fixed = TRUE)
  })
  aud_case("NNS.FSD.uni rejects an invalid type instead of running discrete FSD", {
    expect_error(NNS.FSD.uni(aud_x_om, aud_y_om, type = "discret"), aud_type_msg, fixed = TRUE)
  })
  aud_case("NNS.SD.efficient.set rejects an invalid type", {
    expect_error(NNS.SD.efficient.set(cbind(x = aud_x_om, y = aud_y_om), 1, type = "foo", status = FALSE),
                 aud_type_msg, fixed = TRUE)
  })
  aud_case("NNS.SD.efficient.set rejects an invalid type at degree 2, where type is unused", {
    expect_error(NNS.SD.efficient.set(cbind(x = aud_x_om, y = aud_y_om), 2, type = "foo", status = FALSE),
                 aud_type_msg, fixed = TRUE)
  })
  aud_case("NNS.SD.cluster rejects an invalid type", {
    expect_error(NNS.SD.cluster(cbind(x = aud_x_om, y = aud_y_om), 1, type = "foo"), aud_type_msg, fixed = TRUE)
  })
  aud_case("type must be a single non-missing string", {
    aud_errors(aud_type_msg,
      NNS.FSD(aud_x_om, aud_y_om, type = NA, plot = FALSE),
      NNS.FSD(aud_x_om, aud_y_om, type = c("discrete", "continuous"), plot = FALSE),
      NNS.FSD.uni(aud_x_om, aud_y_om, type = 1),
      NNS.FSD.uni(aud_x_om, aud_y_om, type = "Discrete "))
  })
  aud_case("the C++ entry points reject an invalid type", {
    aud_errors(aud_type_msg,
      .Call(NNS:::`_NNS_NNS_FSD_uni_cpp`, aud_x_om, aud_y_om, "foo"),
      NNS:::sd_dom_matrix_prefix_parallel(cbind(aud_x_om, aud_y_om), 1L, "foo"))
  })
  aud_case("valid types are case-insensitive and keep their verdicts", {
    expect_identical(NNS.FSD(aud_x_om, aud_y_om, type = "Continuous", plot = FALSE), "X OMEGA Y")
    expect_identical(NNS.FSD(aud_x_om, aud_y_om, type = "DISCRETE", plot = FALSE), "NO FSD EXISTS")
    expect_identical(NNS.FSD.uni(aud_x_om, aud_y_om, "CONTINUOUS"), 1L)
    expect_identical(NNS.FSD.uni(aud_x_om, aud_y_om, "discrete"), 0L)
    expect_identical(NNS.SD.efficient.set(cbind(x = aud_x_om, y = aud_y_om), 1, "Continuous", status = FALSE), "x")
  })
})

# Lean: none | BUGS: #24, #28
test_that("I10: stochastic dominance rejects infinite observations by name and keeps the NA message (BUGS #24, #28)", {
  withr::local_seed(20261108)
  aud_case("the uni-directional tests reject an infinite observation in x", {
    expect_error(NNS.FSD.uni(c(1, Inf), c(0, 1)), "`x` must be finite.", fixed = TRUE)
  })
  aud_case("the uni-directional tests reject an infinite observation in y", {
    aud_errors("`y` must be finite.",
      NNS.SSD.uni(c(0, 1), c(-Inf, 1)),
      NNS.TSD.uni(c(0, 1), c(-Inf, 1)),
      NNS.FSD.uni(c(0, 1), c(-Inf, 1), "continuous"))
  })
  aud_case("a total-loss log return makes the efficient set an error, not a superset of the FSD set", {
    # The bankrupt asset used to be FSD-dominated but SSD/TSD-efficient.
    x <- diff(log(c(100, 101, 99, 102, 103)))
    y <- diff(log(c(100, 60, 30, 10, 0)))  # ends in -Inf
    m <- cbind(x = x, y = y)
    aud_errors("`x` must be finite.",
      NNS.SD.efficient.set(m, 1, status = FALSE),
      NNS.SD.efficient.set(m, 2, status = FALSE),
      NNS.SD.efficient.set(m, 3, status = FALSE))
  })
  aud_case("NNS.SD.cluster names `data` for an infinite observation", {
    m <- cbind(a = c(1, Inf, 2), b = c(0, 1, 1), c = c(-1, 0, 0))
    expect_error(NNS.SD.cluster(m, 2), "`data` must be finite.", fixed = TRUE)
  })
  aud_case("the dominance-matrix kernel rejects an infinite observation", {
    m <- cbind(a = c(1, Inf, 2), b = c(0, 1, 1))
    expect_error(NNS:::sd_dom_matrix_prefix_parallel(m, 2L, "discrete"), "`x` must be finite.", fixed = TRUE)
  })
  aud_case("NA and NaN keep the missing-values message", {
    msg <- "You have some missing values, please address."
    aud_errors(msg,
      NNS.TSD.uni(c(1, NaN), 1:3),
      NNS.SSD(c(1, NA, Inf), 1:3, plot = FALSE),
      NNS.SD.efficient.set(cbind(a = c(1, NA), b = c(Inf, 2)), 2, status = FALSE),
      NNS.SD.cluster(cbind(a = c(1, NA), b = c(Inf, 2)), 2))
  })
  aud_case("a large finite loss is still dominated at every order", {
    x <- c(0, 1, 3); y <- c(0, 1, -50)
    expect_identical(c(NNS.FSD.uni(x, y), NNS.SSD.uni(x, y), NNS.TSD.uni(x, y)), c(1L, 1L, 1L))
  })
  aud_case("the bidirectional tests name x for an infinite observation", {
    aud_errors("`x` must be finite.",
      NNS.FSD(c(1, Inf), c(0, 1), plot = FALSE),
      NNS.FSD(c(1, Inf), c(0, 1), type = "continuous", plot = FALSE),
      NNS.SSD(c(1, Inf), c(0, 1), plot = FALSE),
      NNS.TSD(c(1, Inf), c(0, 1), plot = FALSE))
  })
  aud_case("the bidirectional tests name y for an infinite observation", {
    aud_errors("`y` must be finite.",
      NNS.FSD(c(0, 1), c(1, -Inf), plot = FALSE),
      NNS.SSD(c(0, 1), c(1, -Inf), plot = FALSE),
      NNS.TSD(c(0, 1), c(1, -Inf), plot = FALSE))
  })
  aud_case("the bidirectional tests reject infinite input the same way with plot = TRUE", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    expect_error(NNS.SSD(c(1, Inf), c(0, 1), plot = TRUE), "`x` must be finite.", fixed = TRUE)
  })
})

# Lean: none | BUGS: ac35c82, 4007104
test_that("I11: stochastic dominance set and cluster degree must be 1, 2 or 3", {
  withr::local_seed(20261108)
  aud_case("NNS.SD.efficient.set rejects a fractional degree instead of truncating it", {
    m <- cbind(x = aud_x_om, y = aud_y_om)
    for (d in list(1.5, 2.9, 0.5, 3.0000001)) {
      expect_error(NNS.SD.efficient.set(m, d, status = FALSE), aud_degree_msg, fixed = TRUE)
    }
  })
  aud_case("NNS.SD.efficient.set rejects other invalid degrees with the same message", {
    m <- cbind(x = aud_x_om, y = aud_y_om)
    for (d in list(0, 4, -1, NA, NA_real_, NaN, Inf, c(1, 2), numeric(0), "2", TRUE, factor(2), NULL)) {
      expect_error(NNS.SD.efficient.set(m, d, status = FALSE), aud_degree_msg, fixed = TRUE)
    }
    # Beyond the integer range: an error, not a coercion warning.
    expect_no_warning(expect_error(NNS.SD.efficient.set(m, 1e10, status = FALSE), aud_degree_msg, fixed = TRUE))
  })
  aud_case("NNS.SD.efficient.set accepts 1, 2, 3 as double or integer with unchanged results", {
    withr::local_seed(7)
    A <- cbind(a = rnorm(40), b = rnorm(40, 0.3), c = rexp(40) - 1)
    for (d in 1:3) {
      ref <- .Call(NNS:::`_NNS_NNS_SD_efficient_set_parallel_cpp`, A, d, "discrete", FALSE)
      expect_identical(NNS.SD.efficient.set(A, d, status = FALSE), ref)
      expect_identical(NNS.SD.efficient.set(A, as.numeric(d), status = FALSE), ref)
    }
  })
  aud_case("NNS.SD.cluster rejects an invalid degree even when no efficient set is computed", {
    m <- cbind(x = aud_x_om, y = aud_y_om)
    expect_error(NNS.SD.cluster(m, 1.5), aud_degree_msg, fixed = TRUE)
    # min_cluster >= ncol: the efficient set is never called.
    expect_error(NNS.SD.cluster(m, 1.5, min_cluster = 2), aud_degree_msg, fixed = TRUE)
    expect_identical(NNS.SD.cluster(m, 1, min_cluster = 2)$Clusters,
                     NNS.SD.cluster(m, 1L, min_cluster = 2)$Clusters)
  })
  aud_case("the compiled efficient set checks the degree before its LPM ordering pass", {
    m <- cbind(x = aud_x_om, y = aud_y_om)
    # Before the fix, the ordering pass ran 2147483647 multiplications per
    # element before the degree was checked (minutes before the error).
    expect_error(.Call(NNS:::`_NNS_NNS_SD_efficient_set_parallel_cpp`, m, .Machine$integer.max,
                       "discrete", FALSE),
                 "degree must be 1, 2, or 3", fixed = TRUE)
    # No wall-clock bound: the zero-column case below fails on the old order by
    # itself (the zero-column early return, which precedes the ordering pass,
    # came before the check and gave character(0)).
    expect_error(.Call(NNS:::`_NNS_NNS_SD_efficient_set_parallel_cpp`, m[, 0, drop = FALSE], 4L,
                       "discrete", FALSE),
                 "degree must be 1, 2, or 3", fixed = TRUE)
  })
})

test_that("I12: integer targets and integer data are used as given by every partial-moment wrapper and NNS.CDF", {
  withr::local_seed(20261108)
  pm_ratio_targets <- fx_pm_ratio_targets()
  aud_case("LPM.ratio and UPM.ratio use integer targets as given", env = pm_ratio_targets, {
    expect_equal(LPM.ratio(0, 2L, x), 0.5)
    expect_equal(LPM.ratio(0, 1:3, x), c(0.25, 0.5, 0.75))
    expect_equal(UPM.ratio(0, 1:3, x), c(0.75, 0.5, 0.25))
    expect_equal(LPM.ratio(1, 1:3, x), LPM.ratio(1, c(1, 2, 3), x))
  })
  aud_case("Co.LPM, Co.UPM, D.LPM and D.UPM use integer targets as given", env = pm_ratio_targets, {
    expect_equal(Co.LPM(1, x, x, 2L, 2L), 0.25)
    expect_equal(Co.UPM(1, x, x, 2L, 2L), Co.UPM(1, x, x, 2, 2))
    expect_equal(D.LPM(1, 1, x, rev(x), 0L, 0L), 0)
    expect_equal(D.UPM(1, 1, x, rev(x), 3L, 2L), D.UPM(1, 1, x, rev(x), 3, 2))
  })
  aud_case("PM.matrix uses an integer target as given", env = pm_ratio_targets, {
    A <- cbind(x, rev(x))
    expect_equal(PM.matrix(1, 1, c(0L, 0L), A, FALSE)$cupm[1, 1], mean(x^2))
    expect_equal(PM.matrix(1, 1, c(2L, 3L), A, FALSE), PM.matrix(1, 1, c(2, 3), A, FALSE))
  })
  aud_case("NNS.CDF of an integer vector is the ECDF", {
    expect_equal(NNS.CDF(1:10, plot = FALSE)$Function$CDF, (1:10) / 10)
  })
  aud_case("LPM.ratio and UPM.ratio use integer targets", {
    x <- c(0, 1, 2, 3, 10)
    expect_equal(LPM.ratio(0, 1L, x), 0.4)
    expect_equal(UPM.ratio(1, c(1L, 2L), x), UPM.ratio(1, c(1, 2), x))
  })
  aud_case("NNS.CDF of integer data matches the same data as double", {
    x <- rep(0:2, each = 20)
    expect_equal(NNS.CDF(x, plot = FALSE)$Function[[2]], rep(c(1, 2, 3) / 3, each = 20))
    expect_equal(NNS.CDF(x, type = "hazard", plot = FALSE)$Function[[2]],
                 NNS.CDF(as.double(x), type = "hazard", plot = FALSE)$Function[[2]])
  })
})

test_that("I13: multivariate NNS.CDF target: one value per column, bounds-checked on every coordinate", {
  withr::local_seed(20261108)
  A <- cbind(c(0, 1, 2), c(0, 1, 2), c(0, 1, 2))
  aud_errors("target out of bounds", NNS.CDF(A, 0, target = c(0, 0, 1e6), type = "survival", plot = FALSE))
  aud_errors("one value per column", NNS.CDF(A, 0, target = c(1, 1), plot = FALSE))
})

test_that("I14: co-moments reject unequal lengths, empty input and recycled targets; non-finite input is rejected by the ratio and co-moment wrappers", {
  withr::local_seed(20261108)
  aud_errors("must be finite",
             LPM.ratio(0, 2, c(1, 2, NA)),
             UPM.ratio(0, Inf, 1:3),
             Co.LPM(0, c(1, 2, NaN), 1:3, 2, 2),
             Co.UPM(0, 1:3, 1:3, 2, -Inf),
             D.LPM(0, 0, 1:3, c(1, NA, 3), 2, 2),
             D.UPM(0, 0, 1:3, 1:3, NA_real_, 2),
             PM.matrix(1, 1, c(2, 2), cbind(1:3, c(1, 2, NA)), FALSE))
  aud_case("co-moments reject unequal lengths, empty input and recycled targets", {
    set.seed(1)
    x <- rnorm(100); y <- rnorm(99)
    fns <- list(
      function(x, y, tx, ty) Co.LPM(1, x, y, tx, ty),
      function(x, y, tx, ty) Co.UPM(1, x, y, tx, ty),
      function(x, y, tx, ty) D.LPM(1, 1, x, y, tx, ty),
      function(x, y, tx, ty) D.UPM(1, 1, x, y, tx, ty)
    )
    for (f in fns) {
      expect_error(f(x, y, 0, 0), "same length")
      expect_error(f(x, y, NULL, NULL), "same length")
      expect_error(f(x, y, seq(-1, 1, length.out = 200), 0), "same length")
      expect_error(f(numeric(0), numeric(0), 0, 0), "must not be empty")
      expect_error(f(numeric(0), numeric(0), NULL, NULL), "must not be empty")
      expect_error(f(1:5, 1:5, numeric(0), 0), "must not be empty")
      expect_error(f(1:5, 1:5, 0, numeric(0)), "must not be empty")
      expect_error(f(1:5, 1:5, numeric(0), numeric(0)), "must not be empty")
      expect_error(f(x, x, c(0, 1, 2), c(0, 1)), "length 1")
      expect_equal(f(x, x, c(-1, 0, 1), 0.5),
                   c(f(x, x, -1, 0.5), f(x, x, 0, 0.5), f(x, x, 1, 0.5)))
    }
  })
})

test_that("I15: negative, NA and NaN degrees are rejected for scalar and vector targets and by every wrapper", {
  withr::local_seed(20261108)
  aud_errors("non-negative",
             LPM(-1, 0, c(-2, 1)),
             UPM(NaN, 0, c(-2, 1)),
             LPM.ratio(-0.5, 0, c(-2, 1)),
             Co.LPM(-1, c(-2, 1), c(-2, 1), 0, 0),
             D.UPM(1, -2, c(-2, 1), c(-2, 1), 0, 0),
             PM.matrix(-1, 1, NULL, matrix(rnorm(20), 10), TRUE),
             Co.LPM_nD(matrix(rnorm(20), 10), c(0, 0), -1))
  aud_case("a fractional degree is accepted", {
    expect_equal(LPM(0.5, 0, c(-4, 1)), 1)
  })
  aud_case("vector-target LPM / UPM reject negative, NA and NaN degrees", {
    x <- c(-2, 1)
    for (t in list(c(0, 1), rep(0, 33))) {
      expect_error(LPM(-1, t, x), "non-negative")
      expect_error(UPM(-1, t, x), "non-negative")
      expect_error(LPM(NaN, t, x), "non-negative")
      expect_error(UPM(NA_real_, t, x), "non-negative")
    }
  })
})
