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
  # The solver resolves the threshold itself (test-var-root-tolerance.R checks
  # it recovers these order statistics); this checks the absolute ratio
  # residual at p near 1e-21, which used to be the documented contract.
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

# --- PM.matrix(norm = TRUE): zero and unsafe cell totals ----------------------

pm_quadrants <- c("cupm", "dupm", "dlpm", "clpm")
pm_cell_sum <- function(P) P$cupm + P$dupm + P$dlpm + P$clpm
pm_expect_close <- function(actual, expected, tol) {
  expect_lte(max(abs(unname(unlist(actual)) - unname(unlist(expected)))), tol)
}

test_that("PM.matrix(norm = TRUE) gives a column constant at its target the degree-0 cells", {
  M <- cbind(a = c(1, 2, 3), b = c(5, 5, 5))
  P <- PM.matrix(1, 1, NULL, M, FALSE, TRUE)
  # LPM.ratio convention: a constant variable is all lower.
  expect_identical(P$clpm["b", "b"], 1)
  expect_identical(c(P$cupm["b", "b"], P$dupm["b", "b"], P$dlpm["b", "b"]), c(0, 0, 0))
  expect_equal(unname(P$cov.matrix), matrix(c(1, 1 / 3, 1 / 3, 1), 2))
  P0 <- PM.matrix(0, 0, NULL, M, FALSE, TRUE)
  for (q in c(pm_quadrants, "cov.matrix")) {
    expect_equal(P[[q]]["b", ], P0[[q]]["b", ])
    expect_equal(P[[q]][, "b"], P0[[q]][, "b"])
  }
  expect_equal(pm_cell_sum(P), matrix(1, 2, 2, dimnames = dimnames(P$cupm)))
})

test_that("PM.matrix(norm = TRUE) uses the degree-0 cells when every row ties in one of the two columns", {
  # Neither column is constant, but each row ties the target in a or b.
  M <- cbind(a = c(0, 0, 1, -2), b = c(3, -1, 0, 0))
  P <- PM.matrix(2, 1, c(0, 0), M, TRUE, TRUE)
  P0 <- PM.matrix(0, 0, c(0, 0), M, FALSE, FALSE)
  expect_equal(c(P$cupm["a", "b"], P$dupm["a", "b"], P$dlpm["a", "b"], P$clpm["a", "b"]),
               c(0, 0.25, 0.25, 0.5))
  for (q in c(pm_quadrants, "cov.matrix")) expect_equal(P[[q]]["a", "b"], P0[[q]]["a", "b"])
  # The crossed mirror holds for the fallback cell too.
  expect_identical(P$dupm["b", "a"], P$dlpm["a", "b"])
  expect_identical(P$dlpm["b", "a"], P$dupm["a", "b"])
  expect_equal(pm_cell_sum(P), matrix(1, 2, 2, dimnames = dimnames(P$cupm)))
})

test_that("PM.matrix(norm = TRUE) at LPM_degree = 0 counts a tie as lower, unchanged", {
  M <- cbind(a = c(1, 2, 3), b = c(5, 5, 5))
  P <- PM.matrix(0, 1, NULL, M, FALSE, TRUE)
  expect_equal(unname(P$cov.matrix), matrix(c(1, 1 / 3, 1 / 3, 1), 2))
})

test_that("PM.matrix(norm = TRUE) is scale-free at equal degrees when the products underflow or overflow", {
  set.seed(2)
  X <- cbind(rnorm(30), rexp(30), runif(30))
  for (d in 1:2) {
    ref <- unlist(PM.matrix(d, d, NULL, X, FALSE, TRUE))
    for (s in c(1e-200, 2^-700, 1e200)) {
      P <- PM.matrix(d, d, NULL, X * s, FALSE, TRUE)
      pm_expect_close(unlist(P), ref, 16 * .Machine$double.eps)
    }
  }
  # Old behaviour: every product underflowed, all cells 0.
  P <- PM.matrix(1, 1, NULL, X * 1e-200, FALSE, TRUE)
  pm_expect_close(diag(P$cov.matrix), rep(1, 3), 4 * .Machine$double.eps)
})

test_that("PM.matrix(norm = TRUE) at degree 200 does not return 0 when the products overflow", {
  set.seed(2)
  X <- cbind(rnorm(30), rexp(30), runif(30))
  P <- PM.matrix(200, 200, NULL, X * 1e3, FALSE, TRUE)
  expect_true(all(is.finite(unlist(P))))
  pm_expect_close(pm_cell_sum(P), 1, 8 * .Machine$double.eps)
  pm_expect_close(unlist(P), unlist(PM.matrix(200, 200, NULL, X, FALSE, TRUE)), 1e-14)
})

test_that("PM.matrix(norm = TRUE) with unequal degrees matches the exact ratio when the total overflows", {
  # Exact values from rational arithmetic (nns-proofs oracle pm_matrix), in
  # unlist() order: cupm, dupm, dlpm, clpm, each column-major. The plain sums
  # overflow (2^2002), so every cell goes through the rescue.
  X <- cbind(a = c(-1.75, -0.5, 0.25, 1.5, 2, -1.25, 0.75, 1.25),
             b = c(1.5, -2, 0.5, -0.75, 1.25, 0.25, -1.5, 1.75))
  exact <- list(
    "1000 1001" = c(1, 1.0778285412673956e-21, 1.0778285412673956e-21, 3.1782067229759299e-116,
                    0, 0, 1, 0, 0, 1, 0, 0, 2.5944544677354528e-117, 0, 0, 1),
    "700 650" = c(1, 1.1886477434465246e-26, 1.1886477434465246e-26, 3.2171355712192349e-106,
                  0, 2.9359596119441445e-243, 1, 0, 0, 1, 2.9359596119441445e-243, 0,
                  8.2086889175245743e-52, 2.594463309812418e-285, 2.594463309812418e-285, 1))
  for (key in names(exact)) {
    d <- as.numeric(strsplit(key, " ")[[1]])
    got <- unname(unlist(PM.matrix(d[1], d[2], c(0, 0), X, TRUE, TRUE)[pm_quadrants]))
    e <- exact[[key]]
    expect_identical(got == 0, e == 0)
    expect_lte(max(abs(got - e)[e != 0] / e[e != 0]), 1e-14)
  }
})

test_that("PM.matrix(norm = TRUE) log-domain rescue keeps a term whose ratio to the column scale is denormal", {
  # Column a's deviations span more than 2^1074, so |d| / s underflows to 0 for
  # the first row; that row is the only one contributing to cell (a, b) (rows 2
  # and 3 tie in b), whose exact value is cupm = 1. Was NaN (0 / 0).
  M <- cbind(a = c(5e-324, 1e300, -1e300), b = c(1, 0, 0))
  P <- PM.matrix(1, 1, c(0, 0), M, FALSE, TRUE)
  expect_identical(c(P$cupm["a", "b"], P$dupm["a", "b"], P$dlpm["a", "b"], P$clpm["a", "b"]),
                   c(1, 0, 0, 0))
  expect_identical(c(P$cupm["a", "a"], P$clpm["a", "a"]), c(0.5, 0.5))
  expect_true(all(is.finite(unlist(P))))
})

test_that("PM.matrix(norm = TRUE) results with a safe total are unchanged", {
  # Bit pattern of the plain division path (no rescue involved).
  set.seed(123)
  A <- cbind(rnorm(100), rnorm(100), rnorm(100))
  for (d in list(c(1, 1), c(2, 1), c(0, 2))) {
    raw <- PM.matrix(d[1], d[2], NULL, A, TRUE, FALSE)
    tot <- raw$cupm + raw$dupm + raw$dlpm + raw$clpm
    P <- PM.matrix(d[1], d[2], NULL, A, TRUE, TRUE)
    # Upper triangle: the lower one is the (crossed) mirror of it.
    ut <- upper.tri(tot, diag = TRUE)
    for (q in pm_quadrants) expect_identical(P[[q]][ut], (raw[[q]] / tot)[ut])
  }
})
