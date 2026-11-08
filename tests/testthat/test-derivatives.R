# dy.dx, dy.d_ and NNS.diff.

test_that("DER01: dy.dx: missing eval.point, row order, repeated and all five bandwidth positions", {
  withr::local_seed(20261108)
  dydx_duplicate_bands <- fx_dydx_duplicate_bands()
  aud_case("dy.dx stops with a clear message when eval.point is missing", {
    x <- seq(-2, 2, length.out = 101)
    for (ep in list(NULL, numeric(0))) {
      expect_error(dy.dx(x, x^2, eval.point = ep),
                   "eval.point must be provided", fixed = TRUE)
    }
    expect_error(dy.dx(x, x^2), "eval.point must be provided", fixed = TRUE)
  })
  aud_case("dy.dx does not depend on the row order of x", {
    set.seed(17)
    x <- runif(300, -2, 2)
    y <- x^3
    o <- order(x)
    ep <- c(-1.5, 0.5, 1.5)
    expect_equal(dy.dx(x, y, ep), dy.dx(x[o], y[o], ep))
    expect_equal(dy.dx(x, y, ep)$first.derivative, 3 * ep^2, tolerance = 0.1)
  })
  aud_case("dy.dx works when small samples repeat a bandwidth", {
    for (n in c(10, 12)) {
      x <- seq(0, 3, length.out = n)
      d <- dy.dx(x, x^2, eval.point = 1.5)
      expect_true(all(is.finite(c(d$first.derivative, d$second.derivative))))
      expect_equal(d$first.derivative, 3, tolerance = 0.05, info = paste("n =", n))
    }
  })
  aud_case("dy.dx averages over all five bandwidth positions", env = dydx_duplicate_bands, {
    for (n in c(20, 30, 50, 10)) {
      x <- seq(0, 3, length.out = n)
      d <- dy.dx(x, x^2, eval.point = 1.5)
      expect_equal(c(d$first.derivative, d$second.derivative), band_means(x, x^2, 1.5),
                   tolerance = 1e-12, info = paste("n =", n))
    }
  })
})

test_that("DER02: dy.d_: stencils in support, discrete and constant conventions, multivariate boundaries, irrelevant regressor, scalar eval points, one-column frames, row thinning", {
  withr::local_seed(20261108)
  aud_case("local derivative stencils stay in support at boundaries and interior points", {
    x <- seq(-2, 3, length.out = 101)
    for (at in c(-2, -1.99, 0, 2.99, 3)) for (width in c(.03, .1, .25)) {
      s <- NNS:::.nns_dy_stencil(x, at, width)
      expect_true(all(s$nodes >= min(x) & s$nodes <= max(x)))
      # These are derivative identities for a polynomial, including unequal
      # spacing and one-sided endpoint formulas, not comparisons to the code.
      expect_equal(sum(s$first * (s$nodes^2 + 3 * s$nodes + 4)), 2 * at + 3,
                   tolerance = 1e-9)
      expect_equal(sum(s$second * (s$nodes^2 + 3 * s$nodes + 4)), 2,
                   tolerance = 1e-8)
    }
  })
  aud_case("discrete and constant coordinates have explicit derivative conventions", {
    two <- NNS:::.nns_dy_stencil(c(0, 0, 1, 1), .5, .1)
    expect_equal(sum(two$first * c(1, 4)), 3)
    expect_true(all(is.na(two$second)))
    constant <- NNS:::.nns_dy_stencil(rep(2, 10), 2, .1)
    expect_true(is.na(constant$first))
    expect_error(NNS:::.nns_dy_stencil(1:5, 6, .1), "inside the observed")
  })
  aud_case("multivariate derivative boundaries remain finite and reject extrapolation", {
    x <- as.matrix(expand.grid(a = seq(-1, 1, length.out = 15),
                              b = seq(-1, 1, length.out = 15)))
    y <- 2 * x[, 1] + x[, 2]^2
    z <- dy.d_(x, y, wrt = 1:2, eval.points = rbind(c(-1, 0), c(1, 0)), messages = FALSE)
    expect_true(all(is.finite(unlist(z))))
    expect_error(dy.d_(x, y, 1, eval.points = matrix(c(0, 2), 1), messages = FALSE),
                 "inside the observed")
  })
  aud_case("dy.d_ gives an irrelevant balanced regressor a slope near zero", {
    # A complete Cartesian design avoids incidental sample association and
    # sparse-neighbour variation at the evaluated tuples.
    x <- as.matrix(expand.grid(seq(-2, 2, length.out = 21),
                              seq(-2, 2, length.out = 21)))
    ep <- rbind(c(-1, 0), c(0, 1), c(1, -1))
    r <- dy.d_(x, 2 * x[, 1], wrt = 1:2, eval.points = ep, messages = FALSE)
    expect_true(all(r["First", 1][[1]] > 1))
    expect_true(all(abs(r["First", 2][[1]]) < 0.1))
  })
  aud_case("dy.d_ averages a scalar eval point over the observed rows (#65)", {
    # With three regressors the old grid set the other two to the same
    # quantile, so the average of d(x1 x2 x3)/dx1 = x2 x3 came out positive.
    set.seed(611); n <- 600
    x <- cbind(runif(n, -1, 1), runif(n, -1, 1), runif(n, -1, 1)); y <- x[, 1] * x[, 2] * x[, 3]
    ame <- mean(x[, 2] * x[, 3])
    r <- dy.d_(x, y, wrt = 1, eval.points = 0, messages = FALSE)
    expect_lt(abs(unlist(r["First", ])[1] - ame), 0.05)
  })
  aud_case("dy.d_ thins rows without a fixed stride", {
    r <- NNS:::.nns_spread_rows(20000, 10000)
    expect_length(r, 10000)
    expect_false(is.unsorted(r)); expect_false(anyDuplicated(r) > 0)
    # Every residue class is visited evenly, throughout the data.
    for (q in 2:5) expect_lt(max(abs(tabulate(r %% q + 1, q) - 10000 / q)), 0.01 * 10000 / q)
    expect_lt(abs(mean(r[r <= 10000] %% 2) - 0.5), 0.01)
    expect_identical(NNS:::.nns_spread_rows(50, 80), 1:50)
  })
  aud_case("dy.d_ accepts a one-column data frame of evaluation points (#150)", {
    set.seed(1)
    x <- cbind(a = runif(200, -1, 1), b = runif(200, -1, 1))
    y <- x[, 1]^2 + x[, 2]
    a <- dy.d_(x, y, 1, eval.points = data.frame(a = c(-.5, .5)), messages = FALSE)
    b <- dy.d_(x, y, 1, eval.points = c(-.5, .5), messages = FALSE)
    expect_equal(a, b)
  })
})

test_that("DER03: NNS.diff converges to the derivative and ignores an added constant", {
  withr::local_seed(20261108)
  aud_case("NNS.diff converges to the derivative within 1e-7", {
    cases <- list(list(log, 5, 0.2),
                  list(tanh, 0.5, 1 - tanh(0.5)^2),
                  list(function(x) x^4 - 3 * x, -2, -35))
    for (cs in cases) {
      r <- NNS.diff(cs[[1]], cs[[2]], digits = 16)
      expect_equal(unname(r["converged", 1]), 1)
      expect_lt(abs(r["DERIVATIVE", 1] - cs[[3]]), 1e-7)
    }
  })
  aud_case("NNS.diff is unaffected by a constant added to f", {
    d <- function(k, ...) NNS.diff(function(x) sin(x) + k, 1, digits = 16, ...)
    for (k in c(0, 1000)) {
      r <- d(k)
      expect_equal(unname(r["converged", 1]), 1)
      expect_lt(abs(r["DERIVATIVE", 1] - cos(1)), 1e-6)
    }
    expect_equal(d(1000, tol = 1e-4)["DERIVATIVE", 1],
                 d(0, tol = 1e-4)["DERIVATIVE", 1], tolerance = 1e-8)
  })
})
