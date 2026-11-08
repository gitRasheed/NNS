var_root_data <- function() {
  set.seed(42)
  list(
    symmetric = c(-.003, -.002, -.001, .001, .002, .003),
    skewed = c(rexp(40) - 1, 5),
    returns = rnorm(500, 0, .01)
  )
}

test_that("integer-degree LPM.VaR finds the symmetric median at small scale", {
  x <- c(-.003, -.002, -.001, .001, .002, .003)
  expect_equal(LPM.VaR(.5, 4, x), 0, tolerance = 1e-12)
  expect_equal(UPM.VaR(.5, 4, x), 0, tolerance = 1e-12)
})

test_that("LPM.VaR and UPM.VaR scale with the data", {
  p <- c(.01, .05, .5, .95)
  for (x in var_root_data()) for (d in 1:5) {
    spread <- diff(range(x))
    lpm <- LPM.VaR(p, d, x)
    upm <- UPM.VaR(p, d, x)
    for (k in c(1e-4, 1, 1e4)) {
      expect_lt(max(abs(LPM.VaR(p, d, k * x) - k * lpm)), 1e-10 * k * spread)
      expect_lt(max(abs(UPM.VaR(p, d, k * x) - k * upm)), 1e-10 * k * spread)
    }
  }
})

test_that("LPM.VaR and UPM.VaR invert their ratios at any scale", {
  p <- c(.01, .05, .5, .95)
  for (x in var_root_data()) for (d in 1:5) for (k in c(1e-4, 1, 1e4)) {
    y <- k * x
    lr <- vapply(p, function(q) LPM.ratio(d, LPM.VaR(q, d, y), y), numeric(1))
    ur <- vapply(p, function(q) UPM.ratio(d, UPM.VaR(q, d, y), y), numeric(1))
    # Degree 5 uses the general solver, which the fractional-degree tests hold to 1e-7.
    tol <- if (d <= 4) 1e-8 else 1e-7
    expect_lt(max(abs(lr - p)), tol)
    expect_lt(max(abs(ur - p)), tol)
  }
})

test_that("integer-degree VaR inverts the ratio where it is steep", {
  x <- c(rep(0, 99999), 1)
  for (d in 1:4) {
    expect_lt(abs(LPM.ratio(d, LPM.VaR(.036, d, x), x) - .036), 1e-10)
    expect_lt(abs(UPM.ratio(d, UPM.VaR(.964, d, x), x) - .964), 1e-10)
  }
})

# --- Integer-degree VaR is the root in t, not a break near it -----------------
# (VaRRoot.lean: intLPMVaR_snap_counterexample / intUPMVaR_snap_counterexample.)
# A tolerance on the ratio is not a tolerance on t where the ratio is flat, so
# the solver used to return a break whose ratio was within 2^-39 of p.

test_that("integer-degree VaR on c(0, 1) is the closed-form root, also in the far tails", {
  # On c(0, 1), LPM.ratio(d, t) = t^d / (t^d + (1 - t)^d), so the LPM root is
  # r / (1 + r) with r = (p / (1 - p))^(1 / d), and the UPM root is 1 minus it.
  # p = 2^-996 (1 - p rounds to 1) has the exact r = 2^(-996 / d); pow() with
  # an inexact exponent 1/3 would be off by |log(p)| * eps / 3 there.
  x <- c(0, 1)
  for (d in 1:4) for (p in c(2^-996, 1e-12, 1e-6, 0.3)) {
    r <- if (p == 2^-996) 2^(-996 / d) else (p / (1 - p))^(1 / d)
    root <- r / (1 + r)
    expect_no_warning(v <- LPM.VaR(p, d, x))
    expect_no_warning(w <- UPM.VaR(p, d, x))
    expect_lt(abs(v / root - 1), 1e-14)
    expect_lt(abs(w - (1 - root)), 2 * .Machine$double.eps)
  }
  # The former answers were the sample ends, whose ratios are 0 and 1.
  expect_gt(LPM.VaR(1e-12, 4, x), 9.99e-4)
  expect_lt(UPM.VaR(1e-12, 4, x), 1 - 9.99e-4)
})

test_that("integer-degree VaR resolves a percentile just off a break's ratio", {
  # On c(0, 1, 2) at d = 2 the ratio is 1/2 at t = 1; for t = 1 + s the root
  # solves (2 - 3p) s^2 + 2 s + (1 - 2p) = 0 (stable form of the positive root).
  x <- c(0, 1, 2)
  p <- 0.5 + 1e-13
  s <- -(1 - 2 * p) / (1 + sqrt(1 - (2 - 3 * p) * (1 - 2 * p)))
  expect_gt(s, 0)
  expect_lte(abs(LPM.VaR(p, 2, x) - (1 + s)), 4 * .Machine$double.eps)
  # UPM.ratio = q is LPM.ratio = 1 - q (exact here: 1 - q = p).
  expect_lte(abs(UPM.VaR(1 - p, 2, x) - (1 + s)), 4 * .Machine$double.eps)
})

test_that("integer-degree VaR recovers the order statistic of a far cluster at tiny ratios", {
  # The ratios here are 1e-28 to 1e-17: the former median-centred power sums
  # cancelled and returned the wrong order statistic.
  set.seed(11)
  x <- c(rnorm(40, -1e5, 1), rnorm(61, 0, 1))
  target <- sort(x)[c(5, 20, 35)]
  for (d in 3:4) {
    expect_lte(max(abs(LPM.VaR(LPM.ratio(d, target, x), d, x) - target)), 64 * .Machine$double.eps * 1e5)
    expect_lte(max(abs(UPM.VaR(UPM.ratio(d, -target, -x), d, -x) + target)), 64 * .Machine$double.eps * 1e5)
  }
})

test_that("integer-degree VaR ratio residual at tiny percentiles is the rounding of t", {
  # Near min(x) the LPM ratio's relative slope is about d / (t - min(x)), so a
  # threshold correct to a few ulps has relative ratio residual of at most
  # about d * eps * |t| / (t - min(x)) (mirror image for UPM near max(x)).
  set.seed(29)
  x <- rnorm(200)
  eps <- .Machine$double.eps
  for (d in 1:4) {
    p <- c(1e-10, 1e-13)
    v <- LPM.VaR(p, d, x)
    w <- UPM.VaR(p, d, x)
    expect_true(all(v > min(x)))
    expect_true(all(w < max(x)))
    expect_true(all(abs(LPM.ratio(d, v, x) / p - 1) <= 8 * d * eps * abs(v) / (v - min(x))), info = d)
    expect_true(all(abs(UPM.ratio(d, w, x) / p - 1) <= 8 * d * eps * abs(w) / (max(x) - w)), info = d)
  }
})
