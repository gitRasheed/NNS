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

test_that("integer-degree VaR at the smallest normal percentile converges in few steps", {
  # The gap [min x, next order statistic] is hundreds of binades wider than the
  # root offset here; bracketing it by the whole half gap took zeroin about
  # 1000 iterations at d = 2. With min(x) = 0 exactly, the offset is
  # representable: on the first gap LPM = s^d / n and UPM = mean((x - s)^d),
  # which is mean(x^d) to relative 1e-60, so s^d = p * n * mean(x^d) / (1 - p).
  # Its d-th root is refined by one Newton step (pow() with the inexact
  # exponent 1/3 is off by |log(y)| * eps / 3, about 60 eps here).
  set.seed(5)
  x <- rnorm(1e5)
  x <- x - min(x)
  p <- .Machine$double.xmin
  e <- new.env()
  e$iter <- integer(0)
  suppressMessages(trace(
    "uniroot",
    exit = bquote(assign("iter", c(get("iter", envir = .(e)), iter), envir = .(e))),
    where = asNamespace("stats"), print = FALSE
  ))
  on.exit(suppressMessages(untrace("uniroot", where = asNamespace("stats"))), add = TRUE)
  for (d in 1:4) {
    y <- p * length(x) * mean(x^d) / (1 - p)
    root <- y^(1 / d)
    root <- root - (root^d - y) / (d * root^(d - 1))
    expect_no_warning(v <- LPM.VaR(p, d, x))
    expect_no_warning(w <- UPM.VaR(p, d, -x))
    expect_lt(abs(v / root - 1), 8 * .Machine$double.eps)
    expect_lt(abs(w / -root - 1), 8 * .Machine$double.eps)
  }
  expect_length(e$iter, 8L)
  # Bisection through the whole half gap takes about 1000 iterations at d = 2.
  expect_lt(max(e$iter), 100L)
})

test_that("integer-degree VaR solves two points whose spread overflows", {
  # z = x - median was c(0, Inf) here and the solve stopped with "missing value
  # where TRUE/FALSE needed". On c(a, b) the LPM root is a + (b - a) r / (1 + r)
  # with r = (p / (1 - p))^(1 / d); for b = -a = 1e308 that is 1e308 (2 s - 1),
  # s = r / (1 + r), and the UPM root is its negative.
  x <- c(-1e308, 1e308)
  p <- c(0.1, 0.5, 0.9)
  for (d in 1:4) {
    r <- (p / (1 - p))^(1 / d)
    root <- 1e308 * (2 * r / (1 + r) - 1)
    v <- LPM.VaR(p, d, x)
    w <- UPM.VaR(p, d, x)
    expect_true(all(is.finite(v)) && all(is.finite(w)))
    expect_lt(max(abs(v - root)), 1e308 * 1e-14)
    expect_lt(max(abs(w + root)), 1e308 * 1e-14)
  }
})

# --- A root within rounding of an order statistic is that order statistic -----
# The breaks are tested by the sign of h = (1 - p) LPM - p UPM. On tied decimal
# data whose decimal ratio at an order statistic is exactly p, rounding leaves
# h a few eps off 0 at the break, and an exact sign test returned the order
# statistic +/- 1 ulp (NNS.boost then dropped every learner tied at the cutoff).

test_that("integer-degree VaR returns a tied decimal order statistic whose ratio is p", {
  s <- c(0.1, 0.3, 0.5, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.9)
  # Decimal LPM and UPM at 0.7 are 1.2 and 0.3: the ratio is 0.8.
  expect_identical(LPM.VaR(0.8, 1, s), 0.7)
  expect_identical(UPM.VaR(0.2, 1, s), 0.7)
  expect_identical(UPM.VaR(0.8, 1, -s), -0.7)
  expect_identical(LPM.VaR(0.2, 1, -s), -0.7)
  # NNS.boost keeps the learners at or above LPM.VaR(0.8, 1, scores).
  expect_identical(sum(s >= LPM.VaR(0.8, 1, s)), 7L)
})

test_that("integer-degree VaR returns the order statistic wherever the decimal ratio is p", {
  # Scores k / den + offset with integer k: the decimal partial moments at a
  # tie are integers, so "the decimal ratio is num / pden" is an exact integer
  # test. The double data and percentile are within rounding of those
  # decimals, and so is the root.
  probs <- list(c(4, 5), c(1, 5), c(1, 2), c(3, 4), c(1, 4), c(2, 5), c(9, 10), c(1, 10))
  set.seed(123)
  checks <- 0L
  for (i in 1:5000) {
    den <- sample(c(10, 100, 1000), 1)
    offset <- sample(c(0, 0, 1, -3, 17), 1)
    k <- sample(sample(0:den, sample(2:4, 1)), sample(4:40, 1), replace = TRUE)
    if (length(unique(k)) < 2L) next
    s <- k / den + offset
    d <- sample(1:4, 1)
    for (kk in unique(k)) {
      L <- sum(pmax(kk - k, 0)^d)
      U <- sum(pmax(k - kk, 0)^d)
      if (L == 0 || U == 0) next
      for (pr in probs) if (L * pr[2] == pr[1] * (L + U)) {
        checks <- checks + 1L
        t <- kk / den + offset
        p <- pr[1] / pr[2]
        q <- (pr[2] - pr[1]) / pr[2]
        info <- sprintf("d = %d, p = %d/%d, t = %s, s = %s", d, pr[1], pr[2], t,
                        paste(s, collapse = ", "))
        expect_identical(LPM.VaR(p, d, s), t, info = info)
        expect_identical(UPM.VaR(q, d, s), t, info = info)
        expect_identical(UPM.VaR(p, d, -s), -t, info = info)
      }
    }
  }
  expect_gt(checks, 50L)
})

test_that("integer-degree VaR does not return a break for a root resolved off it", {
  # The rounding allowance at a break is far below these offsets (#29).
  expect_gt(LPM.VaR(1e-12, 4, c(0, 1)), 9.99e-4)
  expect_lt(UPM.VaR(1e-12, 4, c(0, 1)), 1 - 9.99e-4)
  # Root 1 + 1.33e-13 (about 600 eps from 1) at d = 1, as at d = 2 above.
  x <- c(0, 1, 2)
  p <- 0.5 + 1e-13
  expect_gt(LPM.VaR(p, 1, x), 1)
  expect_lt(abs(LPM.VaR(p, 1, x) - (1 + 4e-13 / 3)), 8 * .Machine$double.eps)
  expect_lt(UPM.VaR(p, 1, x), 1)
})

test_that("integer-degree VaR resolves a root a few ulps off a break far from 0", {
  # The order-statistic allowance is 2 eps |x_k|, at most a quarter of the
  # half gap; an allowance of 8 eps |x_k| returned these breaks.
  # Root 1e6 + 1 + 2e-9 * (4 / 3), 11 ulps above 1e6 + 1.
  v <- LPM.VaR(0.5 + 1e-9, 1, 1e6 + 0:2)
  expect_gt(v, 1e6 + 1)
  expect_lt(abs(v - (1e6 + 1 + 4e-9 / 3)), 2 * 2^-33)
  w <- UPM.VaR(0.5 + 1e-9, 1, -(1e6 + 0:2))
  expect_identical(w, -v)
  # Root 1e12 + 1.001, 8 ulps above 1e12 + 1.
  x <- 1e12 + 0:3
  t <- 1e12 + 1.001
  expect_lte(abs(LPM.VaR(LPM.ratio(1, t, x), 1, x) - t), 2^-13)
  # The ratio at 1e15 + 1.25 is exactly 0.375: (2 * 0.25 + 1) / 4.
  expect_identical(LPM.VaR(0.375, 1, 1e15 + 0:3), 1e15 + 1.25)
  expect_identical(UPM.VaR(0.375, 1, -(1e15 + 0:3)), -(1e15 + 1.25))
})
