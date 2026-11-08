# LPM.VaR and UPM.VaR. Lean: nns-proofs VaRRoot.lean.

# Lean: LPMratio_strictMonoOn, UPMratio_strictAntiOn | BUGS: #15
test_that("V01: LPM.VaR and UPM.VaR invert their ratios on real and tied data, integer and fractional degrees", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_strictMonoOn (VaR root): LPM.ratio(d, LPM.VaR(p, d, x), x) = p, d = 1..4", {
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      v <- LPM.VaR(aud_var_p, d, x)
      aud_close(LPM.ratio(d, v, x), aud_var_p, aud_ratio_tol(x), info = paste(nm, d))
    }
  })
  aud_case("UPMratio_strictAntiOn (VaR root): UPM.ratio(d, UPM.VaR(p, d, x), x) = p, d = 1..4", {
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      v <- UPM.VaR(aud_var_p, d, x)
      aud_close(UPM.ratio(d, v, x), aud_var_p, aud_ratio_tol(x), info = paste(nm, d))
    }
  })
  aud_case("LPMratio_strictMonoOn (VaR root): fallback path d = 2.5, 5 within the uniroot tolerance", {
    # Fallback solves to sqrt(eps) * range in t (BUGS.md #15, design choice).
    for (nm in names(fr$kept)) for (d in c(2.5, 5)) {
      x <- fr$kept[[nm]]
      v <- LPM.VaR(aud_var_p, d, x)
      aud_close(LPM.ratio(d, v, x), aud_var_p, 1e-6, info = paste(nm, d))
    }
  })
  aud_case("LPMratio_strictMonoOn (VaR root): exact dyadic probe with ties, z = (1,1,1,2,5,5)", {
    z <- c(1, 1, 1, 2, 5, 5)
    v <- LPM.VaR(c(0.1, 0.5, 0.9), 1, z)
    aud_close(v, c(1.3, 2.5, 4.318181818181818), 4 * aud_eps * 5)
    aud_close(LPM.ratio(1, v, z), c(0.1, 0.5, 0.9), 4 * aud_eps)
  })
})

# Lean: LPMratio_monotone, LPMratio_eq_one_of_forall_le, LPMratio_eq_zero_of_forall_ge, LPMratio_of_forall_eq | BUGS: none
test_that("V02: VaR is monotone in p, gives the sample ends at p = 0 and 1 and the value of a constant sample", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("LPMratio_monotone (VaR): LPM.VaR is nondecreasing in p, UPM.VaR nonincreasing", {
    p <- seq(0, 1, by = 0.05)
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      expect_true(all(diff(LPM.VaR(p, d, x)) >= 0), info = paste(nm, d))
      expect_true(all(diff(UPM.VaR(p, d, x)) <= 0), info = paste(nm, d))
    }
  })
  aud_case("LPMratio_eq_zero_of_forall_ge / LPMratio_eq_one_of_forall_le (VaR ends): p = 0 -> min, p = 1 -> max", {
    for (nm in names(fr$kept)) for (d in 1:4) {
      x <- fr$kept[[nm]]
      expect_identical(LPM.VaR(c(0, 1), d, x), c(min(x), max(x)), info = paste(nm, d))
      expect_identical(UPM.VaR(c(1, 0), d, x), c(min(x), max(x)), info = paste(nm, d))
    }
  })
  aud_case("LPMratio_of_forall_eq (VaR): constant sample and n = 1 return the common value", {
    for (d in c(0, 1, 2, 2.5)) {
      expect_identical(unname(LPM.VaR(c(0.3, 0.6), d, rep(4, 3))), c(4, 4), info = d)
      expect_identical(unname(LPM.VaR(0.3, d, 7)), 7, info = d)
    }
  })
})

# Lean: UPMVaR0_eq_LPMVaR0, quantile7_eq, floor_add_one_le_card_le_quantile7, card_lt_quantile7_le_ceil, LPMratio_zero_quantile7_gt, LPMVaR0_counterexample, LPMVaR0_not_minimal | BUGS: none
test_that("V03: degree-0 VaR is the type-7 sample quantile, not the ECDF inverse", {
  withr::local_seed(20261108)
  fr <- aud_fx_ratio()
  aud_case("quantile7_eq / UPMVaR0_eq_LPMVaR0 (VaRRoot.lean): degree-0 VaR is the type-7 sample quantile", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      expect_identical(LPM.VaR(aud_var_p, 0, x), quantile(x, aud_var_p, type = 7), info = nm)
      expect_identical(unname(UPM.VaR(aud_var_p, 0, x)), unname(quantile(x, 1 - aud_var_p, type = 7)), info = nm)
    }
  })
  aud_case("floor_add_one_le_card_le_quantile7 / card_lt_quantile7_le_ceil (d = 0): at p = k/(n-1) VaR is the minimal t with ECDF(t) >= p", {
    # At p = k/(n-1) the type-7 index (n-1)p = k is an integer: floor_add_one_le_card_le_quantile7
    # gives ECDF(VaR) >= (k+1)/n >= p and card_lt_quantile7_le_ceil gives #{x < VaR} <= k, so the
    # ECDF just below VaR is <= k/n < p: VaR is the minimal t with ECDF(t) >= p there.
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      n <- length(x)
      p <- unique(round(seq(1, n - 2, length.out = 9))) / (n - 1)
      p <- p[(n - 1) * p == round((n - 1) * p)]  # type-7 index exactly an integer in double
      v <- unname(LPM.VaR(p, 0, x))
      below <- vapply(v, function(s) max(c(min(x) - 1, x[x < s])), numeric(1))
      expect_true(all(LPM.ratio(0, v, x) >= p), info = nm)
      expect_true(all(LPM.ratio(0, below, x) < p), info = nm)
    }
  })
  aud_case("LPMratio_zero_quantile7_gt / floor_add_one_le_card_le_quantile7 / card_lt_quantile7_le_ceil (d = 0): type-7 VaR count bounds", {
    for (nm in names(fr$kept)) {
      x <- fr$kept[[nm]]
      n <- length(x)
      v <- unname(LPM.VaR(aud_var_p, 0, x))
      h <- (n - 1) * aud_var_p  # the type-7 index, computed in double exactly as quantile() does
      n_le <- vapply(v, function(s) sum(x <= s), numeric(1))
      n_lt <- vapply(v, function(s) sum(x < s), numeric(1))
      expect_true(all(n_le >= floor(h) + 1), info = nm)                  # floor_add_one_le_card_le_quantile7
      expect_true(all(LPM.ratio(0, v, x) > aud_var_p - aud_var_p / n), info = nm) # LPMratio_zero_quantile7_gt
      expect_true(all(n_lt <= ceiling(h)), info = nm)                      # card_lt_quantile7_le_ceil
    }
  })
  aud_case("LPMVaR0_counterexample / LPMVaR0_not_minimal: type 7 is not the ECDF inverse on (0, 1, 2)", {
    x <- c(0, 1, 2)
    v <- unname(LPM.VaR(c(0.1, 0.9), 0, x))
    aud_close(v, c(0.2, 1.8), 4 * aud_eps)
    expect_identical(LPM.ratio(0, v[2], x), 2 / 3)   # 2/3 < 0.9: ECDF(VaR) >= p fails
    expect_identical(LPM.ratio(0, 0, x), 1 / 3)      # 0 < VaR(0.1) yet ECDF(0) >= 0.1
  })
})

# Lean: none | BUGS: #12
test_that("V04: BUGS #12: the VaR solvers stay accurate at extreme degree and data scale", {
  withr::local_seed(20261108)
  aud_case("the non-integer-degree VaR solver finds the symmetric median at degree 3e9 without warnings", {
    expect_no_warning(v <- LPM.VaR(0.5, 3e9, c(0, 1)))
    expect_equal(v, 0.5, tolerance = 1e-6)
  })
  aud_case("UPM.VaR at degree 1e4 finds the symmetric median", {
    expect_equal(UPM.VaR(0.5, 1e4, c(0, 1)), 0.5, tolerance = 1e-6)
  })
  aud_case("the integer-degree VaR solver is scale invariant on data of scale 1e-200", {
    x <- c(-3, 1, 2, 2)
    p <- c(0.1, 0.5, 0.9)
    expect_equal(LPM.VaR(p, 2, 1e-200 * x) * 1e200, LPM.VaR(p, 2, x), tolerance = 1e-12)
    expect_equal(UPM.VaR(p, 4, 1e200 * x) * 1e-200, UPM.VaR(p, 4, x), tolerance = 1e-12)
  })
  aud_case("the non-integer-degree VaR solver is scale invariant on data of scale 1e-200", {
    x <- c(-3, 1, 2, 2)
    p <- c(0.1, 0.5, 0.9)
    expect_equal(LPM.VaR(p, 2.5, 1e-200 * x) * 1e200, LPM.VaR(p, 2.5, x), tolerance = 1e-6)
  })
})

# Lean: intLPMVaR_snap_counterexample, intUPMVaR_snap_counterexample | BUGS: #29, tie-snap, 09b94e6
# A tolerance on the ratio is not a tolerance on t where the ratio is flat, so the
# solver used to return a break whose ratio was within 2^-39 of p.
test_that("V05: BUGS #29: integer-degree VaR is the root in t, not a break near it", {
  withr::local_seed(20261108)
  aud_case("integer-degree VaR on c(0, 1) is the closed-form root, also in the far tails", {
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
  aud_case("integer-degree VaR resolves a percentile just off a break's ratio", {
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
  aud_case("integer-degree VaR does not return a break for a root resolved off it", {
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
  aud_case("integer-degree VaR resolves a root a few ulps off a break far from 0", {
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
})

# Lean: none | BUGS: #8, #10, #29, 2887cd4, 4217b25
test_that("V06: BUGS #8, #10: integer-degree VaR at tiny percentiles, far clusters and overflowing spreads", {
  withr::local_seed(20261108)
  aud_case("integer-degree VaR recovers the order statistic of a far cluster at tiny ratios", {
    # The ratios here are 1e-28 to 1e-17: the former median-centred power sums
    # cancelled and returned the wrong order statistic.
    withr::local_seed(11)
    x <- c(rnorm(40, -1e5, 1), rnorm(61, 0, 1))
    target <- sort(x)[c(5, 20, 35)]
    for (d in 3:4) {
      expect_lte(max(abs(LPM.VaR(LPM.ratio(d, target, x), d, x) - target)), 64 * .Machine$double.eps * 1e5)
      expect_lte(max(abs(UPM.VaR(UPM.ratio(d, -target, -x), d, -x) + target)), 64 * .Machine$double.eps * 1e5)
    }
  })
  aud_case("integer-degree VaR ratio residual at tiny percentiles is the rounding of t", {
    # Near min(x) the LPM ratio's relative slope is about d / (t - min(x)), so a
    # threshold correct to a few ulps has relative ratio residual of at most
    # about d * eps * |t| / (t - min(x)) (mirror image for UPM near max(x)).
    withr::local_seed(29)
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
  aud_case("integer-degree VaR at the smallest normal percentile converges in few steps", {
    # The gap [min x, next order statistic] is hundreds of binades wider than the
    # root offset here; bracketing it by the whole half gap took zeroin about
    # 1000 iterations at d = 2. With min(x) = 0 exactly, the offset is
    # representable: on the first gap LPM = s^d / n and UPM = mean((x - s)^d),
    # which is mean(x^d) to relative 1e-60, so s^d = p * n * mean(x^d) / (1 - p).
    # Its d-th root is refined by one Newton step (pow() with the inexact
    # exponent 1/3 is off by |log(y)| * eps / 3, about 60 eps here).
    withr::local_seed(5)
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
  aud_case("integer-degree VaR solves two points whose spread overflows", {
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
})

# Lean: none | BUGS: tie-snap, 9c65004, d7059dd
# The breaks are tested by the sign of h = (1 - p) LPM - p UPM. On tied decimal data
# whose decimal ratio at an order statistic is exactly p, rounding leaves h a few eps
# off 0 at the break, and an exact sign test returned the order statistic +/- 1 ulp
# (NNS.boost then dropped every learner tied at the cutoff). These cases pass on the
# pre-audit release too; V05 holds the tie-snap failures there.
test_that("V07: tie-snap: integer-degree VaR returns an order statistic whose ratio is p, including 0", {
  withr::local_seed(20261108)
  aud_case("integer-degree VaR returns a tied decimal order statistic whose ratio is p", {
    s <- c(0.1, 0.3, 0.5, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.9)
    # Decimal LPM and UPM at 0.7 are 1.2 and 0.3: the ratio is 0.8.
    expect_identical(LPM.VaR(0.8, 1, s), 0.7)
    expect_identical(UPM.VaR(0.2, 1, s), 0.7)
    expect_identical(UPM.VaR(0.8, 1, -s), -0.7)
    expect_identical(LPM.VaR(0.2, 1, -s), -0.7)
    # NNS.boost keeps the learners at or above LPM.VaR(0.8, 1, scores).
    expect_identical(sum(s >= LPM.VaR(0.8, 1, s)), 7L)
  })
  aud_case("integer-degree VaR returns the order statistic wherever the decimal ratio is p", {
    # Scores k / den + offset with integer k: the decimal partial moments at a
    # tie are integers, so "the decimal ratio is num / pden" is an exact integer
    # test. The double data and percentile (1 / 3 included) are within rounding
    # of those values, and so is the root. Two or three levels among the first 13
    # multiples of 1 / den, repeated 3 to 8 times, make exact ratios common:
    # with this seed there are 354 checks (169, 90, 45 and 50 at d = 1 to 4).
    probs <- list(c(4, 5), c(1, 5), c(1, 2), c(3, 4), c(1, 4), c(2, 5), c(9, 10),
                  c(1, 10), c(1, 3), c(2, 3))
    withr::local_seed(123)
    checks <- 0L
    for (i in 1:5000) {
      den <- sample(c(10, 100, 1000), 1)
      offset <- sample(c(0, 0, 1, -3, 17), 1)
      k <- sample(sample(0:min(den, 12), sample(2:3, 1)), sample(3:8, 1), replace = TRUE)
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
    expect_gt(checks, 200L)
  })
  aud_case("the order statistic 0 is returned as 0 (kills M42: the snap allowance at a break)", {
    # pin: tie-snap (9c65004) at an order statistic equal to 0, no ties; an allowance
    # relative to |x_k| is 0 there, so the break must be matched exactly.
    x <- c(0, 0.8, 0.2, -0.5, -0.2)
    expect_identical(LPM.VaR(LPM.ratio(4, 0, x), 4, x), 0)
    expect_identical(UPM.VaR(UPM.ratio(4, 0, x), 4, x), 0)
    expect_identical(LPM.VaR(LPM.ratio(1, 0, x), 1, x), 0)
    expect_identical(UPM.VaR(UPM.ratio(1, 0, -x), 1, -x), 0)
  })
})
