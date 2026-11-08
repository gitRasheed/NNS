# Univariate LPM and UPM: identities, ECDF limit, reflection, translation, scaling,
# monotonicity, the prefix backend (> 32 targets) and degree handling. Lean:
# nns-proofs PartialMoments.lean, RealDegree.lean, Integration.lean, PrefixBackend.lean.
#
# On the direct scan (<= 32 targets) every integer degree up to 32 is the same
# left-fold product, so those loops run degree 0, 1 and their largest integer
# degree (U01's master identity also runs 2, 3 and 4, the central-moment
# statements); the prefix backend recombines d + 1 columns, so its loops keep every degree.

# Lean: UPM_add_neg_one_pow_mul_LPM, LPM_zero_add_UPM_zero, UPM_one_sub_LPM_one, LPM_add_UPM_of_even, LPM_two_add_UPM_two_mean, rLPM_add_rUPM, centralMoment_two_eq, centralMoment_three_eq, centralMoment_four_eq, UPM_one_mean_eq_LPM_one_mean | BUGS: none
test_that("U01: LPM and UPM decompose the central and absolute moments on dyadic, random and real data", {
  withr::local_seed(20261108)
  fu <- aud_fx_uni()
  aud_case("UPM_add_neg_one_pow_mul_LPM: master identity is exact on dyadic data, d = 0 to 4 and 6", {
    for (d in c(0, 1, 2, 3, 4, 6)) {
      lhs <- aud_upm_each(d, fu$t_exact, fu$x_exact) + (-1)^d * aud_lpm_each(d, fu$t_exact, fu$x_exact)
      rhs <- vapply(fu$t_exact, function(t) sum((fu$x_exact - t)^d) / 8, numeric(1))
      expect_identical(lhs, rhs, label = paste("degree", d))
    }
  })
  aud_case("UPM_add_neg_one_pow_mul_LPM: master identity on random and real data, d = 0 to 4 and 6", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      for (d in c(0, 1, 2, 3, 4, 6)) {
        lhs <- aud_upm_each(d, ts, x) + (-1)^d * aud_lpm_each(d, ts, x)
        rhs <- vapply(ts, function(t) sum((x - t)^d) / length(x), numeric(1))
        aud_within_max(lhs, rhs, aud_pm_tol_v(d, ts, x), label = sprintf("%s d=%d", nm, d))
      }
    }
  })
  aud_case("UPM_add_neg_one_pow_mul_LPM: master identity holds for n = 1", {
    for (d in c(0, 1, 6)) {
      for (t in c(-1, 0.5, 1, 2.5)) {
        expect_identical(UPM(d, t, 1) + (-1)^d * LPM(d, t, 1), (1 - t)^d)
      }
    }
  })
  aud_case("LPM_zero_add_UPM_zero: LPM(0) + UPM(0) is exactly 1, ties included", {
    expect_identical(aud_lpm_each(0, fu$t_exact, fu$x_exact) + aud_upm_each(0, fu$t_exact, fu$x_exact),
                     rep(1, length(fu$t_exact)))
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      expect_identical(aud_lpm_each(0, ts, x) + aud_upm_each(0, ts, x), rep(1, length(ts)),
                       label = nm)
    }
  })
  aud_case("LPM_zero_add_UPM_zero: degenerate samples (n = 1, constant) at the tie", {
    expect_identical(LPM(0, 5, 5) + UPM(0, 5, 5), 1)
    expect_identical(LPM(0, 5, rep(5, 7)), 1)
    expect_identical(UPM(0, 5, rep(5, 7)), 0)
  })
  aud_case("UPM_one_sub_LPM_one: UPM(1) - LPM(1) = mean - t, exact on dyadic data", {
    expect_identical(aud_upm_each(1, fu$t_exact, fu$x_exact) - aud_lpm_each(1, fu$t_exact, fu$x_exact),
                     sum(fu$x_exact) / 8 - fu$t_exact)
  })
  aud_case("UPM_one_sub_LPM_one: UPM(1) - LPM(1) = mean - t on random and real data", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      aud_within_max(aud_upm_each(1, ts, x) - aud_lpm_each(1, ts, x), sum(x) / length(x) - ts,
                     aud_pm_tol_v(1, ts, x), label = nm)
    }
  })
  aud_case("LPM_add_UPM_of_even: LPM(d) + UPM(d) is the d-th central moment for even d", {
    for (d in c(2, 6)) {
      expect_identical(aud_lpm_each(d, fu$t_exact, fu$x_exact) + aud_upm_each(d, fu$t_exact, fu$x_exact),
                       vapply(fu$t_exact, function(t) sum((fu$x_exact - t)^d) / 8, numeric(1)))
      for (nm in names(fu$kept)) {
        x <- fu$kept[[nm]]
        ts <- aud_targets(x)
        aud_within_max(aud_lpm_each(d, ts, x) + aud_upm_each(d, ts, x),
                       vapply(ts, function(t) sum((x - t)^d) / length(x), numeric(1)),
                       aud_pm_tol_v(d, ts, x), label = sprintf("%s d=%d", nm, d))
      }
    }
  })
  aud_case("LPM_two_add_UPM_two_mean: LPM(2, mean) + UPM(2, mean) is the population variance", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      n <- length(x)
      mu <- mean(x)
      expect_lte(abs(LPM(2, mu, x) + UPM(2, mu, x) - var(x) * (n - 1) / n),
                 aud_pm_tol(2, mu, x), label = nm)
    }
  })
  aud_case("rLPM_add_rUPM: LPM(d) + UPM(d) = mean |x - t|^d for every real d, dyadic data", {
    # |x - 1| takes values in {0, 0.25, 1, 2.25, 4}: squares of dyadics, so the
    # d = 0.5 and d = 1.5 powers are exact too; two observations tie at t = 1.
    x <- c(-3, -1.25, 0, 0.75, 1, 1, 2, 3.25)
    for (d in c(0, 0.5, 1, 1.5, 3)) {
      expect_identical(LPM(d, 1, x) + UPM(d, 1, x), sum(abs(x - 1)^d) / 8,
                       label = paste("degree", d))
    }
  })
  aud_case("rLPM_add_rUPM: LPM(d) + UPM(d) = mean |x - t|^d for real d on random and real data", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      for (d in c(0, 0.25, 0.5, 1, 1.5, 2.5, 3, 3.7)) {
        aud_within_max(aud_lpm_each(d, ts, x) + aud_upm_each(d, ts, x),
                       vapply(ts, function(t) sum(abs(x - t)^d) / length(x), numeric(1)),
                       aud_pm_tol_v(d, ts, x), label = sprintf("%s d=%g", nm, d))
      }
    }
  })
})

# Lean: LPM_zero_eq_ecdf, LPM_zero_monotone, LPM_one_monotone, rLPM_eq_strict, rLPM_continuous, rLPM_tendsto_zero', rLPM_tendsto_zero, rLPM_limit_zero_lt_LPM_zero, rUPM_tendsto_zero | BUGS: none
test_that("U02: degree-0 LPM is the right-continuous ECDF, and the d -> 0+ limit drops the tie share", {
  withr::local_seed(20261108)
  fu <- aud_fx_uni()
  aud_case("LPM_zero_eq_ecdf: LPM(0, t, x) is the right-continuous ECDF, scalar targets", {
    expect_identical(aud_lpm_each(0, fu$t_exact, fu$x_exact),
                     vapply(fu$t_exact, function(t) sum(fu$x_exact <= t) / 8, numeric(1)))
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- c(aud_targets(x), head(unique(x), 10))
      expect_identical(aud_lpm_each(0, ts, x), ecdf(x)(ts), label = nm)
    }
  })
  aud_case("LPM_zero_eq_ecdf: LPM(0) with more than 32 targets is the ECDF at every data value", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- rep_len(sort(unique(x)), 60)
      expect_identical(LPM(0, ts, x), ecdf(x)(ts), label = nm)
    }
  })
  aud_case("LPM_zero_monotone / LPM_one_monotone: degree 0 and 1 LPM nondecreasing in t", {
    for (x in c(list(fu$x_exact), fu$kept)) {
      ts <- sort(c(aud_targets(x), head(unique(x), 20)))
      expect_true(all(diff(aud_lpm_each(0, ts, x)) >= 0))
      expect_true(all(diff(aud_lpm_each(1, ts, x)) >= 0))
    }
  })
  aud_case("rLPM_eq_strict: for d > 0 a tie at the target contributes nothing to LPM", {
    x <- c(1, 2, 2, 3)
    expect_identical(LPM(0, 2, x), 0.75)
    for (d in c(1e-9, 0.5, 1, 2)) {
      expect_identical(LPM(d, 2, x), 0.25, label = paste("degree", d))
    }
  })
  aud_case("rLPM_continuous needs d > 0: degree-0 LPM jumps by the tie share at a data point", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      t <- x[1]
      h <- 1e-6 * diff(range(x))
      tie_share <- sum(x == t) / length(x)
      # k/n - j/n is (k - j)/n up to one rounding
      expect_gte(LPM(0, t, x) - LPM(0, t - h, x), tie_share - 2 * .Machine$double.eps,
                 label = nm)
    }
  })
  aud_case("rLPM_tendsto_zero': as d -> 0+, LPM tends to the strict ECDF #{x < t} / n", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      aud_within_max(aud_lpm_each(1e-9, ts, x),
                     vapply(ts, function(t) sum(x < t) / length(x), numeric(1)), 1e-6,
                     label = nm)
    }
  })
  aud_case("rLPM_tendsto_zero / rLPM_limit_zero_lt_LPM_zero: the d -> 0+ limit drops the tie share", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      t <- x[1]
      tie_share <- sum(x == t) / length(x)
      expect_gt(tie_share, 0)
      expect_lte(abs(LPM(0, t, x) - LPM(1e-9, t, x) - tie_share), 1e-6, label = nm)
    }
  })
  aud_case("rUPM_tendsto_zero: UPM is right-continuous in the degree at 0", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      aud_within_max(aud_upm_each(1e-9, ts, x), aud_upm_each(0, ts, x), 1e-6, label = nm)
    }
  })
  aud_case("degree-0 ties with excess_ret = TRUE", {
    expect_equal(LPM(0, 2, 1:3, excess_ret = TRUE), 2 / 3)
    expect_equal(LPM(0, 2, 1:3, TRUE) + UPM(0, 2, 1:3, TRUE), 1)
  })
})

# Lean: LPM_nonneg, UPM_nonneg, rLPM_natCast, rUPM_natCast, rLPM_eq_max, rUPM_eq_max, rLPM_const, rUPM_const | BUGS: none
test_that("U03: real-degree LPM and UPM are means of positive-part powers, also on constant samples", {
  withr::local_seed(20261108)
  fu <- aud_fx_uni()
  aud_case("LPM_nonneg / UPM_nonneg: partial moments are nonnegative at every degree", {
    for (x in fu$kept) {
      ts <- aud_targets(x)
      for (d in c(0, 0.3, 1, 1.5, 6)) {
        expect_true(all(aud_lpm_each(d, ts, x) >= 0))
        expect_true(all(aud_upm_each(d, ts, x) >= 0))
      }
    }
  })
  aud_case("rLPM_natCast / rUPM_natCast: natural-degree kernel equals the real-power definition", {
    # rLPM_natCast: the real-power definition rLPM (rLowerDev, (t - v)^d as
    # rpow) at a natural degree k equals the natural-power LPM. In the kernel the
    # natural degree goes through repeatMultiplication; the oracle below is the
    # real-power route (R's `^` with a double exponent calls C pow; d = 0 keeps
    # the tie, as rpow 0 0 = 1). k = 1, 2 are identical on both routes (pow(v, 1)
    # = v, R squares as v * v), so k >= 3 is where the check bites.
    rpow_lpm <- function(k, t, x) sum(ifelse(x <= t, (t - x)^(k + 0), 0)) / length(x)
    rpow_upm <- function(k, t, x) sum(ifelse(x > t, (x - t)^(k + 0), 0)) / length(x)
    for (nm in names(c(list(exact = fu$x_exact), fu$kept))) {
      x <- c(list(exact = fu$x_exact), fu$kept)[[nm]]
      ts <- c(aud_targets(x), head(x, 5))
      for (k in c(0, 1, 6)) {
        tol <- aud_pm_tol_v(k, ts, x)
        aud_within_max(aud_lpm_each(k, ts, x), vapply(ts, rpow_lpm, numeric(1), k = k, x = x), tol,
                       label = sprintf("LPM %s k=%d", nm, k))
        aud_within_max(aud_upm_each(k, ts, x), vapply(ts, rpow_upm, numeric(1), k = k, x = x), tol,
                       label = sprintf("UPM %s k=%d", nm, k))
      }
    }
  })
  aud_case("rLPM_eq_max / rUPM_eq_max: for d > 0, LPM and UPM are means of positive parts", {
    for (d in c(1, 5)) {
      expect_identical(aud_lpm_each(d, fu$t_exact, fu$x_exact),
                       vapply(fu$t_exact, function(t) sum(pmax(t - fu$x_exact, 0)^d) / 8, numeric(1)))
      expect_identical(aud_upm_each(d, fu$t_exact, fu$x_exact),
                       vapply(fu$t_exact, function(t) sum(pmax(fu$x_exact - t, 0)^d) / 8, numeric(1)))
    }
  })
  aud_case("rLPM_eq_max / rUPM_eq_max: positive-part form for real d > 0 on random and real data", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      for (d in c(1e-9, 0.5, 1.5, 2.5)) {
        aud_within_max(aud_lpm_each(d, ts, x),
                       vapply(ts, function(t) sum(pmax(t - x, 0)^d) / length(x), numeric(1)),
                       aud_pm_tol_v(d, ts, x), label = sprintf("LPM %s d=%g", nm, d))
        aud_within_max(aud_upm_each(d, ts, x),
                       vapply(ts, function(t) sum(pmax(x - t, 0)^d) / length(x), numeric(1)),
                       aud_pm_tol_v(d, ts, x), label = sprintf("UPM %s d=%g", nm, d))
      }
    }
  })
  aud_case("rLPM_const / rUPM_const: constant sample gives [c <= t](t - c)^d and [c < t](c - t)^d", {
    # 0 and 1 results are exact at every n (sums of zeros; n / n). A sum of n
    # copies of v divided by n is exact for n in {1, 2} (v + v = 2v, /2 exact)
    # but not in general for n >= 3, so there the (t - c)^d value is compared
    # to n * eps relative (pm_tol: the scale is exactly (t - c)^d).
    for (n in c(1, 2, 5, 8)) {
      x <- rep(2, n)
      expect_identical(LPM(0, 2, x), 1)
      expect_identical(UPM(0, 2, x), 0)
      for (d in c(0.5, 1, 3, 3.7)) {
        expect_identical(LPM(d, 2, x), 0)
        expect_identical(UPM(d, 2, x), 0)
        expect_identical(UPM(d, 3.5, x), 0)
        expect_identical(LPM(d, 0.5, x), 0)
        if (n <= 2) {
          expect_identical(LPM(d, 3.5, x), 1.5^d, label = sprintf("LPM n=%d d=%g", n, d))
          expect_identical(UPM(d, 0.5, x), 1.5^d, label = sprintf("UPM n=%d d=%g", n, d))
        } else {
          expect_lte(abs(LPM(d, 3.5, x) - 1.5^d), aud_pm_tol(d, 3.5, x),
                     label = sprintf("LPM n=%d d=%g", n, d))
          expect_lte(abs(UPM(d, 0.5, x) - 1.5^d), aud_pm_tol(d, 0.5, x),
                     label = sprintf("UPM n=%d d=%g", n, d))
        }
      }
    }
  })
})

# Lean: rUPM_eq_rLPM_neg, rLPM_add_const, rUPM_add_const, rLPM_mul_const, rUPM_mul_const | BUGS: none
test_that("U04: LPM and UPM follow reflection, translation and scaling of sample and target", {
  withr::local_seed(20261108)
  fu <- aud_fx_uni()
  aud_case("rUPM_eq_rLPM_neg: for d > 0, UPM(d, t, x) = LPM(d, -t, -x) (direct path, structural)", {
    # On the direct (<= 32 target) kernel this is a structural identity, not a
    # numerical check: UPM_C forms x_i - t and LPM_C on (-t, -x) forms
    # (-t) - (-x_i), the same IEEE operation in the same loop order. The
    # numerical content is in the > 32 target test below.
    for (d in c(0.5, 1, 1.5, 3)) {
      expect_identical(aud_upm_each(d, fu$t_exact, fu$x_exact), aud_lpm_each(d, -fu$t_exact, -fu$x_exact))
      for (x in fu$kept) {
        ts <- aud_targets(x)
        expect_identical(aud_upm_each(d, ts, x), aud_lpm_each(d, -ts, -x))
      }
    }
  })
  aud_case("rUPM_eq_rLPM_neg: for integer d > 0, UPM(d, t, x) = LPM(d, -t, -x) on the > 32 target path", {
    # The prefix backend computes LPM and UPM with different recurrences, so
    # here the reflection is a real numerical check. Only integer d in 1..32
    # take the prefix backend; non-integer d falls back to LPM_C / UPM_C and
    # would be the structural identity again.
    for (nm in names(c(list(exact = fu$x_exact), fu$kept))) {
      x <- c(list(exact = fu$x_exact), fu$kept)[[nm]]
      ts <- rep_len(c(aud_targets(x), head(x, 10)), 40)
      for (d in 1:6) {
        aud_within_max(UPM(d, ts, x), LPM(d, -ts, -x), aud_pm_tol_v(d, ts, x),
                       label = sprintf("%s d=%d", nm, d))
      }
    }
  })
  aud_case("rUPM_eq_rLPM_neg: reflection fails at d = 0 when the target is tied", {
    x <- c(1, 2, 2, 3)
    expect_identical(UPM(0, 2, x), 0.25)
    expect_identical(LPM(0, -2, -x), 0.75)
  })
  aud_case("rLPM_add_const / rUPM_add_const: translation of sample and target is exact on dyadic data", {
    # data on the 1/8 grid (random, with ties) shifted by 1024 stays exact
    withr::local_seed(7)
    x <- c(fu$x_exact, round(rnorm(200) * 8) / 8)
    ts <- c(fu$t_exact, x[9:20], round(mean(x) * 8) / 8)
    for (d in c(0, 0.5, 1, 1.5, 4)) {
      expect_identical(aud_lpm_each(d, ts + 1024, x + 1024), aud_lpm_each(d, ts, x))
      expect_identical(aud_upm_each(d, ts + 1024, x + 1024), aud_upm_each(d, ts, x))
    }
  })
  aud_case("rLPM_add_const / rUPM_add_const: translation is exact on integer-valued real data", {
    for (x in list(faithful$waiting, as.numeric(na.omit(airquality$Ozone)), mtcars$hp)) {
      ts <- c(median(x), x[1:5], quantile(x, c(0.25, 0.75), type = 1, names = FALSE))
      for (d in c(0, 0.5, 1, 3)) {
        expect_identical(aud_lpm_each(d, ts - 37, x - 37), aud_lpm_each(d, ts, x))
        expect_identical(aud_upm_each(d, ts - 37, x - 37), aud_upm_each(d, ts, x))
      }
    }
  })
  aud_case("rLPM_mul_const / rUPM_mul_const: scaling by 4 multiplies integer-degree moments by 4^d exactly", {
    for (x in c(list(fu$x_exact), fu$kept)) {
      ts <- aud_targets(x)
      for (d in c(0, 1, 5)) {
        expect_identical(aud_lpm_each(d, 4 * ts, 4 * x), 4^d * aud_lpm_each(d, ts, x))
        expect_identical(aud_upm_each(d, 4 * ts, 4 * x), 4^d * aud_upm_each(d, ts, x))
      }
    }
  })
  aud_case("rLPM_mul_const / rUPM_mul_const: scaling by a > 0 multiplies by a^d for real d", {
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- aud_targets(x)
      for (a in c(0.37, 4, 1000)) {
        for (d in c(0.5, 1.5, 2.5, 3)) {
          tol <- (d + 2) * a^d * aud_pm_tol_v(d, ts, x)
          aud_within_max(aud_lpm_each(d, a * ts, a * x), a^d * aud_lpm_each(d, ts, x), tol,
                         label = sprintf("LPM %s a=%g d=%g", nm, a, d))
          aud_within_max(aud_upm_each(d, a * ts, a * x), a^d * aud_upm_each(d, ts, x), tol,
                         label = sprintf("UPM %s a=%g d=%g", nm, a, d))
        }
      }
    }
  })
})

# Lean: LPM_monotone_target, UPM_antitone_target, rLPM_monotone, rUPM_antitone, rLPM_continuous, rUPM_continuous, rLPM_convexOn, rUPM_convexOn, not_rLPM_convexOn_of_lt_one | BUGS: none
test_that("U05: LPM is nondecreasing, continuous and (d >= 1) convex in the target; UPM mirrors it", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  fu <- aud_fx_uni()
  aud_case("LPM_monotone_target / UPM_antitone_target: monotone in the target at degree 0 on dyadic data", {
    expect_true(all(diff(LPM(0, fm$t_dy, fm$x_dy)) >= 0))
    expect_true(all(diff(UPM(0, fm$t_dy, fm$x_dy)) <= 0))
  })
  aud_case("rLPM_monotone: LPM is nondecreasing in t for real d >= 0 (direct path, exact)", {
    for (x in c(list(fu$x_exact), fu$kept)) {
      ts <- sort(c(aud_targets(x), head(unique(x), 20)))
      for (d in c(0, 0.3, 0.5, 1, 1.5, 3)) {
        expect_true(all(diff(aud_lpm_each(d, ts, x)) >= 0), label = paste("degree", d))
      }
    }
  })
  aud_case("rUPM_antitone: UPM is nonincreasing in t for real d >= 0 (direct path, exact)", {
    for (x in c(list(fu$x_exact), fu$kept)) {
      ts <- sort(c(aud_targets(x), head(unique(x), 20)))
      for (d in c(0, 0.3, 0.5, 1, 1.5, 3)) {
        expect_true(all(diff(aud_upm_each(d, ts, x)) <= 0), label = paste("degree", d))
      }
    }
  })
  aud_case("rLPM_monotone / rUPM_antitone (LPM_monotone_target / UPM_antitone_target): monotone in t on the > 32 target path", {
    for (x in fu$kept) {
      ts <- seq(min(x) - 1, max(x) + 1, length.out = 80)
      for (d in c(0, 0.5, 1, 2, 3)) {
        tol <- aud_pm_tol(d, max(ts), x) + aud_pm_tol(d, min(ts), x)
        expect_true(all(diff(LPM(d, ts, x)) >= -tol), label = paste("LPM degree", d))
        expect_true(all(diff(UPM(d, ts, x)) <= tol), label = paste("UPM degree", d))
      }
    }
  })
  aud_case("rLPM_continuous / rUPM_continuous: for d > 0 a perturbation h moves the moment by O(h)", {
    # |a^d - b^d| <= d |a - b| max(a, b)^(d - 1) for d >= 1, and <= |a - b|^d
    # for 0 < d < 1, applied per observation with a, b the positive parts.
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      h <- 1e-6 * diff(range(x))
      for (t in c(x[1], median(x))) {
        M <- max(abs(x - t)) + h
        for (d in c(0.25, 0.5, 1, 1.5, 3)) {
          bound <- if (d >= 1) d * h * M^(d - 1) else h^d
          bound <- bound + aud_pm_tol(d, t, x) + aud_pm_tol(d, t + h, x)
          expect_lte(abs(LPM(d, t + h, x) - LPM(d, t, x)), bound,
                     label = sprintf("LPM %s d=%g +h", nm, d))
          expect_lte(abs(LPM(d, t - h, x) - LPM(d, t, x)), bound,
                     label = sprintf("LPM %s d=%g -h", nm, d))
          expect_lte(abs(UPM(d, t + h, x) - UPM(d, t, x)), bound,
                     label = sprintf("UPM %s d=%g +h", nm, d))
          expect_lte(abs(UPM(d, t - h, x) - UPM(d, t, x)), bound,
                     label = sprintf("UPM %s d=%g -h", nm, d))
        }
      }
    }
  })
  aud_case("rLPM_convexOn / rUPM_convexOn: midpoint inequality in t for d >= 1", {
    for (nm in names(c(list(exact = fu$x_exact), fu$kept))) {
      x <- c(list(exact = fu$x_exact), fu$kept)[[nm]]
      ts <- aud_targets(x)
      pairs <- utils::combn(length(ts), 2)
      a <- ts[pairs[1, ]]
      b <- ts[pairs[2, ]]
      for (d in c(1, 1.5, 4)) {
        tol <- vapply(seq_along(a), function(i) aud_pm_tol(d, a[i], x) + aud_pm_tol(d, b[i], x),
                      numeric(1))
        expect_true(all(aud_lpm_each(d, (a + b) / 2, x) <=
                          (aud_lpm_each(d, a, x) + aud_lpm_each(d, b, x)) / 2 + tol),
                    label = sprintf("LPM %s d=%g", nm, d))
        expect_true(all(aud_upm_each(d, (a + b) / 2, x) <=
                          (aud_upm_each(d, a, x) + aud_upm_each(d, b, x)) / 2 + tol),
                    label = sprintf("UPM %s d=%g", nm, d))
      }
    }
  })
  aud_case("not_rLPM_convexOn_of_lt_one: for 0 < d < 1 the one-point LPM violates the midpoint inequality", {
    for (d in c(0.25, 0.5, 0.9)) {
      expect_gt(LPM(d, 0.5, 0), (LPM(d, 0, 0) + LPM(d, 1, 0)) / 2)
      expect_identical(LPM(d, 0.5, 0), 0.5^d)
    }
  })
})

# Lean: prefixLPM_eq_LPM, prefixUPM_eq_UPM, workerLPM_eq_rLPMVec, workerUPM_eq_rUPMVec, usesPrefix_iff, rLPMVec_sub_rLPM_tendsto, rLPMVec_sub_rLPM_tendsto_natCast, rLPMVecExact_eq_rLPM, rUPMVecExact_eq_rUPM, prefixDegreeOK_natCast_iff | BUGS: none
test_that("U06: more than 32 targets (prefix backend) give the direct-scan values at integer degrees", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  fu <- aud_fx_uni()
  aud_case("usesPrefix_iff: 32 targets stay on the direct kernel (identical to scalar calls), 33 switch to the prefix backend", {
    # NNS_DIRECT_PATH_MAX_TARGETS = 32: a 32-target call runs the same scalar
    # kernel per target, so it is bit-identical; at 33 targets the prefix backend
    # takes over (a probe shows ~1e-16 differences at d >= 1), held to the
    # backend's rounding bound.
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      for (nt in c(32L, 33L)) {
        ts <- seq(min(x) - 0.5, max(x) + 0.5, length.out = nt)
        ts[c(2, nt - 1)] <- c(x[1], median(x))   # include ties with observations
        for (d in 0:3) {
          lab <- sprintf("%s nt=%d d=%d", nm, nt, d)
          if (nt == 32L) {
            expect_identical(LPM(d, ts, x), aud_lpm_each(d, ts, x), label = paste(lab, "LPM"))
            expect_identical(UPM(d, ts, x), aud_upm_each(d, ts, x), label = paste(lab, "UPM"))
          } else {
            tol <- 8 * (length(x) + d) * aud_eps * max(aud_abs_mom(x, d, min(ts)), aud_abs_mom(x, d, max(ts)))
            aud_within_max(LPM(d, ts, x), aud_lpm_each(d, ts, x), tol, label = paste(lab, "LPM"))
            aud_within_max(UPM(d, ts, x), aud_upm_each(d, ts, x), tol, label = paste(lab, "UPM"))
          }
        }
      }
    }
  })
  aud_case("rLPMVecExact_eq_rLPM (prefixLPM_eq_LPM / prefixUPM_eq_UPM): degree 0 is bit-identical on both paths, ties included", {
    for (x in c(list(fu$x_exact), fu$kept)) {
      ts <- rep_len(c(aud_targets(x), head(x, 10)), 40)
      expect_identical(LPM(0, ts, x), aud_lpm_each(0, ts, x))
      expect_identical(UPM(0, ts, x), aud_upm_each(0, ts, x))
    }
  })
  aud_case("rLPMVecExact_eq_rLPM (workerLPM_eq_rLPMVec / workerUPM_eq_rUPMVec): integer degrees 1..6, 16, 32 agree to n * eps on both paths", {
    # d = 32 is the largest degree the prefix backend accepts
    # (prefixDegreeOK_natCast_iff: d <= 32); targets outside the data range hit
    # the k = 0 / k = n branches of the prefix lpm() / upm().
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- rep_len(c(aud_targets(x), min(x) - 1, max(x) + 1, head(x, 10)), 40)
      for (d in c(1:6, 16, 32)) {
        tol <- aud_pm_tol_v(d, ts, x)
        aud_within_max(LPM(d, ts, x), aud_lpm_each(d, ts, x), tol,
                       label = sprintf("LPM %s d=%d", nm, d))
        aud_within_max(UPM(d, ts, x), aud_upm_each(d, ts, x), tol,
                       label = sprintf("UPM %s d=%d", nm, d))
      }
    }
  })
  aud_case("rLPMVecExact_eq_rLPM: 32 vs 33 targets and degree > 32 give the direct value", {
    x <- c(1, 2, 2, 3)
    expect_identical(LPM(0, rep(2, 32), x), rep(0.75, 32))
    expect_identical(LPM(0, rep(2, 33), x), rep(0.75, 33))
    expect_identical(LPM(1, rep(2, 33), x), rep(0.25, 33))
    z <- fu$random_samples$normal
    ts <- rep_len(aud_targets(z), 40)
    expect_identical(LPM(33, ts, z), aud_lpm_each(33, ts, z))
  })
  aud_case("rLPMVecExact_eq_rLPM: n = 1 and constant samples on the > 32 target path", {
    # every gap is in {-1, 0, 1} and every power sum a small integer, so both
    # paths are exact and bit-identity is justified
    for (xx in list(2, rep(2, 7))) {
      ts <- rep_len(c(1, 2, 3), 40)
      for (d in 0:3) {
        expect_identical(LPM(d, ts, xx), aud_lpm_each(d, ts, xx),
                         label = sprintf("LPM n=%d d=%d", length(xx), d))
        expect_identical(UPM(d, ts, xx), aud_upm_each(d, ts, xx),
                         label = sprintf("UPM n=%d d=%d", length(xx), d))
      }
    }
  })
  native_kernel_fixes <- fx_native_kernel_fixes()
  aud_case("price and t(3) return series, degrees 0..4, 8, 12, > 32 targets, plus UPM(2) just below the maximum", env = native_kernel_fixes, {
    set.seed(7)
    price <- 400 * cumprod(1 + rnorm(1500, 4e-4, 0.012))
    ret <- rt(1000, df = 3) * 0.01
    for (x in list(price, ret)) {
      t <- c(seq(min(x) - 0.1 * diff(range(x)), max(x) + 0.1 * diff(range(x)), length.out = 80),
             sort(x)[seq(1, length(x), by = 25)])
      for (d in c(0, 1, 2, 3, 4, 8, 12)) {
        l <- LPM(d, t, x); u <- UPM(d, t, x)
        expect_true(all(l >= 0) && all(u >= 0))
        expect_equal(l, direct_lpm(d, t, x), tolerance = 1e-10)
        expect_equal(u, direct_upm(d, t, x), tolerance = 1e-10)
        r <- LPM.ratio(d, t, x)
        expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
      }
    }
    # Degree 2 just below the maximum used to come out negative.
    x <- rnorm(1e4, 100, 1)
    t <- sort(x)[(1e4 - 99):1e4] - 1e-9
    expect_true(all(UPM(2, t, x) >= 0))
  })
  native_kernel_fixes <- fx_native_kernel_fixes()
  aud_case("UPM and UPM.ratio relative accuracy just below the maximum (> 32 targets)", env = native_kernel_fixes, {
    set.seed(11)
    x <- rnorm(1e4, 100, 1)
    t <- sort(x)[(1e4 - 39):1e4] - 1e-6
    expect_equal(UPM(2, t, x) / direct_upm(2, t, x), rep(1, 40), tolerance = 1e-10)
    direct_ratio <- vapply(t, function(tt) UPM.ratio(2, tt, x), numeric(1))
    expect_equal(UPM.ratio(2, t, x) / direct_ratio, rep(1, 40), tolerance = 1e-10)
  })
})

# Lean: rLPMVecExact_eq_rLPM, rUPMVecExact_eq_rUPM | BUGS: #0
# The value for (d, t, x) must not depend on how many targets share the call. At
# 8898017 the prefix backend snapped degrees within 1e-12 of an integer (modelled by
# rLPMVec; rLPMVec_sub_rLPM_tendsto and rLPMVec_sub_rLPM_tendsto_natCast describe
# that bug, so no test asserts them).
test_that("U07: BUGS #0: non-integer and near-zero degrees use pow on both dispatch paths", {
  withr::local_seed(20261108)
  fu <- aud_fx_uni()
  aud_case("a degree just above 0 gives the same ratios on both dispatch paths", {
    x <- c(1, 2, 2, 3)
    expect_identical(LPM.ratio(1e-13, c(2, aud_pad), x)[1], LPM.ratio(1e-13, 2, x))
    expect_identical(UPM.ratio(1e-13, c(2, aud_pad), x)[1], UPM.ratio(1e-13, 2, x))
  })
  aud_case("near-integer degrees are evaluated with pow on every path, not rounded", {
    withr::local_seed(1)
    r <- round(rnorm(250, 0, 0.01), 4)
    t <- sort(r)
    expect_identical(LPM.ratio(1e-13, t, r),
                     vapply(t, function(tt) LPM.ratio(1e-13, tt, r), numeric(1)))
    expect_identical(LPM(2 + 1e-13, t, r),
                     vapply(t, function(tt) LPM(2 + 1e-13, tt, r), numeric(1)))
  })
  aud_case("NNS.CDF at a tiny positive degree is not the degree-0 ECDF on tied data", {
    withr::local_seed(1)
    r <- round(rnorm(250, 0, 0.01), 4)
    expect_false(identical(NNS.CDF(r, degree = 1e-13, plot = FALSE)$Function$CDF,
                           NNS.CDF(r, degree = 0, plot = FALSE)$Function$CDF))
  })
  aud_case("a near-integer UPM degree is not rounded on the prefix path", {
    x <- c(0, 0.5, 1, 2)
    expect_identical(UPM(2 + 1e-13, c(0.5, aud_pad), x)[1], UPM(2 + 1e-13, 0.5, x))
    expect_false(identical(UPM(2 + 1e-13, c(0.5, aud_pad), x)[1], UPM(2, 0.5, x)))
  })
  aud_case("rLPMVecExact_eq_rLPM: non-integer degrees are bit-identical on both paths", {
    for (x in c(list(fu$x_exact), fu$kept)) {
      ts <- rep_len(c(aud_targets(x), head(x, 10)), 40)
      for (d in c(0.5, 1.5, 2.5, 2 + 2e-12)) {
        expect_identical(LPM(d, ts, x), aud_lpm_each(d, ts, x), label = paste("LPM degree", d))
        expect_identical(UPM(d, ts, x), aud_upm_each(d, ts, x), label = paste("UPM degree", d))
      }
    }
  })
  aud_case("rLPMVecExact_eq_rLPM: degree 1e-13 at a tied target gives the same LPM on both paths", {
    # fails on 8898017: see BUGS.md #0 (prefix path snaps 1e-13 to 0 and counts the tie)
    x <- c(1, 2, 2, 3)
    expect_identical(LPM(1e-13, 2, x), 0.25)
    expect_identical(LPM(1e-13, rep(2, 33), x), rep(0.25, 33))
  })
  aud_case("rLPMVecExact_eq_rLPM: computed near-zero degree 0.1 + 0.2 - 0.3 gives the same LPM on both paths", {
    # fails on 8898017: see BUGS.md #0
    d <- 0.1 + 0.2 - 0.3
    for (x in c(list(c(1, 2, 2, 3)), fu$real_samples)) {
      ts <- rep_len(c(x[1], median(x)), 40)
      expect_identical(LPM(d, ts, x), aud_lpm_each(d, ts, x))
    }
  })
  aud_case("rLPMVecExact_eq_rLPM / rUPMVecExact_eq_rUPM: near-integer degrees agree to n * eps on both paths", {
    # fails on 8898017: see BUGS.md #0 (degrees within 1e-12 of an integer are
    # evaluated at the integer on the prefix path; gap ~1e-12, about 1e4 * tol)
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- rep_len(c(aud_targets(x), head(x, 10)), 40)
      for (d in c(1 + 1e-13, 2 - 5e-13, 2 + 5e-13, 3 + 9e-13)) {
        tol <- aud_pm_tol_v(d, ts, x)
        aud_within_max(LPM(d, ts, x), aud_lpm_each(d, ts, x), tol,
                       label = sprintf("LPM %s d=%.15g", nm, d))
        aud_within_max(UPM(d, ts, x), aud_upm_each(d, ts, x), tol,
                       label = sprintf("UPM %s d=%.15g", nm, d))
      }
    }
  })
  aud_case("rUPMVecExact_eq_rUPM: degree 1e-13 gives the same UPM on both paths", {
    # fails on 8898017: see BUGS.md #0 (prefix path evaluates UPM(0); direct
    # path gives mean((x - t)_+^1e-13), which differs by ~1e-13 * mean log gap)
    for (nm in names(fu$kept)) {
      x <- fu$kept[[nm]]
      ts <- rep_len(aud_targets(x), 40)
      aud_within_max(UPM(1e-13, ts, x), aud_upm_each(1e-13, ts, x), aud_pm_tol_v(1e-13, ts, x),
                     label = nm)
    }
  })
})

# Lean: LPM_one_eq_integral_LPM_zero, LPM_succ_eq_integral, LPM_two_eq_iterated_integral, UPM_succ_eq_integral, UPM_one_eq_integral_UPM_zero, mean_sub_eq_integral | BUGS: none
test_that("U08: LPM of degree d + 1 is (d + 1) times the integral of LPM of degree d, ties included", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("LPM_one_eq_integral_LPM_zero: LPM(1, t) = integral of the ECDF LPM(0, .), exact with ties", {
    k <- aud_knot_grid(fm$x_int)
    expect_identical(LPM(1, k, fm$x_int), aud_cum_step(k, function(s) LPM(0, s, fm$x_int)))
  })
  aud_case("LPM_succ_eq_integral d = 1: LPM(2, t) = 2 * integral of LPM(1, .), exact", {
    k <- aud_knot_grid(fm$x_int)
    expect_identical(LPM(2, k, fm$x_int), 2 * aud_cum_trap(k, function(s) LPM(1, s, fm$x_int)))
  })
  aud_case("LPM_succ_eq_integral d = 2: LPM(3, t) = 3 * integral of LPM(2, .), exact (Simpson)", {
    k <- aud_knot_grid(fm$x_int)
    expect_identical(LPM(3, k, fm$x_int), aud_cum_simpson_half(k, function(s) LPM(2, s, fm$x_int)))
  })
  aud_case("LPM_two_eq_iterated_integral: LPM(2, t) = 2 * double integral of the ECDF (factor 2)", {
    k <- aud_knot_grid(fm$x_int)
    inner <- aud_cum_step(k, function(s) LPM(0, s, fm$x_int))  # piecewise linear, exact at knots
    outer <- c(0, cumsum((head(inner, -1) + inner[-1]) / 2 * diff(k)))
    expect_identical(LPM(2, k, fm$x_int), 2 * outer)
  })
  aud_case("UPM_one_eq_integral_UPM_zero / UPM_succ_eq_integral: upper analogues, exact", {
    k <- aud_knot_grid(fm$x_int)
    rev_cum <- function(v) rev(cumsum(rev(v)))
    h <- diff(k)
    u0 <- UPM(0, head(k, -1), fm$x_int)          # strict upper side: constant on [k_i, k_i+1)
    expect_identical(UPM(1, k, fm$x_int), c(rev_cum(u0 * h), 0))
    u1 <- UPM(1, k, fm$x_int)
    expect_identical(UPM(2, k, fm$x_int), 2 * c(rev_cum((head(u1, -1) + u1[-1]) / 2 * h), 0))
  })
  aud_case("LPM_succ_eq_integral d = 0, 1, 2: random and real data via the >32-target prefix path", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      k <- aud_knot_grid(x, refine = TRUE)
      expect_gt(length(k), 32)
      scale <- max(abs(x - mean(x)))
      tol <- function(d) 8 * (length(x) + length(k)) * aud_eps * (scale + diff(range(k)))^(d + 1)
      aud_within_max(LPM(1, k, x), aud_cum_step(k, function(s) LPM(0, s, x)), tol(0), label = nm)
      aud_within_max(LPM(2, k, x), 2 * aud_cum_trap(k, function(s) LPM(1, s, x)), tol(1), label = nm)
      aud_within_max(LPM(3, k, x), aud_cum_simpson_half(k, function(s) LPM(2, s, x)), tol(2), label = nm)
    }
  })
  aud_case("mean_sub_eq_integral: mean - t = integral_t^b (1 - F) - integral_a^t F, exact", {
    k <- aud_knot_grid(fm$x_int)
    F_int <- aud_cum_step(k, function(s) LPM(0, s, fm$x_int))
    S_int <- rev(c(0, cumsum(rev((1 - LPM(0, head(k, -1), fm$x_int)) * diff(k)))))
    expect_identical(mean(fm$x_int) - k, S_int - F_int)
  })
})

# Lean: LPM_eq_zero_of_lt, UPM_eq_zero_of_le, LPM_pos_of_lt, LPM_one_of_forall_le | BUGS: none
test_that("U09: LPM vanishes below the sample, is positive above its minimum and is t - mean above its maximum", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("LPM_eq_zero_of_lt / UPM_eq_zero_of_le: tails vanish; degree 0 tie at max is UPM = 0", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      for (d in c(0, 1, 3)) {
        expect_identical(LPM(d, min(x) - 1, x), 0, label = nm)
        expect_identical(UPM(d, max(x), x), 0, label = nm)
      }
      expect_gt(LPM(0, min(x), x), 0)  # lower side counts the tie at min
    }
  })
  aud_case("LPM_pos_of_lt: strictly above some observation the LPM is positive, every degree", {
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      t <- min(x) + diff(range(x)) / 100
      for (d in c(0, 1, 3)) expect_gt(LPM(d, t, x), 0, label = sprintf("%s d=%d", nm, d))
    }
  })
  aud_case("LPM_one_of_forall_le: at or above every observation LPM(1, t) = t - mean", {
    expect_identical(LPM(1, c(7, 8, 100), fm$x_int), c(7, 8, 100) - 3)
    for (nm in names(fm$kept)) {
      x <- fm$kept[[nm]]
      ts <- max(x) + c(0, 0.5, 10)
      aud_within_max(LPM(1, ts, x), ts - mean(x), aud_mom_tol(x, 1, max(x)) + 4 * aud_eps * max(abs(ts)),
                     label = nm)
    }
  })
})

# Lean: none | BUGS: #1, #2, #17
test_that("U10: BUGS #1, #2, #17: integer degrees keep the left-fold product up to 32 and use pow beyond, at any size", {
  withr::local_seed(20261108)
  aud_case("finite degrees above the int range are accepted and use pow", {
    x <- c(1, 1.5, 2)
    y <- c(1, 2, 2, 3)
    expect_identical(LPM(3e9, 3, y), Inf)
    expect_identical(LPM(2^31, 2, x), 1 / 3)
    expect_identical(LPM(2^31, 1.2, x), 0)
    expect_identical(UPM(3e9, 0, x), Inf)
    expect_identical(LPM(3e9, 3, y, excess_ret = TRUE), Inf)
    expect_identical(Co.LPM(3e9, x, x, 2, 2), 1 / 3)
  })
  aud_case("integer degrees up to 32 keep the left-fold product bit for bit", {
    # Deliberate compatibility pin: earlier releases evaluated integer degrees as
    # the product 1 * v * v * ..., and results for degrees <= 32 must not move.
    expect_identical(LPM(32, 1, -0.3), Reduce(`*`, rep(1 - (-0.3), 32), 1))
    expect_identical(UPM(5, 0, 1.1), Reduce(`*`, rep(1.1, 5), 1))
    expect_identical(LPM(3, 1, -0.3, excess_ret = TRUE), Reduce(`*`, rep(1.3, 3), 1))
  })
  aud_case("integer degrees above 32 agree with the power to a few ulps", {
    tol <- 8 * .Machine$double.eps
    expect_equal(LPM(33, 1, -0.3), 1.3^33, tolerance = tol)
    expect_equal(UPM(40, 0, 1.1), 1.1^40, tolerance = tol)
    expect_equal(Co.LPM(40, -0.3, -0.3, 1, 1, degree_y = 1), 1.3^40 * 1.3, tolerance = tol)
  })
  aud_case("the cost of a huge integer degree does not grow with the degree", {
    # The left fold is a chain of dependent multiplies, about 1 ns each: 11 s for
    # degree 1e7 at n = 1000 and minutes for 2^31 - 1 at n = 100, so a fold build
    # fails the bound rather than hanging. pow takes microseconds, and the bound
    # is loose enough for slow or instrumented builds.
    el <- system.time(r1 <- LPM(2^31 - 1, 2, rep(1, 100)))[["elapsed"]]
    el <- el + system.time(r2 <- UPM(1e7, 0, rep(1, 1000)))[["elapsed"]]
    el <- el + system.time(r3 <- Co.LPM(2e8, c(0, 0.5, 2), c(0, 0.5, 2), 1, 1))[["elapsed"]]
    expect_lt(el, 60)
    expect_identical(c(r1, r2, r3), c(1, 1, 1 / 3))
  })
})

# Lean: none | BUGS: #5
test_that("U11: BUGS #5: the prefix path gives Inf, not NaN, when anchored sums overflow or gap powers underflow", {
  withr::local_seed(20261108)
  aud_case("the prefix path returns Inf, not NaN, when an anchored sum overflows", {
    x <- c(-1e308, 1e308)
    expect_identical(LPM(2, c(1e308, aud_pad), x)[1], LPM(2, 1e308, x))
    expect_identical(LPM(2, c(1e308, aud_pad), x)[1], Inf)
  })
  aud_case("the prefix recurrence returns Inf, not NaN, when a gap power underflows", {
    w <- c(-1e200, 0, 1e-170, 1, 2)
    expect_identical(LPM(4, c(2, 1.5, aud_pad), w)[1:2], c(LPM(4, 2, w), LPM(4, 1.5, w)))
    expect_identical(UPM(4, c(-1e-170, aud_pad), c(0, 1e200))[1], UPM(4, -1e-170, c(0, 1e200)))
  })
  aud_case("the prefix path has no NaN at a high degree on widely spread data", {
    withr::local_seed(2)
    z <- runif(40) * 1e10
    out <- LPM(32, z, z)
    expect_false(any(is.nan(out)))
    expect_identical(is.infinite(out), is.infinite(vapply(z, function(t) LPM(32, t, z), numeric(1))))
  })
})

# pin: NNS-Python reference values
test_that("U12: NNS-Python reference values: LPM and UPM at the mean, degrees 0 to 2", {
  withr::local_seed(20261108)
  x <- nns_py$x
  x_df <- as.data.frame(x)
  ref <- data.frame(fn = rep(c("LPM", "UPM"), each = 3), degree = rep(c(0, 1, 2), 2),
                    expected = c(0.49, 0.1032933, 0.02993767, 0.51, 0.1032933, 0.03027411))
  ref$label <- sprintf("%s(%g, mean(x), x)", ref$fn, ref$degree)
  aud_cases(ref, function(r) {
    expect_equal(match.fun(r$fn)(r$degree, mean(x), x), r$expected, tolerance = 1e-5)
  })
  # Input derived from a one-column frame: a named target and a named sample.
  aud_case("LPM(0, colMeans(x_df), unlist(x_df))", {
    expect_equal(LPM(0, colMeans(x_df), unlist(x_df)), 0.49, tolerance = 1e-5)
  })
  aud_case("UPM(0, colMeans(x_df), unlist(x_df))", {
    expect_equal(UPM(0, colMeans(x_df), unlist(x_df)), 0.51, tolerance = 1e-5)
  })
})
