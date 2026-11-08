# Property tests for the univariate partial moments LPM / UPM.
#
# Each test_that names the Lean theorem it mirrors in
# nns-proofs/NNSProofs/PartialMoments.lean (natural degrees) or
# nns-proofs/NNSProofs/RealDegree.lean (real degrees, prefix backend).
#
# Data regimes:
#   exact  : small dyadic data (n a power of 2) where double arithmetic is
#            exact, compared with expect_identical.
#   random : fixed seed, rnorm / rexp / rt(df = 3), rounded to force ties,
#            compared with an absolute tolerance n * eps * scale.
#   real   : datasets shipped with R, targets at mean, median, quantiles,
#            and exact data points (ties).
#
# Tests that fail on the baseline build are kept as the specification and
# carry a "# fails on 8898017: see BUGS.md #<id>" comment.

# ---- helpers ---------------------------------------------------------------

# Absolute deviation moment sum(|x - t|^d) / n: the natural scale of every
# partial-moment identity at (d, t). Naive summation of n terms, each with
# O(d) rounding, has absolute error below (n + d) * eps * scale; the factor 2
# covers the R-side oracle as well.
pm_tol <- function(d, t, x) {
  scale <- sum(abs(x - t)^d) / length(x)
  2 * (length(x) + d + 2) * .Machine$double.eps * scale
}

# The R wrappers vectorise over target; sapply forces the direct (<= 32
# targets) kernel one target at a time.
lpm_each <- function(d, ts, x) vapply(ts, function(t) LPM(d, t, x), numeric(1))
upm_each <- function(d, ts, x) vapply(ts, function(t) UPM(d, t, x), numeric(1))

pm_tol_v <- function(d, ts, x) vapply(ts, function(t) pm_tol(d, t, x), numeric(1))

# Largest excess of |lhs - rhs| over the tolerance; must be <= 0.
expect_within <- function(lhs, rhs, tol, label) {
  expect_lte(max(abs(lhs - rhs) - tol), 0, label = label)
}

targets_of <- function(x) {
  c(mean(x), median(x), quantile(x, c(0.05, 0.25, 0.75, 0.95), names = FALSE),
    x[1], min(x), max(x))
}

# exact regime: n = 8, dyadic values, a repeated value (tie) at 0.5
x_exact <- c(-2, -1, -0.5, 0, 0.5, 0.5, 1.5, 3)
t_exact <- c(-3, -1, 0, 0.25, 0.5, 1, 3, 4)

# random regime: ties forced by rounding
set.seed(20261108)
random_samples <- list(
  normal = round(rnorm(400), 1),
  expo   = round(rexp(300, 2), 2),
  t3     = round(rt(500, df = 3), 1)
)

# real regime
real_samples <- list(
  dax      = as.numeric(diff(log(EuStockMarkets[, "DAX"]))),
  waiting  = faithful$waiting,
  ozone    = as.numeric(na.omit(airquality$Ozone)),
  mag      = quakes$mag,
  mpg      = mtcars$mpg
)

all_samples <- c(random_samples, real_samples)

# ---- PartialMoments.lean ---------------------------------------------------

test_that("UPM_add_neg_one_pow_mul_LPM: master identity is exact on dyadic data, d = 0..6", {
  for (d in 0:6) {
    lhs <- upm_each(d, t_exact, x_exact) + (-1)^d * lpm_each(d, t_exact, x_exact)
    rhs <- vapply(t_exact, function(t) sum((x_exact - t)^d) / 8, numeric(1))
    expect_identical(lhs, rhs, label = paste("degree", d))
  }
})

test_that("UPM_add_neg_one_pow_mul_LPM: master identity on random and real data, d = 0..6", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    for (d in 0:6) {
      lhs <- upm_each(d, ts, x) + (-1)^d * lpm_each(d, ts, x)
      rhs <- vapply(ts, function(t) sum((x - t)^d) / length(x), numeric(1))
      expect_within(lhs, rhs, pm_tol_v(d, ts, x), label = sprintf("%s d=%d", nm, d))
    }
  }
})

test_that("UPM_add_neg_one_pow_mul_LPM: master identity holds for n = 1", {
  for (d in 0:6) {
    for (t in c(-1, 0.5, 1, 2.5)) {
      expect_identical(UPM(d, t, 1) + (-1)^d * LPM(d, t, 1), (1 - t)^d)
    }
  }
})

test_that("LPM_zero_add_UPM_zero: LPM(0) + UPM(0) is exactly 1, ties included", {
  expect_identical(lpm_each(0, t_exact, x_exact) + upm_each(0, t_exact, x_exact),
                   rep(1, length(t_exact)))
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    expect_identical(lpm_each(0, ts, x) + upm_each(0, ts, x), rep(1, length(ts)),
                     label = nm)
  }
})

test_that("LPM_zero_add_UPM_zero: degenerate samples (n = 1, constant) at the tie", {
  expect_identical(LPM(0, 5, 5) + UPM(0, 5, 5), 1)
  expect_identical(LPM(0, 5, rep(5, 7)), 1)
  expect_identical(UPM(0, 5, rep(5, 7)), 0)
})

test_that("UPM_one_sub_LPM_one: UPM(1) - LPM(1) = mean - t, exact on dyadic data", {
  expect_identical(upm_each(1, t_exact, x_exact) - lpm_each(1, t_exact, x_exact),
                   sum(x_exact) / 8 - t_exact)
})

test_that("UPM_one_sub_LPM_one: UPM(1) - LPM(1) = mean - t on random and real data", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    expect_within(upm_each(1, ts, x) - lpm_each(1, ts, x), sum(x) / length(x) - ts,
                  pm_tol_v(1, ts, x), label = nm)
  }
})

test_that("LPM_add_UPM_of_even: LPM(d) + UPM(d) is the d-th central moment for even d", {
  for (d in c(2, 4, 6)) {
    expect_identical(lpm_each(d, t_exact, x_exact) + upm_each(d, t_exact, x_exact),
                     vapply(t_exact, function(t) sum((x_exact - t)^d) / 8, numeric(1)))
    for (nm in names(all_samples)) {
      x <- all_samples[[nm]]
      ts <- targets_of(x)
      expect_within(lpm_each(d, ts, x) + upm_each(d, ts, x),
                    vapply(ts, function(t) sum((x - t)^d) / length(x), numeric(1)),
                    pm_tol_v(d, ts, x), label = sprintf("%s d=%d", nm, d))
    }
  }
})

test_that("LPM_two_add_UPM_two_mean: LPM(2, mean) + UPM(2, mean) is the population variance", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    n <- length(x)
    mu <- mean(x)
    expect_lte(abs(LPM(2, mu, x) + UPM(2, mu, x) - var(x) * (n - 1) / n),
               pm_tol(2, mu, x), label = nm)
  }
})

test_that("LPM_nonneg / UPM_nonneg: partial moments are nonnegative at every degree", {
  for (x in all_samples) {
    ts <- targets_of(x)
    for (d in c(0, 0.3, 1, 1.5, 2, 3, 6)) {
      expect_true(all(lpm_each(d, ts, x) >= 0))
      expect_true(all(upm_each(d, ts, x) >= 0))
    }
  }
})

test_that("LPM_zero_eq_ecdf: LPM(0, t, x) is the right-continuous ECDF, scalar targets", {
  expect_identical(lpm_each(0, t_exact, x_exact),
                   vapply(t_exact, function(t) sum(x_exact <= t) / 8, numeric(1)))
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- c(targets_of(x), head(unique(x), 10))
    expect_identical(lpm_each(0, ts, x), ecdf(x)(ts), label = nm)
  }
})

test_that("LPM_zero_eq_ecdf: LPM(0) with more than 32 targets is the ECDF at every data value", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- rep_len(sort(unique(x)), 60)
    expect_identical(LPM(0, ts, x), ecdf(x)(ts), label = nm)
  }
})

test_that("LPM_zero_monotone / LPM_one_monotone: degree 0 and 1 LPM nondecreasing in t", {
  for (x in c(list(x_exact), all_samples)) {
    ts <- sort(c(targets_of(x), head(unique(x), 20)))
    expect_true(all(diff(lpm_each(0, ts, x)) >= 0))
    expect_true(all(diff(lpm_each(1, ts, x)) >= 0))
  }
})

# ---- RealDegree.lean -------------------------------------------------------

test_that("rLPM_natCast / rUPM_natCast: integer-typed degree equals double degree", {
  for (k in 0:4) {
    expect_identical(lpm_each(k, t_exact, x_exact), lpm_each(as.integer(k), t_exact, x_exact))
    expect_identical(upm_each(k, t_exact, x_exact), upm_each(as.integer(k), t_exact, x_exact))
  }
})

test_that("rLPM_add_rUPM: LPM(d) + UPM(d) = mean |x - t|^d for every real d, dyadic data", {
  # |x - 1| takes values in {0, 0.25, 1, 2.25, 4}: squares of dyadics, so the
  # d = 0.5 and d = 1.5 powers are exact too; two observations tie at t = 1.
  x <- c(-3, -1.25, 0, 0.75, 1, 1, 2, 3.25)
  for (d in c(0, 0.5, 1, 1.5, 2, 3)) {
    expect_identical(LPM(d, 1, x) + UPM(d, 1, x), sum(abs(x - 1)^d) / 8,
                     label = paste("degree", d))
  }
})

test_that("rLPM_add_rUPM: LPM(d) + UPM(d) = mean |x - t|^d for real d on random and real data", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    for (d in c(0, 0.25, 0.5, 1, 1.5, 2.5, 3, 3.7)) {
      expect_within(lpm_each(d, ts, x) + upm_each(d, ts, x),
                    vapply(ts, function(t) sum(abs(x - t)^d) / length(x), numeric(1)),
                    pm_tol_v(d, ts, x), label = sprintf("%s d=%g", nm, d))
    }
  }
})

test_that("rLPM_eq_max / rUPM_eq_max: for d > 0, LPM and UPM are means of positive parts", {
  for (d in c(1, 2, 3, 5)) {
    expect_identical(lpm_each(d, t_exact, x_exact),
                     vapply(t_exact, function(t) sum(pmax(t - x_exact, 0)^d) / 8, numeric(1)))
    expect_identical(upm_each(d, t_exact, x_exact),
                     vapply(t_exact, function(t) sum(pmax(x_exact - t, 0)^d) / 8, numeric(1)))
  }
})

test_that("rLPM_eq_max / rUPM_eq_max: positive-part form for real d > 0 on random and real data", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    for (d in c(1e-9, 0.5, 1.5, 2.5)) {
      expect_within(lpm_each(d, ts, x),
                    vapply(ts, function(t) sum(pmax(t - x, 0)^d) / length(x), numeric(1)),
                    pm_tol_v(d, ts, x), label = sprintf("LPM %s d=%g", nm, d))
      expect_within(upm_each(d, ts, x),
                    vapply(ts, function(t) sum(pmax(x - t, 0)^d) / length(x), numeric(1)),
                    pm_tol_v(d, ts, x), label = sprintf("UPM %s d=%g", nm, d))
    }
  }
})

test_that("rLPM_eq_strict: for d > 0 a tie at the target contributes nothing to LPM", {
  x <- c(1, 2, 2, 3)
  expect_identical(LPM(0, 2, x), 0.75)
  for (d in c(1e-9, 0.5, 1, 2)) {
    expect_identical(LPM(d, 2, x), 0.25, label = paste("degree", d))
  }
})

test_that("rUPM_eq_rLPM_neg: for d > 0, UPM(d, t, x) = LPM(d, -t, -x)", {
  for (d in c(0.5, 1, 1.5, 2, 3)) {
    expect_identical(upm_each(d, t_exact, x_exact), lpm_each(d, -t_exact, -x_exact))
    for (x in all_samples) {
      ts <- targets_of(x)
      expect_identical(upm_each(d, ts, x), lpm_each(d, -ts, -x))
    }
  }
})

test_that("rUPM_eq_rLPM_neg: reflection fails at d = 0 when the target is tied", {
  x <- c(1, 2, 2, 3)
  expect_identical(UPM(0, 2, x), 0.25)
  expect_identical(LPM(0, -2, -x), 0.75)
})

test_that("rLPM_add_const / rUPM_add_const: translation of sample and target is exact on dyadic data", {
  # data on the 1/8 grid (random, with ties) shifted by 1024 stays exact
  set.seed(7)
  x <- c(x_exact, round(rnorm(200) * 8) / 8)
  ts <- c(t_exact, x[9:20], round(mean(x) * 8) / 8)
  for (d in c(0, 0.5, 1, 1.5, 2, 3, 4)) {
    expect_identical(lpm_each(d, ts + 1024, x + 1024), lpm_each(d, ts, x))
    expect_identical(upm_each(d, ts + 1024, x + 1024), upm_each(d, ts, x))
  }
})

test_that("rLPM_add_const / rUPM_add_const: translation is exact on integer-valued real data", {
  for (x in list(faithful$waiting, as.numeric(na.omit(airquality$Ozone)), mtcars$hp)) {
    ts <- c(median(x), x[1:5], quantile(x, c(0.25, 0.75), type = 1, names = FALSE))
    for (d in c(0, 0.5, 1, 2, 3)) {
      expect_identical(lpm_each(d, ts - 37, x - 37), lpm_each(d, ts, x))
      expect_identical(upm_each(d, ts - 37, x - 37), upm_each(d, ts, x))
    }
  }
})

test_that("rLPM_mul_const / rUPM_mul_const: scaling by 4 multiplies integer-degree moments by 4^d exactly", {
  for (x in c(list(x_exact), all_samples)) {
    ts <- targets_of(x)
    for (d in 0:5) {
      expect_identical(lpm_each(d, 4 * ts, 4 * x), 4^d * lpm_each(d, ts, x))
      expect_identical(upm_each(d, 4 * ts, 4 * x), 4^d * upm_each(d, ts, x))
    }
  }
})

test_that("rLPM_mul_const / rUPM_mul_const: scaling by a > 0 multiplies by a^d for real d", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    for (a in c(0.37, 4, 1000)) {
      for (d in c(0.5, 1.5, 2.5, 3)) {
        tol <- (d + 2) * a^d * pm_tol_v(d, ts, x)
        expect_within(lpm_each(d, a * ts, a * x), a^d * lpm_each(d, ts, x), tol,
                      label = sprintf("LPM %s a=%g d=%g", nm, a, d))
        expect_within(upm_each(d, a * ts, a * x), a^d * upm_each(d, ts, x), tol,
                      label = sprintf("UPM %s a=%g d=%g", nm, a, d))
      }
    }
  }
})

test_that("rLPM_const / rUPM_const: constant sample gives [c <= t](t - c)^d and [c < t](c - t)^d", {
  for (n in c(1, 5)) {
    x <- rep(2, n)
    expect_identical(LPM(0, 2, x), 1)
    expect_identical(UPM(0, 2, x), 0)
    for (d in c(0.5, 1, 2, 3)) {
      expect_identical(LPM(d, 2, x), 0)
      expect_identical(UPM(d, 2, x), 0)
      expect_identical(LPM(d, 3.5, x), 1.5^d)
      expect_identical(UPM(d, 3.5, x), 0)
      expect_identical(UPM(d, 0.5, x), 1.5^d)
      expect_identical(LPM(d, 0.5, x), 0)
    }
  }
})

test_that("rLPM_monotone: LPM is nondecreasing in t for real d >= 0 (direct path, exact)", {
  for (x in c(list(x_exact), all_samples)) {
    ts <- sort(c(targets_of(x), head(unique(x), 20)))
    for (d in c(0, 0.3, 0.5, 1, 1.5, 2, 3)) {
      expect_true(all(diff(lpm_each(d, ts, x)) >= 0), label = paste("degree", d))
    }
  }
})

test_that("rUPM_antitone: UPM is nonincreasing in t for real d >= 0 (direct path, exact)", {
  for (x in c(list(x_exact), all_samples)) {
    ts <- sort(c(targets_of(x), head(unique(x), 20)))
    for (d in c(0, 0.3, 0.5, 1, 1.5, 2, 3)) {
      expect_true(all(diff(upm_each(d, ts, x)) <= 0), label = paste("degree", d))
    }
  }
})

test_that("rLPM_monotone / rUPM_antitone: monotone in t on the > 32 target path", {
  for (x in all_samples) {
    ts <- seq(min(x) - 1, max(x) + 1, length.out = 80)
    for (d in c(0, 0.5, 1, 2, 3)) {
      tol <- pm_tol(d, max(ts), x) + pm_tol(d, min(ts), x)
      expect_true(all(diff(LPM(d, ts, x)) >= -tol), label = paste("LPM degree", d))
      expect_true(all(diff(UPM(d, ts, x)) <= tol), label = paste("UPM degree", d))
    }
  }
})

test_that("rLPM_continuous / rUPM_continuous: for d > 0 a perturbation h moves the moment by O(h)", {
  # |a^d - b^d| <= d |a - b| max(a, b)^(d - 1) for d >= 1, and <= |a - b|^d
  # for 0 < d < 1, applied per observation with a, b the positive parts.
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    h <- 1e-6 * diff(range(x))
    for (t in c(x[1], median(x))) {
      M <- max(abs(x - t)) + h
      for (d in c(0.25, 0.5, 1, 1.5, 2, 3)) {
        bound <- if (d >= 1) d * h * M^(d - 1) else h^d
        bound <- bound + pm_tol(d, t, x) + pm_tol(d, t + h, x)
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

test_that("rLPM_continuous needs d > 0: degree-0 LPM jumps by the tie share at a data point", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    t <- x[1]
    h <- 1e-6 * diff(range(x))
    tie_share <- sum(x == t) / length(x)
    # k/n - j/n is (k - j)/n up to one rounding
    expect_gte(LPM(0, t, x) - LPM(0, t - h, x), tie_share - 2 * .Machine$double.eps,
               label = nm)
  }
})

test_that("rLPM_convexOn / rUPM_convexOn: midpoint inequality in t for d >= 1", {
  for (nm in names(c(list(exact = x_exact), all_samples))) {
    x <- c(list(exact = x_exact), all_samples)[[nm]]
    ts <- targets_of(x)
    pairs <- utils::combn(length(ts), 2)
    a <- ts[pairs[1, ]]
    b <- ts[pairs[2, ]]
    for (d in c(1, 1.5, 2, 3, 4)) {
      tol <- vapply(seq_along(a), function(i) pm_tol(d, a[i], x) + pm_tol(d, b[i], x),
                    numeric(1))
      expect_true(all(lpm_each(d, (a + b) / 2, x) <=
                        (lpm_each(d, a, x) + lpm_each(d, b, x)) / 2 + tol),
                  label = sprintf("LPM %s d=%g", nm, d))
      expect_true(all(upm_each(d, (a + b) / 2, x) <=
                        (upm_each(d, a, x) + upm_each(d, b, x)) / 2 + tol),
                  label = sprintf("UPM %s d=%g", nm, d))
    }
  }
})

test_that("not_rLPM_convexOn_of_lt_one: for 0 < d < 1 the one-point LPM violates the midpoint inequality", {
  for (d in c(0.25, 0.5, 0.9)) {
    expect_gt(LPM(d, 0.5, 0), (LPM(d, 0, 0) + LPM(d, 1, 0)) / 2)
    expect_identical(LPM(d, 0.5, 0), 0.5^d)
  }
})

test_that("rLPM_tendsto_zero': as d -> 0+, LPM tends to the strict ECDF #{x < t} / n", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    expect_within(lpm_each(1e-9, ts, x),
                  vapply(ts, function(t) sum(x < t) / length(x), numeric(1)), 1e-6,
                  label = nm)
  }
})

test_that("rLPM_tendsto_zero / rLPM_limit_zero_lt_LPM_zero: the d -> 0+ limit drops the tie share", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    t <- x[1]
    tie_share <- sum(x == t) / length(x)
    expect_gt(tie_share, 0)
    expect_lte(abs(LPM(0, t, x) - LPM(1e-9, t, x) - tie_share), 1e-6, label = nm)
  }
})

test_that("rUPM_tendsto_zero: UPM is right-continuous in the degree at 0", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- targets_of(x)
    expect_within(upm_each(1e-9, ts, x), upm_each(0, ts, x), 1e-6, label = nm)
  }
})

# ---- RealDegree.lean: prefix backend (> 32 targets) vs direct path ---------
#
# The specification is rLPMVecExact_eq_rLPM / rUPMVecExact_eq_rUPM: the value
# for (d, t, x) must not depend on how many targets share the call. At
# 8898017 the prefix backend snaps degrees within 1e-12 of an integer
# (modelled by rLPMVec; rLPMVec_sub_rLPM_tendsto proves the gap tends to the
# tie share as d -> 0+, rLPMVec_sub_rLPM_tendsto_natCast that it is only
# o(1) near k >= 1). Those theorems describe the bug, so no test asserts them.

test_that("rLPMVecExact_eq_rLPM: degree 0 is bit-identical on both paths, ties included", {
  for (x in c(list(x_exact), all_samples)) {
    ts <- rep_len(c(targets_of(x), head(x, 10)), 40)
    expect_identical(LPM(0, ts, x), lpm_each(0, ts, x))
    expect_identical(UPM(0, ts, x), upm_each(0, ts, x))
  }
})

test_that("rLPMVecExact_eq_rLPM: integer degrees 1..6 agree to n * eps on both paths", {
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- rep_len(c(targets_of(x), head(x, 10)), 40)
    for (d in 1:6) {
      tol <- pm_tol_v(d, ts, x)
      expect_within(LPM(d, ts, x), lpm_each(d, ts, x), tol,
                    label = sprintf("LPM %s d=%d", nm, d))
      expect_within(UPM(d, ts, x), upm_each(d, ts, x), tol,
                    label = sprintf("UPM %s d=%d", nm, d))
    }
  }
})

test_that("rLPMVecExact_eq_rLPM: 32 vs 33 targets and degree > 32 give the direct value", {
  x <- c(1, 2, 2, 3)
  expect_identical(LPM(0, rep(2, 32), x), rep(0.75, 32))
  expect_identical(LPM(0, rep(2, 33), x), rep(0.75, 33))
  expect_identical(LPM(1, rep(2, 33), x), rep(0.25, 33))
  z <- random_samples$normal
  ts <- rep_len(targets_of(z), 40)
  expect_identical(LPM(33, ts, z), lpm_each(33, ts, z))
})

test_that("rLPMVecExact_eq_rLPM: non-integer degrees are bit-identical on both paths", {
  for (x in c(list(x_exact), all_samples)) {
    ts <- rep_len(c(targets_of(x), head(x, 10)), 40)
    for (d in c(0.5, 1.5, 2.5, 2 + 2e-12)) {
      expect_identical(LPM(d, ts, x), lpm_each(d, ts, x), label = paste("LPM degree", d))
      expect_identical(UPM(d, ts, x), upm_each(d, ts, x), label = paste("UPM degree", d))
    }
  }
})

test_that("rLPMVecExact_eq_rLPM: degree 1e-13 at a tied target gives the same LPM on both paths", {
  # fails on 8898017: see BUGS.md #0 (prefix path snaps 1e-13 to 0 and counts the tie)
  x <- c(1, 2, 2, 3)
  expect_identical(LPM(1e-13, 2, x), 0.25)
  expect_identical(LPM(1e-13, rep(2, 33), x), rep(0.25, 33))
})

test_that("rLPMVecExact_eq_rLPM: computed near-zero degree 0.1 + 0.2 - 0.3 gives the same LPM on both paths", {
  # fails on 8898017: see BUGS.md #0
  d <- 0.1 + 0.2 - 0.3
  for (x in c(list(c(1, 2, 2, 3)), real_samples)) {
    ts <- rep_len(c(x[1], median(x)), 40)
    expect_identical(LPM(d, ts, x), lpm_each(d, ts, x))
  }
})

test_that("rLPMVecExact_eq_rLPM / rUPMVecExact_eq_rUPM: near-integer degrees agree to n * eps on both paths", {
  # fails on 8898017: see BUGS.md #0 (degrees within 1e-12 of an integer are
  # evaluated at the integer on the prefix path; gap ~1e-12, about 1e4 * tol)
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- rep_len(c(targets_of(x), head(x, 10)), 40)
    for (d in c(1 + 1e-13, 2 - 5e-13, 2 + 5e-13, 3 + 9e-13)) {
      tol <- pm_tol_v(d, ts, x)
      expect_within(LPM(d, ts, x), lpm_each(d, ts, x), tol,
                    label = sprintf("LPM %s d=%.15g", nm, d))
      expect_within(UPM(d, ts, x), upm_each(d, ts, x), tol,
                    label = sprintf("UPM %s d=%.15g", nm, d))
    }
  }
})

test_that("rUPMVecExact_eq_rUPM: degree 1e-13 gives the same UPM on both paths", {
  # fails on 8898017: see BUGS.md #0 (prefix path evaluates UPM(0); direct
  # path gives mean((x - t)_+^1e-13), which differs by ~1e-13 * mean log gap)
  for (nm in names(all_samples)) {
    x <- all_samples[[nm]]
    ts <- rep_len(targets_of(x), 40)
    expect_within(UPM(1e-13, ts, x), upm_each(1e-13, ts, x), pm_tol_v(1e-13, ts, x),
                  label = nm)
  }
})
