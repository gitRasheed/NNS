# Property tests for the co- and divergent partial moments Co.LPM, Co.UPM,
# D.LPM, D.UPM, the PM.matrix quadrant matrices and the n-dimensional
# Co.LPM_nD / Co.UPM_nD / DPM_nD.
#
# Each test_that names the Lean theorem it mirrors in
# nns-proofs/NNSProofs/CoPartialMoments.lean (bivariate co-moments) or
# nns-proofs/NNSProofs/CoMomentsMatrix.lean (degree_y, PM.matrix). The nD
# functions have no Lean statement; their tests cite the SPEC.md section
# instead, so the trace stays honest.
#
# Data regimes:
#   exact  : n = 8 dyadic data (so PM.matrix's multiply-by-1/n equals the
#            co-moments' divide-by-n bit for bit), compared with
#            expect_identical.
#   random : fixed seed, rnorm / rexp / rt(df = 3), rounded to force ties,
#            compared with an absolute tolerance n * eps * scale.
#   real   : datasets shipped with R, targets at the mean, median, quartiles
#            and at an exact observation (a tie in both coordinates).
#
# Targets are always passed as numbers when results of different functions
# are compared: PM.matrix (Rcpp colMeans), Co.LPM (Rcpp sugar mean) and R's
# mean() can differ by an ulp, which moves tie membership at degree 0. The one
# exception is the covariance_decomposition test with target_x = target_y =
# NULL: it exercises the bivariate mean-default path, compares each function
# only with others resolving the same NULL target, and checks the degree-0
# quadrants through their sum, which is 1 whatever the tie membership.
#
# Every case and degree of a property is collected into one vector and
# checked with a single expectation; a failure names the offending
# case/degree entries.
#
# Tests that fail on the baseline build are kept as the specification and
# carry a "# fails on 8898017: see BUGS.md #<id>" comment.

# ---- helpers ---------------------------------------------------------------

eps <- .Machine$double.eps

# Natural scale of a co-moment identity at degrees (dx, dy): the mean of
# |x - tx|^dx |y - ty|^dy. Naive in-order summation of n terms with O(d)
# rounding per term stays below (n + d) * eps * scale; k covers the number
# of moments combined plus the R-side oracle.
co_tol <- function(c, dx, dy, k = 4) {
  k * (length(c$x) + dx + dy + 2) * eps * mean(abs(c$x - c$tx)^dx * abs(c$y - c$ty)^dy)
}

# Matrix scale: entries are compared against the largest entry.
mat_tol <- function(M, n, k = 4) k * (n + ncol(M)) * eps * max(abs(M))

# Elementwise |actual - expected| <= tol, one expectation, failing entries named.
# An empty comparison fails (a case list that collapsed to nothing would
# otherwise pass vacuously), and tol must be a scalar or match actual.
expect_within <- function(actual, expected, tol) {
  non_empty <- length(actual) > 0
  same_shape <- length(expected) %in% c(1L, length(actual))
  tol_shape <- length(tol) %in% c(1L, length(actual))
  if (!(non_empty && same_shape && tol_shape)) {
    return(expect(FALSE, paste(c(if (!non_empty) "empty comparison",
                                 if (!same_shape) "length mismatch (expected)",
                                 if (!tol_shape) "length mismatch (tol)"), collapse = "; ")))
  }
  err <- abs(actual - expected)
  tol <- rep_len(tol, length(err))
  bad <- which(!(is.finite(err) & err <= tol))
  where <- names(actual)
  if (is.null(where)) where <- as.character(seq_along(err))
  msg <- sprintf("%s: |diff| = %.3g > tol %.3g", where[bad], err[bad], tol[bad])
  expect(length(bad) == 0, paste(c("not within tolerance:", head(msg, 8)), collapse = "\n"))
}

# Collect f over named cases / over degree pairs into one named vector.
over <- function(cases, f) unlist(lapply(cases, f))
at_degrees <- function(pairs, f) {
  unlist(lapply(pairs, function(p) setNames(f(p[1], p[2]), paste0("d", p[1], p[2]))))
}

min_eigen <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

degree_pairs <- list(c(0, 0), c(1, 1), c(2, 1), c(1, 0), c(0, 2), c(2, 2))
pairs_03 <- lapply(0:15, function(k) c(k %/% 4, k %% 4))
pairs_02 <- lapply(0:8, function(k) c(k %/% 3, k %% 3))

# One case per target pair: mean, median, quartiles and an observation
# (x[i], y[i]) so that both coordinates tie their target.
co_cases <- function(x, y, label) {
  i <- ceiling(length(x) / 3)
  q <- function(v, p) quantile(v, p, names = FALSE)
  targets <- list(
    mean   = c(mean(x), mean(y)),
    median = c(median(x), median(y)),
    q25    = c(q(x, 0.25), q(y, 0.25)),
    q75    = c(q(x, 0.75), q(y, 0.75)),
    point  = c(x[i], y[i])
  )
  lapply(setNames(targets, paste(label, names(targets))), function(t)
    list(x = x, y = y, tx = t[1], ty = t[2]))
}

# exact regime: n = 8 dyadic values, ties at the targets (0.5, 0.5)
x_exact <- c(-2, -1, -0.5, 0, 0.5, 0.5, 1.5, 3)
y_exact <- c(1, -1, 0.5, 0.5, -2, 0.5, 2, -0.25)
exact_cases <- list(
  "exact ties" = list(x = x_exact, y = y_exact, tx = 0.5, ty = 0.5),
  "exact off"  = list(x = x_exact, y = y_exact, tx = 0.25, ty = -1.5),
  "exact out"  = list(x = x_exact, y = y_exact, tx = -4, ty = 4)
)

# random regime: ties forced by rounding
set.seed(20261108)
z_norm <- round(rnorm(400), 1)
z_t3 <- round(rt(300, df = 3), 1)
random_cases <- c(
  co_cases(z_norm, round(0.6 * z_norm + rnorm(400), 1), "normal"),
  co_cases(round(rexp(250, 2), 2), round(rexp(250, 1), 2), "expo"),
  co_cases(z_t3, round(-0.4 * z_t3 + rt(300, df = 3), 1), "t3")
)

# real regime
eu <- diff(log(EuStockMarkets))
aq <- na.omit(airquality[, c("Ozone", "Temp")])
real_cases <- c(
  co_cases(as.numeric(eu[, "DAX"]), as.numeric(eu[, "FTSE"]), "eu"),
  co_cases(faithful$eruptions, faithful$waiting, "faithful"),
  co_cases(as.numeric(aq$Ozone), as.numeric(aq$Temp), "airquality"),
  co_cases(quakes$mag, quakes$depth, "quakes"),
  co_cases(mtcars$mpg, mtcars$hp, "mtcars")
)

sample_cases <- c(random_cases, real_cases)
all_cases <- c(exact_cases, sample_cases)

# Matrices for PM.matrix: exact (n = 8, dyadic), random with ties, real.
X_exact <- cbind(a = x_exact, b = y_exact, c = c(0, 0.5, -1, 2, 2, -0.5, 1, 0.5))
t_exact <- c(0.5, 0.5, 0.5)
X_random <- cbind(u = round(rnorm(60), 1), v = round(rexp(60), 1),
                  w = round(rt(60, df = 3), 1))
X_random <- cbind(X_random, dup = X_random[, "u"])  # rank deficient on purpose
pm_samples <- list(
  random = X_random,
  eu     = eu,
  mtcars = as.matrix(mtcars),
  quakes = as.matrix(quakes[, c("lat", "long", "depth", "mag")])
)
pm_all <- c(list(exact = X_exact), pm_samples)

# The four quadrants with the Lean argument conventions.
co_lpm <- function(c, dx, dy) Co.LPM(dx, c$x, c$y, c$tx, c$ty, dy)
co_upm <- function(c, dx, dy) Co.UPM(dx, c$x, c$y, c$tx, c$ty, dy)
d_lpm  <- function(c, dl, du) D.LPM(dl, du, c$x, c$y, c$tx, c$ty)
d_upm  <- function(c, dl, du) D.UPM(dl, du, c$x, c$y, c$tx, c$ty)

co_moment <- function(c, d = 1, e = 1) mean((c$x - c$tx)^d * (c$y - c$ty)^e)
decomposition <- function(c) {
  co_upm(c, 1, 1) + co_lpm(c, 1, 1) - d_lpm(c, 1, 1) - d_upm(c, 1, 1)
}
quadrant_sum0 <- function(c) {
  co_lpm(c, 0, 0) + co_upm(c, 0, 0) + d_lpm(c, 0, 0) + d_upm(c, 0, 0)
}
mixed_master_lhs <- function(c, d, e) {
  co_upm(c, d, e) + (-1)^d * d_upm(c, d, e) + (-1)^e * d_lpm(c, e, d) +
    (-1)^(d + e) * co_lpm(c, d, e)
}

# ---- CoPartialMoments.lean -------------------------------------------------

test_that("coMoment_decomposition: exact dyadic data, CoUPM + CoLPM - DLPM - DUPM is identical to the co-moment", {
  expect_identical(over(exact_cases, decomposition), over(exact_cases, co_moment))
})

test_that("coMoment_decomposition: random and real data, every target, within n * eps", {
  expect_within(over(sample_cases, decomposition), over(sample_cases, co_moment),
                over(sample_cases, function(c) co_tol(c, 1, 1)))
})

test_that("coMoment_decomposition: n = 1 is exact", {
  expect_identical(decomposition(list(x = 3, y = 5, tx = 1, ty = 2)), 6)
  expect_identical(decomposition(list(x = 3, y = -5, tx = 1, ty = 2)), -14)
})

test_that("covariance_decomposition: targets at the means give the population covariance", {
  mean_cases <- sample_cases[grep("mean$", names(sample_cases))]
  pop_cov <- function(c) cov(c$x, c$y) * (length(c$x) - 1) / length(c$x)
  expect_within(over(mean_cases, decomposition), over(mean_cases, pop_cov),
                over(mean_cases, function(c) co_tol(c, 1, 1)))
})

test_that("covariance_decomposition: target_x = target_y = NULL resolves to the means (bivariate mean-default path)", {
  # one case per data set; its numeric targets are mean(x), mean(y), so co_tol
  # is evaluated at the targets NULL resolves to
  mean_cases <- sample_cases[grep("mean$", names(sample_cases))]
  null_t <- function(c) list(x = c$x, y = c$y, tx = NULL, ty = NULL)
  pop_cov <- function(c) cov(c$x, c$y) * (length(c$x) - 1) / length(c$x)
  expect_within(over(mean_cases, function(c) decomposition(null_t(c))), over(mean_cases, pop_cov),
                over(mean_cases, function(c) co_tol(c, 1, 1)))
  expect_within(over(mean_cases, function(c) quadrant_sum0(null_t(c))), 1, 4 * eps)
})

test_that("quadrants_sum_eq_one: exact dyadic data, the degree-0 quadrants sum to exactly 1", {
  expect_identical(unname(over(exact_cases, quadrant_sum0)), rep(1, length(exact_cases)))
})

test_that("quadrants_sum_eq_one: random and real data, ties at the target, within 4 eps", {
  expect_within(over(sample_cases, quadrant_sum0), 1, 4 * eps)
})

test_that("quadrants_sum_eq_one: n = 1, and a constant sample at its target sits in the co-lower quadrant", {
  expect_identical(quadrant_sum0(list(x = 3, y = 5, tx = 3, ty = 5)), 1)
  expect_identical(quadrant_sum0(list(x = 3, y = 5, tx = 4, ty = 4)), 1)
  const <- list(x = rep(2, 7), y = rep(-1, 7), tx = 2, ty = -1)
  expect_identical(c(co_lpm(const, 0, 0), co_upm(const, 0, 0),
                     d_lpm(const, 0, 0), d_upm(const, 0, 0)), c(1, 0, 0, 0))
})

test_that("quadrants_sum_eq_one: degree 0 vs degree > 0 at full ties (0^0 = 1 but 0^d = 0)", {
  const <- list(x = rep(2, 7), y = rep(-1, 7), tx = 2, ty = -1)
  expect_identical(c(co_lpm(const, 1, 1), co_upm(const, 1, 1), d_lpm(const, 1, 1),
                     d_upm(const, 1, 1), co_lpm(const, 0, 1), co_lpm(const, 1, 0)),
                   numeric(6))
})

test_that("CoLPM_nonneg, CoUPM_nonneg, DLPM_nonneg, DUPM_nonneg: every quadrant is >= 0", {
  quads <- function(c) at_degrees(degree_pairs, function(d, e) {
    c(co_lpm(c, d, e), co_upm(c, d, e), d_lpm(c, d, e), d_upm(c, d, e))
  })
  expect_true(all(over(all_cases, quads) >= 0))
})

test_that("DLPM_swap: D.LPM(dl, du, x, y, tx, ty) is identical to D.UPM(dl, du, y, x, ty, tx)", {
  lhs <- function(c) at_degrees(degree_pairs, function(l, u) D.LPM(l, u, c$x, c$y, c$tx, c$ty))
  rhs <- function(c) at_degrees(degree_pairs, function(l, u) D.UPM(l, u, c$y, c$x, c$ty, c$tx))
  expect_identical(over(all_cases, lhs), over(all_cases, rhs))
})

test_that("CoLPM_swap: Co.LPM(dx, x, y, tx, ty, dy) is identical to Co.LPM(dy, y, x, ty, tx, dx)", {
  lhs <- function(c) at_degrees(degree_pairs, function(d, e) Co.LPM(d, c$x, c$y, c$tx, c$ty, e))
  rhs <- function(c) at_degrees(degree_pairs, function(d, e) Co.LPM(e, c$y, c$x, c$ty, c$tx, d))
  expect_identical(over(all_cases, lhs), over(all_cases, rhs))
})

test_that("CoUPM_swap: Co.UPM(dx, x, y, tx, ty, dy) is identical to Co.UPM(dy, y, x, ty, tx, dx)", {
  lhs <- function(c) at_degrees(degree_pairs, function(d, e) Co.UPM(d, c$x, c$y, c$tx, c$ty, e))
  rhs <- function(c) at_degrees(degree_pairs, function(d, e) Co.UPM(e, c$y, c$x, c$ty, c$tx, d))
  expect_identical(over(all_cases, lhs), over(all_cases, rhs))
})

self_cases <- function(x, ts) lapply(ts, function(t) list(x = x, y = x, tx = t, ty = t))
exact_self <- self_cases(x_exact, c(-1, 0.5, 0.25, 4))

test_that("CoLPM_self: exact dyadic data, Co.LPM(d, x, x, t, t, e) is identical to LPM(d + e, t, x)", {
  lhs <- function(c) at_degrees(pairs_02, function(d, e) co_lpm(c, d, e))
  rhs <- function(c) at_degrees(pairs_02, function(d, e) LPM(d + e, c$tx, c$x))
  expect_identical(over(exact_self, lhs), over(exact_self, rhs))
})

test_that("CoUPM_self: exact dyadic data, Co.UPM(d, x, x, t, t, e) is identical to UPM(d + e, t, x)", {
  lhs <- function(c) at_degrees(pairs_02, function(d, e) co_upm(c, d, e))
  rhs <- function(c) at_degrees(pairs_02, function(d, e) UPM(d + e, c$tx, c$x))
  expect_identical(over(exact_self, lhs), over(exact_self, rhs))
})

sample_self <- lapply(sample_cases, function(c) list(x = c$x, y = c$x, tx = c$tx, ty = c$tx))

test_that("CoLPM_self: random and real data at tied targets, within n * eps", {
  expect_within(
    over(sample_self, function(c) at_degrees(pairs_02, function(d, e) co_lpm(c, d, e))),
    over(sample_self, function(c) at_degrees(pairs_02, function(d, e) LPM(d + e, c$tx, c$x))),
    over(sample_self, function(c) at_degrees(pairs_02, function(d, e) co_tol(c, d, e))))
})

test_that("CoUPM_self: random and real data at tied targets, within n * eps", {
  expect_within(
    over(sample_self, function(c) at_degrees(pairs_02, function(d, e) co_upm(c, d, e))),
    over(sample_self, function(c) at_degrees(pairs_02, function(d, e) UPM(d + e, c$tx, c$x))),
    over(sample_self, function(c) at_degrees(pairs_02, function(d, e) co_tol(c, d, e))))
})

eu_self_40 <- local({
  x <- as.numeric(eu[, "DAX"])
  ts <- c(quantile(x, seq(0.01, 0.99, length.out = 30), names = FALSE), x[1:10])
  list(x = x, y = x, tx = ts, ty = ts)
})
self_tol_40 <- function(c, d, e) {
  2 * (length(c$x) + d + e + 2) * eps *
    vapply(c$tx, function(t) mean(abs(c$x - t)^(d + e)), numeric(1))
}

test_that("CoLPM_self: 40 targets (LPM prefix backend) agree with the co-moment kernel", {
  c <- eu_self_40
  expect_within(at_degrees(pairs_02, function(d, e) co_lpm(c, d, e)),
                at_degrees(pairs_02, function(d, e) LPM(d + e, c$tx, c$x)),
                at_degrees(pairs_02, function(d, e) self_tol_40(c, d, e)))
})

test_that("CoUPM_self: 40 targets (UPM prefix backend) agree with the co-moment kernel", {
  c <- eu_self_40
  expect_within(at_degrees(pairs_02, function(d, e) co_upm(c, d, e)),
                at_degrees(pairs_02, function(d, e) UPM(d + e, c$tx, c$x)),
                at_degrees(pairs_02, function(d, e) self_tol_40(c, d, e)))
})

test_that("CoLPM_self / CoUPM_self: n = 1 and constant samples", {
  edge <- list(n1_low = list(x = 3, y = 3, tx = 5, ty = 5),
               n1_high = list(x = 3, y = 3, tx = 1, ty = 1),
               const = list(x = rep(2, 5), y = rep(2, 5), tx = 2, ty = 2))
  both <- function(c) at_degrees(pairs_02, function(d, e) c(co_lpm(c, d, e), co_upm(c, d, e)))
  univ <- function(c) at_degrees(pairs_02, function(d, e) c(LPM(d + e, c$tx, c$x), UPM(d + e, c$tx, c$x)))
  expect_identical(over(edge, both), over(edge, univ))
})

test_that("DLPM_self: a sample never diverges from itself, D.LPM and D.UPM are exactly 0", {
  self_all <- c(exact_self, sample_self)
  divergent <- function(c) at_degrees(degree_pairs, function(l, u) c(d_lpm(c, l, u), d_upm(c, l, u)))
  expect_true(all(over(self_all, divergent) == 0))
})

# ---- CoMomentsMatrix.lean: degree_y ----------------------------------------

test_that("coLPM_R_none / coUPM_R_none: omitting degree_y is identical to degree_y = degree", {
  dflt <- function(c) unlist(lapply(0:3, function(d) {
    c(Co.LPM(d, c$x, c$y, c$tx, c$ty), Co.UPM(d, c$x, c$y, c$tx, c$ty))
  }))
  expl <- function(c) unlist(lapply(0:3, function(d) c(co_lpm(c, d, d), co_upm(c, d, d))))
  expect_identical(over(all_cases, dflt), over(all_cases, expl))
})

test_that("coLPM_R_some / coUPM_R_some: degree_y applies to y only (exact dyadic data)", {
  lower <- function(d, t, v) ifelse(v <= t, (t - v)^d, 0)
  upper <- function(d, t, v) ifelse(v > t, (v - t)^d, 0)
  kern <- function(c) at_degrees(pairs_03, function(d, e) c(co_lpm(c, d, e), co_upm(c, d, e)))
  oracle <- function(c) at_degrees(pairs_03, function(d, e) {
    c(mean(lower(d, c$tx, c$x) * lower(e, c$ty, c$y)),
      mean(upper(d, c$tx, c$x) * upper(e, c$ty, c$y)))
  })
  expect_identical(over(exact_cases, kern), over(exact_cases, oracle))
})

test_that("mixed_master: exact dyadic data, signed quadrants are identical to mean((x - tx)^d (y - ty)^e)", {
  lhs <- function(c) at_degrees(pairs_03, function(d, e) mixed_master_lhs(c, d, e))
  rhs <- function(c) at_degrees(pairs_03, function(d, e) co_moment(c, d, e))
  expect_identical(over(exact_cases, lhs), over(exact_cases, rhs))
})

test_that("mixed_master: random and real data, degrees 0..3, within n * eps", {
  expect_within(
    over(sample_cases, function(c) at_degrees(pairs_03, function(d, e) mixed_master_lhs(c, d, e))),
    over(sample_cases, function(c) at_degrees(pairs_03, function(d, e) co_moment(c, d, e))),
    over(sample_cases, function(c) at_degrees(pairs_03, function(d, e) co_tol(c, d, e))))
})

# Marginal identities: lhs(c, d), rhs(c, d), tol degrees (dx, dy) as functions of d.
marginal <- list(
  CoLPM_zero_add_DUPM = list(
    lhs = function(c, d) co_lpm(c, d, 0) + d_upm(c, d, 0),
    rhs = function(c, d) LPM(d, c$tx, c$x), tol = function(c, d) co_tol(c, d, 0)),
  CoUPM_zero_add_DLPM = list(
    lhs = function(c, d) co_upm(c, d, 0) + d_lpm(c, 0, d),
    rhs = function(c, d) UPM(d, c$tx, c$x), tol = function(c, d) co_tol(c, d, 0)),
  CoLPM_zero_left_add_DLPM = list(
    lhs = function(c, d) co_lpm(c, 0, d) + d_lpm(c, d, 0),
    rhs = function(c, d) LPM(d, c$ty, c$y), tol = function(c, d) co_tol(c, 0, d)),
  CoUPM_zero_left_add_DUPM = list(
    lhs = function(c, d) co_upm(c, 0, d) + d_upm(c, 0, d),
    rhs = function(c, d) UPM(d, c$ty, c$y), tol = function(c, d) co_tol(c, 0, d))
)
over_d <- function(cases, f) over(cases, function(c) vapply(0:3, function(d) f(c, d), numeric(1)))

test_that("CoLPM_zero_add_DUPM: Co.LPM(dx, ..., degree_y = 0) + D.UPM(dx, 0, ...) = LPM(dx, tx, x)", {
  m <- marginal$CoLPM_zero_add_DUPM
  expect_identical(over_d(exact_cases, m$lhs), over_d(exact_cases, m$rhs))
  expect_within(over_d(sample_cases, m$lhs), over_d(sample_cases, m$rhs), over_d(sample_cases, m$tol))
})

test_that("CoUPM_zero_add_DLPM: Co.UPM(dx, ..., degree_y = 0) + D.LPM(0, dx, ...) = UPM(dx, tx, x)", {
  m <- marginal$CoUPM_zero_add_DLPM
  expect_identical(over_d(exact_cases, m$lhs), over_d(exact_cases, m$rhs))
  expect_within(over_d(sample_cases, m$lhs), over_d(sample_cases, m$rhs), over_d(sample_cases, m$tol))
})

test_that("CoLPM_zero_left_add_DLPM: Co.LPM(0, ..., degree_y = dy) + D.LPM(dy, 0, ...) = LPM(dy, ty, y)", {
  m <- marginal$CoLPM_zero_left_add_DLPM
  expect_identical(over_d(exact_cases, m$lhs), over_d(exact_cases, m$rhs))
  expect_within(over_d(sample_cases, m$lhs), over_d(sample_cases, m$rhs), over_d(sample_cases, m$tol))
})

test_that("CoUPM_zero_left_add_DUPM: Co.UPM(0, ..., degree_y = dy) + D.UPM(0, dy, ...) = UPM(dy, ty, y)", {
  m <- marginal$CoUPM_zero_left_add_DUPM
  expect_identical(over_d(exact_cases, m$lhs), over_d(exact_cases, m$rhs))
  expect_within(over_d(sample_cases, m$lhs), over_d(sample_cases, m$rhs), over_d(sample_cases, m$tol))
})

test_that("CoLPM_zero_le_LPM / CoUPM_zero_le_UPM: degree_y = 0 is bounded by the univariate moment", {
  slack <- function(c, d) {
    c(LPM(d, c$tx, c$x) - co_lpm(c, d, 0), UPM(d, c$tx, c$x) - co_upm(c, d, 0)) + co_tol(c, d, 0)
  }
  expect_true(all(over(sample_cases, function(c) sapply(0:3, function(d) slack(c, d))) >= 0))
})

# Cauchy-Schwarz: bound^2 - value^2 >= 0 up to a relative rounding margin.
cs_slack <- function(c, d, e) {
  rel <- 1 + 8 * (length(c$x) + 2 * d + 2 * e + 2) * eps
  c(CoLPM = LPM(2 * d, c$tx, c$x) * LPM(2 * e, c$ty, c$y) * rel - co_lpm(c, d, e)^2,
    CoUPM = UPM(2 * d, c$tx, c$x) * UPM(2 * e, c$ty, c$y) * rel - co_upm(c, d, e)^2,
    DLPM  = UPM(2 * e, c$tx, c$x) * LPM(2 * d, c$ty, c$y) * rel - d_lpm(c, d, e)^2,
    DUPM  = LPM(2 * d, c$tx, c$x) * UPM(2 * e, c$ty, c$y) * rel - d_upm(c, d, e)^2)
}

test_that("CoLPM_sq_le, CoUPM_sq_le, DLPM_sq_le, DUPM_sq_le: Cauchy-Schwarz bounds", {
  slack <- over(all_cases, function(c) at_degrees(degree_pairs, function(d, e) cs_slack(c, d, e)))
  expect_true(all(slack >= 0), info = paste(names(slack)[slack < 0], collapse = ", "))
})

test_that("CoLPM_sq_le: equality at x = y (Co.LPM(d, x, x)^2 = LPM(2d)^2)", {
  x <- as.numeric(eu[, "SMI"])
  t <- median(x)
  expect_equal(vapply(1:2, function(d) Co.LPM(d, x, x, t, t)^2, numeric(1)),
               vapply(1:2, function(d) LPM(2 * d, t, x)^2, numeric(1)), tolerance = 1e-12)
})

# ---- vector targets --------------------------------------------------------

test_that("SPEC vector targets: 100 targets (parallel path) are identical to scalar calls", {
  c <- random_cases[["normal mean"]]
  ts <- seq(-2, 2, length.out = 100)
  scalar <- function(f) vapply(seq_along(ts), f, numeric(1))
  expect_identical(Co.LPM(1, c$x, c$y, ts, 0.3),
                   scalar(function(k) Co.LPM(1, c$x, c$y, ts[k], 0.3)))
  expect_identical(Co.UPM(2, c$x, c$y, 0.1, ts, 1),
                   scalar(function(k) Co.UPM(2, c$x, c$y, 0.1, ts[k], 1)))
  expect_identical(D.LPM(1, 2, c$x, c$y, ts, rev(ts)),
                   scalar(function(k) D.LPM(1, 2, c$x, c$y, ts[k], rev(ts)[k])))
  expect_identical(D.UPM(0, 1, c$x, c$y, ts, ts),
                   scalar(function(k) D.UPM(0, 1, c$x, c$y, ts[k], ts[k])))
})

eu_vec_40 <- local({
  x <- as.numeric(eu[, "CAC"]); y <- as.numeric(eu[, "SMI"])
  list(x = x, y = y,
       tx = c(quantile(x, seq(0.02, 0.98, length.out = 30), names = FALSE), x[1:10]),
       ty = c(quantile(y, seq(0.98, 0.02, length.out = 30), names = FALSE), y[1:10]))
})
split_targets <- function(c) {
  lapply(seq_along(c$tx), function(k) list(x = c$x, y = c$y, tx = c$tx[k], ty = c$ty[k]))
}

test_that("quadrants_sum_eq_one: holds elementwise for 40 vector targets", {
  expect_within(quadrant_sum0(eu_vec_40), 1, 4 * eps)
})

test_that("coMoment_decomposition: holds elementwise for 40 vector targets", {
  per_target <- split_targets(eu_vec_40)
  expect_within(decomposition(eu_vec_40), unname(over(per_target, co_moment)),
                unname(over(per_target, function(c) co_tol(c, 1, 1))))
})

# ---- CoMomentsMatrix.lean: PM.matrix ---------------------------------------

pm_scale <- function(pop_adj, n, dL, dU) {
  if (pop_adj && n > 1 && dL > 0 && dU > 0) n / (n - 1) else 1
}

pairwise <- function(f, X, tg) {
  k <- ncol(X)
  M <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) M[i, j] <- f(X[, i], X[, j], tg[i], tg[j])
  M
}

# PM.matrix quadrants rebuilt from the pairwise bivariate functions.
pm_reference <- function(dL, dU, tg, X, pop_adj) {
  s <- pm_scale(pop_adj, nrow(X), dL, dU)
  list(
    cupm = s * pairwise(function(a, b, u, v) Co.UPM(dU, a, b, u, v), X, tg),
    dupm = s * pairwise(function(a, b, u, v) D.UPM(dL, dU, a, b, u, v), X, tg),
    dlpm = s * pairwise(function(a, b, u, v) D.LPM(dL, dU, a, b, u, v), X, tg),
    clpm = s * pairwise(function(a, b, u, v) Co.LPM(dL, a, b, u, v), X, tg)
  )
}
pm_quadrants <- function(...) lapply(PM.matrix(...)[c("cupm", "dupm", "dlpm", "clpm")], unname)

test_that("clpmM / cupmM / dupmM / dlpmM (definitions): exact dyadic data, entries identical to the pairwise functions", {
  got <- lapply(degree_pairs, function(d) pm_quadrants(d[1], d[2], t_exact, X_exact, FALSE))
  ref <- lapply(degree_pairs, function(d) pm_reference(d[1], d[2], t_exact, X_exact, FALSE))
  expect_identical(got, ref)
})

test_that("clpmM / cupmM / dupmM / dlpmM (definitions): random and real data, mixed degrees, pop_adj scale", {
  grid <- expand.grid(sample = names(pm_samples), deg = seq_along(degree_pairs),
                      pop_adj = c(FALSE, TRUE), stringsAsFactors = FALSE)
  errs <- unlist(Map(function(s, k, pa) {
    X <- pm_samples[[s]]; tg <- colMeans(X); d <- degree_pairs[[k]]
    got <- pm_quadrants(d[1], d[2], tg, X, pa)
    ref <- pm_reference(d[1], d[2], tg, X, pa)
    vapply(names(ref), function(q) max(abs(got[[q]] - ref[[q]])) / mat_tol(ref[[q]], nrow(X)),
           numeric(1))
  }, grid$sample, grid$deg, grid$pop_adj))
  expect_within(errs, 0, 1)
})

test_that("PMmatrix_dupm_transpose: dupm is identical to t(dlpm), normalised or not", {
  for (X in pm_all) for (d in degree_pairs) for (nm in c(FALSE, TRUE)) {
    P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm)
    expect_identical(P$dupm, t(P$dlpm))
  }
})

test_that("PMmatrix_clpm_symm / PMmatrix_cupm_symm / PMmatrix_cov_symm: exactly symmetric", {
  for (X in pm_all) for (d in degree_pairs) for (nm in c(FALSE, TRUE)) {
    P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm)
    expect_identical(P[c("clpm", "cupm", "cov.matrix")],
                     lapply(P[c("clpm", "cupm", "cov.matrix")], t))
  }
})

test_that("dupmM_diag / dlpmM_diag: divergent diagonals are exactly 0", {
  diags <- unlist(lapply(pm_all, function(X) lapply(degree_pairs, function(d) {
    P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE)
    c(diag(P$dupm), diag(P$dlpm))
  })))
  expect_true(all(diags == 0))
})

test_that("clpmM_diag / cupmM_diag: diagonals are s * LPM(2 dL) and s * UPM(2 dU)", {
  for (X in pm_samples) for (d in degree_pairs) for (pa in c(FALSE, TRUE)) {
    tg <- colMeans(X); n <- nrow(X); s <- pm_scale(pa, n, d[1], d[2])
    P <- PM.matrix(d[1], d[2], tg, X, pa)
    lpm2 <- s * vapply(seq_len(ncol(X)), function(i) LPM(2 * d[1], tg[i], X[, i]), numeric(1))
    upm2 <- s * vapply(seq_len(ncol(X)), function(i) UPM(2 * d[2], tg[i], X[, i]), numeric(1))
    expect_within(c(unname(diag(P$clpm)), unname(diag(P$cupm))), c(lpm2, upm2),
                  rep(c(mat_tol(P$clpm, n), mat_tol(P$cupm, n)), each = ncol(X)))
  }
})

test_that("PMmatrix_cov_eq: cov.matrix = cupm + clpm - dupm - dlpm, normalised or not", {
  for (X in pm_all) for (d in degree_pairs) for (nm in c(FALSE, TRUE)) {
    P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm)
    expect_within(P$cov.matrix, P$cupm + P$clpm - P$dupm - P$dlpm,
                  8 * eps * max(abs(P$cupm) + abs(P$clpm) + abs(P$dupm) + abs(P$dlpm)))
  }
})

test_that("covM_one_apply: exact dyadic data, degree-1 cov.matrix is identical to the co-moment about the targets", {
  P <- unname(PM.matrix(1, 1, t_exact, X_exact, FALSE)$cov.matrix)
  expect_identical(P, pairwise(function(a, b, u, v) mean((a - u) * (b - v)), X_exact, t_exact))
})

test_that("covM_one_apply: random and real data, targets at the 30% quantiles", {
  for (X in pm_samples) {
    n <- nrow(X)
    tg <- apply(X, 2, quantile, probs = 0.3, names = FALSE)
    P <- unname(PM.matrix(1, 1, tg, X, FALSE)$cov.matrix)
    co <- pairwise(function(a, b, u, v) mean((a - u) * (b - v)), X, tg)
    tol <- 4 * (n + 4) * eps * pairwise(function(a, b, u, v) mean(abs((a - u) * (b - v))), X, tg)
    expect_within(P, co, tol)
  }
})

cov_tol <- function(X) 4 * (nrow(X) + 4) * eps * outer(apply(X, 2, sd), apply(X, 2, sd))

test_that("PMmatrix_cov_eq_popCov: PM.matrix(1, 1, NULL, X, FALSE)$cov.matrix is the population covariance", {
  for (X in pm_samples) {
    n <- nrow(X)
    expect_within(PM.matrix(1, 1, NULL, X, FALSE)$cov.matrix, cov(X) * (n - 1) / n, cov_tol(X))
  }
})

test_that("PMmatrix_cov_eq_sampleCov: PM.matrix(1, 1, NULL, X, TRUE)$cov.matrix is cov(X)", {
  for (X in pm_samples) {
    P <- PM.matrix(1, 1, NULL, X, TRUE)$cov.matrix
    expect_within(P, cov(X), cov_tol(X))
    expect_identical(dimnames(P), dimnames(cov(X)))
  }
})

test_that("pmScale_true / PMmatrix_cov_popAdj: pop_adj scales all quadrants by n / (n - 1) when both degrees > 0", {
  for (X in pm_samples) for (d in list(c(1, 1), c(2, 1), c(1, 3))) {
    tg <- colMeans(X); n <- nrow(X)
    P1 <- PM.matrix(d[1], d[2], tg, X, TRUE)
    P0 <- PM.matrix(d[1], d[2], tg, X, FALSE)
    expect_within(unlist(P1), unlist(P0) * n / (n - 1), 4 * eps * max(abs(unlist(P1))))
  }
})

test_that("PMmatrix_popAdj_ignored: pop_adj changes nothing when either degree is 0", {
  for (X in pm_all) for (d in list(c(0, 0), c(1, 0), c(0, 1), c(2, 0), c(0, 2))) {
    for (nm in c(FALSE, TRUE)) {
      expect_identical(PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm),
                       PM.matrix(d[1], d[2], colMeans(X), X, FALSE, nm))
    }
  }
})

test_that("PMmatrix_popAdj_ignored: pop_adj changes nothing when n = 1", {
  X1 <- matrix(c(3, 5, -1), 1, dimnames = list(NULL, c("a", "b", "c")))
  for (d in degree_pairs) for (nm in c(FALSE, TRUE)) {
    expect_identical(PM.matrix(d[1], d[2], c(1, 2, 0), X1, TRUE, nm),
                     PM.matrix(d[1], d[2], c(1, 2, 0), X1, FALSE, nm))
  }
})

test_that("clpmM_posSemidef / cupmM_posSemidef / covM_posSemidef: unnormalised matrices are PSD", {
  grid <- expand.grid(sample = names(pm_samples), target = 1:3, deg = seq_along(degree_pairs),
                      pop_adj = c(FALSE, TRUE), stringsAsFactors = FALSE)
  margins <- unlist(Map(function(s, k, g, pa) {
    X <- pm_samples[[s]]; d <- degree_pairs[[g]]
    tg <- list(colMeans(X), apply(X, 2, median), X[2, ])[[k]]
    P <- PM.matrix(d[1], d[2], tg, X, pa)
    vapply(P[c("clpm", "cupm", "cov.matrix")], function(M) {
      min_eigen(M) / (8 * (nrow(X) + ncol(M)) * eps * max(abs(M)))
    }, numeric(1))
  }, grid$sample, grid$target, grid$deg, grid$pop_adj))
  expect_true(all(margins >= -1), info = paste(names(margins)[margins < -1], collapse = ", "))
})

test_that("clpmM_quadForm / cupmM_quadForm / covM_quadForm: exact dyadic quadratic forms are sums of squares", {
  v <- c(1, -2, 0.5)
  for (d in list(c(1, 1), c(2, 1), c(0, 2))) {
    P <- PM.matrix(d[1], d[2], t_exact, X_exact, FALSE)
    Lm <- sweep(X_exact, 2, t_exact, function(x, t) ifelse(x <= t, (t - x)^d[1], 0))
    Um <- sweep(X_exact, 2, t_exact, function(x, t) ifelse(x > t, (x - t)^d[2], 0))
    expect_equal(
      c(v %*% P$clpm %*% v, v %*% P$cupm %*% v, v %*% P$cov.matrix %*% v),
      c(mean((Lm %*% v)^2), mean((Um %*% v)^2), mean(((Um - Lm) %*% v)^2)),
      tolerance = 1e-14)
  }
})

# ---- CoMomentsMatrix.lean: norm = TRUE -------------------------------------

norm_sum <- function(P) P$cupm + P$dupm + P$dlpm + P$clpm

test_that("PMmatrix_norm_sum_eq_one: normalised quadrants sum to 1 in every cell with positive total", {
  for (X in pm_samples) for (d in degree_pairs) {
    tg <- colMeans(X)
    expect_true(all(norm_sum(PM.matrix(d[1], d[2], tg, X, FALSE)) > 0))
    expect_within(norm_sum(PM.matrix(d[1], d[2], tg, X, FALSE, TRUE)), 1, 8 * eps)
  }
})

test_that("PMmatrix_norm_cov_mem_Icc / normQ_mem_Icc: normalised quadrants lie in [0, 1], cov in [-1, 1]", {
  for (X in pm_all) for (d in degree_pairs) for (pa in c(FALSE, TRUE)) {
    P <- PM.matrix(d[1], d[2], colMeans(X), X, pa, TRUE)
    q <- unlist(P[c("cupm", "dupm", "dlpm", "clpm")])
    expect_true(all(q >= 0 & q <= 1 + 2 * eps) && all(abs(P$cov.matrix) <= 1 + 4 * eps))
  }
})

test_that("PMmatrix_norm_popAdj: pop_adj cancels under norm = TRUE", {
  for (X in pm_samples) for (d in degree_pairs) {
    tg <- colMeans(X)
    expect_within(unlist(PM.matrix(d[1], d[2], tg, X, TRUE, TRUE)),
                  unlist(PM.matrix(d[1], d[2], tg, X, FALSE, TRUE)), 8 * eps)
  }
})

test_that("PMmatrix_zero_zero_sum: exact dyadic data, degree-0 quadrants sum to exactly 1, normalised or not", {
  for (nm in c(FALSE, TRUE)) {
    expect_identical(unname(norm_sum(PM.matrix(0, 0, t_exact, X_exact, FALSE, nm))), matrix(1, 3, 3))
  }
})

test_that("PMmatrix_zero_zero_sum: random and real data, degree-0 quadrants sum to 1 within 8 eps", {
  for (X in pm_samples) for (nm in c(FALSE, TRUE)) {
    expect_within(norm_sum(PM.matrix(0, 0, colMeans(X), X, TRUE, nm)), 1, 8 * eps)
  }
})

test_that("PMmatrix_zero_zero_norm: at degree 0 norm = TRUE changes nothing", {
  expect_identical(PM.matrix(0, 0, t_exact, X_exact, FALSE, TRUE),
                   PM.matrix(0, 0, t_exact, X_exact, FALSE, FALSE))
  for (X in pm_samples) {
    expect_within(unlist(PM.matrix(0, 0, colMeans(X), X, TRUE, TRUE)),
                  unlist(PM.matrix(0, 0, colMeans(X), X, TRUE, FALSE)), 8 * eps)
  }
})

X_tied <- cbind(a = c(1, 2, 3, 4), b = c(5, 5, 5, 5), c = c(0, 2, -1, 3))
t_tied <- c(2.5, 5, 1)

# A zero cell total (dL > 0 and every row ties the target in column i or j)
# takes the degree-0 cells, as LPM.ratio and the normalised nD functions do.
# CoMomentsMatrix.lean's PMmatrix_norm_of_total_eq_zero, PMmatrix_norm_sum_of_tied
# and PMmatrix_norm_sum_eq_zero_example model the earlier all-zero convention.

test_that("PMmatrix_norm_of_total_eq_zero (degree-0 convention): a tied column takes the degree-0 cells (dL > 0)", {
  P0 <- PM.matrix(0, 0, t_tied, X_tied, FALSE, FALSE)
  for (d in list(c(1, 1), c(2, 1), c(1, 0))) {
    P <- PM.matrix(d[1], d[2], t_tied, X_tied, FALSE, TRUE)
    for (q in names(P)) {
      expect_within(P[[q]]["b", ], P0[[q]]["b", ], 4 * eps)
      expect_within(P[[q]][, "b"], P0[[q]][, "b"], 4 * eps)
    }
    expect_within(norm_sum(P), 1, 8 * eps)
  }
})

test_that("PMmatrix_norm_cov_diag: dL > 0 gives diagonal 1, also for a tied column", {
  P <- PM.matrix(1, 1, NULL, cbind(a = c(1, 2, 3), b = c(5, 5, 5)), FALSE, TRUE)
  expect_within(unname(P$cov.matrix), matrix(c(1, 1 / 3, 1 / 3, 1), 2), 4 * eps)
  P <- PM.matrix(2, 1, t_tied, X_tied, TRUE, TRUE)
  expect_within(unname(diag(P$cov.matrix)), c(1, 1, 1), 4 * eps)
})

test_that("PMmatrix_norm_sum_eq_zero_example (degree-0 convention): PM.matrix(1, 1, 0, matrix(0), TRUE, TRUE) is all lower", {
  P <- PM.matrix(1, 1, 0, matrix(0), TRUE, TRUE)
  expect_identical(vapply(P, function(M) M[1, 1], numeric(1)),
                   c(cupm = 0, dupm = 0, dlpm = 0, clpm = 1, cov.matrix = 1))
})

test_that("PMmatrix_norm_cov_diag_zero_left: at dL = 0 the normalised cov diagonal is 1, even for a tied column", {
  P <- PM.matrix(0, 1, NULL, cbind(a = 1:3, b = c(5, 5, 5)), FALSE, TRUE)
  expect_within(unname(P$cov.matrix), matrix(c(1, 1 / 3, 1 / 3, 1), 2), 4 * eps)
  for (dU in 0:2) {
    expect_within(unname(diag(PM.matrix(0, dU, t_tied, X_tied, TRUE, TRUE)$cov.matrix)), 1, 4 * eps)
  }
})

test_that("PMmatrix_norm_sum_eq_one_zero_left / pmTotal_zero_left_pos: at dL = 0 every cell sums to 1", {
  for (dU in 0:2) {
    expect_within(norm_sum(PM.matrix(0, dU, t_tied, X_tied, TRUE, TRUE)), 1, 8 * eps)
  }
})

test_that("pmTotal_of_tied_zero_left: a tied column has diagonal total exactly 1 at dL = 0", {
  totals <- vapply(0:2, function(dU) norm_sum(PM.matrix(0, dU, t_tied, X_tied, TRUE, FALSE))["b", "b"],
                   numeric(1))
  expect_identical(totals, c(1, 1, 1))
})

test_that("PMmatrix_norm_cov_not_posSemidef: the normalised cov.matrix is not PSD (target 0)", {
  X <- rbind(c(1, 0, 3), c(3, 0, 1), c(2, 0, 0), c(0, 3, 3))
  M <- PM.matrix(1, 1, rep(0, 3), X, FALSE, TRUE)$cov.matrix
  # Cell (1, 2) has total 0 (every row ties in column 1 or 2) and takes the
  # degree-0 value -1; the Lean counterexample uses the earlier 0. Not PSD
  # either way.
  expect_within(M, matrix(c(1, -1, 1, -1, 1, 1, 1, 1, 1), 3), 4 * eps)
  v <- c(1, 1, -1)
  expect_lt(drop(v %*% M %*% v), 0)
})

test_that("PMmatrix_norm_cov_not_posSemidef_mean: not PSD with column-mean targets and positive totals", {
  Z <- rbind(c(-1, -1, 2), c(2, 0, 2), c(1, 2, -2), c(-2, 1, -1))
  M <- PM.matrix(1, 1, NULL, Z, FALSE, TRUE)$cov.matrix
  expect_within(M, matrix(c(1, .2, .2, .2, 1, -1, .2, -1, 1), 3), 1e-14)
  v <- c(-2, 5, 5)
  expect_within(drop(v %*% M %*% v), -4, 1e-12)
  expect_lt(min_eigen(M), -0.07)
})

test_that("PMmatrix_norm_clpm_not_posSemidef: the normalised clpm is not PSD", {
  W <- rbind(c(-2, -2, 0), c(-2, 1, -2), c(-2, 0, 1))
  M <- PM.matrix(1, 1, rep(0, 3), W, FALSE, TRUE)$clpm
  expect_within(M, matrix(c(1, 2 / 3, 2 / 3, 2 / 3, .8, 0, 2 / 3, 0, .8), 3), 1e-14)
  v <- c(-4, 3, 3)
  expect_within(drop(v %*% M %*% v), -8 / 5, 1e-12)
})

# ---- nD functions (no Lean theorem; SPEC.md "Co.LPM_nD", "DPM_nD") ---------

nd_sum <- function(X, t, d, norm) {
  Co.LPM_nD(X, t, d, norm) + Co.UPM_nD(X, t, d, norm) + DPM_nD(X, t, d, norm)
}

test_that("SPEC DPM_nD (nD analogue of quadrants_sum_eq_one): exact dyadic data, degree-0 cells sum to exactly 1", {
  expect_identical(c(nd_sum(X_exact, t_exact, 0, FALSE), nd_sum(X_exact[, 1:2], t_exact[1:2], 0, FALSE)),
                   c(1, 1))
})

test_that("SPEC DPM_nD (nD analogue of quadrants_sum_eq_one): random and real data, degree-0 cells partition the sample", {
  sums <- unlist(lapply(pm_samples, function(X) {
    lapply(list(colMeans(X), apply(X, 2, median), X[3, ]), function(tg) {
      c(nd_sum(X, tg, 0, FALSE), nd_sum(X[, 1:2], tg[1:2], 0, TRUE))
    })
  }))
  expect_within(sums, 1, 4 * eps)
})

# Independent oracle for the k-column functions (SPEC.md "Co.LPM_nD",
# "DPM_nD"): with D = X - t row-wise, Co.LPM_nD = mean prod_j (t_j - x_ij)_+^d,
# Co.UPM_nD = mean prod_j (x_ij - t_j)_+^d, DPM_nD = mean over rows that are
# neither strictly all-below nor strictly all-above of prod_j |x_ij - t_j|^d.
# At degree 0 the cells are the indicators all(D <= 0), all(D > 0) and
# neither. ts attributes are dropped first (pmax on an mts fails).
nd_oracle <- function(X, t, d) {
  X <- matrix(as.numeric(X), nrow(X), ncol(X))
  D <- sweep(X, 2, t)
  rowprod <- function(M) apply(M, 1, prod)
  if (d == 0) {
    lo <- apply(D <= 0, 1, all)
    hi <- apply(D > 0, 1, all)
    return(c(clpm = mean(lo), cupm = mean(hi), dpm = mean(!lo & !hi)))
  }
  strict <- apply(D < 0, 1, all) | apply(D > 0, 1, all)
  c(clpm = mean(rowprod(pmax(-D, 0)^d)), cupm = mean(rowprod(pmax(D, 0)^d)),
    dpm = mean(ifelse(strict, 0, rowprod(abs(D)^d))))
}
nd_triple <- function(X, t, d) {
  c(clpm = Co.LPM_nD(X, t, d, FALSE), cupm = Co.UPM_nD(X, t, d, FALSE), dpm = DPM_nD(X, t, d, FALSE))
}
# Scale of a k-column moment: mean prod_j |x_ij - t_j|^d (the three cells' total).
nd_tol <- function(X, t, d, k = 4) {
  X <- matrix(as.numeric(X), nrow(X), ncol(X))
  k * (nrow(X) + d * ncol(X) + 2) * eps * mean(apply(abs(sweep(X, 2, t))^d, 1, prod))
}

# t_exact ties two rows in column c, so Co.LPM_nD is 0 for d > 0 there (zero
# branch only); c(0.5, 1, 0.5) makes all three cells positive at every degree
# with five ties; c(1, 1, 1) and c(0, 0, 0) put one of the co-cells at 0.
nd_exact_targets <- list(t_exact = t_exact, mixed = c(0.5, 1, 0.5),
                         ones = c(1, 1, 1), zeros = c(0, 0, 0))

test_that("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: k = 3 exact dyadic data, identical to the k-column product oracle", {
  per_target <- function(f) unlist(lapply(nd_exact_targets, function(t) {
    unlist(lapply(0:2, function(d) setNames(f(X_exact, t, d), paste0("d", d, ".", c("clpm", "cupm", "dpm")))))
  }))
  got <- per_target(nd_triple)
  ref <- per_target(nd_oracle)
  expect_identical(got, ref)
  # the mixed target really exercises all three cells at d > 0
  mixed_pos <- got[grep("^mixed\\.d[12]\\.", names(got))]
  expect_true(length(mixed_pos) == 6 && all(mixed_pos > 0))
})

test_that("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: k = 3..4, random and real data, within n * eps of the oracle", {
  vals <- function(f) unlist(lapply(pm_samples, function(X) {
    lapply(list(mean = colMeans(X), obs = X[3, ]), function(tg) {
      unlist(lapply(0:2, function(d) setNames(f(X, unname(tg), d), paste0("d", d, ".", c("clpm", "cupm", "dpm")))))
    })
  }))
  expect_within(vals(nd_triple), vals(nd_oracle), vals(function(X, t, d) rep(nd_tol(X, t, d), 3)))
})

test_that("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: n = 1 rows (k = 3)", {
  t <- c(2, 2, 2)
  expect_identical(c(Co.LPM_nD(rbind(c(1, 2, 3)), t, 1, FALSE),
                     Co.UPM_nD(rbind(c(3, 4, 5)), t, 1, FALSE),
                     DPM_nD(rbind(c(1, 3, 5)), t, 1, FALSE),
                     Co.LPM_nD(rbind(c(0, 1, 1.5)), t, 2, FALSE)),
                   c(0, 6, 3, 2 * 2 * 1 * 1 * 0.5 * 0.5))
  # a tie in one coordinate: 0^0 = 1 at degree 0 but 0^d = 0 for d > 0
  expect_identical(c(Co.LPM_nD(rbind(c(1, 2, 3)), t, 0, FALSE),
                     DPM_nD(rbind(c(1, 2, 3)), t, 0, FALSE),
                     DPM_nD(rbind(c(1, 2, 3)), t, 1, FALSE),
                     Co.LPM_nD(rbind(c(1, 2, 1)), t, 0, FALSE),
                     Co.LPM_nD(rbind(c(1, 2, 1)), t, 1, FALSE)),
                   c(0, 1, 0, 1, 0))
})

test_that("SPEC Co.LPM_nD / DPM_nD: a sample constant at its target, raw cells", {
  X <- matrix(1, 3, 2)
  expect_identical(c(nd_triple(X, c(1, 1), 0), nd_triple(X, c(1, 1), 1)),
                   c(clpm = 1, cupm = 0, dpm = 0, clpm = 0, cupm = 0, dpm = 0))
})

test_that("SPEC Co.LPM_nD: norm = TRUE cells sum to 1 for a sample constant at its target (degree > 0)", {
  # fails on 8898017: see BUGS.md #18 (normaliser 0 -> every norm = TRUE cell
  # returns 0, so the sum is 0 instead of 1)
  sums <- vapply(1:2, function(d) nd_sum(matrix(1, 3, 2), c(1, 1), d, TRUE), numeric(1))
  expect_within(sums, 1, 8 * eps)
})

nd_cases <- function(f) over(all_cases, function(c) vapply(0:2, function(d) f(c, d), numeric(1)))

test_that("SPEC Co.LPM_nD: k = 2 is identical to Co.LPM", {
  expect_identical(nd_cases(function(c, d) Co.LPM_nD(cbind(c$x, c$y), c(c$tx, c$ty), d, FALSE)),
                   nd_cases(function(c, d) co_lpm(c, d, d)))
})

test_that("SPEC Co.UPM_nD: k = 2 reduces to Co.UPM", {
  expect_within(nd_cases(function(c, d) Co.UPM_nD(cbind(c$x, c$y), c(c$tx, c$ty), d, FALSE)),
                nd_cases(function(c, d) co_upm(c, d, d)),
                nd_cases(function(c, d) co_tol(c, d, d, k = 2)))
})

test_that("SPEC DPM_nD: k = 2 reduces to D.LPM(d, d) + D.UPM(d, d)", {
  expect_within(nd_cases(function(c, d) DPM_nD(cbind(c$x, c$y), c(c$tx, c$ty), d, FALSE)),
                nd_cases(function(c, d) d_lpm(c, d, d) + d_upm(c, d, d)),
                nd_cases(function(c, d) co_tol(c, d, d)))
})

test_that("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: one column reduces to LPM / UPM / 0", {
  one_col <- function(c, d) {
    c(Co.LPM_nD(matrix(c$x), c$tx, d, FALSE), Co.UPM_nD(matrix(c$x), c$tx, d, FALSE),
      DPM_nD(matrix(c$x), c$tx, d, FALSE))
  }
  univ <- function(c, d) c(LPM(d, c$tx, c$x), UPM(d, c$tx, c$x), 0)
  tol <- function(c, d) rep(co_tol(list(x = c$x, y = c$x, tx = c$tx, ty = c$tx), d, 0), 3)
  over3 <- function(f) over(all_cases, function(c) unlist(lapply(0:2, function(d) f(c, d))))
  expect_within(over3(one_col), over3(univ), over3(tol))
})

test_that("SPEC Co.LPM_nD: norm = TRUE cells sum to 1 when the total is positive (degree > 0)", {
  sums <- unlist(lapply(pm_samples, function(X) vapply(1:2, function(d) nd_sum(X, colMeans(X), d, TRUE),
                                                       numeric(1))))
  expect_within(sums, 1, 8 * eps)
})
