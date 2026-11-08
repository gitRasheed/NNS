# Co- and divergent partial moments, PM.matrix and the nD partial moments. Lean:
# nns-proofs CoPartialMoments.lean, CoMomentsMatrix.lean; the nD functions cite SPEC.md.

# Lean: coMoment_decomposition, covariance_decomposition | BUGS: none
test_that("C01: co-moment quadrants decompose the co-moment and the covariance, for scalar and vector targets", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("coMoment_decomposition: exact dyadic data, CoUPM + CoLPM - DLPM - DUPM is identical to the co-moment", {
    expect_identical(aud_over(fc$exact_cases, aud_decomposition), aud_over(fc$exact_cases, aud_co_moment))
  })
  aud_case("coMoment_decomposition: random and real data, every target, within n * eps", {
    aud_within_cells(aud_over(fc$sample_cases, aud_decomposition), aud_over(fc$sample_cases, aud_co_moment),
                     aud_over(fc$sample_cases, function(c) aud_co_tol(c, 1, 1)))
  })
  aud_case("coMoment_decomposition: n = 1 is exact", {
    expect_identical(aud_decomposition(list(x = 3, y = 5, tx = 1, ty = 2)), 6)
    expect_identical(aud_decomposition(list(x = 3, y = -5, tx = 1, ty = 2)), -14)
  })
  aud_case("covariance_decomposition: targets at the means give the population covariance", {
    mean_cases <- fc$sample_cases[grep("mean$", names(fc$sample_cases))]
    pop_cov <- function(c) cov(c$x, c$y) * (length(c$x) - 1) / length(c$x)
    aud_within_cells(aud_over(mean_cases, aud_decomposition), aud_over(mean_cases, pop_cov),
                     aud_over(mean_cases, function(c) aud_co_tol(c, 1, 1)))
  })
  aud_case("covariance_decomposition: target_x = target_y = NULL resolves to the means (bivariate mean-default path)", {
    # one case per data set; its numeric targets are mean(x), mean(y), so co_tol
    # is evaluated at the targets NULL resolves to
    mean_cases <- fc$sample_cases[grep("mean$", names(fc$sample_cases))]
    null_t <- function(c) list(x = c$x, y = c$y, tx = NULL, ty = NULL)
    pop_cov <- function(c) cov(c$x, c$y) * (length(c$x) - 1) / length(c$x)
    aud_within_cells(aud_over(mean_cases, function(c) aud_decomposition(null_t(c))), aud_over(mean_cases, pop_cov),
                     aud_over(mean_cases, function(c) aud_co_tol(c, 1, 1)))
    aud_within_cells(aud_over(mean_cases, function(c) aud_quadrant_sum0(null_t(c))), 1, 4 * aud_eps)
  })
  aud_case("SPEC vector targets: 100 targets (parallel path) are identical to scalar calls", {
    c <- fc$random_cases[["normal mean"]]
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
  aud_case("coMoment_decomposition: holds elementwise for 40 vector targets", {
    per_target <- aud_split_targets(fc$eu_vec_40)
    aud_within_cells(aud_decomposition(fc$eu_vec_40), unname(aud_over(per_target, aud_co_moment)),
                     unname(aud_over(per_target, function(c) aud_co_tol(c, 1, 1))))
  })
})

# Lean: quadrants_sum_eq_one | BUGS: a08d04d
test_that("C02: degree-0 quadrants partition the sample, ties in the co-lower quadrant", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("D.LPM and D.UPM count degree-0 ties on the documented side", {
    x <- c(-1, 0, 0, 1, 2, 0)
    y <- c(0, 0, 1, -1, 0, 2)
    expect_identical(D.LPM(0, 0, x, y, 0, 0), mean(x > 0 & y <= 0))
    expect_identical(D.UPM(0, 0, x, y, 0, 0), mean(x <= 0 & y > 0))
  })
  aud_case("quadrants_sum_eq_one: exact dyadic data, the degree-0 quadrants sum to exactly 1", {
    expect_identical(unname(aud_over(fc$exact_cases, aud_quadrant_sum0)), rep(1, length(fc$exact_cases)))
  })
  aud_case("quadrants_sum_eq_one: random and real data, ties at the target, within 4 eps", {
    aud_within_cells(aud_over(fc$sample_cases, aud_quadrant_sum0), 1, 4 * aud_eps)
  })
  aud_case("quadrants_sum_eq_one: n = 1, and a constant sample at its target sits in the co-lower quadrant", {
    expect_identical(aud_quadrant_sum0(list(x = 3, y = 5, tx = 3, ty = 5)), 1)
    expect_identical(aud_quadrant_sum0(list(x = 3, y = 5, tx = 4, ty = 4)), 1)
    const <- list(x = rep(2, 7), y = rep(-1, 7), tx = 2, ty = -1)
    expect_identical(c(aud_co_lpm(const, 0, 0), aud_co_upm(const, 0, 0),
                       aud_d_lpm(const, 0, 0), aud_d_upm(const, 0, 0)), c(1, 0, 0, 0))
  })
  aud_case("quadrants_sum_eq_one: degree 0 vs degree > 0 at full ties (0^0 = 1 but 0^d = 0)", {
    const <- list(x = rep(2, 7), y = rep(-1, 7), tx = 2, ty = -1)
    expect_identical(c(aud_co_lpm(const, 1, 1), aud_co_upm(const, 1, 1), aud_d_lpm(const, 1, 1),
                       aud_d_upm(const, 1, 1), aud_co_lpm(const, 0, 1), aud_co_lpm(const, 1, 0)),
                     numeric(6))
  })
  aud_case("quadrants_sum_eq_one: holds elementwise for 40 vector targets", {
    aud_within_cells(aud_quadrant_sum0(fc$eu_vec_40), 1, 4 * aud_eps)
  })
})

# Lean: mixed_master | BUGS: none
test_that("C03: the mixed-degree master identity holds and degree_y applies to y only", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("coLPM_R_none / coUPM_R_none: omitting degree_y is identical to degree_y = degree", {
    dflt <- function(c) unlist(lapply(0:3, function(d) {
      c(Co.LPM(d, c$x, c$y, c$tx, c$ty), Co.UPM(d, c$x, c$y, c$tx, c$ty))
    }))
    expl <- function(c) unlist(lapply(0:3, function(d) c(aud_co_lpm(c, d, d), aud_co_upm(c, d, d))))
    expect_identical(aud_over(fc$all_cases, dflt), aud_over(fc$all_cases, expl))
  })
  aud_case("coLPM_R_some / coUPM_R_some: degree_y applies to y only (exact dyadic data)", {
    lower <- function(d, t, v) ifelse(v <= t, (t - v)^d, 0)
    upper <- function(d, t, v) ifelse(v > t, (v - t)^d, 0)
    kern <- function(c) aud_at_degrees(fc$pairs_03, function(d, e) c(aud_co_lpm(c, d, e), aud_co_upm(c, d, e)))
    oracle <- function(c) aud_at_degrees(fc$pairs_03, function(d, e) {
      c(mean(lower(d, c$tx, c$x) * lower(e, c$ty, c$y)),
        mean(upper(d, c$tx, c$x) * upper(e, c$ty, c$y)))
    })
    expect_identical(aud_over(fc$exact_cases, kern), aud_over(fc$exact_cases, oracle))
  })
  aud_case("mixed_master: exact dyadic data, signed quadrants are identical to mean((x - tx)^d (y - ty)^e)", {
    lhs <- function(c) aud_at_degrees(fc$pairs_03, function(d, e) aud_mixed_master_lhs(c, d, e))
    rhs <- function(c) aud_at_degrees(fc$pairs_03, function(d, e) aud_co_moment(c, d, e))
    expect_identical(aud_over(fc$exact_cases, lhs), aud_over(fc$exact_cases, rhs))
  })
  aud_case("mixed_master: random and real data, degrees 0..3, within n * eps", {
    aud_within_cells(
      aud_over(fc$sample_cases, function(c) aud_at_degrees(fc$pairs_03, function(d, e) aud_mixed_master_lhs(c, d, e))),
      aud_over(fc$sample_cases, function(c) aud_at_degrees(fc$pairs_03, function(d, e) aud_co_moment(c, d, e))),
      aud_over(fc$sample_cases, function(c) aud_at_degrees(fc$pairs_03, function(d, e) aud_co_tol(c, d, e))))
  })
})

# Lean: DLPM_swap, CoLPM_swap, CoUPM_swap, DLPM_self | BUGS: none
test_that("C04: co-moments are symmetric under swapping the variables, and a sample never diverges from itself", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("DLPM_swap: D.LPM(dl, du, x, y, tx, ty) is identical to D.UPM(dl, du, y, x, ty, tx)", {
    lhs <- function(c) aud_at_degrees(fc$degree_pairs, function(l, u) D.LPM(l, u, c$x, c$y, c$tx, c$ty))
    rhs <- function(c) aud_at_degrees(fc$degree_pairs, function(l, u) D.UPM(l, u, c$y, c$x, c$ty, c$tx))
    expect_identical(aud_over(fc$all_cases, lhs), aud_over(fc$all_cases, rhs))
  })
  aud_case("CoLPM_swap: Co.LPM(dx, x, y, tx, ty, dy) is identical to Co.LPM(dy, y, x, ty, tx, dx)", {
    lhs <- function(c) aud_at_degrees(fc$degree_pairs, function(d, e) Co.LPM(d, c$x, c$y, c$tx, c$ty, e))
    rhs <- function(c) aud_at_degrees(fc$degree_pairs, function(d, e) Co.LPM(e, c$y, c$x, c$ty, c$tx, d))
    expect_identical(aud_over(fc$all_cases, lhs), aud_over(fc$all_cases, rhs))
  })
  aud_case("CoUPM_swap: Co.UPM(dx, x, y, tx, ty, dy) is identical to Co.UPM(dy, y, x, ty, tx, dx)", {
    lhs <- function(c) aud_at_degrees(fc$degree_pairs, function(d, e) Co.UPM(d, c$x, c$y, c$tx, c$ty, e))
    rhs <- function(c) aud_at_degrees(fc$degree_pairs, function(d, e) Co.UPM(e, c$y, c$x, c$ty, c$tx, d))
    expect_identical(aud_over(fc$all_cases, lhs), aud_over(fc$all_cases, rhs))
  })
  aud_case("DLPM_self: a sample never diverges from itself, D.LPM and D.UPM are exactly 0", {
    self_all <- c(fc$exact_self, fc$sample_self)
    divergent <- function(c) aud_at_degrees(fc$degree_pairs, function(l, u) c(aud_d_lpm(c, l, u), aud_d_upm(c, l, u)))
    expect_true(all(aud_over(self_all, divergent) == 0))
  })
})

# Lean: CoLPM_self, CoUPM_self | BUGS: none
test_that("C05: the co-moment of a sample with itself is the univariate moment of the summed degree", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("CoLPM_self: exact dyadic data, Co.LPM(d, x, x, t, t, e) is identical to LPM(d + e, t, x)", {
    lhs <- function(c) aud_at_degrees(fc$pairs_02, function(d, e) aud_co_lpm(c, d, e))
    rhs <- function(c) aud_at_degrees(fc$pairs_02, function(d, e) LPM(d + e, c$tx, c$x))
    expect_identical(aud_over(fc$exact_self, lhs), aud_over(fc$exact_self, rhs))
  })
  aud_case("CoUPM_self: exact dyadic data, Co.UPM(d, x, x, t, t, e) is identical to UPM(d + e, t, x)", {
    lhs <- function(c) aud_at_degrees(fc$pairs_02, function(d, e) aud_co_upm(c, d, e))
    rhs <- function(c) aud_at_degrees(fc$pairs_02, function(d, e) UPM(d + e, c$tx, c$x))
    expect_identical(aud_over(fc$exact_self, lhs), aud_over(fc$exact_self, rhs))
  })
  aud_case("CoLPM_self: random and real data at tied targets, within n * eps", {
    aud_within_cells(
      aud_over(fc$sample_self, function(c) aud_at_degrees(fc$pairs_02, function(d, e) aud_co_lpm(c, d, e))),
      aud_over(fc$sample_self, function(c) aud_at_degrees(fc$pairs_02, function(d, e) LPM(d + e, c$tx, c$x))),
      aud_over(fc$sample_self, function(c) aud_at_degrees(fc$pairs_02, function(d, e) aud_co_tol(c, d, e))))
  })
  aud_case("CoUPM_self: random and real data at tied targets, within n * eps", {
    aud_within_cells(
      aud_over(fc$sample_self, function(c) aud_at_degrees(fc$pairs_02, function(d, e) aud_co_upm(c, d, e))),
      aud_over(fc$sample_self, function(c) aud_at_degrees(fc$pairs_02, function(d, e) UPM(d + e, c$tx, c$x))),
      aud_over(fc$sample_self, function(c) aud_at_degrees(fc$pairs_02, function(d, e) aud_co_tol(c, d, e))))
  })
  aud_case("CoLPM_self: 40 targets (LPM prefix backend) agree with the co-moment kernel", {
    c <- fc$eu_self_40
    aud_within_cells(aud_at_degrees(fc$pairs_02, function(d, e) aud_co_lpm(c, d, e)),
                     aud_at_degrees(fc$pairs_02, function(d, e) LPM(d + e, c$tx, c$x)),
                     aud_at_degrees(fc$pairs_02, function(d, e) aud_self_tol_40(c, d, e)))
  })
  aud_case("CoUPM_self: 40 targets (UPM prefix backend) agree with the co-moment kernel", {
    c <- fc$eu_self_40
    aud_within_cells(aud_at_degrees(fc$pairs_02, function(d, e) aud_co_upm(c, d, e)),
                     aud_at_degrees(fc$pairs_02, function(d, e) UPM(d + e, c$tx, c$x)),
                     aud_at_degrees(fc$pairs_02, function(d, e) aud_self_tol_40(c, d, e)))
  })
  aud_case("CoLPM_self / CoUPM_self: n = 1 and constant samples", {
    edge <- list(n1_low = list(x = 3, y = 3, tx = 5, ty = 5),
                 n1_high = list(x = 3, y = 3, tx = 1, ty = 1),
                 const = list(x = rep(2, 5), y = rep(2, 5), tx = 2, ty = 2))
    both <- function(c) aud_at_degrees(fc$pairs_02, function(d, e) c(aud_co_lpm(c, d, e), aud_co_upm(c, d, e)))
    univ <- function(c) aud_at_degrees(fc$pairs_02, function(d, e) c(LPM(d + e, c$tx, c$x), UPM(d + e, c$tx, c$x)))
    expect_identical(aud_over(edge, both), aud_over(edge, univ))
  })
})

# Lean: DLPM_nonneg, DUPM_nonneg, CoLPM_nonneg, CoUPM_nonneg, CoLPM_zero_add_DUPM, CoUPM_zero_add_DLPM, CoLPM_zero_left_add_DLPM, CoUPM_zero_left_add_DUPM, CoLPM_zero_le_LPM, CoUPM_zero_le_UPM, CoLPM_sq_le, CoUPM_sq_le, DLPM_sq_le, DUPM_sq_le | BUGS: none
test_that("C06: quadrants are nonnegative, marginalise at degree 0 and obey Cauchy-Schwarz", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("CoLPM_nonneg, CoUPM_nonneg, DLPM_nonneg, DUPM_nonneg: every quadrant is >= 0", {
    quads <- function(c) aud_at_degrees(fc$degree_pairs, function(d, e) {
      c(aud_co_lpm(c, d, e), aud_co_upm(c, d, e), aud_d_lpm(c, d, e), aud_d_upm(c, d, e))
    })
    expect_true(all(aud_over(fc$all_cases, quads) >= 0))
  })
  aud_case("CoLPM_zero_add_DUPM: Co.LPM(dx, ..., degree_y = 0) + D.UPM(dx, 0, ...) = LPM(dx, tx, x)", {
    m <- fc$marginal$CoLPM_zero_add_DUPM
    expect_identical(aud_over_d(fc$exact_cases, m$lhs), aud_over_d(fc$exact_cases, m$rhs))
    aud_within_cells(aud_over_d(fc$sample_cases, m$lhs), aud_over_d(fc$sample_cases, m$rhs), aud_over_d(fc$sample_cases, m$tol))
  })
  aud_case("CoUPM_zero_add_DLPM: Co.UPM(dx, ..., degree_y = 0) + D.LPM(0, dx, ...) = UPM(dx, tx, x)", {
    m <- fc$marginal$CoUPM_zero_add_DLPM
    expect_identical(aud_over_d(fc$exact_cases, m$lhs), aud_over_d(fc$exact_cases, m$rhs))
    aud_within_cells(aud_over_d(fc$sample_cases, m$lhs), aud_over_d(fc$sample_cases, m$rhs), aud_over_d(fc$sample_cases, m$tol))
  })
  aud_case("CoLPM_zero_left_add_DLPM: Co.LPM(0, ..., degree_y = dy) + D.LPM(dy, 0, ...) = LPM(dy, ty, y)", {
    m <- fc$marginal$CoLPM_zero_left_add_DLPM
    expect_identical(aud_over_d(fc$exact_cases, m$lhs), aud_over_d(fc$exact_cases, m$rhs))
    aud_within_cells(aud_over_d(fc$sample_cases, m$lhs), aud_over_d(fc$sample_cases, m$rhs), aud_over_d(fc$sample_cases, m$tol))
  })
  aud_case("CoUPM_zero_left_add_DUPM: Co.UPM(0, ..., degree_y = dy) + D.UPM(0, dy, ...) = UPM(dy, ty, y)", {
    m <- fc$marginal$CoUPM_zero_left_add_DUPM
    expect_identical(aud_over_d(fc$exact_cases, m$lhs), aud_over_d(fc$exact_cases, m$rhs))
    aud_within_cells(aud_over_d(fc$sample_cases, m$lhs), aud_over_d(fc$sample_cases, m$rhs), aud_over_d(fc$sample_cases, m$tol))
  })
  aud_case("CoLPM_zero_le_LPM / CoUPM_zero_le_UPM: degree_y = 0 is bounded by the univariate moment", {
    slack <- function(c, d) {
      c(LPM(d, c$tx, c$x) - aud_co_lpm(c, d, 0), UPM(d, c$tx, c$x) - aud_co_upm(c, d, 0)) + aud_co_tol(c, d, 0)
    }
    expect_true(all(aud_over(fc$sample_cases, function(c) sapply(0:3, function(d) slack(c, d))) >= 0))
  })
  aud_case("CoLPM_sq_le, CoUPM_sq_le, DLPM_sq_le, DUPM_sq_le: Cauchy-Schwarz bounds", {
    slack <- aud_over(fc$all_cases, function(c) aud_at_degrees(fc$degree_pairs, function(d, e) aud_cs_slack(c, d, e)))
    expect_true(all(slack >= 0), info = paste(names(slack)[slack < 0], collapse = ", "))
  })
  aud_case("CoLPM_sq_le: equality at x = y (Co.LPM(d, x, x)^2 = LPM(2d)^2)", {
    x <- as.numeric(fc$eu[, "SMI"])
    t <- median(x)
    expect_equal(vapply(1:2, function(d) Co.LPM(d, x, x, t, t)^2, numeric(1)),
                 vapply(1:2, function(d) LPM(2 * d, t, x)^2, numeric(1)), tolerance = 1e-12)
  })
})

# Lean: PMmatrix_dupm_transpose, PMmatrix_clpm_symm, PMmatrix_cupm_symm, PMmatrix_cov_symm, dupmM_diag, dlpmM_diag, clpmM_diag, cupmM_diag | BUGS: none
test_that("C07: PM.matrix entries are the pairwise co-moments, with crossed divergent triangles", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("clpmM / cupmM / dupmM / dlpmM (definitions): exact dyadic data, entries identical to the pairwise functions", {
    got <- lapply(fc$degree_pairs, function(d) aud_pm_quadrants(d[1], d[2], fc$t_exact, fc$X_exact, FALSE))
    ref <- lapply(fc$degree_pairs, function(d) aud_pm_reference(d[1], d[2], fc$t_exact, fc$X_exact, FALSE))
    expect_identical(got, ref)
  })
  aud_case("clpmM / cupmM / dupmM / dlpmM (definitions): random and real data, mixed degrees, pop_adj scale", {
    grid <- expand.grid(sample = names(fc$pm_samples), deg = seq_along(fc$degree_pairs),
                        pop_adj = c(FALSE, TRUE), stringsAsFactors = FALSE)
    errs <- unlist(Map(function(s, k, pa) {
      X <- fc$pm_samples[[s]]; tg <- colMeans(X); d <- fc$degree_pairs[[k]]
      got <- aud_pm_quadrants(d[1], d[2], tg, X, pa)
      ref <- aud_pm_reference(d[1], d[2], tg, X, pa)
      vapply(names(ref), function(q) max(abs(got[[q]] - ref[[q]])) / aud_mat_tol(ref[[q]], nrow(X)),
             numeric(1))
    }, grid$sample, grid$deg, grid$pop_adj))
    aud_within_cells(errs, 0, 1)
  })
  aud_case("PMmatrix_dupm_transpose: dupm is identical to t(dlpm), normalised or not", {
    for (X in fc$pm_all) for (d in fc$degree_pairs) for (nm in c(FALSE, TRUE)) {
      P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm)
      expect_identical(P$dupm, t(P$dlpm))
    }
  })
  aud_case("PMmatrix_clpm_symm / PMmatrix_cupm_symm / PMmatrix_cov_symm: exactly symmetric", {
    for (X in fc$pm_all) for (d in fc$degree_pairs) for (nm in c(FALSE, TRUE)) {
      P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm)
      expect_identical(P[c("clpm", "cupm", "cov.matrix")],
                       lapply(P[c("clpm", "cupm", "cov.matrix")], t))
    }
  })
  aud_case("dupmM_diag / dlpmM_diag: divergent diagonals are exactly 0", {
    diags <- unlist(lapply(fc$pm_all, function(X) lapply(fc$degree_pairs, function(d) {
      P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE)
      c(diag(P$dupm), diag(P$dlpm))
    })))
    expect_true(all(diags == 0))
  })
  aud_case("clpmM_diag / cupmM_diag: diagonals are s * LPM(2 dL) and s * UPM(2 dU)", {
    for (X in fc$pm_samples) for (d in fc$degree_pairs) for (pa in c(FALSE, TRUE)) {
      tg <- colMeans(X); n <- nrow(X); s <- aud_pm_scale(pa, n, d[1], d[2])
      P <- PM.matrix(d[1], d[2], tg, X, pa)
      lpm2 <- s * vapply(seq_len(ncol(X)), function(i) LPM(2 * d[1], tg[i], X[, i]), numeric(1))
      upm2 <- s * vapply(seq_len(ncol(X)), function(i) UPM(2 * d[2], tg[i], X[, i]), numeric(1))
      aud_within_cells(c(unname(diag(P$clpm)), unname(diag(P$cupm))), c(lpm2, upm2),
                       rep(c(aud_mat_tol(P$clpm, n), aud_mat_tol(P$cupm, n)), each = ncol(X)))
    }
  })
})

# Lean: PMmatrix_cov_eq, covM_one_apply, PMmatrix_cov_eq_popCov, PMmatrix_cov_eq_sampleCov, pmScale_true, PMmatrix_cov_popAdj, PMmatrix_popAdj_ignored | BUGS: #22
test_that("C08: PM.matrix cov.matrix is the population or sample covariance decomposition; pop_adj acts only at positive degrees", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("PMmatrix_cov_eq: cov.matrix = cupm + clpm - dupm - dlpm, normalised or not", {
    for (X in fc$pm_all) for (d in fc$degree_pairs) for (nm in c(FALSE, TRUE)) {
      P <- PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm)
      aud_within_cells(P$cov.matrix, P$cupm + P$clpm - P$dupm - P$dlpm,
                       8 * aud_eps * max(abs(P$cupm) + abs(P$clpm) + abs(P$dupm) + abs(P$dlpm)))
    }
  })
  aud_case("covM_one_apply: exact dyadic data, degree-1 cov.matrix is identical to the co-moment about the targets", {
    P <- unname(PM.matrix(1, 1, fc$t_exact, fc$X_exact, FALSE)$cov.matrix)
    expect_identical(P, aud_pairwise(function(a, b, u, v) mean((a - u) * (b - v)), fc$X_exact, fc$t_exact))
  })
  aud_case("covM_one_apply: random and real data, targets at the 30% quantiles", {
    for (X in fc$pm_samples) {
      n <- nrow(X)
      tg <- apply(X, 2, quantile, probs = 0.3, names = FALSE)
      P <- unname(PM.matrix(1, 1, tg, X, FALSE)$cov.matrix)
      co <- aud_pairwise(function(a, b, u, v) mean((a - u) * (b - v)), X, tg)
      tol <- 4 * (n + 4) * aud_eps * aud_pairwise(function(a, b, u, v) mean(abs((a - u) * (b - v))), X, tg)
      aud_within_cells(P, co, tol)
    }
  })
  aud_case("PMmatrix_cov_eq_popCov: PM.matrix(1, 1, NULL, X, FALSE)$cov.matrix is the population covariance", {
    for (X in fc$pm_samples) {
      n <- nrow(X)
      aud_within_cells(PM.matrix(1, 1, NULL, X, FALSE)$cov.matrix, cov(X) * (n - 1) / n, aud_cov_tol(X))
    }
  })
  aud_case("PMmatrix_cov_eq_sampleCov: PM.matrix(1, 1, NULL, X, TRUE)$cov.matrix is cov(X)", {
    for (X in fc$pm_samples) {
      P <- PM.matrix(1, 1, NULL, X, TRUE)$cov.matrix
      aud_within_cells(P, cov(X), aud_cov_tol(X))
      expect_identical(dimnames(P), dimnames(cov(X)))
    }
  })
  aud_case("pmScale_true / PMmatrix_cov_popAdj: pop_adj scales all quadrants by n / (n - 1) when both degrees > 0", {
    for (X in fc$pm_samples) for (d in list(c(1, 1), c(2, 1), c(1, 3))) {
      tg <- colMeans(X); n <- nrow(X)
      P1 <- PM.matrix(d[1], d[2], tg, X, TRUE)
      P0 <- PM.matrix(d[1], d[2], tg, X, FALSE)
      aud_within_cells(unlist(P1), unlist(P0) * n / (n - 1), 4 * aud_eps * max(abs(unlist(P1))))
    }
  })
  aud_case("PMmatrix_popAdj_ignored: pop_adj changes nothing when either degree is 0", {
    for (X in fc$pm_all) for (d in list(c(0, 0), c(1, 0), c(0, 1), c(2, 0), c(0, 2))) {
      for (nm in c(FALSE, TRUE)) {
        expect_identical(PM.matrix(d[1], d[2], colMeans(X), X, TRUE, nm),
                         PM.matrix(d[1], d[2], colMeans(X), X, FALSE, nm))
      }
    }
  })
  aud_case("PMmatrix_popAdj_ignored: pop_adj changes nothing when n = 1", {
    X1 <- matrix(c(3, 5, -1), 1, dimnames = list(NULL, c("a", "b", "c")))
    for (d in fc$degree_pairs) for (nm in c(FALSE, TRUE)) {
      expect_identical(PM.matrix(d[1], d[2], c(1, 2, 0), X1, TRUE, nm),
                       PM.matrix(d[1], d[2], c(1, 2, 0), X1, FALSE, nm))
    }
  })
})

# Lean: clpmM_posSemidef, cupmM_posSemidef, covM_posSemidef, clpmM_quadForm, cupmM_quadForm, covM_quadForm, PMmatrix_norm_cov_not_posSemidef, PMmatrix_norm_cov_not_posSemidef_mean, PMmatrix_norm_clpm_not_posSemidef | BUGS: none
test_that("C09: unnormalised PM.matrix quadrants are PSD sums of squares; normalised ones are not", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("clpmM_posSemidef / cupmM_posSemidef / covM_posSemidef: unnormalised matrices are PSD", {
    grid <- expand.grid(sample = names(fc$pm_samples), target = 1:3, deg = seq_along(fc$degree_pairs),
                        pop_adj = c(FALSE, TRUE), stringsAsFactors = FALSE)
    margins <- unlist(Map(function(s, k, g, pa) {
      X <- fc$pm_samples[[s]]; d <- fc$degree_pairs[[g]]
      tg <- list(colMeans(X), apply(X, 2, median), X[2, ])[[k]]
      P <- PM.matrix(d[1], d[2], tg, X, pa)
      vapply(P[c("clpm", "cupm", "cov.matrix")], function(M) {
        aud_min_eigen(M) / (8 * (nrow(X) + ncol(M)) * aud_eps * max(abs(M)))
      }, numeric(1))
    }, grid$sample, grid$target, grid$deg, grid$pop_adj))
    expect_true(all(margins >= -1), info = paste(names(margins)[margins < -1], collapse = ", "))
  })
  aud_case("clpmM_quadForm / cupmM_quadForm / covM_quadForm: exact dyadic quadratic forms are sums of squares", {
    v <- c(1, -2, 0.5)
    for (d in list(c(1, 1), c(2, 1), c(0, 2))) {
      P <- PM.matrix(d[1], d[2], fc$t_exact, fc$X_exact, FALSE)
      Lm <- sweep(fc$X_exact, 2, fc$t_exact, function(x, t) ifelse(x <= t, (t - x)^d[1], 0))
      Um <- sweep(fc$X_exact, 2, fc$t_exact, function(x, t) ifelse(x > t, (x - t)^d[2], 0))
      expect_equal(
        c(v %*% P$clpm %*% v, v %*% P$cupm %*% v, v %*% P$cov.matrix %*% v),
        c(mean((Lm %*% v)^2), mean((Um %*% v)^2), mean(((Um - Lm) %*% v)^2)),
        tolerance = 1e-14)
    }
  })
  aud_case("PMmatrix_norm_cov_not_posSemidef: the normalised cov.matrix is not PSD (target 0)", {
    X <- rbind(c(1, 0, 3), c(3, 0, 1), c(2, 0, 0), c(0, 3, 3))
    M <- PM.matrix(1, 1, rep(0, 3), X, FALSE, TRUE)$cov.matrix
    # Cell (1, 2) has total 0 (every row ties in column 1 or 2) and takes the
    # degree-0 value -1; the Lean counterexample uses the earlier 0. Not PSD
    # either way.
    aud_within_cells(M, matrix(c(1, -1, 1, -1, 1, 1, 1, 1, 1), 3), 4 * aud_eps)
    v <- c(1, 1, -1)
    expect_lt(drop(v %*% M %*% v), 0)
  })
  aud_case("PMmatrix_norm_cov_not_posSemidef_mean: not PSD with column-mean targets and positive totals", {
    Z <- rbind(c(-1, -1, 2), c(2, 0, 2), c(1, 2, -2), c(-2, 1, -1))
    M <- PM.matrix(1, 1, NULL, Z, FALSE, TRUE)$cov.matrix
    aud_within_cells(M, matrix(c(1, .2, .2, .2, 1, -1, .2, -1, 1), 3), 1e-14)
    v <- c(-2, 5, 5)
    aud_within_cells(drop(v %*% M %*% v), -4, 1e-12)
    expect_lt(aud_min_eigen(M), -0.07)
  })
  aud_case("PMmatrix_norm_clpm_not_posSemidef: the normalised clpm is not PSD", {
    W <- rbind(c(-2, -2, 0), c(-2, 1, -2), c(-2, 0, 1))
    M <- PM.matrix(1, 1, rep(0, 3), W, FALSE, TRUE)$clpm
    aud_within_cells(M, matrix(c(1, 2 / 3, 2 / 3, 2 / 3, .8, 0, 2 / 3, 0, .8), 3), 1e-14)
    v <- c(-4, 3, 3)
    aud_within_cells(drop(v %*% M %*% v), -8 / 5, 1e-12)
  })
})

# Lean: PMmatrix_norm_sum_eq_one, normQ_mem_Icc, PMmatrix_norm_cov_mem_Icc, PMmatrix_norm_popAdj, PMmatrix_zero_zero_sum, PMmatrix_zero_zero_norm | BUGS: none
test_that("C10: PM.matrix(norm = TRUE) cells sum to 1, lie in [0, 1] and cancel pop_adj", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("PMmatrix_norm_sum_eq_one: normalised quadrants sum to 1 in every cell with positive total", {
    for (X in fc$pm_samples) for (d in fc$degree_pairs) {
      tg <- colMeans(X)
      expect_true(all(aud_cell_sum(PM.matrix(d[1], d[2], tg, X, FALSE)) > 0))
      aud_within_cells(aud_cell_sum(PM.matrix(d[1], d[2], tg, X, FALSE, TRUE)), 1, 8 * aud_eps)
    }
  })
  aud_case("PMmatrix_norm_cov_mem_Icc / normQ_mem_Icc: normalised quadrants lie in [0, 1], cov in [-1, 1]", {
    for (X in fc$pm_all) for (d in fc$degree_pairs) for (pa in c(FALSE, TRUE)) {
      P <- PM.matrix(d[1], d[2], colMeans(X), X, pa, TRUE)
      q <- unlist(P[c("cupm", "dupm", "dlpm", "clpm")])
      expect_true(all(q >= 0 & q <= 1 + 2 * aud_eps) && all(abs(P$cov.matrix) <= 1 + 4 * aud_eps))
    }
  })
  aud_case("PMmatrix_norm_popAdj: pop_adj cancels under norm = TRUE", {
    for (X in fc$pm_samples) for (d in fc$degree_pairs) {
      tg <- colMeans(X)
      aud_within_cells(unlist(PM.matrix(d[1], d[2], tg, X, TRUE, TRUE)),
                       unlist(PM.matrix(d[1], d[2], tg, X, FALSE, TRUE)), 8 * aud_eps)
    }
  })
  aud_case("PMmatrix_zero_zero_sum: exact dyadic data, degree-0 quadrants sum to exactly 1, normalised or not", {
    for (nm in c(FALSE, TRUE)) {
      expect_identical(unname(aud_cell_sum(PM.matrix(0, 0, fc$t_exact, fc$X_exact, FALSE, nm))), matrix(1, 3, 3))
    }
  })
  aud_case("PMmatrix_zero_zero_sum: random and real data, degree-0 quadrants sum to 1 within 8 eps", {
    for (X in fc$pm_samples) for (nm in c(FALSE, TRUE)) {
      aud_within_cells(aud_cell_sum(PM.matrix(0, 0, colMeans(X), X, TRUE, nm)), 1, 8 * aud_eps)
    }
  })
  aud_case("PMmatrix_zero_zero_norm: at degree 0 norm = TRUE changes nothing", {
    expect_identical(PM.matrix(0, 0, fc$t_exact, fc$X_exact, FALSE, TRUE),
                     PM.matrix(0, 0, fc$t_exact, fc$X_exact, FALSE, FALSE))
    for (X in fc$pm_samples) {
      aud_within_cells(unlist(PM.matrix(0, 0, colMeans(X), X, TRUE, TRUE)),
                       unlist(PM.matrix(0, 0, colMeans(X), X, TRUE, FALSE)), 8 * aud_eps)
    }
  })
})

# Lean: PMmatrix_norm_of_total_eq_zero, PMmatrix_norm_sum_of_tied, PMmatrix_norm_sum_eq_zero_example, PMmatrix_norm_cov_diag, PMmatrix_norm_cov_diag_zero_left, pmTotal_zero_left_pos, PMmatrix_norm_sum_eq_one_zero_left, pmTotal_of_tied_zero_left | BUGS: #19, e89cc99
test_that("C11: BUGS #19: PM.matrix(norm = TRUE) zero-total cells take the degree-0 values", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("PM.matrix(norm = TRUE) gives a column constant at its target the degree-0 cells", {
    M <- cbind(a = c(1, 2, 3), b = c(5, 5, 5))
    P <- PM.matrix(1, 1, NULL, M, FALSE, TRUE)
    # LPM.ratio convention: a constant variable is all lower.
    expect_identical(P$clpm["b", "b"], 1)
    expect_identical(c(P$cupm["b", "b"], P$dupm["b", "b"], P$dlpm["b", "b"]), c(0, 0, 0))
  })
  aud_case("PM.matrix(target = NULL) ties a constant column to its value, not to its rounded mean", {
    # The compiled mean of three 0.1 values is 0.1 plus one ulp, so the
    # deviations used to be about -1.4e-17 instead of ties, and the zero-total
    # convention depended on how the mean rounded (b = 5 took it, b = 0.1 not).
    for (v in c(0.1, 5, -1 / 3, 1e300)) {
      M <- cbind(a = c(1, 2, 3), b = v)
      expect_identical(PM.matrix(1, 1, NULL, M, FALSE, TRUE), PM.matrix(1, 1, c(2, v), M, FALSE, TRUE))
      expect_identical(PM.matrix(1, 1, NULL, M, FALSE, TRUE)$cov.matrix[1, 2], 1 / 3)
      P <- PM.matrix(1, 1, NULL, M, FALSE, FALSE)
      expect_identical(c(P$clpm["b", "b"], P$cupm["b", "b"], P$clpm["a", "b"]), c(0, 0, 0))
      # Degree 0: every row of b is a tie, counted as lower.
      expect_identical(PM.matrix(0, 0, NULL, M, FALSE, FALSE)$clpm["b", "b"], 1)
    }
    # A non-constant column keeps the compiled mean.
    M <- cbind(a = c(0.1, 0.2, 0.4), b = 0.1)
    expect_identical(PM.matrix(1, 1, NULL, M, FALSE, FALSE),
                     PM.matrix(1, 1, c((0.1 + 0.2 + 0.4) / 3, 0.1), M, FALSE, FALSE))
  })
  aud_case("PM.matrix(norm = TRUE) uses the degree-0 cells when every row ties in one of the two columns", {
    # Neither column is constant, but each row ties the target in a or b.
    M <- cbind(a = c(0, 0, 1, -2), b = c(3, -1, 0, 0))
    P <- PM.matrix(2, 1, c(0, 0), M, TRUE, TRUE)
    P0 <- PM.matrix(0, 0, c(0, 0), M, FALSE, FALSE)
    expect_equal(c(P$cupm["a", "b"], P$dupm["a", "b"], P$dlpm["a", "b"], P$clpm["a", "b"]),
                 c(0, 0.25, 0.25, 0.5))
    for (q in c(aud_quadrant_names, "cov.matrix")) expect_equal(P[[q]]["a", "b"], P0[[q]]["a", "b"])
    # The crossed mirror holds for the fallback cell too.
    expect_identical(P$dupm["b", "a"], P$dlpm["a", "b"])
    expect_identical(P$dlpm["b", "a"], P$dupm["a", "b"])
    expect_equal(aud_cell_sum(P), matrix(1, 2, 2, dimnames = dimnames(P$cupm)))
  })
  aud_case("PMmatrix_norm_of_total_eq_zero (degree-0 convention): a tied column takes the degree-0 cells (dL > 0)", {
    P0 <- PM.matrix(0, 0, fc$t_tied, fc$X_tied, FALSE, FALSE)
    for (d in list(c(1, 1), c(2, 1), c(1, 0))) {
      P <- PM.matrix(d[1], d[2], fc$t_tied, fc$X_tied, FALSE, TRUE)
      for (q in names(P)) {
        aud_within_cells(P[[q]]["b", ], P0[[q]]["b", ], 4 * aud_eps)
        aud_within_cells(P[[q]][, "b"], P0[[q]][, "b"], 4 * aud_eps)
      }
      aud_within_cells(aud_cell_sum(P), 1, 8 * aud_eps)
    }
  })
  aud_case("PMmatrix_norm_cov_diag: dL > 0 gives diagonal 1, also for a tied column", {
    P <- PM.matrix(1, 1, NULL, cbind(a = c(1, 2, 3), b = c(5, 5, 5)), FALSE, TRUE)
    aud_within_cells(unname(P$cov.matrix), matrix(c(1, 1 / 3, 1 / 3, 1), 2), 4 * aud_eps)
    P <- PM.matrix(2, 1, fc$t_tied, fc$X_tied, TRUE, TRUE)
    aud_within_cells(unname(diag(P$cov.matrix)), c(1, 1, 1), 4 * aud_eps)
  })
  aud_case("PMmatrix_norm_sum_eq_zero_example (degree-0 convention): PM.matrix(1, 1, 0, matrix(0), TRUE, TRUE) is all lower", {
    P <- PM.matrix(1, 1, 0, matrix(0), TRUE, TRUE)
    expect_identical(vapply(P, function(M) M[1, 1], numeric(1)),
                     c(cupm = 0, dupm = 0, dlpm = 0, clpm = 1, cov.matrix = 1))
  })
  aud_case("PMmatrix_norm_cov_diag_zero_left: at dL = 0 the normalised cov diagonal is 1, even for a tied column", {
    P <- PM.matrix(0, 1, NULL, cbind(a = 1:3, b = c(5, 5, 5)), FALSE, TRUE)
    aud_within_cells(unname(P$cov.matrix), matrix(c(1, 1 / 3, 1 / 3, 1), 2), 4 * aud_eps)
    for (dU in 0:2) {
      aud_within_cells(unname(diag(PM.matrix(0, dU, fc$t_tied, fc$X_tied, TRUE, TRUE)$cov.matrix)), 1, 4 * aud_eps)
    }
  })
  aud_case("PMmatrix_norm_sum_eq_one_zero_left / pmTotal_zero_left_pos: at dL = 0 every cell sums to 1", {
    for (dU in 0:2) {
      aud_within_cells(aud_cell_sum(PM.matrix(0, dU, fc$t_tied, fc$X_tied, TRUE, TRUE)), 1, 8 * aud_eps)
    }
  })
  aud_case("pmTotal_of_tied_zero_left: a tied column has diagonal total exactly 1 at dL = 0", {
    totals <- vapply(0:2, function(dU) aud_cell_sum(PM.matrix(0, dU, fc$t_tied, fc$X_tied, TRUE, FALSE))["b", "b"],
                     numeric(1))
    expect_identical(totals, c(1, 1, 1))
  })
})

# Lean: none | BUGS: #12, #19, 59e155f, 82ffbe2
test_that("C12: PM.matrix(norm = TRUE) is scale-free when cell totals underflow or overflow", {
  withr::local_seed(20261108)
  aud_case("PM.matrix(norm = TRUE) is scale-free at equal degrees when the products underflow or overflow", {
    withr::local_seed(2)
    X <- cbind(rnorm(30), rexp(30), runif(30))
    for (d in 1:2) {
      ref <- unlist(PM.matrix(d, d, NULL, X, FALSE, TRUE))
      for (s in c(1e-200, 2^-700, 1e200)) {
        P <- PM.matrix(d, d, NULL, X * s, FALSE, TRUE)
        aud_close_flat(unlist(P), ref, 16 * .Machine$double.eps)
      }
    }
    # Old behaviour: every product underflowed, all cells 0.
    P <- PM.matrix(1, 1, NULL, X * 1e-200, FALSE, TRUE)
    aud_close_flat(diag(P$cov.matrix), rep(1, 3), 4 * .Machine$double.eps)
  })
  aud_case("PM.matrix(norm = TRUE) at degree 200 does not return 0 when the products overflow", {
    withr::local_seed(2)
    X <- cbind(rnorm(30), rexp(30), runif(30))
    P <- PM.matrix(200, 200, NULL, X * 1e3, FALSE, TRUE)
    expect_true(all(is.finite(unlist(P))))
    aud_close_flat(aud_cell_sum(P), 1, 8 * .Machine$double.eps)
    aud_close_flat(unlist(P), unlist(PM.matrix(200, 200, NULL, X, FALSE, TRUE)), 1e-14)
  })
  aud_case("PM.matrix(norm = TRUE) with unequal degrees matches the exact ratio when the total overflows", {
    # Exact values from rational arithmetic on the double inputs (nns-proofs
    # oracle pm_matrix), rounded to double, in unlist() order: cupm, dupm, dlpm,
    # clpm, each column-major. A 0 is an exact 0 (ties) or a value below the
    # smallest denormal.
    pm_check_exact <- function(X, exact, tol) {
      for (key in names(exact)) {
        d <- as.numeric(strsplit(key, " ")[[1]])
        got <- unname(unlist(PM.matrix(d[1], d[2], c(0, 0), X, TRUE, TRUE)[aud_quadrant_names]))
        e <- exact[[key]]
        expect_identical(got == 0, e == 0)
        expect_lte(max(abs(got - e)[e != 0] / e[e != 0]), tol(d))
      }
    }
    # The rescue's relative error is of order degree * eps: the quadrant exponent
    # a log2 s_i + b log2 s_j is rounded to about |exponent| * eps before exp2,
    # and the rounding of |d| / s is raised to the degree.
    bound <- function(d) 64 * max(d) * .Machine$double.eps

    # Dyadic data: s = 2 and every |d| / s is exact, so the exponents are
    # integers and the only rounding is in the sums. This exercises the
    # exact-exponent path only, hence the tighter tolerance. The plain sums
    # overflow (2^2002), so every cell goes through the rescue.
    X <- cbind(a = c(-1.75, -0.5, 0.25, 1.5, 2, -1.25, 0.75, 1.25),
               b = c(1.5, -2, 0.5, -0.75, 1.25, 0.25, -1.5, 1.75))
    pm_check_exact(X, list(
      "1000 1001" = c(1, 1.0778285412673956e-21, 1.0778285412673956e-21, 3.1782067229759299e-116,
                      0, 0, 1, 0, 0, 1, 0, 0, 2.5944544677354528e-117, 0, 0, 1),
      "700 650" = c(1, 1.1886477434465246e-26, 1.1886477434465246e-26, 3.2171355712192349e-106,
                    0, 2.9359596119441445e-243, 1, 0, 0, 1, 2.9359596119441445e-243, 0,
                    8.2086889175245743e-52, 2.594463309812418e-285, 2.594463309812418e-285, 1)),
      function(d) 1e-14)

    # Non-dyadic scales: the same data times 1.3 and 0.7 (relative error about
    # 2.7e-14 at degrees (1000, 1001)).
    Y <- cbind(a = X[, "a"] * 1.3, b = X[, "b"] * 0.7)
    pm_check_exact(Y, list(
      "1000 1001" = c(1, 1.4011771036479542e-21, 1.4011771036479542e-21, 1.5573212942580642e-116,
                      0, 0, 1, 0, 0, 1, 0, 0, 1.53518015842312e-117, 0, 0, 1),
      "700 650" = c(1, 2.3871821306951924e-32, 2.3871821306951924e-32, 9.946387187063837e-91,
                    0, 1.0604359987379094e-256, 1, 0, 0, 1, 1.0604359987379094e-256, 0,
                    2.0352090648213927e-40, 4.666051567129481e-293, 4.666051567129481e-293, 1)),
      bound)
    Z <- cbind(a = c(2, 0, 1.9, -1.7, 1.1) * 1.3, b = c(0, 2, 1.1, -1.9, 1.7) * 0.7)
    pm_check_exact(Z, list(
      "2000 2001" = c(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      7.010615702072899e-284, 1, 1, 4.000968056720675e-90)),
      bound)
    # Cell (a, b) has two comparable terms whose rescaled products are denormal
    # (about 3e-313), so it is summed in the log domain (relative error about
    # 2.9e-13 in clpm).
    W <- cbind(a = c(4, 0, 3, -2.95) * 1.3, b = c(0, 4, 3, -3.05) * 0.7)
    pm_check_exact(W, list(
      "1250 1251" = c(1, 0.9205749312463112, 0.9205749312463112, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0.07942506875368878, 0.07942506875368878, 5.073457005883048e-296)),
      bound)
  })
  aud_case("PM.matrix(norm = TRUE) log-domain rescue keeps a term whose ratio to the column scale is denormal", {
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
  aud_case("PM.matrix(norm = TRUE) results with a safe total are unchanged", {
    # Bit pattern of the plain division path (no rescue involved).
    withr::local_seed(123)
    A <- cbind(rnorm(100), rnorm(100), rnorm(100))
    for (d in list(c(1, 1), c(2, 1), c(0, 2))) {
      raw <- PM.matrix(d[1], d[2], NULL, A, TRUE, FALSE)
      tot <- raw$cupm + raw$dupm + raw$dlpm + raw$clpm
      P <- PM.matrix(d[1], d[2], NULL, A, TRUE, TRUE)
      # Upper triangle: the lower one is the (crossed) mirror of it.
      ut <- upper.tri(tot, diag = TRUE)
      for (q in aud_quadrant_names) expect_identical(P[[q]][ut], (raw[[q]] / tot)[ut])
    }
  })
})

# Lean: quadrants_sum_eq_one | BUGS: none
test_that("C13: nD co-moments match the k-column product oracle, including n = 1 and constant samples", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("SPEC DPM_nD (nD analogue of quadrants_sum_eq_one): exact dyadic data, degree-0 cells sum to exactly 1", {
    expect_identical(c(aud_nd_sum(fc$X_exact, fc$t_exact, 0, FALSE), aud_nd_sum(fc$X_exact[, 1:2], fc$t_exact[1:2], 0, FALSE)),
                     c(1, 1))
  })
  aud_case("SPEC DPM_nD (nD analogue of quadrants_sum_eq_one): random and real data, degree-0 cells partition the sample", {
    sums <- unlist(lapply(fc$pm_samples, function(X) {
      lapply(list(colMeans(X), apply(X, 2, median), X[3, ]), function(tg) {
        c(aud_nd_sum(X, tg, 0, FALSE), aud_nd_sum(X[, 1:2], tg[1:2], 0, TRUE))
      })
    }))
    aud_within_cells(sums, 1, 4 * aud_eps)
  })
  aud_case("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: k = 3 exact dyadic data, identical to the k-column product oracle", {
    per_target <- function(f) unlist(lapply(fc$nd_exact_targets, function(t) {
      unlist(lapply(0:2, function(d) setNames(f(fc$X_exact, t, d), paste0("d", d, ".", c("clpm", "cupm", "dpm")))))
    }))
    got <- per_target(aud_nd_triple)
    ref <- per_target(aud_nd_oracle)
    expect_identical(got, ref)
    # the mixed target really exercises all three cells at d > 0
    mixed_pos <- got[grep("^mixed\\.d[12]\\.", names(got))]
    expect_true(length(mixed_pos) == 6 && all(mixed_pos > 0))
  })
  aud_case("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: k = 3..4, random and real data, within n * eps of the oracle", {
    vals <- function(f) unlist(lapply(fc$pm_samples, function(X) {
      lapply(list(mean = colMeans(X), obs = X[3, ]), function(tg) {
        unlist(lapply(0:2, function(d) setNames(f(X, unname(tg), d), paste0("d", d, ".", c("clpm", "cupm", "dpm")))))
      })
    }))
    aud_within_cells(vals(aud_nd_triple), vals(aud_nd_oracle), vals(function(X, t, d) rep(aud_nd_tol(X, t, d), 3)))
  })
  aud_case("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: n = 1 rows (k = 3)", {
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
  aud_case("SPEC Co.LPM_nD / DPM_nD: a sample constant at its target, raw cells", {
    X <- matrix(1, 3, 2)
    expect_identical(c(aud_nd_triple(X, c(1, 1), 0), aud_nd_triple(X, c(1, 1), 1)),
                     c(clpm = 1, cupm = 0, dpm = 0, clpm = 0, cupm = 0, dpm = 0))
  })
  aud_case("degree-0 nD co-moments partition a sample tied at the target", {
    X <- cbind(c(1, 2, 3, 3, 4, 5), c(3, 3, 1, 5, 3, 2))
    t <- c(3, 3)
    expect_equal(c(Co.LPM_nD(X, t), Co.UPM_nD(X, t), DPM_nD(X, t)), c(0.5, 0, 0.5))
    G <- as.matrix(expand.grid(-1:1, -1:1, -1:1))
    t <- c(0, 0, 0)
    expect_equal(c(Co.LPM_nD(G, t), Co.UPM_nD(G, t), DPM_nD(G, t)), c(8, 1, 18) / 27)
  })
})

# Lean: none | BUGS: none
test_that("C14: nD co-moments reduce to the bivariate and univariate functions at k = 2 and k = 1", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("SPEC Co.LPM_nD: k = 2 is identical to Co.LPM", {
    expect_identical(aud_nd_cases(function(c, d) Co.LPM_nD(cbind(c$x, c$y), c(c$tx, c$ty), d, FALSE)),
                     aud_nd_cases(function(c, d) aud_co_lpm(c, d, d)))
  })
  aud_case("SPEC Co.UPM_nD: k = 2 reduces to Co.UPM", {
    aud_within_cells(aud_nd_cases(function(c, d) Co.UPM_nD(cbind(c$x, c$y), c(c$tx, c$ty), d, FALSE)),
                     aud_nd_cases(function(c, d) aud_co_upm(c, d, d)),
                     aud_nd_cases(function(c, d) aud_co_tol(c, d, d, k = 2)))
  })
  aud_case("SPEC DPM_nD: k = 2 reduces to D.LPM(d, d) + D.UPM(d, d)", {
    aud_within_cells(aud_nd_cases(function(c, d) DPM_nD(cbind(c$x, c$y), c(c$tx, c$ty), d, FALSE)),
                     aud_nd_cases(function(c, d) aud_d_lpm(c, d, d) + aud_d_upm(c, d, d)),
                     aud_nd_cases(function(c, d) aud_co_tol(c, d, d)))
  })
  aud_case("SPEC Co.LPM_nD / Co.UPM_nD / DPM_nD: one column reduces to LPM / UPM / 0", {
    one_col <- function(c, d) {
      c(Co.LPM_nD(matrix(c$x), c$tx, d, FALSE), Co.UPM_nD(matrix(c$x), c$tx, d, FALSE),
        DPM_nD(matrix(c$x), c$tx, d, FALSE))
    }
    univ <- function(c, d) c(LPM(d, c$tx, c$x), UPM(d, c$tx, c$x), 0)
    tol <- function(c, d) rep(aud_co_tol(list(x = c$x, y = c$x, tx = c$tx, ty = c$tx), d, 0), 3)
    over3 <- function(f) aud_over(fc$all_cases, function(c) unlist(lapply(0:2, function(d) f(c, d))))
    aud_within_cells(over3(one_col), over3(univ), over3(tol))
  })
})

# Lean: none | BUGS: #12, #18
test_that("C15: BUGS #18: normalised nD co-moments take the degree-0 values on zero totals and survive underflow", {
  withr::local_seed(20261108)
  fc <- aud_fx_co()
  aud_case("normalised nD partial moments with every row tied take the degree-0 values", {
    X <- matrix(1, 3, 2)
    expect_identical(c(Co.LPM_nD(X, c(1, 1), 1), Co.UPM_nD(X, c(1, 1), 1), DPM_nD(X, c(1, 1), 1)),
                     c(1, 0, 0))
  })
  aud_case("normalised nD partial moments with a tied column sum to 1", {
    X <- cbind(1:5, rep(1, 5))
    cells <- c(Co.LPM_nD(X, c(3, 1), 2), Co.UPM_nD(X, c(3, 1), 2), DPM_nD(X, c(3, 1), 2))
    expect_identical(cells, c(Co.LPM_nD(X, c(3, 1), 0), Co.UPM_nD(X, c(3, 1), 0), DPM_nD(X, c(3, 1), 0)))
    expect_equal(sum(cells), 1)
  })
  aud_case("one-column normalised Co.LPM_nD on constant data agrees with LPM.ratio", {
    expect_identical(Co.LPM_nD(matrix(1, 3, 1), 1, 1), LPM.ratio(1, 1, rep(1, 3)))
    expect_identical(Co.UPM_nD(matrix(1, 3, 1), 1, 1), UPM.ratio(1, 1, rep(1, 3)))
  })
  aud_case("the batched nD kernel takes the degree-0 value on a zero total", {
    withr::local_seed(1)
    Z <- cbind(rnorm(6), 0)
    expect_identical(NNS:::Co.LPM_nD.batch(Z, Z, 1, TRUE), NNS:::Co.LPM_nD.batch(Z, Z, 0, TRUE))
  })
  aud_case("the batched upper kernel is the reflected lower kernel on ordinary data", {
    withr::local_seed(3)
    Y <- matrix(rnorm(60), 20)
    expect_identical(NNS:::Co.LPM_nD.batch(Y, Y, 2, TRUE, upper = TRUE),
                     NNS:::Co.LPM_nD.batch(-Y, -Y, 2, TRUE))
    expect_equal(NNS:::Co.LPM_nD.batch(Y, Y, 2, TRUE, upper = TRUE),
                 apply(Y, 1, function(t) Co.UPM_nD(Y, t, 2)))
  })
  aud_case("normalised Co.LPM_nD is not 0 when every product underflows", {
    y <- matrix(c(-1, 1, -2, 2), 2)
    expect_equal(Co.LPM_nD(1e-100 * y, c(0, 0), 4), Co.LPM_nD(y, c(0, 0), 4))
    expect_equal(Co.LPM_nD(1e-100 * y, c(0, 0), 4), 0.5)
  })
  aud_case("normalised nD partial moments are invariant to per-column scale at degree 200", {
    withr::local_seed(1)
    M <- cbind(rnorm(250, 0, 0.01), rnorm(250, 0, 0.01))
    cells <- function(D) c(Co.LPM_nD(D, c(0, 0), 200), Co.UPM_nD(D, c(0, 0), 200), DPM_nD(D, c(0, 0), 200))
    expect_equal(cells(M %*% diag(c(1e3, 1e-3))), cells(M), tolerance = 1e-10)
    expect_equal(sum(cells(M)), 1)
  })
  aud_case("an underflowed nD total with one tied row uses the untied rows, not the degree-0 counts", {
    Z <- rbind(c(1, 0), c(-1e-200, -1e-200), c(2e-200, 1e-200))
    expect_equal(c(Co.LPM_nD(Z, c(0, 0), 2), Co.UPM_nD(Z, c(0, 0), 2), DPM_nD(Z, c(0, 0), 2)),
                 c(0.2, 0.8, 0))
  })
  aud_case("SPEC Co.LPM_nD: norm = TRUE cells sum to 1 when the total is positive (degree > 0)", {
    sums <- unlist(lapply(fc$pm_samples, function(X) vapply(1:2, function(d) aud_nd_sum(X, colMeans(X), d, TRUE),
                                                            numeric(1))))
    aud_within_cells(sums, 1, 8 * aud_eps)
  })
})

# Lean: none | BUGS: #11
test_that("C16: BUGS #11: multivariate NNS.CDF with a constant column", {
  withr::local_seed(20261108)
  aud_case("multivariate NNS.CDF at degree 1 with a constant column equals the degree-0 CDF", {
    X <- cbind(1:5, rep(1, 5))
    expect_identical(NNS.CDF(X, 1, plot = FALSE)$Function$CDF, NNS.CDF(X, 0, plot = FALSE)$Function$CDF)
  })
  aud_case("multivariate NNS.CDF target.value agrees with the curve on a constant column", {
    X <- cbind(1:5, rep(1, 5))
    expect_identical(NNS.CDF(X, 1, target = c(3, 1), plot = FALSE)$target.value,
                     NNS.CDF(X, 1, plot = FALSE)$Function$CDF[3])
  })
  aud_case("multivariate NNS.CDF survival at degree 1 is Co.UPM_nD on a constant column", {
    X <- cbind(1:5, rep(1, 5))
    S <- NNS.CDF(X, 1, type = "survival", plot = FALSE)$Function$CDF
    expect_identical(S, apply(X, 1, function(t) Co.UPM_nD(X, t, 1)))
    expect_identical(S, rep(0, 5))
  })
})

# pin: NNS-Python reference values
test_that("C17: NNS-Python reference values: Co.LPM, Co.UPM, D.LPM and D.UPM for vectors and one-column frames", {
  withr::local_seed(20261108)
  x <- nns_py$x
  y <- nns_py$y
  x_df <- as.data.frame(x)
  y_df <- as.data.frame(y)
  # Co.LPM and Co.UPM take one degree (d2 is NA); D.LPM and D.UPM take two.
  ref <- read.table(header = TRUE, text = "
    fn     d1 d2 target expected
    Co.UPM 0  NA NULL   0.28
    Co.UPM 0  NA mean   0.28
    Co.UPM 1  NA mean   0.01204606
    Co.UPM 2  NA mean   0.0009799173
    Co.LPM 0  NA NULL   0.24
    Co.LPM 0  NA mean   0.24
    Co.LPM 1  NA mean   0.01058035
    Co.LPM 2  NA mean   0.0008940764
    D.LPM  0  0  NULL   0.23
    D.LPM  0  0  mean   0.23
    D.LPM  1  0  mean   0.06404049
    D.LPM  0  1  mean   0.05311669
    D.LPM  1  1  mean   0.01513793
    D.LPM  2  0  mean   0.02248309
    D.LPM  0  2  mean   0.01727327
    D.LPM  2  2  mean   0.001554909
    D.UPM  0  0  NULL   0.25
    D.UPM  0  0  mean   0.25
    D.UPM  0  1  mean   0.05488706
    D.UPM  1  0  mean   0.05843498
    D.UPM  1  1  mean   0.01199175
    D.UPM  0  2  mean   0.01512857
    D.UPM  2  0  mean   0.01926167
    D.UPM  2  2  mean   0.0009941733")
  ref$d1 <- as.numeric(ref$d1)
  ref$d2 <- as.numeric(ref$d2)
  ref <- rbind(cbind(ref, input = "vector"), cbind(ref, input = "frame"))
  ref$label <- sprintf("%s(%s) at the %s target, %s input", ref$fn,
                       ifelse(is.na(ref$d2), ref$d1, paste(ref$d1, ref$d2, sep = ", ")), ref$target, ref$input)
  aud_cases(ref, function(r) {
    vec <- r$input == "vector"
    xx <- if (vec) x else x_df
    yy <- if (vec) y else y_df
    tx <- if (r$target == "NULL") NULL else if (vec) mean(x) else colMeans(x_df)
    ty <- if (r$target == "NULL") NULL else if (vec) mean(y) else colMeans(y_df)
    got <- if (is.na(r$d2)) match.fun(r$fn)(r$d1, xx, yy, tx, ty) else match.fun(r$fn)(r$d1, r$d2, xx, yy, tx, ty)
    expect_equal(got, r$expected, tolerance = 1e-5)
  })
})

# pin: NNS-Python reference values
test_that("C18: PM.matrix reference covariance of a 3 x 2 matrix and frame, mean and NULL target, and the norm = TRUE identity", {
  withr::local_seed(20261108)
  A <- matrix(c(1,1,3,2,2,3), ncol = 2)
  T1 <- matrix(c(1.3333333, 0.6666667, 0.6666667, 0.3333333), ncol=2)
  T2 <- matrix(c(0.8888889, 0.4444444, 0.4444444, 0.2222222), ncol=2)
  aud_case("the reference T1 is cov(A)", {
    expect_equal(T1, cov(A), tolerance=1e-5)
  })
  # pop_adj = TRUE gives T1 and FALSE gives T2; a frame's result carries its V1, V2 names.
  ref <- expand.grid(input = c("matrix", "frame"), target = c("mean", "NULL"), pop_adj = c(TRUE, FALSE),
                     stringsAsFactors = FALSE)
  ref$label <- sprintf("PM.matrix(1, 1) of a %s, %s target, pop_adj = %s", ref$input, ref$target, ref$pop_adj)
  aud_cases(ref, function(r) {
    a <- if (r$input == "matrix") A else as.data.frame(A)
    target <- if (r$target == "mean") colMeans(a) else NULL
    want <- if (r$pop_adj) T1 else T2
    if (r$input == "frame") dimnames(want) <- list(c("V1", "V2"), c("V1", "V2"))
    expect_equal(NNS::PM.matrix(1, 1, target, a, pop_adj = r$pop_adj)$cov.matrix, want, tolerance = 1e-5)
  })
  aud_case("NNS::PM.matrix - norm TRUE returns signed normalized covariance decomposition", {
    x <- nns_py$x
    y <- nns_py$y
    z <- nns_py$z
    A <- cbind(x, y, z)
    pm <- NNS::PM.matrix(1, 1, NULL, A, pop_adj = TRUE, norm = TRUE)

    expect_equal(
      pm$cov.matrix,
      pm$cupm + pm$clpm - pm$dlpm - pm$dupm,
      tolerance = 1e-10
    )
    expect_equal(unname(diag(pm$cov.matrix)), rep(1, ncol(A)), tolerance = 1e-10)
  })
})

test_that("C19: co-moments apply each degree to its own variable when one is fractional", {
  withr::local_seed(20261108)
  lpm_leg <- function(v, t, d) if (d == 0) as.numeric(v <= t) else pmax(t - v, 0)^d
  upm_leg <- function(v, t, d) if (d == 0) as.numeric(v > t) else pmax(v - t, 0)^d
  expect_equal(D.UPM(.5, 1, c(-4, -9), c(2, 3), 0, 0), 6.5)
  expect_equal(D.UPM(1, .5, c(-2, -3), c(4, 9), 0, 0), 6.5)

  set.seed(11)
  x <- rnorm(60); y <- rnorm(60)
  degrees <- c(0, .5, 1, 1.5, 2, 3)
  # One comparison per function and target pair over the 36 degree pairs. As
  # lists, each cell keeps the scalar tolerance; a plain vector comparison
  # averages the error over the cells that differ.
  ab <- expand.grid(a = degrees, b = degrees)
  cells <- function(f) as.list(mapply(f, ab$a, ab$b))
  for (tx in c(0.1, -0.3)) for (ty in c(-0.2, 0.4)) {
    expect_equal(cells(function(a, b) Co.LPM(a, x, y, tx, ty, b)),
                 cells(function(a, b) mean(lpm_leg(x, tx, a) * lpm_leg(y, ty, b))))
    expect_equal(cells(function(a, b) Co.UPM(a, x, y, tx, ty, b)),
                 cells(function(a, b) mean(upm_leg(x, tx, a) * upm_leg(y, ty, b))))
    expect_equal(cells(function(a, b) D.LPM(a, b, x, y, tx, ty)),
                 cells(function(a, b) mean(upm_leg(x, tx, b) * lpm_leg(y, ty, a))))
    expect_equal(cells(function(a, b) D.UPM(a, b, x, y, tx, ty)),
                 cells(function(a, b) mean(lpm_leg(x, tx, a) * upm_leg(y, ty, b))))
  }
})
