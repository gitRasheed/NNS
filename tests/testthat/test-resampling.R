# NNS.meboot, NNS.MC and NNS.SS.

test_that("RES01: NNS.MC rho grid: inside the bounds, towards upper_rho, default full range, invalid grids, residual-scale switch", {
  withr::local_seed(20261108)
  mc_rho_grid <- fx_mc_rho_grid()
  aud_case("NNS.MC exp keeps the rho grid inside the bounds", env = mc_rho_grid, {
    for (b in list(c(0.2, 0.8), c(-0.8, -0.2), c(-0.5, 0.5))) {
      r <- mc_rhos(lower_rho = b[1], upper_rho = b[2], by = 0.2, exp = 2)
      expect_equal(range(r), b, info = paste(b, collapse = ", "))
      expect_true(all(r >= b[1] & r <= b[2]))
    }
  })
  aud_case("NNS.MC exp moves the grid towards upper_rho", env = mc_rho_grid, {
    r <- mc_rhos(lower_rho = 0.2, upper_rho = 0.8, by = 0.2, exp = 2)
    expect_true(all(r[2:3] > c(0.6, 0.4)))
    r <- mc_rhos(lower_rho = -0.8, upper_rho = -0.2, by = 0.2, exp = 2)
    expect_true(all(r[2:3] > c(-0.4, -0.6)))
  })
  aud_case("NNS.MC keeps its default full-range grid", {
    r <- seq(-1, 1, 0.01)
    for (e in c(1, 2)) {
      old <- rev(c(-(abs(r[r <= 0])^e), r[r > 0]^(1 / e)))
      set.seed(3)
      m <- NNS.MC(rnorm(20), reps = 1, drift = FALSE, exp = e)
      expect_identical(names(m$replicates), paste0("rho = ", old))
    }
  })
  aud_case("Monte Carlo exponent and interval validation rejects invalid grids", {
    for (e in list(0, -1, Inf, NA_real_, c(1, 2)))
      expect_error(NNS.MC(1:10, exp = e), "exp must")
    expect_error(NNS.MC(1:10, lower_rho = .8, upper_rho = .2), "lower_rho")
    expect_error(NNS.MC(1:10, by = 0), "by must")
  })
  aud_case("residual-scale expansion can be disabled independently", {
    set.seed(5)
    x <- rnorm(60)
    run <- function(expand) {
      set.seed(8)
      NNS.meboot(x, reps = 8, rho = 0, drift = FALSE,
                 expand.sd = FALSE, force.clt = FALSE,
                 expand.residual.sd = expand)["replicates", ]$replicates
    }
    expect_gt(mean(apply(run(TRUE), 2, sd)), mean(apply(run(FALSE), 2, sd)))
  })
})

test_that("RES02: NNS.meboot drift: flat baseline, target drift and scale, kappa, force.clt spread, mean level, rounding tolerance", {
  withr::local_seed(20261108)
  meboot_drift <- fx_meboot_drift()
  aud_case("drift = FALSE replicates are a flat baseline plus residual-scale noise", env = meboot_drift, {
    R <- meboot_reps(drift = FALSE)
    expect_lt(abs(mean_slope(R)), 0.1)
    expect_lt(mean_detrended_sd(R), 1.1 * sd(res))
    expect_lt(mean(apply(R, 2, sd)), 1.1 * sd(res))
  })
  aud_case("target_drift and target_drift_scale set the replicate slope", env = meboot_drift, {
    expect_equal(mean_slope(meboot_reps(target_drift = 1)), 1, tolerance = 0.05)
    expect_equal(mean_slope(meboot_reps(target_drift_scale = 0.5)), 0.5 * coef(lm(x ~ t))[[2]], tolerance = 0.05)
    expect_lt(mean_detrended_sd(meboot_reps(target_drift = 1)), 1.1 * sd(res))
  })
  aud_case("scl.adjustment kappa uses the residual mean and sd", env = meboot_drift, {
    set.seed(1)
    b <- NNS.meboot(x, reps = 50, rho = 1, scl.adjustment = TRUE)
    zz <- c(b["xmin", ]$xmin, b["z", ]$z, b["xmax", ]$xmax)
    uv <- (sum((b["desintxb", ]$desintxb - mean(res))^2) + sum(diff(zz^2) / 12)) / length(x)
    expect_equal(b["kappa", ]$kappa, sd(res) / sqrt(uv) - 1, tolerance = 1e-10)
    expect_equal(mean(apply(b["replicates", ]$replicates, 2, sd)),
                 mean(apply(meboot_reps(), 2, sd)), tolerance = 0.01)
  })
  aud_case("force.clt spreads changed-drift replicate means on the residual scale", {
    set.seed(7)
    tt <- 1:100
    y <- 1000 + 10 * tt + rnorm(100)
    fit <- lm(y ~ tt)
    clt_reps <- function(...) {
      set.seed(1)
      NNS.meboot(y, reps = 10, rho = 1, ...)["replicates", ]$replicates
    }

    # The spread is that of a mean of n = 100 observations (#101).
    flat <- colMeans(clt_reps(drift = FALSE))
    expect_equal(mean(flat), mean(y), tolerance = 1e-6)
    expect_equal(sd(flat), sd(resid(fit)) / sqrt(100), tolerance = 1e-5)

    sloped <- colMeans(clt_reps(target_drift = 5))
    ref <- 5 * tt + resid(fit)
    expect_equal(mean(sloped), mean(y), tolerance = 1e-6)
    expect_equal(sd(sloped), sd(ref) / sqrt(100), tolerance = 1e-5)
  })
  aud_case("drift = FALSE keeps mean(x) as the replicate level", env = meboot_drift, {
    for (clt in c(TRUE, FALSE)) {
      R <- meboot_reps(drift = FALSE, force.clt = clt)
      expect_equal(mean(R), mean(x), tolerance = 1e-3)
      expect_gt(min(R), 0)
    }
  })
  aud_case("a target_drift equal to the original slope up to rounding keeps the default output", {
    ap <- as.numeric(AirPassengers); tt <- seq_along(ap)
    d_lm <- unname(coef(lm(ap ~ tt))[2])   # differs from fast_lm's slope by ~1e-15
    reps <- function(...) { set.seed(1); NNS.meboot(ap, reps = 20, rho = 1, ...)["replicates", ][[1]] }
    base <- reps()
    expect_equal(reps(target_drift = d_lm), base)
    expect_equal(reps(target_drift_scale = 1), base)
    # A real change of drift still takes the changed-drift path.
    expect_false(isTRUE(all.equal(reps(target_drift = 1.01 * d_lm), base)))
    # A target strictly inside the tolerance reproduces the default exactly,
    # which pins the snap to the original slope.
    d <- NNS:::fast_lm(tt, ap)$coef[2]
    tol <- sqrt(.Machine$double.eps) * max(abs(d), sd(ap) / length(ap))
    expect_identical(reps(target_drift = d + tol / 2), base)
  })
  aud_case("the slope tolerance allows for rounding at a large series level", env = meboot_drift, {
    t <- 1:20; set.seed(7); x <- 1e6 + .0001 * t + .01 * rnorm(20)
    d_lm <- unname(coef(lm(x ~ t))[2])   # same OLS fit, ~1e-11 from fast_lm's slope
    reps <- function(...) { set.seed(1); NNS.meboot(x, reps = 5, rho = .5, ...)["replicates", ][[1]] }
    expect_identical(reps(target_drift = d_lm), reps())
    # A real change of drift still takes the changed-drift path.
    expect_false(isTRUE(all.equal(reps(target_drift = .1), reps())))
  })
  aud_case("force.clt spreads replicate means like means of n observations (#101)", env = meboot_drift, {
    set.seed(1); x <- cumsum(rnorm(100))
    spread <- function(r) { set.seed(2); sd(colMeans(NNS.meboot(x, reps = r, rho = 1)["replicates", ]$replicates)) }
    s20 <- spread(20); s80 <- spread(80)
    # Before: 0.841 at 20 replicates and 0.420 at 80 (sd(x) / sqrt(reps)).
    expect_lt(abs(s80 / s20 - 1), 0.2)
    expect_equal(s80, sd(x) / sqrt(length(x)), tolerance = 0.2)
  })
})

test_that("RES03: degenerate series: exactly linear and constant series bootstrap, and NNS.SS and NNS.MC accept them; NNS.SS with constant groups", {
  withr::local_seed(20261108)
  meboot_zero_residuals <- fx_meboot_zero_residuals()
  aud_case("an exactly linear series bootstraps like a nearly linear one at any rho", env = meboot_zero_residuals, {
    x <- as.numeric(1:10)
    set.seed(5)
    near <- x + 1e-9 * rnorm(10)
    q <- qnorm(1:5 / 6)
    clt_means <- 5.5 + (q - mean(q)) / sd(q) * sd(x) / sqrt(10)
    for (type in c("spearman", "pearson", "nnsdep")) {
      for (rho in c(-1, 0, 0.5, 1)) {
        R <- reps_of(x, rho = rho, type = type)
        expect_lt(max_resid(R), 1e-6)
        expect_true(all(slopes(R) > 1 - 1e-6 & slopes(R) < 1.05 + 1e-6))
        expect_equal(unname(sort(colMeans(R))), clt_means, tolerance = 1e-6)
        expect_equal(unname(sort(colMeans(R))),
                     unname(sort(colMeans(reps_of(near, rho = rho, type = type)))), tolerance = 1e-6)
      }
    }
    expect_true(all(is.finite(reps_of(x, rho = 1, scl.adjustment = TRUE))))
    expect_replicate_matrix(reps_of(x, rho = 1, reps = 1), matrix(x, 10, 1))
    expect_equal(unname(colMeans(reps_of(x, rho = 1, drift = FALSE))), rep(5.5, 5), tolerance = 1e-6)

    r <- reps_of(ts(x, frequency = 4), rho = 0.5)
    expect_equal(attr(r, "tsp"), tsp(ts(x, frequency = 4)))
    expect_lt(max_resid(r), 1e-6)
  })
  aud_case("a constant series bootstraps to constant columns", env = meboot_zero_residuals, {
    expect_replicate_matrix(reps_of(rep(3, 10), rho = 0.5), matrix(3, 10, 5))
    expect_replicate_matrix(reps_of(rep(3, 10), rho = 0.5, type = "nnsdep"), matrix(3, 10, 5))
    expect_replicate_matrix(reps_of(rep(3, 10), rho = 0.5, scl.adjustment = TRUE), matrix(3, 10, 5))
    expect_replicate_matrix(reps_of(rep(3, 10), target_drift = 0), matrix(3, 10, 5))
  })
  aud_case("NNS.SS and NNS.MC work on an exactly linear series", {
    set.seed(1)
    x <- rnorm(10, 5)
    set.seed(2)
    r <- NNS.SS(x, 1:10, confidence.interval = TRUE, reps = 49)
    expect_true(all(is.finite(r$boot_vals)))
    expect_true(r$lower <= r$p_star && r$p_star <= r$upper)

    set.seed(1)
    mc <- NNS.MC(1:20, reps = 3, by = 0.5)
    expect_length(mc$replicates, 5)
    expect_true(all(is.finite(unlist(mc$replicates))))
  })
  aud_case("NNS.SS gives an interval when one group is constant", {
    x <- rep(1, 30)
    y <- rep(c(0, 1), 15)
    set.seed(1)
    r <- NNS.SS(x, y, confidence.interval = TRUE, reps = 99)
    expect_equal(r$p_star, 0.75)
    expect_true(all(is.finite(r$boot_vals)))
    expect_true(r$lower <= 0.75 && 0.75 <= r$upper)

    set.seed(1)
    r <- NNS.SS(y, x, confidence.interval = TRUE, reps = 99)
    expect_equal(r$p_star, 0.25)
    expect_true(r$lower <= 0.25 && 0.25 <= r$upper)
  })
  aud_case("NNS.SS interval for two constant groups is the point estimate", {
    r <- NNS.SS(rep(2, 10), rep(1, 15), confidence.interval = TRUE, reps = 9)
    expect_equal(c(r$lower, r$upper), c(1, 1))
    r <- NNS.SS(rep(1, 10), rep(1, 15), confidence.interval = TRUE, reps = 9)
    expect_equal(c(r$lower, r$upper), c(0.5, 0.5))
  })
  aud_case("NNS.SS bootstraps non-constant groups with NNS.meboot as before", {
    set.seed(42)
    x <- rnorm(40, 0.3)
    y <- rnorm(30)
    set.seed(1)
    r <- NNS.SS(x, y, confidence.interval = TRUE, reps = 49)
    set.seed(1)
    xb <- NNS.meboot(x, reps = 49, rho = 1)["replicates", ]$replicates
    yb <- NNS.meboot(y, reps = 49, rho = 1)["replicates", ]$replicates
    vals <- vapply(1:49, function(i) NNS.SS(xb[, i], yb[, i])$p_star, numeric(1))
    expect_identical(r$boot_vals, vals)
  })
})
