# File-level objects of the upstream tests, one constructor per former test file.
# A block calls the constructor and runs that file's cases in the environment it
# returns (aud_case(env = )), so none of these objects is a top-level name.

# ---- test-anova.R ------------------------------------------------------------

# former test-anova-grand-mean-medians.R
fx_anova_grand_mean_medians <- function() {
  a <- c(-1.2, 0.3, 0.8, 2.5, 1.1, -0.4)
  b <- c(0.1, 1.4, 2.2, 3.0, 3.9, 4.6, 0.5, 2.7, 1.9, -0.3)
  x <- 1:9
  y <- c(4, 7, 15, 20)
  environment()
}

# former test-anova-robust.R
fx_anova_robust <- function() {
  anova <- function(..., seed = 1) {
    set.seed(seed)
    tryCatch(NNS.ANOVA(..., plot = FALSE), error = conditionMessage)
  }
  x_na <- c(-2, NA, 0, NaN, 1, 3)
  y_na <- c(NA, -1, 1, 2)
  x_ok <- c(-2, 0, 1, 3)
  y_ok <- c(-1, 1, 2)
  x <- c(-2, -1, 0, 2, 6, 8, 10)
  y <- c(-3, 1, 4, 5)
  boot_stats <- function(seed) {
    set.seed(seed)
    cb <- matrix(sample(x, 1000 * length(x), replace = TRUE), nrow = length(x))
    tb <- matrix(sample(y, 1000 * length(y), replace = TRUE), nrow = length(y))
    list(c = colMeans(cb), t = colMeans(tb))
  }
  bounds <- function(r) {
    if (is.list(r)) unname(c(r$Effect_Size_LB, r$Effect_Size_UB)) else r
  }
  environment()
}

# former test_ANOVA.R
fx_anova_ref <- function() {
  x <- nns_py$x
  y <- nns_py$y
  z <- nns_py$z
  A <- data.frame(cbind(x,y,z))
  R1 <- c("Certainty" = 0.7642063)
  R2 <- matrix(c(
    1.0000000,
    0.7776676,
    0.7790700,
    0.7776676,
    1.0000000,
    0.9487158,
    0.7790700,
    0.9487158,
    1.0000000
  ),ncol=3)
  colnames(R2) <- c("x", "y", "z")
  rownames(R2) <- c("x", "y", "z")
  B <- NNS::NNS.ANOVA(cbind(x,y,z))
  C <- NNS::NNS.ANOVA(cbind(x,y,z), pairwise=T)
  environment()
}

# ---- test-arma.R -------------------------------------------------------------

# former test-arma-forecasts.R
fx_arma_forecasts <- function() {
  ap <- as.numeric(AirPassengers)
  arma <- function(x, ...) NNS.ARMA(x, plot = FALSE, ...)
  one_step_path <- function(x, h, ...) {
    out <- numeric(h)
    for (i in seq_len(h)) {
      out[i] <- arma(x, h = 1, ...)
      x <- c(x, out[i])
    }
    out
  }
  environment()
}

# former test-arma-nonlin-shrink.R
fx_arma_nonlin_shrink <- function() {
  arma <- function(method, shrink = FALSE, h = 1) {
    NNS.ARMA(AirPassengers, h = h, seasonal.factor = 12, method = method, shrink = shrink, plot = FALSE)
  }
  environment()
}

# former test-arma-optim-objective.R
fx_arma_optim_objective <- function() {
  mse <- expression(mean((predicted - actual)^2))
  run_optim <- function(...) {
    NNS.ARMA.optim(AirPassengers[1:132], training.set = 105, seasonal.factor = c(6, 12, 24),
                   obj.fn = mse, ncores = 1, print.trace = FALSE, ...)
  }
  environment()
}

# former test-arma-optim-scoring.R
fx_arma_optim_scoring <- function() {
  mse <- expression(mean((predicted - actual)^2))
  default_score <- function(predicted, actual) eval(eval(formals(NNS.ARMA.optim)$obj.fn))
  random_walk <- function(seed) {
    set.seed(seed)
    30 + cumsum(rnorm(100))
  }
  two_cycles <- function() {
    set.seed(9)
    t <- 1:160
    50 + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
  }
  run_optim <- function(...) NNS.ARMA.optim(..., ncores = 1, print.trace = FALSE)
  validation_forecast <- function(o, actual, negative.values = FALSE) {
    p <- o$errors + actual + o$bias.shift
    if (negative.values) p else pmax(0, p)
  }

  # A deterministic branch fixture tests smoothing data flow independently of
  # changes to gravity or the empirically selected model on one random seed.
  smoothing_fixture <- function(x, h = NULL) {
    calls <- integer()
    testthat::local_mocked_bindings(
      NNS.reg = function(x, y, ...) {
        calls <<- c(calls, length(y))
        list(Fitted.xy = data.frame(y.hat = rep(mean(y), length(y))))
      },
      NNS.ARMA = function(variable, h, training.set = NULL, ...) {
        if (!is.null(training.set)) variable <- head(variable, training.set)
        rep(tail(variable, 1), h)
      },
      .package = "NNS"
    )
    result <- run_optim(x, h = h, training.set = 80, seasonal.factor = 1,
                        obj.fn = mse, bias.correction = "none")
    list(result = result, calls = calls)
  }
  trend_cycles <- function() {
    set.seed(1)
    t <- 1:160
    50 + 0.05 * t + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
  }
  shrink_vs_equal_weights <- function(objective, obj.fn) {
    x <- trend_cycles()
    o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), lin.only = TRUE, objective = objective, obj.fn = obj.fn)
    p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = "lin",
                  weights = c(0.5, 0.5), plot = FALSE)
    expect_equal(o$weights, c(0.5, 0.5))
    expect_false(o$shrink)
    expect_equal(o$errors, p - tail(x, 32))
  }
  environment()
}

# former test-arma-pi-lower-bound.R
fx_arma_pi_lower_bound <- function() {
  arma_pi <- function(x, ...) {
    set.seed(123)
    NNS.ARMA(x, plot = FALSE, pred.int = 0.95, h = 6, ...)
  }
  environment()
}

# former test-arma-single-lag-shrink.R
fx_arma_single_lag_shrink <- function() {
  shrunk_lin <- function(h, ...) {
    NNS.ARMA(AirPassengers, h = h, method = "lin", shrink = TRUE, plot = FALSE, ...)
  }
  shrunk_recursion <- function(x, lag, h) {
    for (j in seq_len(h)) {
      y <- rev(x[seq(length(x) + 1 - lag, 1, -lag)])
      i <- seq_along(y)
      fit <- coef(lm(y ~ i))
      x <- c(x, max(0, (fit[[1]] + fit[[2]] * (length(y) + 1) + mean(y)) / 2))
    }
    tail(x, h)
  }
  environment()
}

# former test-arma-train-support.R
fx_arma_train_support <- function() {
  lin_arma <- function(x, ...) NNS.ARMA(x, h = 3, seasonal.factor = 1, method = "lin", plot = FALSE, ...)
  environment()
}

# former test-optim-pi-from-returned-model.R
fx_optim_pi_from_returned_model <- function() {
  mse <- expression(mean((predicted - actual)^2))
  two_cycles <- function() {
    set.seed(9)
    t <- 1:160
    50 + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
  }
  run_optim <- function(...) NNS.ARMA.optim(..., ncores = 1, print.trace = FALSE)
  environment()
}

# ---- test-boost.R ------------------------------------------------------------

# former test-boost-plots.R
fx_boost_plots <- function() {
  .boost_plot_fixture <- function() {
    set.seed(42)
    list(
      learner.results = runif(40, 0.5, 1),
      threshold.probability = 0.80,
      learner.threshold = 0.9,
      feature.frequency = c(a = 5L, b = 3L, c = 1L)
    )
  }
  environment()
}

# ---- test-classification.R ---------------------------------------------------

# former test-dimred-threshold.R
fx_dimred_threshold <- function() {
  kept <- function(equation) {
    vapply(equation, function(eq) {
      sum(abs(eq$Coefficient[eq$Variable != "DENOMINATOR"]) > 0)
    }, integer(1))
  }
  environment()
}

# ---- test-dependence.R -------------------------------------------------------

# former test-copula-marginal-anchor.R
fx_copula_marginal_anchor <- function() {
  grid_X <- function() {
    as.matrix(expand.grid(c(0.1, 0.2, 0.3, 0.5, 1.7, 4.1),
                          c(1, 1.5, 2, 9),
                          c(0.05, 0.4, 0.45, 3.3, 7)))
  }
  environment()
}

# former test_Copula.R
fx_copula_ref <- function() {
  x <- nns_py$x
  y <- nns_py$y
  z <- nns_py$z
  A <- data.frame(x,y)
  Z <- data.frame(x,y,z)
  B <- NNS.copula(A, continuous=T, plot=F)
  C <- NNS.copula(A, continuous=F, plot=F)
  D <- NNS.copula(Z, continuous=T, plot=F)
  E <- NNS.copula(Z, continuous=F, plot=F)
  environment()
}

# ---- test-derivatives.R ------------------------------------------------------

# former test-dydx-duplicate-bands.R
fx_dydx_duplicate_bands <- function() {
  band_means <- function(x, y, p) {
    h_s <- round(exp(seq(log(2), log(floor(sqrt(length(x)))), length.out = 5)))
    h_step <- gravity(diff(sort(unique(x)))) * h_s
    lo <- pmax(min(x), p - h_step)
    hi <- pmin(max(x), p + h_step)
    est <- matrix(unlist(NNS.reg(x, y, plot = FALSE, point.est = c(lo, rep(p, 5), hi),
                                 point.only = TRUE, ncores = 1, smooth = TRUE)$Point.est), ncol = 3)
    run_1 <- hi - p; run_2 <- p - lo
    rise_1 <- est[, 3] - est[, 2]; rise_2 <- est[, 2] - est[, 1]
    c(mean((rise_1 + rise_2) / (run_1 + run_2)),
      mean((rise_1 / run_1 - rise_2 / run_2) / ((run_1 + run_2) / 2)))
  }
  environment()
}

# ---- test-normalization.R ----------------------------------------------------

# former test_Normalization.R
fx_normalization_ref <- function() {
  x <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
  y <- c(10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0)
  z <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)
  A <- cbind(x, y, z)
  environment()
}

# ---- test-partial-moments-ratios-cdf.R ---------------------------------------

# former test-cdf-moments-var.R
fx_cdf_moments_var <- function() {
  strict_joint_survival <- function(A, t) mean(rowSums(sweep(A, 2, t, ">")) == ncol(A))
  environment()
}

# former test-cdf-target-value.R
fx_cdf_target_value <- function() {
  cdf_types <- c("CDF", "survival", "hazard", "cumulative hazard")
  cdf_target_data <- function() {
    set.seed(11)
    list(grid = 1:20, skewed = sort(round(rlnorm(30), 2)))
  }
  cdf_at <- function(x, degree, type, target) {
    NNS.CDF(x, degree, target = target, type = type, plot = FALSE)$target.value
  }
  environment()
}

# former test-pm-ratio-targets.R
fx_pm_ratio_targets <- function() {
  x <- c(1, 2, 3, 10)
  environment()
}

# ---- test-partial-moments-univariate.R ---------------------------------------

# former test-native-kernel-fixes.R
fx_native_kernel_fixes <- function() {
  direct_lpm <- function(d, t, x) vapply(t, function(tt) LPM(d, tt, x), numeric(1))
  direct_upm <- function(d, t, x) vapply(t, function(tt) UPM(d, tt, x), numeric(1))
  environment()
}

# ---- test-partial-moments-var.R ----------------------------------------------

# former test-var-root-tolerance.R
fx_var_root_tolerance <- function() {
  var_root_data <- function() {
    set.seed(42)
    list(
      symmetric = c(-.003, -.002, -.001, .001, .002, .003),
      skewed = c(rexp(40) - 1, 5),
      returns = rnorm(500, 0, .01)
    )
  }
  environment()
}

# ---- test-regression.R -------------------------------------------------------

# former test_Partition_Map.R
fx_partition_map_ref <- function() {
  x <- nns_py$x
  y <- nns_py$y
  T_ORDER <- 2
  T_DT <- data.frame(x, y, quadrant = "q", prior.quadrant = "pq", stringsAsFactors = FALSE)
  T_DT$quadrant <- c("q11","q44","q44","q12","q34","q23","q31","q13","q41","q21","q23","q13","q43","q42","q43","q13","q24","q22","q32","q12","q12","q11","q34","q34","q34","q41","q41","q42","q42","q12","q22","q23","q41",
                     "q41","q23","q44","q21","q31","q11","q41","q32","q42","q43","q21","q31","q44","q23","q31","q32","q14","q22","q11","q12","q14","q21","q24","q41","q34","q33","q32","q14","q34","q32","q34","q33","q42",
                     "q13","q22","q44","q12","q24","q11","q34","q34","q42","q12","q14","q22","q22","q13","q43","q32","q14","q41","q11","q33","q43","q24","q42","q21","q13","q33","q41","q11","q12","q11","q44","q43","q44",
                     "q23")
  T_DT$prior.quadrant <- c("q1","q4","q4","q1","q3","q2","q3","q1","q4","q2",
                           "q2","q1","q4","q4","q4","q1","q2","q2","q3","q1",
                           "q1","q1","q3","q3","q3","q4","q4","q4","q4","q1",
                           "q2","q2","q4","q4","q2","q4","q2","q3","q1","q4",
                           "q3","q4","q4","q2","q3","q4","q2","q3","q3","q1",
                           "q2","q1","q1","q1","q2","q2","q4","q3","q3","q3",
                           "q1","q3","q3","q3","q3","q4","q1","q2","q4","q1",
                           "q2","q1","q3","q3","q4","q1","q1","q2","q2","q1",
                           "q4","q3","q1","q4","q1","q3","q4","q2","q4","q2",
                           "q1","q3","q4","q1","q1","q1","q4","q4","q4","q2")
  T_regression_points <- data.frame(
    "quadrant"= c("q1", "q2", "q3", "q4"),
    "x"=c( 0.6826845, 0.3294027, 0.7310836, 0.3137827),
    "y"=c( 0.7416389, 0.8025298, 0.2608277, 0.3627194),
    stringsAsFactors = FALSE
  )
  R1 <- NNS.part(x,y,Voronoi=FALSE,min.obs.stop=TRUE)
  environment()
}

# ---- test-resampling.R -------------------------------------------------------

# former test-mc-rho-grid.R
fx_mc_rho_grid <- function() {
  mc_rhos <- function(...) {
    set.seed(3)
    m <- NNS.MC(rnorm(20), reps = 2, drift = FALSE, ...)
    as.numeric(sub("rho = ", "", names(m$replicates)))
  }
  environment()
}

# former test-meboot-drift.R
fx_meboot_drift <- function() {
  x <- as.numeric(AirPassengers)
  t <- seq_along(x)
  res <- resid(lm(x ~ t))
  meboot_reps <- function(...) {
    set.seed(1)
    NNS.meboot(x, reps = 50, rho = 1, ...)["replicates", ]$replicates
  }
  mean_slope <- function(R) mean(apply(R, 2, function(r) coef(lm(r ~ t))[2]))
  mean_detrended_sd <- function(R) mean(apply(R, 2, function(r) sd(resid(lm(r ~ t)))))
  environment()
}

# former test-meboot-zero-residuals.R
fx_meboot_zero_residuals <- function() {
  reps_of <- function(x, reps = 5, ...) {
    set.seed(1)
    NNS.meboot(x, reps = reps, ...)["replicates", ]$replicates
  }
  expect_replicate_matrix <- function(actual, expected) {
    expect_s3_class(actual, "NNS.matrix")
    expect_identical(dim(actual), dim(expected))
    expect_equal(as.numeric(actual), as.numeric(expected))
  }
  slopes <- function(R) apply(R, 2, function(r) coef(lm(r ~ seq_along(r)))[[2]])
  max_resid <- function(R) max(abs(apply(R, 2, function(r) resid(lm(r ~ seq_along(r))))))
  environment()
}

# ---- test-seasonality.R ------------------------------------------------------

# former test-seas-all-phases.R
fx_seas_all_phases <- function() {
  gmd <- function(v) {
    k <- length(v)
    if (k < 2 || min(v) == max(v)) return(0)
    s <- sort(v)
    2 * sum((2 * seq_len(k) - k - 1) * s) / (k * (k - 1))
  }
  detrend <- function(x) { t <- seq_along(x); x - unname(coef(lm(x ~ t))[2]) * (t - mean(t)) }
  score <- function(x, m) {
    n <- length(x); use_cv <- all(x > 0) || all(x < 0)
    sum(sapply(1:m, function(j) {
      v <- x[seq(j, n, m)]
      length(v) * (if (use_cv) gmd(v) / abs(mean(v)) else gmd(v))
    })) / n
  }
  environment()
}

# ---- test-stack.R ------------------------------------------------------------

# former test-stack-boost.R
fx_stack_boost <- function() {
  capture_frame <- function(fun, call) {
    captured <- new.env()
    suppressMessages(trace(
      fun, where = asNamespace("NNS"), print = FALSE,
      exit = bquote(assign("frame", environment(), envir = .(captured)))
    ))
    on.exit(suppressMessages(untrace(fun, where = asNamespace("NNS"))))
    force(call)
    captured$frame
  }
  boost_binary <- function(lo, hi) {
    y <- ifelse(iris$Species[1:140] == "setosa", hi, lo)
    NNS.boost(iris[1:140, 1:4], y, IVs.test = iris[c(1:3, 141:143), 1:4],
              epochs = 2, learner.trials = 3, feature.importance = FALSE,
              status = FALSE, seed = 1, pred.int = 0.9)
  }
  balance_y <- c(1, 1, 1, 2, 2, 2, 3)
  environment()
}

# former test-stack-pipeline-validation.R
fx_stack_pipeline_validation <- function() {
  .capture_pipeline_stack <- function(...) {
    captured <- new.env()
    suppressMessages(trace("NNS.stack", where = asNamespace("NNS"), print = FALSE,
      exit = bquote(assign("frame", environment(), envir = .(captured)))))
    on.exit(suppressMessages(untrace("NNS.stack", where = asNamespace("NNS"))))
    result <- NNS::NNS.stack(...)
    list(result = result, frame = captured$frame)
  }
  environment()
}

# former test-stack-univariate-k.R
fx_stack_univariate_k <- function() {
  capture_frame <- function(fun, call) {
    captured <- new.env()
    suppressMessages(trace(
      fun, where = asNamespace("NNS"), print = FALSE,
      exit = bquote(assign("frame", environment(), envir = .(captured)))
    ))
    on.exit(suppressMessages(untrace(fun, where = asNamespace("NNS"))))
    force(call)
    captured$frame
  }
  wt <- mtcars[, "wt", drop = FALSE]
  environment()
}

# ---- test-var.R --------------------------------------------------------------

# former test-var-forecast-rows.R
fx_var_forecast_rows <- function() {
  var_fixture <- function(n = 40) {
    set.seed(251)
    t <- seq_len(n)
    cbind(x = sin(t / 5) + rnorm(n, sd = .2),
          y = cos(t / 7) + rnorm(n, sd = .2))
  }

  # Runs NNS.VAR and records the IVs.test passed to each multivariate NNS.stack.
  var_with_test_rows <- function(variables, h, tau) {
    tests <- list()
    stack <- NNS.stack
    local_mocked_bindings(
      NNS.stack = function(IVs.train, DV.train, IVs.test = NULL, ...) {
        tests[[length(tests) + 1L]] <<- IVs.test
        stack(IVs.train, DV.train, IVs.test, ...)
      },
      .package = "NNS"
    )
    res <- NNS.VAR(variables, h = h, tau = tau, ncores = 1, status = FALSE)
    list(result = res, tests = tests)
  }

  # Lag-matrix rows n+1..n+h: column "<var>_tau_<k>" at period t is the value at t - k.
  forecast_rows <- function(res, cols, h) {
    extended <- rbind(as.data.frame(res$interpolated_and_extrapolated),
                      as.data.frame(res$univariate))
    n <- nrow(extended) - h
    vapply(cols, function(col) {
      k <- as.integer(sub(".*_tau_", "", col))
      extended[n + seq_len(h) - k, sub("_tau_[0-9]+$", "", col)]
    }, numeric(h))
  }
  expect_forecast_rows <- function(run, h, n_stacks) {
    expect_length(run$tests, n_stacks)
    for (iv in run$tests) {
      expected <- forecast_rows(run$result, colnames(iv), h)
      expect_equal(unname(as.matrix(iv)), unname(matrix(expected, nrow = h)))
    }
  }
  environment()
}
