# NNS.ARMA and NNS.ARMA.optim: seasonal means and weights, dynamic forecasts,
# shrink, the optimiser's scores and stages, and prediction intervals.
#
# NNS.ARMA.optim run time: the same models with less repeated work.

test_that("ARMA01: NNS.ARMA argument contract: weights, h and the warn option", {
  withr::local_seed(20261108)
  arma_forecasts <- fx_arma_forecasts()
  aud_case("weights need one value per seasonal lag", env = arma_forecasts, {
    expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = 1),
                 "one value per seasonal lag")
  })
  aud_case("weights must be finite and must not sum to 0", env = arma_forecasts, {
    expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = c(0, 0)), "finite")
    expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = c(NA, 1)), "finite")
  })
  aud_case("dynamic = TRUE rejects numeric weights", env = arma_forecasts, {
    expect_error(arma(ap, h = 2, seasonal.factor = FALSE, best.periods = 2, dynamic = TRUE,
                      weights = c(1, 0)),
                 "dynamic")
  })
  aud_case("h must be a positive integer", env = arma_forecasts, {
    for (h in list(0, 2.5, -1, NA, c(1, 2))) {
      expect_error(arma(ap, h = h, seasonal.factor = 12), "positive integer")
    }
    expect_length(arma(ap, h = 3L, seasonal.factor = 12), 3)
  })
  aud_case("the warn option is restored when NNS.ARMA fails", env = arma_forecasts, {
    old <- getOption("warn")
    on.exit(options(warn = old), add = TRUE)
    options(warn = 1)
    expect_error(arma(ap[1:3], h = 2, seasonal.factor = 12))
    expect_equal(getOption("warn"), 1)
    expect_error(arma(ap, h = 2, seasonal.factor = c(12, 24), weights = 1))
    expect_equal(getOption("warn"), 1)
  })
})

test_that("ARMA02: method = 'means' and seasonal weights: subset means, relative and scale-free weights, equal and zero weights, modulo periods, intervals", {
  withr::local_seed(20261108)
  arma_forecasts <- fx_arma_forecasts()
  aud_case("method = 'means' forecasts the seasonal subset means", env = arma_forecasts, {
    monthly <- as.numeric(tapply(ap, rep(1:12, 12), mean))
    expect_equal(arma(ap, h = 14, seasonal.factor = 12, method = "means"),
                 rep(monthly, length.out = 14))
  })
  aud_case("method = 'means' with several lags weights each lag's subset mean", env = arma_forecasts, {
    expected <- 0.25 * mean(ap[seq(1, 144, 12)]) + 0.75 * mean(ap[seq(1, 144, 24)])
    expect_equal(arma(ap, h = 1, seasonal.factor = c(12, 24), weights = c(1, 3),
                      method = "means"),
                 expected)
  })
  aud_case("forecasts do not depend on the scale of the weights", env = arma_forecasts, {
    for (m in c("lin", "nonlin", "both")) {
      expect_equal(arma(ap, h = 3, seasonal.factor = c(12, 24), weights = c(1, 1), method = m),
                   arma(ap, h = 3, seasonal.factor = c(12, 24), weights = c(0.5, 0.5), method = m))
    }
  })
  aud_case("a single-lag weight is relative", env = arma_forecasts, {
    expect_equal(arma(ap, h = 3, seasonal.factor = 12, weights = 0.5),
                 arma(ap, h = 3, seasonal.factor = 12))
  })
  aud_case("weights = 'equal' works with a numeric seasonal.factor", env = arma_forecasts, {
    expect_equal(arma(ap, h = 3, seasonal.factor = c(12, 24), weights = "equal"),
                 arma(ap, h = 3, seasonal.factor = c(12, 24), weights = c(1, 1)))
  })
  aud_case("all-zero seasonal weights fall back to equal weights", env = arma_forecasts, {
    set.seed(1)
    x <- rnorm(100, 10)
    expect_equal(arma(x, h = 2, modulo = 7),
                 arma(x, h = 2, modulo = 7, weights = "equal"))
  })
  aud_case("method = 'means' returns prediction intervals", env = arma_forecasts, {
    set.seed(123)
    r <- arma(ap, h = 4, seasonal.factor = 12, method = "means", pred.int = 0.95)
    expect_equal(dim(r), c(4L, 3L))
    expect_true(all(is.finite(as.matrix(r))))
    expect_true(all(r[[2]] <= r[[1]] & r[[1]] <= r[[3]]))
  })
  aud_case("a modulo period beyond the n/5 scan keeps its ARMA weight", {
    set.seed(1)
    x <- 10 + sin(2 * pi * (1:96) / 24) + rnorm(96, sd = .1)
    s <- NNS.seas(x, modulo = 24, mod.only = FALSE, plot = FALSE)
    w <- NNS:::ARMA.seas.weighting(FALSE, s$all.periods)
    expect_gt(w$Weights[match(24, w$lag)], 0)
    expect_lt(s$all.periods$Coefficient.of.Variation[s$all.periods$Period == 24],
              s$all.periods$Variable.Coefficient.of.Variation[1])
  })
})

test_that("ARMA03: constant, trending and negative series: constant forecasts, best period, nonseasonal trend from lag 1, negative values kept", {
  withr::local_seed(20261108)
  arma_forecasts <- fx_arma_forecasts()
  aud_case("nonlinear forecasts keep negative values when the series has them", env = arma_forecasts, {
    set.seed(1)
    s <- 3 * sin(2 * pi * (1:240) / 12) + rnorm(240, 0, 0.5)
    f <- arma(s, h = 12, seasonal.factor = 12)
    expect_lt(max(f[7:11]), -0.5)
  })
  aud_case("seasonal.factor = TRUE with best.periods = NULL uses the best period", env = arma_forecasts, {
    best <- NNS.seas(ap, plot = FALSE)$best.period
    expect_equal(arma(ap, h = 3, seasonal.factor = TRUE, best.periods = NULL),
                 arma(ap, h = 3, seasonal.factor = best))
  })
  aud_case("a constant series forecasts the constant with automatic seasonality", env = arma_forecasts, {
    expect_equal(arma(rep(5, 60), h = 3), c(5, 5, 5))
  })
  aud_case("a constant series forecasts the constant with a numeric seasonal.factor", env = arma_forecasts, {
    for (m in c("lin", "nonlin", "both", "means")) {
      expect_equal(arma(rep(5, 60), h = 3, seasonal.factor = 12, method = m), c(5, 5, 5))
      expect_equal(arma(rep(5, 60), h = 3, seasonal.factor = c(12, 24), method = m), c(5, 5, 5))
      expect_equal(arma(rep(-2, 30), h = 2, seasonal.factor = 7, method = m), c(-2, -2))
    }
  })
  aud_case("NNS.ARMA forecasts a nonseasonal trend from lag 1", {
    f <- NNS.ARMA(1:100, h = 2, plot = FALSE, seasonal.plot = FALSE)
    expect_equal(f, NNS.ARMA(1:100, h = 2, seasonal.factor = 1, plot = FALSE))
    expect_equal(f, c(101, 102), tolerance = 1e-3)
  })
})

test_that("ARMA04: dynamic = TRUE: repeated one-step forecasts, modulo, later periods and equal weights at every step", {
  withr::local_seed(20261108)
  arma_forecasts <- fx_arma_forecasts()
  aud_case("dynamic = TRUE matches repeated one-step forecasts with best.periods", env = arma_forecasts, {
    expect_equal(arma(ap, h = 3, seasonal.factor = FALSE, best.periods = 2, dynamic = TRUE),
                 one_step_path(ap, 3, seasonal.factor = FALSE, best.periods = 2))
  })
  aud_case("dynamic = TRUE honours modulo", env = arma_forecasts, {
    expect_equal(arma(ap, h = 3, seasonal.factor = FALSE, best.periods = 2, modulo = 12,
                      dynamic = TRUE),
                 one_step_path(ap, 3, seasonal.factor = FALSE, best.periods = 2, modulo = 12))
  })
  aud_case("dynamic = TRUE uses periods that appear later in the forecast path", env = arma_forecasts, {
    set.seed(1)
    x <- rnorm(20, 10)
    expect_equal(arma(x, h = 6, seasonal.factor = FALSE, best.periods = 4, method = "lin",
                      dynamic = TRUE),
                 one_step_path(x, 6, seasonal.factor = FALSE, best.periods = 4, method = "lin"))
  })
  aud_case("dynamic = TRUE applies weights = 'equal' at every step", env = arma_forecasts, {
    expect_equal(arma(ap, h = 3, seasonal.factor = FALSE, best.periods = 2, dynamic = TRUE,
                      weights = "equal"),
                 one_step_path(ap, 3, seasonal.factor = FALSE, best.periods = 2, weights = "equal"))
  })
})

test_that("ARMA05: shrink: blend with the means estimate, both components, every forecast moves, affine level, single shrunk lag", {
  withr::local_seed(20261108)
  arma_nonlin_shrink <- fx_arma_nonlin_shrink()
  arma_single_lag_shrink <- fx_arma_single_lag_shrink()
  aud_case("shrink blends the nonlin estimate with the means estimate", env = arma_nonlin_shrink, {
    expect_equal(arma("nonlin", TRUE), (arma("nonlin") + arma("means")) / 2)
  })
  aud_case("shrink with method = \"both\" shrinks both components", env = arma_nonlin_shrink, {
    expect_equal(arma("both", TRUE), (arma("lin", TRUE) + arma("nonlin", TRUE)) / 2)
  })
  aud_case("shrink moves every nonlin forecast towards the means forecast", env = arma_nonlin_shrink, {
    nonlin <- arma("nonlin", h = 12)
    shrunk <- arma("nonlin", TRUE, h = 12)
    means <- arma("means", h = 12)
    expect_true(all(nonlin != means))
    expect_true(all(shrunk != nonlin))
    expect_true(all(shrunk > pmin(nonlin, means) & shrunk < pmax(nonlin, means)))
  })
  aud_case("nonlinear shrink preserves positive affine changes of level", {
    f <- function(x) NNS.ARMA(x, h = 1, seasonal.factor = 12,
                              method = "nonlin", shrink = TRUE, plot = FALSE)
    expect_equal(f(2 * AirPassengers + 7), 2 * f(AirPassengers) + 7, tolerance = 1e-6)
  })
  aud_case("a single shrunk lag refits on its own forecasts past the lag", env = arma_single_lag_shrink, {
    single <- shrunk_lin(30, seasonal.factor = 12)
    expect_equal(single, shrunk_lin(30, seasonal.factor = c(12, 24), weights = c(1, 0)))
    expect_equal(single, shrunk_recursion(as.numeric(AirPassengers), 12, 30))
  })
  aud_case("a single shrunk lag is unchanged up to the lag", env = arma_single_lag_shrink, {
    expect_equal(shrunk_lin(1, seasonal.factor = 12), 332.0568182, tolerance = 1e-6)
    expect_equal(shrunk_lin(12, seasonal.factor = 12), shrunk_recursion(as.numeric(AirPassengers), 12, 12))
  })
  aud_case("a single shrunk lag longer than the history forecasts NaN, not 0", {
    expect_equal(NNS.ARMA(rep(100, 11), seasonal.factor = 12, h = 3, method = "lin", shrink = TRUE, plot = FALSE),
                 c(NaN, 100, 100))
  })
})

test_that("ARMA06: NNS.ARMA.optim reports the score of the validation forecast it returns: bias shift, split-window correction, no default shift, numeric responses", {
  withr::local_seed(20261108)
  arma_optim_objective <- fx_arma_optim_objective()
  aud_case("obj.fn scores the returned validation forecast", env = arma_optim_objective, {
    o <- run_optim()
    expect_equal(o$validation.forecast, o$results)
    expect_equal(o$obj.fn, mean((o$results - AirPassengers[106:132])^2))
  })
  aud_case("obj.fn scores the bias-shifted validation forecast when h is set", env = arma_optim_objective, {
    o <- run_optim(h = 12)
    actual <- AirPassengers[106:132]
    validation <- pmax(0, o$errors + actual + o$bias.shift)
    expect_equal(o$obj.fn, mean((validation - actual)^2))
  })
  # Force a persistent bias, rather than relying on a random fixture selecting a
  # particular gravity shift after other forecasting changes.
  aud_case("an accepted split-window correction is included in the returned score", env = arma_optim_objective, {
    testthat::local_mocked_bindings(
      NNS.ARMA = function(variable, h, ...) rep(200, h),
      .package = "NNS"
    )
    o <- run_optim(bias.correction = "split")
    actual <- AirPassengers[106:132]
    expect_true(o$bias.validation$accepted)
    expect_true(o$bias.shift != 0)
    expect_equal(o$raw.results, rep(200, length(actual)))
    expect_equal(o$results, pmax(0, o$raw.results + o$bias.shift))
    expect_equal(o$obj.fn, mean((o$results - actual)^2))
  })
  aud_case("bias is estimated before the rows used to accept it", {
    score <- expression(mean((predicted - actual)^2))
    p <- rep(10, 8)
    choose <- function(y, policy = "split", objective = "min", obj = score)
      NNS:::.nns_arma_bias_choice(p, y, obj, objective, TRUE, policy)
    persistent <- choose(rep(7, 8))
    transient <- choose(c(rep(7, 4), rep(13, 4)))
    expect_equal(persistent$estimate.rows, 1:4)
    expect_equal(persistent$acceptance.rows, 5:8)
    expect_equal(persistent$estimated.shift, -3)
    expect_equal(transient$estimated.shift, -3)
    expect_true(persistent$accepted)
    expect_equal(persistent$applied.shift, -3)
    expect_false(transient$accepted)
    expect_equal(transient$applied.shift, 0)
    expect_equal(choose(rep(7, 8), objective = "max",
                        obj = expression(-mean((predicted - actual)^2)))$applied.shift, -3)
    expect_equal(choose(rep(7, 8), "none")$status, "disabled")
  })
  aud_case("short validation windows do not claim separate bias validation", {
    x <- NNS:::.nns_arma_bias_choice(4:6, 1:3, expression(mean((predicted - actual)^2)),
                                    "min", TRUE, "split")
    expect_equal(x$applied.shift, 0)
    expect_equal(x$status, "insufficient_validation")
  })
  aud_case("NNS.ARMA.optim adds no bias shift by default (#150)", {
    expect_identical(eval(formals(NNS.ARMA.optim)$bias.correction), c("none", "split"))
    # The same-window rule (estimate and accept on one window) is not offered.
    expect_error(NNS.ARMA.optim(as.numeric(AirPassengers[1:96]), h = 12, seasonal.factor = 12,
      print.trace = FALSE, ncores = 1, pred.int = NULL, bias.correction = "legacy"))
    fit <- NNS.ARMA.optim(as.numeric(AirPassengers[1:96]), h = 12,
      seasonal.factor = c(6, 12), print.trace = FALSE, ncores = 1, pred.int = NULL)
    expect_identical(fit$bias.shift, 0)
    expect_identical(fit$bias.validation$policy, "none")
    expect_equal(fit$results, fit$raw.results)
  })
  aud_case("ARMA nonlinear fits and optimizer smoothing retain numeric responses", {
    for (levels in list(c(0, 1), c(-5, 5))) {
      y <- rep(levels, 15)
      fit <- NNS.ARMA(y, h = 2, seasonal.factor = 3, method = "nonlin",
                      plot = FALSE, negative.values = TRUE)
      expect_equal(fit, levels[1] + diff(levels) * c(10 / 9, -1 / 9),
                   tolerance = 1e-10)
      optimized <- NNS.ARMA.optim(y, h = 2, seasonal.factor = 3,
        negative.values = TRUE, ncores = 1, plot = FALSE, print.trace = FALSE,
        obj.fn = expression(mean((predicted - actual)^2)), bias.correction = "none")
      expect_equal(optimized$results,
                   levels[1] + diff(levels) * c(7 / 12, 5 / 12), tolerance = 1e-10)
    }
  })
})

test_that("ARMA07: NNS.ARMA.optim stages: lin with negative values, nonlin, both, equal weights, and shrink against the equal-weights score", {
  withr::local_seed(20261108)
  arma_optim_scoring <- fx_arma_optim_scoring()
  aud_case("the lin model forecasts with negative.values", env = arma_optim_scoring, {
    t <- 1:160
    x <- pmax(0.1, 20 - 0.15 * t + sin(2 * pi * t / 12))
    o <- run_optim(x, training.set = 128, seasonal.factor = 12, lin.only = TRUE,
               negative.values = TRUE, obj.fn = mse)
    p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = 12, method = "lin",
                  negative.values = TRUE, plot = FALSE)
    expect_equal(o$errors, p - tail(x, 32))
    expect_equal(o$obj.fn, mean((validation_forecast(o, tail(x, 32), TRUE) - tail(x, 32))^2))
  })
  aud_case("linear.approximation = FALSE evaluates the nonlin method", env = arma_optim_scoring, {
    x <- AirPassengers[1:96]
    o <- run_optim(x, training.set = 84, seasonal.factor = c(6, 12), linear.approximation = FALSE)
    p <- NNS.ARMA(x, training.set = 84, h = 12, seasonal.factor = 12, method = "nonlin", plot = FALSE)
    expect_equal(o$method, "nonlin")
    expect_equal(as.numeric(o$periods), 12)
    expect_equal(o$errors, p - tail(x, 12))
    expect_equal(o$obj.fn, default_score(validation_forecast(o, tail(x, 12)), tail(x, 12)))
  })
  aud_case("the both stage is scored with NNS.ARMA(method = \"both\")", env = arma_optim_scoring, {
    x <- random_walk(5)
    o <- run_optim(x, training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
    p <- NNS.ARMA(x, training.set = 80, h = 20, seasonal.factor = 2, method = "both",
                  shrink = o$shrink, plot = FALSE)
    expect_equal(o$method, "both")
    expect_equal(o$errors, p - tail(x, 20))
    expect_equal(o$obj.fn, mean((validation_forecast(o, tail(x, 20)) - tail(x, 20))^2))
  })
  aud_case("equal weights are selected when they score better", env = arma_optim_scoring, {
    o <- run_optim(two_cycles(), h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
    expect_equal(as.numeric(o$periods), c(7, 12))
    expect_equal(o$weights, c(0.5, 0.5))
  })
  aud_case("obj.fn is the equal-weights score when equal weights win (min)", env = arma_optim_scoring, {
    x <- two_cycles()
    o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
    p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = o$method,
                  weights = c(0.5, 0.5), shrink = o$shrink, plot = FALSE)
    expect_equal(o$errors, p - tail(x, 32))
    expect_equal(o$obj.fn, mean((validation_forecast(o, tail(x, 32)) - tail(x, 32))^2))
  })
  aud_case("obj.fn is the equal-weights score when equal weights win (max)", env = arma_optim_scoring, {
    x <- two_cycles()
    o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), objective = "max",
               obj.fn = expression(-mean((predicted - actual)^2)))
    p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = o$method,
                  weights = c(0.5, 0.5), shrink = o$shrink, plot = FALSE)
    expect_equal(o$errors, p - tail(x, 32))
    expect_equal(o$obj.fn, -mean((validation_forecast(o, tail(x, 32)) - tail(x, 32))^2))
  })
  # Equal weights score 12.04 and no weights 12.67, so shrink at 12.40 must lose to the equal-weights score.
  aud_case("shrink is tested against the equal-weights score (min)", env = arma_optim_scoring, {
    shrink_vs_equal_weights("min", mse)
  })
  aud_case("shrink is tested against the equal-weights score (max)", env = arma_optim_scoring, {
    shrink_vs_equal_weights("max", expression(-mean((predicted - actual)^2)))
  })
})

test_that("ARMA08: NNS.ARMA.optim nns.regress: full-series forecast, training-window smooth, must beat shrink", {
  withr::local_seed(20261108)
  arma_optim_scoring <- fx_arma_optim_scoring()
  aud_case("the out-of-sample forecast applies nns.regress to the full series", env = arma_optim_scoring, {
    x <- c(rep(10, 79), 100, rep(12, 20))
    fit <- smoothing_fixture(x, h = 20)
    expect_true(fit$result$nns.regress)
    expect_true(fit$result$shrink)
    expect_equal(fit$calls, c(80L, 100L))
    expect_equal(fit$result$results, rep(mean(x), 20))
  })
  aud_case("nns.regress is not selected by a smooth that saw the validation window", env = arma_optim_scoring, {
    x <- c(rep(0, 79), 14.5, rep(8, 20))
    # A full-series mean would beat the last-value forecast; the prefix mean loses.
    expect_lt((mean(x) - 8)^2, (14.5 - 8)^2)
    expect_gt((mean(x[1:80]) - 8)^2, (14.5 - 8)^2)
    fit <- smoothing_fixture(x)
    expect_false(fit$result$nns.regress)
    expect_equal(fit$calls, 80L)
    expect_equal(fit$result$results, rep(14.5, 20))
  })
  aud_case("validation results with nns.regress use a smooth of the training window", env = arma_optim_scoring, {
    x <- c(rep(10, 79), 100, rep(12, 20))
    fit <- smoothing_fixture(x)
    expect_true(fit$result$nns.regress)
    expect_equal(fit$calls, 80L)
    expect_equal(fit$result$results, rep(mean(x[1:80]), 20))
  })
  aud_case("nns.regress must beat the accepted shrink score", env = arma_optim_scoring, {
    # Smoothed 1.670 beats raw 3.809 but loses to raw with shrink, 1.447.
    o <- run_optim(random_walk(5), training.set = 80, seasonal.factor = c(1, 2, 4, 6, 12), obj.fn = mse)
    expect_true(o$shrink)
    expect_false(o$nns.regress)
  })
})

test_that("ARMA09: NNS.ARMA.optim internals match their reference constructions: one forecast per validation, nonlinear step, partition frames, threaded PM, lags, scan slots", {
  withr::local_seed(20261108)
  aud_case("NNS.ARMA.optim computes each validation forecast once", {
    orig <- NNS::NNS.ARMA
    keys <- character()
    testthat::local_mocked_bindings(
      NNS.ARMA = function(variable, ...) {
        a <- list(...)
        keys <<- c(keys, paste(length(variable), deparse(a$training.set), a$h, a$method,
                               paste(a$seasonal.factor, collapse = ","),
                               paste(a$weights, collapse = ","), isTRUE(a$shrink),
                               sep = "|"))
        orig(variable, ...)
      },
      .package = "NNS"
    )
    set.seed(9)
    t <- 1:160
    x <- 50 + 5 * sin(2 * pi * t / 12) + 4 * sin(2 * pi * t / 7) + rnorm(160, sd = 1.5)
    fit <- NNS.ARMA.optim(x, h = 10, seasonal.factor = c(7, 12), ncores = 1,
                          print.trace = FALSE,
                          obj.fn = expression(mean((predicted - actual)^2)))
    # The selected multi-period model and its equal-weight candidate were both
    # requested twice before.
    expect_equal(as.numeric(fit$periods), c(7, 12))
    expect_equal(fit$weights, c(0.5, 0.5))
    expect_false(anyDuplicated(keys) > 0)
  })
  aud_case("the nonlinear step uses the squared-index weighted rise and run", {
    x <- as.numeric(AirPassengers)[1:120]
    for (lag in c(1, 12)) {
      gv <- NNS:::generate.vectors(x, lag)
      ci <- gv$Component.index[[1]]
      cs <- gv$Component.series[[1]]
      rp <- NNS.reg(ci, cs, type = "REG", return.values = FALSE, plot = FALSE,
                    multivariate.call = TRUE)
      xs <- rp$x[nrow(rp)] - rp$x
      ys <- rp$y[nrow(rp)] - rp$y
      xs <- xs[-length(xs)]
      ys <- ys[-length(ys)]
      expected <- tail(cs, 1) +
        mean(rep(ys, seq_along(ys)^2)) / mean(rep(xs, seq_along(xs)^2))
      expect_equal(NNS.ARMA(x, h = 1, seasonal.factor = lag, method = "nonlin",
                            plot = FALSE), expected, tolerance = 1e-12)
    }
  })
  aud_case("partition data frames are the ones as.data.frame() builds", {
    set.seed(3)
    x <- rnorm(60)
    y <- x + rnorm(60)
    for (type in list(NULL, "XONLY")) {
      out <- NNS:::NNS_part_cpp(x, y, type = type, order_in = 3L, obs_req = 8L,
                                min_obs_stop = TRUE, noise_reduction = "off")
      for (nm in c("dt", "regression.points", "segments_h", "segments_v")) {
        expect_identical(out[[nm]],
                         as.data.frame(as.list(out[[nm]]), stringsAsFactors = FALSE))
      }
    }
  })
  aud_case("small partial-moment inputs agree with the threaded path", {
    # 30 x 3 runs on the calling thread; 200 stacked copies are dispatched.
    set.seed(4)
    m <- matrix(rnorm(90), 30, 3)
    big <- m[rep(seq_len(30), 200), ]
    tg <- colMeans(m)
    for (deg in c(0, 1)) {
      expect_equal(DPM_nD(m, tg, deg, TRUE), DPM_nD(big, tg, deg, TRUE), tolerance = 1e-12)
      expect_equal(Co.LPM_nD(m, tg, deg, FALSE), Co.LPM_nD(big, tg, deg, FALSE),
                   tolerance = 1e-12)
      expect_equal(PM.matrix(deg, deg, tg, m, FALSE), PM.matrix(deg, deg, tg, big, FALSE),
                   tolerance = 1e-12)
    }
  })
  aud_case("a plain numeric predictor is encoded as its one-column data frame", {
    enc <- NNS:::.nns_reg_encode_predictors
    x <- c(a = 3, b = 1, c = 2)
    expect_identical(enc(x), enc(data.frame(x = x)))
    expect_identical(enc(1:5), enc(data.frame(x = 1:5)))
    expect_error(enc(c(1, NA)), "finite numeric")
  })
  aud_case("regression points match the data frame construction", {
    # The data frame version the vector code replaced.
    reference <- function(x, y, order, noise.reduction = "off") {
      reducer <- function(z) NNS:::.nns_reg_reduce_value(z, noise.reduction, FALSE)
      if (identical(order, "max")) {
        s <- split(y, x)
        rp <- data.frame(x = as.numeric(names(s)), y = vapply(s, reducer, numeric(1L)))
      } else {
        part <- NNS.part(x, y, type = "XONLY", order = order, obs.req = 0,
                         min.obs.stop = TRUE, noise.reduction = noise.reduction)
        rp <- as.data.frame(part$regression.points[, c("x", "y"), drop = FALSE])
      }
      rp <- rp[is.finite(rp$x) & is.finite(rp$y), , drop = FALSE]
      rp <- rp[order(rp$x, method = "radix"), , drop = FALSE]
      merge <- function(rp) {
        if (!anyDuplicated(rp$x)) return(rp)
        g <- split(rp$y, rp$x)
        rp <- data.frame(x = as.numeric(names(g)), y = vapply(g, reducer, numeric(1L)))
        rp[order(rp$x, method = "radix"), , drop = FALSE]
      }
      rp <- merge(rp)
      rp <- rbind(rp, data.frame(x = c(min(x), max(x)),
                                 y = c(reducer(y[x == min(x)]), reducer(y[x == max(x)]))))
      rp <- merge(rp[order(rp$x, method = "radix"), , drop = FALSE])
      rownames(rp) <- NULL
      rp
    }
    set.seed(5)
    cases <- list(list(rnorm(40), rnorm(40)),
                  list(round(runif(50) * 4), rnorm(50)),
                  list(1:30, cumsum(rnorm(30))))
    for (cs in cases) for (ord in list(2L, 4L, "max")) for (nr in c("off", "median")) {
      got <- NNS:::.nns_reg_build_points(cs[[1]], cs[[2]], ord, nr, FALSE)
      expect_identical(got[c("x", "y")], reference(cs[[1]], cs[[2]], ord, nr))
      expect_false(is.null(attr(got, "nns.order")))
    }
  })
  aud_case("linear.approximation defaults to TRUE", {
    # FALSE scores every candidate set with the nonlinear methods; the default
    # passes the linear search's selection to the nonlinear comparison.
    expect_identical(eval(formals(NNS.ARMA.optim)$linear.approximation), TRUE)
  })
  aud_case("several lags of at least h give the step-by-step result", {
    # The loop: each step regresses every lag's component series on its index,
    # combines the lags with sum(), and appends the forecast to the series.
    loop <- function(x, lags, h, w) {
      out <- numeric(h)
      for (j in seq_len(h)) {
        gv <- NNS:::generate.vectors(x, lags)
        est <- vapply(seq_along(lags), function(i) {
          ci <- gv$Component.index[[i]]
          cf <- NNS:::fast_lm(ci, gv$Component.series[[i]])$coef
          as.numeric(cf[1] + cf[2] * (ci[length(ci)] + 1))
        }, numeric(1))
        out[j] <- pmax(0, sum(est * w))
        x <- c(x, out[j])
      }
      out
    }
    x <- as.numeric(AirPassengers)
    # c(6, 24) has a lag shorter than h = 12 and keeps the loop.
    for (lags in list(c(12, 24), c(24, 36, 30), c(6, 24))) for (h in c(1, 6, 12)) {
      w <- rep(1 / length(lags), length(lags))
      expect_identical(NNS.ARMA(x, h = h, seasonal.factor = lags, method = "lin",
                                weights = "equal", plot = FALSE),
                       loop(x, lags, h, w))
    }
  })
  aud_case("the seasonality scan scores every lag in its own slot", {
    # 1200 lags of 6000 values are scored on several threads.
    set.seed(8)
    t <- 1:6000
    x <- 10 + sin(2 * pi * t / 24) + rnorm(6000, sd = .3)
    M <- NNS:::NNS_seas_cpp(variable = x, modulo = NULL, mod_only = TRUE,
                            max_period = NULL)$all.periods
    direct <- function(m) {
      comp <- split(x, (seq_along(x) - 1L) %% m)
      d <- vapply(comp, function(z) {
        k <- length(z)
        s <- sort(z)
        sum((2 * seq_len(k) - k - 1) * s) * 2 / (k * (k - 1)) / abs(mean(z))
      }, numeric(1))
      sum(lengths(comp) * d) / length(x)
    }
    for (m in c(24L, 48L, 600L, 1199L)) {
      expect_equal(M$Coefficient.of.Variation[M$Period == m], direct(m), tolerance = 1e-10)
    }
  })
})

test_that("ARMA10: prediction intervals: residuals of the returned model, LPM/UPM VaR bands, bounds around the forecast, one-point and flat paths, clamped nonnegative series", {
  withr::local_seed(20261108)
  optim_pi_from_returned_model <- fx_optim_pi_from_returned_model()
  arma_pi_lower_bound <- fx_arma_pi_lower_bound()
  arma_train_support <- fx_arma_train_support()
  aud_case("errors and bias.shift come from the returned shrink model", env = optim_pi_from_returned_model, {
    x <- two_cycles()
    o <- run_optim(x, training.set = 128, seasonal.factor = c(7, 12), obj.fn = mse,
                   bias.correction = "split")
    p <- NNS.ARMA(x, training.set = 128, h = 32, seasonal.factor = c(7, 12), method = "lin",
                  weights = c(0.5, 0.5), shrink = TRUE, plot = FALSE)
    expect_true(o$shrink)
    expect_equal(o$errors, p - tail(x, 32))
    expect_equal(o$bias.validation$estimated.shift,
                 -NNS.gravity(o$errors[o$bias.validation$estimate.rows]))
    expect_equal(o$bias.shift, o$bias.validation$applied.shift)
    expect_equal(o$results, p + o$bias.shift)
  })
  aud_case("errors come from the returned nns.regress model", env = optim_pi_from_returned_model, {
    x <- as.numeric(AirPassengers[1:132])
    o <- run_optim(x, training.set = 105, seasonal.factor = c(6, 12, 24), objective = "max")
    smoothed <- NNS.reg(1:105, x[1:105], plot = FALSE, smooth = TRUE)$Fitted.xy$y.hat
    p <- NNS.ARMA(smoothed, training.set = 105, h = 27, seasonal.factor = o$periods,
                  method = o$method, weights = o$weights, shrink = TRUE, plot = FALSE)
    expect_true(o$nns.regress)
    expect_equal(o$errors, p - tail(x, 27))
    expect_equal(o$results, p + o$bias.shift)
  })
  aud_case("the interval adds the LPM and UPM VaR of actual - results", env = optim_pi_from_returned_model, {
    x <- two_cycles()
    o <- run_optim(x, training.set = 128, seasonal.factor = c(7, 12), obj.fn = mse)
    r <- tail(x, 32) - o$results
    expect_equal(o$lower.pred.int, o$results + LPM.VaR(0.025, 0, r))
    expect_equal(o$upper.pred.int, o$results + UPM.VaR(0.025, 0, r))
  })
  aud_case("the out-of-sample interval uses the validation residuals", env = optim_pi_from_returned_model, {
    x <- two_cycles()
    o <- run_optim(x, h = 10, seasonal.factor = c(7, 12), obj.fn = mse)
    r <- -(o$errors + o$bias.shift)
    expect_equal(o$lower.pred.int, o$results + LPM.VaR(0.025, 0, r))
    expect_equal(o$upper.pred.int, o$results + UPM.VaR(0.025, 0, r))
  })
  aud_case("a nonnegative series takes residuals from the clamped forecasts it returns", {
    # A persistent positive residual centre shifts alternating raw forecasts
    # downward; the lower member must be clamped before the bands are computed.
    testthat::local_mocked_bindings(
      NNS.ARMA = function(variable, h, ...) rep(c(1, 9), length.out = h),
      .package = "NNS"
    )
    v <- rep(1, 80)
    o <- NNS.ARMA.optim(v, seasonal.factor = c(2, 4), training.set = 64, negative.values = FALSE,
                        obj.fn = expression(mean(predicted)), objective = "min",
                        print.trace = FALSE, ncores = 1, plot = FALSE,
                        bias.correction = "split")
    expect_true(any(o$results == 0))
    res <- tail(v, 16) - o$results
    expect_equal(o$upper.pred.int, pmax(0, o$results + UPM.VaR(.025, 0, res)))
    expect_equal(o$lower.pred.int, pmax(0, o$results + LPM.VaR(.025, 0, res)))
    # Residuals of the unclamped forecasts move the bands on this series.
    unclamped <- -(o$errors + o$bias.shift)
    expect_gt(max(abs(c(pmax(0, o$results + UPM.VaR(.025, 0, unclamped)) - o$upper.pred.int,
                        pmax(0, o$results + LPM.VaR(.025, 0, unclamped)) - o$lower.pred.int))), 0.01)
  })
  aud_case("lower prediction bound sits below the forecast on a negative-valued series", env = arma_pi_lower_bound, {
    set.seed(1)
    x <- sin(1:200 / 2) * 5 - 10 + rnorm(200, sd = 0.5)
    r <- arma_pi(x, seasonal.factor = 13)
    expect_true(all(r[[1]] < 0))
    expect_true(all(r[[2]] < r[[1]]))
    expect_true(all(r[[1]] < r[[3]]))
  })
  aud_case("lower prediction bound sits below the forecast on a mean-zero series", env = arma_pi_lower_bound, {
    set.seed(2)
    x <- rnorm(250, 0, 0.01)
    r <- arma_pi(x, seasonal.factor = FALSE)
    expect_true(all(r[[2]] < r[[1]]))
    expect_true(all(r[[1]] < r[[3]]))
  })
  aud_case("pred.int works for a one-point forecast path (#98)", {
    set.seed(123)
    r <- NNS.ARMA(AirPassengers, h = 1, seasonal.factor = 12, pred.int = .95, plot = FALSE)
    expect_equal(nrow(r), 1L)
    expect_lte(r[[2]], r[[1]]); expect_gte(r[[3]], r[[1]])
    set.seed(123)
    r <- NNS.ARMA(AirPassengers, h = 1, seasonal.factor = c(12, 24), method = "lin", pred.int = .95, plot = FALSE)
    expect_lt(r[[2]], r[[3]])   # the per-lag spread widens it
  })
  aud_case("pred.int works for a flat forecast path (#98)", {
    set.seed(123)
    r <- NNS.ARMA(rep(2, 30), h = 12, seasonal.factor = 2, method = "means", pred.int = .95, plot = FALSE)
    expect_equal(unname(unlist(r[, 1])), rep(2, 12))
    expect_equal(unname(unlist(r[, 2])), rep(2, 12))
    expect_equal(unname(unlist(r[, 3])), rep(2, 12))
  })
  aud_case("held-out negative values do not change the clamp", env = arma_train_support, {
    expect_equal(lin_arma(c(20:1, -1), training.set = 20), lin_arma(20:1))
    expect_equal(lin_arma(20:1), c(0, 0, 0))
    expect_true(all(lin_arma(c(20:1, -1)) < 0))
  })
  aud_case("NNS.ARMA.optim clamps validation forecasts from a nonnegative training window", {
    x <- c(seq(58, 2, -2), 1, 0.5, 0.2, 0.1, -0.1)
    run <- function(...) NNS.ARMA.optim(x, training.set = 29, seasonal.factor = c(1, 2), ncores = 1,
                                        print.trace = FALSE, obj.fn = expression(mean((predicted - actual)^2)), ...)
    o <- run()
    expect_true(all(o$errors + tail(x, 5) >= 0))
    expect_true(all(o$results >= 0))
    # The out-of-sample forecast is fitted on the whole series, which is negative.
    expect_true(any(run(h = 3)$results < 0))
  })
  aud_case("ARMA.optim out-of-sample bands sit on the clamped forecasts", {
    # A fixed positive bias in the validation prefix persists in its acceptance
    # half. The same shift must clamp the lower out-of-sample forecast to zero.
    local_mocked_bindings(
      NNS.ARMA = function(variable, h, ...) rep(c(1, 9), length.out = h),
      .package = "NNS"
    )
    x <- rep(1, 80)
    o <- NNS.ARMA.optim(x, h = 6, training.set = 64, seasonal.factor = c(2, 4),
                        obj.fn = expression(mean((predicted - actual)^2)),
                        ncores = 1, print.trace = FALSE, plot = FALSE,
                        bias.correction = "split")
    expect_true(o$bias.validation$accepted)
    expect_true(any(o$results == 0))
    expect_true(all(o$upper.pred.int >= o$results))
    expect_true(all(o$upper.pred.int[o$results == 0] > 0))
    expect_equal(o$results, pmax(0, o$raw.results + o$bias.shift))
  })
})
