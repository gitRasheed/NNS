# NNS.VAR: forecast validation, forecast rows and interpolation.
#
# NNS.VAR's interpolated_and_extrapolated values for mixed-frequency panels.

test_that("VAR01: forecast validation sees only the past: covariates, mixed-frequency prefix, relevance coefficients, two-level targets", {
  withr::local_seed(20261108)
  aud_case("VAR validation covariates cannot see observations after their origin", {
    captured <- list()
    testthat::local_mocked_bindings(
      NNS.ARMA.optim = function(variable, h, ...) {
        list(results = tail(variable, 1) + seq_len(h) / 10)
      },
      .nns_var_relevant = function(IVs.train, DV.train, ..., .validation = NULL) {
        colnames(IVs.train)
      },
      NNS.stack = function(IVs.train, DV.train, IVs.test, ..., .validation = NULL) {
        captured[[length(captured) + 1L]] <<- .validation
        list(stack = rep(mean(DV.train), nrow(IVs.test)))
      },
      .package = "NNS"
    )
    t <- seq_len(40)
    x <- cbind(a = sin(t / 5), b = cos(t / 7))
    NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
    first <- captured
    x[33:40, ] <- x[33:40, ] + 100
    captured <- list()
    NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
    expect_length(first, 2)
    for (i in 1:2) {
      expect_equal(captured[[i]]$IVs.train, first[[i]]$IVs.train)
      expect_equal(captured[[i]]$IVs.test, first[[i]]$IVs.test)
      expect_equal(captured[[i]]$DV.train, first[[i]]$DV.train)
      expect_equal(captured[[i]]$DV.test - first[[i]]$DV.test, rep(100, 8))
    }
    # Later validation lags must also use forecasts, not realized heldout values.
    expect_equal(first[[1]]$IVs.test$b_tau_0, x[32, "b"] + (1:8) / 10)
    expect_equal(tail(first[[1]]$IVs.test$b_tau_2, 6), x[32, "b"] + (1:6) / 10)
  })
  aud_case("VAR completes a mixed-frequency prefix before seeing later releases", {
    captured <- list()
    testthat::local_mocked_bindings(
      NNS.ARMA.optim = function(variable, h, ...) {
        list(results = tail(variable, 1) + seq_len(h) / 10)
      },
      # Missing releases are completed from the observed prefix alone.
      .nns_var_index_fit = function(x, y, x.new, ...) rep(mean(y), length(x.new)),
      NNS.reg = function(x, y, point.est, ...) {
        list(Point.est = rep(mean(y), length(point.est)))
      },
      .nns_var_relevant = function(IVs.train, DV.train, ..., .validation = NULL) {
        colnames(IVs.train)
      },
      NNS.stack = function(IVs.train, DV.train, IVs.test, ..., .validation = NULL) {
        if (!is.null(.validation)) captured[[length(captured) + 1L]] <<- .validation
        list(stack = rep(mean(DV.train), nrow(IVs.test)))
      },
      .package = "NNS"
    )
    t <- seq_len(40)
    x <- cbind(a = 20 + sin(t / 5), b = 20 + cos(t / 7))
    x[32, "b"] <- NA_real_  # trailing missing release at the validation origin
    x[34, "a"] <- NA_real_
    x[c(35, 38), "b"] <- NA_real_
    out <- NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
    first <- captured
    expect_equal(unname(out$validation$observed.targets), c(7L, 6L))
    expect_length(first[[1]]$DV.test, 7L)
    expect_length(first[[2]]$DV.test, 6L)
    expect_equal(tail(first[[2]]$DV.train, 1), mean(x[1:31, "b"]))
    x[33:40, ] <- x[33:40, ] + 100
    captured <- list()
    NNS.VAR(x, h = 4, tau = 2, ncores = 1, status = FALSE)
    for (i in 1:2) {
      expect_equal(captured[[i]]$IVs.train, first[[i]]$IVs.train)
      expect_equal(captured[[i]]$DV.train, first[[i]]$DV.train)
      expect_equal(captured[[i]]$IVs.test, first[[i]]$IVs.test)
    }
  })
  aud_case("VAR relevance coefficients never fit heldout rows in pipeline validation", {
    trained_rows <- integer()
    testthat::local_mocked_bindings(
      NNS.reg = function(x, y, point.est = NULL, ...) {
        trained_rows <<- c(trained_rows, length(y))
        list(equation = data.frame(Coefficient = rep(1, ncol(x))),
             Point.est = rep(mean(y), NROW(point.est)))
      },
      .package = "NNS"
    )
    t <- seq_len(40)
    x <- data.frame(a = sin(t), b = cos(t), c = t)
    y <- sin(t / 5)
    fold <- list(IVs.train = x[1:32, ], DV.train = y[1:32],
                 IVs.test = x[33:40, ] * 0, DV.test = y[33:40])
    NNS:::.nns_var_relevant(x, y, 8,
      expression(mean((predicted - actual)^2)), "min", "cor", .validation = fold)
    expect_true(length(trained_rows) > 0)
    expect_true(all(trained_rows == 32L))
  })
  aud_case("numeric two-level VAR targets remain forecast regressions", {
    captured <- new.env()
    captured$frames <- list()
    suppressMessages(trace(".nns_var_relevant", where = asNamespace("NNS"),
      print = FALSE, exit = bquote(assign("frames",
        c(get("frames", envir = .(captured)), list(environment())),
        envir = .(captured)))))
    withr::defer(suppressMessages(untrace(".nns_var_relevant",
                                        where = asNamespace("NNS"))))
    local_mocked_bindings(
      NNS.ARMA.optim = function(variable, h, ...) {
        list(results = rep(mean(variable), h))
      },
      .package = "NNS"
    )
    prediction.reference <- list()
    for (levels in list(c(0, 1), c(-5, 5))) for (method in c("cor", "NNS.dep")) {
      captured$frames <- list()
      t <- seq_len(40)
      x <- cbind(a = rep(levels, 20), b = sin(t / 5))
      out <- NNS.VAR(x, h = 2, tau = 1, ncores = 1, status = FALSE,
                     dim.red.method = method,
                     obj.fn = expression(mean((predicted - actual)^2)))
      expect_identical(dim(out$ensemble), c(2L, 2L))
      expect_true(all(is.finite(as.matrix(out$ensemble))))
      expect_equal(unname(out$validation$observed.targets), c(8L, 8L))
      # Relevance must extract real regression coefficients, not NULL from a
      # list of nominal class equations and then silently retain every column.
      fr <- captured$frames[[1L]]
      expect_length(fr$coef, 3L)
      expect_true(all(is.finite(fr$coef)))
      if (method == "NNS.dep") {
        expect_equal(fr$coef, c(0, 1, 0), tolerance = 1e-12)
        expect_equal(na.omit(out$relevant_variables$a), "a_tau_1",
                     ignore_attr = TRUE)
        # The first lag is observed at the origin; the seven later lags are
        # prefix-mean forecasts. Their midpoint errors yield 7/8 * range^2/4.
        expect_equal(fr$nns.ord, diff(levels)^2 * 7 / 32, tolerance = 1e-12)
        expect_equal(out$multivariate$a, c(levels[1], mean(levels)),
                     tolerance = 1e-10)
      } else {
        expected <- apply(fr$CV.IVs.train, 2, stats::cor, y = fr$CV.DV.train,
                           method = "spearman")
        expect_equal(fr$coef, unname(expected), tolerance = 1e-12)
        expect_equal(na.omit(out$relevant_variables$a), c("b_tau_0", "a_tau_1"),
                     ignore_attr = TRUE)
        expect_gt(fr$nns.ord / diff(levels)^2, .21)
        expect_lt(fr$nns.ord / diff(levels)^2, .25)
        expect_true(all(out$multivariate$a > min(levels) &
                          out$multivariate$a < max(levels)))
        expect_lt(out$multivariate$a[1], mean(levels))
      }
      if (identical(levels, c(0, 1))) prediction.reference[[method]] <- out$multivariate$a else
        expect_equal(out$multivariate$a,
                     levels[1] + diff(levels) * prediction.reference[[method]],
                     tolerance = 1e-10)
    }
  })
})

test_that("VAR02: forecast rows: multi-step, per-variable lags, long list lags, a single series, response to the univariate forecasts", {
  withr::local_seed(20261108)
  var_forecast_rows <- fx_var_forecast_rows()
  aud_case("multi-step NNS.VAR evaluates the multivariate stage on the forecast rows", env = var_forecast_rows, {
    run <- var_with_test_rows(var_fixture(), h = 4, tau = 1)
    expect_forecast_rows(run, h = 4, n_stacks = 2)
  })
  aud_case("one-step NNS.VAR with per-variable lags uses the forecast row", env = var_forecast_rows, {
    run <- var_with_test_rows(var_fixture(), h = 1, tau = c(1, 3))
    expect_forecast_rows(run, h = 1, n_stacks = 2)
  })
  aud_case("list lags longer than h mix observed and forecast values in the forecast rows", env = var_forecast_rows, {
    run <- var_with_test_rows(var_fixture(), h = 2, tau = list(c(1, 3), 2))
    expect_forecast_rows(run, h = 2, n_stacks = 2)
  })
  aud_case("a single series uses its own extended lags as forecast rows", env = var_forecast_rows, {
    run <- var_with_test_rows(var_fixture()[, "x", drop = FALSE], h = 3, tau = 2)
    expect_forecast_rows(run, h = 3, n_stacks = 1)
  })
  aud_case("the multivariate forecast responds to the univariate forecasts", env = var_forecast_rows, {
    shift <- 0
    local_mocked_bindings(
      NNS.ARMA.optim = function(variable, h, ...) {
        list(results = tail(variable, 1) + shift + seq_len(h) / 10)
      },
      .package = "NNS"
    )
    x <- var_fixture()
    low <- NNS.VAR(x, h = 3, tau = 1, ncores = 1, status = FALSE)
    shift <- 5
    high <- NNS.VAR(x, h = 3, tau = 1, ncores = 1, status = FALSE)
    expect_equal(as.matrix(high$univariate) - as.matrix(low$univariate),
                 matrix(5, 3, 2), ignore_attr = TRUE)
    expect_gt(max(abs(as.matrix(high$multivariate) - as.matrix(low$multivariate))), 0.01)
  })
})

test_that("VAR03: interpolation: the index fit, trailing and interior NAs, one or two observations", {
  withr::local_seed(20261108)
  aud_case("the index fit is NNS 13.0's one-neighbour fit on the duplicated index", {
    # Regression points 1, 7 and 13 give three cells. The newest observation
    # (13) is a cell of its own; the others are summarised by the gravity (here
    # the median) of their index and values: (2.5, 11) and (8.5, 15.5).
    x <- c(1, 4, 7, 10, 13)
    y <- c(10, 12, 15, 16, 20)
    fit <- .nns_var_index_fit(x, y, c(0, 2, 5, 8, 11, 14, 16), points = c(1, 7, 13))
    # Inside the data: the value of the nearest regression point.
    expect_equal(fit[2:5], c(11, 11, 15.5, 20))
    # After the data: from 20, the 3:2:1 average of the slopes to the centre of
    # the regression points (8.5) and to 10.75 and 11.875, whose nearest points
    # are 8.5 (the first of two tied) and 13: (3 * 1 + 2 * 2 + 0) / 6 per period.
    expect_equal(fit[6:7], 20 + c(1, 3) * 7 / 6, tolerance = 1e-12)
    # Before the data: from the first cell's 11, the same rule gives -0.3 per
    # period.
    expect_equal(fit[1], 10.7, tolerance = 1e-12)
  })
  aud_case("with trailing NAs, every NA comes from the index fit", {
    # A quarterly series in a monthly panel whose last two months are not yet
    # observed.
    set.seed(11)
    n <- 120
    t <- seq_len(n)
    q <- 100 + 0.5 * t + 2 * sin(2 * pi * t / 12) + rnorm(n, sd = 0.3)
    m <- 50 + cumsum(rnorm(n, 0.1, 0.5))
    keep <- seq(1, n - 2, by = 3)
    x <- rep(NA_real_, n)
    x[keep] <- q[keep]
    out <- NNS.VAR(data.frame(q = x, m = m), h = 0, tau = 3, status = FALSE, ncores = 1)

    miss <- setdiff(t, keep)
    expect_identical(out$q[keep], x[keep])
    expect_identical(out$m, m)
    expect_equal(out$q[miss], .nns_var_index_fit(keep, x[keep], miss), tolerance = 1e-12)
    # The newest observation keeps its own regression point, so the trailing
    # values lie on a line through it.
    last <- max(keep)
    expect_equal(2 * out$q[last + 1] - out$q[last + 2], x[last], tolerance = 1e-10)
  })
  aud_case("without trailing NAs, interior NAs are interpolated linearly", {
    set.seed(13)
    m <- 50 + cumsum(rnorm(30))
    z <- 10 + 0.3 * seq_len(30) + rnorm(30)
    z[c(5, 6, 17)] <- NA
    out <- NNS.VAR(data.frame(z = z, m = m), h = 0, status = FALSE, ncores = 1)
    obs <- which(!is.na(z))
    expect_equal(out$z[c(5, 6, 17)], approx(obs, z[obs], c(5, 6, 17))$y, tolerance = 1e-12)
    expect_identical(out$z[obs], z[obs])
  })
  aud_case("a series with one or two observations is extended from them", {
    set.seed(12)
    m <- 50 + cumsum(rnorm(40))
    out <- NNS.VAR(data.frame(z = c(3, 4, rep(NA, 38)), m = m), h = 0, status = FALSE, ncores = 1)
    expect_equal(out$z, c(3, 4, 5:42))
    out <- NNS.VAR(data.frame(z = c(3, rep(NA, 39)), m = m), h = 0, status = FALSE, ncores = 1)
    expect_identical(out$z, rep(3, 40))
  })
})
