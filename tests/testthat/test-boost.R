# NNS.boost: the iris run, the learner threshold, diagnostics plots, response
# scale and labels, and input handling.
#
# Full iris integration run mirroring the documented NNS.boost example. Heavy:
# 100 learner trials (exhaustive: 15 subsets), 100 epochs, final NNS.stack.

# pin: shaped like the documented NNS.boost iris example (epochs, learner.trials, type, threshold)
test_that("BST01: NNS.boost iris run has probability threshold, epochs, and a replicated multivariate final stack", {
  skip_on_cran()
  withr::local_seed(20261108)

  stack_args <- new.env(parent = emptyenv())
  suppressMessages(trace(
    NNS.stack,
    where = asNamespace("NNS"),
    print = FALSE,
    tracer = bquote({
      assign("captured", list(
        method = method,
        folds = folds,
        stack = stack,
        n_predictors = ncol(IVs.train),
        predictor_names = names(IVs.train)
      ), envir = .(stack_args))
    })
  ))
  on.exit(suppressMessages(untrace(NNS.stack, where = asNamespace("NNS"))),
          add = TRUE)

  pdf_file <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)

  set.seed(123)
  msgs <- capture_messages(
    a <- NNS.boost(
      iris[1:140, 1:4],
      iris[1:140, 5],
      IVs.test = iris[141:150, 1:4],
      epochs = 100,
      learner.trials = 100,
      type = "CLASS",
      depth = NULL,
      balance = TRUE,
      feature.importance = TRUE,
      status = TRUE
    )
  )

  grDevices::dev.off()

  # 1-3: the console states the probability and the distinct LPM.VaR cutoff,
  # and never prints the obsolete status string.
  expect_true(any(grepl("Learner threshold probability = 0.80", msgs,
                        fixed = TRUE)))
  cutoff_line <- grep("objective cutoff = ", msgs, fixed = TRUE, value = TRUE)
  expect_length(cutoff_line, 1L)
  cutoff <- as.numeric(sub(".*objective cutoff = ([0-9.]+).*", "\\1",
                           cutoff_line))
  expect_true(is.finite(cutoff))
  expect_false(isTRUE(all.equal(cutoff, 0.80)))
  expect_false(any(grepl("Learner Accuracy Threshold", msgs, fixed = TRUE)))

  # 4-5: both diagnostic panels were drawn (nonempty device output).
  expect_true(file.exists(pdf_file))
  expect_gt(file.size(pdf_file), 0)

  # 6: iris has 15 unique feature subsets < 100 trials (exhaustive), yet
  # epochs still ran.
  expect_true(any(grepl("% of epochs", msgs, fixed = TRUE)))

  # 7-9: results and feature outputs. Classification results keep the
  # historical numeric coding (integer codes, base category 1) and
  # $class.levels recovers the label for each code.
  expect_length(a$results, 10L)
  expect_true(is.numeric(a$results))
  expect_true(all(a$results %in% seq_along(levels(iris$Species))))
  expect_identical(a$class.levels, levels(iris$Species))
  expect_gte(mean(a$results == as.numeric(iris[141:150, 5])), 0.8)
  expect_gte(mean(a$class.levels[a$results] ==
                    as.character(iris[141:150, 5])), 0.8)
  expect_equal(sum(a$feature.weights), 1)
  expect_true(all(names(a$feature.frequency) %in% names(iris)[1:4]))
  expect_true(length(a$feature.frequency) >= 1L)

  # 10: the final call used replicated multivariate predictors through
  # NNS.stack(method = 1, folds = 5, stack = FALSE).
  captured <- get("captured", envir = stack_args)
  expect_equal(as.numeric(captured$method), 1)
  expect_equal(as.numeric(captured$folds), 5)
  expect_false(captured$stack)
  expect_gt(captured$n_predictors, 1L)
  base_names <- sub("[.][0-9]+$", "", captured$predictor_names)
  expect_true(all(base_names %in% names(iris)[1:4]))
  expect_false(any(grepl("xstar", captured$predictor_names, fixed = TRUE)))

  # Replication counts follow the historical scaling rule.
  expected_replication <- pmax(1L, as.integer(round(
    as.numeric(a$feature.frequency) / min(as.numeric(a$feature.frequency))
  )))
  expect_identical(captured$n_predictors, sum(expected_replication))
})

test_that("BST02: NNS.boost passes a supplied threshold to LPM.VaR as a probability (iris run)", {
  skip_on_cran()
  withr::local_seed(20261108)

  # 15 trials still enumerate all 15 feature subsets, and none of the checks
  # depends on the epoch count.
  set.seed(123)
  msgs <- capture_messages(
    b <- NNS.boost(
      iris[1:140, 1:4],
      iris[1:140, 5],
      IVs.test = iris[141:150, 1:4],
      epochs = 10,
      learner.trials = 15,
      type = "CLASS",
      depth = NULL,
      balance = TRUE,
      threshold = 0.65,
      feature.importance = FALSE,
      status = TRUE
    )
  )

  expect_true(any(grepl("Learner threshold probability = 0.65", msgs,
                        fixed = TRUE)))

  # The supplied probability produces a distinct objective cutoff: 0.65 is
  # not applied as a literal accuracy cutoff.
  cutoff_line <- grep("objective cutoff = ", msgs, fixed = TRUE, value = TRUE)
  expect_length(cutoff_line, 1L)
  cutoff <- as.numeric(sub(".*objective cutoff = ([0-9.]+).*", "\\1",
                           cutoff_line))
  expect_true(is.finite(cutoff))
  expect_false(isTRUE(all.equal(cutoff, 0.65)))

  expect_length(b$results, 10L)
})

test_that("BST03: learner threshold: a probability mapped through LPM.VaR, supplied or default, never fivenum or empirical quantiles", {
  withr::local_seed(20261108)
  aud_case("learner threshold maps probabilities through LPM.VaR", {
    scores <- c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

    default_max <- .nns_boost_threshold(NULL, "max", FALSE, scores)
    expect_identical(default_max$probability, 0.80)
    expect_equal(default_max$cutoff, as.numeric(LPM.VaR(0.80, 1, scores)))

    default_min <- .nns_boost_threshold(NULL, "min", FALSE, scores)
    expect_identical(default_min$probability, 0.20)
    expect_equal(default_min$cutoff, as.numeric(LPM.VaR(0.20, 1, scores)))

    extreme_max <- .nns_boost_threshold(NULL, "max", TRUE, scores)
    expect_identical(extreme_max$probability, 1)
    expect_equal(extreme_max$cutoff, as.numeric(LPM.VaR(1, 1, scores)))

    extreme_min <- .nns_boost_threshold(NULL, "min", TRUE, scores)
    expect_identical(extreme_min$probability, 0)
    expect_equal(extreme_min$cutoff, as.numeric(LPM.VaR(0, 1, scores)))
  })
  aud_case("a supplied threshold is used as the LPM.VaR probability", {
    scores <- c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00)

    supplied <- .nns_boost_threshold(0.65, "max", FALSE, scores)
    expect_identical(supplied$probability, 0.65)
    expect_equal(supplied$cutoff, as.numeric(LPM.VaR(0.65, 1, scores)))

    # The probability is not treated as a literal objective-score cutoff.
    expect_false(isTRUE(all.equal(supplied$cutoff, 0.65)))

    # extreme overrides a supplied probability.
    overridden <- .nns_boost_threshold(0.65, "max", TRUE, scores)
    expect_identical(overridden$probability, 1)
  })
  aud_case("threshold selection does not use fivenum or empirical quantiles", {
    boost_source <- paste(deparse(body(NNS.boost)), collapse = "\n")
    helper_source <- paste(deparse(body(.nns_boost_threshold)), collapse = "\n")

    for (src in c(boost_source, helper_source)) {
      expect_false(grepl("fivenum", src, fixed = TRUE))
      expect_false(grepl("quantile(", src, fixed = TRUE))
      expect_false(grepl("stats::quantile", src, fixed = TRUE))
    }

    # The cutoff comes from LPM.VaR, and epochs reuse the same cutoff variable.
    expect_true(grepl("LPM.VaR", helper_source, fixed = TRUE))
    expect_true(grepl("learner.threshold", boost_source, fixed = TRUE))
    expect_false(grepl("Learner Accuracy Threshold", boost_source, fixed = TRUE))
  })
})

test_that("BST04: diagnostics: two panels, the cutoff line and label, drawn only with feature.importance, on a PDF device", {
  withr::local_seed(20261108)
  boost_plots <- fx_boost_plots()
  aud_case("diagnostics draws both panels in one two-panel layout", env = boost_plots, {
    fx <- .boost_plot_fixture()
    seen <- new.env(parent = emptyenv())

    local_mocked_bindings(
      .nns_boost_plot_learner_distribution = function(learner.results,
                                                      threshold.probability,
                                                      learner.threshold) {
        seen$dist_layout <- graphics::par("mfrow")
        seen$dist_threshold <- learner.threshold
        invisible(NULL)
      },
      .nns_boost_plot_feature_frequency = function(feature.frequency) {
        seen$freq_layout <- graphics::par("mfrow")
        invisible(NULL)
      },
      .package = "NNS"
    )

    pdf_file <- withr::local_tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_file)
    on.exit(grDevices::dev.off(), add = TRUE)

    .nns_boost_plot_diagnostics(
      learner.results = fx$learner.results,
      threshold.probability = fx$threshold.probability,
      learner.threshold = fx$learner.threshold,
      feature.frequency = fx$feature.frequency
    )

    # Both helpers were invoked, under the same undisturbed two-panel layout.
    expect_identical(seen$dist_layout, c(2L, 1L))
    expect_identical(seen$freq_layout, c(2L, 1L))

    # The threshold line receives learner.threshold, not the probability.
    expect_identical(seen$dist_threshold, fx$learner.threshold)

    # The original graphics parameters are restored after both panels.
    expect_identical(graphics::par("mfrow"), c(1L, 1L))
  })
  aud_case("the learner-distribution panel draws the cutoff line and label", {
    src <- paste(deparse(body(.nns_boost_plot_learner_distribution)),
                 collapse = "\n")
    expect_true(grepl("abline", src, fixed = TRUE))
    expect_true(grepl("v = learner.threshold", src, fixed = TRUE))
    expect_true(grepl("LPM.VaR(p = %.2f)", src, fixed = TRUE))
    expect_true(grepl("hist", src, fixed = TRUE))

    freq_src <- paste(deparse(body(.nns_boost_plot_feature_frequency)),
                      collapse = "\n")
    expect_true(grepl("barplot", freq_src, fixed = TRUE))
    expect_true(grepl("horiz = TRUE", freq_src, fixed = TRUE))
  })
  aud_case("NNS.boost invokes both plot helpers only when feature.importance = TRUE", {
    calls <- new.env(parent = emptyenv())
    calls$dist <- 0L
    calls$freq <- 0L

    local_mocked_bindings(
      .nns_boost_plot_learner_distribution = function(...) {
        calls$dist <- calls$dist + 1L
        invisible(NULL)
      },
      .nns_boost_plot_feature_frequency = function(...) {
        calls$freq <- calls$freq + 1L
        invisible(NULL)
      },
      .package = "NNS"
    )

    pdf_file <- withr::local_tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_file)
    on.exit(grDevices::dev.off(), add = TRUE)

    run_boost <- function(feature.importance) {
      NNS.boost(
        iris[1:60, 1:4],
        iris[1:60, 5],
        IVs.test = iris[61:70, 1:4],
        epochs = 3,
        learner.trials = 5,
        type = "CLASS",
        feature.importance = feature.importance,
        status = FALSE,
        seed = 123L
      )
    }

    invisible(run_boost(feature.importance = TRUE))
    expect_identical(calls$dist, 1L)
    expect_identical(calls$freq, 1L)

    invisible(run_boost(feature.importance = FALSE))
    expect_identical(calls$dist, 1L)
    expect_identical(calls$freq, 1L)
  })
  aud_case("NNS.boost completes its diagnostics on a PDF device without errors", {
    pdf_file <- withr::local_tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_file)

    expect_no_error(
      invisible(NNS.boost(
        iris[1:60, 1:4],
        iris[1:60, 5],
        IVs.test = iris[61:70, 1:4],
        epochs = 3,
        learner.trials = 5,
        type = "CLASS",
        feature.importance = TRUE,
        status = FALSE,
        seed = 123L
      ))
    )

    grDevices::dev.off()
    expect_true(file.exists(pdf_file))
    expect_gt(file.size(pdf_file), 0)
  })
})

test_that("BST05: responses keep their scale and labels: 0/1, -5/5, two-valued objectives, numeric binary and permuted factor levels", {
  withr::local_seed(20261108)
  stack_boost <- fx_stack_boost()
  aud_case("NNS.boost returns 0/1 response values for a 0/1 numeric response", env = stack_boost, {
    f <- boost_binary(0, 1)
    expect_equal(f$results, c(1, 1, 1, 0, 0, 0))
    expect_equal(unname(as.matrix(f$pred.int)), cbind(c(1, 1, 1, 0, 0, 0), c(1, 1, 1, 0, 0, 0)))
  })
  aud_case("NNS.boost returns -5/5 response values for a -5/5 numeric response", env = stack_boost, {
    # Near-zero epoch scores can trigger the best-trial fallback by roundoff.
    # Accept only that documented warning, without requiring it on every platform.
    fallback_warning <- paste0(
      "No feature combination re-passed the objective cutoff during epochs; ",
      "using the best learner-trial combination. ",
      "Consider a lower [threshold] probability."
    )
    f <- withCallingHandlers(
      boost_binary(-5, 5),
      warning = function(w) {
        expect_identical(conditionMessage(w), fallback_warning)
        invokeRestart("muffleWarning")
      }
    )
    expect_equal(f$results, c(5, 5, 5, -5, -5, -5))
    expect_equal(unname(as.matrix(f$pred.int)), cbind(c(5, 5, 5, -5, -5, -5), c(5, 5, 5, -5, -5, -5)))
  })
  aud_case("boost's objective sees response values for a two-valued numeric response", {
    set.seed(1)
    x <- matrix(rnorm(240), 80, 3)
    y <- ifelse(x[, 1] + rnorm(80, sd = .5) > 0, 5, -5)
    seen <- numeric(0)
    record <- function(a) seen <<- c(seen, a)
    loss <- as.expression(bquote({
      .(record)(actual)
      sum((predicted - actual)^2)
    }))
    NNS.boost(x, y, x[1:10, ], obj.fn = loss, learner.trials = 3, epochs = 0,
              threshold = 1, CV.size = .25, seed = 1, status = FALSE,
              feature.importance = FALSE)
    expect_setequal(unique(seen), c(-5, 5))
    seen <- numeric(0)
    NNS.boost(x, ifelse(y > 0, .2, -.1), x[1:10, ], obj.fn = loss, learner.trials = 3,
              epochs = 0, threshold = 1, CV.size = .25, seed = 1, status = FALSE,
              feature.importance = FALSE)
    expect_true(all(seen %in% c(-.1, .2)))
  })
  aud_case("boost does not decode numeric binary stack results twice", {
    x <- iris[1:100, 1:2]
    query <- x[c(1, 100), ]
    for (values in list(c(-5, 8), c(-.1, .2))) {
      y <- rep(values, each = 50)
      seen <- numeric()
      record <- function(actual) seen <<- c(seen, actual)
      loss <- as.expression(bquote({ .(record)(actual); sum((predicted - actual)^2) }))
      fit <- NNS.boost(x, y, query, obj.fn = loss, learner.trials = 3,
        epochs = 0, threshold = 1, seed = 7, folds = 2, status = FALSE,
        feature.importance = FALSE)
      expect_setequal(unique(seen), values)
      expect_equal(fit$results, values)
    }
  })
  aud_case("boost keeps factor label decisions when class levels are permuted", {
    run <- function(y) NNS.boost(iris[, 1:3], y, iris[c(1, 55, 110), 1:3],
      learner.trials = 3, epochs = 0, balance = TRUE, threshold = 1, type = "CLASS",
      seed = 13, folds = 2, status = FALSE, feature.importance = FALSE)
    a <- run(iris$Species)
    b <- run(factor(iris$Species, levels = rev(levels(iris$Species))))
    expect_equal(a$class.levels[a$results], b$class.levels[b$results])
    expect_equal(a$feature.weights, b$feature.weights)
  })
})

test_that("BST06: inputs: duplicate predictors, test-column names, NULL distance, positional objective, single-member class", {
  withr::local_seed(20261108)
  stack_boost <- fx_stack_boost()
  aud_case("NNS.boost accepts duplicate predictor columns like NNS.reg (cbind(x, x))", {
    set.seed(123)
    x <- rnorm(60)
    y <- as.numeric(x > 0)
    train <- cbind(x, x)

    # Same cbind(x, x) dimension trick as NNS.reg; NNS.boost must not error with
    # "[IVs.train] predictor names must be unique." any more, and duplicate-named
    # test columns must align rather than error.
    expect_error(
      NNS.boost(IVs.train = train, DV.train = y, IVs.test = train[1:5, ],
                learner.trials = 5, epochs = 5, CV.size = 0.25,
                status = FALSE),
      NA
    )
  })
  aud_case("NNS.boost applies the same test-column rule", {
    set.seed(9)
    x <- data.frame(a = rnorm(80), b = rnorm(80), c = rnorm(80))
    y <- as.numeric(x$a + x$b > 0)

    foreign <- setNames(x[1:6, ], c("t.a", "t.b", "t.c"))
    expect_error(
      NNS.boost(IVs.train = x, DV.train = y, IVs.test = foreign,
                epochs = 8, learner.trials = 5, folds = 1, status = FALSE),
      NA
    )

    partial <- setNames(x[1:6, ], c("a", "t.b", "t.c"))
    expect_error(
      NNS.boost(IVs.train = x, DV.train = y, IVs.test = partial,
                epochs = 8, learner.trials = 5, folds = 1, status = FALSE),
      "exactly match"
    )

    by.name <- NNS.boost(IVs.train = x, DV.train = y,
                         IVs.test = x[1:6, c("c", "a", "b")],
                         epochs = 8, learner.trials = 5, folds = 1,
                         status = FALSE, seed = 7)$results
    in.order <- NNS.boost(IVs.train = x, DV.train = y,
                          IVs.test = x[1:6, c("a", "b", "c")],
                          epochs = 8, learner.trials = 5, folds = 1,
                          status = FALSE, seed = 7)$results
    expect_equal(by.name, in.order, tolerance = 1e-12)
  })
  aud_case("NNS.boost NULL distance matches explicit NNS and propagates L2", {
    x <- data.frame(a = c(0, 1, 0.2, 0.9, 0.4, 0.8, 0.3, 0.7),
                    b = c(0, 0.2, 1, 0.9, 0.7, 0.1, 0.6, 0.4))
    y <- c(0, 2, 4, 8, 5, 3, 6, 7)
    args <- list(IVs.train = x, DV.train = y, IVs.test = x[1:2, ], learner.trials = 2,
                 epochs = 0, status = FALSE, seed = 42)
    a <- do.call(NNS.boost, c(args, list(dist = NULL)))
    b <- do.call(NNS.boost, c(args, list(dist = "NNS")))
    expect_equal(a$results, b$results, tolerance = 1e-12)
    expect_no_error(do.call(NNS.boost, c(args, list(dist = "L2"))))
  })
  aud_case("NNS.boost preserves old positional objective argument order", {
    x <- data.frame(a = c(0, 1, 0.2, 0.9, 0.4, 0.8, 0.3, 0.7),
                    b = c(0, 0.2, 1, 0.9, 0.7, 0.1, 0.6, 0.4))
    y <- c(0, 2, 4, 8, 5, 3, 6, 7)
    test <- x[1:2, ]
    obj <- expression(sum(abs(predicted - actual)))

    positional <- NNS.boost(x, y, test, NULL, NULL, 2, 0, 0.25, FALSE,
                            NULL, NULL, obj, "min", FALSE, FALSE, FALSE,
                            NULL, FALSE, 42)
    named <- NNS.boost(IVs.train = x, DV.train = y, IVs.test = test,
                       learner.trials = 2, epochs = 0, CV.size = 0.25,
                       balance = FALSE, ts.test = NULL, threshold = NULL,
                       obj.fn = obj, objective = "min", extreme = FALSE,
                       features.only = FALSE, feature.importance = FALSE,
                       pred.int = NULL, status = FALSE, seed = 42)

    expect_equal(positional$results, named$results, tolerance = 1e-12)
  })
  aud_case("NNS.boost balancing keeps a single-member class", env = stack_boost, {
    set.seed(9)
    x <- data.frame(a = rnorm(60), b = rnorm(60))
    frame <- capture_frame("NNS.boost", NNS::NNS.boost(
      x, factor(rep(1:3, c(29, 29, 2))), IVs.test = x[1:3, ], type = "CLASS",
      balance = TRUE, epochs = 0, learner.trials = 3, feature.importance = FALSE,
      status = FALSE, seed = 1
    ))
    set.seed(1)
    z <- frame$.balance_training(data.frame(i = seq_along(balance_y)), balance_y)
    expect_identical(as.integer(table(z$y)), c(4L, 4L, 4L))
    expect_identical(balance_y[z$x$i], z$y)
  })
})
