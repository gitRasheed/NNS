# NNS.reg and NNS.part fits: reducers, test columns, R2, regression or
# classification, the argument contract and predictions.

# pin: NNS-Python reference values
test_that("REG01: fits: not pulled below noise, mode reducer on any scale, integer-valued data, mode reducers agree, NNS.part map and infinities", {
  withr::local_seed(20261108)
  partition_map_ref <- fx_partition_map_ref()
  aud_case("default NNS.reg fits are not pulled below the data on pure noise", {
    # The default reducer is gravity: on noise its fits sat 0.047 below the
    # noise.reduction = "mean" fits.
    d <- vapply(1:40, function(s) {
      set.seed(s); x <- rnorm(200); y <- rnorm(200)
      mean(NNS.reg(x, y, plot = FALSE, residual.plot = FALSE)$Fitted.xy$y.hat -
             NNS.reg(x, y, noise.reduction = "mean", plot = FALSE, residual.plot = FALSE)$Fitted.xy$y.hat)
    }, numeric(1))
    expect_lt(abs(mean(d)), 0.02)
  })
  aud_case("the regression mode reducer works on any scale", {
    # It used the discrete mode, which rounds to integers: on x in [0, 1] the
    # fit had 3 regression points and R2 0.72.
    set.seed(1); x <- runif(300); y <- sin(6 * x) + rnorm(300, sd = .1)
    f <- NNS.reg(x, y, noise.reduction = "mode", plot = FALSE, residual.plot = FALSE)
    expect_gt(nrow(f$regression.points), 10)
    expect_gt(f$R2, 0.9)
    g <- NNS.reg(100 * x, y, noise.reduction = "mode", plot = FALSE, residual.plot = FALSE)
    expect_equal(g$Fitted.xy$y.hat, f$Fitted.xy$y.hat, tolerance = 1e-8)
    p <- NNS.part(x, y, order = 2, noise.reduction = "mode")$regression.points
    expect_false(all(p$x == round(p$x)))
  })
  aud_case("NNS.part keeps the x centroid of integer-valued x", {
    p <- NNS.part(1:20, 2 * (1:20), type = "XONLY", order = 2,
                  noise.reduction = "mean", obs.req = 0)
    expect_equal(p$regression.points$x, c(5.5, 15.5))
    expect_equal(p$regression.points$y, c(11, 31))
  })
  aud_case("NNS.reg fits an integer-valued line exactly", {
    x <- 1:10
    fit <- NNS.reg(x, 3 * x + 1, plot = FALSE, residual.plot = FALSE, ncores = 1)
    expect_equal(fit$R2, 1)
    expect_equal(fit$Fitted.xy$y.hat, 3 * x + 1)
  })
  aud_case("native and pure-R multivariate mode reducers agree", {
    set.seed(4)
    X <- matrix(round(rnorm(120), 1), ncol = 2)
    y <- round(X[, 1] + X[, 2]^2 + rnorm(60, sd = 0.3), 1)
    fit <- function(native) {
      op <- options(NNS.native.mreg = native)
      on.exit(options(op))
      NNS.reg(X, y, noise.reduction = "mode", point.est = X[1:4, ],
              plot = FALSE, residual.plot = FALSE, ncores = 1)
    }
    a <- fit(TRUE)
    b <- fit(FALSE)
    expect_equal(a$RPM, b$RPM)
    expect_equal(a$Fitted.xy$y.hat, b$Fitted.xy$y.hat)
    expect_equal(a$Point.est, b$Point.est)
  })
  aud_case("NNS.part", env = partition_map_ref, {
    expect_equal(R1$order, T_ORDER, tolerance=1e-5)
    expect_equal(as.data.frame(R1$dt), T_DT, tolerance=1e-5)
    expect_equal(as.data.frame(R1$regression.points), T_regression_points, tolerance=1e-5)
  })
  aud_case("native partition complete-case handling keeps infinities", {
    part <- NNS.part(
      x = c(Inf, Inf, NA_real_, NaN),
      y = c(1, 2, 3, 4),
      order = 1,
      obs.req = 0,
      noise.reduction = "median"
    )

    expect_true(any(is.infinite(part$regression.points$x)))
    expect_false(any(is.na(part$regression.points$x)))
  })
})

test_that("REG02: test columns: positional when no names match, never silently reordered, aligned by name, partial overlap and schema errors", {
  withr::local_seed(20261108)
  aud_case("test columns sharing no names with training are matched positionally", {
    # cbind() names its columns after the supplied expressions, so a test set
    # built as cbind(test.x_1, test.x_2) carries names that describe nothing
    # about the training predictors x_1 / x_2. Those names cannot express an
    # alignment, so the columns are taken in the order supplied.
    set.seed(123)
    x_1 <- rnorm(120); x_2 <- rnorm(120)
    y <- 10 * x_1 + 10 * x_2 + rnorm(120)

    set.seed(321)
    test.x_1 <- rnorm(20); test.x_2 <- rnorm(20)

    named <- NNS.stack(cbind(x_1, x_2), y, IVs.test = cbind(test.x_1, test.x_2),
                       method = 1, folds = 2, ncores = 1, status = FALSE)$stack
    bare <- NNS.stack(cbind(x_1, x_2), y,
                      IVs.test = unname(cbind(test.x_1, test.x_2)),
                      method = 1, folds = 2, ncores = 1, status = FALSE)$stack

    expect_length(named, 20L)
    expect_equal(named, bare, tolerance = 1e-12)
  })
  aud_case("foreign test names do not silently reorder columns", {
    set.seed(4)
    x <- data.frame(a = rnorm(80), b = rnorm(80))
    y <- x$a - 2 * x$b

    ordered <- x[1:5, c("a", "b")]
    reversed <- x[1:5, c("b", "a")]

    foreign <- setNames(ordered, c("zz", "yy"))
    foreign.reversed <- setNames(reversed, c("zz", "yy"))

    base <- NNS.reg(x, y, point.est = ordered, plot = FALSE,
                    residual.plot = FALSE, ncores = 1)$Point.est

    # Positional, so the same values in the same order agree ...
    expect_equal(
      NNS.reg(x, y, point.est = foreign, plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      base, tolerance = 1e-12
    )
    # ... and swapping the columns changes the answer, since the names are not
    # consulted at all.
    expect_false(isTRUE(all.equal(
      NNS.reg(x, y, point.est = foreign.reversed, plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      base
    )))
  })
  aud_case("matching names are still aligned by name, not by position", {
    set.seed(5)
    x <- data.frame(a = rnorm(80), b = rnorm(80))
    y <- x$a - 2 * x$b

    expect_equal(
      NNS.reg(x, y, point.est = x[1:5, c("a", "b")], plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      NNS.reg(x, y, point.est = x[1:5, c("b", "a")], plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      tolerance = 1e-12
    )
  })
  aud_case("partially overlapping names remain an error", {
    # A test set that names some training predictors and not others is
    # ambiguous, and in practice a wrong or misspelled column.
    set.seed(6)
    x <- data.frame(a = rnorm(60), b = rnorm(60))
    y <- x$a + x$b

    bad <- setNames(x[1:5, ], c("a", "c"))

    expect_error(
      NNS.reg(x, y, point.est = bad, plot = FALSE, residual.plot = FALSE),
      "exactly match"
    )
    expect_error(
      NNS.stack(x, y, IVs.test = bad, method = 1, folds = 2, ncores = 1,
                status = FALSE),
      "exactly match"
    )
  })
  aud_case("named test vectors follow the same rule", {
    set.seed(7)
    x <- data.frame(a = rnorm(60), b = rnorm(60))
    y <- x$a - 2 * x$b

    base <- NNS.reg(x, y, point.est = c(a = 0.5, b = -0.5), plot = FALSE,
                    residual.plot = FALSE, ncores = 1)$Point.est

    # Permutation of the training names is aligned by name.
    expect_equal(
      NNS.reg(x, y, point.est = c(b = -0.5, a = 0.5), plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      base, tolerance = 1e-12
    )
    # Names shared with nothing in training fall back to position.
    expect_equal(
      NNS.reg(x, y, point.est = c(q = 0.5, z = -0.5), plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      base, tolerance = 1e-12
    )
    expect_error(
      NNS.reg(x, y, point.est = c(a = 0.5, z = -0.5), plot = FALSE,
              residual.plot = FALSE),
      "exactly match"
    )
  })
  aud_case("a single predictor accepts any test column name", {
    set.seed(8)
    x <- rnorm(60); y <- x^2 + rnorm(60)

    expect_equal(
      NNS.reg(x, y, point.est = cbind(test.x = x[1:4]), plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      NNS.reg(x, y, point.est = x[1:4], plot = FALSE,
              residual.plot = FALSE, ncores = 1)$Point.est,
      tolerance = 1e-12
    )
  })
  aud_case("named prediction columns are reordered and schema checked", {
    set.seed(1)
    x <- data.frame(a = rnorm(40), b = rnorm(40))
    y <- x$a - 2 * x$b
    p1 <- x[1:5, c("a", "b")]
    p2 <- x[1:5, c("b", "a")]

    f1 <- NNS.reg(x, y, point.est = p1, plot = FALSE,
                  residual.plot = FALSE, ncores = 1)
    f2 <- NNS.reg(x, y, point.est = p2, plot = FALSE,
                  residual.plot = FALSE, ncores = 1)
    expect_equal(f1$Point.est, f2$Point.est, tolerance = 1e-12)

    bad <- p1
    names(bad) <- c("a", "c")
    expect_error(
      NNS.reg(x, y, point.est = bad, plot = FALSE,
              residual.plot = FALSE),
      "exactly match"
    )
  })
})

test_that("REG03: R2: Racine-Hastie squared correlation, bounded, degenerate series, reflection geometry", {
  withr::local_seed(20261108)
  aud_case("Racine-Hastie R2 is the squared observed-fitted correlation", {
    actual <- c(-2, -1, 0, 1, 2)
    predicted <- c(8, 4, 0, -4, -8)

    expected <- stats::cor(actual, predicted)^2
    expect_equal(NNS:::.nns_reg_r2(actual, predicted), expected, tolerance = 1e-14)
    expect_equal(expected, 1)
  })
  aud_case("Racine-Hastie R2 remains bounded when predictive R2 is negative", {
    actual <- c(-2, -1, 0, 1, 2)
    predicted <- c(20, 10, 0, -10, -20)

    predictive_r2 <- 1 - sum((actual - predicted)^2) /
      sum((actual - mean(actual))^2)
    nns_r2 <- NNS:::.nns_reg_r2(actual, predicted)

    expect_lt(predictive_r2, 0)
    expect_gte(nns_r2, 0)
    expect_lte(nns_r2, 1)
    expect_equal(nns_r2, stats::cor(actual, predicted)^2, tolerance = 1e-14)
  })
  aud_case("Racine-Hastie R2 handles constant-series degeneracy explicitly", {
    expect_equal(NNS:::.nns_reg_r2(rep(3, 5), rep(3, 5)), 1)
    expect_equal(NNS:::.nns_reg_r2(rep(3, 5), rep(4, 5)), 0)
    expect_equal(NNS:::.nns_reg_r2(1:5, rep(3, 5)), 0)
  })
  aud_case("reported NNS.reg and NNS.M.reg R2 values are bounded", {
    x <- seq(-2, 2, length.out = 40)
    y <- sin(4 * x) + x^2
    univariate <- NNS.reg(x, y, order = 1, plot = FALSE, residual.plot = FALSE)

    X <- cbind(x, x^2)
    multivariate <- NNS.M.reg(
      X, y, order = 1, n.best = 1,
      plot = FALSE, residual.plot = FALSE
    )

    expect_gte(univariate$R2, 0)
    expect_lte(univariate$R2, 1)
    expect_gte(multivariate$R2, 0)
    expect_lte(multivariate$R2, 1)
  })
  aud_case("response reflection preserves default mean-regression geometry", {
    set.seed(128)
    x <- data.frame(a = rnorm(140), b = rnorm(140))
    y <- x$a^2 + 2 * x$b + rnorm(140, sd = .3)
    fit <- function(response) NNS.reg(x, response, point.est = x[1:10, ],
      noise.reduction = "mean", plot = FALSE, ncores = 1)
    a <- fit(y); b <- fit(-y)
    expect_equal(a$n.best, b$n.best)
    expect_equal(a$rhs.partitions, b$rhs.partitions)
    expect_equal(a$Point.est, -b$Point.est, tolerance = 1e-10)
    expect_equal(a$R2, b$R2, tolerance = 1e-10)
  })
  aud_case("Racine-Hastie fit is squared correlation rather than one minus relative SSE", {
    y <- seq(-3, 4, length.out = 30)
    prediction <- 2 * y + 5
    expect_equal(NNS:::.nns_reg_r2(y, prediction), 1)
    expect_lt(1 - sum((y - prediction)^2) / sum((y - mean(y))^2), 0)
  })
})

test_that("REG04: regression or classification: REG keeps numeric binary responses, counts and XONLY stay regression, CLASS ties and numeric class values, explicit n.best", {
  withr::local_seed(20261108)
  aud_case("REG keeps numeric binary responses on the regression surface", {
    x <- c(0, 1, 2, 3)
    y <- c(0, 0, 1, 1)
    fit <- NNS.reg(x, y, type = "REG", order = "max",
                   point.est = 1.5, plot = FALSE)
    expect_equal(fit$Point.est, .5)
    expect_null(fit$Prediction.Accuracy)
    expect_null(fit$class.levels)
    expect_equal(NNS.reg(x, y, type = "reg", order = "max",
                         point.est = 1.5, plot = FALSE)$Point.est, .5)
    expect_type(NNS.reg(x, y, plot = FALSE)$Prediction.Accuracy, "double")
    expect_error(NNS.reg(x, factor(y), type = "REG", plot = FALSE),
                 "requires a numeric response")
  })
  aud_case("REG reaches full multivariate and dimension reduction fits", {
    t <- seq_len(40)
    x <- data.frame(a = t, b = sin(t / 5))
    y <- rep(c(0, 1), 20)
    for (method in list(NULL, "cor", "NNS.dep")) {
      fit <- NNS.reg(x, y, type = "REG", dim.red.method = method,
                     point.est = x[1:3, ], plot = FALSE, ncores = 1)
      expect_null(fit$Prediction.Accuracy)
      expect_null(fit$class.levels)
      expect_true(all(is.finite(fit$Point.est)))
      if (!is.null(method)) {
        expect_s3_class(fit$equation, "data.frame")
        expect_length(fit$equation$Coefficient, ncol(x) + 1L)
        expect_true(all(is.finite(fit$equation$Coefficient)))
      }
    }
  })
  aud_case("Stack REG preserves numeric response scale in both components", {
    t <- seq_len(40)
    x <- data.frame(a = t, b = sin(t / 5))
    run <- function(y) NNS.stack(x, y, IVs.test = x[1:3, ], type = "REG",
      method = c(1, 2), folds = 2, ncores = 1, seed = 17, status = FALSE,
      pred.int = .8)
    binary <- run(rep(c(0, 1), 20))
    signed <- run(rep(c(-5, 5), 20))
    expect_null(binary$class.levels)
    expect_null(signed$class.levels)
    for (component in c("reg", "dim.red", "stack"))
      expect_equal(signed[[component]], 10 * binary[[component]] - 5,
                   tolerance = 1e-10)
    for (component in c("reg.pred.int", "dim.red.pred.int", "pred.int"))
      expect_equal(as.matrix(signed[[component]]),
                   10 * as.matrix(binary[[component]]) - 5, tolerance = 1e-10)
  })
  aud_case("numeric derivative responses use regression even when binary", {
    y <- rep(c(0, 1), 15)
    expect_equal(dy.dx(seq_along(y), y, eval.point = "overall"), 1 / 29)
    expect_equal(dy.dx(seq_along(y), 10 * y - 5, eval.point = "overall"), 10 / 29)
  })
  aud_case("classification returns numeric observed class values", {
    train <- iris[-(141:150), ]
    test <- iris[141:150, 1:4]
    fit <- NNS.reg(train[, 1:4], train[, 5], point.est = test,
                   type = "CLASS", plot = FALSE,
                   residual.plot = FALSE, ncores = 1)
    expect_true(is.numeric(fit$Point.est))
    expect_true(all(fit$Point.est %in% 1:3))
    expect_identical(fit$class.levels, levels(iris$Species))
  })
  aud_case("count responses remain regression unless CLASS is explicit", {
    x <- seq_len(60)
    y <- rep(1:3, each = 20)
    regression <- NNS.reg(x, y, plot = FALSE)
    expect_null(regression$Prediction.Accuracy)

    classification <- NNS.reg(x, y, type = "CLASS", plot = FALSE)
    expect_true(is.numeric(classification$Prediction.Accuracy))
  })
  aud_case("XONLY does not activate classification", {
    x <- seq_len(40)
    y <- sin(x / 4)
    fit <- NNS.reg(x, y, type = "XONLY", plot = FALSE)
    expect_null(fit$Prediction.Accuracy)
  })
  aud_case("NNS.reg CLASS predicts one of a cell's tied classes", {
    # The cell x = 1:4 holds classes 1, 1, 3, 3; the mean of the tied modes made it class 2.
    x <- 1:8; y <- c(1, 1, 3, 3, 2, 2, 2, 2)
    f <- NNS.reg(x, y, type = "CLASS", order = 2, plot = FALSE)
    expect_true(all(f$Fitted.xy$y.hat[1:4] %in% c(1, 3)))
  })
  aud_case("categorical NNS.reg honours an explicit n.best", {
    x <- data.frame(a = factor(c("a", "a", "b")), b = factor(c("u", "v", "u")))
    q <- data.frame(a = "b", b = "v")   # an unobserved joint combination
    est <- function(nb) NNS.reg(x, c(0, 10, 20), point.est = q, n.best = nb,
      noise.reduction = "mean", plot = FALSE, residual.plot = FALSE)$Point.est
    expect_equal(est(NULL), 15)         # the two nearest cells, tied
    expect_equal(est("all"), est(3))
    expect_false(isTRUE(all.equal(est("all"), est(NULL))))
  })
})

test_that("REG05: argument contract: warning settings, dimensions, factor encoding, missing predictions, distance modes, part order max", {
  withr::local_seed(20261108)
  aud_case("NNS.reg preserves warning settings on success and failure", {
    old <- getOption("warn")
    on.exit(options(warn = old), add = TRUE)
    options(warn = 1)

    expect_silent(NNS.reg(1:20, (1:20)^2, plot = FALSE))
    expect_equal(getOption("warn"), 1)

    expect_error(NNS.reg(1:20, 1:10, plot = FALSE), "rows")
    expect_equal(getOption("warn"), 1)
  })
  aud_case("training and response dimensions are checked before fitting", {
    expect_error(
      NNS.reg(matrix(1:20, ncol = 2), 1:5, plot = FALSE),
      "rows"
    )
  })
  aud_case("factor encoding is training-fitted and unseen levels fail", {
    x <- data.frame(group = factor(rep(c("a", "b"), each = 15)),
                    z = seq_len(30))
    y <- rep(c(1, 2), each = 15)
    ok <- data.frame(group = factor(c("a", "b"), levels = c("a", "b")),
                     z = c(4, 20))
    expect_silent(
      NNS.reg(x, y, point.est = ok, type = "CLASS",
              plot = FALSE, residual.plot = FALSE)
    )

    bad <- data.frame(group = "c", z = 10)
    expect_error(
      NNS.reg(x, y, point.est = bad, type = "CLASS",
              plot = FALSE, residual.plot = FALSE),
      "unseen level"
    )
  })
  aud_case("distance modes are validated and recorded", {
    set.seed(6)
    x <- cbind(a = rnorm(30), b = rnorm(30))
    y <- x[, 1] + x[, 2]
    for (d in c("NNS", "L1", "L2", "FACTOR")) {
      fit <- NNS.reg(x, y, point.est = x[1:2, ], dist = d,
                     plot = FALSE, residual.plot = FALSE)
      expect_identical(fit$dist, d)
    }
    expect_error(NNS.reg(x, y, dist = "DTW", plot = FALSE), "dist")
  })
  aud_case("missing and nonfinite prediction values fail", {
    x <- cbind(a = 1:20, b = 21:40)
    y <- 1:20
    expect_error(
      NNS.reg(x, y, point.est = matrix(c(NA, 2), nrow = 1),
              plot = FALSE, residual.plot = FALSE),
      "finite"
    )
  })
  aud_case("NNS.part order max never passes a character value to C++", {
    fit <- NNS.part(c(1, 1, 2, 3), c(1, 3, 4, 9),
                    order = "max", type = "XONLY")
    expect_equal(nrow(fit$regression.points), 3)
    expect_identical(fit$order, 3L)
  })
})

test_that("REG06: predictions agree with their fit: dim.red rows, x* equation, zero weights, order max, one row, n.best one, interval rows", {
  withr::local_seed(20261108)
  aud_case("dimension reduction is invariant to unrelated prediction rows", {
    set.seed(2)
    x <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
    y <- x$a + x$b^2 - x$c
    base <- x[1:3, ]
    extra <- rbind(base, data.frame(a = 1e6, b = -1e6, c = 1e6))

    f1 <- NNS.reg(x, y, point.est = base, dim.red.method = "cor",
                  plot = FALSE, residual.plot = FALSE)
    f2 <- NNS.reg(x, y, point.est = extra, dim.red.method = "cor",
                  plot = FALSE, residual.plot = FALSE)
    expect_equal(f1$Point.est, f2$Point.est[1:3], tolerance = 1e-12)
  })
  aud_case("numeric dimension-reduction equation reproduces X star", {
    set.seed(3)
    x <- cbind(a = rnorm(40), b = rnorm(40), c = rnorm(40))
    y <- rowSums(x)
    coef <- c(2, -1, 0.5)
    fit <- NNS.reg(x, y, dim.red.method = coef, plot = FALSE)

    mins <- apply(x, 2, min)
    ranges <- apply(x, 2, max) - mins
    norm <- matrix(0.5, nrow(x), ncol(x))
    active <- ranges > 0
    norm[, active] <- sweep(sweep(x[, active, drop = FALSE], 2,
                                  mins[active], "-"), 2, ranges[active], "/")
    expected <- as.numeric(norm %*% coef / sum(abs(coef) > 0))
    expect_equal(as.numeric(fit$x.star$x), expected, tolerance = 1e-12)
    expect_equal(tail(fit$equation$Coefficient, 1), sum(abs(coef) > 0))
  })
  aud_case("all-zero dimension weights fall back to a defined model", {
    x <- cbind(a = rep(1, 30), b = rep(2, 30))
    y <- seq_len(30)
    fit <- NNS.reg(x, y, dim.red.method = c(0, 0), plot = FALSE)
    expect_true(all(is.finite(fit$x.star$x)))
    expect_gt(tail(fit$equation$Coefficient, 1), 0)
  })
  aud_case("order max fitted values use the same prediction function", {
    x <- c(1, 1, 2, 3, 4, 5)
    y <- c(1, 3, 4, 9, 16, 25)
    fit <- NNS.reg(x, y, order = "max", point.est = x, plot = FALSE)
    expect_equal(fit$Fitted.xy$y.hat, fit$Point.est, tolerance = 1e-12)
  })
  aud_case("one-row matrix and vector predictions are identical", {
    set.seed(4)
    x <- cbind(a = rnorm(50), b = rnorm(50))
    y <- x[, 1] - x[, 2]
    v <- unname(x[5, ])
    m <- matrix(v, nrow = 1)

    fv <- NNS.reg(x, y, point.est = v, plot = FALSE,
                  residual.plot = FALSE, ncores = 1)
    fm <- NNS.reg(x, y, point.est = m, plot = FALSE,
                  residual.plot = FALSE, ncores = 1)
    expect_equal(fv$Point.est, fm$Point.est, tolerance = 1e-12)
  })
  aud_case("n.best one fitted values equal predictions on training rows", {
    set.seed(5)
    x <- cbind(a = rnorm(40), b = rnorm(40))
    y <- x[, 1]^2 + x[, 2]
    fit <- NNS.reg(x, y, point.est = x, n.best = 1,
                   plot = FALSE, residual.plot = FALSE, ncores = 1)
    expect_equal(fit$Fitted.xy$y.hat, fit$Point.est, tolerance = 1e-12)
  })
  aud_case("prediction intervals preserve prediction row count", {
    x <- seq_len(40)
    y <- x + sin(x)
    points <- c(-10, 1, 20, 50)
    fit <- NNS.reg(x, y, point.est = points,
                   confidence.interval = 0.9, plot = FALSE)
    expect_equal(nrow(fit$pred.int), length(points))
    expect_true(all(fit$pred.int$pred.int.neg <= fit$pred.int$pred.int.pos))
  })
})
