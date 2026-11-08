# NNS.stack: distance and backend, native searches, inputs, scoring and the
# validation pipeline.
#
# NNS.stack classification weights its two models as NNS 13.0 did: by their
# squared objective scores, or inverse squares for a minimized objective.

test_that("STK01: distance and backend: NULL is NNS, method 1 OOF search, native backend options", {
  withr::local_seed(20261108)
  aud_case("NNS.stack NULL distance matches explicit NNS", {
    x <- data.frame(a = c(0, 1, 0.2, 0.9, 0.4, 0.8, 0.3, 0.7),
                    b = c(0, 0.2, 1, 0.9, 0.7, 0.1, 0.6, 0.4))
    y <- c(0, 2, 4, 8, 5, 3, 6, 7)
    args <- list(IVs.train = x, DV.train = y, IVs.test = x[1:2, ], method = 1,
                 folds = 2, ncores = 1, status = FALSE, stack = FALSE, seed = 42)
    a <- do.call(NNS.stack, c(args, list(dist = NULL)))
    b <- do.call(NNS.stack, c(args, list(dist = "NNS")))
    expect_equal(a$reg, b$reg, tolerance = 1e-12)
  })
  aud_case("NNS.stack Method 1 OOF search uses NNS distance for NULL and NNS", {
    x <- data.frame(a = c(0, 0.2, 0.9, 1, 0.15, 0.85, 0.35, 0.65, 0.45, 0.55),
                    b = c(0, 1, 0.2, 1, 0.85, 0.15, 0.7, 0.3, 0.55, 0.45))
    y <- c(0, 10, 20, 30, 11, 19, 12, 18, 14, 16)
    args <- list(IVs.train = x, DV.train = y, IVs.test = x[1:3, ], method = 1,
                 folds = 2, ncores = 1, status = FALSE, stack = FALSE, seed = 7)

    nns_null <- do.call(NNS.stack, c(args, list(dist = NULL)))
    nns_alias <- do.call(NNS.stack, c(args, list(dist = "NNS")))
    l2 <- do.call(NNS.stack, c(args, list(dist = "L2")))

    expect_identical(nns_null$NNS.reg.n.best, nns_alias$NNS.reg.n.best)
    expect_equal(nns_null$OBJfn.reg, nns_alias$OBJfn.reg, tolerance = 1e-12)
    expect_equal(nns_null$reg, nns_alias$reg, tolerance = 1e-12)
    expect_false(isTRUE(all.equal(nns_null$NNS.reg.n.best, l2$NNS.reg.n.best)) &&
                   isTRUE(all.equal(nns_null$OBJfn.reg, l2$OBJfn.reg)) &&
                   isTRUE(all.equal(nns_null$reg, l2$reg)))
  })
  aud_case("native backend options are accepted", {
    old <- options(NNS.native.stack = TRUE, NNS.native.mreg = TRUE, NNS.native.univariate = TRUE)
    on.exit(options(old), add = TRUE)
    expect_true(isTRUE(getOption("NNS.native.stack")))
  })
})

test_that("STK02: native method 1 and 2 searches select what the references select", {
  withr::local_seed(20261108)
  aud_case("method 1 native and reference select identical n.best on regression smoke case", {
    set.seed(11); x <- data.frame(a = rnorm(24), b = rnorm(24)); y <- x$a - x$b^2
    ref <- stack_with_native(FALSE, NNS.stack(x, y, IVs.test = x[1:3,], method = 1, folds = 2, ncores = 1, status = FALSE))
    nat <- stack_with_native(TRUE, NNS.stack(x, y, IVs.test = x[1:3,], method = 1, folds = 2, ncores = 1, status = FALSE))
    expect_identical(nat$NNS.reg.n.best, ref$NNS.reg.n.best)
    expect_numeric_close(nat$reg, ref$reg)
  })
  aud_case("method 2 native and reference select identical dimension count on smoke case", {
    set.seed(12); x <- data.frame(a = rnorm(24), b = rnorm(24), c = rnorm(24)); y <- x$a + x$b
    ref <- stack_with_native(FALSE, NNS.stack(x, y, IVs.test = x[1:3,], method = 2, folds = 2, ncores = 1, status = FALSE))
    nat <- stack_with_native(TRUE, NNS.stack(x, y, IVs.test = x[1:3,], method = 2, folds = 2, ncores = 1, status = FALSE))
    expect_identical(nat$NNS.dim.red.threshold, ref$NNS.dim.red.threshold)
    expect_numeric_close(nat$dim.red, ref$dim.red)
  })
})

test_that("STK03: inputs: character responses and predictors, duplicate predictors and names, empty factor levels, one test row, constant fold predictor, single-member class", {
  withr::local_seed(20261108)
  stack_boost <- fx_stack_boost()
  aud_case("NNS.stack codes a character response like the equivalent factor", {
    idx <- c(101:140, 1:100)
    x <- iris[idx, 1:4]
    y <- as.character(iris$Species[idx])
    run <- function(dv) NNS.stack(x, dv, IVs.test = iris[141:150, 1:4], method = 1,
                                  folds = 2, status = FALSE, ncores = 1, seed = 1)
    chr <- run(y)
    expect_identical(chr$class.levels, c("setosa", "versicolor", "virginica"))
    expect_identical(chr, run(factor(y)))
  })
  aud_case("NNS.stack encodes a character predictor like the equivalent factor", {
    set.seed(4)
    x <- data.frame(x = rnorm(60), g = rep(c("zeta", "alpha", "mid"), 20))
    y <- x$x + rep(c(2, 3, 0), 20) + rnorm(60)
    xf <- x
    xf$g <- factor(xf$g)
    run <- function(iv) NNS.stack(iv, y, IVs.test = iv[1:5, ], method = 2, folds = 2,
                                  status = FALSE, ncores = 1, seed = 1)
    chr <- run(x)
    expect_identical(names(chr$NNS.dim.red.coefficients), c("x", "g__alpha", "g__mid", "g__zeta"))
    expect_identical(chr, run(xf))
  })
  aud_case("NNS.stack balancing keeps a single-member class", env = stack_boost, {
    set.seed(9)
    x <- data.frame(a = rnorm(60), b = rnorm(60))
    frame <- capture_frame("NNS.stack", NNS::NNS.stack(
      x, factor(rep(1:3, c(29, 29, 2))), IVs.test = x[1:3, ], method = 1,
      balance = TRUE, folds = 2, seed = 1, ncores = 1, status = FALSE
    ))
    set.seed(1)
    idx <- frame$.balance_indices(balance_y)
    expect_identical(as.integer(table(balance_y[idx])), c(4L, 4L, 4L))
  })
  aud_case("NNS.stack accepts duplicate predictor columns like NNS.reg (cbind(x, x))", {
    set.seed(123)
    x <- rnorm(100); y <- rnorm(100)

    # The cbind(x, x) dimension trick works for NNS.reg; NNS.stack must not
    # error with "[IVs.train] predictor names must be unique." any more.
    expect_error(NNS.reg(x = cbind(x, x), y = y, plot = FALSE), NA)
    res <- NNS.stack(IVs.train = cbind(x, x), DV.train = y,
                     method = 1, ncores = 1, folds = 2, status = FALSE)
    expect_type(res, "list")
    expect_true("reg" %in% names(res))
  })
  aud_case("NNS.stack aligns duplicate-named test columns (cbind(x, x)) instead of erroring", {
    set.seed(123)
    x <- rnorm(60); y <- rnorm(60)
    train <- cbind(x, x)

    expect_error(
      NNS.stack(IVs.train = train, DV.train = y, IVs.test = train[1:5, ],
                method = 1, ncores = 1, folds = 2, status = FALSE),
      NA
    )
  })
  aud_case("a declared factor level with no training rows is a valid test value (#150)", {
    x <- data.frame(g = factor(c(rep(c("a", "b"), 19), "a", "rare")))
    y <- factor(rep(c("L", "R"), 20))
    fit <- NNS.boost(x, y, ts.test = 8, learner.trials = 1, epochs = 0,
                     features.only = TRUE, feature.importance = FALSE,
                     status = FALSE, seed = 9)
    expect_equal(unname(unlist(fit$feature.weights)), 1)

    g <- factor(c("a", "b", "a", "b"), levels = c("a", "b", "c"))
    at <- function(level) NNS.reg(data.frame(g = g), c(1, 2, 1, 2),
      point.est = data.frame(g = factor(level, levels = levels(g))),
      plot = FALSE, residual.plot = FALSE)$Point.est
    expect_equal(at("a"), 1)
    expect_equal(at("b"), 2)
    expect_true(is.finite(at("c")))
    # Level order still does not matter.
    g2 <- factor(g, levels = c("c", "b", "a"))
    expect_equal(NNS.reg(data.frame(g = g2), c(1, 2, 1, 2),
      point.est = data.frame(g = factor("c", levels = levels(g2))),
      plot = FALSE, residual.plot = FALSE)$Point.est, at("c"))

    set.seed(1)
    gg <- factor(sample(c("a", "b"), 60, TRUE), levels = c("a", "b", "c"))
    zz <- rnorm(60)
    s <- NNS.stack(data.frame(g = gg, z = zz), zz + (gg == "b"),
      data.frame(g = factor(c("a", "c"), levels = levels(gg)), z = c(0, 0)),
      folds = 2, ncores = 1, status = FALSE)
    expect_true(all(is.finite(s$stack)))
    # An undeclared character value is still unseen.
    expect_error(NNS.reg(data.frame(g = c("a", "b", "a", "b")), c(1, 2, 1, 2),
      point.est = data.frame(g = "c"), plot = FALSE), "unseen level")
  })
  aud_case("stack accepts a single test row with a factor predictor (#150)", {
    # vapply() returned the one-row dummy block as a vector, which was then
    # reshaped into one column per row instead of one row.
    g <- factor(rep(c("a", "b"), 20))
    z <- seq(0, 1, length.out = 40)
    fit <- NNS.stack(data.frame(g = g, z = z), z + (g == "b"),
      data.frame(g = factor("b", levels = c("a", "b")), z = .5),
      folds = 2, ncores = 1, status = FALSE)
    expect_length(fit$stack, 1)
    expect_true(abs(fit$stack - 1.5) < .25)
  })
  aud_case("stack keeps the same method 1 columns when a fold has a constant predictor (#150)", {
    # b is constant in the first chronological training fold, so its fold
    # coefficient is 0. Dropping it gave that fold a univariate design with no
    # neighbour counts to tune, and the stack stopped.
    x <- data.frame(a = rep(1:10, 6), b = c(rep(0, 40), rep(c(-10, 10), 10)))
    y <- x$a + 2 * x$b
    fit <- NNS.stack(x, y, x[51:60, ], ts.test = 10, folds = 2, order = "max",
                     dim.red.method = "cor", ncores = 1, status = FALSE)
    expect_length(fit$stack, 10)
    expect_true(all(is.finite(fit$stack)))
    expect_lt(max(abs(fit$stack - y[51:60])), 2)
  })
})

test_that("STK04: scoring: ALL candidates, single-predictor n.best and OBJfn.reg, NNS 13.0 class weights, near-equal blend, inverse squares, neighbour-count probabilities", {
  withr::local_seed(20261108)
  stack_boost <- fx_stack_boost()
  stack_univariate_k <- fx_stack_univariate_k()
  aud_case("NNS.stack scores the multivariate ALL candidate with NNS.reg(n.best = 'all')", env = stack_boost, {
    set.seed(5)
    x <- data.frame(a = rnorm(160), b = rnorm(160))
    y <- rnorm(160)
    frame <- capture_frame("NNS.stack", NNS::NNS.stack(
      x, y, IVs.test = x[1:3, ], method = 1, seed = 1, ncores = 1, status = FALSE
    ))
    oof <- rep(NA_real_, length(y))
    for (s in frame$splits) {
      design <- frame$.method1_design_for_split(s$train, s$validation)
      fit <- NNS.reg(design$train, frame$y_fit_all[s$train], point.est = design$test,
                     n.best = "all", order = frame$order, dist = frame$dist, type = "REG",
                     factor.2.dummy = FALSE, plot = FALSE, residual.plot = FALSE,
                     ncores = 1, point.only = TRUE)
      oof[s$validation] <- fit$Point.est
    }
    expect_equal(frame$candidate_raws$all, oof, tolerance = 1e-12)
    expect_equal(unname(frame$candidate_scores["all"]), sum((oof - y)^2), tolerance = 1e-12)
  })
  aud_case("NNS.stack scores the univariate ALL candidate with the deployed NNS.reg fit", env = stack_boost, {
    set.seed(5)
    x <- data.frame(a = rnorm(160))
    y <- rnorm(160)
    frame <- capture_frame("NNS.stack", NNS::NNS.stack(
      x, y, IVs.test = x[1:3, , drop = FALSE], method = 1, seed = 1, ncores = 1, status = FALSE
    ))
    oof <- rep(NA_real_, length(y))
    for (s in frame$splits) {
      design <- frame$.method1_design_for_split(s$train, s$validation)
      fit <- NNS.reg(design$train[[1L]], frame$y_fit_all[s$train], point.est = design$test[[1L]],
                     n.best = "all", order = frame$order, dist = frame$dist, type = "REG",
                     factor.2.dummy = FALSE, plot = FALSE, residual.plot = FALSE,
                     ncores = 1, point.only = TRUE)
      oof[s$validation] <- fit$Point.est
    }
    expect_equal(frame$candidate_raws$all, oof, tolerance = 1e-12)
    expect_equal(unname(frame$candidate_scores["all"]), sum((oof - y)^2), tolerance = 1e-12)
  })
  aud_case("NNS.stack reports no n.best for a single predictor", env = stack_univariate_k, {
    s <- NNS.stack(wt, mtcars$mpg, method = 1, seed = 1, ncores = 1, status = FALSE)
    expect_identical(s$NNS.reg.n.best, NA_integer_)
  })
  aud_case("NNS.stack OBJfn.reg for a single predictor scores the deployed NNS.reg fit", env = stack_univariate_k, {
    y <- mtcars$mpg
    frame <- capture_frame("NNS.stack", s <- NNS::NNS.stack(
      wt, y, method = 1, seed = 1, ncores = 1, status = FALSE
    ))
    univariate_fit <- function(train_x, train_y, test_x) {
      NNS.reg(train_x, train_y, point.est = test_x, order = frame$order, type = "REG",
              factor.2.dummy = FALSE, plot = FALSE, residual.plot = FALSE,
              ncores = 1)$Point.est
    }
    oof <- rep(NA_real_, length(y))
    for (sp in frame$splits) {
      design <- frame$.method1_design_for_split(sp$train, sp$validation)
      oof[sp$validation] <- univariate_fit(design$train[[1L]], frame$y_fit_all[sp$train],
                                           design$test[[1L]])
    }
    expect_equal(s$OBJfn.reg, sum((oof - y)^2), tolerance = 1e-12)
    expect_equal(s$reg, univariate_fit(frame$full_design$train[[1L]], frame$y_fit_all,
                                       frame$full_design$test[[1L]]), tolerance = 1e-12)
  })
  aud_case("the weights follow NNS 13.0's rule", {
    expect_equal(.nns_stack_objective_weights(.9, .6, "max"),
                 c(reg = .81, dim.red = .36) / 1.17)
    expect_equal(.nns_stack_objective_weights(2, 1, "min"), c(reg = .2, dim.red = .8))
    # A zero loss counts as 1e-10, so that model takes essentially all the weight.
    expect_equal(.nns_stack_objective_weights(0, 1, "min"), c(reg = 1, dim.red = 0))
    expect_equal(.nns_stack_objective_weights(0, 0, "max"), c(reg = .5, dim.red = .5))
    expect_equal(.nns_stack_objective_weights(NA, .8, "max"), c(reg = 0, dim.red = 1))
  })
  aud_case("nearly equal objective scores give a nearly equal blend", {
    # Out-of-fold accuracies of 130 and 131 of 138 rows. Choosing the blend
    # weight by its own out-of-fold accuracy gave dim.red all the weight, and the
    # stack then misclassified the sixth test row, which reg classifies correctly.
    set.seed(123)
    test.set <- sort(unlist(lapply(split(seq_len(nrow(iris)), iris$Species),
                                   function(rows) sample(rows, 4))))
    b <- NNS.stack(IVs.train = iris[-test.set, 1:4], DV.train = iris[-test.set, 5],
                   IVs.test = iris[test.set, 1:4], type = "CLASS", balance = TRUE,
                   ncores = 1, folds = 5, status = FALSE)
    expect_equal(c(b$OBJfn.reg, b$OBJfn.dim.red), c(130, 131) / 138)
    expect_equal(b$weights, c(reg = 130^2, dim.red = 131^2) / (130^2 + 131^2))
    expect_equal(unname(b$probabilities),
                 unname(b$weights[["reg"]] * b$reg.probabilities +
                          b$weights[["dim.red"]] * b$dim.red.probabilities))
    expect_equal(b$stack, as.numeric(iris[test.set, 5]))
  })
  aud_case("a minimized objective weights the models by inverse squares", {
    set.seed(5)
    test.set <- sample(nrow(iris), 30)
    b <- NNS.stack(IVs.train = iris[-test.set, 1:4], DV.train = iris[-test.set, 5],
                   IVs.test = iris[test.set, 1:4], type = "CLASS",
                   obj.fn = expression(mean(predicted != actual)), objective = "min",
                   ncores = 1, folds = 5, status = FALSE)
    s <- pmax(c(b$OBJfn.reg, b$OBJfn.dim.red), 1e-10)
    expect_equal(unname(b$weights), (1 / s^2) / sum(1 / s^2))
  })
  aud_case("class probabilities over a neighbour-count path equal separate fits (#150)", {
    set.seed(11)
    X <- cbind(u = round(runif(120), 1), v = round(runif(120), 1))  # tied distances
    cls <- sample(1:3, 120, TRUE)
    P <- outer(cls, 1:3, "==") * 1
    Q <- cbind(u = round(runif(25), 1), v = round(runif(25), 1))
    ks <- list(1L, 2L, 3L, 7L, "all")
    path <- NNS:::.nns_class_base(X, P, Q, NULL, NULL, "NNS", 1L,
                                  need.fitted = FALSE, k.path = ks)$predicted.path
    for (a in seq_along(ks)) {
      one <- NNS:::.nns_class_base(X, P, Q, NULL, ks[[a]], "NNS", 1L)$predicted
      expect_identical(path[[a]], one)
    }
    # Point-only classification skips the fitted classes but not the predictions.
    full <- NNS.reg(X, factor(cls), point.est = Q, plot = FALSE, residual.plot = FALSE)
    po <- NNS.reg(X, factor(cls), point.est = Q, plot = FALSE, residual.plot = FALSE,
                  point.only = TRUE)
    expect_identical(po$Point.est, full$Point.est)
    expect_identical(po$Point.probabilities, full$Point.probabilities)
    expect_null(po$Fitted.xy)
  })
})

test_that("STK05: pipeline validation: supplied prefix only, selected coordinates, split modes, two-level regression", {
  withr::local_seed(20261108)
  stack_pipeline_validation <- fx_stack_pipeline_validation()
  aud_case("forecast validation uses only the supplied prefix and generated predictors", env = stack_pipeline_validation, {
    set.seed(131)
    x <- data.frame(a = rnorm(45), b = rnorm(45), c = rnorm(45))
    y <- x$a - 2 * x$b + .1 * rnorm(45)
    validation <- list(IVs.train = x[1:40, ], DV.train = y[1:40],
                       IVs.test = x[41:43, ] + 1, DV.test = y[41:43])
    run <- function(xx, yy) .capture_pipeline_stack(xx, yy, IVs.test = x[41:45, ],
      .validation = validation, ts.test = 5, folds = 1, ncores = 1, status = FALSE)
    a <- run(x, y)
    changed <- x; changed[41:45, ] <- changed[41:45, ] + 100
    yy <- y; yy[41:45] <- yy[41:45] - 200
    b <- run(changed, yy)
    expect_equal(a$result$OBJfn.reg, b$result$OBJfn.reg)
    expect_equal(a$result$OBJfn.dim.red, b$result$OBJfn.dim.red)
    expect_equal(a$frame$dim_best_count, b$frame$dim_best_count)
    expect_equal(a$frame$candidate_raws, b$frame$candidate_raws)
    expect_equal(a$frame$splits[[1]]$train, 1:40)
    expect_equal(a$frame$splits[[1]]$validation, 41:43)
    expect_equal(a$frame$x.cv[41:43, ], validation$IVs.test)
    # The final fit still receives all original training rows.
    expect_equal(b$frame$x, changed)
    expect_equal(b$frame$y_fit_all, yy)
    expect_false(isTRUE(all.equal(a$result$stack, b$result$stack)))
  })
  aud_case("stack method one retains selected coordinates instead of duplicate indices", env = stack_pipeline_validation, {
    set.seed(128)
    x <- data.frame(a = rnorm(80), b = rnorm(80), c = rnorm(80))
    y <- x$a * x$b + .5 * x$c
    fit <- .capture_pipeline_stack(x, y, IVs.test = x[1:3, ], folds = 2,
      ncores = 1, status = FALSE, dim.red.method = c(a = 1, b = -1, c = .1))
    fr <- fit$frame
    selected <- which(abs(fr$dim_full_coef) > 0)
    expect_equal(fr$.method1_full_design()$train,
                 fr$full_design$train[, selected, drop = FALSE])
    if (length(selected) > 1) expect_false(isTRUE(all.equal(
      fr$.method1_full_design()$train[[1]], fr$.method1_full_design()$train[[2]])))
  })
  aud_case("pipeline validation rejects incompatible split modes", {
    x <- data.frame(a = 1:10, b = (1:10)^2)
    y <- seq(.1, 1, .1)
    v <- list(IVs.train = x[1:8, ], DV.train = y[1:8],
              IVs.test = x[9:10, ], DV.test = y[9:10])
    expect_error(NNS.stack(x, y, .validation = v, ts.test = 2, folds = 2,
      ncores = 1, status = FALSE), "folds = 1")
    v$DV.test <- NA_real_
    expect_error(NNS.stack(x, y, .validation = v, ts.test = 2, folds = 1,
      ncores = 1, status = FALSE), "finite held-out")
  })
  aud_case("numeric two-level forecast responses remain regression with a validation fold", {
    x <- data.frame(a = seq_len(30), b = sin(seq_len(30)))
    y <- rep(c(2, 5), 15)
    v <- list(IVs.train = x[1:25, ], DV.train = y[1:25],
              IVs.test = x[26:30, ], DV.test = y[26:30])
    run <- function(response) NNS.stack(x, response, IVs.test = x[26:30, ],
      .validation = v, ts.test = 5, folds = 1, method = 1,
      ncores = 1, status = FALSE)
    a <- run(y)
    changed <- y; changed[26:30] <- 7
    b <- run(changed)
    expect_null(a$class.levels)
    expect_equal(a$OBJfn.reg, b$OBJfn.reg)
    expect_error(NNS.stack(x, y, .validation = v, ts.test = 5, folds = 1,
      type = "CLASS", ncores = 1, status = FALSE), "requires regression")
  })
})
