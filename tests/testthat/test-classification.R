# Classification: class cells, label and level invariance, named weights, and the
# dimension-reduction thresholds NNS.stack reports and NNS.reg accepts.
#
# NNS.stack reports the dimension-reduction cutoff that keeps its tuned count,
# and NNS.reg accepts what NNS.stack reports, NA included.

test_that("CLS01: class cells keep class frequencies, counts and numeric labels, including two-valued coordinates", {
  withr::local_seed(20261108)
  aud_case("classification interpolates class frequencies, not arbitrary class codes", {
    x <- rep(0:2, each = 4)
    labels <- rep(c("A", "C", "B"), each = 4)
    run <- function(levels) NNS.reg(x, factor(labels, levels = levels),
      order = "max", point.est = c(.25, .5, .75, 1.25, 1.5), plot = FALSE)
    a <- run(c("A", "B", "C")); b <- run(c("C", "A", "B"))
    expect_equal(a$class.levels[a$Point.est], b$class.levels[b$Point.est])
    expect_equal(unclass(a$Point.probabilities), b$Point.probabilities[, c("A", "B", "C")])
    expect_equal(a$Point.probabilities[2, ], c(A = .5, B = 0, C = .5))
    expect_equal(a$class.levels[a$Point.est[2]], "A")
    expect_equal(rowSums(a$Point.probabilities), rep(1, 5))
  })
  aud_case("class cells retain counts and numeric labels keep their identity", {
    x <- factor(rep(c("left", "right"), each = 4))
    y <- c(10, 10, 10, -7, -7, -7, -7, 10)
    fit <- NNS.reg(x, y, type = "CLASS", point.est = c("left", "right"), plot = FALSE)
    expect_equal(fit$Point.est, c(10, -7))
    expect_equal(unname(unclass(fit$Point.probabilities)), rbind(c(.25, .75), c(.75, .25)))
    relabel <- ifelse(y == 10, -100, 8)
    other <- NNS.reg(x, relabel, type = "CLASS", point.est = c("left", "right"), plot = FALSE)
    expect_equal(other$Point.est, c(-100, 8))
  })
  aud_case("two-valued numeric coordinates retain distinct class cells", {
    x <- rep(c(0, 1), 20)
    y <- factor(ifelse(x == 0, "zero", "one"))
    fit <- NNS.reg(x, y, point.est = c(0, .5, 1), plot = FALSE)
    expect_equal(nrow(fit$RPM), 2L)
    expect_equal(fit$Prediction.Accuracy, 1)
    expect_equal(unname(unclass(fit$Point.probabilities)),
                 rbind(c(0, 1), c(.5, .5), c(1, 0)))
    mixed <- data.frame(group = factor(x), continuous = rep(1:20, each = 2))
    cells <- NNS.reg(mixed, y, point.est = mixed[1:2, ], n.best = 1,
                     plot = FALSE, ncores = 1)
    expect_equal(cells$class.levels[cells$Point.est], c("zero", "one"))
    expect_true(all(apply(cells$RPM.probabilities, 1, max) == 1))
  })
})

test_that("CLS02: nominal label and level order do not change stack probabilities, decisions or categorical regression", {
  withr::local_seed(20261108)
  aud_case("stack probabilities and decoded decisions are invariant to nominal label order", {
    x <- iris[, 1:3]
    run <- function(y) NNS.stack(x, y, IVs.test = x[c(1, 52, 103), ],
      folds = 2, seed = 17, ncores = 1, status = FALSE)
    a <- run(iris$Species)
    b <- run(factor(iris$Species, levels = rev(levels(iris$Species))))
    expect_equal(a$class.levels[a$stack], b$class.levels[b$stack])
    expect_equal(unclass(a$probabilities), b$probabilities[, colnames(a$probabilities)])
    expect_equal(a$weights, b$weights)
    expect_equal(rowSums(a$probabilities), rep(1, 3))
    expect_true(all(is.finite(a$NNS.dim.red.threshold)))
    expect_true(all(is.finite(b$NNS.dim.red.threshold)))
    expect_true(a$NNS.dim.red.count %in% 1:3)
  })
  aud_case("numeric stack outputs and custom objectives see original class values", {
    x <- iris[1:100, 1:2]
    y <- rep(c(-5, 8), each = 50)
    observed <- numeric()
    record <- function(actual) observed <<- c(observed, actual)
    loss <- as.expression(bquote({ .(record)(actual); mean(predicted == actual) }))
    fit <- NNS.stack(x, y, x[c(1, 100), ], method = 1, obj.fn = loss,
      objective = "max", folds = 2, ncores = 1, status = FALSE, seed = 9)
    expect_setequal(unique(observed), c(-5, 8))
    expect_equal(fit$stack, c(-5, 8))
  })
  aud_case("nominal predictor level order does not alter categorical group regression", {
    y <- rep(c(2.1, 5.3, -1.8), each = 6) + rep(c(-.2, 0, .2), 6)
    labels <- rep(c("a", "b", "c"), each = 6)
    run <- function(levels) NNS.reg(factor(labels, levels = levels), y,
      point.est = c("a", "b", "c"), noise.reduction = "mean", plot = FALSE)
    a <- run(c("a", "b", "c")); b <- run(c("c", "a", "b"))
    expect_equal(a$Point.est, b$Point.est)
    expect_equal(a$Point.est, c(2.1, 5.3, -1.8))
  })
  aud_case("both stack components canonicalize unordered predictor levels", {
    set.seed(128)
    x <- data.frame(a = rnorm(90), group = factor(rep(c("z", "a", "m"), 30)))
    reordered <- x; reordered$group <- factor(x$group, levels = c("m", "z", "a"))
    response <- x$a + rep(c(1, 2, -1), 30) + rnorm(90, sd = .1)
    run <- function(data, y) NNS.stack(data, y, IVs.test = data[1:6, ],
      method = c(1, 2), folds = 2, seed = 19, ncores = 1, status = FALSE)
    expect_equal(run(x, response), run(reordered, response))
    target <- factor(rep(c("B", "A", "C"), 30))
    expect_equal(run(x, target), run(reordered, target))
  })
})

test_that("CLS03: named dimension weights and factor coefficient names align in regression, classification and Stack", {
  withr::local_seed(20261108)
  aud_case("named dimension weights align in regression and nominal classification", {
    set.seed(12)
    x <- data.frame(a = rep(c(-1, 1), 20), b = rnorm(40))
    labels <- factor(ifelse(x$a > 0, "A", "B"))
    for (response in list(labels, seq_len(nrow(x)) / 10)) {
      fit <- function(w) NNS.reg(x, response, point.est = x[1:4, ],
        dim.red.method = w, plot = FALSE)
      a <- fit(c(a = 1, b = 0)); b <- fit(c(b = 0, a = 1))
      expect_equal(a$equation, b$equation)
      expect_equal(a$Point.est, b$Point.est)
      expect_equal(a$Point.probabilities, b$Point.probabilities)
      expect_error(fit(c(a = 1, unknown = 0)), "exactly match")
    }
    fit <- function(w) NNS.stack(x, labels, x[1:4, ], dim.red.method = w,
      method = c(1, 2), folds = 2, ncores = 1, seed = 1, status = FALSE)
    expect_equal(fit(c(a = 1, b = 0)), fit(c(b = 0, a = 1)))
    expect_error(fit(c(a = 1, unknown = 0)), "exactly match")
  })
  aud_case("nominal Stack shares the factor coefficient names and fold columns", {
    x <- data.frame(group = factor(rep(c("a", "b"), 20)),
                    value = rep(1:20, each = 2))
    y <- factor(ifelse(x$group == "a", "L", "R"))
    fit <- function(w) NNS.stack(x, y, x[1:4, ], dim.red.method = w,
      method = 2, folds = 2, ncores = 1, seed = 9, status = FALSE)
    weights <- c(group__a = 1, group__b = -1, value = 0)
    a <- fit(weights); b <- fit(rev(weights))
    expect_equal(a, b)
    expect_equal(a$class.levels[a$stack], as.character(y[1:4]))
    expect_equal(a$NNS.dim.red.coefficients[[1]]$Variable,
                 c(names(weights), "DENOMINATOR"))
    # A level present only in one validation row keeps its design column even
    # when absent from that fold's fitted observations.
    x$group <- as.character(x$group); x$group[1] <- "rare"
    x$group <- factor(x$group)
    rare <- NNS.stack(x, y, x[1:4, ], method = c(1, 2), folds = 2,
                      ncores = 1, seed = 9, status = FALSE)
    expect_equal(rowSums(rare$probabilities), rep(1, 4))
  })
})

# pin: the logistic-regression example
test_that("CLS04: the dim.red threshold NNS.stack reports keeps its tuned count in NNS.reg: two classes, more classes, the logistic example", {
  withr::local_seed(20261108)
  dimred_threshold <- fx_dimred_threshold()
  aud_case("a two-class stack reports one cutoff that keeps its tuned count", env = dimred_threshold, {
    # As in the logistic regression example: a 0/1 response, stacked as 1/2,
    # whose cutoff is then passed to NNS.reg.
    set.seed(21)
    n <- 400
    x <- data.frame(matrix(rnorm(n * 6), n, dimnames = list(NULL, paste0("v", 1:6))))
    y <- rbinom(n, 1, plogis(1.5 * x$v1 - x$v2 + .5 * x$v3))
    s <- NNS.stack(IVs.train = x, DV.train = y + 1, type = "CLASS", status = FALSE,
                   ncores = 1, folds = 3)
    expect_length(s$NNS.dim.red.threshold, 1)
    expect_true(is.finite(s$NNS.dim.red.threshold))
    expect_lt(s$NNS.dim.red.count, 6)
    r <- NNS.reg(x, y, dim.red.method = "cor", threshold = s$NNS.dim.red.threshold,
                 plot = FALSE)
    expect_equal(unname(kept(r$equation)), rep(s$NNS.dim.red.count, 2))
  })
  aud_case("with more classes, the stack reports one cutoff per class", env = dimred_threshold, {
    set.seed(9)
    n <- 600
    x <- data.frame(matrix(rnorm(n * 6), n))
    cls <- factor(max.col(cbind(x$X1, x$X2, -x$X1 - x$X2 + rnorm(n, sd = .5))))
    s <- NNS.stack(IVs.train = x, DV.train = cls, type = "CLASS", status = FALSE,
                   ncores = 1, folds = 5)
    expect_named(s$NNS.dim.red.threshold, levels(cls))
    expect_lt(s$NNS.dim.red.count, 6)
    r <- NNS.reg(x, cls, dim.red.method = "cor", threshold = s$NNS.dim.red.threshold,
                 plot = FALSE)
    expect_equal(unname(kept(r$equation)), rep(s$NNS.dim.red.count, 3))
  })
  aud_case("the logistic regression example runs with a renamed test column", env = dimred_threshold, {
    # read.csv of an Excel "CSV UTF-8" test file names its first column X...Male.
    set.seed(23)
    n <- 300
    x <- data.frame(Male = rbinom(n, 1, .4), Pregnancy.Test = rbinom(n, 1, .2),
                    Folic.Acid = rbinom(n, 1, .2), Wine = rbinom(n, 1, .3),
                    Ginger.Ale = rbinom(n, 1, .2))
    y <- rbinom(n, 1, plogis(-1 + 2 * x$Pregnancy.Test + 1.5 * x$Folic.Acid - x$Wine))
    test <- setNames(x[1:40, ], c("X...Male", names(x)[-1]))
    s <- NNS.stack(IVs.train = x, DV.train = y + 1, IVs.test = test, balance = TRUE,
                   order = "max", type = "CLASS", status = FALSE, ncores = 1, folds = 3)
    expect_length(s$stack, 40)
    expect_true(is.finite(s$NNS.dim.red.threshold))
    r <- NNS.reg(x = x, y = y, point.est = test, plot = FALSE, dim.red.method = "cor",
                 threshold = s$NNS.dim.red.threshold)
    expect_length(r$Point.est, 40)
    expect_true(all(kept(r$equation) >= s$NNS.dim.red.count))
  })
})

test_that("CLS05: NNS.reg thresholds: per class by name or order, NA and NULL defaults, the all-cells neighbour count", {
  withr::local_seed(20261108)
  dimred_threshold <- fx_dimred_threshold()
  aud_case("NNS.reg applies a per-class threshold by class name or in class order", env = dimred_threshold, {
    # Each cutoff keeps two different predictors for its class; 0.8 applied to
    # versicolor would keep none.
    th <- c(setosa = .8, versicolor = .1, virginica = .7)
    fit <- function(threshold) {
      NNS.reg(iris[, 1:4], iris[, 5], dim.red.method = "cor", threshold = threshold,
              plot = FALSE)$equation
    }
    in.order <- fit(unname(th))
    expect_equal(fit(th[c(3, 1, 2)]), in.order)
    for (k in names(th)) expect_equal(in.order[[k]], fit(th[[k]])[[k]])
    expect_equal(unname(kept(in.order)), c(2L, 2L, 2L))
    expect_error(fit(c(.1, .2)), "one cutoff per class")
    expect_error(
      NNS.reg(iris[, 2:4], iris[, 1], dim.red.method = "cor", threshold = c(.1, .2),
              plot = FALSE),
      "single number for regression"
    )
  })
  aud_case("NA thresholds and neighbour counts select the defaults", {
    x <- iris[, 2:4]
    y <- iris[, 1]
    expect_identical(
      NNS.reg(x, y, dim.red.method = "cor", threshold = NA, plot = FALSE)$equation,
      NNS.reg(x, y, dim.red.method = "cor", threshold = 0, plot = FALSE)$equation
    )
    expect_identical(
      NNS.reg(iris[, 1:4], iris[, 5], dim.red.method = "cor", threshold = NULL,
              plot = FALSE)$equation,
      NNS.reg(iris[, 1:4], iris[, 5], dim.red.method = "cor", threshold = 0,
              plot = FALSE)$equation
    )
    expect_identical(
      NNS.reg(x, y, n.best = NA, point.est = x[1:5, ], plot = FALSE, ncores = 1)$Point.est,
      NNS.reg(x, y, point.est = x[1:5, ], plot = FALSE, ncores = 1)$Point.est
    )
    # Without a dimension-reduction component the stack reports NA, which
    # NNS.reg accepts.
    s <- NNS.stack(IVs.train = x, DV.train = y, method = 1, status = FALSE,
                   ncores = 1, folds = 2, seed = 1)
    expect_true(is.na(s$NNS.dim.red.threshold))
    expect_error(
      NNS.reg(x, y, dim.red.method = "cor", threshold = s$NNS.dim.red.threshold,
              plot = FALSE),
      NA
    )
  })
  aud_case("a classification stack reports the cell count when all cells win", {
    # A single predictor's only neighbour candidate is all cells.
    s <- suppressWarnings(NNS.stack(IVs.train = iris[, 3], DV.train = iris[, 5],
                                    type = "CLASS", status = FALSE, ncores = 1,
                                    folds = 3))
    expect_true(is.finite(s$NNS.reg.n.best))
    expect_gte(s$NNS.reg.n.best, 1)
  })
})
