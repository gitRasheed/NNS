# NNS.dep and NNS.copula: p-values, independence, range, cells and levels.
#
# A product grid is exactly independent, so the copula is 0.
#
# NNS.dep measures each cell's concordance against the within-cell
# permutation null, so independence no longer leaves a floor (#70).
#
# NNS.dep reads the debiased cell dependence at every level of the recursive
# mean split, on x, on y and on the joint quadrants, and reports the strongest
# level (#95).

test_that("DEP01: NNS.dep p.value: orientation, graphics, upper-tail and two-sided permutation p-values, ties, asym and infinite values", {
  withr::local_seed(20261108)
  aud_case("NNS.dep p.value without print.map opens no graphics device", {
    graphics.off()
    set.seed(11); x <- rnorm(50); y <- rnorm(50)
    expect_silent(NNS.dep(x, y, p.value = TRUE))
    expect_null(dev.list())
  })
  aud_case("NNS.dep Dependence p.value is the upper-tail permutation p-value", {
    set.seed(5); x <- rnorm(50); y <- x + rnorm(50)
    obs <- NNS.dep(x, y)$Dependence
    set.seed(6); perms <- replicate(100L, sample.int(50))
    null <- apply(perms, 2, function(i) NNS.dep(x, y[i])$Dependence)
    set.seed(6); p <- NNS.dep(x, y, p.value = TRUE)
    expect_equal(p$`Dependence p.value`, (1 + sum(null >= obs)) / 101)
    # 7 / 101 before the per-cell debiasing (#70) removed the independence floor.
    expect_equal(p$`Dependence p.value`, 1 / 101)
  })
  aud_case("NNS.dep Dependence p.value is at least 1 / 101", {
    set.seed(7); x <- rnorm(50)
    expect_equal(NNS.dep(x, x, p.value = TRUE)$`Dependence p.value`, 1 / 101)
  })
  aud_case("NNS.dep Correlation p.value is the two-sided permutation p-value", {
    set.seed(5); x <- rnorm(50); y <- x + rnorm(50)
    obs <- NNS.dep(x, y)$Correlation
    set.seed(6); perms <- replicate(100L, sample.int(50))
    null <- apply(perms, 2, function(i) NNS.dep(x, y[i])$Correlation)
    set.seed(6); p <- NNS.dep(x, y, p.value = TRUE)
    expect_equal(p$`Correlation p.value`,
                 min(1, 2 * min(1 + sum(null <= obs), 1 + sum(null >= obs)) / 101))
    expect_equal(p$`Correlation p.value`, 2 / 101)
  })
  aud_case("NNS.dep Correlation p.value counts tied null correlations", {
    # 99 of the 100 positions of y's event give Correlation 0; it was 0 before.
    x <- c(1, rep(0, 99)); y <- c(0, 1, rep(0, 98))
    set.seed(871); p <- NNS.dep(x, y, p.value = TRUE)
    expect_equal(p$Correlation, 0)
    expect_equal(p$`Correlation p.value`, 1)
  })
  aud_case("NNS.dep p.value = TRUE honours asym", {
    set.seed(3); x <- rnorm(100); y <- x^2 + rnorm(100, sd = 0.1)
    direct <- NNS.dep(y, x, asym = TRUE)
    set.seed(4); p <- NNS.dep(y, x, p.value = TRUE, asym = TRUE)
    expect_equal(p$Correlation, direct$Correlation)
    expect_equal(p$Dependence, direct$Dependence)
    # x is not a function of y = x^2: 0.6005389 before the per-cell debiasing
    # (#70), 0.1027830 before the multi-level measure (#95), 0.1156500 before
    # the gravity splits were corrected (#117).
    expect_equal(p$Dependence, 0.1045294, tolerance = 1e-6)
  })
  aud_case("NNS.dep rejects infinite values", {
    set.seed(8); x <- rnorm(30); y <- x + rnorm(30)
    expect_error(NNS.dep(replace(x, 3, Inf), y), "x must be finite")
    expect_error(NNS.dep(replace(x, 3, -Inf), y), "x must be finite")
    expect_error(NNS.dep(x, replace(y, 5, Inf)), "y must be finite")
    expect_error(NNS.dep(cbind(replace(x, 3, Inf), y)), "x must be finite")
    expect_error(NNS.dep(replace(x, 3, NaN), y), "missing values")
  })
})

test_that("DEP02: NNS.copula is 0 under independence: skewed margins, non-mean and boundary targets, margins on one side", {
  withr::local_seed(20261108)
  copula_marginal_anchor <- fx_copula_marginal_anchor()
  aud_case("NNS.copula is 0 on independent skewed margins at the default target", env = copula_marginal_anchor, {
    X <- grid_X()
    expect_equal(NNS.copula(X), 0, tolerance = 1e-6)
    expect_equal(NNS.copula(X, continuous = FALSE), 0, tolerance = 1e-6)
    expect_equal(NNS.copula(X[, 1:2]), 0, tolerance = 1e-6)
    x <- X[, 1]
    expect_equal(NNS.copula(cbind(x, x, -x)), 1)
  })
  aud_case("NNS.copula is 0 on independent margins at a non-mean target", env = copula_marginal_anchor, {
    X <- grid_X()
    t <- c(0.25, 1.75, 0.42)
    expect_equal(NNS.copula(X, target = t), 0, tolerance = 1e-6)
    expect_equal(NNS.copula(X, target = t, continuous = FALSE), 0, tolerance = 1e-6)
    Y <- rbind(c(-10, 0.1, 0.1), c(0.1, -10, 0.1), c(0.1, 0.1, -10))
    expect_lte(NNS.copula(Y, target = rep(0, 3)), 1)
  })
  aud_case("NNS.copula is 0 on an independent grid at a boundary target", {
    expect_equal(NNS.copula(expand.grid(0:1, 0:1), target = c(0, 0)), 0)
    expect_equal(NNS.copula(expand.grid(0:1, 0:1), target = c(0, 1)), 0)
    expect_equal(NNS.copula(expand.grid(0:1, 0:1, 0:2), target = c(0, 0, 0)), 0)
  })
  aud_case("NNS.copula is 0 when every margin lies on one side of its target", {
    # Each pair is then trivially concordant or discordant; these returned 1.
    g <- expand.grid(x = 1:10, y = 1:10)
    expect_identical(NNS.copula(g, target = c(0, 0)), 0)
    expect_identical(NNS.copula(g, target = c(0, 0), continuous = FALSE), 0)
    expect_identical(NNS.copula(g, target = c(10, 10)), 0)
    expect_identical(NNS.copula(g, target = c(0, 11)), 0)
    expect_identical(NNS.copula(expand.grid(1:4, 1:4, 1:4), target = c(0, 0, 0)), 0)
    set.seed(8); z <- cbind(rexp(1000), rexp(1000))
    expect_identical(NNS.copula(z, target = c(0, 0)), 0)
    x <- 1:10
    expect_identical(NNS.copula(cbind(x, x), target = c(0, 0)), 0)
    b <- cbind(rep(c(0, 1), 50), rep(c(1, 0), 50))
    expect_identical(NNS.copula(b, target = c(1, 1)), 0)
    expect_identical(NNS.copula(b, target = c(1, 1), continuous = FALSE), 0)
  })
})

# pin: NNS-Python reference values
test_that("DEP03: NNS.copula drops constant margins, is 1 for perfect dependence (inside targets, mixed signs), is not 1 for an interior pair sum, and the reference value", {
  withr::local_seed(20261108)
  copula_ref <- fx_copula_ref()
  aud_case("NNS.copula drops constant margins", {
    expect_equal(NNS.copula(cbind(-2:2, 1)), 0)
    expect_equal(NNS.copula(cbind(-2:2, 1), continuous = FALSE), 0)
    x <- 1:5
    expect_equal(NNS.copula(cbind(x, x, 3)), 1)
  })
  aud_case("A constant margin leaves the copula of the varying margins unchanged", {
    # p = 0.5 in both margins, where the old fixed anchors were exact; values match the old cbind(x, y)
    x <- c(1, 2, 3, 4, 5, 6)
    y <- c(3, 1, 6, 2, 5, 4)
    expect_equal(NNS.copula(cbind(x, y, 5)), 0.690065559342, tolerance = 1e-10)
    expect_equal(NNS.copula(cbind(x, y, 5), continuous = FALSE), sqrt(1/3), tolerance = 1e-10)
    Z <- iris[, 1:3]
    expect_equal(NNS.copula(cbind(Z, k = 1)), NNS.copula(Z))
    expect_equal(NNS.copula(cbind(Z, k = 1), target = c(5, 3, 4, 0)), NNS.copula(Z, target = c(5, 3, 4)))
  })
  aud_case("NNS.copula is 1 for perfect dependence at targets inside the data", {
    x <- 1:10
    for (cont in c(TRUE, FALSE)) {
      expect_equal(NNS.copula(cbind(x, x), continuous = cont), 1)
      expect_equal(NNS.copula(cbind(x, -x), continuous = cont), 1)
      expect_equal(NNS.copula(cbind(x, 2 * x + 1), target = c(3, 7), continuous = cont), 1)
      expect_equal(NNS.copula(cbind(x, -x), target = c(3.5, -3.5), continuous = cont), 1)
      expect_equal(NNS.copula(cbind(x, x, -x), target = c(3.5, 3.5, -3.5), continuous = cont), 1)
    }
  })
  aud_case("NNS.copula does not return 1 for an interior pair sum of 1", {
    X <- cbind(c(5, 3, 4, 5, 3, 3, 2, 5, 2, 2),
               c(2, 5, 1, 1, 5, 5, 5, 2, 2, 1),
               c(1, 1, 4, 2, 4, 5, 2, 4, 4, 5))
    expect_equal(NNS.copula(X), 0.577972677175, tolerance = 1e-10)
  })
  aud_case("NNS.copula is 1 for perfect dependence with mixed signs", {
    x <- 1:100
    expect_equal(NNS.copula(cbind(x, x, -x)), 1)
    expect_equal(NNS.copula(cbind(x, x, -x), continuous = FALSE), 1)
    expect_equal(NNS.copula(cbind(x, x, x, -x)), 1)
  })
  aud_case("Copula", env = copula_ref, {
    expect_equal(B, 0.2542705, tolerance=1e-5)
    expect_equal(C, 0.1969773, tolerance=1e-5)
    expect_equal(D, 0.2506271, tolerance=1e-5)
    expect_equal(E, 0.2088075, tolerance=1e-5)
  })
})

test_that("DEP04: NNS.dep under independence: no floor, weak dependence separated, Correlation p-value size", {
  withr::local_seed(20261108)
  aud_case("NNS.dep has no independence floor", {
    set.seed(70)
    dep <- replicate(20, NNS.dep(rnorm(1000), rnorm(1000))$Dependence)
    # Before the debiasing the mean was about 0.32 at n = 1000.
    expect_lt(mean(dep), 0.06)
    set.seed(71)
    dep100 <- replicate(40, NNS.dep(rnorm(100), rnorm(100))$Dependence)
    # About 0.60 before.
    expect_lt(mean(dep100), 0.15)
  })
  aud_case("NNS.dep separates weak dependence from independence", {
    set.seed(74)
    pairs <- replicate(10, {
      x <- rnorm(1000)
      c(NNS.dep(x, x + rnorm(1000))$Dependence, NNS.dep(x, rnorm(1000))$Dependence)
    })
    expect_gt(mean(pairs[1, ]), mean(pairs[2, ]) + 0.1)
    # Nonlinear dependence is still detected.
    x <- rnorm(1000)
    expect_gt(NNS.dep(x, x^2 + rnorm(1000, sd = 0.1))$Dependence, 0.6)
  })
  aud_case("NNS.dep Correlation p.value holds its size under independence", {
    # 28 / 200 rejections at 0.05 when it was the smaller one-sided tail.
    set.seed(5)
    p <- replicate(200L, NNS.dep(rnorm(30), rnorm(30), p.value = TRUE)$`Correlation p.value`)
    expect_lte(mean(p <= 0.05), 0.10)
  })
})

test_that("DEP05: NNS.dep range: Dependence in [0, 1] bounds |Correlation|; +-1 for x, -x, binary, shifted binary and decimal-discrete pairs", {
  withr::local_seed(20261108)
  aud_case("NNS.dep keeps perfect dependence at 1", {
    set.seed(72); x <- rnorm(500)
    expect_equal(unname(unlist(NNS.dep(x, x))), c(1, 1))
    expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1))
  })
  aud_case("NNS.dep Dependence stays in [0, 1] and bounds |Correlation|", {
    set.seed(73)
    for (i in 1:20) {
      n <- sample(c(20, 60, 200), 1)
      x <- rnorm(n); y <- 0.5 * x + rnorm(n)
      for (asym in c(FALSE, TRUE)) {
        r <- NNS.dep(x, y, asym = asym)
        expect_gte(r$Dependence, 0)
        expect_lte(r$Dependence, 1)
        expect_lte(abs(r$Correlation), r$Dependence + 1e-12)
      }
    }
  })
  aud_case("NNS.dep is 1 for identical binary data at any group size", {
    b <- rep(0:1, c(51, 49))
    expect_equal(unname(unlist(NNS.dep(b, b))), c(1, 1))
  })
  aud_case("NNS.dep keeps the sign of a fully discordant cell", {
    x <- c(-2, -1, 1, 2)
    expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1))
    expect_equal(unname(unlist(NNS.dep(x, -x, asym = TRUE))), c(-1, 1))
    expect_equal(NNS.dep(cbind(x, -x))$Correlation[1, 2], -1)
  })
  aud_case("a shifted binary pair keeps perfect dependence", {
    b <- rep(0:1, each = 50)
    expect_equal(unlist(NNS.dep(b + .1, b + .1)), c(Correlation = 1, Dependence = 1))
  })
  aud_case("y determined by a decimal-valued discrete x has dependence 1", {
    x <- rep(c(.1, .2, .3, .7), 25)
    expect_equal(NNS.dep(x, x^2)$Dependence, 1)
  })
})

test_that("DEP06: NNS.dep cells: constant cells, ties with the cell mean, tied grids, argument and row order, the order heuristic on tied data", {
  withr::local_seed(20261108)
  aud_case("a cell with constant x counts only when y is determined", {
    set.seed(75)
    xb <- rbinom(500, 1, 0.5); yb <- rbinom(500, 1, 0.5)
    expect_lt(NNS.dep(xb, yb)$Dependence, 0.1)
    expect_equal(NNS.dep(xb, xb)$Dependence, 1)
  })
  aud_case("NNS.dep degree-0 discordant count excludes ties with the cell mean", {
    # The middle point ties both cell means, so it counts as co-lower and the
    # cell is not all-discordant. Three points then carry no evidence beyond
    # the within-cell permutation null (#70), so the result is 0, not -1.
    x <- -1:1
    expect_equal(unname(unlist(NNS.dep(x, -x))), c(0, 0))
    # With five points the countermonotone signal clears the null: the debiased
    # degree-0 and degree-1 terms average to 2/3 (no square root since #95).
    x <- -2:2
    expect_equal(unname(unlist(NNS.dep(x, -x))), c(-1, 1) * 2 / 3)
  })
  aud_case("NNS.dep reads an exactly independent tied grid as 0", {
    # The discrete blend used to add 1/6 here (#95).
    G <- as.matrix(expand.grid(-1:1, -1:1))[rep(1:9, each = 100), ]
    expect_equal(unname(unlist(NNS.dep(G[, 1], G[, 2]))), c(0, 0))
    expect_equal(unname(unlist(NNS.dep(G[, 1], G[, 2], asym = TRUE))), c(0, 0))
  })
  aud_case("symmetric NNS.dep on tied data ignores argument and row order", {
    # Zero covariances and values equal to a cell mean rounded either way
    # depending on the order the cell was summed in, more so at a high level.
    set.seed(487); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
    expect_equal(unlist(NNS.dep(y, x)), unlist(NNS.dep(x, y)), tolerance = 1e-12)
    for (lev in c(0, 1e9)) for (div in c(1, 10)) {
      d <- vapply(1:300, function(s) {
        set.seed(s); n <- sample(10:80, 1); x <- lev + sample(0:5, n, TRUE) / div; y <- lev + sample(0:5, n, TRUE) / div
        o <- sample(n)
        a <- unlist(NNS.dep(x, y))
        max(abs(a - unlist(NNS.dep(y, x))), abs(a - unlist(NNS.dep(x[o], y[o]))))
      }, numeric(1))
      expect_lt(max(d), 1e-12)
    }
    # Decimal ties give the result of their exact integer image.
    set.seed(208); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
    expect_equal(NNS.dep(x / 10, y / 10), NNS.dep(x, y), tolerance = 1e-12)   # 0.167 vs -0.189 before
    set.seed(78); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
    # 1e9 + 0.1 is stored to 1.2e-7, so the shifted data matches its image to about 1e-6.
    expect_equal(NNS.dep(1e9 + x / 10, 1e9 + y / 10), NNS.dep(x, y), tolerance = 1e-6)   # 0.110 vs 0.191
    # Distinct values at a high level are not ties.
    ts <- 1.7e9 + seq_len(100) / 1e6
    expect_equal(unname(unlist(NNS.dep(ts, -ts))), c(-1, 1))
    expect_equal(NNS.dep(iris[, 1], iris[, 2]), NNS.dep(iris[, 2], iris[, 1]), tolerance = 1e-12)
  })
  aud_case("the local dependence NNS.reg's default order reads is unchanged on tied data", {
    set.seed(388); n <- sample(10:80, 1); x <- sample(0:5, n, TRUE); y <- sample(0:5, n, TRUE)
    # Pinned again after #117 corrected the gravity splits (0.054, order 1 before).
    expect_equal(NNS:::.nns_dep_local(x, y), 0.281843356364165, tolerance = 1e-12)
    expect_identical(NNS:::.nns_reg_default_order(x, y), 2L)
    expect_equal(NNS.reg(x, y, point.est = 4, plot = FALSE, residual.plot = FALSE)$Point.est,
                 2.4009920634920636, tolerance = 1e-12)
  })
  aud_case("constant cells keep the direction of a discrete relationship", {
    b <- rep(0:1, each = 50)
    expect_equal(unlist(NNS.dep(b, -b)), c(Correlation = -1, Dependence = 1))
    expect_equal(NNS.dep(b + .1, -(b + .1))$Correlation, -1)
    u <- rep(1:3, 30)
    expect_equal(NNS.dep(u, (u - 2)^2)$Correlation, 0)
    g <- rep(1:3, each = 30)
    expect_equal(NNS.dep(g, (g - 2)^2)$Correlation, 0)
    expect_equal(NNS.dep(rev(g), rev((g - 2)^2))$Correlation, 0)
    x <- 1e6 + rep(c(-.01, 0, .01), each = 1000); y <- 1e6 + rep(c(.01, 0, .01), each = 1000)
    expect_equal(c(NNS.dep(x, y)$Correlation, NNS.dep(rev(x), rev(y))$Correlation), c(0, 0))
    expect_equal(NNS.dep(u, u)$Correlation, 1)
    expect_equal(NNS.dep(rep(c(0, 1, 2, 10), 30), rep(0:3, 30))$Correlation, 1)
  })
})

test_that("DEP07: NNS.dep levels: noisy predictors, non-functional relations, asym, discrete margins, exact null, native vs partition map, finest-cell order", {
  withr::local_seed(20261108)
  aud_case("a strong but noisy predictor is not scored 0", {
    # Before #95 the iris signal scored 0 and the noise column 0.065.
    set.seed(3)
    X <- cbind(signal = iris$Sepal.Length, noise = rnorm(150))
    y <- iris$Petal.Width
    expect_gt(NNS.dep(X[, 1], y, asym = TRUE)$Dependence, 0.75)
    expect_lt(NNS.dep(X[, 2], y, asym = TRUE)$Dependence, 0.3)
    f <- NNS.reg(X, y, dim.red.method = "NNS.dep", plot = FALSE, residual.plot = FALSE)
    expect_gt(f$equation$Coefficient[1], 4 * f$equation$Coefficient[2])
    expect_gt(f$R2, 0.6)

    # y = x + e (r = 0.71) ranged from 0.03 to 0.32 across seeds.
    dep <- vapply(1:5, function(s) {
      set.seed(s); x <- rnorm(1000)
      NNS.dep(x, x + rnorm(1000))$Dependence
    }, numeric(1))
    expect_true(all(dep > 0.55 & dep < 0.8))
  })
  aud_case("relations without a function are detected by the joint quadrants", {
    set.seed(123)
    df <- data.frame(x = runif(10000, -1, 1), y = runif(10000, -1, 1))
    df <- subset(df, x^2 + y^2 <= 1 & x^2 + y^2 >= 0.95)
    circle <- NNS.dep(df$x, df$y)
    expect_gt(circle$Dependence, 0.8)
    expect_lt(abs(circle$Correlation), 0.2)

    set.seed(5); x <- rnorm(1000)
    expect_gt(NNS.dep(x, x * rnorm(1000))$Dependence, 0.4)
  })
  aud_case("asym = TRUE measures the dependence of y on x only", {
    set.seed(6); x <- rnorm(500); y <- x^2 + rnorm(500, sd = 0.1)
    expect_gt(NNS.dep(x, y, asym = TRUE)$Dependence, 0.8)
    expect_lt(NNS.dep(y, x, asym = TRUE)$Dependence, 0.3)
    # The symmetric measure includes both conditional partitions.
    sym <- NNS.dep(x, y)$Dependence
    expect_gte(sym, NNS.dep(x, y, asym = TRUE)$Dependence)
    expect_gte(sym, NNS.dep(y, x, asym = TRUE)$Dependence)
  })
  aud_case("discrete margins: no blend floor, and |Correlation| <= Dependence", {
    g <- expand.grid(x = c(rep(0, 9), 1), y = c(rep(0, 9), 1))   # exactly independent
    expect_equal(unname(unlist(NNS.dep(g$x, g$y))), c(0, 0))
    set.seed(4); x <- sample(0:4, 100, TRUE); y <- c(0, 3, 0, 4, 3)[x + 1]
    for (asym in c(FALSE, TRUE)) {
      r <- NNS.dep(x, y, asym = asym)
      expect_equal(r$Dependence, 1)
      expect_lte(abs(r$Correlation), r$Dependence)
    }
    b <- rep(0:1, each = 50)
    expect_equal(unname(unlist(NNS.dep(b, 1 - b))), c(-1, 1))
  })
  aud_case("small tied cells use the exact permutation null", {
    # Seven points form a single cell. x has six tied values, so the null of the
    # degree-1 term is exact over the 7 distinct arrangements (the normal
    # approximation used for larger cells gives 0.659 here instead of 0.483).
    x <- c(rep(0, 6), 1); y <- c(1, 1, 1, 2, 2, 3, 9)
    dx <- x - mean(x); dy <- y - mean(y)
    t1 <- function(v) abs(sum(v * dy)) / sum(abs(v * dy))
    e1 <- mean(vapply(1:7, function(k) t1(replace(rep(dx[1], 7), k, dx[7])), numeric(1)))
    a <- sum(dx <= 0); b <- sum(dy <= 0); conc <- sum((dx <= 0) == (dy <= 0))
    k <- max(0, a + b - 7):min(a, b)
    e0 <- sum(dhyper(k, a, 7 - a, b) * abs(2 * (7 - a - b + 2 * k) / 7 - 1))
    t0 <- abs(2 * conc / 7 - 1)
    expected <- max(0, ((t0 - e0) / (1 - e0) + (t1(dx) - e1) / (1 - e1)) / 2)
    expect_equal(e1, 0.48295, tolerance = 1e-5)
    expect_equal(NNS.dep(x, y)$Dependence, expected, tolerance = 1e-12)
    expect_equal(NNS.dep(x, y, asym = TRUE)$Dependence, expected, tolerance = 1e-12)
  })
  aud_case("the native and partition-map paths agree", {
    set.seed(9); x <- rnorm(300); y <- x^2 + rnorm(300)
    native <- c(unlist(NNS.dep(x, y)), unlist(NNS.dep(x, y, asym = TRUE)))
    old <- options(NNS.native = FALSE); on.exit(options(old))
    expect_identical(c(unlist(NNS.dep(x, y)), unlist(NNS.dep(x, y, asym = TRUE))), native)
  })
  aud_case("NNS.reg's order heuristic keeps the local (finest-cell) dependence", {
    # The default order reads dependence as a fine-scale signal-to-noise
    # measure. The multi-level NNS.dep also counts dependence between coarse
    # cells (0.64 for y = x + e below) and would over-partition noisy data, so
    # the heuristic keeps the finest-cell measure. Its copula component now
    # averages both response orientations, preserving the sign-invariance
    # expected of a smoothing-complexity choice.
    set.seed(1); x <- rnorm(500); y <- x + rnorm(500)
    expect_equal(NNS:::.nns_dep_local(x, y), 0)
    expect_gt(NNS.dep(x, y, asym = TRUE)$Dependence, 0.55)
    expect_equal(NNS:::.nns_reg_default_order(x, y), NNS:::.nns_reg_default_order(x, -y))
    set.seed(2); x <- runif(400); y <- sin(6 * x) + rnorm(400, sd = .2)
    expect_equal(NNS:::.nns_dep_local(x, y), 0.394552637766, tolerance = 1e-10)
    expect_equal(NNS:::.nns_reg_default_order(x, y), NNS:::.nns_reg_default_order(x, -y))
    set.seed(3); u <- sample(1:4, 200, TRUE); v <- u^2 + rnorm(200, sd = .5)
    expect_equal(NNS:::.nns_dep_local(u, v, asym = FALSE), 0.0905745380689, tolerance = 1e-10)
  })
})
