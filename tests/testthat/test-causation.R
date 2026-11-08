# NNS.caus and Uni.caus: input contract, the [0, 1] index, lags, null
# calibration, p-values and the causation matrix.
#
# NNS.caus follows book chapter 14 (#68): explicitly lagged pairs on the copula
# scale, lagged co-partial-moment dependence, and a [0, 1] index with 0.5 neutral.
#
# NNS.caus after #96: aligned dim.red weights, a null centred for any margins,
# one lag per permutation test, and circular shifts for time series.

test_that("CAU01: NNS.caus input contract: frames, factors, shapes, lengths, finiteness, constants and tau", {
  withr::local_seed(20261108)
  aud_case("NNS.caus treats single-column data frames like vectors", {
    x <- as.double(1:20); y <- x^2
    v <- NNS.caus(x, y)
    expect_equal(NNS.caus(data.frame(x), y), v)
    expect_equal(NNS.caus(x, data.frame(y)), v)
    expect_equal(NNS.caus(data.frame(x), data.frame(y)), v)
    expect_error(NNS.caus(data.frame(x = replace(x, 2, Inf)), y), "must be finite")
  })
  aud_case("NNS.caus keeps a factor in a single-column data frame", {
    f <- factor(rep(c("a", "b"), 30))
    y <- as.double(1:60)^2
    v <- NNS.caus(f, y, factor.2.dummy = TRUE)
    expect_equal(unname(v[["C(x--->y)"]] + v[["C(y--->x)"]]), 1)
    expect_equal(NNS.caus(data.frame(f), y, factor.2.dummy = TRUE), v)
  })
  aud_case("NNS.caus rejects multi-column x or y when y is supplied", {
    set.seed(1); x <- rnorm(60); y <- x^2 + rnorm(60, sd = 0.1)
    expect_error(NNS.caus(data.frame(x, x), y), "single variable")
    expect_error(NNS.caus(x, cbind(y, y)), "single variable")
    expect_error(NNS.caus(factor(rep(1:3, 20)), y, factor.2.dummy = TRUE), "single variable")
    expect_length(NNS.caus(factor(rep(1:2, 30)), y, factor.2.dummy = TRUE), 2)
  })
  aud_case("NNS.caus rejects unequal lengths", {
    set.seed(1); x <- rnorm(60); y <- rnorm(60)
    expect_error(NNS.caus(x, y[-1]), "same length")
  })
  aud_case("NNS.caus rejects infinite values", {
    set.seed(1); x <- rnorm(60); y <- rnorm(60)
    expect_error(NNS.caus(replace(x, 4, Inf), y), "must be finite")
    expect_error(NNS.caus(x, replace(y, 4, -Inf)), "must be finite")
  })
  aud_case("NNS.caus rejects constant inputs", {
    set.seed(1); x <- rnorm(60)
    expect_error(NNS.caus(x, rep(1, 60)), "must not be constant")
    expect_error(NNS.caus(rep(1, 60), x, tau = 3), "must not be constant")
  })
  aud_case("NNS.caus rejects tau that is not \"cs\", \"ts\" or a non-negative integer", {
    set.seed(1); x <- rnorm(60); y <- x^2 + rnorm(60, sd = 0.1)
    for (tau in list(1.5, -1, "3", "foo", c(1, 2), NA, TRUE))
      expect_error(NNS.caus(x, y, tau = tau), "tau must be")
    expect_equal(NNS.caus(x, y, tau = 3L), NNS.caus(x, y, tau = 3))
    expect_equal(NNS.caus(x, y, tau = "cs"), NNS.caus(x, y, tau = 0))
  })
  aud_case("NNS.caus rejects tau that leaves fewer than three lagged rows", {
    set.seed(1); x <- rnorm(20); y <- rnorm(20)
    expect_length(NNS.caus(x, y, tau = 17), 2)
    expect_error(NNS.caus(x, y, tau = 18), "tau must be")
    expect_error(NNS.caus(c(1, 3), c(2, 1)), "tau must be")
  })
  aud_case("NNS.caus tau = \"ts\" stops when no seasonal period qualifies", {
    ar1 <- function(seed) {
      set.seed(seed); e <- rnorm(50); a <- numeric(50)
      for (i in 2:50) a[i] <- 0.9 * a[i - 1] + e[i]
      list(a = a, b = rnorm(50))
    }
    # NNS.seas reports only lags with >= 5 points per component (#69), so an
    # AR(1) usually returns Period 1 (no seasonality), which qualifies as
    # tau = 1. For seed 6 every period NNS.seas reports for x exceeds sqrt(50).
    s <- ar1(6)
    expect_true(all(NNS.seas(s$a, plot = FALSE)$periods > sqrt(50)))
    expect_error(NNS.caus(s$a, s$b, tau = "ts"), "no seasonal period")
    s <- ar1(8)
    expect_length(NNS.caus(s$a, s$b, tau = "ts"), 2)
  })
})

test_that("CAU02: the index is in [0, 1], sums to 1, is neutral for x = y, swaps with the arguments and is invariant to monotone maps", {
  withr::local_seed(20261108)
  aud_case("the index is in [0, 1], sums to 1 and is neutral for x = y", {
    set.seed(1); x <- rnorm(300); y <- x^2 + rnorm(300, sd = 0.2)
    v <- NNS.caus(x, y)
    expect_named(v, c("C(x--->y)", "C(y--->x)"))
    expect_true(all(v >= 0 & v <= 1))
    expect_equal(unname(sum(v)), 1)
    expect_equal(unname(NNS.caus(x, x)), c(0.5, 0.5))
    expect_equal(unname(NNS.caus(x, x, tau = 2)), c(0.5, 0.5))
  })
  aud_case("swapping the arguments swaps the directions", {
    set.seed(2); x <- rnorm(300); y <- sin(2 * x) + rnorm(300, sd = 0.2)
    expect_equal(unname(NNS.caus(y, x)), rev(unname(NNS.caus(x, y))))
    expect_equal(unname(NNS.caus(y, x, tau = 1)), rev(unname(NNS.caus(x, y, tau = 1))))
  })
  aud_case("a non-invertible function is attributed to its argument", {
    set.seed(3); x <- rnorm(1000)
    expect_gt(NNS.caus(x, x^2)[["C(x--->y)"]], 0.9)
    expect_lt(NNS.caus(x^2, x)[["C(x--->y)"]], 0.1)
  })
  aud_case("an invertible relationship has no cross-sectional direction", {
    set.seed(4); x <- rnorm(500)
    expect_equal(unname(NNS.caus(x, 2 * x + 1)), c(0.5, 0.5))
  })
  aud_case("the result is invariant to shifts, rescaling and increasing transforms", {
    set.seed(5); x <- rnorm(500, 1); y <- x - 1.5 + rnorm(500, sd = 0.3)
    v <- NNS.caus(x, y)
    expect_equal(NNS.caus(x, y + 10), v)
    expect_equal(NNS.caus(3 * x - 7, exp(y)), v)
    expect_equal(NNS.caus(x, y, tau = 2), NNS.caus(x + 100, 0.01 * y, tau = 2))
  })
})

test_that("CAU03: lags: tau identifies lagged causation, tau = \"ts\" resolves once, dim.red with ts.test, weights drive dimension reduction", {
  withr::local_seed(20261108)
  aud_case("tau lags x against y and identifies lagged causation", {
    set.seed(6)
    lagged <- replicate(10, {
      x <- rnorm(400); y <- c(0, 0.8 * x[-400]) + 0.3 * rnorm(400)
      NNS.caus(x, y, tau = 1)[["C(x--->y)"]]
    })
    expect_gt(min(lagged), 0.7)
    # A contemporaneous comparison cannot see a pure lag-1 effect.
    set.seed(6); x <- rnorm(400); y <- c(0, 0.8 * x[-400]) + 0.3 * rnorm(400)
    expect_lt(abs(NNS.caus(x, y, tau = 0)[["C(x--->y)"]] - 0.5), 0.2)
  })
  aud_case("tau = \"ts\" with p.value resolves the lag once", {
    set.seed(1); x <- as.numeric(arima.sim(list(ar = .9), 50)); y <- rnorm(50)
    r <- NNS.caus(x, y, tau = "ts", p.value = TRUE, nperm = 20, seed = 1)
    expect_equal(r$causation, NNS.caus(x, y, tau = "ts"))
    expect_identical(r$p.value$permute, "shift")
    expect_length(r$p.value$null_distribution, 20)
  })
  aud_case("dim.red with ts.test weights the aligned pairs the model fits", {
    set.seed(9)
    X <- cbind(signal = rnorm(200), noise = rnorm(200))
    y <- X[, 1]^2 + rnorm(200, sd = .1)
    fit <- NNS.stack(X[1:160, ], y[1:160], IVs.test = X[161:200, ], ts.test = 32,
                     folds = 1, method = 2, stack = FALSE, dim.red.method = "NNS.caus",
                     order = 5, optimize.threshold = FALSE, ncores = 1, status = FALSE)
    w <- fit$NNS.dim.red.coefficients
    expect_gt(w[1], 4 * w[2])                               # 0 / 0.0186 before
    expect_lt(sqrt(mean((fit$dim.red - y[161:200])^2)), 0.5)  # 1.618 before
    # cs and ts give the same weight: both score aligned pairs.
    expect_identical(NNS:::Uni.caus(y, X[, 1], "ts"), NNS:::Uni.caus(y, X[, 1], "cs"))
  })
  aud_case("NNS.caus weights still drive dimension reduction", {
    set.seed(9); X <- cbind(x1 = rnorm(200), x2 = rnorm(200)); y <- X[, 1]^2 + rnorm(200, sd = 0.1)
    fit <- NNS.reg(X, y, dim.red.method = "NNS.caus", plot = FALSE, ncores = 1)
    expect_true(all(is.finite(fit$Fitted.xy$y.hat)))
    w <- fit$equation$Coefficient
    expect_gt(abs(w[1]), abs(w[2]))
  })
})

test_that("CAU04: null calibration: independent series centre at 0.5 whatever their ties, and P is unbiased", {
  withr::local_seed(20261108)
  aud_case("independent series stay near the neutral 0.5", {
    set.seed(7)
    v <- replicate(20, NNS.caus(rnorm(500), rnorm(500), tau = 1)[["C(x--->y)"]])
    expect_lt(abs(mean(v) - 0.5), 0.1)
    expect_lt(max(abs(v - 0.5)), 0.4)
  })
  aud_case("independent series with different tie structures centre at 0.5", {
    set.seed(71)
    v <- replicate(20, NNS.caus(rbinom(1000, 1, .05), rnorm(1000))[[1]])
    expect_lt(abs(mean(v) - 0.5), 0.05)                     # 0.408 before
    set.seed(72)
    w <- replicate(20, NNS.caus(rexp(1000), sample(1:5, 1000, TRUE))[[1]])
    expect_lt(abs(mean(w) - 0.5), 0.05)
  })
  aud_case("P is unbiased under independence for a binary effect", {
    # The degree-2 deviation with the m / (m - 1) factor is a U-statistic, so
    # the within-cell average is unbiased whatever the effect's distribution.
    set.seed(73)
    p <- replicate(200, {
      cause <- NNS:::.nns_caus_u(rnorm(500)); effect <- NNS:::.nns_caus_u(rbinom(500, 1, .1))
      n <- length(cause); k <- floor(sqrt(n)); cell <- pmax(1L, ceiling(cause * k))
      d <- vapply(split(effect, cell), NNS:::.nns_caus_dev, numeric(1)); w <- lengths(split(effect, cell))
      1 - sum(w * d) / sum(w) / NNS:::.nns_caus_dev(effect)
    })
    expect_lt(abs(mean(p)), 0.01)
    v <- c(3, 1, 4, 1, 5)
    expect_equal(NNS:::.nns_caus_dev(v), var(v))
    # About its centre, an unresolved sign is not resolved: knowing y = x^2
    # leaves x's deviation in place.
    set.seed(74); x <- rnorm(1000)
    r <- NNS:::.nns_caus_raw(x, x^2, 0L)
    expect_gt(r$P[1], 0.8)
    expect_lt(r$P[2], 0.1)
  })
})

test_that("CAU05: p.value: x -> y orientation, circular-shift null for time series, shuffle for cross-sections, each offset once, rows of a one-column frame", {
  withr::local_seed(20261108)
  aud_case("NNS.caus p.value statistic is oriented x -> y", {
    set.seed(1); x <- rnorm(100); y <- x^2 + rnorm(100, sd = 0.1)
    xy <- NNS.caus(x, y, p.value = TRUE, nperm = 10, seed = 2)
    yx <- NNS.caus(y, x, p.value = TRUE, nperm = 10, seed = 2)
    expect_identical(names(yx$causation), c("C(x--->y)", "C(y--->x)"))
    # The statistic is C(x--->y) - C(y--->x): y = x^2 is driven by x.
    expect_gt(xy$p.value$observed_signed, 0.5)
    expect_equal(yx$p.value$observed_signed, -xy$p.value$observed_signed)
  })
  aud_case("time series use a circular-shift null, cross-sections a shuffle", {
    set.seed(2); x <- rnorm(100); y <- c(0, x[-100]) + rnorm(100)
    expect_identical(NNS.caus(x, y, tau = 1, p.value = TRUE, nperm = 5, seed = 1)$p.value$permute, "shift")
    expect_identical(NNS.caus(x, y, p.value = TRUE, nperm = 5, seed = 1)$p.value$permute, "y")
    # A shift is a rotation at least max(tau + 1, n / 10) steps from alignment.
    expect_identical(NNS:::.nns_caus_shift_offsets(100, 2L), 10:90)
    expect_identical(NNS:::.nns_caus_shift(1:100, 30L), c(31:100, 1:30))
  })
  aud_case("the shift null uses each admissible offset at most once", {
    # n = 100, tau = 49 admits only offset 50. Drawn 19 times with replacement it
    # gave p = 1 / 20 whenever the observed statistic beat it.
    p <- vapply(1:10, function(s) {
      set.seed(s); x <- rnorm(100); y <- rnorm(100)
      NNS.caus(x, y, tau = 49, p.value = TRUE, nperm = 19, seed = s)$p.value$two.sided
    }, numeric(1))
    expect_true(all(p >= 0.5))
    set.seed(1); M <- matrix(rnorm(120), 60, 2)
    r <- NNS.caus(M, tau = 29, p.value = TRUE, nperm = 19, seed = 1)
    expect_gte(r$p.value_matrix[1, 2], 0.5)
    set.seed(2); d <- NNS:::.nns_caus_shift_draw(100, 2L, 19)
    expect_length(d, 19)
    expect_false(anyDuplicated(d) > 0)
    expect_true(all(d %in% 10:90))
  })
  aud_case("p.value permutes the rows of a single-column data frame", {
    set.seed(1); x <- rnorm(80); y <- x^2 + rnorm(80, sd = 0.2)
    v <- NNS.caus(x, y, p.value = TRUE, nperm = 20, seed = 3)
    d <- NNS.caus(x, data.frame(y = y), p.value = TRUE, nperm = 20, seed = 3)
    expect_equal(d$p.value$null_distribution, v$p.value$null_distribution)
    expect_gt(length(unique(d$p.value$null_distribution)), 1)
    expect_equal(NNS.caus(data.frame(x = x), y, p.value = TRUE, nperm = 20, seed = 3, permute = "x")$p.value$two.sided,
                 NNS.caus(x, y, p.value = TRUE, nperm = 20, seed = 3, permute = "x")$p.value$two.sided)
  })
})

test_that("CAU06: causation matrix: C(row -> column) with a 0 diagonal; its p-values keep factors numeric, shift one series per pair, honour permute, use the signed statistic", {
  withr::local_seed(20261108)
  aud_case("the causation matrix holds C(row ---> column) with a 0 diagonal", {
    set.seed(8); a <- rnorm(300); b <- a^2 + rnorm(300, sd = 0.1); cc <- rnorm(300)
    m <- NNS.caus(cbind(a, b, cc))
    expect_equal(unname(diag(m)), c(0, 0, 0))
    expect_equal(m["a", "b"], NNS.caus(a, b)[["C(x--->y)"]])
    expect_equal(m["b", "a"], NNS.caus(a, b)[["C(y--->x)"]])
    expect_equal(m["a", "b"] + m["b", "a"], 1)
    expect_gt(m["a", "b"], 0.8)
  })
  aud_case("the matrix p-value path keeps factor columns numeric", {
    r <- NNS.caus(iris, p.value = TRUE, nperm = 2, seed = 1)
    expect_equal(dim(r$causality), c(5, 5))
    expect_true(all(is.finite(r$p.value_matrix)))
  })
  aud_case("the matrix time-series null shifts one series per pair", {
    # Shifting every column independently can leave a pair aligned when the two
    # shifts nearly coincide; each pair now keeps its first column and shifts
    # the second, so the shift helper runs once per pair and permutation.
    calls <- 0L
    testthat::local_mocked_bindings(
      .nns_caus_shift = function(v, tau = 0L) { calls <<- calls + 1L; n <- length(v); c(v[31:n], v[1:30]) },
      .package = "NNS")
    set.seed(5); M <- matrix(rnorm(240), 60, 4)
    r <- NNS.caus(M, tau = 1, p.value = TRUE, nperm = 2, seed = 1)
    expect_identical(calls, 2L * ncol(utils::combn(4, 2)))   # 12; 8 with column shifts
    expect_true(all(is.finite(r$p.value_matrix)))
  })
  aud_case("the matrix p-value honours permute like the vector path", {
    set.seed(107); M <- matrix(rnorm(120), 60, 2)
    for (tau in c(0, 1)) for (m in c("x", "y", "both", "shift")) {
      mat <- NNS.caus(M, tau = tau, p.value = TRUE, nperm = 19, permute = m, seed = 7)$p.value_matrix[1, 2]
      vec <- NNS.caus(M[, 1], M[, 2], tau = tau, p.value = TRUE, nperm = 19, permute = m, seed = 7)$p.value$two.sided
      expect_equal(mat, vec, info = paste(tau, m))
    }
  })
  aud_case("the matrix p-value uses the same signed statistic as the vector path on tied data", {
    set.seed(17); invisible(rnorm(120))
    M <- cbind(sample(0:1, 60, TRUE), sample(0:1, 60, TRUE))
    mat <- NNS.caus(M, p.value = TRUE, nperm = 19, permute = "both", seed = 7)$p.value_matrix[1, 2]
    vec <- NNS.caus(M[, 1], M[, 2], p.value = TRUE, nperm = 19, permute = "both", seed = 7)$p.value$two.sided
    expect_equal(mat, vec)
  })
})
