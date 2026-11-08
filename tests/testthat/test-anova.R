# NNS.ANOVA: reference values, grand statistics, median ties, constant groups,
# missing values and the robust bootstrap.

test_that("AOV01: NNS.ANOVA reference values and the size-weighted grand mean and pooled median", {
  withr::local_seed(20261108)
  anova_ref <- fx_anova_ref()
  anova_grand_mean_medians <- fx_anova_grand_mean_medians()
  aud_case("NNS.ANOVA", env = anova_ref, {
    expect_equal(B, R1, tolerance=1e-4)
  })
  aud_case("NNS.ANOVA - pairwise", env = anova_ref, {
    expect_equal(C, R2, tolerance=1e-4)
  })
  aud_case("list and two-sample paths share the size-weighted grand mean", env = anova_grand_mean_medians, {
    two <- NNS.ANOVA(a, b, plot = FALSE, confidence.interval = NULL)
    expect_equal(two$Grand_Statistic, mean(c(a, b)))
    expect_equal(two$Certainty, 0.02320132, tolerance = 1e-6)
    expect_equal(NNS.ANOVA(list(a, b), plot = FALSE), c(Certainty = two$Certainty))
  })
  aud_case("grand median is the median of the pooled observations", env = anova_grand_mean_medians, {
    r <- NNS.ANOVA(x, y, medians = TRUE, means.only = TRUE, plot = FALSE,
                   confidence.interval = NULL)
    expect_identical(r$Grand_Statistic, 6)
  })
  aud_case("multi-group medians evaluate at the pooled median", env = anova_grand_mean_medians, {
    expected <- NNS:::NNS.ANOVA.bin(x, y, medians = TRUE, means.only = TRUE,
                                    mean.of.means = 6, plot = FALSE)$Certainty
    expect_equal(expected, 0.1789941, tolerance = 1e-6)
    expect_equal(NNS.ANOVA(list(x, y), medians = TRUE, means.only = TRUE, plot = FALSE),
                 c(Certainty = expected))
  })
})

test_that("AOV02: median agreement uses the whole CDF jump at a tied median, in multi-group and pairwise paths, and only there", {
  withr::local_seed(20261108)
  aud_case("median agreement uses the entire CDF jump at a tied median", {
    x <- c(rep(0, 80), rep(1, 20))
    fit <- NNS.ANOVA(x, x, means.only = TRUE, medians = TRUE,
                     confidence.interval = NULL, plot = FALSE)
    maximum <- ((length(x) * 2 - 2) / (length(x) * 2))^2
    expect_equal(fit$Grand_Statistic, 0)
    expect_equal(fit$Control_CDF, 0.8)
    expect_equal(fit$Treatment_CDF, 0.8)
    expect_equal(fit$Certainty, maximum)

    # A common median can have different jump heights in the two groups.
    y <- c(rep(0, 60), rep(2, 40))
    expect_equal(NNS.ANOVA(x, y, means.only = TRUE, medians = TRUE,
                           confidence.interval = NULL, plot = FALSE)$Certainty,
                 maximum)
    # The full distribution profile still detects differing tail structure.
    expect_lt(NNS.ANOVA(x, y, means.only = FALSE, medians = TRUE,
                        confidence.interval = NULL, plot = FALSE)$Certainty,
              maximum)
  })
  aud_case("median ties are handled in multi-group and pairwise paths", {
    x <- c(rep(0, 80), rep(1, 20))
    maximum <- (198 / 200)^2
    expect_equal(unname(NNS.ANOVA(list(x, x, x), means.only = TRUE,
                                  medians = TRUE, plot = FALSE)), maximum)
    pairwise <- NNS.ANOVA(cbind(x, x, x), means.only = TRUE,
                          medians = TRUE, pairwise = TRUE, plot = FALSE)
    expect_equal(as.numeric(pairwise[lower.tri(pairwise)]), rep(maximum, 3))
  })
  aud_case("median distance agrees with the original score away from atoms", {
    x <- c(-4, -2, 0, 2, 4)
    y <- c(1, 3, 5, 7, 9)
    target <- median(c(x, y))
    expect_false(target %in% c(x, y))
    distance <- max(abs(mean(x <= target) - 0.5), abs(mean(y <= target) - 0.5))
    expected <- (0.5 - distance)^2 / 0.25 * (8 / 10)^2
    expect_equal(NNS.ANOVA(x, y, means.only = TRUE, medians = TRUE,
                           confidence.interval = NULL, plot = FALSE)$Certainty,
                 expected)
    expect_equal(NNS.ANOVA(rep(0, 10), rep(1, 10), means.only = TRUE,
                           medians = TRUE, confidence.interval = NULL,
                           plot = FALSE)$Certainty, 0)
  })
  aud_case("mean agreement retains its degree-one target", {
    x <- c(-2, -1, 1, 2)
    expected <- (6 / 8)^2
    expect_equal(NNS.ANOVA(x, 4 * x, means.only = TRUE,
                           confidence.interval = NULL, plot = FALSE)$Certainty,
                 expected)
  })
})

test_that("AOV03: constant groups: the zero-area expectile balance, and agreement in mean but not in profile", {
  withr::local_seed(20261108)
  aud_case("constant groups use their zero-area expectile balance", {
    x <- rep(0.25, 10)
    maximum <- (18 / 20)^2
    for (central.only in c(TRUE, FALSE)) {
      for (use.medians in c(TRUE, FALSE)) {
        fit <- NNS.ANOVA(x, x, means.only = central.only,
                         medians = use.medians, confidence.interval = NULL,
                         plot = FALSE)
        expect_equal(fit$Certainty, maximum)
        # Agreement does not change the public CDF step convention.
        expect_equal(fit$Control_CDF, 1)
        expect_equal(fit$Treatment_CDF, 1)
        expect_equal(NNS.ANOVA(x, x + 1, means.only = central.only,
                               medians = use.medians,
                               confidence.interval = NULL,
                               plot = FALSE)$Certainty, 0)
      }
    }

    expect_equal(unname(NNS.ANOVA(list(x, x, x), plot = FALSE)), maximum)
    pairwise <- NNS.ANOVA(cbind(x, x, x + 1), pairwise = TRUE, plot = FALSE)
    expect_equal(pairwise[1, 2], maximum)
    expect_equal(pairwise[1, 3], 0)

    # The weighted formula must not move a common decimal constant by one ULP.
    unequal <- NNS.ANOVA(rep(0.1, 2), rep(0.1, 4),
                         confidence.interval = NULL, plot = FALSE)
    expect_identical(unequal$Grand_Statistic, 0.1)
    expect_equal(unequal$Certainty, (4 / 6)^2)
  })
  aud_case("a constant and a variable group can agree in mean but not profile", {
    x <- rep(0, 10)
    y <- rep(c(-1, 1), 5)
    maximum <- (18 / 20)^2
    means <- NNS.ANOVA(x, y, means.only = TRUE,
                       confidence.interval = NULL, plot = FALSE)
    profile <- NNS.ANOVA(x, y, confidence.interval = NULL, plot = FALSE)
    expect_equal(means$Certainty, maximum)
    expect_equal(means$Control_CDF, 1)
    expect_equal(means$Treatment_CDF, 0.5)
    expect_gt(profile$Certainty, 0)
    expect_lt(profile$Certainty, means$Certainty)
  })
})

test_that("AOV04: two-sample NNS.ANOVA omits NA and NaN with means, medians and robust = TRUE", {
  withr::local_seed(20261108)
  anova_robust <- fx_anova_robust()
  aud_case("two-sample NNS.ANOVA omits NA and NaN", env = anova_robust, {
    expect_equal(anova(x_na, y_na), anova(x_ok, y_ok))
  })
  aud_case("two-sample NNS.ANOVA omits NA and NaN with medians = TRUE", env = anova_robust, {
    expect_equal(anova(x_na, y_na, medians = TRUE), anova(x_ok, y_ok, medians = TRUE))
  })
  aud_case("two-sample NNS.ANOVA omits NA and NaN with robust = TRUE", env = anova_robust, {
    expect_equal(anova(x_na, y_na, robust = TRUE), anova(x_ok, y_ok, robust = TRUE))
  })
})

test_that("AOV05: robust mode bootstraps each group from its full sample, past min(n1, n2)", {
  withr::local_seed(20261108)
  anova_robust <- fx_anova_robust()
  aud_case("robust mode bootstraps each group independently from its full sample", env = anova_robust, {
    set.seed(314)
    cx <- replicate(100, x[sample.int(length(x), replace = TRUE)])
    ty <- replicate(100, y[sample.int(length(y), replace = TRUE)])
    cs <- vapply(seq_len(100), function(i) {
      NNS:::NNS.ANOVA.bin(cx[, i], ty[, i], plot = FALSE)$Certainty
    }, numeric(1))
    expected <- c(NNS.gravity(cs), unname(quantile(cs, c(0.025, 0.975))))

    r <- anova(x, y, robust = TRUE, seed = 314)
    expect_type(r, "list")
    expect_equal(
      unname(unlist(r[c("Robust Certainty Estimate", "Lower 95% CI", "Upper 95% CI")])),
      expected
    )
  })
  aud_case("robust mode resamples past the first min(n1, n2) values of the longer group", env = anova_robust, {
    set.seed(1)
    long <- c(rnorm(10), rnorm(90, 5))
    short <- rnorm(10)
    r <- anova(long, short, robust = TRUE, seed = 2)
    expect_equal(r$Certainty, 0)
    expect_equal(r[["Robust Certainty Estimate"]], 0)
    expect_equal(r[["Upper 95% CI"]], 0)
  })
})

test_that("AOV06: one-sided tails give one finite effect bound, with and without robust", {
  withr::local_seed(20261108)
  anova_robust <- fx_anova_robust()
  aud_case("tails = 'left' gives a finite lower effect bound and no upper bound", env = anova_robust, {
    b <- boot_stats(32)
    expect_equal(
      bounds(anova(x, y, tails = "left", seed = 32)),
      c(unname(quantile(b$t, 0.05) - quantile(b$c, 0.95)), Inf)
    )
  })
  aud_case("tails = 'right' gives a finite upper effect bound and no lower bound", env = anova_robust, {
    b <- boot_stats(32)
    expect_equal(
      bounds(anova(x, y, tails = "right", seed = 32)),
      c(-Inf, unname(quantile(b$t, 0.95) - quantile(b$c, 0.05)))
    )
  })
  aud_case("one-sided tails also work with robust = TRUE", env = anova_robust, {
    left <- bounds(anova(x, y, tails = "left", robust = TRUE))
    right <- bounds(anova(x, y, tails = "right", robust = TRUE))
    expect_true(is.numeric(left) && is.finite(left[1]) && left[2] == Inf)
    expect_true(is.numeric(right) && right[1] == -Inf && is.finite(right[2]))
  })
})
