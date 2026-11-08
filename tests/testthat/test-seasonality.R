# NNS.seas and the lag-vector generators.
#
# Reference for the scan (#97): the linear trend is removed (keeping the mean),
# each component's dispersion is the Gini mean difference, the coefficient of
# variation is used when the series stays on one side of zero, and lags are
# scanned only up to n / 5 so every component has at least 5 points (#69).
#
# NNS.seas reports the fundamental period (#129): a multiple of the true period
# scores as well in expectation, and the lowest score often fell on a multiple.

test_that("SEAS01: generate.vectors and generate.lin.vectors bounds: empty input, len = 0, lag > n and l = 0", {
  withr::local_seed(20261108)
  aud_case("vector generators do not read out of bounds or divide by zero", {
    g <- NNS:::generate.vectors(c(1, 2, 3), 7L)
    expect_length(g$Component.series[[1]], 0)
    expect_equal(NNS:::generate.vectors(as.numeric(1:10), 3L)$Component.series[[1]], c(2, 5, 8))
    expect_equal(NNS:::generate.vectors(as.numeric(1:10), 3L, 7L)$Component.series[[1]], c(2, 5))
    expect_error(NNS:::generate.vectors(c(1, 2, 3), 2L, 5L), "len")
    lv <- NNS:::generate.lin.vectors(c(1, 2, 3), 7L, 2L)
    expect_length(lv$Component.series[[1]], 0)
  })
  aud_case("generate.vectors gives empty components for empty input, len = 0 and lag > n", {
    for (g in list(NNS:::generate.vectors(numeric(0), 2L),
                   NNS:::generate.vectors(c(1, 2, 3, 4), 2L, 0L),
                   NNS:::generate.vectors(c(1, 2, 3), 7L))) {
      expect_identical(g$Component.series[[1]], numeric(0))
      expect_identical(g$Component.index[[1]], integer(0))
    }
  })
  aud_case("generate.lin.vectors(x, 0) errors instead of killing the R session", {
    # Run in a child process: on the unfixed kernel this call raises SIGFPE.
    script <- tempfile(fileext = ".R")
    on.exit(unlink(script))
    writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
                 "tryCatch(NNS:::generate.lin.vectors(1:4, 0L),",
                 "         error = function(e) cat(conditionMessage(e)))"), script)
    out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                                    stdout = TRUE, stderr = TRUE))
    expect_null(attr(out, "status"))
    expect_match(paste(out, collapse = "\n"), "'l' must be a positive integer")
  })
})

# pin: the ch25 toy example of the book
test_that("SEAS02: NNS.seas scan: AirPassengers period 12, the ch25 toy example, >= 5 points per component, length-weighted components, max.period", {
  withr::local_seed(20261108)
  seas_all_phases <- fx_seas_all_phases()
  aud_case("NNS.seas ranks period 12 first on AirPassengers", env = seas_all_phases, {
    r <- NNS.seas(as.numeric(AirPassengers), plot = FALSE)
    expect_equal(r$best.period, 12L)
    x <- detrend(as.numeric(AirPassengers))
    expect_equal(r$all.periods$Coefficient.of.Variation[1], score(x, 12), tolerance = 1e-10)
    expect_equal(r$all.periods$Variable.Coefficient.of.Variation[1], score(x, 1), tolerance = 1e-10)
  })
  aud_case("NNS.seas reproduces the ch25 lag-2 toy example once components have 5 points", env = seas_all_phases, {
    # The 8-point book series has 4 points per lag-2 component: below the minimum.
    expect_equal(NNS.seas(c(10, 18, 11, 21, 12, 20, 13, 23), plot = FALSE)$periods, 1L)
    x <- c(10, 18, 11, 21, 12, 20, 13, 23, 14, 22)
    r <- NNS.seas(x, plot = FALSE)
    expect_equal(r$periods, 2L)
    expect_equal(r$all.periods$Coefficient.of.Variation, score(detrend(x), 2), tolerance = 1e-10)
    # Without detrending the components are the book's (10, ..., 14) and (18, ..., 22).
    r0 <- NNS.seas(x, plot = FALSE, detrend = FALSE)
    expect_equal(r0$all.periods$Coefficient.of.Variation, score(x, 2), tolerance = 1e-12)
  })
  aud_case("NNS.seas scans only lags with at least 5 points per component", {
    set.seed(69)
    for (n in c(60, 200, 1000)) {
      r <- NNS.seas(rnorm(n, 10), plot = FALSE)
      expect_true(all(r$all.periods$Period <= n %/% 5))
    }
    # Pure noise no longer flags nearly every lag (493 of 500 before).
    set.seed(1)
    flagged <- replicate(5, sum(NNS.seas(rnorm(1000, 10), plot = FALSE)$all.periods$Period > 1))
    expect_lt(mean(flagged), 120)
  })
  aud_case("each lag m scores all m component series, weighted by length", env = seas_all_phases, {
    x <- detrend(as.numeric(lynx))
    n <- length(x)
    ref <- sapply(1:(n %/% 5), function(m) score(x, m))
    pass <- which(ref < ref[1])
    r <- NNS.seas(as.numeric(lynx), plot = FALSE)$all.periods
    expect_equal(sort(r$Period), pass)
    expect_equal(r$Coefficient.of.Variation[order(r$Period)], ref[pass], tolerance = 1e-10)
  })
  aud_case("max.period caps the scan", {
    set.seed(4); x <- 10 + rnorm(500)
    expect_true(all(NNS.seas(x, plot = FALSE, max.period = 7)$all.periods$Period <= 7))
    expect_error(NNS.seas(x, plot = FALSE, max.period = 0), "max.period")
  })
})

test_that("SEAS03: dispersion rules: zero-crossing series, strict CV, the no-seasonality weight, constant phases and series, interval-scale offsets", {
  withr::local_seed(20261108)
  seas_all_phases <- fx_seas_all_phases()
  aud_case("a series that crosses zero compares dispersions, not CVs", env = seas_all_phases, {
    x <- sin(2 * pi * (1:96) / 12)
    x <- x - mean(x)
    r <- NNS.seas(x, plot = FALSE, detrend = FALSE)$all.periods
    expect_equal(r$Variable.Coefficient.of.Variation[1], gmd(x))
    # 24 and 36 have fewer than 5 points per component at n = 96.
    expect_equal(r$Period[1], 12L)
    expect_lt(r$Coefficient.of.Variation[1], 1e-12)
  })
  aud_case("no lag is reported unless its CV is strictly below CV(X)", {
    expect_equal(NNS.seas(rep(5, 20), plot = FALSE)$periods, 1L)
  })
  aud_case("the no-seasonality row gets weight 1, not 0", {
    M <- data.frame(Period = 1L, Coefficient.of.Variation = 0.5, Variable.Coefficient.of.Variation = 0.5)
    expect_equal(NNS:::ARMA.seas.weighting(FALSE, M)$Weights, 1)
  })
  aud_case("a constant zero phase scores 0 and keeps its lag", env = seas_all_phases, {
    x <- rep(c(0, 10), 10)
    r <- NNS.seas(x, plot = FALSE, detrend = FALSE)$all.periods
    expect_equal(sort(r$Period), c(2L, 4L))
    expect_equal(r$Coefficient.of.Variation, rep(0, 2))
    expect_equal(NNS.ARMA(x, h = 2, plot = FALSE, seasonal.plot = FALSE), c(0, 10))
  })
  aud_case("a constant decimal series has no seasonal period", {
    expect_equal(NNS.seas(rep(.1, 20), plot = FALSE)$best.period, 1)
    expect_equal(NNS.seas(rep(c(.1, 1), 10), plot = FALSE)$best.period, 2)
  })
  aud_case("absolute seasonality scoring is invariant to interval-scale offsets", {
    t <- seq_len(120)
    x <- 15 + sin(2 * pi * t / 12) + .2 * cos(t)
    a <- NNS.seas(x, plot = FALSE, dispersion = "absolute")
    b <- NNS.seas(x + 273.15, plot = FALSE, dispersion = "absolute")
    expect_equal(a$periods, b$periods)
    expect_equal(a$all.periods, b$all.periods, tolerance = 1e-10)
  })
})

test_that("SEAS04: trend and noise: noisy periods found, linear trends removed, exactly linear series constant, noise flags level-free", {
  withr::local_seed(20261108)
  aud_case("NNS.seas finds the true period of a noisy seasonal series", {
    # The best period is the true period or one of its multiples (which share
    # its seasonal structure); before #69 it was a long non-multiple lag.
    for (p in c(7, 12, 24)) for (s in 1:5) {
      set.seed(100 * p + s); x <- 10 + 2 * sin(2 * pi * (1:240) / p) + rnorm(240)
      r <- NNS.seas(x, plot = FALSE)
      expect_equal(r$best.period %% p, 0)
      expect_true(p %in% r$periods)
    }
  })
  aud_case("a linear trend does not hide the seasonal period (#97)", {
    bp <- sapply(1201:1210, function(s) {
      set.seed(s); t <- 1:240
      NNS.seas(10 + .07 * t + 2 * sin(2 * pi * t / 12) + rnorm(240), plot = FALSE)$best.period
    })
    expect_true(all(bp %% 12 == 0))                 # 0 of 10 before
    expect_equal(NNS.seas(1:50 + 0.5, plot = FALSE)$best.period, 1)
  })
  aud_case("an exactly linear series detrends to a constant, not a period", {
    for (a in c(0, 10, -50, 1e6, 1e9)) for (b in c(.3, -2, 1e-3)) for (n in c(12L, 50L, 144L)) {
      r <- NNS.seas(a + b * seq_len(n), plot = FALSE)
      expect_equal(r$best.period, 1, info = paste(a, b, n))
      expect_true(all(r$all.periods$Coefficient.of.Variation >= 0), info = paste(a, b, n))
    }
    r <- NNS.seas(10 + .3 * (1:50), plot = FALSE)$all.periods
    expect_equal(c(r$Coefficient.of.Variation, r$Variable.Coefficient.of.Variation), c(0, 0))
    t <- 1:120
    expect_equal(NNS.seas(1e6 + .3 * t + sin(2 * pi * t / 12), plot = FALSE)$best.period, 12)
    # Small or slow periods on a high, long series are not rounding noise.
    t <- seq_len(12000)
    expect_equal(NNS.seas(1e9 + .3 * t + .01 * sin(2 * pi * t / 12), plot = FALSE, max.period = 12)$best.period, 12)
    x <- 1e9 + .3 * t + .1 * cos(2 * pi * (t - .5) / 1200)
    expect_equal(NNS.seas(x, plot = FALSE, max.period = 1200)$best.period, 1200)
    expect_equal(NNS.seas(1e9 + cumsum(rep(.3, 1e5)), plot = FALSE, max.period = 24)$best.period, 1)
  })
  aud_case("noise flags do not depend on the level or the tail shape (#97)", {
    flagged <- function(gen) mean(sapply(1:5, function(s) {
      set.seed(s); sum(NNS.seas(gen(), plot = FALSE)$periods > 1)
    }))
    # Before: 170 of 200 at mean 0, 0 at mean 1, 198 for t3 noise at mean 10.
    for (v in c(flagged(function() rnorm(1000)), flagged(function() 1 + rnorm(1000)),
                flagged(function() 10 + rt(1000, 3))))
      expect_true(v > 60 && v < 140)
  })
})

test_that("SEAS05: p.value filter: shuffled copies, modulo periods beyond the scan and the fundamental period", {
  withr::local_seed(20261108)
  aud_case("p.value keeps periods that beat the best lag of shuffled copies (#97)", {
    fp <- sapply(1:10, function(s) {
      set.seed(s); NNS.seas(10 + rnorm(300), plot = FALSE, p.value = TRUE, nperm = 40)$best.period
    })
    expect_lte(sum(fp > 1), 2)
    set.seed(3); t <- 1:240
    r <- NNS.seas(10 + 2 * sin(2 * pi * t / 12) + rnorm(240), plot = FALSE, p.value = TRUE, nperm = 40)
    expect_equal(r$best.period %% 12, 0)
    expect_true(all(r$all.periods$p.value[match(r$periods, r$all.periods$Period)] <= 0.05))
  })
  aud_case("p.value tests modulo periods beyond the scan against a null that scores them", {
    # Each shuffled copy scores the same modulo-added periods as the observed
    # series; with a null that only scanned up to n / 5, noise passed in 27 of
    # 40 series here.
    hits <- sapply(1:20, function(s) {
      set.seed(s)
      a <- NNS.seas(10 + rnorm(96), modulo = 20:47, p.value = TRUE, nperm = 40, plot = FALSE)
      any(a$periods > 1)
    })
    expect_lte(sum(hits), 3)
    set.seed(3)
    x <- 10 + sin(2 * pi * (1:96) / 24) + rnorm(96, sd = .3)
    a <- NNS.seas(x, modulo = 24, mod.only = FALSE, p.value = TRUE, nperm = 40, plot = FALSE)
    expect_equal(a$best.period, 24L)
  })
  aud_case("the permutation filter applies to the fundamental period", {
    set.seed(1); x <- 10 + sin(2 * pi * (1:540) / 12) + rnorm(540, sd = .3)
    set.seed(3); out <- NNS.seas(x, plot = FALSE, p.value = TRUE, nperm = 30)
    expect_equal(out$best.period, 12)
    set.seed(3); w <- rnorm(300, 50, 5)
    set.seed(1); expect_equal(NNS.seas(w, plot = FALSE, p.value = TRUE, nperm = 30)$best.period, 1L)
  })
})

test_that("SEAS06: fundamental period: delete-one scores, 12 not a multiple, a longer period not a divisor, the full table", {
  withr::local_seed(20261108)
  aud_case("delete-one scores match a direct recomputation", {
    set.seed(5)
    x <- 10 + rexp(57)
    direct <- function(x, m, use_cv) {
      grp <- (seq_along(x) - 1L) %% m
      score <- function(keep) {
        comp <- split(x[keep], grp[keep])
        d <- vapply(comp, function(z) {
          k <- length(z); g <- sum(abs(outer(z, z, "-"))) / (k * (k - 1))
          if (use_cv) g / abs(mean(z)) else g
        }, numeric(1))
        sum(lengths(comp) * d) / sum(keep)
      }
      vapply(seq_along(x), function(i) score(seq_along(x) != i), numeric(1))
    }
    for (m in c(3L, 7L)) for (cv in c(TRUE, FALSE))
      expect_equal(NNS:::.nns_seas_delete_one(x, m, cv), direct(x, m, cv))
    expect_null(NNS:::.nns_seas_delete_one(x, 20L, TRUE))   # components of 2 or 3: too short
  })
  aud_case("a monthly cycle is reported at 12, not at a multiple", {
    t <- 1:540
    p <- vapply(1:6, function(s) {
      set.seed(s); x <- 10 + sin(2 * pi * t / 12) + rnorm(540, sd = .3)
      out <- NNS.seas(x, plot = FALSE)
      c(best = out$best.period, first = out$all.periods$Period[1L], periods = out$periods[1L],
        lowest = with(out$all.periods, Period[which.min(Coefficient.of.Variation)]))
    }, numeric(4))
    expect_true(all(p["best", ] == 12))
    expect_true(all(p["first", ] == 12) && all(p["periods", ] == 12))
    expect_true(all(p["lowest", ] %% 12 == 0 & p["lowest", ] > 12))   # the lowest score sat on a multiple
    expect_equal(NNS.seas(nottem, plot = FALSE)$best.period, 12)
  })
  aud_case("a longer true period is not reduced to a divisor", {
    t <- 1:540
    for (s in 1:3) {
      set.seed(s)
      expect_equal(NNS.seas(10 + sin(2 * pi * t / 24) + rnorm(540, sd = .3), plot = FALSE)$best.period, 24)
      set.seed(s)
      mix <- 10 + sin(2 * pi * t / 12) + .5 * sin(2 * pi * t / 24) + rnorm(540, sd = .3)
      expect_equal(NNS.seas(mix, plot = FALSE)$best.period, 24)
    }
  })
  aud_case("the table keeps every scored period and its score", {
    set.seed(2); x <- 10 + sin(2 * pi * (1:540) / 12) + rnorm(540, sd = .3)
    out <- NNS.seas(x, plot = FALSE)
    raw <- NNS:::NNS_seas_cpp(variable = x - coef(lm(x ~ seq_along(x)))[2] * (seq_along(x) - mean(seq_along(x))),
                              modulo = NULL, mod_only = TRUE, max_period = NULL)$all.periods
    expect_setequal(out$all.periods$Period, raw$Period)
    expect_equal(out$all.periods[order(out$all.periods$Period), "Coefficient.of.Variation"],
                 raw[order(raw$Period), "Coefficient.of.Variation"])
    expect_equal(out$all.periods$Period[-1L], raw$Period[raw$Period != out$best.period])
  })
})
