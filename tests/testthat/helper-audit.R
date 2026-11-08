# Shared helpers and fixtures of the audit tests (test-partial-moments-*.R,
# test-moments.R, test-dominance-verdicts.R, test-input-contract.R). Every
# name starts with aud_ so nothing here shadows an upstream test object.

# ---- case runner -----------------------------------------------------------

# Runs one labelled case in its own environment: a failed expectation is
# reported with "[label] " in front, and an error ends this case only, as one
# labelled failure. Signalling after the restart has unwound keeps an error
# from being reported twice (tryCatch would swallow failed expectations).
aud_case <- function(label, code) {
  expr <- substitute(code)
  env <- new.env(parent = parent.frame())
  tag <- function(e) paste0("[", label, "] ", conditionMessage(e))
  withRestarts(
    withCallingHandlers(eval(expr, env),
      expectation_failure = function(e) {
        e$message <- tag(e)
        testthat::exp_signal(e)
        invokeRestart("continue_test")
      },
      error = function(e) {
        if (!inherits(e, "expectation")) invokeRestart("aud_case_error", e)
      }),
    aud_case_error = function(e) testthat::fail(paste0(tag(e), " (unexpected error)")))
  invisible()
}

# Runs f(row) for every row of the data frame `cases`, labelled by cases$label.
aud_cases <- function(cases, f) {
  for (i in seq_len(nrow(cases))) {
    row <- cases[i, , drop = FALSE]
    aud_case(row$label, f(row))
  }
}

# A table of error cases sharing one message: each further argument is a call
# that must raise `msg` (fixed match), run as its own case labelled by the call.
aud_errors <- function(msg, ...) {
  calls <- as.list(substitute(list(...)))[-1]
  env <- parent.frame()
  for (cl in calls) {
    aud_case(paste(deparse(cl, width.cutoff = 500L), collapse = " "),
             expect_error(eval(cl, env), msg, fixed = TRUE))
  }
}

aud_memo <- function(f) {
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- f()
    cache
  }
}

aud_eps <- .Machine$double.eps

# ---- tolerances and comparisons --------------------------------------------

# Largest excess of |lhs - rhs| over the tolerance; must be <= 0.
aud_within_max <- function(lhs, rhs, tol, label = NULL) {
  expect_lte(max(abs(lhs - rhs) - tol), 0, label = label)
}

aud_close_flat <- function(actual, expected, tol) {
  expect_lte(max(abs(unname(unlist(actual)) - unname(unlist(expected)))), tol)
}

# Elementwise |actual - expected| <= tol, one expectation, failing entries
# named. An empty comparison fails (a case list that collapsed to nothing
# would otherwise pass vacuously); tol is a scalar or matches actual.
aud_within_cells <- function(actual, expected, tol) {
  non_empty <- length(actual) > 0
  same_shape <- length(expected) %in% c(1L, length(actual))
  tol_shape <- length(tol) %in% c(1L, length(actual))
  if (!(non_empty && same_shape && tol_shape)) {
    return(expect(FALSE, paste(c(if (!non_empty) "empty comparison",
                                 if (!same_shape) "length mismatch (expected)",
                                 if (!tol_shape) "length mismatch (tol)"), collapse = "; ")))
  }
  err <- abs(actual - expected)
  tol <- rep_len(tol, length(err))
  bad <- which(!(is.finite(err) & err <= tol))
  where <- names(actual)
  if (is.null(where)) where <- as.character(seq_along(err))
  msg <- sprintf("%s: |diff| = %.3g > tol %.3g", where[bad], err[bad], tol[bad])
  expect(length(bad) == 0, paste(c("not within tolerance:", head(msg, 8)), collapse = "\n"))
}

aud_close <- function(object, expected, tol, info = NULL) {
  expect_equal(length(object), length(expected), info = info)
  # NaN-safe: a NaN difference is a failure, not an error.
  dev <- max(abs(object - expected))
  expect_true(is.finite(dev) && dev <= tol,
              info = paste("max |diff| =", format(dev, digits = 3), "> tol =", format(tol, digits = 3), info))
}

# ---- univariate partial moments ---------------------------------------------

# Absolute deviation moment sum(|x - t|^d) / n: the natural scale of every
# partial-moment identity at (d, t). Naive summation of n terms, each with
# O(d) rounding, has absolute error below (n + d) * eps * scale; the factor 2
# covers the R-side oracle as well.
aud_pm_tol <- function(d, t, x) {
  scale <- sum(abs(x - t)^d) / length(x)
  2 * (length(x) + d + 2) * .Machine$double.eps * scale
}
aud_pm_tol_v <- function(d, ts, x) vapply(ts, function(t) aud_pm_tol(d, t, x), numeric(1))

# One call per target keeps LPM/UPM on the direct (<= 32 targets) kernel.
aud_lpm_each <- function(d, ts, x) vapply(ts, function(t) LPM(d, t, x), numeric(1))
aud_upm_each <- function(d, ts, x) vapply(ts, function(t) UPM(d, t, x), numeric(1))

aud_targets <- function(x) {
  c(mean(x), median(x), quantile(x, c(0.05, 0.25, 0.75, 0.95), names = FALSE),
    x[1], min(x), max(x))
}

# 32 extra targets force the prefix backend.
aud_pad <- seq(-1, 1, length.out = 32)

# Samples: exact (n = 8, dyadic, a tie at 0.5), random with ties forced by
# rounding, real. `kept` is the covering design of the sample axis: one random
# sample with ties, one large-n small-scale real sample, one small-n real
# sample with ties. The cut samples are still drawn so the seeded data match.
aud_fx_uni <- aud_memo(function() withr::with_seed(20261108, {
  random_samples <- list(
    normal = round(rnorm(400), 1),
    expo   = round(rexp(300, 2), 2),
    t3     = round(rt(500, df = 3), 1)
  )
  real_samples <- list(
    dax      = as.numeric(diff(log(EuStockMarkets[, "DAX"]))),
    waiting  = faithful$waiting,
    ozone    = as.numeric(na.omit(airquality$Ozone)),
    mag      = quakes$mag,
    mpg      = mtcars$mpg
  )
  all_samples <- c(random_samples, real_samples)
  list(x_exact = c(-2, -1, -0.5, 0, 0.5, 0.5, 1.5, 3),
       t_exact = c(-3, -1, 0, 0.25, 0.5, 1, 3, 4),
       random_samples = random_samples, real_samples = real_samples,
       all_samples = all_samples, kept = all_samples[c("normal", "dax", "mpg")])
}))

# ---- ratios, CDF and VaR ----------------------------------------------------

# Absolute tolerance for a ratio in [0, 1] computed from n terms.
aud_ratio_tol <- function(x) 8 * length(x) * aud_eps

aud_targets_for <- function(x) {
  unique(c(mean(x), median(x), quantile(x, c(.05, .25, .75, .95), names = FALSE),
           x[1], x[length(x) %/% 2], min(x), max(x)))
}

aud_ecdf_count <- function(t, x) vapply(t, function(s) sum(x <= s), numeric(1)) / length(x)
aud_surv_count <- function(t, x) vapply(t, function(s) sum(x > s), numeric(1)) / length(x)

# More than 32 targets: the prefix backend. The margin is relative to the
# range so that most grid points are strictly inside (min x, max x) for every
# sample (dax has range ~0.15: an absolute margin of 1 left it 4 interior points).
aud_grid_for <- function(x, k = 64) {
  m <- 0.1 * diff(range(x))
  sort(c(seq(min(x) - m, max(x) + m, length.out = k), x[1:min(8, length(x))]))
}

aud_cdf_fun <- function(x, d = 0, type = "CDF", target = NULL) {
  NNS.CDF(x, d, target = target, type = type, plot = FALSE)
}

aud_var_p <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

# exact: n = 8, dyadic, ties at 0 and 1, dyadic mean 0.34375.
aud_fx_ratio <- aud_memo(function() withr::with_seed(20261108, {
  x_exact <- c(-2, -0.5, 0, 0, 0.25, 1, 1, 3)
  samples <- list(
    exact   = x_exact,
    norm    = round(rnorm(200), 2),
    exp     = round(rexp(150), 1),
    t3      = round(rt(300, df = 3), 1),
    dax     = as.numeric(diff(log(EuStockMarkets))[, "DAX"]),
    waiting = faithful$waiting,
    ozone   = as.numeric(na.omit(airquality$Ozone)),
    mag     = quakes$mag,
    mpg     = mtcars$mpg
  )
  list(x_exact = x_exact, samples = samples,
       kept = samples[c("exact", "norm", "dax", "mpg")])
}))

# ---- moments, integration and dominance ------------------------------------

# Absolute scale of a degree-k moment about t: mean(|x - t|^k).
aud_abs_mom <- function(x, k, t = mean(x)) mean(abs(x - t)^k)

# Naive summation of n terms of degree k: error below (n + k) * eps * scale.
aud_mom_tol <- function(x, k, t = mean(x)) 4 * (length(x) + k + 2) * aud_eps * aud_abs_mom(x, k, t)

# Partial-moment assembly of NNS.moments (momentsPM in Moments.lean).
aud_moments_pm <- function(x, population = TRUE) {
  n <- length(x)
  mu <- mean(x)
  V <- LPM(2, mu, x) + UPM(2, mu, x)
  S <- UPM(3, mu, x) - LPM(3, mu, x)
  K <- UPM(4, mu, x) + LPM(4, mu, x)
  if (population) {
    return(list(mean = mu, variance = V, skewness = S / V^1.5, kurtosis = K / V^2 - 3))
  }
  s2 <- V * n / (n - 1)
  list(mean = mu, variance = s2,
       skewness = n / ((n - 1) * (n - 2)) * (n * S) / s2^1.5,
       kurtosis = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * (n * K) / s2^2 -
         3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}

# Absolute tolerances for the standardized moments (both sides carry the
# summation error of the degree-k sum, relative to V^(k/2)).
aud_moments_tol <- function(x) {
  V <- aud_abs_mom(x, 2)
  c(mean = 4 * aud_eps * max(abs(x)),
    variance = aud_mom_tol(x, 2),
    skewness = 4 * aud_mom_tol(x, 3) / V^1.5 + 8 * length(x) * aud_eps * aud_abs_mom(x, 3) / V^1.5,
    kurtosis = 4 * aud_mom_tol(x, 4) / V^2 + 8 * length(x) * aud_eps * aud_abs_mom(x, 4) / V^2)
}

aud_moments_equal <- function(got, want, tol, label) {
  for (nm in c("mean", "variance", "skewness", "kurtosis")) {
    aud_within_max(got[[nm]], want[[nm]], tol[[nm]], label = paste(label, nm))
  }
}

# Exact quadrature of a piecewise polynomial on consecutive knots (the LPM of
# degree d is a polynomial of degree d between observations); cumulative
# integrals from k[1] to each knot.
aud_cum_step <- function(k, f) c(0, cumsum(head(f(k), -1) * diff(k)))            # d = 0, right-continuous
aud_cum_trap <- function(k, f) { v <- f(k); c(0, cumsum((head(v, -1) + v[-1]) / 2 * diff(k))) }
aud_cum_simpson_half <- function(k, f) {  # 3 * integral, kept dyadic: h * (fa + 4 fm + fb) / 2
  v <- f(k); m <- f((head(k, -1) + k[-1]) / 2)
  c(0, cumsum(diff(k) * (head(v, -1) + 4 * m + v[-1]) / 2))
}

# One point below every observation, the sorted distinct observations, one
# point above. refine = TRUE adds the midpoints (still a partition on whose
# cells LPM_d is one polynomial), which also pushes small real samples past 32
# targets onto the prefix backend.
aud_knot_grid <- function(x, refine = FALSE) {
  u <- sort(unique(x))
  k <- c(u[1] - 1, u, u[length(u)] + 1)
  if (refine) k <- sort(c(k, (head(k, -1) + k[-1]) / 2))
  k
}

# The R verdict for SSD on the merged knots (SPEC "sd_dom_pair, SSD"), as a
# function of the knot differences d = LPM1_x - LPM1_y. NA = inside the band.
aud_ssd_expected <- function(d, tol) {
  if (max(d) > tol) return(0L)
  if (all(d <= 0) && any(d < -tol)) return(1L)
  NA_integer_
}

# Implementation-free LPM(1, t, x) = mean(max(t - x, 0)), the SSD oracle.
aud_lpm1_ref <- function(ts, x) vapply(ts, function(t) mean(pmax(t - x, 0)), numeric(1))

# Exact TSD criterion of sdDomPair_tsd_iff for integer samples. Sufficiency is
# LPM_two_le_of_le_at_knots (knots, interior vertices s* = b / a, mean gate);
# strictness is checked at the knots only, which LPM_two_eq_of_eq_at_knots
# makes exact. All quantities are integers scaled by mx * my.
aud_tsd_parts <- function(x, y) {
  mx <- length(x); my <- length(y)
  k <- sort(unique(c(x, y)))
  A <- vapply(k, function(t) sum(x <= t) * my - sum(y <= t) * mx, numeric(1))
  B <- vapply(k, function(t) sum(x[x <= t]) * my - sum(y[y <= t]) * mx, numeric(1))
  C <- vapply(k, function(t) sum(x[x <= t]^2) * my - sum(y[y <= t]^2) * mx, numeric(1))
  q_knot <- A * k^2 - 2 * B * k + C
  lo <- head(k, -1); hi <- k[-1]; a <- head(A, -1); b <- head(B, -1); cc <- head(C, -1)
  inside <- a < 0 & lo * a > b & b > hi * a
  list(knots_ok = all(q_knot <= 0),
       vertex_bad = any(inside & (cc * a - b^2 < 0)),
       mean_ok = sum(x) * my >= sum(y) * mx,
       differs = any(q_knot != 0))
}
aud_tsd_exact <- function(x, y) {
  p <- aud_tsd_parts(x, y)
  as.integer(p$knots_ok && !p$vertex_bad && p$mean_ok && p$differs)
}

# Exact continuous-FSD (Omega) criterion of sdDomPair_omega_iff for integer
# samples: R_x(t) <= R_y(t) for all t, R = LPM.ratio(1, t, .) = L1 / (L1 + U1)
# (0/0 -> k/m), and the curves differ somewhere. Counts and prefix sums only.
# On a knot cell R_x <= R_y iff g(t) = L1x * U1y - L1y * U1x <= 0, a quadratic
# checked at the knots and at an interior vertex; a quadratic that vanishes at
# both knots and at the midpoint is zero, so "differs" is knots + midpoints.
aud_omega_parts <- function(x, y) {
  mx <- length(x); my <- length(y); sx <- sum(x); sy <- sum(y)
  k <- sort(unique(c(x, y)))
  ax <- vapply(k, function(t) sum(x <= t), numeric(1))
  px <- vapply(k, function(t) sum(x[x <= t]), numeric(1))
  ay <- vapply(k, function(t) sum(y <= t), numeric(1))
  py <- vapply(k, function(t) sum(y[y <= t]), numeric(1))
  lx <- ax * k - px; dx <- sx - 2 * px + (2 * ax - mx) * k   # m * L1, m * (L1 + U1)
  ly <- ay * k - py; dy <- sy - 2 * py + (2 * ay - my) * k
  nx <- ifelse(dx == 0, ax, lx); dx <- ifelse(dx == 0, mx, dx)
  ny <- ifelse(dy == 0, ay, ly); dy <- ifelse(dy == 0, my, dy)
  g_knot <- nx * dy - ny * dx                                 # sign(R_x - R_y)
  i <- seq_len(length(k) - 1)
  a1 <- ax[i]; p1 <- px[i]; a2 <- ay[i]; p2 <- py[i]; lo <- k[i]; hi <- k[i + 1]
  al <- -a1 * (my - a2) + a2 * (mx - a1)
  be <- a1 * (sy - p2) + p1 * (my - a2) - a2 * (sx - p1) - p2 * (mx - a1)
  ga <- -p1 * (sy - p2) + p2 * (sx - p1)
  inside <- al < 0 & -be > 2 * al * hi & -be < 2 * al * lo   # -be / (2 al) in (lo, hi)
  mid <- lo + hi
  list(knots_ok = all(g_knot <= 0),
       vertex_bad = any(inside & (4 * al * ga - be^2 < 0)),
       knot_differs = any(g_knot != 0),
       mid_differs = any(al * mid^2 + 2 * be * mid + 4 * ga != 0))
}
aud_omega_exact <- function(x, y) {
  p <- aud_omega_parts(x, y)
  as.integer(p$knots_ok && !p$vertex_bad && (p$knot_differs || p$mid_differs))
}

aud_verdicts <- function(p) c(
  fsd = NNS.FSD.uni(p$x, p$y, type = "discrete"),
  ssd = NNS.SSD.uni(p$x, p$y),
  tsd = NNS.TSD.uni(p$x, p$y)
)

# Samples and dominance pairs. The four seeds are the ones the data have
# always been drawn with, so the pairs (and their verdict tables) are fixed.
aud_fx_md <- aud_memo(function() {
  fx <- list(
    # exact regime: n = 8, integers with ties at 1 and 4; mean 3
    x_int = c(0, 1, 1, 3, 4, 4, 4, 7),
    # dyadic sample with mean exactly 0 and a tie at 0.25
    x_dy = c(-2, -1, -0.5, 0, 0.5, 0.5, 1.5, 3) - 0.25,
    t_dy = c(-3, -1.25, -0.25, 0, 0.25, 1, 2.75, 4)
  )
  fx$random_samples <- withr::with_seed(20261108, list(
    normal = round(rnorm(400), 1),
    expo   = round(rexp(300, 2), 2),
    t3     = round(rt(500, df = 3), 1)
  ))
  fx$real_samples <- list(
    dax       = as.numeric(diff(log(EuStockMarkets[, "DAX"]))),
    eruptions = faithful$eruptions,
    ozone     = as.numeric(na.omit(airquality$Ozone)),
    mag       = quakes$mag,
    mpg       = mtcars$mpg
  )
  fx$all_samples <- c(fx$random_samples, fx$real_samples)
  fx$kept <- fx$all_samples[c("normal", "dax", "mpg")]
  # Random small integer pairs (values 0..6, sizes 1..6).
  fx$int_pairs <- withr::with_seed(7, replicate(400, list(
    x = sample(0:6, sample(1:6, 1), replace = TRUE),
    y = sample(0:6, sample(1:6, 1), replace = TRUE)), simplify = FALSE))
  # Pairs that are FSD by construction: y_i <= x_i, some strictly, ties kept.
  fx$fsd_pairs <- withr::with_seed(11, replicate(100, {
    x <- round(rt(sample(5:40, 1), df = 3), 1)
    list(x = x, y = round(x - rexp(length(x)) * rbinom(length(x), 1, 0.5), 1))
  }, simplify = FALSE))
  fx$int_verdicts <- t(vapply(fx$int_pairs, aud_verdicts, integer(3)))
  fx$int_verdicts_rev <- t(vapply(fx$int_pairs, function(p) aud_verdicts(list(x = p$y, y = p$x)), integer(3)))
  fx$fsd_verdicts <- t(vapply(fx$fsd_pairs, aud_verdicts, integer(3)))
  fx$int_omega <- vapply(fx$int_pairs, function(p) NNS.FSD.uni(p$x, p$y, type = "continuous"), integer(1))
  fx$int_omega_rev <- vapply(fx$int_pairs, function(p) NNS.FSD.uni(p$y, p$x, type = "continuous"), integer(1))
  # Pairs on which the interior TSD vertex check fires. On int_pairs the vertex
  # never decides a verdict, so reject-sample a wider regime (values 0..5, sizes
  # 2..15) until 25 pairs pass the knots, the min gate and the mean gate and are
  # rejected by the vertex alone; keep up to 200 vertex-firing pairs overall.
  # The 1e5 iteration cap is a safety net only: seed 23 reaches 25 in about 2.5e4.
  withr::with_seed(23, {
    vertex_pairs <- list(); n_vertex_decided <- 0L; iter <- 0L
    while (n_vertex_decided < 25L && iter < 100000L) {
      iter <- iter + 1L
      x <- sample(0:5, sample(2:15, 1), replace = TRUE)
      y <- sample(0:5, sample(2:15, 1), replace = TRUE)
      pt <- aud_tsd_parts(x, y)
      if (!pt$vertex_bad) next
      decided <- pt$knots_ok && pt$mean_ok && pt$differs && min(x) >= min(y)
      if (decided || length(vertex_pairs) < 200L) {
        vertex_pairs[[length(vertex_pairs) + 1L]] <- list(x = x, y = y, decided = decided)
        n_vertex_decided <- n_vertex_decided + decided
      }
    }
  })
  fx$vertex_pairs <- vertex_pairs
  fx$n_vertex_decided <- n_vertex_decided
  fx
})

# ---- co- and divergent partial moments --------------------------------------

# Natural scale of a co-moment identity at degrees (dx, dy): the mean of
# |x - tx|^dx |y - ty|^dy. Naive in-order summation of n terms with O(d)
# rounding per term stays below (n + d) * eps * scale; k covers the number of
# moments combined plus the R-side oracle.
aud_co_tol <- function(c, dx, dy, k = 4) {
  k * (length(c$x) + dx + dy + 2) * aud_eps * mean(abs(c$x - c$tx)^dx * abs(c$y - c$ty)^dy)
}

# Matrix scale: entries are compared against the largest entry.
aud_mat_tol <- function(M, n, k = 4) k * (n + ncol(M)) * aud_eps * max(abs(M))

# Collect f over named cases / over degree pairs into one named vector.
aud_over <- function(cases, f) unlist(lapply(cases, f))
aud_at_degrees <- function(pairs, f) {
  unlist(lapply(pairs, function(p) setNames(f(p[1], p[2]), paste0("d", p[1], p[2]))))
}

aud_min_eigen <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

# One case per target pair: mean, median, quartiles and an observation
# (x[i], y[i]) so that both coordinates tie their target.
aud_co_cases <- function(x, y, label) {
  i <- ceiling(length(x) / 3)
  q <- function(v, p) quantile(v, p, names = FALSE)
  targets <- list(
    mean   = c(mean(x), mean(y)),
    median = c(median(x), median(y)),
    q25    = c(q(x, 0.25), q(y, 0.25)),
    q75    = c(q(x, 0.75), q(y, 0.75)),
    point  = c(x[i], y[i])
  )
  lapply(setNames(targets, paste(label, names(targets))), function(t)
    list(x = x, y = y, tx = t[1], ty = t[2]))
}

aud_co_lpm <- function(c, dx, dy) Co.LPM(dx, c$x, c$y, c$tx, c$ty, dy)
aud_co_upm <- function(c, dx, dy) Co.UPM(dx, c$x, c$y, c$tx, c$ty, dy)
aud_d_lpm  <- function(c, dl, du) D.LPM(dl, du, c$x, c$y, c$tx, c$ty)
aud_d_upm  <- function(c, dl, du) D.UPM(dl, du, c$x, c$y, c$tx, c$ty)

aud_co_moment <- function(c, d = 1, e = 1) mean((c$x - c$tx)^d * (c$y - c$ty)^e)
aud_decomposition <- function(c) {
  aud_co_upm(c, 1, 1) + aud_co_lpm(c, 1, 1) - aud_d_lpm(c, 1, 1) - aud_d_upm(c, 1, 1)
}
aud_quadrant_sum0 <- function(c) {
  aud_co_lpm(c, 0, 0) + aud_co_upm(c, 0, 0) + aud_d_lpm(c, 0, 0) + aud_d_upm(c, 0, 0)
}
aud_mixed_master_lhs <- function(c, d, e) {
  aud_co_upm(c, d, e) + (-1)^d * aud_d_upm(c, d, e) + (-1)^e * aud_d_lpm(c, e, d) +
    (-1)^(d + e) * aud_co_lpm(c, d, e)
}

aud_self_cases <- function(x, ts) lapply(ts, function(t) list(x = x, y = x, tx = t, ty = t))
aud_self_tol_40 <- function(c, d, e) {
  2 * (length(c$x) + d + e + 2) * aud_eps *
    vapply(c$tx, function(t) mean(abs(c$x - t)^(d + e)), numeric(1))
}
aud_over_d <- function(cases, f) aud_over(cases, function(c) vapply(0:3, function(d) f(c, d), numeric(1)))
aud_cs_slack <- function(c, d, e) {
  rel <- 1 + 8 * (length(c$x) + 2 * d + 2 * e + 2) * aud_eps
  c(CoLPM = LPM(2 * d, c$tx, c$x) * LPM(2 * e, c$ty, c$y) * rel - aud_co_lpm(c, d, e)^2,
    CoUPM = UPM(2 * d, c$tx, c$x) * UPM(2 * e, c$ty, c$y) * rel - aud_co_upm(c, d, e)^2,
    DLPM  = UPM(2 * e, c$tx, c$x) * LPM(2 * d, c$ty, c$y) * rel - aud_d_lpm(c, d, e)^2,
    DUPM  = LPM(2 * d, c$tx, c$x) * UPM(2 * e, c$ty, c$y) * rel - aud_d_upm(c, d, e)^2)
}
aud_split_targets <- function(c) {
  lapply(seq_along(c$tx), function(k) list(x = c$x, y = c$y, tx = c$tx[k], ty = c$ty[k]))
}

aud_pm_scale <- function(pop_adj, n, dL, dU) {
  if (pop_adj && n > 1 && dL > 0 && dU > 0) n / (n - 1) else 1
}
aud_pairwise <- function(f, X, tg) {
  k <- ncol(X)
  M <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) M[i, j] <- f(X[, i], X[, j], tg[i], tg[j])
  M
}
aud_pm_reference <- function(dL, dU, tg, X, pop_adj) {
  s <- aud_pm_scale(pop_adj, nrow(X), dL, dU)
  list(
    cupm = s * aud_pairwise(function(a, b, u, v) Co.UPM(dU, a, b, u, v), X, tg),
    dupm = s * aud_pairwise(function(a, b, u, v) D.UPM(dL, dU, a, b, u, v), X, tg),
    dlpm = s * aud_pairwise(function(a, b, u, v) D.LPM(dL, dU, a, b, u, v), X, tg),
    clpm = s * aud_pairwise(function(a, b, u, v) Co.LPM(dL, a, b, u, v), X, tg)
  )
}
aud_pm_quadrants <- function(...) lapply(PM.matrix(...)[c("cupm", "dupm", "dlpm", "clpm")], unname)
aud_quadrant_names <- c("cupm", "dupm", "dlpm", "clpm")
aud_cell_sum <- function(P) P$cupm + P$dupm + P$dlpm + P$clpm
aud_cov_tol <- function(X) 4 * (nrow(X) + 4) * aud_eps * outer(apply(X, 2, sd), apply(X, 2, sd))

aud_nd_sum <- function(X, t, d, norm) {
  Co.LPM_nD(X, t, d, norm) + Co.UPM_nD(X, t, d, norm) + DPM_nD(X, t, d, norm)
}
aud_nd_oracle <- function(X, t, d) {
  X <- matrix(as.numeric(X), nrow(X), ncol(X))
  D <- sweep(X, 2, t)
  rowprod <- function(M) apply(M, 1, prod)
  if (d == 0) {
    lo <- apply(D <= 0, 1, all)
    hi <- apply(D > 0, 1, all)
    return(c(clpm = mean(lo), cupm = mean(hi), dpm = mean(!lo & !hi)))
  }
  strict <- apply(D < 0, 1, all) | apply(D > 0, 1, all)
  c(clpm = mean(rowprod(pmax(-D, 0)^d)), cupm = mean(rowprod(pmax(D, 0)^d)),
    dpm = mean(ifelse(strict, 0, rowprod(abs(D)^d))))
}
aud_nd_triple <- function(X, t, d) {
  c(clpm = Co.LPM_nD(X, t, d, FALSE), cupm = Co.UPM_nD(X, t, d, FALSE), dpm = DPM_nD(X, t, d, FALSE))
}
aud_nd_tol <- function(X, t, d, k = 4) {
  X <- matrix(as.numeric(X), nrow(X), ncol(X))
  k * (nrow(X) + d * ncol(X) + 2) * aud_eps * mean(apply(abs(sweep(X, 2, t))^d, 1, prod))
}
aud_nd_cases <- function(f) {
  aud_over(aud_fx_co()$all_cases, function(c) vapply(0:2, function(d) f(c, d), numeric(1)))
}

# Bivariate cases and PM.matrix samples. One seed for the whole sequence of
# draws (z_norm, the random cases, X_random), in the order they were made.
aud_fx_co <- aud_memo(function() withr::with_seed(20261108, {
  # exact regime: n = 8 dyadic values, ties at the targets (0.5, 0.5)
  x_exact <- c(-2, -1, -0.5, 0, 0.5, 0.5, 1.5, 3)
  y_exact <- c(1, -1, 0.5, 0.5, -2, 0.5, 2, -0.25)
  exact_cases <- list(
    "exact ties" = list(x = x_exact, y = y_exact, tx = 0.5, ty = 0.5),
    "exact off"  = list(x = x_exact, y = y_exact, tx = 0.25, ty = -1.5),
    "exact out"  = list(x = x_exact, y = y_exact, tx = -4, ty = 4)
  )
  z_norm <- round(rnorm(400), 1)
  z_t3 <- round(rt(300, df = 3), 1)
  random_cases <- c(
    aud_co_cases(z_norm, round(0.6 * z_norm + rnorm(400), 1), "normal"),
    aud_co_cases(round(rexp(250, 2), 2), round(rexp(250, 1), 2), "expo"),
    aud_co_cases(z_t3, round(-0.4 * z_t3 + rt(300, df = 3), 1), "t3")
  )
  eu <- diff(log(EuStockMarkets))
  aq <- na.omit(airquality[, c("Ozone", "Temp")])
  real_cases <- c(
    aud_co_cases(as.numeric(eu[, "DAX"]), as.numeric(eu[, "FTSE"]), "eu"),
    aud_co_cases(faithful$eruptions, faithful$waiting, "faithful"),
    aud_co_cases(as.numeric(aq$Ozone), as.numeric(aq$Temp), "airquality"),
    aud_co_cases(quakes$mag, quakes$depth, "quakes"),
    aud_co_cases(mtcars$mpg, mtcars$hp, "mtcars")
  )
  sample_cases <- c(random_cases, real_cases)
  # Matrices for PM.matrix: exact (n = 8, dyadic), random with ties, real.
  X_exact <- cbind(a = x_exact, b = y_exact, c = c(0, 0.5, -1, 2, 2, -0.5, 1, 0.5))
  t_exact <- c(0.5, 0.5, 0.5)
  X_random <- cbind(u = round(rnorm(60), 1), v = round(rexp(60), 1),
                    w = round(rt(60, df = 3), 1))
  X_random <- cbind(X_random, dup = X_random[, "u"])  # rank deficient on purpose
  pm_samples <- list(
    random = X_random,
    eu     = eu,
    mtcars = as.matrix(mtcars),
    quakes = as.matrix(quakes[, c("lat", "long", "depth", "mag")])
  )
  list(
    x_exact = x_exact, y_exact = y_exact, exact_cases = exact_cases,
    z_norm = z_norm, z_t3 = z_t3, random_cases = random_cases, eu = eu, aq = aq,
    real_cases = real_cases, sample_cases = sample_cases,
    all_cases = c(exact_cases, sample_cases),
    X_exact = X_exact, t_exact = t_exact, X_random = X_random,
    pm_samples = pm_samples, pm_all = c(list(exact = X_exact), pm_samples),
    degree_pairs = list(c(0, 0), c(1, 1), c(2, 1), c(1, 0), c(0, 2), c(2, 2)),
    pairs_03 = lapply(0:15, function(k) c(k %/% 4, k %% 4)),
    pairs_02 = lapply(0:8, function(k) c(k %/% 3, k %% 3)),
    exact_self = aud_self_cases(x_exact, c(-1, 0.5, 0.25, 4)),
    sample_self = lapply(sample_cases, function(c) list(x = c$x, y = c$x, tx = c$tx, ty = c$tx)),
    eu_self_40 = local({
      x <- as.numeric(eu[, "DAX"])
      ts <- c(quantile(x, seq(0.01, 0.99, length.out = 30), names = FALSE), x[1:10])
      list(x = x, y = x, tx = ts, ty = ts)
    }),
    eu_vec_40 = local({
      x <- as.numeric(eu[, "CAC"]); y <- as.numeric(eu[, "SMI"])
      list(x = x, y = y,
           tx = c(quantile(x, seq(0.02, 0.98, length.out = 30), names = FALSE), x[1:10]),
           ty = c(quantile(y, seq(0.98, 0.02, length.out = 30), names = FALSE), y[1:10]))
    }),
    X_tied = cbind(a = c(1, 2, 3, 4), b = c(5, 5, 5, 5), c = c(0, 2, -1, 3)),
    t_tied = c(2.5, 5, 1),
    nd_exact_targets = list(t_exact = t_exact, mixed = c(0.5, 1, 0.5),
                            ones = c(1, 1, 1), zeros = c(0, 0, 0)),
    marginal = list(
      CoLPM_zero_add_DUPM = list(
        lhs = function(c, d) aud_co_lpm(c, d, 0) + aud_d_upm(c, d, 0),
        rhs = function(c, d) LPM(d, c$tx, c$x), tol = function(c, d) aud_co_tol(c, d, 0)),
      CoUPM_zero_add_DLPM = list(
        lhs = function(c, d) aud_co_upm(c, d, 0) + aud_d_lpm(c, 0, d),
        rhs = function(c, d) UPM(d, c$tx, c$x), tol = function(c, d) aud_co_tol(c, d, 0)),
      CoLPM_zero_left_add_DLPM = list(
        lhs = function(c, d) aud_co_lpm(c, 0, d) + aud_d_lpm(c, d, 0),
        rhs = function(c, d) LPM(d, c$ty, c$y), tol = function(c, d) aud_co_tol(c, 0, d)),
      CoUPM_zero_left_add_DUPM = list(
        lhs = function(c, d) aud_co_upm(c, 0, d) + aud_d_upm(c, 0, d),
        rhs = function(c, d) UPM(d, c$ty, c$y), tol = function(c, d) aud_co_tol(c, 0, d))
    )
  )
}))

# ---- input contract and dominance arguments ---------------------------------

aud_x4 <- c(1, 2, 2, 3)
aud_x5 <- c(1, 2, 2, 3, 10)
aud_x_om <- c(0, 1, 4, 5)
aud_y_om <- c(0, 2, 2.5, 3)  # discrete: no FSD either way; continuous: X OMEGA Y
aud_type_msg <- "`type` must be \"discrete\" or \"continuous\"."
aud_degree_msg <- "`degree` must be 1, 2 or 3."
