#' LPM VaR
#'
#' Generates a value at risk (VaR) quantile based on the Lower Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for left-tail VaR. Values outside [0, 1] are clamped; \code{NA} is an error.
#' @param degree integer; \code{(degree = 0)} gives the interpolated empirical quantile (\code{quantile(x, percentile, type = 7)}), which on tied data need not be an exact inverse of \code{LPM.ratio(0, ., x)}. Use it for probability statements, such as a VaR exceeded with a given frequency.  \code{(degree = 1)} gives an area-weighted (severity-weighted) threshold, the \code{percentile}-expectile of \code{x} (Newey and Powell 1987): \code{percentile} of the area of deviations, not of the observations, lies below it, so it is not a quantile of \code{x} (for exponential data the degree-1 5\% point has about a quarter of the observations below it).  Higher degrees weight large deviations more.
#' @param x a numeric vector. \code{NA} and \code{NaN} are dropped; infinite values, factors and character data are an error. A data frame (including a tibble or data.table) is flattened to one vector.
#' @return Returns the type-7 empirical quantile (\code{degree = 0}), whose strict or inclusive sample tail fraction may differ from \code{percentile} because of interpolation or ties, or the threshold with \code{percentile} of the directional deviation area below it (\code{degree > 0}). For integer \code{degree} 1 to 4 the threshold is solved for itself, not to a tolerance on the ratio: it is the root to rounding error on the scale of the data, also for percentiles near 0 or 1, where the ratio is flat. A root within rounding of an order statistic returns that order statistic (for example \code{min(x)} at a tiny \code{percentile} when \code{min(x)} is not 0); where the offset from it is representable, as on data with minimum 0, the root is resolved to relative precision for percentiles down to the smallest normal double, about 2.2e-308.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, left-tail
#' LPM.VaR(0.05, 0, x)
#' }
#' @export

LPM.VaR <- function(percentile, degree, x) {
  x <- .NNS_prepare_VaR_x(x)
  percentile <- .NNS_prepare_VaR_percentile(percentile)
  
  if (degree == 0) {
    return(stats::quantile(x, percentile, na.rm = TRUE))
  }
  
  if (.NNS_is_supported_integer_degree(degree)) {
    return(.NNS_LPM_VaR_integer(percentile, as.integer(degree), x))
  }
  
  .NNS_LPM_VaR_optimize(percentile, degree, x)
}



#' UPM VaR
#'
#' Generates an upside value at risk (VaR) quantile based on the Upper Partial Moment ratio.
#'
#' @param percentile numeric [0, 1]; The percentile for right-tail VaR. Values outside [0, 1] are clamped; \code{NA} is an error.
#' @param degree integer; \code{(degree = 0)} gives the interpolated empirical quantile (\code{quantile(x, 1 - percentile, type = 7)}), which on tied data need not be an exact inverse of \code{UPM.ratio(0, ., x)}. Use it for probability statements.  \code{(degree = 1)} gives an area-weighted (severity-weighted) threshold, the (1 - \code{percentile})-expectile of \code{x} (Newey and Powell 1987): \code{percentile} of the area of deviations, not of the observations, lies above it, so it is not a quantile of \code{x}.  Higher degrees weight large deviations more.
#' @param x a numeric vector. \code{NA} and \code{NaN} are dropped; infinite values, factors and character data are an error. A data frame (including a tibble or data.table) is flattened to one vector.
#' @return Returns the type-7 empirical quantile (\code{degree = 0}), whose strict or inclusive sample tail fraction may differ from \code{percentile} because of interpolation or ties, or the threshold with \code{percentile} of the directional deviation area above it (\code{degree > 0}). For integer \code{degree} 1 to 4 the threshold is solved for itself, not to a tolerance on the ratio: it is the root to rounding error on the scale of the data, also for percentiles near 0 or 1, where the ratio is flat. A root within rounding of an order statistic returns that order statistic (for example \code{max(x)} at a tiny \code{percentile} when \code{max(x)} is not 0); where the offset from it is representable, as on data with maximum 0, the root is resolved to relative precision for percentiles down to the smallest normal double, about 2.2e-308.
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995, 2nd edition: \url{https://ovvo-financial.github.io/NNS/book/})
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' ## For 5th percentile, right-tail
#' UPM.VaR(0.05, 0, x)
#' @export

UPM.VaR <- function(percentile, degree, x) {
  x <- .NNS_prepare_VaR_x(x)
  percentile <- .NNS_prepare_VaR_percentile(percentile)
  
  if (degree == 0) {
    return(stats::quantile(x, 1 - percentile, na.rm = TRUE))
  }
  
  if (.NNS_is_supported_integer_degree(degree)) {
    return(.NNS_UPM_VaR_integer(percentile, as.integer(degree), x))
  }
  
  .NNS_UPM_VaR_optimize(percentile, degree, x)
}



# ==============================================================================
# Internal helpers
# ==============================================================================

.NNS_prepare_VaR_x <- function(x) {
  # Any data frame (tbl and data.table included) is flattened, as NNS.CDF
  # and the stochastic dominance functions do for one column.
  if (is.data.frame(x)) {
    .nns_stop_if_factor_columns(x, "x", call = sys.call(-1L))
    x <- unlist(x, use.names = FALSE)
  }
  
  x <- .nns_as_numeric(x, "x", call = sys.call(-1L))
  x <- x[!is.na(x)]
  
  if (!length(x)) {
    stop("x must contain at least one non-NA numeric value.")
  }
  # NA/NaN are dropped above; +/-Inf would reach the degree > 0 solvers as
  # Inf - Inf, so every degree rejects it the same way.
  if (!all(is.finite(x))) .nns_stop_arg("`x` must be finite.", sys.call(-1L))
  
  x
}



# Percentiles outside [0, 1] are clamped; a missing one has no threshold and
# is an error at every degree, as a missing target is for the ratios.
.NNS_prepare_VaR_percentile <- function(percentile) {
  percentile <- .nns_as_numeric(percentile, "percentile", call = sys.call(-1L))
  if (anyNA(percentile)) .nns_stop_arg("`percentile` must not be NA.", sys.call(-1L))
  pmin(pmax(percentile, 0), 1)
}



.NNS_is_supported_integer_degree <- function(degree) {
  # Range first: as.integer() warns on a degree beyond the integer range.
  length(degree) == 1L &&
    is.finite(degree) &&
    degree >= 1 &&
    degree <= 4 &&
    degree == as.integer(degree)
}



# ------------------------------------------------------------------------------
# Integer-degree VaR inversion for degrees 1:4
# ------------------------------------------------------------------------------
#
# For integer degree d,
#
#   LPM_d(t) = mean((t - x)^d [x <= t]),   UPM_d(t) = mean((x - t)^d [x > t]).
#
# Between adjacent distinct observations u_k < u_{k+1} (gap g) the below/above
# sets are fixed. With s = t - u_k and r = u_{k+1} - t = g - s,
#
#   LPM_d(t) = sum_j C(d,j) s^(d-j) lower_j(k),  lower_j(k) = mean((u_k - x)^j [x <= u_k])
#   UPM_d(t) = sum_j C(d,j) r^(d-j) upper_j(k),  upper_j(k) = mean((x - u_{k+1})^j [x >= u_{k+1}])
#
# Each sum is anchored at its own adjacent order statistic, so every term is
# >= 0 and nothing cancels (an expansion around one fixed centre cancels from
# degree 3 up when a cluster sits far from that centre). LPM.ratio(t) = p reads
#
#   h(t) = w_L * LPM_d(t) - w_U * UPM_d(t) = 0,   w_L = 1 - p, w_U = p,
#
# which is strictly increasing in t on [min x, max x] for non-constant data, so
# the root exists and is unique. UPM.VaR passes the weights the other way round
# (w_L = q, w_U = 1 - q), so a tail weight is never formed as 1 - (1 - q).
#
# The solver
#   1. locates the gap by the sign of h at the breaks (binary search);
#   2. solves h = 0 on that gap with uniroot() from the nearer anchor (in s on
#      the lower half, in r on the upper half), so the root is resolved to a
#      few ulps of t, also where the ratio is flat near min x and max x; a
#      root very close to the anchor is bracketed by a bound on its offset, so
#      tail percentiles do not cost hundreds of bisection steps;
#   3. returns a break only where the sign of h puts the root there. A tolerance
#      on the ratio is not a tolerance on t: the former snap to a break within
#      2^-39 of its ratio gave LPM.VaR(1e-12, 4, c(0, 1)) = 0, root about 1e-3.
# ------------------------------------------------------------------------------

.NNS_LPM_VaR_integer <- function(percentile, degree, x) {
  p <- pmin(pmax(as.numeric(percentile), 0), 1)
  .NNS_VaR_integer_solve(
    w_L = 1 - p, w_U = p,
    at_min = p <= 0, at_max = p >= 1,
    degree = degree, x = x
  )
}



.NNS_UPM_VaR_integer <- function(percentile, degree, x) {
  q <- pmin(pmax(as.numeric(percentile), 0), 1)
  # UPM.ratio(t) = q is LPM.ratio(t) = 1 - q, with the weights kept exact.
  .NNS_VaR_integer_solve(
    w_L = q, w_U = 1 - q,
    at_min = q >= 1, at_max = q <= 0,
    degree = degree, x = x
  )
}



.NNS_VaR_integer_solve <- function(w_L, w_U, at_min, at_max, degree, x) {
  n <- length(x)
  x_sorted <- sort(x)
  x_min <- x_sorted[1L]
  x_max <- x_sorted[n]
  
  if (n == 1L || x_min == x_max) {
    return(rep(x_min, length(w_L)))
  }
  
  # The ratio is translation- and scale-free, but deviations of extreme
  # magnitude under- or overflow in the degree-d powers: build the anchored sums
  # on data centred at the median and scaled by an exact power of two (1 for
  # ordinary data). Answers are mapped back from the original order statistics.
  z <- x_sorted - x_sorted[(n + 1L) %/% 2L]
  if (!is.finite(z[1L]) || !is.finite(z[n])) {
    # The spread from the median overflows (data near both ends of the double
    # range): centre at the midrange, formed without overflow, instead.
    z <- x_sorted - (x_min / 2 + x_max / 2)
  }
  scale <- .nns_pow2_rescale(max(-z[1L], z[n]))
  z <- z * scale
  
  runs <- rle(z)
  u <- runs$values
  m <- length(u)
  k_break <- cumsum(runs$lengths)
  u_orig <- x_sorted[k_break]
  
  out <- numeric(length(w_L))
  out[at_min] <- x_min
  out[at_max] <- x_max
  middle <- !(at_min | at_max)
  
  if (!any(middle)) return(out)
  
  gap <- diff(u)
  d1 <- degree + 1L
  bin <- choose(degree, 0:degree)
  
  # lower[k, j + 1] = lower_j(k), upper[k, j + 1] = upper_j(k), k = 1, ..., m - 1
  # (src/partial_moments_rcpp.cpp, accumulated in long double).
  sums <- .Call("_NNS_VaR_anchor_sums_RCPP", u, as.numeric(runs$lengths), degree)
  lower <- sums$lower
  upper <- sums$upper
  
  # sum_j C(d,j) h^(d-j) coef_j for each row of coef; all terms are >= 0.
  combine <- function(coef, h) {
    value <- 0
    for (j in 0:degree) {
      value <- value + bin[j + 1L] * h^(degree - j) * coef[, j + 1L]
    }
    value
  }
  
  # LPM and UPM at the left end of every gap (t = u_k), in the solver's arithmetic.
  lpm_break <- lower[, d1]
  upm_break <- combine(upper, gap)
  
  # Offset o in [0, half] from the anchor with f(o) = 0, given f(0) = f0 and
  # f(half) = f_half of opposite signs. In o the near moment is
  # sum_j C(d,j) o^(d-j) coef_j; its j < d terms are >= 0 and at the root they
  # sum to at most excess = |f0| / w (the far moment only shrinks as o grows),
  # so each term bounds o (the j = 0 term always does: coef_0 is the fraction
  # of the sample at or beyond the anchor). The tightest bound is within about
  # a factor d of the root when the root is close to the anchor, where
  # [0, half] can be hundreds of binades wide and zeroin would bisect through
  # them (about 1000 iterations at d = 2 for a percentile of 1e-307). The
  # narrow bracket is used only where it saves more than 30 binades, so other
  # roots keep the [0, half] bracket. The bracket starts at twice the bound,
  # where h is well away from 0 (at d = 1 in the tails the bound is the root
  # itself, up to rounding), and is doubled while rounding leaves it short of the root.
  solve_offset <- function(f, coef, excess, f0, half, f_half) {
    j <- 0:(degree - 1L)
    term <- bin[j + 1L] * coef[1L, j + 1L]
    pos <- term > 0
    hi <- 2 * min(exp((log(excess) - log(term[pos])) / (degree - j[pos])))
    f_hi <- f_half
    if (is.finite(hi) && hi > 0 && hi < 2^-30 * half) {
      repeat {
        f_hi <- f(hi)
        if (f_hi == 0) return(hi)
        if (sign(f_hi) != sign(f0)) break
        hi <- 2 * hi
        if (hi >= half) {
          hi <- half
          f_hi <- f_half
          break
        }
      }
    } else {
      hi <- half
    }
    o <- stats::uniroot(
      f, c(0, hi), f.lower = f0, f.upper = f_hi,
      tol = .Machine$double.xmin, maxiter = 2000L
    )$root
    min(max(o, 0), hi)
  }
  
  solve_one <- function(wl, wu) {
    # h(u_1) = -wu * UPM(u_1) < 0 and h(u_m) = wl * LPM(u_m) > 0: find k with
    # h(u_k) <= 0 < h(u_{k+1}) (h(u_m) is taken as positive).
    lo <- 1L
    hi <- m
    while (hi - lo > 1L) {
      mid <- (lo + hi) %/% 2L
      if (wl * lpm_break[mid] - wu * upm_break[mid] <= 0) lo <- mid else hi <- mid
    }
    k <- lo
    g <- gap[k]
    lo_coef <- lower[k, , drop = FALSE]
    up_coef <- upper[k, , drop = FALSE]
    
    h_s <- function(s) wl * combine(lo_coef, s) - wu * combine(up_coef, g - s)
    h_r <- function(r) wl * combine(lo_coef, g - r) - wu * combine(up_coef, r)
    
    h_left <- wl * lpm_break[k] - wu * upm_break[k]
    if (h_left >= 0) return(u_orig[k])
    h_right <- h_r(0)
    if (h_right <= 0) return(u_orig[k + 1L])
    
    half <- g / 2
    h_half <- h_s(half)
    if (h_half == 0) return(u_orig[k] + half / scale)
    
    # Solve from the nearer anchor; zeroin stops at about 4 * eps * |s| (or |r|),
    # so the offset from the anchor is resolved to relative precision.
    if (h_half > 0) {
      s <- solve_offset(h_s, lo_coef, -h_left / wl, h_left, half, h_half)
      t <- u_orig[k] + s / scale
    } else {
      r <- solve_offset(h_r, up_coef, h_right / wu, h_right, half, h_half)
      t <- u_orig[k + 1L] - r / scale
    }
    min(max(t, u_orig[k]), u_orig[k + 1L])
  }
  
  idx <- which(middle)
  out[idx] <- vapply(idx, function(i) solve_one(w_L[i], w_U[i]), numeric(1))
  out
}



# ------------------------------------------------------------------------------
# Old-method fallback for unsupported degrees
# ------------------------------------------------------------------------------

.NNS_LPM_VaR_optimize <- function(percentile, degree, x) {
  p <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  x_min <- min(x)
  x_max <- max(x)
  
  if (x_min == x_max) {
    return(rep(x_min, length(p)))
  }
  
  # The ratio is monotone in the target: solve ratio(b) = p by bracketing, to a
  # tolerance relative to the range of x.
  tol <- .Machine$double.eps^0.5 * (x_max - x_min)
  vapply(
    p,
    function(pp) {
      g <- function(b) as.numeric(.Call("_NNS_LPM_ratio_RCPP", degree, b, x)) - pp
      g_lo <- g(x_min); g_hi <- g(x_max)
      if (g_lo == 0) return(x_min)
      if (g_hi == 0) return(x_max)
      if (sign(g_lo) == sign(g_hi)) {
        return(if (abs(g_lo) <= abs(g_hi)) x_min else x_max)
      }
      stats::uniroot(g, c(x_min, x_max), f.lower = g_lo, f.upper = g_hi, tol = tol)$root
    },
    numeric(1)
  )
}



.NNS_UPM_VaR_optimize <- function(percentile, degree, x) {
  p <- pmin(pmax(as.numeric(percentile), 0), 1)
  
  x_min <- min(x)
  x_max <- max(x)
  
  if (x_min == x_max) {
    return(rep(x_min, length(p)))
  }
  
  # The ratio is monotone in the target: solve ratio(b) = p by bracketing, to a
  # tolerance relative to the range of x.
  tol <- .Machine$double.eps^0.5 * (x_max - x_min)
  vapply(
    p,
    function(pp) {
      g <- function(b) as.numeric(.Call("_NNS_UPM_ratio_RCPP", degree, b, x)) - pp
      g_lo <- g(x_min); g_hi <- g(x_max)
      if (g_lo == 0) return(x_min)
      if (g_hi == 0) return(x_max)
      if (sign(g_lo) == sign(g_hi)) {
        return(if (abs(g_lo) <= abs(g_hi)) x_min else x_max)
      }
      stats::uniroot(g, c(x_min, x_max), f.lower = g_lo, f.upper = g_hi, tol = tol)$root
    },
    numeric(1)
  )
}