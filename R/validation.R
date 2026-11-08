# Shared input validation for the partial-moment family.
#
# A factor is an integer vector underneath, but its level codes are not its
# values, so it is rejected rather than silently coerced. Other numeric data
# (double, integer, logical, and classes built on them such as Date or xts) is
# accepted; character, list and other types are rejected.

.nns_stop_arg <- function(msg, call) {
  stop(simpleError(msg, call = call))
}

# `null_ok` only changes the message, for arguments where NULL has a meaning
# (resolved by the caller before this is reached).
.nns_as_numeric <- function(x, name, call = sys.call(-1L), null_ok = FALSE) {
  if (is.factor(x)) {
    .nns_stop_arg(sprintf("`%s` must be numeric, not a factor.", name), call)
  }
  if (!typeof(x) %in% c("double", "integer", "logical")) {
    .nns_stop_arg(sprintf("`%s` must be numeric%s.", name, if (null_ok) " or NULL" else ""), call)
  }
  as.numeric(x)
}

# Data frames and matrices go through as.matrix()/unlist(), which turn a
# factor column into strings or level codes; reject it while the columns are
# still visible.
.nns_stop_if_factor_columns <- function(x, name, call = sys.call(-1L)) {
  if (is.factor(x) || (is.data.frame(x) && any(vapply(x, is.factor, logical(1L))))) {
    .nns_stop_arg(sprintf("`%s` must be numeric, not a factor.", name), call)
  }
  invisible(x)
}

# A single TRUE/FALSE flag. Numeric 0/1 keep working; NA, vectors and other
# types are an error instead of depending on how C++ or `&&` coerce them.
.nns_as_flag <- function(x, name, call = sys.call(-1L)) {
  if (!(is.logical(x) || is.numeric(x)) || is.factor(x) ||
      length(x) != 1L || is.na(x)) {
    .nns_stop_arg(sprintf("`%s` must be TRUE or FALSE.", name), call)
  }
  as.logical(x)
}

# Single-column data (vector, or a one-column data.frame/matrix/tbl/
# data.table/xts/zoo) as a plain numeric vector; used by the stochastic
# dominance functions.
.nns_single_column_numeric <- function(v, arg_name, call = sys.call(-1L)) {
  if (inherits(v, c("tbl", "data.table", "xts", "zoo")) || is.data.frame(v) || is.matrix(v)) {
    if (!is.null(dim(v)) && ncol(v) > 1) {
      .nns_stop_arg(sprintf("%s must be a single-column object or numeric vector.", arg_name), call)
    }
    v <- unlist(v, use.names = FALSE)
    if (!is.factor(v)) v <- as.vector(v)
  }
  .nns_as_numeric(v, arg_name, call = call)
}

# The `type` of the stochastic dominance functions: one of "discrete" or
# "continuous" in any letter case, returned lower-cased. Anything else is an
# error rather than a silent fallback, which differed between entry points.
.nns_sd_type <- function(type, call = sys.call(-1L)) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !tolower(type) %in% c("discrete", "continuous")) {
    .nns_stop_arg("`type` must be \"discrete\" or \"continuous\".", call)
  }
  tolower(type)
}

# One sample of a stochastic dominance comparison: single-column numeric data
# without missing or infinite values. The tolerances and the midrange shift of
# the dominance tests are defined for finite data only.
.nns_sd_sample <- function(v, name, call = sys.call(-1L)) {
  v <- .nns_single_column_numeric(v, name, call = call)
  if (anyNA(v)) .nns_stop_arg("You have some missing values, please address.", call)
  if (!all(is.finite(v))) .nns_stop_arg(sprintf("`%s` must be finite.", name), call)
  v
}
