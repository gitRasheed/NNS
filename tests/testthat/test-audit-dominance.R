# Stochastic dominance audit: the `type` argument, non-finite observations and
# the plot curves of the bidirectional tests.

x_om <- c(0, 1, 4, 5)
y_om <- c(0, 2, 2.5, 3)  # discrete: no FSD either way; continuous: X OMEGA Y
type_msg <- "`type` must be \"discrete\" or \"continuous\"."

# ---- type (#23) -------------------------------------------------------------

test_that("NNS.FSD rejects an invalid type instead of running the Omega order", {
  expect_error(NNS.FSD(x_om, y_om, type = "foo", plot = FALSE), type_msg, fixed = TRUE)
})

test_that("NNS.FSD.uni rejects an invalid type instead of running discrete FSD", {
  expect_error(NNS.FSD.uni(x_om, y_om, type = "discret"), type_msg, fixed = TRUE)
})

test_that("NNS.SD.efficient.set rejects an invalid type", {
  expect_error(NNS.SD.efficient.set(cbind(x = x_om, y = y_om), 1, type = "foo", status = FALSE),
               type_msg, fixed = TRUE)
})

test_that("NNS.SD.efficient.set rejects an invalid type at degree 2, where type is unused", {
  expect_error(NNS.SD.efficient.set(cbind(x = x_om, y = y_om), 2, type = "foo", status = FALSE),
               type_msg, fixed = TRUE)
})

test_that("NNS.SD.cluster rejects an invalid type", {
  expect_error(NNS.SD.cluster(cbind(x = x_om, y = y_om), 1, type = "foo"), type_msg, fixed = TRUE)
})

test_that("type must be a single non-missing string", {
  expect_error(NNS.FSD(x_om, y_om, type = NA, plot = FALSE), type_msg, fixed = TRUE)
  expect_error(NNS.FSD(x_om, y_om, type = c("discrete", "continuous"), plot = FALSE), type_msg, fixed = TRUE)
  expect_error(NNS.FSD.uni(x_om, y_om, type = 1), type_msg, fixed = TRUE)
  expect_error(NNS.FSD.uni(x_om, y_om, type = "Discrete "), type_msg, fixed = TRUE)
})

test_that("the C++ entry points reject an invalid type", {
  expect_error(.Call(NNS:::`_NNS_NNS_FSD_uni_cpp`, x_om, y_om, "foo"), type_msg, fixed = TRUE)
  expect_error(NNS:::sd_dom_matrix_prefix_parallel(cbind(x_om, y_om), 1L, "foo"), type_msg, fixed = TRUE)
})

test_that("valid types are case-insensitive and keep their verdicts", {
  expect_identical(NNS.FSD(x_om, y_om, type = "Continuous", plot = FALSE), "X OMEGA Y")
  expect_identical(NNS.FSD(x_om, y_om, type = "DISCRETE", plot = FALSE), "NO FSD EXISTS")
  expect_identical(NNS.FSD.uni(x_om, y_om, "CONTINUOUS"), 1L)
  expect_identical(NNS.FSD.uni(x_om, y_om, "discrete"), 0L)
  expect_identical(NNS.SD.efficient.set(cbind(x = x_om, y = y_om), 1, "Continuous", status = FALSE), "x")
})

# ---- non-finite observations (#24) ------------------------------------------

test_that("the uni-directional tests reject an infinite observation in x", {
  expect_error(NNS.SSD.uni(c(1, Inf), c(0, 1)), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.TSD.uni(c(1, Inf), c(0, 1)), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.FSD.uni(c(1, Inf), c(0, 1)), "`x` must be finite.", fixed = TRUE)
})

test_that("the uni-directional tests reject an infinite observation in y", {
  expect_error(NNS.SSD.uni(c(0, 1), c(-Inf, 1)), "`y` must be finite.", fixed = TRUE)
  expect_error(NNS.TSD.uni(c(0, 1), c(-Inf, 1)), "`y` must be finite.", fixed = TRUE)
  expect_error(NNS.FSD.uni(c(0, 1), c(-Inf, 1), "continuous"), "`y` must be finite.", fixed = TRUE)
})

test_that("a total-loss log return makes the efficient set an error, not a superset of the FSD set", {
  # The bankrupt asset used to be FSD-dominated but SSD/TSD-efficient.
  x <- diff(log(c(100, 101, 99, 102, 103)))
  y <- diff(log(c(100, 60, 30, 10, 0)))  # ends in -Inf
  m <- cbind(x = x, y = y)
  expect_error(NNS.SD.efficient.set(m, 1, status = FALSE), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.SD.efficient.set(m, 2, status = FALSE), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.SD.efficient.set(m, 3, status = FALSE), "`x` must be finite.", fixed = TRUE)
})

test_that("NNS.SD.cluster names `data` for an infinite observation", {
  m <- cbind(a = c(1, Inf, 2), b = c(0, 1, 1), c = c(-1, 0, 0))
  expect_error(NNS.SD.cluster(m, 2), "`data` must be finite.", fixed = TRUE)
})

test_that("the dominance-matrix kernel rejects an infinite observation", {
  m <- cbind(a = c(1, Inf, 2), b = c(0, 1, 1))
  expect_error(NNS:::sd_dom_matrix_prefix_parallel(m, 2L, "discrete"), "`x` must be finite.", fixed = TRUE)
})

test_that("NA and NaN keep the missing-values message", {
  msg <- "You have some missing values, please address."
  expect_error(NNS.TSD.uni(c(1, NaN), 1:3), msg, fixed = TRUE)
  expect_error(NNS.SSD(c(1, NA, Inf), 1:3, plot = FALSE), msg, fixed = TRUE)
  expect_error(NNS.SD.efficient.set(cbind(a = c(1, NA), b = c(Inf, 2)), 2, status = FALSE), msg, fixed = TRUE)
  expect_error(NNS.SD.cluster(cbind(a = c(1, NA), b = c(Inf, 2)), 2), msg, fixed = TRUE)
})

test_that("a large finite loss is still dominated at every order", {
  x <- c(0, 1, 3); y <- c(0, 1, -50)
  expect_identical(c(NNS.FSD.uni(x, y), NNS.SSD.uni(x, y), NNS.TSD.uni(x, y)), c(1L, 1L, 1L))
})

# ---- bidirectional tests and their plot curves (#28) ------------------------

test_that("the bidirectional tests name x for an infinite observation", {
  expect_error(NNS.FSD(c(1, Inf), c(0, 1), plot = FALSE), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.FSD(c(1, Inf), c(0, 1), type = "continuous", plot = FALSE), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.SSD(c(1, Inf), c(0, 1), plot = FALSE), "`x` must be finite.", fixed = TRUE)
  expect_error(NNS.TSD(c(1, Inf), c(0, 1), plot = FALSE), "`x` must be finite.", fixed = TRUE)
})

test_that("the bidirectional tests name y for an infinite observation", {
  expect_error(NNS.FSD(c(0, 1), c(1, -Inf), plot = FALSE), "`y` must be finite.", fixed = TRUE)
  expect_error(NNS.SSD(c(0, 1), c(1, -Inf), plot = FALSE), "`y` must be finite.", fixed = TRUE)
  expect_error(NNS.TSD(c(0, 1), c(1, -Inf), plot = FALSE), "`y` must be finite.", fixed = TRUE)
})

test_that("the bidirectional tests reject infinite input the same way with plot = TRUE", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(NNS.SSD(c(1, Inf), c(0, 1), plot = TRUE), "`x` must be finite.", fixed = TRUE)
})

test_that("plot = FALSE does not compute the plot curves", {
  local_mocked_bindings(LPM = function(...) stop("plot curve computed"),
                        LPM.ratio = function(...) stop("plot curve computed"))
  expect_identical(NNS.FSD(x_om, y_om, plot = FALSE), "NO FSD EXISTS")
  expect_identical(NNS.FSD(x_om, y_om, type = "continuous", plot = FALSE), "X OMEGA Y")
  expect_identical(NNS.SSD(c(1, 1), c(0, 2), plot = FALSE), "X SSD Y")
  expect_identical(NNS.TSD(c(0, 2), c(1, 1), plot = FALSE), "Y TSD X")
})

test_that("the verdict does not depend on plot", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(28)
  x <- round(rnorm(30), 2); y <- round(rnorm(25, 0.3), 2)
  expect_identical(NNS.FSD(x, y, type = "continuous", plot = TRUE), NNS.FSD(x, y, type = "continuous", plot = FALSE))
  expect_identical(NNS.SSD(x, y, plot = TRUE), NNS.SSD(x, y, plot = FALSE))
  expect_identical(NNS.TSD(x, y, plot = TRUE), NNS.TSD(x, y, plot = FALSE))
})

# ---- TSD plot axis (#27) ----------------------------------------------------

test_that("NNS.TSD plots against the thresholds, as NNS.SSD does", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- c(0, 0.01, 0.02, 10); y <- c(1, 2, 3, 4)
  NNS.TSD(x, y, plot = TRUE)
  usr <- graphics::par("usr")[1:2]
  # the x-axis spans the data range [0, 10], not the index range [1, 8]
  expect_lte(usr[1], 0)
  expect_gte(usr[2], 10)
})

# ---- degree of the efficient set and clusters ---------------------------------

degree_msg <- "`degree` must be 1, 2 or 3."

test_that("NNS.SD.efficient.set rejects a fractional degree instead of truncating it", {
  m <- cbind(x = x_om, y = y_om)
  for (d in list(1.5, 2.9, 0.5, 3.0000001)) {
    expect_error(NNS.SD.efficient.set(m, d, status = FALSE), degree_msg, fixed = TRUE)
  }
})

test_that("NNS.SD.efficient.set rejects other invalid degrees with the same message", {
  m <- cbind(x = x_om, y = y_om)
  for (d in list(0, 4, -1, NA, NA_real_, NaN, Inf, c(1, 2), numeric(0), "2", TRUE, factor(2), NULL)) {
    expect_error(NNS.SD.efficient.set(m, d, status = FALSE), degree_msg, fixed = TRUE)
  }
  # Beyond the integer range: an error, not a coercion warning.
  expect_no_warning(expect_error(NNS.SD.efficient.set(m, 1e10, status = FALSE), degree_msg, fixed = TRUE))
})

test_that("NNS.SD.efficient.set accepts 1, 2, 3 as double or integer with unchanged results", {
  set.seed(7)
  A <- cbind(a = rnorm(40), b = rnorm(40, 0.3), c = rexp(40) - 1)
  for (d in 1:3) {
    ref <- .Call(NNS:::`_NNS_NNS_SD_efficient_set_parallel_cpp`, A, d, "discrete", FALSE)
    expect_identical(NNS.SD.efficient.set(A, d, status = FALSE), ref)
    expect_identical(NNS.SD.efficient.set(A, as.numeric(d), status = FALSE), ref)
  }
})

test_that("NNS.SD.cluster rejects an invalid degree even when no efficient set is computed", {
  m <- cbind(x = x_om, y = y_om)
  expect_error(NNS.SD.cluster(m, 1.5), degree_msg, fixed = TRUE)
  # min_cluster >= ncol: the efficient set is never called.
  expect_error(NNS.SD.cluster(m, 1.5, min_cluster = 2), degree_msg, fixed = TRUE)
  expect_identical(NNS.SD.cluster(m, 1, min_cluster = 2)$Clusters,
                   NNS.SD.cluster(m, 1L, min_cluster = 2)$Clusters)
})

test_that("the compiled efficient set checks the degree before its LPM ordering pass", {
  m <- cbind(x = x_om, y = y_om)
  # Before the fix, the ordering pass ran 2147483647 multiplications per
  # element before the degree was checked.
  elapsed <- system.time(
    expect_error(.Call(NNS:::`_NNS_NNS_SD_efficient_set_parallel_cpp`, m, .Machine$integer.max,
                       "discrete", FALSE),
                 "degree must be 1, 2, or 3", fixed = TRUE)
  )[["elapsed"]]
  expect_lt(elapsed, 2)
  expect_error(.Call(NNS:::`_NNS_NNS_SD_efficient_set_parallel_cpp`, m[, 0, drop = FALSE], 4L,
                     "discrete", FALSE),
               "degree must be 1, 2, or 3", fixed = TRUE)
})
