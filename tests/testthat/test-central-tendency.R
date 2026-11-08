# NNS.mode and NNS.gravity.
#
# NNS.gravity averages the median, a histogram mode, the mean and the midpoint
# of the quartiles. It used floor-position order statistics (the lower middle
# value for even n) and named histogram bins by their left edges, so it sat
# below the centre of symmetric data and did not reflect with the data (#117).
#
# The continuous NNS.mode returns data values on discrete-valued data, finds a
# mode at the edge of the data, and backs the regression "mode" reducer (#114).

test_that("CT01: discrete NNS.mode beyond the int range, tied modes and the gravity tie-break", {
  withr::local_seed(20261108)
  aud_case("discrete NNS.mode handles values beyond the int range", {
    expect_equal(NNS.mode(c(3e9, 3e9, 3e9, 1, 2), discrete = TRUE), 3e9)
    expect_equal(NNS.mode(c(-5e9, -5e9, -5e9, 1, 2), discrete = TRUE), -5e9)
    expect_equal(NNS.mode(c(1, 2, 2, 3, 3, 4, 4, 5), discrete = TRUE, multi = TRUE), c(2, 3, 4))
  })
  aud_case("discrete NNS.mode keeps tied modes beyond the int range apart", {
    x <- c(3e9, 3e9, -5e9, -5e9, 1)
    expect_identical(NNS.mode(x, discrete = TRUE, multi = TRUE), c(-5e9, 3e9))
    # The tie goes to the mode nearest the gravity, -6.5e8 (the median is 1).
    expect_identical(NNS.mode(x, discrete = TRUE, multi = FALSE), 3e9)
  })
  aud_case("discrete NNS.mode returns the tied mode nearest the gravity, not their mean", {
    nearest <- function(x) {
      m <- NNS.mode(x, discrete = TRUE, multi = TRUE)
      m[which.min(abs(m - NNS.gravity(x)))]                  # ties: the lower
    }
    for (x in list(c(1, 1, 3, 3), c(1, 1, 2, 2), rep(0:19, each = 400), c(1, 3), c(1, 2))) {
      r <- NNS:::mode_class(x)
      expect_true(r %in% x)
      expect_identical(r, nearest(x))
    }
    expect_identical(NNS:::mode_class(c(1, 1, 3, 3)), 1)     # 2 before
    expect_identical(NNS.mode(rep(0:19, each = 400), discrete = TRUE, multi = FALSE), 9)
    expect_identical(NNS.mode(c(1, 3), discrete = TRUE, multi = TRUE), c(1, 3))
  })
  aud_case("equally high mode peaks resolve to the one nearest the gravity", {
    x <- rep(0:19, each = 400)
    m <- NNS.mode(x)
    r <- NNS.mode(x, multi = FALSE)
    expect_true(r %in% m)
    expect_equal(r, m[which.min(abs(m - NNS.gravity(x)))])
    # Tied clusters: always a real peak, never a value in the gap between them
    # (averaging them gave 8.91 for this sample).
    y <- c(rep(1:3, each = 300), rep(10:12, each = 300), rep(20:22, each = 150))
    expect_true(NNS.mode(y, multi = FALSE) %in% NNS.mode(y))
  })
})

test_that("CT02: NNS.gravity is centred, reflects, moves a quarter as far as the mean for an outlier, and stays fast with a tight IQR", {
  withr::local_seed(20261108)
  aud_case("gravity is centred on symmetric data and reflects with it", {
    expect_identical(NNS.gravity(1:4), 2.5)                 # 2.25 before
    expect_identical(-NNS.gravity(-(1:4)), 2.5)             # 2.75 before
    g <- vapply(1:500, function(s) { set.seed(s); x <- rnorm(10); c(NNS.gravity(x), -NNS.gravity(-x)) }, numeric(2))
    expect_lt(abs(mean(g[1, ])), 0.05)                      # -0.176 before
    expect_lt(mean(abs(g[1, ] - g[2, ])), 0.1)              # 0.35 before
  })
  aud_case("one outlier moves gravity a quarter as far as the mean", {
    expect_equal(NNS.gravity(c(rep(0, 99), 1e6)), 2500)
    expect_identical(NNS.gravity(c(1, 5, 9)), 5)            # up to three values: the median
  })
  aud_case("NNS.gravity / NNS.mode stay fast with a tight IQR and an outlier", {
    set.seed(1)
    x <- c(5 + 1e-7 * rnorm(1000), 6)
    elapsed <- system.time({ g <- NNS.gravity(x); m <- NNS.mode(x, multi = FALSE) })[["elapsed"]]
    expect_lt(elapsed, 2)
    expect_equal(g, 5, tolerance = 1e-3)
    expect_equal(m, 5, tolerance = 1e-5)
    # Tiny histograms (<= 3 bins) are well defined; tied modes resolve to the
    # one nearest the gravity, then the lower value, not their average. The
    # gravity is 0.5, halfway (0.375 before #117), so the tie goes to 0.
    expect_equal(NNS.mode(c(0, 0, 1, 1), multi = FALSE), 0)
    expect_equal(NNS.gravity(c(0, 0, 1, 1)), 0.5)
  })
})

test_that("CT03: NNS.mode peaks: data values at the ends, a density highest at its minimum, peaks centred on observations", {
  withr::local_seed(20261108)
  aud_case("discrete-valued data return their own values, end values included", {
    expect_identical(NNS.mode(c(1, 1, 1, 5, 5, 5, 9, 9, 9)), c(1, 5, 9))   # 1, 3.33, 8 before
    expect_identical(NNS.mode(rep(0:4, each = 50)), as.numeric(0:4))       # 0.885 1.897 2.909
    expect_identical(NNS.mode(rep(0:19, each = 400)), as.numeric(0:19))    # 18 centres, no 0 or 19
    # Tied modes: the one nearest the gravity, so never a value in a gap.
    expect_identical(NNS.mode(c(1, 1, 1, 5, 5, 5, 9, 9, 9), multi = FALSE), 5)
    expect_identical(NNS.mode(c(0, 0, 1, 1), multi = FALSE), 0)
    # A mode is a value more frequent than its neighbours: two clusters, two modes.
    x <- rep(c(1, 2, 3, 7, 8, 9), c(10, 30, 10, 10, 25, 10))
    expect_identical(NNS.mode(x), c(2, 8))
    expect_identical(NNS.mode(x, multi = FALSE), 2)
    # Values need not be integers.
    expect_identical(NNS.mode(rep(c(0.1, 0.25, 0.4), c(20, 50, 20))), 0.25)
  })
  aud_case("a density highest at its minimum has its mode there", {
    e <- vapply(1:20, function(s) { set.seed(s); NNS.mode(rexp(1000), multi = FALSE) }, numeric(1))
    expect_lt(mean(e), 0.3)                                                # 0.64 before
    # Isolated extremes of a symmetric sample are not modes.
    k <- vapply(1:20, function(s) { set.seed(s); x <- rnorm(100); m <- NNS.mode(x); sum(m == min(x) | m == max(x)) }, numeric(1))
    expect_identical(sum(k), 0)
  })
  aud_case("peaks are centred on the observations, not on bin edges", {
    # A sharp spike at 0.5 on a uniform background: the left-edge bin names put
    # the mode up to a bin (about 0.11) below it (0.418 to 0.493 before).
    m <- vapply(1:10, function(s) {
      set.seed(s); NNS.mode(c(0.5 + rnorm(1000, sd = 0.001), runif(1000, 0, 10)), multi = FALSE)
    }, numeric(1))
    expect_lt(max(abs(m - 0.5)), 0.01)
  })
})
