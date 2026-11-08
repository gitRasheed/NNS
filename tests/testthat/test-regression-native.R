# Native multivariate regression kernels and distances against their R references.

test_that("MREG01: NNS.distance exponential rank weights use rate 1/k on every path and thread count", {
  withr::local_seed(20261108)
  aud_case("NNS.distance weights ranks by dexp(rank, rate = 1/k)", {
    set.seed(3)
    n <- 300
    X <- cbind(a = runif(n), b = runif(n))
    y <- sin(4 * X[, 1]) + X[, 2]^2 + rnorm(n, sd = .1)
    rpm <- NNS.reg(X, y, plot = FALSE, order = 2)$RPM
    pt <- c(a = .37, b = .61)

    mins <- pmin(vapply(rpm[, 1:2], min, 0), pt)
    maxs <- pmax(vapply(rpm[, 1:2], max, 0), pt)
    Z <- sweep(sweep(as.matrix(rpm[, 1:2]), 2, pt), 2, maxs - mins, "/")
    d <- rowSums(abs(Z) + Z^2)
    norm <- function(v) if (all(is.finite(v)) && sum(v) > 0) v / sum(v) else rep(0, length(v))
    reference <- function(k) {
      ids <- order(d)[seq_len(k)]
      s <- d[ids]; r <- seq_len(k)
      w <- norm(rep(1 / k, k) + norm(dt(s, k)) + norm(1 / s) +
        norm(dexp(r, rate = 1 / k)) + norm(abs(rev(dlnorm(r, 0, sd(r), log = TRUE)))) +
        norm(r^-2) + norm(dnorm(s, 0, sd(s))) + norm(exp(-s / (2 * var(s)))))
      sum(w * rpm$y.hat[ids])
    }

    for (k in c(3, 5, 10)) {
      expect_equal(NNS.distance(rpm, pt, k = k), reference(min(k, nrow(rpm))), tolerance = 1e-10)
    }
  })
  aud_case("corrected exponential ranks agree across distance paths and threads", {
    set.seed(133)
    rpm <- data.frame(a = runif(18), b = runif(18), y.hat = rnorm(18))
    query <- data.frame(a = c(.23, .41, .71), b = c(.37, .63, .19))
    for (cores in c(1L, 2L)) {
      path <- NNS:::NNS.distance.path.bulk(rpm, query, kmax = 8L, ncores = cores)
      for (k in c(1L, 3L, 5L, 8L)) {
        scalar <- vapply(seq_len(nrow(query)), function(i) {
          NNS.distance(rpm, query[i, ], k = k)
        }, numeric(1))
        one <- NNS:::NNS.distance.path.single.bulk(rpm, query, k, ncores = cores)
        bulk <- NNS:::NNS.distance.bulk(rpm, query, k)
        expect_equal(as.numeric(path[, k]), scalar, tolerance = 1e-12)
        expect_equal(as.numeric(one), scalar, tolerance = 1e-12)
        expect_equal(as.numeric(bulk), scalar, tolerance = 1e-12)
      }
    }
  })
})

test_that("MREG02: native prediction kernels equal their references: 1-D fast path, partial selection, all-k path, xstar, univariate helper, thread determinism", {
  withr::local_seed(20261108)
  aud_case("1-D fast path matches generic distances with unequal ranges", {
    x <- c(-2, -1, 0, 0, 1, 2)
    rpm_x <- cbind(x, x)
    y_reg <- c(4, 1, 2, 6, 1, 4)
    y_cls <- c(1, 1, 2, 2, 3, 3)
    xt_scalar <- c(-0.5, 0, 0.5, 3)
    Xtest <- cbind(xt_scalar, xt_scalar)

    # Deliberately inconsistent supplied ranges for duplicated coordinates.
    # The 1-D shortcut must still reproduce the generic per-column scaling.
    mins <- c(-2, -4)
    maxs <- c(2, 8)

    rpm_reg <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_reg)
    rpm_cls <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_cls)
    dist_names <- c("NNS", "L2", "L1")

    for (dist_code in 0:2) {
      dist_name <- dist_names[dist_code + 1L]

      for (k in c(1L, 4L, nrow(rpm_x), nrow(rpm_x) + 2L)) {
        kk <- min(k, nrow(rpm_x))

        reference <- .nns_mreg_predict_reference(
          Xtest, rpm_reg, kk, dist_name, mins, maxs, FALSE
        )
        legacy <- NNS_mreg_predict_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE
        )
        single_one <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
        )
        single_many <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
        )
        path_one <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
        )
        path_many <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
        )

        expect_equal(legacy, reference, tolerance = 1e-12)
        expect_equal(single_one, reference, tolerance = 1e-12)
        expect_equal(single_many, single_one, tolerance = 1e-12)
        expect_equal(path_one[, kk], reference, tolerance = 1e-12)
        expect_equal(path_many, path_one, tolerance = 1e-12)

        reference_cls <- .nns_mreg_predict_reference(
          Xtest, rpm_cls, kk, dist_name, mins, maxs, TRUE
        )
        legacy_cls <- NNS_mreg_predict_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE
        )
        single_cls_one <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
        )
        single_cls_many <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
        )
        path_cls_one <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
        )
        path_cls_many <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
        )

        expect_equal(legacy_cls, reference_cls, tolerance = 1e-12)
        expect_equal(single_cls_one, reference_cls, tolerance = 1e-12)
        expect_equal(single_cls_many, single_cls_one, tolerance = 1e-12)
        expect_equal(path_cls_one[, kk], reference_cls, tolerance = 1e-12)
        expect_equal(path_cls_many, path_cls_one, tolerance = 1e-12)
      }
    }
  })
  aud_case("partial-selection and single-k kernels match references", {
    rpm_x <- rbind(
      c(0, 0), c(0, 2), c(2, 0), c(2, 2),
      c(1, 1), c(1, 1), c(3, 1), c(-1, 1)
    )
    y_reg <- c(0, 2, 4, 6, 3, 5, 7, -1)
    y_cls <- c(1, 1, 2, 2, 3, 3, 4, 4)
    Xtest <- rbind(c(1, 0), c(1, 1), c(1, 2), c(0.5, 1.5))
    mins <- apply(rpm_x, 2, min)
    maxs <- apply(rpm_x, 2, max)

    rpm_reg <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_reg)
    rpm_cls <- data.frame(V1 = rpm_x[, 1], V2 = rpm_x[, 2], y.hat = y_cls)
    dist_names <- c("NNS", "L2", "L1")

    for (dist_code in 0:2) {
      dist_name <- dist_names[dist_code + 1L]

      for (k in c(1L, 2L, 4L, nrow(rpm_x), nrow(rpm_x) + 3L)) {
        kk <- min(k, nrow(rpm_x))

        ref_reg <- .nns_mreg_predict_reference(
          Xtest, rpm_reg, kk, dist_name, mins, maxs, FALSE
        )
        legacy_reg <- NNS_mreg_predict_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE
        )
        single_reg_1 <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
        )
        single_reg_2 <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
        )
        path_reg_1 <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 1L
        )
        path_reg_2 <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_reg, Xtest, k, dist_code, mins, maxs, FALSE, 2L
        )

        expect_equal(legacy_reg, ref_reg, tolerance = 1e-12)
        expect_equal(single_reg_1, ref_reg, tolerance = 1e-12)
        expect_equal(single_reg_2, single_reg_1, tolerance = 1e-12)
        expect_equal(path_reg_1[, kk], ref_reg, tolerance = 1e-12)
        expect_equal(path_reg_2, path_reg_1, tolerance = 1e-12)

        ref_cls <- .nns_mreg_predict_reference(
          Xtest, rpm_cls, kk, dist_name, mins, maxs, TRUE
        )
        legacy_cls <- NNS_mreg_predict_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE
        )
        single_cls_1 <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
        )
        single_cls_2 <- NNS_mreg_predict_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
        )
        path_cls_1 <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 1L
        )
        path_cls_2 <- NNS_mreg_predict_path_v2_cpp(
          rpm_x, y_cls, Xtest, k, dist_code, mins, maxs, TRUE, 2L
        )

        expect_equal(legacy_cls, ref_cls, tolerance = 1e-12)
        expect_equal(single_cls_1, ref_cls, tolerance = 1e-12)
        expect_equal(single_cls_2, single_cls_1, tolerance = 1e-12)
        expect_equal(path_cls_1[, kk], ref_cls, tolerance = 1e-12)
        expect_equal(path_cls_2, path_cls_1, tolerance = 1e-12)
      }
    }
  })
  aud_case("v2 all-k path matches legacy path on smoke data", {
    rpm <- cbind(c(0,1,2), c(0,1,4)); x <- cbind(c(.5,1.5)); mins <- 0; maxs <- 2
    old <- NNS_mreg_predict_path_cpp(rpm[,1,drop=FALSE], rpm[,2], x, 3L, 0L, mins, maxs, FALSE)
    nat <- NNS_mreg_predict_path_v2_cpp(rpm[,1,drop=FALSE], rpm[,2], x, 3L, 0L, mins, maxs, FALSE, 1L)
    expect_equal(nat, old, tolerance = 1e-12)
  })
  aud_case("v2 all-k path is deterministic with multiple native threads", {
    rpm_x <- matrix(seq(0, 1, length.out = 40), ncol = 2)
    y <- seq_len(nrow(rpm_x)) / 10
    xt <- rpm_x[c(1, 5, 10), , drop = FALSE]
    one <- NNS_mreg_predict_path_v2_cpp(rpm_x, y, xt, 6L, 0L,
                                        apply(rpm_x, 2, min), apply(rpm_x, 2, max),
                                        FALSE, 1L)
    many <- NNS_mreg_predict_path_v2_cpp(rpm_x, y, xt, 6L, 0L,
                                         apply(rpm_x, 2, min), apply(rpm_x, 2, max),
                                         FALSE, 2L)
    expect_equal(many, one, tolerance = 1e-12)
  })
  aud_case("native xstar path matches active coefficient projections including zero fallback", {
    train <- matrix(1:12, 4); test <- matrix(13:18, 2); coef <- c(2, 0, -1)
    ord <- order(abs(coef), decreasing = TRUE, na.last = NA, method = "radix")
    out <- NNS_xstar_path_cpp(train, test, coef, ord, 1L)
    expect_equal(dim(out$train), c(4L, 3L)); expect_equal(out$column_order, ord)
    z <- NNS_xstar_path_cpp(train, test, c(0,0,0), seq_len(3), 1L)
    expect_equal(z$train[,2], rowMeans(train[,1:2]))
  })
  aud_case("native xstar path is deterministic with multiple threads", {
    train <- matrix(seq_len(60), 20, 3); test <- matrix(seq_len(15), 5, 3)
    coef <- c(3, 0, -2); ord <- order(abs(coef), decreasing = TRUE, na.last = NA, method = "radix")
    one <- NNS_xstar_path_cpp(train, test, coef, ord, 1L)
    many <- NNS_xstar_path_cpp(train, test, coef, ord, 2L)
    expect_equal(many$train, one$train, tolerance = 1e-14)
    expect_equal(many$test, one$test, tolerance = 1e-14)
    expect_identical(many$representative, one$representative)
  })
  aud_case("univariate fast helper matches point estimates from public path", {
    set.seed(13); x <- rnorm(20); y <- sin(x); test <- x[1:5]
    fast <- .nns_reg_univariate_fast(x, y, test)
    pub <- NNS.reg(x, y, point.est = test, plot = FALSE, residual.plot = FALSE, point.only = TRUE)
    expect_equal(fast$prediction, as.numeric(pub$Point.est), tolerance = 1e-12)
  })
  aud_case("native kernels equal the reference when the k nearest regression points are equidistant", {
    # Four regression points on the corners of a square around the test points:
    # the k nearest are at one distance, so the normal and rbf weights hit their
    # sd == 0 and var == 0 branches. With range 4 every scaled coordinate and
    # distance is exact; on near-ties R and C++ can pick different neighbours.
    rpm <- data.frame(x1 = c(0, 2, 0, 2, 4), x2 = c(0, 0, 2, 2, 4), y.hat = c(1, 2, 3, 4, 10))
    X <- as.matrix(rpm[, 1:2]); Xtest <- rbind(c(1, 1), c(1, 1.5))
    for (dist in c("L2", "L1", "NNS")) for (k in 2:4) {
      code <- match(dist, c("NNS", "L2", "L1", "FACTOR")) - 1L
      expect_equal(NNS:::NNS_mreg_predict_v2_cpp(X, rpm$y.hat, Xtest, k, code, c(0, 0), c(4, 4), FALSE, 1L),
                   NNS:::.nns_mreg_predict_reference(Xtest, rpm, k, dist, c(0, 0), c(4, 4), FALSE), tolerance = 1e-12)
    }
  })
})

test_that("MREG03: native setup: RPM and duplicate map, one-based rightmost-closed interval IDs on every boundary", {
  withr::local_seed(20261108)
  aud_case("native mreg setup returns RPM and duplicate map smoke outputs", {
    X <- cbind(a = c(1,1,2,2), b = c(1,1,2,2)); y <- c(1,2,3,4)
    prep <- .nns_mreg_prepare_model(X, y, order = "max", use.native = TRUE)
    expect_true(is.data.frame(prep$RPM)); expect_identical(as.integer(NNS_duplicate_column_map_cpp(X)), c(1L,1L))
  })
  aud_case("native mreg setup preserves one-based rightmost-closed findInterval IDs", {
    X <- matrix(as.numeric(1:12), ncol = 1, dimnames = list(NULL, "x"))
    y <- as.numeric(1:12)
    boundaries <- list(as.numeric(1:12))
    native <- NNS_mreg_setup_cpp(X, y, boundaries, 3L, FALSE)
    reference_ids <- as.character(findInterval(X[, 1], boundaries[[1]],
                                               left.open = FALSE,
                                               rightmost.closed = TRUE))
    expect_identical(as.character(native$ids), reference_ids)
    expect_identical(tail(as.character(native$ids), 2), c("11", "11"))

    native_rpm <- as.data.frame(native$RPM)
    names(native_rpm) <- c("x", "y.hat")
    reference_rpm <- .nns_mreg_build_rpm(as.data.frame(X), y, reference_ids,
                                         "off", FALSE)
    expect_identical(native_rpm$x, reference_rpm$x)
    expect_identical(native_rpm$y.hat, reference_rpm$y.hat)

    reference_groups <- split(seq_len(nrow(X)), reference_ids)
    expect_identical(
      names(reference_groups),
      c("1", "10", "11", "2", "3", "4", "5", "6", "7", "8", "9")
    )
    merged_group_position <- match("11", names(reference_groups))
    expect_false(is.na(merged_group_position))
    expect_equal(native_rpm$x[merged_group_position], 11.5)
    expect_equal(native_rpm$y.hat[merged_group_position], 11.5)
    expect_equal(tail(native_rpm$x, 1), 9)
    expect_equal(tail(native_rpm$y.hat, 1), 9)
  })
  aud_case("native interval IDs match R for below first, exact boundaries, and above final", {
    boundaries <- list(as.numeric(c(1, 2, 3, 4)))
    X <- matrix(as.numeric(c(0, 1, 1.5, 2, 3.5, 4, 5)), ncol = 1)
    native <- NNS_mreg_setup_cpp(X, seq_len(nrow(X)), boundaries, 3L, FALSE)
    reference_ids <- as.character(findInterval(X[, 1], boundaries[[1]],
                                               left.open = FALSE,
                                               rightmost.closed = TRUE))
    expect_identical(as.character(native$ids), reference_ids)
    expect_identical(as.character(native$ids), as.character(c(0, 1, 1, 2, 3, 3, 4)))
  })
  aud_case("native interval IDs match R across multiple predictors and boundary counts", {
    X <- cbind(
      a = as.numeric(c(0, 1, 1.5, 2, 3.5, 4, 5, 4, 4)),
      b = as.numeric(c(-1, 0, 5, 10, 20, 20, 25, 20, 20)),
      c = rep(7, 9)
    )
    boundaries <- list(
      as.numeric(c(1, 2, 3, 4)),
      as.numeric(c(0, 10, 20)),
      as.numeric(7)
    )
    native <- NNS_mreg_setup_cpp(X, seq_len(nrow(X)), boundaries, 3L, FALSE)
    id_parts <- lapply(seq_along(boundaries), function(j) {
      findInterval(X[, j], boundaries[[j]], left.open = FALSE,
                   rightmost.closed = TRUE)
    })
    reference_ids <- do.call(paste, c(as.data.frame(id_parts), sep = "."))
    expect_identical(as.character(native$ids), reference_ids)

    native_rpm <- as.data.frame(native$RPM)
    names(native_rpm) <- c(colnames(X), "y.hat")
    reference_rpm <- .nns_mreg_build_rpm(as.data.frame(X), seq_len(nrow(X)),
                                         reference_ids, "off", FALSE)
    expect_identical(native_rpm$a, reference_rpm$a)
    expect_identical(native_rpm$b, reference_rpm$b)
    expect_identical(native_rpm$c, reference_rpm$c)
    expect_identical(native_rpm$y.hat, reference_rpm$y.hat)
  })
})

test_that("MREG04: distance modes: NULL is NNS, NNS differs from L1 and L2, R and C++ agree, explicit modes, code dispatch, bulk wrappers by position", {
  withr::local_seed(20261108)
  aud_case("NNS.reg NULL distance is native NNS alias", {
    x <- cbind(a = c(0, 1, 0.2, 0.9, 0.4, 0.8), b = c(0, 0.2, 1, 0.9, 0.7, 0.1))
    y <- c(0, 2, 4, 8, 5, 3)
    pts <- cbind(a = c(0.35, 0.75), b = c(0.85, 0.15))

    fit_null <- NNS.reg(x, y, point.est = pts, n.best = 3, dist = NULL,
                        plot = FALSE, residual.plot = FALSE, ncores = 1)
    fit_nns <- NNS.reg(x, y, point.est = pts, n.best = 3, dist = "nns",
                       plot = FALSE, residual.plot = FALSE, ncores = 1)

    expect_identical(fit_null$dist, "NNS")
    expect_equal(fit_null$Point.est, fit_nns$Point.est, tolerance = 1e-12)
  })
  aud_case("native NNS distance differs from L1 and L2 when rankings differ", {
    x <- cbind(a = c(0, 0.2, 0.9, 1), b = c(0, 1, 0.2, 1))
    y <- c(0, 10, 20, 30)
    pts <- cbind(a = 0.55, b = 0.55)

    nns <- NNS.reg(x, y, point.est = pts, n.best = 2, dist = NULL,
                   plot = FALSE, residual.plot = FALSE, ncores = 1)$Point.est
    l1 <- NNS.reg(x, y, point.est = pts, n.best = 2, dist = "L1",
                  plot = FALSE, residual.plot = FALSE, ncores = 1)$Point.est
    l2 <- NNS.reg(x, y, point.est = pts, n.best = 2, dist = "L2",
                  plot = FALSE, residual.plot = FALSE, ncores = 1)$Point.est

    expect_false(isTRUE(all.equal(nns, l1)))
    expect_false(isTRUE(all.equal(nns, l2)))
  })
  aud_case("native R and C++ multivariate prediction paths agree", {
    x <- cbind(a = c(0, 1, 0.2, 0.9, 0.4, 0.8), b = c(0, 0.2, 1, 0.9, 0.7, 0.1))
    y <- c(0, 2, 4, 8, 5, 3)
    pts <- cbind(a = c(0.35, 0.75), b = c(0.85, 0.15))

    fit <- NNS.reg(x, y, point.est = pts, n.best = 3, dist = NULL,
                   plot = FALSE, residual.plot = FALSE, ncores = 1)
    ref <- .nns_mreg_predict_reference(pts, fit$RPM, 3, "NNS",
                                       apply(x, 2, min), apply(x, 2, max), FALSE)
    expect_equal(fit$Point.est, ref, tolerance = 1e-12)
  })
  aud_case("explicit distance modes and invalid values are handled", {
    withr::local_seed(20261108)
    x <- cbind(a = rnorm(12), b = rep(1:3, 4))
    y <- seq_len(12)
    for (d in c("NNS", "L1", "L2", "FACTOR")) {
      expect_no_error(NNS.reg(x, y, point.est = x[1:2, ], n.best = 2, dist = d,
                              plot = FALSE, residual.plot = FALSE, ncores = 1))
    }
  })
  aud_case("native distance codes dispatch to their intended calculations", {
    rpm_x <- matrix(c(0, 0,
                      0.6, 0.6,
                      1, 0,
                      0.2, 1), ncol = 2, byrow = TRUE)
    yhat <- c(1, 2, 3, 4)
    # Keep the dispatch fixture away from mathematically tied distances.
    # Symmetric coordinates such as c(0.45, 0.15) give the first two RPM rows
    # identical theoretical NNS/L1/L2 distances, whose last-bit ordering can
    # differ across compilers (notably Apple ARM64) and change rank weights.
    xtest <- matrix(c(0.47, 0.11,
                      0.2, 0.8), ncol = 2, byrow = TRUE)
    mins <- c(0, 0)
    maxs <- c(1, 1)
    rpm <- data.frame(a = rpm_x[, 1], b = rpm_x[, 2], y.hat = yhat)
    code_map <- c(NNS = 0L, L2 = 1L, L1 = 2L, FACTOR = 3L)

    for (d in names(code_map)) {
      ref <- .nns_mreg_predict_reference(xtest, rpm, 2L, d, mins, maxs, FALSE)
      got <- NNS_mreg_predict_cpp(rpm_x, yhat, xtest, 2L, code_map[[d]], mins, maxs, FALSE)
      expect_equal(got, ref, tolerance = 1e-12, info = d)
    }

    expect_error(
      NNS_mreg_predict_cpp(rpm_x, yhat, xtest, 2L, 99L, mins, maxs, FALSE),
      "distance code"
    )
  })
  aud_case("bulk distance wrappers match an unnamed matrix by position (#150)", {
    rpm <- data.frame(X1 = c(0, 1, 0, 1), X2 = c(0, 0, 1, 1), y.hat = c(1, 2, 3, 4))
    q <- matrix(c(.2, .9, .1, .8), 2)
    expect_equal(NNS:::NNS.distance.bulk(rpm, q, 2),
                 NNS:::NNS.distance.bulk(rpm, setNames(as.data.frame(q), c("X1", "X2")), 2))
    expect_equal(NNS:::NNS.distance.path.bulk(rpm, q, 3),
                 NNS:::NNS.distance.path.bulk(rpm, setNames(as.data.frame(q), c("X1", "X2")), 3))
  })
})
