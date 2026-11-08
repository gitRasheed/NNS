# Efficient sets and dominance clusters.

# pin: NNS-Python reference values
test_that("S01: NNS-Python reference efficient sets at orders 1 to 3, column order and caller matrix untouched", {
  withr::local_seed(20261108)
  x <- nns_py$x
  y <- nns_py$y
  z <- nns_py$z
  xx <- x+10
  yy <- y+10
  zz <- z+10
  Z <- matrix(c(x,y,z,xx,yy,zz),ncol=6)
  colnames(Z) <- c("x", "y", "z", "xx", "yy", "zz")
  # type is left at its default for degrees 2 and 3.
  sets <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, text = "
    degree | type       | expected
    1      | discrete   | yy zz xx
    1      | continuous | yy zz xx
    2      |            | yy xx
    3      |            | yy xx")
  sets$degree <- as.numeric(sets$degree)
  sets$label <- sprintf("NNS.SD.efficient.set(Z, degree = %g%s)", sets$degree,
                        ifelse(sets$type == "", "", paste0(", type = ", sets$type)))
  aud_cases(sets, function(r) {
    got <- if (r$type == "") NNS.SD.efficient.set(x=Z, degree=r$degree, status=F) else
      NNS.SD.efficient.set(x=Z, degree=r$degree, type=r$type, status=F)
    expect_equal(got, strsplit(r$expected, " ")[[1]])
  })
  aud_case("NNS.SD.efficient.set leaves the caller's matrix unnamed", {
    set.seed(1)
    m <- matrix(rnorm(300), 100, 3)
    m2 <- m
    named <- m + 0
    colnames(named) <- c("X_1", "X_2", "X_3")

    out <- NNS.SD.efficient.set(m, 1, status = FALSE)

    expect_null(colnames(m))
    expect_null(colnames(m2))
    expect_identical(out, NNS.SD.efficient.set(named, 1, status = FALSE))
  })
  aud_case("the efficient set does not depend on column order", {
    x <- c(-1, -1, 2); y <- c(-2, 1, 1)
    expect_equal(NNS.FSD.uni(x, y, "continuous"), 1)
    expect_identical(NNS.SD.efficient.set(cbind(Y = y, X = x), 1, "continuous", status = FALSE), "X")
    expect_identical(NNS.SD.efficient.set(cbind(X = x, Y = y), 1, "continuous", status = FALSE), "X")
    expect_identical(NNS.SD.cluster(cbind(Y = y, X = x), 1, "continuous")$Clusters,
                     NNS.SD.cluster(cbind(X = x, Y = y), 1, "continuous")$Clusters)

    a <- 0:5; h <- 3 * 2^-44
    M <- cbind(A = a, B = a + h * c(0, -2, 3, 0, -1, 0), C = a + h * c(0, -1, 0, 3, -2, 0))
    D <- NNS:::sd_dom_matrix_prefix_parallel(M, 2L, "discrete")
    expect_identical(D[cbind(c(1, 2, 3), c(2, 3, 1))], rep(1L, 3))
    for (p in list(1:3, c(2, 3, 1), c(3, 2, 1), c(2, 1, 3)))
      expect_setequal(NNS.SD.efficient.set(M[, p], 2, status = FALSE), c("A", "B", "C"))

    differ <- character(0)
    for (s in 1:200) {
      set.seed(s)
      nc <- sample(4:6, 1)
      M <- matrix(sample(-3:3, 5 * nc, replace = TRUE), 5, nc, dimnames = list(NULL, paste0("C", 1:nc)))
      p <- sample(nc)
      for (d in 1:3) for (type in if (d == 1) c("discrete", "continuous") else "discrete") {
        if (!setequal(NNS.SD.efficient.set(M[, p], d, type, status = FALSE),
                      NNS.SD.efficient.set(M, d, type, status = FALSE)))
          differ <- c(differ, paste(s, d, type))
      }
    }
    expect_identical(differ, character(0))
  })
})

test_that("S02: NNS.SD.cluster keeps the dendrogram leaf order and the final peeled cluster", {
  withr::local_seed(20261108)
  aud_case("the returned dendrogram keeps hclust's leaf order", {
    set.seed(123)
    A <- matrix(rnorm(200), 10, 20)
    colnames(A) <- paste0("V", 1:20)
    pdf(NULL); on.exit(dev.off())
    res <- NNS.SD.cluster(A, degree = 1, min_cluster = 2, dendrogram = TRUE)

    cl <- res$Clusters
    expect_gt(length(cl), 1)
    v <- unlist(cl, use.names = FALSE)
    lab <- rep(seq_along(cl), lengths(cl))
    k <- seq_along(v)
    d <- as.dist(outer(lab, lab, function(a, b) 20 * abs(a - b)) + 1e-3 * outer(k, k, function(i, j) abs(i - j)))
    attr(d, "Labels") <- v
    hc <- hclust(d, method = "complete")

    expect_identical(res$Dendrogram$order, hc$order)
    expect_identical(order.dendrogram(as.dendrogram(res$Dendrogram)), res$Dendrogram$order)
    expect_identical(unname(labels(as.dendrogram(res$Dendrogram))), v[hc$order])
  })
  aud_case("variables left when peeling stops at once form the final cluster", {
    A <- cbind(A = 1:5, B = 2:6)
    expect_identical(NNS.SD.cluster(A, min_cluster = 2)$Clusters, list(Cluster_1 = c("A", "B")))
    expect_identical(NNS.SD.cluster(A, min_cluster = 3)$Clusters, list(Cluster_1 = c("A", "B")))
    expect_identical(NNS.SD.cluster(A[, "A", drop = FALSE])$Clusters, list(Cluster_1 = "A"))

    pdf(NULL); on.exit(dev.off())
    res <- NNS.SD.cluster(A, min_cluster = 2, dendrogram = TRUE)
    expect_identical(res$Clusters, list(Cluster_1 = c("A", "B")))
    expect_identical(unname(res$Dendrogram$labels), c("A", "B"))
    res <- NNS.SD.cluster(cbind(A = 1:5, B = 11:15), dendrogram = TRUE)
    expect_identical(res$Clusters, list(Cluster_1 = "B", Cluster_2 = "A"))
    expect_warning(res <- NNS.SD.cluster(A[, "A", drop = FALSE], dendrogram = TRUE), "Not enough variables")
    expect_identical(res$Clusters, list(Cluster_1 = "A"))
  })
})
