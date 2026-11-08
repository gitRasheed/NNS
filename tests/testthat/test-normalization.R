# NNS.norm.

test_that("NORM01: NNS.norm: linear and nonlinear matrix paths and equal- and unequal-length lists", {
  withr::local_seed(20261108)
  normalization_ref <- fx_normalization_ref()
  aud_case("NNS.norm - linear matrix path equalizes means", env = normalization_ref, {
    N <- NNS.norm(A, linear = TRUE)
    expect_equal(dim(N), dim(A))
    expect_equal(as.numeric(colMeans(N)), rep(mean(colMeans(A)), 3), tolerance = 1e-12)
  })
  aud_case("NNS.norm - nonlinear matrix path returns finite matrix", env = normalization_ref, {
    N <- NNS.norm(A)
    expect_equal(dim(N), dim(A))
    expect_true(all(is.finite(N)))
  })
  aud_case("NNS.norm - equal-length list matches matrix path with linear = FALSE", env = normalization_ref, {
    N_matrix <- NNS.norm(A)
    N_list <- NNS.norm(list(x, y, z))
    expect_equal(unname(N_list), unname(N_matrix), tolerance = 1e-12)
  })
  aud_case("NNS.norm - unequal-length list forces linear scaling", {
    vec1 <- c(1, 2, 3, 4, 5, 6, 7)
    vec2 <- c(10, 20, 30, 40, 50, 60)
    vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
    N <- NNS.norm(list(vec1, vec2, vec3))
    expect_true(is.list(N))
    expect_equal(lengths(N, use.names = FALSE), c(7L, 6L, 9L))
    m <- c(mean(vec1), mean(vec2), mean(vec3))
    expect_equal(unname(sapply(N, mean)), rep(mean(m), 3), tolerance = 1e-12)
  })
  aud_case("NNS.norm - zero-sum length diffs still detected as unequal", {
    vec1 <- c(1, 2, 3, 4, 5)
    vec2 <- c(10, 20, 30, 40, 50, 60)
    vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    N <- NNS.norm(list(vec1, vec2, vec3))
    expect_true(is.list(N))
    expect_equal(lengths(N, use.names = FALSE), c(5L, 6L, 5L))
    m <- c(mean(vec1), mean(vec2), mean(vec3))
    expect_equal(unname(sapply(N, mean)), rep(mean(m), 3), tolerance = 1e-12)
  })
})
