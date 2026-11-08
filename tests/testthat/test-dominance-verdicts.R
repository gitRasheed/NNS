# Stochastic dominance verdicts. Lean: nns-proofs Integration.lean, DominanceAlgo.lean.

# Lean: LPM_one_le_of_LPM_zero_le, LPM_two_le_of_LPM_zero_le, LPM_two_le_of_LPM_one_le, LPM_add_le_of_forall_le, mean_le_of_LPM_two_le, exists_le_of_forall_LPM_le, mean_le_of_LPM_one_le, sdDomPair_omega_iff, UPM_succ_le_of_forall_le | BUGS: none
test_that("D01: FSD implies SSD implies TSD, and dominance orders the minima and means", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("LPM_one_le_of_LPM_zero_le: FSD.uni == 1 implies SSD.uni == 1 (constructed FSD pairs)", {
    expect_gte(sum(fm$fsd_verdicts[, "fsd"]), 90)
    expect_identical(fm$fsd_verdicts[fm$fsd_verdicts[, "fsd"] == 1L, "ssd"],
                     rep(1L, sum(fm$fsd_verdicts[, "fsd"])))
  })
  aud_case("LPM_two_le_of_LPM_zero_le: FSD.uni == 1 implies TSD.uni == 1 (constructed FSD pairs)", {
    expect_identical(fm$fsd_verdicts[fm$fsd_verdicts[, "fsd"] == 1L, "tsd"],
                     rep(1L, sum(fm$fsd_verdicts[, "fsd"])))
  })
  aud_case("LPM_one_le_of_LPM_zero_le / LPM_two_le_of_LPM_one_le: FSD => SSD => TSD on random integer pairs", {
    v <- fm$int_verdicts
    expect_gt(sum(v[, "fsd"]), 20)
    expect_true(all(v[v[, "fsd"] == 1L, "ssd"] == 1L))
    expect_true(all(v[v[, "ssd"] == 1L, "tsd"] == 1L))
  })
  aud_case("LPM_add_le_of_forall_le: FSD => SSD => TSD on real data (mtcars, faithful)", {
    pairs <- list(
      mpg = list(x = mtcars$mpg[mtcars$cyl == 4], y = mtcars$mpg[mtcars$cyl == 8]),
      eruptions = list(x = faithful$eruptions[faithful$waiting > 70],
                       y = faithful$eruptions[faithful$waiting <= 70]),
      dax_shift = list(x = fm$real_samples$dax, y = fm$real_samples$dax - 0.001)
    )
    for (nm in names(pairs)) {
      expect_identical(aud_verdicts(pairs[[nm]]), c(fsd = 1L, ssd = 1L, tsd = 1L), label = nm)
    }
  })
  aud_case("exists_le_of_forall_LPM_le / mean_le_of_LPM_one_le (SSD) / mean_le_of_LPM_two_le (TSD): dominance forces min and mean gates", {
    for (i in seq_along(fm$int_pairs)) {
      p <- fm$int_pairs[[i]]
      if (any(fm$int_verdicts[i, ] == 1L) || fm$int_omega[i] == 1L) expect_gte(min(p$x), min(p$y))
      if (fm$int_verdicts[i, "ssd"] == 1L) {                       # mean_le_of_LPM_one_le
        expect_gte(sum(p$x) * length(p$y), sum(p$y) * length(p$x))
      }
      if (fm$int_verdicts[i, "tsd"] == 1L) {                       # mean_le_of_LPM_two_le
        expect_gte(sum(p$x) * length(p$y), sum(p$y) * length(p$x))
      }
    }
  })
  aud_case("LPM_one_le_of_LPM_zero_le + UPM_succ_le_of_forall_le (d = 0) => sdDomPair_omega_iff: discrete FSD == 1 implies Omega == 1", {
    # F_x <= F_y gives L1_x <= L1_y (LPM_one_le_of_LPM_zero_le) and U1_x >= U1_y
    # (UPM_succ_le_of_forall_le at d = 0, since UPM0 = 1 - F), hence R_x <= R_y
    expect_gt(sum(fm$int_verdicts[, "fsd"]), 20)
    expect_true(all(fm$int_omega[fm$int_verdicts[, "fsd"] == 1L] == 1L))
    fsd_omega <- vapply(fm$fsd_pairs, function(p) NNS.FSD.uni(p$x, p$y, type = "continuous"), integer(1))
    expect_true(all(fsd_omega[fm$fsd_verdicts[, "fsd"] == 1L] == 1L))
  })
})

# Lean: LPM_one_eq_of_eq_at_knots, sdDomPair_ssd_iff, LPM_one_le_of_le_at_knots, sdDomPair_tsd_iff, LPM_two_le_of_le_at_knots, tsdKnotOnly_vertexCex, sdDomPair_tsd_vertexCex, tsdVertexCex_le_at_knots, tsdVertexCex_mean_le, tsdVertexCex_not_le_at_three, tsdCex_le_on_Iic, tsdCex_mean_lt, tsdCex_not_le_at_three | BUGS: none
test_that("D02: SSD and TSD verdicts equal the exact knot, vertex and mean-gate criteria", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("sdDomPair_ssd_iff: SSD.uni == 1 iff LPM1_x <= LPM1_y at all merged knots, strict at one (base-R oracle)", {
    # LPM_one_le_of_le_at_knots gives the knot sufficiency, LPM_one_eq_of_eq_at_knots
    # makes knot-only strictness exact; the oracle is mean(pmax(t - x, 0)).
    ssd_knot <- vapply(fm$int_pairs, function(p) {
      k <- sort(unique(c(p$x, p$y)))
      d <- aud_lpm1_ref(k, p$x) - aud_lpm1_ref(k, p$y)
      as.integer(all(d <= 0) && any(d < 0))
    }, integer(1))
    expect_gt(sum(ssd_knot), 50)
    expect_identical(fm$int_verdicts[, "ssd"], ssd_knot)
  })
  aud_case("sdDomPair_ssd_iff: SSD knot criterion on random rounded and real data (base-R oracle)", {
    withr::local_seed(3)
    pairs <- list(
      list(x = fm$random_samples$normal[1:200] + 0.3, y = fm$random_samples$normal[201:400]),
      list(x = fm$random_samples$t3[1:250], y = fm$random_samples$t3[251:500]),
      list(x = round(rexp(150), 2), y = round(rexp(120), 2) - 0.2),
      list(x = fm$real_samples$dax, y = as.numeric(diff(log(EuStockMarkets[, "FTSE"])))),
      list(x = as.numeric(diff(log(EuStockMarkets[, "SMI"]))), y = as.numeric(diff(log(EuStockMarkets[, "CAC"])))),
      list(x = fm$real_samples$mag + 0.05, y = fm$real_samples$mag),
      list(x = mtcars$mpg[mtcars$am == 1], y = mtcars$mpg[mtcars$am == 0])
    )
    for (i in seq_along(pairs)) {
      p <- pairs[[i]]
      k <- sort(unique(c(p$x, p$y)))
      d <- aud_lpm1_ref(k, p$x) - aud_lpm1_ref(k, p$y)
      S <- max(abs(c(p$x, p$y)))
      want <- aud_ssd_expected(d, 64 * aud_eps * S + 4 * (length(p$x) + length(p$y)) * aud_eps * S)
      expect_false(is.na(want), label = paste("pair", i, "outside the tolerance band"))
      expect_identical(NNS.SSD.uni(p$x, p$y), want, label = paste("pair", i))
    }
  })
  aud_case("sdDomPair_tsd_iff: TSD.uni equals the exact knot + vertex + mean-gate criterion on random integer pairs", {
    tsd_ref <- vapply(fm$int_pairs, function(p) aud_tsd_exact(p$x, p$y), integer(1))
    expect_gt(sum(tsd_ref), 50)
    # on this regime only the knots and the mean gate decide
    mean_decided <- vapply(fm$int_pairs, function(p) {
      q <- aud_tsd_parts(p$x, p$y); q$knots_ok && q$differs && !q$mean_ok
    }, logical(1))
    expect_gt(sum(mean_decided), 0)
    expect_identical(fm$int_verdicts[, "tsd"], tsd_ref)
  })
  aud_case("sdDomPair_tsd_iff (vertex branch): TSD.uni agrees with the exact criterion where the interior vertex fires", {
    # LPM_two_le_of_le_at_knots needs the vertex hypothesis (src/SD.cpp:264-277);
    # here the knots, min and mean gates all pass and the vertex alone says 0.
    expect_identical(fm$n_vertex_decided, 25L)
    decided <- vapply(fm$vertex_pairs, `[[`, logical(1), "decided")
    expect_gt(sum(!decided), 0)
    for (i in seq_along(fm$vertex_pairs)) {
      p <- fm$vertex_pairs[[i]]
      expect_identical(aud_tsd_exact(p$x, p$y), 0L, label = paste("vertex pair", i))
      expect_identical(NNS.TSD.uni(p$x, p$y), 0L, label = paste("vertex pair", i))
    }
    # the same pairs reversed exercise the vertex walk on the other orientation
    rev_ref <- vapply(fm$vertex_pairs, function(p) aud_tsd_exact(p$y, p$x), integer(1))
    rev_got <- vapply(fm$vertex_pairs, function(p) NNS.TSD.uni(p$y, p$x), integer(1))
    expect_identical(rev_got, rev_ref)
  })
  aud_case("tsdKnotOnly_vertexCex / sdDomPair_tsd_vertexCex: knot-only TSD check is insufficient; x=(0,4,5,1,4), y=(4,2,4,0) is no-TSD", {
    x <- c(0, 4, 5, 1, 4); y <- c(4, 2, 4, 0)
    k <- sort(unique(c(x, y)))
    expect_true(all(LPM(2, k, x) - LPM(2, k, y) <= 0))     # tsdVertexCex_le_at_knots
    expect_gte(mean(x), mean(y))                            # tsdVertexCex_mean_le
    expect_equal(LPM(2, 3, x) - LPM(2, 3, y), 0.1, tolerance = 8 * aud_eps)  # tsdVertexCex_not_le_at_three
    q <- aud_tsd_parts(x, y)
    expect_true(q$knots_ok && q$mean_ok && q$differs)       # tsdKnotOnly_vertexCex: knot-only says 1
    expect_true(q$vertex_bad)                               # the vertex branch decides
    expect_identical(aud_tsd_exact(x, y), 0L)
    expect_identical(NNS.TSD.uni(x, y), 0L)                 # sdDomPair_tsd_vertexCex
    expect_identical(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
  })
  aud_case("tsdCex (Integration.lean): the TSD mean gate is load-bearing; x = 0, y = (-1, 2) is no-TSD", {
    x <- 0; y <- c(-1, 2)
    s <- seq(-2, 2, by = 0.25)
    expect_true(all(LPM(2, s, x) <= LPM(2, s, y)))  # tsdCex_le_on_Iic
    expect_lt(mean(x), mean(y))                      # tsdCex_mean_lt
    expect_lt(LPM(2, 3, y), LPM(2, 3, x))            # tsdCex_not_le_at_three
    expect_identical(aud_tsd_exact(x, y), 0L)            # the mean gate decides
    expect_identical(NNS.TSD.uni(x, y), 0L)
    expect_identical(NNS.TSD(x, y, plot = FALSE), "NO TSD EXISTS")
  })
  aud_case("the TSD mean gate rejects a mean gap that the centred-pair veto lets through (kills M47)", {
    # pin: SSD/TSD mean gate (SD.cpp:239), kills M47. y is a spread of x whose mean is
    # higher by 100 eps * 1.2 = 2.7e-14: below the 128 eps max|x| = 3.4e-14 veto of
    # sd_dom_pair_centred (SD.cpp:145), above the inner gate, so only the gate says 0.
    expect_identical(NNS.TSD.uni(c(-1, -0.5, 0, 0.5, 1), c(-1.2, -0.5, 0, 0.5, 1.2) + 100 * .Machine$double.eps * 1.2), 0L)
  })
})

# Lean: sdDomPair_omega_iff | BUGS: none
test_that("D03: continuous FSD (Omega) equals the exact LPM.ratio(1) criterion", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("sdDomPair_omega_iff: FSD.uni(continuous) equals the exact LPM.ratio(1) criterion on random integer pairs", {
    parts <- lapply(fm$int_pairs, function(p) aud_omega_parts(p$x, p$y))
    ref <- vapply(fm$int_pairs, function(p) aud_omega_exact(p$x, p$y), integer(1))
    expect_gt(sum(ref), 50)
    # both interior branches of the kernel are exercised: the stationary point
    # decides (knots pass, interior crossing) and the midpoint alone shows the
    # curves differ (equal at every knot)
    expect_gt(sum(vapply(parts, function(q) q$knots_ok && q$vertex_bad, logical(1))), 0)
    expect_gt(sum(vapply(parts, function(q) q$knots_ok && !q$vertex_bad && !q$knot_differs &&
                           q$mid_differs, logical(1))), 0)
    expect_identical(fm$int_omega, ref)
    expect_identical(fm$int_omega_rev, vapply(fm$int_pairs, function(p) aud_omega_exact(p$y, p$x), integer(1)))
  })
  aud_case("sdDomPair_omega_iff: book pair x = (0,1,4,5), y = (0,2,2.5,3) is Omega but not FSD", {
    x <- c(0, 1, 4, 5); y <- c(0, 2, 2.5, 3)
    expect_identical(aud_omega_exact(x, y), 1L)
    expect_identical(NNS.FSD.uni(x, y, type = "continuous"), 1L)
    expect_identical(NNS.FSD.uni(x, y, type = "discrete"), 0L)
    expect_identical(NNS.FSD(x, y, type = "continuous", plot = FALSE), "X OMEGA Y")
    expect_identical(NNS.FSD(x, y, type = "discrete", plot = FALSE), "NO FSD EXISTS")
  })
  aud_case("sdDomPair_omega_iff: crossing ratio curves c(1,1) vs c(0,2) give 0 in both directions", {
    expect_identical(aud_omega_exact(c(1, 1), c(0, 2)), 0L)
    expect_identical(aud_omega_exact(c(0, 2), c(1, 1)), 0L)
    expect_identical(NNS.FSD.uni(c(1, 1), c(0, 2), type = "continuous"), 0L)
    expect_identical(NNS.FSD.uni(c(0, 2), c(1, 1), type = "continuous"), 0L)
    expect_identical(NNS.FSD(c(1, 1), c(0, 2), type = "continuous", plot = FALSE),
                     "NO OMEGA DOMINANCE EXISTS")
  })
  aud_case("sdDomPair_omega_iff (0/0 -> k/m fallback): constant samples", {
    # all observations at t: L1 + U1 = 0 and the ratio is k/m (LPM.ratio gives 1)
    expect_identical(LPM.ratio(1, 2, c(2, 2)), 1)
    expect_identical(NNS.FSD.uni(c(2, 2), 2, type = "continuous"), 0L)  # same distribution
    expect_identical(NNS.FSD.uni(2, c(2, 2), type = "continuous"), 0L)
    expect_identical(NNS.FSD.uni(c(3, 3), c(2, 2), type = "continuous"), 1L)
    expect_identical(NNS.FSD.uni(c(2, 2), c(3, 3), type = "continuous"), 0L)
    expect_identical(aud_omega_exact(c(3, 3), c(2, 2)), 1L)
  })
})

# Lean: sdDomPair_eq_zero_or_one, LPM_eq_of_identical, LPMDominates_shift, sdDomPairCentred_eq, sdDomPairCentred_omega_eq | BUGS: none
test_that("D04: verdicts are 0/1, antisymmetric and irreflexive; equal, one-point, constant and shifted samples", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("sdDomPair_eq_zero_or_one: every .uni verdict is the integer 0 or 1", {
    for (v in list(fm$int_verdicts, fm$int_verdicts_rev, fm$fsd_verdicts, fm$int_omega, fm$int_omega_rev)) {
      expect_type(v, "integer")
      expect_true(all(v %in% c(0L, 1L)))
    }
  })
  aud_case("sdDomPair_{fsd,ssd,tsd,omega}_iff: antisymmetry, never both X dominates Y and Y dominates X", {
    # exact dominance is antisymmetric (curves <= both ways means equal curves)
    expect_false(any(fm$int_verdicts == 1L & fm$int_verdicts_rev == 1L))
    expect_false(any(fm$int_omega == 1L & fm$int_omega_rev == 1L))
  })
  aud_case("LPM_eq_of_identical / sdDomPair_*_iff: a sample never dominates itself or a permutation of itself", {
    for (x in c(list(fm$x_int, fm$x_dy, 3), fm$random_samples[1], fm$real_samples["mpg"])) {
      expect_identical(aud_verdicts(list(x = x, y = rev(x))), c(fsd = 0L, ssd = 0L, tsd = 0L))
      expect_identical(NNS.FSD.uni(x, rev(x), type = "continuous"), 0L)
    }
  })
  aud_case("sdDomPair_*_iff (curves never differ): equal distributions of different sizes give 0", {
    x <- c(1, 2); y <- c(1, 1, 2, 2)
    expect_identical(aud_verdicts(list(x = x, y = y)), c(fsd = 0L, ssd = 0L, tsd = 0L))
    expect_identical(aud_verdicts(list(x = y, y = x)), c(fsd = 0L, ssd = 0L, tsd = 0L))
    expect_identical(NNS.FSD.uni(x, y, type = "continuous"), 0L)
    expect_identical(NNS.FSD.uni(y, x, type = "continuous"), 0L)
  })
  aud_case("sdDomPair_*_iff: n = 1 samples, 2 dominates 1 in every order, 1 vs 1 in none", {
    expect_identical(aud_verdicts(list(x = 2, y = 1)), c(fsd = 1L, ssd = 1L, tsd = 1L))
    expect_identical(aud_verdicts(list(x = 1, y = 2)), c(fsd = 0L, ssd = 0L, tsd = 0L))
    expect_identical(aud_verdicts(list(x = 1, y = 1)), c(fsd = 0L, ssd = 0L, tsd = 0L))
  })
  aud_case("sdDomPair_*_iff: constant samples, equal constants give 0, a higher constant dominates", {
    expect_identical(aud_verdicts(list(x = c(2, 2), y = 2)), c(fsd = 0L, ssd = 0L, tsd = 0L))
    expect_identical(aud_verdicts(list(x = rep(3, 4), y = c(1, 3))), c(fsd = 1L, ssd = 1L, tsd = 1L))
    # mean-preserving spread: the constant SSD/TSD-dominates, FSD fails
    expect_identical(aud_verdicts(list(x = c(1, 1), y = c(0, 2))), c(fsd = 0L, ssd = 1L, tsd = 1L))
  })
  aud_case("LPMDominates_shift / sdDomPairCentred_eq / sdDomPairCentred_omega_eq: shifting both samples by 1024 keeps every verdict", {
    shifted <- t(vapply(fm$int_pairs, function(p) aud_verdicts(list(x = p$x + 1024, y = p$y + 1024)),
                        integer(3)))
    expect_identical(shifted, fm$int_verdicts)
    shifted_omega <- vapply(fm$int_pairs, function(p) NNS.FSD.uni(p$x + 1024, p$y + 1024, type = "continuous"),
                            integer(1))
    expect_identical(shifted_omega, fm$int_omega)
  })
})

# Lean: none | BUGS: #27, #28
test_that("D05: bidirectional wrappers agree with the .uni verdicts whatever plot is (BUGS #27, #28)", {
  withr::local_seed(20261108)
  fm <- aud_fx_md()
  aud_case("plot = FALSE does not compute the plot curves", {
    local_mocked_bindings(LPM = function(...) stop("plot curve computed"),
                          LPM.ratio = function(...) stop("plot curve computed"))
    expect_identical(NNS.FSD(aud_x_om, aud_y_om, plot = FALSE), "NO FSD EXISTS")
    expect_identical(NNS.FSD(aud_x_om, aud_y_om, type = "continuous", plot = FALSE), "X OMEGA Y")
    expect_identical(NNS.SSD(c(1, 1), c(0, 2), plot = FALSE), "X SSD Y")
    expect_identical(NNS.TSD(c(0, 2), c(1, 1), plot = FALSE), "Y TSD X")
  })
  aud_case("the verdict does not depend on plot", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    withr::local_seed(28)
    x <- round(rnorm(30), 2); y <- round(rnorm(25, 0.3), 2)
    expect_identical(NNS.FSD(x, y, type = "continuous", plot = TRUE), NNS.FSD(x, y, type = "continuous", plot = FALSE))
    expect_identical(NNS.SSD(x, y, plot = TRUE), NNS.SSD(x, y, plot = FALSE))
    expect_identical(NNS.TSD(x, y, plot = TRUE), NNS.TSD(x, y, plot = FALSE))
  })
  aud_case("NNS.TSD plots against the thresholds, as NNS.SSD does", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    x <- c(0, 0.01, 0.02, 10); y <- c(1, 2, 3, 4)
    NNS.TSD(x, y, plot = TRUE)
    usr <- graphics::par("usr")[1:2]
    # the x-axis spans the data range [0, 10], not the index range [1, 8]
    expect_lte(usr[1], 0)
    expect_gte(usr[2], 10)
  })
  aud_case("bidirectional wrappers agree with the .uni verdicts and swap labels with the arguments (SPEC bidirectional wrappers)", {
    for (i in seq_len(60)) {
      p <- fm$int_pairs[[i]]
      for (ord in c("FSD", "SSD", "TSD")) {
        fn <- switch(ord, FSD = NNS.FSD, SSD = NNS.SSD, TSD = NNS.TSD)
        col <- tolower(ord)
        want <- if (fm$int_verdicts[i, col] == 1L) paste("X", ord, "Y")
                else if (fm$int_verdicts_rev[i, col] == 1L) paste("Y", ord, "X")
                else paste("NO", ord, "EXISTS")
        expect_identical(fn(p$x, p$y, plot = FALSE), want)
      }
      want <- if (fm$int_omega[i] == 1L) "X OMEGA Y"
              else if (fm$int_omega_rev[i] == 1L) "Y OMEGA X"
              else "NO OMEGA DOMINANCE EXISTS"
      expect_identical(NNS.FSD(p$x, p$y, type = "continuous", plot = FALSE), want)
    }
  })
})
