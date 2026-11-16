test_that("find_poisson_cutoff pools counts below N", {
  res <- find_poisson_cutoff(c(2, 4, 60), c(1e-6, 1e-6, 1e-6), N = 50)
  expect_equal(res$r, 6 / 2e-6)
  expect_equal(res$stderr, sqrt((6 / 2e-6)^2 / 6))
  expect_equal(res$type, "poisson_cutoff")
})

test_that("find_poisson_cutoff handles degenerate input", {
  expect_equal(
    find_poisson_cutoff(c(60, 70), c(1e-5, 1e-5), N = 50)$r,
    NA_real_
  )
  zero <- find_poisson_cutoff(c(0, 0), c(1e-5, 1e-5), N = 50)
  expect_equal(zero$r, 0)
  expect_equal(zero$stderr, NA_real_)
})

test_that("find_mle recovers the true rate from expected counts", {
  N <- 100
  r_true <- 3e7
  dil <- c(1e-6, 1e-5, 1e-4)
  counts <- N * (1 - exp(-r_true * dil / N))
  expect_equal(find_mle(counts, dil, N = N)$r, r_true, tolerance = 1e-6)
})

test_that("find_mle approaches the Poisson estimate at low counts", {
  counts <- c(2, 3, 4)
  dil <- rep(1e-6, 3)
  mle <- find_mle(counts, dil, N = 1e6)$r
  expect_equal(
    mle,
    find_poisson_cutoff(counts, dil, N = 50)$r,
    tolerance = 1e-3
  )
})

test_that("find_mle handles edge cases", {
  expect_equal(find_mle(c(5, 200), c(1e-5, 1e-5), N = 100)$r, Inf)
  zero <- find_mle(c(0, 0), c(1e-5, 1e-5), N = 100)
  expect_equal(zero$r, 0)
  expect_equal(zero$stderr, NA_real_)
})

test_that("input validation errors are informative", {
  expect_snapshot(find_poisson_cutoff(1:2, 1, N = 50), error = TRUE)
  expect_snapshot(find_mle(-1, 1e-5, N = 50), error = TRUE)
  expect_snapshot(find_mle("a", 1e-5, N = 50), error = TRUE)
})

read_cfu_example <- function() {
  path <- system.file("extdata", "cfu_example.tsv", package = "CFUestimatoR")
  read.delim(path)
}

test_that("estimate_cfu uses the MLE per group by default", {
  cfu <- read_cfu_example()
  res <- estimate_cfu(
    cfu,
    count = colony_count,
    group = c(measurement, time_hours, growth_media, sample_or_strain),
    N = 5000,
    sampling_vol = vol_plated_ml
  )

  expect_named(
    res,
    c(
      "measurement",
      "time_hours",
      "growth_media",
      "sample_or_strain",
      "estimator",
      "CFU_ml",
      "CFU_ml_stderr"
    )
  )
  expect_equal(nrow(res), 7)
  expect_equal(unique(res$estimator), "ML_estimator")

  e1 <- cfu[cfu$sample_or_strain == "E1" & cfu$growth_media == "PPY", ]
  row <- res[res$sample_or_strain == "E1" & res$growth_media == "PPY", ]
  expect_equal(
    row$CFU_ml,
    find_mle(e1$colony_count, e1$dilution, N = 5000)$r / 0.1
  )
})

test_that("estimate_cfu uses the per-group volume from a column", {
  cfu <- read_cfu_example()
  res <- estimate_cfu(
    cfu,
    count = colony_count,
    group = c(sample_or_strain, growth_media),
    N = 5000,
    sampling_vol = vol_plated_ml
  )

  ctrl <- cfu[cfu$sample_or_strain == "control", ]
  expect_equal(
    res$CFU_ml[res$sample_or_strain == "control"],
    find_mle(ctrl$colony_count, ctrl$dilution, N = 5000)$r / 0.08
  )
  expect_equal(
    estimate_cfu(
      cfu,
      count = colony_count,
      group = c(sample_or_strain, growth_media),
      N = 5000,
      sampling_vol = "vol_plated_ml"
    ),
    res
  )
})

test_that("estimate_cfu can use the Poisson cutoff estimator", {
  cfu <- read_cfu_example()
  res <- estimate_cfu(
    cfu,
    count = colony_count,
    group = c(sample_or_strain, growth_media),
    method = "poisson",
    N = 500,
    sampling_vol = vol_plated_ml
  )

  expect_equal(nrow(res), 7)
  expect_equal(unique(res$estimator), "poisson_cutoff")

  e1 <- cfu[cfu$sample_or_strain == "E1" & cfu$growth_media == "PPY", ]
  row <- res[res$sample_or_strain == "E1" & res$growth_media == "PPY", ]
  expect_equal(
    row$CFU_ml,
    find_poisson_cutoff(e1$colony_count, e1$dilution, N = 500)$r / 0.1
  )
})

test_that("estimate_cfu accepts a scalar volume from a variable", {
  d <- data.frame(count = c(2, 4, 3), dilution = 1e-6)
  vol <- 0.5
  expect_equal(
    estimate_cfu(d, sampling_vol = vol),
    estimate_cfu(d, sampling_vol = 0.5)
  )
})

test_that("estimate_cfu supports a grouping column named type", {
  d <- data.frame(
    count = c(2, 4, 3, 5),
    dilution = 1e-6,
    type = c("a", "a", "b", "b")
  )
  res <- estimate_cfu(d, group = type)
  expect_named(res, c("type", "estimator", "CFU_ml", "CFU_ml_stderr"))
  expect_equal(res$type, c("a", "b"))
})

test_that("estimate_cfu returns Inf when counts exceed N", {
  cfu <- read_cfu_example()
  res <- estimate_cfu(
    cfu,
    count = colony_count,
    group = sample_or_strain,
    N = 100,
    sampling_vol = vol_plated_ml
  )
  expect_equal(unique(res$CFU_ml), Inf)
})

test_that("defaults are N = 5000 (mle), N = 300 (poisson), 0.1 ml", {
  d <- data.frame(
    count = c(5, 8, 4, 62, 77, 51),
    dilution = rep(c(1e-7, 1e-6), each = 3)
  )
  mle <- find_mle(d$count, d$dilution, N = 5000)
  pois <- find_poisson_cutoff(d$count, d$dilution, N = 300)
  expect_equal(find_mle(d$count, d$dilution), mle)
  expect_equal(find_poisson_cutoff(d$count, d$dilution), pois)
  expect_equal(estimate_cfu(d)$CFU_ml, mle$r / 0.1)
  expect_equal(estimate_cfu(d, method = "poisson")$CFU_ml, pois$r / 0.1)
})

test_that("estimate_cfu group accepts tidyselect expressions", {
  d <- data.frame(
    count = c(2, 4, 3, 5, 6, 7, 8, 9),
    dilution = 1e-6,
    strain = rep(c("a", "b"), each = 4),
    rep = rep(c("x", "y"), 4)
  )
  bare <- estimate_cfu(d, group = strain)
  expect_equal(bare$strain, c("a", "b"))
  expect_equal(estimate_cfu(d, group = "strain"), bare)
  expect_equal(estimate_cfu(d, group = all_of("strain")), bare)

  multi <- estimate_cfu(d, group = c(strain, rep))
  expect_equal(multi$strain, c("a", "a", "b", "b"))
  expect_equal(multi$rep, c("x", "y", "x", "y"))
  expect_equal(estimate_cfu(d, group = where(is.character)), multi)
  expect_equal(estimate_cfu(d, group = starts_with("s")), bare)
})

test_that("estimate_cfu accepts bare or string count and dilution columns", {
  d <- data.frame(n = c(2, 4, 3), d = 1e-6)
  res <- estimate_cfu(d, count = n, dilution = d, sampling_vol = 1)
  expect_equal(
    estimate_cfu(d, count = "n", dilution = "d", sampling_vol = 1),
    res
  )
  expect_equal(res$CFU_ml, find_mle(d$n, d$d, N = 5000)$r)
})

test_that("estimate_cfu works with tibbles", {
  d <- tibble::tibble(count = c(2, 4, 3), dilution = 1e-6, g = "a")
  expect_equal(estimate_cfu(d, group = g)$g, "a")
})

test_that("estimate_cfu works without grouping", {
  d <- data.frame(n = c(2, 4, 3), d = 1e-6)
  res <- estimate_cfu(d, count = n, dilution = d, sampling_vol = 1)
  expect_equal(nrow(res), 1)
  expect_equal(res$CFU_ml, find_mle(d$n, d$d, N = 5000)$r)
})

test_that("estimate_cfu validates arguments", {
  d <- data.frame(count = 1, dilution = 1e-5)
  expect_snapshot(estimate_cfu(d, group = nope), error = TRUE)
  expect_snapshot(estimate_cfu(d, count = nope), error = TRUE)
  expect_snapshot(estimate_cfu(d, count = c(count, dilution)), error = TRUE)
  d$CFU_ml <- "a"
  expect_snapshot(estimate_cfu(d, group = CFU_ml), error = TRUE)
  expect_snapshot(estimate_cfu(d, sampling_vol = -1), error = TRUE)
  expect_snapshot(estimate_cfu(d, sampling_vol = nope), error = TRUE)
  varying <- data.frame(count = c(1, 2), dilution = 1e-5, v = c(0.1, 0.2))
  expect_snapshot(estimate_cfu(varying, sampling_vol = v), error = TRUE)
  expect_snapshot(estimate_cfu(d, method = "bogus"), error = TRUE)
})

test_that("results match the original web implementation", {
  # Reference values from https://huggingface.co/spaces/KMichaelMartini/CFUestimator
  # with sampling volume 1, max colony number N = 5000 (MLE) and cutoff
  # M = 300 (Poisson with cutoff). The naive Poisson estimate is the
  # cutoff estimator with no cutoff (N = Inf).
  d <- data.frame(
    dilution = rep(c(1e-7, 1e-6, 1e-5), each = 3),
    count = c(5, 8, 4, 62, 77, 51, 456, 597, 356)
  )

  mle <- find_mle(d$count, d$dilution, N = 5000, V = 1)
  expect_equal(mle$r, 50741066.7853, tolerance = 1e-9)
  expect_equal(mle$stderr, 1262704.1279, tolerance = 1e-9)

  naive <- find_poisson_cutoff(d$count, d$dilution, N = Inf, V = 1)
  expect_equal(naive$r, 48528528.5285, tolerance = 1e-9)
  expect_equal(naive$stderr, 1207192.2668, tolerance = 1e-9)

  cutoff <- find_poisson_cutoff(d$count, d$dilution, N = 300, V = 1)
  expect_equal(cutoff$r, 62727272.7273, tolerance = 1e-9)
  expect_equal(cutoff$stderr, 4359846.8394, tolerance = 1e-9)

  # Also through the high-level interface
  res <- estimate_cfu(d, sampling_vol = 1, N = 5000)
  expect_equal(res$CFU_ml, 50741066.7853, tolerance = 1e-9)
  res <- estimate_cfu(d, method = "poisson", sampling_vol = 1, N = 300)
  expect_equal(res$CFU_ml_stderr, 4359846.8394, tolerance = 1e-9)
})
