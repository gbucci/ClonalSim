test_that("show method works", {
  set.seed(123)
  sim <- simulateTumor()

  # Capture output
  output <- capture.output(show(sim))

  expect_true(any(grepl("ClonalSimData", output)))
  expect_true(any(grepl("Number of mutations", output)))
  expect_true(any(grepl("Clone frequencies", output)))
})

test_that("summary method works", {
  set.seed(123)
  sim <- simulateTumor()

  summ <- summary(sim)

  expect_type(summ, "list")
  expect_s3_class(summ, "summary.ClonalSimData")
  expect_true("n_mutations" %in% names(summ))
  expect_true("n_clones" %in% names(summ))
  expect_true("tumor_purity" %in% names(summ))
  expect_true("vaf_summary" %in% names(summ))
  expect_true("depth_summary" %in% names(summ))

  expect_gt(summ$n_mutations, 0)
  expect_gt(summ$n_clones, 0)
  expect_true(summ$tumor_purity >= 0 & summ$tumor_purity <= 1)
})

test_that("plot method creates ggplot objects", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  sim <- simulateTumor()

  # VAF density plot
  p1 <- plot(sim, type = "vaf_density")
  expect_s3_class(p1, "gg")
  expect_s3_class(p1, "ggplot")

  # VAF scatter plot
  p2 <- plot(sim, type = "vaf_scatter")
  expect_s3_class(p2, "gg")

  # Depth histogram
  p3 <- plot(sim, type = "depth_histogram")
  expect_s3_class(p3, "gg")

  # Clone matrix
  p4 <- plot(sim, type = "clone_matrix")
  expect_s3_class(p4, "gg")

  # Invalid type
  expect_error(plot(sim, type = "invalid_type"))
})
