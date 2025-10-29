test_that("applyBiologicalNoise works correctly", {
  # Basic functionality
  vaf <- applyBiologicalNoise(0.5, 100, concentration = 50)
  expect_length(vaf, 100)
  expect_true(all(vaf >= 0 & vaf <= 1))
  expect_true(abs(mean(vaf) - 0.5) < 0.1)  # Should be close to 0.5

  # Zero mutations
  vaf_zero <- applyBiologicalNoise(0.5, 0, concentration = 50)
  expect_length(vaf_zero, 0)

  # High concentration (low variability)
  vaf_high_conc <- applyBiologicalNoise(0.5, 1000, concentration = 200)
  expect_true(sd(vaf_high_conc) < 0.05)  # Low spread

  # Low concentration (high variability)
  vaf_low_conc <- applyBiologicalNoise(0.5, 1000, concentration = 10)
  expect_true(sd(vaf_low_conc) > 0.05)  # Higher spread

  # Edge cases
  expect_error(applyBiologicalNoise(-0.1, 10))  # Invalid frequency
  expect_error(applyBiologicalNoise(1.5, 10))   # Invalid frequency
  expect_error(applyBiologicalNoise(0.5, 10, concentration = -1))  # Invalid concentration
})

test_that("simulateDepth works correctly", {
  # Negative binomial
  depth_nb <- simulateDepth(100, mean_depth = 100, distribution = "negative_binomial", dispersion = 20)
  expect_length(depth_nb, 100)
  expect_true(all(depth_nb >= 10))  # Minimum depth enforced
  expect_true(abs(mean(depth_nb) - 100) < 20)  # Close to mean

  # Poisson
  depth_pois <- simulateDepth(100, mean_depth = 100, distribution = "poisson")
  expect_length(depth_pois, 100)
  expect_true(all(depth_pois >= 10))

  # Uniform
  depth_uni <- simulateDepth(100, mean_depth = 100, distribution = "uniform")
  expect_length(depth_uni, 100)
  expect_true(all(depth_uni == 100))  # All equal to mean

  # Errors
  expect_error(simulateDepth(10, mean_depth = -1))  # Negative depth
  expect_error(simulateDepth(10, distribution = "invalid"))  # Invalid distribution
})

test_that("simulateSequencingReads works correctly", {
  true_vaf <- rep(0.3, 100)
  depth <- rep(100, 100)

  reads <- simulateSequencingReads(true_vaf, depth, error_rate = 0.001)

  # Check structure
  expect_type(reads, "list")
  expect_named(reads, c("vaf", "alt_reads", "ref_reads"))
  expect_length(reads$vaf, 100)
  expect_length(reads$alt_reads, 100)
  expect_length(reads$ref_reads, 100)

  # Check values
  expect_true(all(reads$vaf >= 0 & reads$vaf <= 1))
  expect_true(all(reads$alt_reads >= 0))
  expect_true(all(reads$ref_reads >= 0))
  expect_true(all(reads$alt_reads <= depth))
  expect_true(all(reads$alt_reads + reads$ref_reads == depth))

  # VAF should be close to true_vaf
  expect_true(abs(mean(reads$vaf) - 0.3) < 0.05)

  # Errors
  expect_error(simulateSequencingReads(c(0.3, 0.4), c(100)))  # Length mismatch
  expect_error(simulateSequencingReads(c(-0.1, 0.5), c(100, 100)))  # Invalid VAF
  expect_error(simulateSequencingReads(c(0.3, 0.4), c(-10, 100)))  # Negative depth
  expect_error(simulateSequencingReads(c(0.3, 0.4), c(100, 100), error_rate = -0.01))  # Negative error rate
})
