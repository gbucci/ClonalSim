test_that("simulateTumor creates valid ClonalSimData object", {
  sim <- simulateTumor(
    subclone_freqs = c(0.3, 0.4, 0.3),
    n_mut_per_clone = c(20, 25, 15),
    seed = 123
  )

  expect_s4_class(sim, "ClonalSimData")
  expect_true(validObject(sim))
})

test_that("simulateTumor with default parameters works", {
  sim <- simulateTumor(seed = 123)

  expect_s4_class(sim, "ClonalSimData")
  mutations <- getMutations(sim)
  expect_gt(nrow(mutations), 0)
  expect_true(all(c("Mutation", "Chromosome", "Position", "Ref", "Alt",
                    "True_VAF", "VAF", "Depth", "Alt_reads") %in% colnames(mutations)))
})

test_that("simulateTumor respects input parameters", {
  sim <- simulateTumor(
    subclone_freqs = c(0.2, 0.3),
    n_mut_per_clone = c(10, 15),
    n_mut_founder = 5,
    seed = 123
  )

  params <- getSimParams(sim)
  expect_equal(params$subclone_freqs, c(Clone1 = 0.2, Clone2 = 0.3))
  expect_equal(params$n_mut_per_clone, c(10, 15))
  expect_equal(params$n_mut_founder, 5)

  mutations <- getMutations(sim)
  # Check total mutations: founder + shared + private
  expect_gt(nrow(mutations), 25)  # At least 5 + 10 + 15
})

test_that("simulateTumor validates input", {
  # Invalid frequencies
  expect_error(simulateTumor(subclone_freqs = c(-0.1, 0.5)))
  expect_error(simulateTumor(subclone_freqs = c(0.6, 0.6)))  # Sum > 1

  # Mismatched lengths
  expect_error(simulateTumor(
    subclone_freqs = c(0.3, 0.4),
    n_mut_per_clone = c(10, 20, 30)
  ))
})

test_that("simulateTumor noise parameters work", {
  # No noise
  sim_no_noise <- simulateTumor(
    subclone_freqs = c(0.5, 0.5),
    biological_noise = list(enabled = FALSE),
    sequencing_noise = list(enabled = FALSE),
    seed = 123
  )

  mutations <- getMutations(sim_no_noise)
  # Without noise, True_VAF should equal VAF
  expect_equal(mutations$True_VAF, mutations$VAF)

  # High biological noise
  sim_high_bio <- simulateTumor(
    subclone_freqs = c(0.5, 0.5),
    biological_noise = list(enabled = TRUE, concentration = 10),
    seed = 123
  )
  mutations_high <- getMutations(sim_high_bio)
  expect_gt(sd(mutations_high$True_VAF), 0.05)

  # Low biological noise
  sim_low_bio <- simulateTumor(
    subclone_freqs = c(0.5, 0.5),
    biological_noise = list(enabled = TRUE, concentration = 200),
    seed = 123
  )
  mutations_low <- getMutations(sim_low_bio)
  expect_lt(sd(mutations_low$True_VAF), sd(mutations_high$True_VAF))
})

test_that("simulateTumor is reproducible with seed", {
  sim1 <- simulateTumor(seed = 42)
  sim2 <- simulateTumor(seed = 42)

  mutations1 <- getMutations(sim1)
  mutations2 <- getMutations(sim2)

  expect_equal(mutations1$VAF, mutations2$VAF)
  expect_equal(mutations1$Depth, mutations2$Depth)
  expect_equal(mutations1$Alt_reads, mutations2$Alt_reads)
})

test_that("simulateTumor creates correct mutation types", {
  sim <- simulateTumor(seed = 123)
  mutations <- getMutations(sim)

  expect_true("Type" %in% colnames(mutations))
  expect_true(all(mutations$Type %in% c("founder", "shared", "private")))

  # Should have all three types with default parameters
  types <- table(mutations$Type)
  expect_gt(types["founder"], 0)
  expect_gt(types["shared"], 0)
  expect_gt(types["private"], 0)
})
