test_that("simulateTumor creates valid ClonalSimData object", {
  set.seed(123)
  sim <- simulateTumor(
    subclone_freqs = c(0.3, 0.4, 0.3),
    n_mut_per_clone = c(20, 25, 15)
  )

  expect_s4_class(sim, "ClonalSimData")
  expect_true(validObject(sim))
})

test_that("simulateTumor with default parameters works", {
  set.seed(123)
  sim <- simulateTumor()

  expect_s4_class(sim, "ClonalSimData")
  mutations <- getMutations(sim)
  expect_gt(nrow(mutations), 0)
  expect_true(all(c("Mutation", "Chromosome", "Position", "Ref", "Alt",
                    "True_VAF", "VAF", "Depth", "Alt_reads") %in% colnames(mutations)))
})

test_that("simulateTumor respects input parameters", {
  set.seed(123)
  sim <- simulateTumor(
    subclone_freqs = c(0.2, 0.3),
    n_mut_per_clone = c(10, 15),
    n_mut_founder = 5
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
  set.seed(123)
  sim_no_noise <- suppressWarnings(simulateTumor(
    subclone_freqs = c(0.5, 0.5),
    n_mut_per_clone = c(30, 30),
    biological_noise = list(enabled = FALSE),
    sequencing_noise = list(enabled = FALSE)
  ))

  mutations <- getMutations(sim_no_noise)
  # Without noise, True_VAF should equal VAF
  expect_equal(mutations$True_VAF, mutations$VAF)

  # High biological noise
  set.seed(123)
  sim_high_bio <- suppressWarnings(simulateTumor(
    subclone_freqs = c(0.5, 0.5),
    n_mut_per_clone = c(30, 30),
    biological_noise = list(enabled = TRUE, concentration = 10)
  ))
  mutations_high <- getMutations(sim_high_bio)
  expect_gt(sd(mutations_high$True_VAF), 0.05)

  # Low biological noise
  set.seed(123)
  sim_low_bio <- suppressWarnings(simulateTumor(
    subclone_freqs = c(0.5, 0.5),
    n_mut_per_clone = c(30, 30),
    biological_noise = list(enabled = TRUE, concentration = 200)
  ))
  mutations_low <- getMutations(sim_low_bio)
  expect_lt(sd(mutations_low$True_VAF), sd(mutations_high$True_VAF))
})

test_that("simulateTumor is reproducible with seed", {
  set.seed(42)
  sim1 <- simulateTumor()
  set.seed(42)
  sim2 <- simulateTumor()

  mutations1 <- getMutations(sim1)
  mutations2 <- getMutations(sim2)

  expect_equal(mutations1$VAF, mutations2$VAF)
  expect_equal(mutations1$Depth, mutations2$Depth)
  expect_equal(mutations1$Alt_reads, mutations2$Alt_reads)
})

test_that("simulateTumor creates correct mutation types", {
  set.seed(123)
  sim <- simulateTumor()
  mutations <- getMutations(sim)

  expect_true("Type" %in% colnames(mutations))
  expect_true(all(mutations$Type %in% c("founder", "shared", "private", "germline")))

  # Should have three somatic types with default parameters (germline disabled by default)
  types <- table(mutations$Type)
  expect_gt(types["founder"], 0)
  expect_gt(types["shared"], 0)
  expect_gt(types["private"], 0)
})

# Germline variant tests
test_that("simulateTumor generates germline variants when enabled", {
  set.seed(123)
  sim <- simulateTumor(
    subclone_freqs = c(0.7),
    n_mut_per_clone = c(50),
    germline_variants = list(enabled = TRUE, n_variants = 100)
  )

  mutations <- getMutations(sim)
  germline_muts <- mutations[mutations$Type == "germline", ]

  # Check germline mutations were generated
  expect_equal(nrow(germline_muts), 100)

  # Check germline VAF is around 0.5 (allowing for noise)
  mean_germline_vaf <- mean(germline_muts$VAF)
  expect_gt(mean_germline_vaf, 0.4)
  expect_lt(mean_germline_vaf, 0.6)
})

test_that("germline VAF stays at ~0.5 regardless of tumor purity", {
  # High purity (70%)
  set.seed(456)
  sim_high <- simulateTumor(
    subclone_freqs = c(0.7),
    n_mut_per_clone = c(30),
    germline_variants = list(enabled = TRUE, n_variants = 100)
  )

  # Low purity (30%)
  set.seed(789)
  sim_low <- suppressWarnings(simulateTumor(
    subclone_freqs = c(0.1, 0.2),
    n_mut_per_clone = c(20, 30),
    germline_variants = list(enabled = TRUE, n_variants = 100)
  ))

  germline_high <- getMutations(sim_high)[getMutations(sim_high)$Type == "germline", ]
  germline_low <- getMutations(sim_low)[getMutations(sim_low)$Type == "germline", ]

  mean_vaf_high <- mean(germline_high$VAF)
  mean_vaf_low <- mean(germline_low$VAF)

  # Both should be around 0.5
  expect_gt(mean_vaf_high, 0.4)
  expect_lt(mean_vaf_high, 0.6)
  expect_gt(mean_vaf_low, 0.4)
  expect_lt(mean_vaf_low, 0.6)

  # Should be similar despite different purity
  expect_lt(abs(mean_vaf_high - mean_vaf_low), 0.1)
})

test_that("germline variants can be disabled", {
  sim <- simulateTumor(
    germline_variants = list(enabled = FALSE)
  )

  mutations <- getMutations(sim)
  germline_muts <- mutations[mutations$Type == "germline", ]

  # Should have no germline mutations
  expect_equal(nrow(germline_muts), 0)
})

test_that("germline variants respect custom parameters", {
  set.seed(123)
  sim <- simulateTumor(
    subclone_freqs = c(0.5),
    n_mut_per_clone = c(30),
    germline_variants = list(
      enabled = TRUE,
      n_variants = 50,
      vaf_expected = 0.5
    )
  )

  mutations <- getMutations(sim)
  germline_muts <- mutations[mutations$Type == "germline", ]

  # Check count
  expect_equal(nrow(germline_muts), 50)

  # Check VAF
  mean_vaf <- mean(germline_muts$VAF)
  expect_gt(mean_vaf, 0.4)
  expect_lt(mean_vaf, 0.6)
})

test_that("shared mutations warning for non-existent clones", {
  # Should warn when default n_mut_shared references non-existent clones
  expect_warning(
    simulateTumor(
      subclone_freqs = c(0.7),  # Only 1 clone
      n_mut_per_clone = c(50)
    ),
    "references non-existent clones"
  )
})
