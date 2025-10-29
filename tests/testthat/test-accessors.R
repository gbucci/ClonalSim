test_that("Accessor functions work correctly", {
  sim <- simulateTumor(seed = 123)

  # getMutations
  mutations <- getMutations(sim)
  expect_s3_class(mutations, "data.frame")
  expect_gt(nrow(mutations), 0)

  # getSimParams
  params <- getSimParams(sim)
  expect_type(params, "list")
  expect_true("subclone_freqs" %in% names(params))

  # getTrueVAF
  true_vaf <- getTrueVAF(sim)
  expect_type(true_vaf, "double")
  expect_equal(length(true_vaf), nrow(mutations))
  expect_true(all(true_vaf >= 0 & true_vaf <= 1))

  # getObservedVAF
  obs_vaf <- getObservedVAF(sim)
  expect_type(obs_vaf, "double")
  expect_equal(length(obs_vaf), nrow(mutations))
  expect_true(all(obs_vaf >= 0 & obs_vaf <= 1))

  # getClonalStructure
  structure <- getClonalStructure(sim)
  expect_s3_class(structure, "data.frame")
  expect_true("Clone" %in% colnames(structure))
  expect_true("Frequency" %in% colnames(structure))

  # getMetadata
  metadata <- getMetadata(sim)
  expect_type(metadata, "list")
})

test_that("Accessors fail on wrong object type", {
  wrong_object <- list(mutations = data.frame())

  expect_error(getMutations(wrong_object))
  expect_error(getSimParams(wrong_object))
  expect_error(getTrueVAF(wrong_object))
  expect_error(getObservedVAF(wrong_object))
  expect_error(getClonalStructure(wrong_object))
})
