test_that("toGRanges creates valid GRanges object", {
  skip_if_not_installed("GenomicRanges")

  set.seed(123)
  sim <- simulateTumor()
  gr <- toGRanges(sim)

  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), nrow(getMutations(sim)))

  # Check metadata columns
  mcols_names <- colnames(GenomicRanges::mcols(gr))
  expect_true("VAF" %in% mcols_names)
  expect_true("Depth" %in% mcols_names)
  expect_true("Type" %in% mcols_names)
})

test_that("toGRanges without metadata works", {
  skip_if_not_installed("GenomicRanges")

  set.seed(123)
  sim <- simulateTumor()
  gr <- toGRanges(sim, include_metadata = FALSE)

  expect_s4_class(gr, "GRanges")
  # Should have fewer or no metadata columns
})

test_that("toVCF creates valid VRanges object", {
  skip_if_not_installed("VariantAnnotation")

  set.seed(123)
  sim <- simulateTumor()
  vr <- suppressWarnings(toVCF(sim, sample_name = "TestSample"))

  expect_s4_class(vr, "VRanges")
  # VRanges filters out mutations where ref == alt
  muts <- getMutations(sim)
  valid_muts <- muts[muts$Ref != muts$Alt, ]
  expect_equal(length(vr), nrow(valid_muts))

  # Check that depths are present
  expect_true("totalDepth" %in% names(S4Vectors::mcols(vr)) ||
                !is.null(VariantAnnotation::totalDepth(vr)))
})

test_that("toDataFrame works correctly", {
  set.seed(123)
  sim <- simulateTumor()
  df <- toDataFrame(sim)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), nrow(getMutations(sim)))
  expect_true("True_VAF" %in% colnames(df))

  # Without true VAF
  df_no_true <- toDataFrame(sim, include_true_vaf = FALSE)
  expect_false("True_VAF" %in% colnames(df_no_true))
})

test_that("Export functions work with file output", {
  set.seed(123)
  sim <- simulateTumor()
  tmp_file <- tempfile(fileext = ".csv")

  # toDataFrame with file - check invisibility
  expect_invisible(toDataFrame(sim, file = tmp_file))
  expect_true(file.exists(tmp_file))

  # Clean up
  unlink(tmp_file)
})

test_that("toPyClone creates correct format", {
  set.seed(123)
  sim <- simulateTumor()
  tmp_file <- tempfile(fileext = ".tsv")

  result <- toPyClone(sim, file = tmp_file, sample_id = "test")
  expect_true(file.exists(tmp_file))

  # Read and check format
  pyclone_data <- read.table(tmp_file, header = TRUE, sep = "\t")
  expect_true("mutation_id" %in% colnames(pyclone_data))
  expect_true("ref_counts" %in% colnames(pyclone_data))
  expect_true("var_counts" %in% colnames(pyclone_data))

  # Clean up
  unlink(tmp_file)
})

test_that("toSciClone creates correct format", {
  set.seed(123)
  sim <- simulateTumor()
  tmp_file <- tempfile(fileext = ".tsv")

  result <- toSciClone(sim, file = tmp_file)
  expect_true(file.exists(tmp_file))

  # Read and check format
  sciclone_data <- read.table(tmp_file, header = TRUE, sep = "\t")
  expect_true("chr" %in% colnames(sciclone_data))
  expect_true("pos" %in% colnames(sciclone_data))
  expect_true("vaf" %in% colnames(sciclone_data))

  # VAF should be in percentage
  expect_true(all(sciclone_data$vaf >= 0 & sciclone_data$vaf <= 100))

  # Clean up
  unlink(tmp_file)
})
