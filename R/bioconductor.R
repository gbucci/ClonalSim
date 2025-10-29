#' Convert ClonalSimData to GRanges object
#'
#' @description Converts mutation data to a GenomicRanges::GRanges object,
#'   the standard Bioconductor data structure for genomic intervals.
#'
#' @param object ClonalSimData object
#' @param include_metadata logical, include all metadata columns. Default: TRUE
#'
#' @return GRanges object with mutation information
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#' gr <- toGRanges(sim)
#' print(gr)
#'
#' # Access metadata columns
#' GenomicRanges::mcols(gr)
#'
#' # Subset by chromosome
#' gr_chr1 <- gr[GenomicRanges::seqnames(gr) == "chr1"]
#'
toGRanges <- function(object, include_metadata = TRUE) {
  stopifnot(inherits(object, "ClonalSimData"))

  mutations <- object@mutations

  # Create GRanges object
  gr <- makeGRangesFromDataFrame(
    mutations,
    seqnames.field = "Chromosome",
    start.field = "Position",
    end.field = "Position",
    strand = "*",
    keep.extra.columns = include_metadata
  )

  # Add additional metadata if not already included
  if (include_metadata) {
    # Ensure all important columns are in metadata
    metadata_cols <- c("Mutation", "Ref", "Alt", "True_VAF", "VAF",
                       "Depth", "Alt_reads", "Clone", "Type", "Clone_IDs")
    existing_cols <- colnames(GenomicRanges::mcols(gr))
    missing_cols <- setdiff(metadata_cols, existing_cols)

    for (col in missing_cols) {
      if (col %in% colnames(mutations)) {
        GenomicRanges::mcols(gr)[[col]] <- mutations[[col]]
      }
    }
  }

  return(gr)
}

#' Convert ClonalSimData to VCF format
#'
#' @description Exports mutation data to Variant Call Format (VCF), the
#'   standard format for variant data. Returns a VariantAnnotation::VRanges
#'   object which can be written to VCF file.
#'
#' @param object ClonalSimData object
#' @param sample_name character, sample name for VCF. Default: "TumorSample"
#' @param output_file character, optional file path to write VCF. If NULL,
#'   returns VRanges object without writing. Default: NULL
#'
#' @return VRanges object (invisibly if output_file is specified)
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#'
#' # Get VRanges object
#' vr <- toVCF(sim)
#' print(vr)
#'
#' # Write to VCF file
#' \dontrun{
#' toVCF(sim, output_file = "simulated_mutations.vcf")
#' }
#'
toVCF <- function(object, sample_name = "TumorSample", output_file = NULL) {
  stopifnot(inherits(object, "ClonalSimData"))

  mutations <- object@mutations

  # Create VRanges object
  vr <- VRanges(
    seqnames = Rle(mutations$Chromosome),
    ranges = IRanges(start = mutations$Position, width = 1),
    ref = mutations$Ref,
    alt = mutations$Alt,
    totalDepth = mutations$Depth,
    refDepth = mutations$Depth - mutations$Alt_reads,
    altDepth = mutations$Alt_reads,
    sampleNames = sample_name
  )

  # Add custom INFO fields as metadata
  S4Vectors::mcols(vr)$TRUE_VAF <- mutations$True_VAF
  S4Vectors::mcols(vr)$VAF <- mutations$VAF
  S4Vectors::mcols(vr)$CLONE <- mutations$Clone
  S4Vectors::mcols(vr)$TYPE <- mutations$Type
  S4Vectors::mcols(vr)$CLONE_IDS <- mutations$Clone_IDs

  # Write to file if specified
  if (!is.null(output_file)) {
    writeVcf(vr, filename = output_file)
    message("VCF file written to: ", output_file)
    return(invisible(vr))
  }

  return(vr)
}

#' Export mutation data to data frame
#'
#' @description Simple export to data.frame (useful for compatibility with
#'   other tools like PyClone, SciClone)
#'
#' @param object ClonalSimData object
#' @param file character, optional file path to write CSV. Default: NULL
#' @param include_true_vaf logical, include True_VAF column. Default: TRUE
#'
#' @return data.frame with mutation data
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#' df <- toDataFrame(sim)
#' head(df)
#'
#' # Export to CSV
#' \dontrun{
#' toDataFrame(sim, file = "mutations.csv")
#' }
#'
toDataFrame <- function(object, file = NULL, include_true_vaf = TRUE) {
  stopifnot(inherits(object, "ClonalSimData"))

  df <- object@mutations

  if (!include_true_vaf) {
    df <- df[, colnames(df) != "True_VAF"]
  }

  if (!is.null(file)) {
    write.csv(df, file = file, row.names = FALSE)
    message("Data frame written to: ", file)
    return(invisible(df))
  }

  return(df)
}

#' Convert to PyClone input format
#'
#' @description Exports data in format suitable for PyClone clonal deconvolution
#'   analysis.
#'
#' @param object ClonalSimData object
#' @param file character, file path to write TSV
#' @param sample_id character, sample identifier. Default: "sample1"
#'
#' @return data.frame formatted for PyClone (invisibly)
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#'
#' \dontrun{
#' toPyClone(sim, file = "pyclone_input.tsv")
#' }
#'
toPyClone <- function(object, file, sample_id = "sample1") {
  stopifnot(inherits(object, "ClonalSimData"))

  mutations <- object@mutations

  # PyClone format
  pyclone_df <- data.frame(
    mutation_id = mutations$Mutation,
    ref_counts = mutations$Depth - mutations$Alt_reads,
    var_counts = mutations$Alt_reads,
    normal_cn = 2,  # Assume diploid
    minor_cn = 0,   # Assume no CNAs
    major_cn = 2,   # Assume no CNAs
    sample_id = sample_id,
    stringsAsFactors = FALSE
  )

  write.table(pyclone_df, file = file, sep = "\t", row.names = FALSE,
              quote = FALSE)
  message("PyClone input written to: ", file)

  return(invisible(pyclone_df))
}

#' Convert to SciClone input format
#'
#' @description Exports data in format suitable for SciClone clonal analysis.
#'
#' @param object ClonalSimData object
#' @param file character, file path to write TSV
#'
#' @return data.frame formatted for SciClone (invisibly)
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#'
#' \dontrun{
#' toSciClone(sim, file = "sciclone_input.tsv")
#' }
#'
toSciClone <- function(object, file) {
  stopifnot(inherits(object, "ClonalSimData"))

  mutations <- object@mutations

  # SciClone format
  sciclone_df <- data.frame(
    chr = mutations$Chromosome,
    pos = mutations$Position,
    ref_reads = mutations$Depth - mutations$Alt_reads,
    var_reads = mutations$Alt_reads,
    vaf = mutations$VAF * 100,  # SciClone uses percentage
    stringsAsFactors = FALSE
  )

  write.table(sciclone_df, file = file, sep = "\t", row.names = FALSE,
              quote = FALSE)
  message("SciClone input written to: ", file)

  return(invisible(sciclone_df))
}
