#' Show method for ClonalSimData
#'
#' @param object ClonalSimData object
#' @return NULL (prints to console)
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#' show(sim)
#'
setMethod("show", "ClonalSimData", function(object) {
  cat("ClonalSimData object\n")
  cat("==========================================\n")
  cat("Number of mutations:", nrow(object@mutations), "\n")

  if (length(object@params$subclone_freqs) > 0) {
    cat("Number of clones:", length(object@params$subclone_freqs), "\n")
    cat("Clone frequencies:", paste(round(object@params$subclone_freqs, 3),
                                     collapse = ", "), "\n")
    cat("Tumor purity:", round(sum(object@params$subclone_freqs), 3), "\n")
  }

  if (nrow(object@mutations) > 0) {
    cat("\nMutation types:\n")
    print(table(object@mutations$Type))

    cat("\nSequencing depth:\n")
    cat("  Mean:", round(mean(object@mutations$Depth), 2), "\n")
    cat("  Range:", min(object@mutations$Depth), "-",
        max(object@mutations$Depth), "\n")

    cat("\nVAF summary (observed):\n")
    print(summary(object@mutations$VAF))
  }

  if (length(object@metadata) > 0) {
    cat("\nMetadata:\n")
    if (!is.null(object@metadata$date)) {
      cat("  Created:", object@metadata$date, "\n")
    }
    if (!is.null(object@metadata$version)) {
      cat("  Package version:", object@metadata$version, "\n")
    }
  }
})

#' Summary method for ClonalSimData
#'
#' @param object ClonalSimData object
#' @return list with summary statistics
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#' summary(sim)
#'
setMethod("summary", "ClonalSimData", function(object) {
  stats <- list(
    n_mutations = nrow(object@mutations),
    n_clones = length(object@params$subclone_freqs),
    tumor_purity = sum(object@params$subclone_freqs),
    mutation_types = table(object@mutations$Type),
    vaf_summary = list(
      true = summary(object@mutations$True_VAF),
      observed = summary(object@mutations$VAF)
    ),
    depth_summary = list(
      mean = mean(object@mutations$Depth),
      median = median(object@mutations$Depth),
      sd = sd(object@mutations$Depth),
      range = range(object@mutations$Depth)
    )
  )

  class(stats) <- c("summary.ClonalSimData", "list")
  stats
})

#' Print summary method
#'
#' @param x summary.ClonalSimData object
#' @param ... additional arguments (unused)
#' @return NULL (prints to console)
#' @export
#'
print.summary.ClonalSimData <- function(x, ...) {
  cat("ClonalSim Summary\n")
  cat("==========================================\n")
  cat("Total mutations:", x$n_mutations, "\n")
  cat("Number of clones:", x$n_clones, "\n")
  cat("Tumor purity:", round(x$tumor_purity, 3), "\n")

  cat("\nMutation type counts:\n")
  print(x$mutation_types)

  cat("\nTrue VAF (biological):\n")
  print(x$vaf_summary$true)

  cat("\nObserved VAF (with sequencing noise):\n")
  print(x$vaf_summary$observed)

  cat("\nSequencing depth:\n")
  cat("  Mean:", round(x$depth_summary$mean, 2), "\n")
  cat("  Median:", x$depth_summary$median, "\n")
  cat("  SD:", round(x$depth_summary$sd, 2), "\n")
  cat("  Range:", x$depth_summary$range[1], "-", x$depth_summary$range[2], "\n")

  invisible(x)
}

#' Plot method for ClonalSimData
#'
#' @param x ClonalSimData object
#' @param y unused (for S4 compatibility)
#' @param type character indicating plot type: "vaf_density", "vaf_scatter",
#'   "depth_histogram", "clone_matrix". Default is "vaf_density"
#' @param ... additional arguments passed to plotting functions
#' @return ggplot2 object
#' @export
#'
#' @examples
#' sim <- simulateTumor(subclone_freqs = c(0.3, 0.4, 0.3))
#' plot(sim, type = "vaf_density")
#' plot(sim, type = "vaf_scatter")
#'
setMethod("plot", signature(x = "ClonalSimData", y = "missing"),
          function(x, y, type = "vaf_density", ...) {
  mutations <- x@mutations
  params <- x@params

  if (type == "vaf_density") {
    p <- ggplot(mutations, aes(x = VAF)) +
      geom_density(fill = "#984EA3", alpha = 0.6, color = "#984EA3", linewidth = 1) +
      geom_rug(alpha = 0.3, color = "#984EA3") +
      theme_minimal() +
      labs(title = "VAF Density Plot - Mixed Tumor Sample",
           subtitle = paste0(nrow(mutations), " mutations from ",
                           length(params$subclone_freqs), " subclones"),
           x = "Variant Allele Frequency (VAF)",
           y = "Density")

    # Add vertical lines for expected clone frequencies
    if (length(params$subclone_freqs) > 0) {
      p <- p +
        geom_vline(xintercept = params$subclone_freqs,
                   linetype = "dashed", color = "red", alpha = 0.5, linewidth = 0.8) +
        geom_vline(xintercept = sum(params$subclone_freqs),
                   linetype = "dashed", color = "darkred", alpha = 0.7, linewidth = 1)
    }

  } else if (type == "vaf_scatter") {
    type_colors <- c("founder" = "#E41A1C", "shared" = "#377EB8", "private" = "#4DAF4A")
    p <- ggplot(mutations, aes(x = seq_len(nrow(mutations)), y = VAF, color = Type)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = type_colors) +
      theme_minimal() +
      labs(title = "Mutational Profile: VAF of All Mutations",
           x = "Mutation Index",
           y = "Variant Allele Frequency (VAF)",
           color = "Type")

    if (length(params$subclone_freqs) > 0) {
      p <- p +
        geom_hline(yintercept = params$subclone_freqs,
                   linetype = "dashed", alpha = 0.3, color = "gray50")
    }

  } else if (type == "depth_histogram") {
    p <- ggplot(mutations, aes(x = Depth)) +
      geom_histogram(bins = 30, fill = "#377EB8", alpha = 0.7) +
      theme_minimal() +
      labs(title = "Sequencing Depth Distribution",
           x = "Depth (read count)",
           y = "Number of mutations")

  } else if (type == "clone_matrix") {
    # Prepare clone matrix data
    if (!"Clone_IDs" %in% colnames(mutations)) {
      stop("Clone_IDs column missing from mutations")
    }

    n_clones <- length(params$subclone_freqs)
    clone_cols <- paste0("Clone", seq_len(n_clones))

    clonal_matrix <- mutations[, c("Mutation", "VAF", "Type", "Clone_IDs")]
    for (i in seq_len(n_clones)) {
      clonal_matrix[[paste0("Clone", i)]] <- grepl(as.character(i),
                                                     clonal_matrix$Clone_IDs)
    }

    clonal_matrix <- clonal_matrix[order(-clonal_matrix$VAF), ]
    clonal_matrix$Index <- seq_len(nrow(clonal_matrix))

    matrix_long <- pivot_longer(clonal_matrix,
                                 cols = all_of(clone_cols),
                                 names_to = "Clone",
                                 values_to = "Present")

    p <- ggplot(matrix_long, aes(x = Clone, y = Index, fill = Present)) +
      geom_tile(color = "white", linewidth = 0.1) +
      scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#756bb1"),
                        labels = c("Absent", "Present")) +
      theme_minimal() +
      labs(title = "Mutation Presence Matrix in Subclones",
           subtitle = "Ordered by decreasing VAF",
           x = "Subclone",
           y = "Mutation (ordered by VAF)",
           fill = "Status") +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank())

  } else {
    stop("Unknown plot type. Choose from: vaf_density, vaf_scatter, ",
         "depth_histogram, clone_matrix")
  }

  return(p)
})
