#!/usr/bin/env Rscript

# Simulation of mutational profile of a tumor sample with subclones
# Bioinformatics script: simulates allelic frequencies from 4 subclones
# Includes shared mutations according to evolutionary hierarchy

library(ggplot2)

set.seed(123)  # For reproducibility

# ===== SIMULATION PARAMETERS =====

# Frequencies of the 4 subclones (must sum to <= 1)
subclone_freqs <- c(0.15, 0.25, 0.30, 0.30)  # Subclone 1, 2, 3, 4
names(subclone_freqs) <- paste0("Clone", 1:4)

# Number of mutations per subclone (private mutations)
n_mut_per_clone <- c(20, 25, 30, 15)  # Unique mutations of each clone

# Number of founder mutations (present in all clones)
n_mut_founder <- 10

# Number of shared mutations between subgroups of clones
# Hierarchical structure: Clone1 -> Clone2,3 -> Clone4
n_mut_shared <- list(
  "2 3 4" = 15,     # Shared by Clone2, Clone3, Clone4
  "3 4" = 12,       # Shared by Clone3 and Clone4
  "1 2" = 8         # Shared by Clone1 and Clone2
)

# Technical noise (sequencing error)
noise_sd <- 0.02  # Standard deviation of gaussian noise

# ===== FUNCTIONS =====

# Function to generate mutations with allelic frequency
generate_mutations <- function(n_mut, base_freq, clone_ids, type = "private", noise_sd = 0.02) {

  if (n_mut == 0) return(NULL)

  # Generate mutation names
  if (type == "founder") {
    mut_names <- paste0("Founder_", 1:n_mut)
    clone_label <- "Founder"
  } else if (type == "shared") {
    clone_str <- paste(clone_ids, collapse = "_")
    mut_names <- paste0("Shared_C", clone_str, "_mut", 1:n_mut)
    clone_label <- paste0("Clone", paste(clone_ids, collapse = "+"))
  } else {
    mut_names <- paste0("Clone", clone_ids, "_mut", 1:n_mut)
    clone_label <- paste0("Clone", clone_ids)
  }

  # Base allelic frequency (with small biological variation)
  vaf <- rnorm(n_mut, mean = base_freq, sd = noise_sd)

  # Limit VAF between 0 and 1
  vaf <- pmax(0.01, pmin(0.99, vaf))

  # Create dataframe
  data.frame(
    Mutation = mut_names,
    VAF = vaf,
    Clone = clone_label,
    Type = type,
    Clone_IDs = paste(clone_ids, collapse = ","),
    stringsAsFactors = FALSE
  )
}

# ===== DATA GENERATION =====

# List to collect all mutations
mutation_list <- list()
idx <- 1

# 1. Founder mutations (present in all clones)
founder_freq <- sum(subclone_freqs)  # Sum of all frequencies
mutation_list[[idx]] <- generate_mutations(
  n_mut_founder,
  founder_freq,
  clone_ids = 1:4,
  type = "founder",
  noise_sd = noise_sd
)
idx <- idx + 1

# 2. Shared mutations between subgroups of clones
for (i in 1:length(n_mut_shared)) {
  # Extract the clones that share the mutations
  shared_clones <- as.numeric(strsplit(names(n_mut_shared)[i], " ")[[1]])

  # Calculate frequency as sum of frequencies of involved clones
  shared_freq <- sum(subclone_freqs[shared_clones])

  mutation_list[[idx]] <- generate_mutations(
    n_mut_shared[[i]],
    shared_freq,
    clone_ids = shared_clones,
    type = "shared",
    noise_sd = noise_sd
  )
  idx <- idx + 1
}

# 3. Private mutations for each subclone
for (i in 1:length(subclone_freqs)) {
  mutation_list[[idx]] <- generate_mutations(
    n_mut_per_clone[i],
    subclone_freqs[i],
    clone_ids = i,
    type = "private",
    noise_sd = noise_sd
  )
  idx <- idx + 1
}

# Combine all data (remove any NULLs)
mutation_list <- mutation_list[!sapply(mutation_list, is.null)]
complete_data <- do.call(rbind, mutation_list)

# Add additional columns for realism
complete_data$Chromosome <- sample(paste0("chr", 1:22), nrow(complete_data), replace = TRUE)
complete_data$Position <- sample(1e6:2e8, nrow(complete_data), replace = TRUE)
complete_data$Ref <- sample(c("A", "T", "C", "G"), nrow(complete_data), replace = TRUE)
complete_data$Alt <- sample(c("A", "T", "C", "G"), nrow(complete_data), replace = TRUE)

# Simulate sequencing coverage (depth)
complete_data$Depth <- rpois(nrow(complete_data), lambda = 100)
complete_data$Alt_reads <- round(complete_data$VAF * complete_data$Depth)

# Reorder columns
complete_data <- complete_data[, c("Mutation", "Chromosome", "Position",
                                    "Ref", "Alt", "VAF", "Depth", "Alt_reads",
                                    "Clone", "Type")]

# ===== OUTPUT =====

# Save to CSV file
write.csv(complete_data, "simulated_mutational_profile.csv", row.names = FALSE)

# Print summary
cat("\n===== SIMULATION SUMMARY =====\n")
cat("Subclone frequencies:\n")
print(subclone_freqs)
cat("\nTotal number of mutations:", nrow(complete_data), "\n")
cat("  - Founder mutations:", n_mut_founder, "\n")
cat("  - Shared mutations:\n")
for (i in 1:length(n_mut_shared)) {
  cat("    Clones", names(n_mut_shared)[i], ":", n_mut_shared[[i]], "\n")
}
cat("  - Private mutations:\n")
for (i in 1:length(n_mut_per_clone)) {
  cat("    Clone", i, ":", n_mut_per_clone[i], "\n")
}

# Statistics per clone
cat("\n===== MEAN VAF PER MUTATION TYPE =====\n")
print(aggregate(VAF ~ Type + Clone, data = complete_data, FUN = mean))

# Count mutations per type
cat("\n===== MUTATION COUNT PER TYPE =====\n")
print(table(complete_data$Type))

# ===== VISUALIZATION =====

# Define colors for mutation types
type_colors <- c("founder" = "#E41A1C",
                 "shared" = "#377EB8",
                 "private" = "#4DAF4A")

# Plot 1: VAF distribution per mutation type
p1 <- ggplot(complete_data, aes(x = VAF, fill = Type)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  facet_wrap(~Type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = type_colors) +
  theme_minimal() +
  labs(title = "VAF Distribution per Mutation Type",
       x = "Variant Allele Frequency (VAF)",
       y = "Number of mutations") +
  theme(legend.position = "none")

ggsave("VAF_distribution_per_type.png", p1, width = 10, height = 8)

# Plot 2: VAF scatter plot colored by type
p2 <- ggplot(complete_data, aes(x = 1:nrow(complete_data), y = VAF, color = Type)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = type_colors) +
  theme_minimal() +
  labs(title = "Mutational Profile: VAF of All Mutations",
       x = "Mutation Index",
       y = "Variant Allele Frequency (VAF)",
       color = "Type") +
  geom_hline(yintercept = subclone_freqs, linetype = "dashed", alpha = 0.3, color = "gray50") +
  annotate("text", x = nrow(complete_data) * 0.95,
           y = subclone_freqs,
           label = paste0("Clone", 1:4),
           hjust = 1, vjust = -0.5, size = 3, color = "gray30")

ggsave("VAF_scatter_plot.png", p2, width = 12, height = 6)

# Plot 3: CUMULATIVE DENSITY PLOT - mixed sample
# This simulates what you would see when sequencing the actual tumor DNA
p3 <- ggplot(complete_data, aes(x = VAF)) +
  geom_density(fill = "#984EA3", alpha = 0.6, color = "#984EA3", linewidth = 1) +
  geom_rug(alpha = 0.3, color = "#984EA3") +
  theme_minimal() +
  labs(title = "Cumulative Density Plot - Mixed Tumor Sample",
       subtitle = paste0("Simulation of ", nrow(complete_data),
                        " mutations from 4 subclones with frequencies: ",
                        paste(subclone_freqs, collapse = ", ")),
       x = "Variant Allele Frequency (VAF)",
       y = "Density") +
  # Add vertical lines for expected clone frequencies
  geom_vline(xintercept = subclone_freqs,
             linetype = "dashed",
             color = "red",
             alpha = 0.5,
             linewidth = 0.8) +
  geom_vline(xintercept = sum(subclone_freqs),
             linetype = "dashed",
             color = "darkred",
             alpha = 0.7,
             linewidth = 1) +
  annotate("text",
           x = c(subclone_freqs, sum(subclone_freqs)),
           y = max(density(complete_data$VAF)$y) * 0.9,
           label = c(paste0("C", 1:4), "Founder"),
           angle = 90,
           vjust = -0.5,
           size = 3.5,
           color = "red")

ggsave("VAF_cumulative_density.png", p3, width = 12, height = 7)

# Plot 4: Violin plot per mutation type
p4 <- ggplot(complete_data, aes(x = Type, y = VAF, fill = Type)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1) +
  scale_fill_manual(values = type_colors) +
  theme_minimal() +
  labs(title = "VAF Distribution per Mutation Type (Violin Plot)",
       x = "Mutation Type",
       y = "Variant Allele Frequency (VAF)") +
  theme(legend.position = "none")

ggsave("VAF_violin_plot.png", p4, width = 10, height = 7)

# Plot 5: Heatmap-like to visualize clonal hierarchy
# Create a presence/absence matrix for each clone
# First ensure Clone_IDs exists
if (!"Clone_IDs" %in% colnames(complete_data)) {
  # If it doesn't exist, create Clone_IDs based on the Clone field
  complete_data$Clone_IDs <- sapply(complete_data$Clone, function(x) {
    if (x == "Founder") return("1,2,3,4")
    # Extract numbers from clone name
    nums <- gsub("[^0-9,+]", "", x)
    nums <- gsub("\\+", ",", nums)
    return(nums)
  })
}

clonal_matrix <- data.frame(
  Mutation = complete_data$Mutation,
  VAF = complete_data$VAF,
  Type = complete_data$Type,
  Clone1 = grepl("1", complete_data$Clone_IDs),
  Clone2 = grepl("2", complete_data$Clone_IDs),
  Clone3 = grepl("3", complete_data$Clone_IDs),
  Clone4 = grepl("4", complete_data$Clone_IDs)
)

# Order by decreasing VAF
clonal_matrix <- clonal_matrix[order(-clonal_matrix$VAF), ]
clonal_matrix$Index <- 1:nrow(clonal_matrix)

# Reshape for ggplot
library(tidyr)
matrix_long <- pivot_longer(clonal_matrix,
                              cols = starts_with("Clone"),
                              names_to = "Clone",
                              values_to = "Present")

p5 <- ggplot(matrix_long, aes(x = Clone, y = Index, fill = Present)) +
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

ggsave("clonal_matrix.png", p5, width = 8, height = 10)

cat("\n===== GENERATED FILES =====\n")
cat("1. simulated_mutational_profile.csv - Complete mutation table\n")
cat("2. VAF_distribution_per_type.png - VAF histograms per mutation type\n")
cat("3. VAF_scatter_plot.png - Scatter plot of all VAFs\n")
cat("4. VAF_cumulative_density.png - DENSITY PLOT of mixed sample\n")
cat("5. VAF_violin_plot.png - Violin plot per mutation type\n")
cat("6. clonal_matrix.png - Heatmap of mutation presence in clones\n\n")

# Show first rows
cat("First rows of dataset:\n")
print(head(complete_data, 15))
