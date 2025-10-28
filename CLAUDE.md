# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ClonalSim is an R-based tumor clonal evolution simulator. It generates realistic mutational profiles of tumor samples with hierarchical clonal structure, simulating:
- Founder mutations (present in all subclones)
- Shared mutations (present in subgroups according to evolutionary hierarchy)
- Private mutations (specific to individual subclones)
- Technical noise (sequencing errors)

The simulation outputs Variant Allele Frequency (VAF) data that mimics real tumor sequencing, useful for testing clonal deconvolution algorithms, benchmarking variant callers, and teaching tumor heterogeneity concepts.

## Key Commands

### Running the Simulator

```bash
# Make executable and run
chmod +x simulate_tumor_clones.R
./simulate_tumor_clones.R

# Or from R
Rscript simulate_tumor_clones.R
```

### Installing Dependencies

```bash
# Install required R packages
R -e "install.packages(c('ggplot2', 'tidyr'), repos='https://cloud.r-project.org')"

# Or from within R
install.packages(c("ggplot2", "tidyr"))
```

## Code Architecture

### Core Simulation Logic

The script follows a structured pipeline:

1. **Parameter Configuration** (lines 11-32): Defines subclone frequencies, mutation counts, hierarchical relationships, and technical noise parameters

2. **Mutation Generation Function** (lines 34-69): `generate_mutations()` creates mutations with specified VAF (Variant Allele Frequency), adding Gaussian noise to simulate biological/technical variation. It handles three mutation types:
   - `founder`: present in all clones
   - `shared`: present in clone subgroups
   - `private`: clone-specific

3. **Data Generation** (lines 71-120): Creates mutations in hierarchical order:
   - Founder mutations (VAF = sum of all clone frequencies)
   - Shared mutations (VAF = sum of involved clone frequencies)
   - Private mutations (VAF = individual clone frequency)

4. **Output Generation** (lines 137-311): Produces CSV data file and 5 visualization plots

### Key Design Patterns

- **Hierarchical Structure**: The default evolutionary hierarchy is Clone1 → Clone2,3 → Clone4, where Clone2/3/4 share mutations that Clone1 doesn't have, and Clone3/4 share additional mutations
- **VAF Calculation**: VAF for mutation groups = sum of frequencies of clones carrying that mutation
- **Additive Model**: Clone frequencies are additive; their sum represents tumor purity (contamination from normal cells if sum < 1.0)

### Important Variables

- `subclone_freqs`: Vector of clone frequencies (default: 0.15, 0.25, 0.30, 0.30)
- `n_mut_shared`: Named list defining shared mutation structure and counts
- `noise_sd`: Standard deviation for Gaussian noise (default: 0.02)
- `Clone_IDs`: Comma-separated clone identifiers for tracking mutation presence

## Output Files

### Generated Artifacts

1. **simulated_mutational_profile.csv**: Complete dataset with columns:
   - Mutation, Chromosome, Position, Ref, Alt
   - VAF (Variant Allele Frequency)
   - Depth (sequencing coverage), Alt_reads
   - Clone (which clone(s) carry it), Type (mutation type)

2. **VAF_cumulative_density.png**: Most important plot - shows VAF distribution of mixed sample (as real sequencing would show), with expected peaks at clone frequencies

3. **VAF_scatter_plot.png**: Scatter plot of all VAFs with reference lines

4. **VAF_distribution_per_type.png**: Histograms faceted by mutation type

5. **VAF_violin_plot.png**: Violin plots showing VAF distributions by type

6. **clonal_matrix.png**: Heatmap showing presence/absence of mutations in each clone, ordered by VAF

7. **fish_plot.png**: Fish plot visualization showing temporal evolution of clonal populations across simulated timepoints (diagnosis, treatment, progression)

## Customization Patterns

### Simulating High Heterogeneity

Increase number of subclones and create more complex hierarchical structure:
```r
subclone_freqs <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.25)
n_mut_shared <- list(
  "2 3 4 5 6" = 15,
  "3 4 5 6" = 12,
  "4 5 6" = 10,
  "5 6" = 8
)
```

### Simulating Low Purity Tumors

Reduce clone frequencies so their sum is << 1.0:
```r
subclone_freqs <- c(0.05, 0.10, 0.12, 0.13)  # Sum = 0.40 (40% tumor cells)
noise_sd <- 0.03  # Increase noise
```

## Biological Model Assumptions

- **Diploidy**: Assumes diploid genome (max theoretical VAF = 0.5 for heterozygous mutations)
- **Clonal Hierarchy**: Tree-like phylogenetic structure (no reticulate evolution)
- **No CNAs**: Does not simulate Copy Number Alterations
- **No LOH**: Does not include Loss of Heterozygosity
- **Heterozygous Mutations**: All mutations assumed heterozygous
- **Gaussian Noise**: Technical errors modeled as Gaussian with configurable standard deviation

## Limitations and Extension Points

The README identifies future extensions:
- Copy Number Variation (CNV) support
- Loss of Heterozygosity (LOH) simulation
- Non-hierarchical (reticulate) evolution models
- Direct VCF format export
- Shiny interactive interface
- RNA-seq allelic expression data

## Development Notes

- The script uses `set.seed(123)` for reproducibility - modify or remove for different random realizations
- VAF values are clamped to [0.01, 0.99] range (line 58)
- Sequencing depth is simulated from Poisson distribution with lambda=100 (line 129)
- Chromosome and position assignments are random and for demonstration purposes only
- All variable names, comments, and output files are in English for consistency
