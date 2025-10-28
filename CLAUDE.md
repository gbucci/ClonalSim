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

**Note**: The `fishplot` package is automatically installed from GitHub if not present (lines 10-16). This requires `devtools` package, which will also be auto-installed if needed.

## Code Architecture

### Core Simulation Logic

The script follows a structured, modular pipeline:

1. **Parameter Configuration** (lines 21-64): Defines:
   - Subclone frequencies and mutation counts
   - Hierarchical evolutionary relationships
   - Biological noise parameters (Beta distribution for intra-tumor heterogeneity)
   - Technical sequencing noise parameters (depth variation, error rates, binomial sampling)

2. **Noise Model Functions** (lines 66-134): Modular functions implementing realistic noise:
   - `apply_biological_noise()`: Beta distribution for tumor heterogeneity (more realistic than Gaussian)
   - `simulate_depth()`: Negative binomial distribution for coverage variation (overdispersion)
   - `simulate_sequencing_reads()`: Binomial sampling for stochastic read counts + sequencing errors

3. **Mutation Generation Function** (lines 136-196): `generate_mutations()` creates mutations with realistic VAF:
   - Handles three mutation types: `founder`, `shared`, `private`
   - Step 1: Apply biological noise (intra-tumor heterogeneity)
   - Step 2: Simulate variable sequencing depth
   - Step 3: Apply technical noise (binomial sampling + base errors)
   - Returns both `True_VAF` (biological) and `VAF` (observed with sequencing noise)

4. **Data Generation** (lines 198-263): Creates mutations in hierarchical order:
   - Founder mutations (VAF = sum of all clone frequencies)
   - Shared mutations (VAF = sum of involved clone frequencies)
   - Private mutations (VAF = individual clone frequency)

5. **Output Generation** (lines 265-onwards): Produces CSV data file and 7 visualization plots, with detailed noise statistics

### Key Design Patterns

- **Hierarchical Structure**: The default evolutionary hierarchy is Clone1 → Clone2,3 → Clone4, where Clone2/3/4 share mutations that Clone1 doesn't have, and Clone3/4 share additional mutations
- **VAF Calculation**: VAF for mutation groups = sum of frequencies of clones carrying that mutation
- **Additive Model**: Clone frequencies are additive; their sum represents tumor purity (contamination from normal cells if sum < 1.0)
- **Two-Stage Noise Model**: Separates biological variation (tumor heterogeneity) from technical variation (sequencing)
- **Modular Functions**: Each noise component is independent and can be enabled/disabled

### Important Variables and Parameters

**Clone Configuration:**
- `subclone_freqs`: Vector of clone frequencies (default: 0.15, 0.25, 0.30, 0.30)
- `n_mut_shared`: Named list defining shared mutation structure and counts
- `Clone_IDs`: Comma-separated clone identifiers for tracking mutation presence

**Biological Noise (simulate_tumor_clones.R:43-49):**
- `biological_noise$enabled`: Enable/disable biological variation (default: TRUE)
- `biological_noise$concentration`: Beta distribution concentration parameter (default: 50)
  - Higher values = less intra-tumor heterogeneity
  - Lower values = more variable VAF within same clone
  - Typical range: 20-100

**Technical Sequencing Noise (simulate_tumor_clones.R:51-64):**
- `sequencing_noise$enabled`: Enable/disable sequencing noise (default: TRUE)
- `sequencing_noise$mean_depth`: Average sequencing coverage (default: 100)
- `sequencing_noise$depth_variation`: Distribution type (default: "negative_binomial")
  - Options: "negative_binomial" (realistic), "poisson" (simpler), "uniform" (testing)
- `sequencing_noise$depth_dispersion`: Negative binomial size parameter (default: 20)
  - Lower values = more variable coverage across positions
  - Higher values = more uniform coverage
- `sequencing_noise$error_rate`: Base miscall rate (default: 0.001 = 0.1%, typical for Illumina)
- `sequencing_noise$binomial_sampling`: Use binomial distribution for read counts (default: TRUE)

## Output Files

### Generated Artifacts

1. **simulated_mutational_profile.csv**: Complete dataset with columns:
   - Mutation, Chromosome, Position, Ref, Alt
   - **True_VAF**: Biological truth with intra-tumor heterogeneity (Beta distribution)
   - **VAF**: Observed VAF after sequencing noise (binomial sampling + errors)
   - Depth: Sequencing coverage (negative binomial distribution with realistic overdispersion)
   - Alt_reads: Alternative allele read count (from binomial sampling, not deterministic)
   - Clone: Which clone(s) carry the mutation
   - Type: Mutation type (founder/shared/private)

2. **VAF_cumulative_density.png**: Most important plot - shows VAF distribution of mixed sample (as real sequencing would show), with expected peaks at clone frequencies

3. **VAF_scatter_plot.png**: Scatter plot of all VAFs with reference lines

4. **VAF_distribution_per_type.png**: Histograms faceted by mutation type

5. **VAF_violin_plot.png**: Violin plots showing VAF distributions by type

6. **clonal_matrix.png**: Heatmap showing presence/absence of mutations in each clone, ordered by VAF

7. **fish_plot.png**: Fish plot visualization showing temporal evolution of clonal populations across simulated timepoints (diagnosis, treatment, progression). The plot shows:
   - **Horizontal axis**: Time progression in months (0, 6, 12, 18, 24)
   - **Vertical axis**: Clonal frequencies (proportion of tumor)
   - **Colored regions**: Different clones, with width representing their frequency at each timepoint
   - **Flow pattern**: How clones emerge, expand, and evolve over time
   - **Scaling**: Clone frequencies are scaled by 50x (scale_factor, line 333) for maximum visibility in the visualization

## Customization Patterns

### Simulating High Intra-Tumor Heterogeneity

```r
# More subclones with complex evolutionary structure
subclone_freqs <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.25)
n_mut_shared <- list(
  "2 3 4 5 6" = 15,
  "3 4 5 6" = 12,
  "4 5 6" = 10,
  "5 6" = 8
)

# Increase biological variability within clones
biological_noise$concentration <- 20  # Lower = more variable
```

### Simulating Low Purity Tumors

```r
# Low tumor content (40% tumor cells, 60% normal contamination)
subclone_freqs <- c(0.05, 0.10, 0.12, 0.13)  # Sum = 0.40
```

### Simulating Low-Coverage Sequencing

```r
# Exome sequencing typically has lower, more variable coverage
sequencing_noise$mean_depth <- 50
sequencing_noise$depth_dispersion <- 10  # More variable coverage
sequencing_noise$error_rate <- 0.002     # Slightly higher error rate
```

### Simulating High-Coverage, High-Quality Sequencing

```r
# Deep targeted sequencing with uniform coverage
sequencing_noise$mean_depth <- 500
sequencing_noise$depth_dispersion <- 100  # Very uniform coverage
sequencing_noise$error_rate <- 0.0005     # Lower error rate
```

### Disabling Noise for Testing

```r
# Turn off all noise to get ideal, deterministic VAFs
biological_noise$enabled <- FALSE
sequencing_noise$enabled <- FALSE
# Or keep sequencing but disable stochastic sampling:
sequencing_noise$binomial_sampling <- FALSE
```

## Use Cases

This simulator is designed for:

### Research Applications
- Testing clonal deconvolution algorithms (e.g., PyClone, SciClone)
- Benchmarking variant callers and mutation detection pipelines
- Development of tumor phylogenetic analysis methods
- Intratumoral heterogeneity studies
- Positive control for analysis pipelines
- Sensitivity testing of detection methods

### Education
- Teaching clonal evolution concepts
- Visualizing tumor heterogeneity
- Next-Generation Sequencing (NGS) data analysis exercises

## Biological Model Assumptions

- **Diploidy**: Assumes diploid genome (max theoretical VAF = 0.5 for heterozygous mutations)
- **Clonal Hierarchy**: Tree-like phylogenetic structure (no reticulate evolution)
- **No CNAs**: Does not simulate Copy Number Alterations (CNAs)
- **No LOH**: Does not include Loss of Heterozygosity (LOH)
- **Heterozygous Mutations**: All mutations assumed heterozygous
- **Sample Purity**: The sum of clonal frequencies represents tumor fraction (if < 1.0, remainder is normal cell contamination)

## Noise Model Details

### Biological Noise (Intra-Tumor Heterogeneity)

**Why Beta Distribution?**
- Bounded between 0 and 1 (appropriate for frequencies)
- More realistic than Gaussian for VAF data
- Can model both symmetric and skewed distributions
- Controlled by concentration parameter (alpha + beta)

**Implementation:** `apply_biological_noise()` (simulate_tumor_clones.R:70-90)
- For a clone with true frequency `f`, samples VAF from Beta(α, β)
- α = f × concentration, β = (1-f) × concentration
- Higher concentration → VAF closer to expected frequency
- Lower concentration → more spread (spatial/temporal heterogeneity)

### Technical Sequencing Noise

**1. Depth Variation (Coverage Overdispersion)**
- **Why Negative Binomial?** Real sequencing shows overdispersion (variance > mean)
- Poisson assumes variance = mean (unrealistic for NGS)
- Negative binomial allows extra variability due to:
  - GC content bias
  - Mappability differences
  - PCR amplification artifacts
- **Implementation:** `simulate_depth()` (simulate_tumor_clones.R:93-112)
- Dispersion parameter controls variability (lower = more variable)

**2. Binomial Sampling (Stochastic Read Counts)**
- **Why Binomial?** Each read is independently sampled from the DNA pool
- Given true VAF and depth D, alt reads ~ Binomial(D, VAF)
- This introduces natural sampling variation
- Low-frequency variants at low depth have high uncertainty
- **Implementation:** `simulate_sequencing_reads()` (simulate_tumor_clones.R:115-134)

**3. Sequencing Errors**
- Base miscall rate (default 0.1% for Illumina)
- Added to true VAF before binomial sampling
- Represents sequencer optical/chemistry errors

## Limitations and Extension Points

The README identifies future extensions:
- Copy Number Variation (CNV) support
- Loss of Heterozygosity (LOH) simulation
- Non-hierarchical (reticulate) evolution models
- Direct VCF format export
- Shiny interactive interface
- RNA-seq allelic expression data

## Troubleshooting

### Common Issues

**Error: "package 'ggplot2' is not available"**
```bash
# Install missing packages
R -e "install.packages(c('ggplot2', 'tidyr'), repos='https://cloud.r-project.org')"
```

**Error: "fishplot installation failed"**
- Ensure `devtools` is installed: `install.packages("devtools")`
- The script auto-installs fishplot from GitHub (chrisamiller/fishplot)
- Manual installation: `devtools::install_github("chrisamiller/fishplot")`

**Warning: "removed rows containing missing values"**
- This is normal - some mutations with very low/high VAFs are truncated to 0.01-0.99 range

**Frequencies don't sum to 1.0**
- This is intentional and represents normal cell contamination
- For 100% tumor purity, ensure `sum(subclone_freqs) = 1.0`

## Development Notes

- The script uses `set.seed(123)` (line 19) for reproducibility - modify or remove for different random realizations
- True VAF from Beta distribution is clamped to [0.01, 0.99] range (line 87)
- Sequencing depth uses negative binomial distribution (not Poisson) for realistic overdispersion (line 99)
- Minimum depth enforced at 10 reads (line 109)
- Binomial sampling ensures alt_reads ≤ depth always (line 124)
- Chromosome and position assignments are random and for demonstration purposes only (lines 253-256)
- All variable names, comments, and output files are in English for consistency
- Recent improvements:
  - Fish plot scaling iteratively improved (1x → 20x → 50x) for clone visibility
  - **Noise model completely refactored**: replaced Gaussian with Beta + negative binomial + binomial sampling
  - **Modular design**: noise components can be independently configured/disabled

## Git Workflow

This repository uses standard git workflow:
- Main branch: `main`
- The repository has a DOI citation (CITATION.cff) for academic use
- When making changes, follow conventional commit messages
- Generated output files (.png, .csv) are tracked in git for demonstration purposes
