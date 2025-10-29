# ClonalSim: Tumor Clonal Evolution Simulator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17465573.svg)](https://doi.org/10.5281/zenodo.17465573)

**ClonalSim** is a Bioconductor-compatible R package for simulating tumor clonal evolution with realistic sequencing noise. It generates mutational profiles of heterogeneous tumor samples with hierarchical clonal structure.

## Features

### Realistic Noise Modeling

ClonalSim implements a **two-stage noise model**:

1. **Biological Noise** (Intra-Tumor Heterogeneity)
   - Beta distribution for VAF variation
   - Models spatial/temporal heterogeneity within tumors
   - Configurable concentration parameter

2. **Technical Sequencing Noise**
   - **Depth overdispersion**: Negative binomial distribution (GC bias, mappability, PCR artifacts)
   - **Stochastic read sampling**: Binomial distribution (realistic sampling variation)
   - **Base calling errors**: Illumina-like error rates (0.1% default)

### Simulated Clonal Structure

- Hierarchical evolutionary relationships (e.g., Clone1 → Clone2,3 → Clone4)
- Three mutation types:
  - **Founder mutations**: present in all subclones
  - **Shared mutations**: present in clone subgroups
  - **Private mutations**: specific to individual subclones
- Configurable clone frequencies and mutation counts

### Generated Output

The package provides:

1. **ClonalSimData S4 object** with:
   - Mutation data (True_VAF and observed VAF)
   - Simulation parameters
   - Clonal structure information
   - Metadata (date, version, seed)

2. **Visualization functions**:
   - VAF density plots
   - VAF scatter plots by mutation type
   - Sequencing depth histograms
   - Clonal matrix heatmaps

3. **Export formats**:
   - GenomicRanges (GRanges objects)
   - VCF format
   - PyClone input format
   - SciClone input format
   - CSV/data.frame

## Installation

### From GitHub (Development Version)

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install ClonalSim
devtools::install_github("gbucci/ClonalSim")
```

### From Bioconductor (When Available)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ClonalSim")
```

### Dependencies

The package requires:
- R >= 4.3.0
- Bioconductor packages: GenomicRanges, IRanges, S4Vectors, VariantAnnotation
- CRAN packages: ggplot2, tidyr

## Quick Start

### Basic Simulation

```r
library(ClonalSim)

# Run simulation with default parameters
sim <- simulateTumor()

# View summary
sim

# Get detailed statistics
summary(sim)

# Visualize results
plot(sim, type = "vaf_density")
```

### Accessing Results

```r
# Get mutation data
mutations <- getMutations(sim)
head(mutations)

# Get true vs observed VAF
true_vaf <- getTrueVAF(sim)      # Biological truth
observed_vaf <- getObservedVAF(sim)  # With sequencing noise

# Get simulation parameters
params <- getSimParams(sim)
params$subclone_freqs

# Get clonal structure
structure <- getClonalStructure(sim)
```

### Visualization

```r
# VAF density plot (most important)
plot(sim, type = "vaf_density")

# VAF scatter plot by mutation type
plot(sim, type = "vaf_scatter")

# Sequencing depth distribution
plot(sim, type = "depth_histogram")

# Clonal matrix heatmap
plot(sim, type = "clone_matrix")
```

## Common Use Cases

### Low Purity Tumor (Normal Cell Contamination)

```r
sim_low_purity <- simulateTumor(
  subclone_freqs = c(0.05, 0.10, 0.15)  # Sum = 0.30 (30% tumor purity)
)
```

### High Coverage Sequencing

```r
sim_deep <- simulateTumor(
  sequencing_noise = list(
    mean_depth = 500,
    depth_dispersion = 100,  # More uniform coverage
    error_rate = 0.0005
  )
)
```

### Low Coverage Exome Sequencing

```r
sim_exome <- simulateTumor(
  sequencing_noise = list(
    mean_depth = 50,
    depth_dispersion = 10,  # More variable coverage
    error_rate = 0.002
  )
)
```

### High Intra-Tumor Heterogeneity

```r
sim_high_het <- simulateTumor(
  subclone_freqs = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.25),
  biological_noise = list(
    enabled = TRUE,
    concentration = 20  # Lower = more heterogeneity
  )
)
```

### Ideal Data (No Noise, for Testing)

```r
sim_ideal <- simulateTumor(
  biological_noise = list(enabled = FALSE),
  sequencing_noise = list(enabled = FALSE)
)
```

## Bioconductor Integration

### Export to GRanges

```r
library(GenomicRanges)

gr <- toGRanges(sim)
gr

# Subset by chromosome
gr_chr1 <- gr[seqnames(gr) == "chr1"]
```

### Export to VCF

```r
# Get VRanges object
vr <- toVCF(sim, sample_name = "TumorSample")

# Write to VCF file
toVCF(sim, output_file = "simulated_mutations.vcf")
```

### Export for Clonal Deconvolution Tools

```r
# PyClone format
toPyClone(sim, file = "pyclone_input.tsv", sample_id = "sample1")

# SciClone format
toSciClone(sim, file = "sciclone_input.tsv")

# Simple CSV
toDataFrame(sim, file = "mutations.csv")
```

## Custom Clonal Structures

```r
# Define custom evolutionary hierarchy
sim_custom <- simulateTumor(
  subclone_freqs = c(0.1, 0.15, 0.2, 0.25, 0.3),
  n_mut_per_clone = c(30, 40, 50, 40, 30),
  n_mut_shared = list(
    "2 3 4 5" = 20,  # Shared by clones 2,3,4,5
    "3 4 5" = 15,    # Shared by clones 3,4,5
    "4 5" = 10,      # Shared by clones 4,5
    "1 2" = 8        # Shared by clones 1,2
  ),
  n_mut_founder = 15
)
```

## Benchmarking Workflow Example

```r
# 1. Generate ground truth
sim <- simulateTumor(
  subclone_freqs = c(0.2, 0.3, 0.5),
  n_mut_per_clone = c(50, 75, 50),
  seed = 42
)

# 2. Get true mutations
true_mutations <- getMutations(sim)

# 3. Export for your variant caller
toVCF(sim, output_file = "ground_truth.vcf")

# 4. Run your variant caller on the VCF

# 5. Compare results and calculate metrics
# - Sensitivity: TP / (TP + FN)
# - Precision: TP / (TP + FP)
# - VAF correlation: cor(true_vaf, called_vaf)
```

## Documentation

### Vignettes

```r
# View introductory vignette
browseVignettes("ClonalSim")
```

### Function Help

```r
# Main function
?simulateTumor

# Noise model functions
?applyBiologicalNoise
?simulateDepth
?simulateSequencingReads

# Export functions
?toGRanges
?toVCF
?toPyClone

# Accessor functions
?getMutations
?getSimParams
```

## Use Cases

ClonalSim is designed for:

- **Benchmarking**: Test variant callers and mutation detection pipelines
- **Algorithm Development**: Develop and test clonal deconvolution algorithms
- **Education**: Teach tumor heterogeneity and clonal evolution concepts
- **Method Validation**: Positive controls for analysis pipelines
- **Research**: Study effects of sequencing parameters on variant detection

## Advanced: Standalone Script

For users who prefer the original standalone script (generates all plots automatically):

```r
# The original script is preserved in inst/scripts/
source(system.file("scripts", "simulate_tumor_clones.R", package = "ClonalSim"))
```

Or directly:

```bash
# From the package directory
Rscript inst/scripts/simulate_tumor_clones.R
```

## Output Data Structure

The `ClonalSimData` object contains:

```r
# Mutations data.frame with columns:
# - Mutation: Unique identifier
# - Chromosome, Position, Ref, Alt: Genomic coordinates
# - True_VAF: Biological truth (with heterogeneity)
# - VAF: Observed VAF (with sequencing noise)
# - Depth: Sequencing coverage
# - Alt_reads: Alternative allele read count
# - Clone: Which clone(s) carry the mutation
# - Type: founder, shared, or private
# - Clone_IDs: Comma-separated clone identifiers
```

## Biological Model

### Assumptions

- **Diploidy**: Assumes diploid genome (max theoretical VAF = 0.5 for heterozygous mutations)
- **Clonal Hierarchy**: Tree-like phylogenetic structure (no reticulate evolution)
- **Heterozygous Mutations**: All mutations assumed heterozygous
- **No CNAs/LOH**: Does not currently simulate Copy Number Alterations or Loss of Heterozygosity

### Noise Model Details

**Biological Noise (Beta Distribution):**
- For clone with frequency `f`, VAF ~ Beta(α, β)
- α = f × concentration, β = (1-f) × concentration
- Higher concentration → less variability

**Technical Noise:**
- Depth ~ NegativeBinomial(mean, dispersion)
- Alt_reads ~ Binomial(depth, VAF + error_rate)
- More realistic than Gaussian/Poisson models

## Citation

If you use ClonalSim in your research, please cite:

```
Bucci, G. (2025). ClonalSim: Simulation of Tumor Clonal Evolution with
Realistic Sequencing Noise. https://github.com/gbucci/ClonalSim
DOI: 10.5281/zenodo.17465573
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

Created for bioinformatics applications in computational oncology.

## Support

- **Issues**: https://github.com/gbucci/ClonalSim/issues
- **Documentation**: Run `browseVignettes("ClonalSim")` in R
- **Email**: gabriele.bucci@example.com

## References

1. McGranahan N, Swanton C. *Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future.* Cell. 2017
2. Dentro SC, Wedge DC, Van Loo P. *Principles of Reconstructing the Subclonal Architecture of Cancers.* Cold Spring Harb Perspect Med. 2017
3. Roth A, et al. *PyClone: statistical inference of clonal population structure in cancer.* Nat Methods. 2014
4. Miller CA, et al. *SciClone: inferring clonal architecture and tracking the spatial and temporal patterns of tumor evolution.* PLoS Comput Biol. 2014

## Related Tools

- **PyClone**: Statistical inference of clonal populations
- **SciClone**: Clonal architecture inference
- **ABSOLUTE**: Tumor purity and ploidy estimation
- **THetA**: Copy number and purity inference
- **PhylogicNDT**: Cancer phylogeny reconstruction

---

**Note**: ClonalSim generates synthetic data for educational, research, and benchmarking purposes.
