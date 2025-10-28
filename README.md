# Tumor Mutational Profile Simulator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17465573.svg)](https://doi.org/10.5281/zenodo.17465573)

R script for simulating mutational profiles of tumor samples with hierarchical clonal structure.

## Description

This tool allows you to simulate a tumor DNA sample composed of a mixture of clones and subclones with specific Variant Allele Frequencies (VAF - Variant Allele Frequency). The script generates realistic data that includes:

- **Founder mutations**: present in all subclones (initial events)
- **Shared mutations**: present in subgroups of clones according to evolutionary hierarchy
- **Private mutations**: specific to individual subclones
- **Technical noise**: sequencing error simulation

## Features

### Simulated Clonal Structure

The script simulates a heterogeneous tumor with:
- 4 subclones with configurable frequencies (default: 0.15, 0.25, 0.30, 0.30)
- Evolutionary hierarchy: Clone1 → Clone2,3 → Clone4
- Mutations distributed according to realistic evolutionary model

### Generated Output

1. **simulated_mutational_profile.csv**: Complete dataset with all mutations
2. **VAF_cumulative_density.png**: Density plot of mixed sample (as in real sequencing)
3. **VAF_scatter_plot.png**: Scatter plot of all VAFs
4. **VAF_distribution_per_type.png**: Histograms per mutation type
5. **VAF_violin_plot.png**: Violin plots to visualize distributions
6. **clonal_matrix.png**: Heatmap of mutation presence/absence in clones

## Requirements

### Software
- R >= 4.0
- R packages:
  - `ggplot2`
  - `tidyr`

### Installing Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install r-base r-cran-ggplot2 r-cran-tidyr

# Or from R
install.packages(c("ggplot2", "tidyr"))
```

## Usage

### Basic Execution

```bash
# Make script executable
chmod +x simulate_tumor_clones.R

# Run
./simulate_tumor_clones.R
```

Or from R:

```r
source("simulate_tumor_clones.R")
```

### Parameter Customization

Modify the parameters at the beginning of the script:

```r
# Frequencies of the 4 subclones (must sum to <= 1)
subclone_freqs <- c(0.15, 0.25, 0.30, 0.30)

# Number of mutations per subclone (private)
n_mut_per_clone <- c(20, 25, 30, 15)

# Number of founder mutations
n_mut_founder <- 10

# Shared mutation structure (evolutionary hierarchy)
n_mut_shared <- list(
  "2 3 4" = 15,    # Shared by Clone2, Clone3, Clone4
  "3 4" = 12,      # Shared by Clone3 and Clone4
  "1 2" = 8        # Shared by Clone1 and Clone2
)

# Technical noise (standard deviation)
noise_sd <- 0.02
```

## Output Data Structure

### CSV File

The `simulated_mutational_profile.csv` file contains the following columns:

| Column | Description |
|---------|-------------|
| `Mutation` | Unique mutation identifier |
| `Chromosome` | Chromosome (chr1-chr22) |
| `Position` | Genomic position |
| `Ref` | Reference allele |
| `Alt` | Alternative allele |
| `VAF` | Variant Allele Frequency (0-1) |
| `Depth` | Sequencing coverage |
| `Alt_reads` | Number of reads with alternative allele |
| `Clone` | Clone(s) carrying the mutation |
| `Type` | Mutation type (founder/shared/private) |
| `Clone_IDs` | Numerical IDs of involved clones |

### Example Output

```
Mutation            Chromosome  Position  Ref  Alt  VAF     Depth  Alt_reads  Clone       Type
Founder_1           chr19       26274770  C    T    0.9888  101    100        Founder     founder
Shared_C2_3_4_mut1  chr9        76062670  C    C    0.8745  85     74         Clone2+3+4  shared
Clone1_mut1         chr15       45123456  A    G    0.1498  95     14         Clone1      private
```

## Interpreting Results

### Cumulative Density Plot

The **VAF_cumulative_density.png** plot represents the distribution of allelic frequencies as if you had sequenced the DNA from the real tumor (mixture of all clones).

**Expected peaks:**
- **~1.0**: Founder mutations (present in all tumor cells)
- **~0.85**: Mutations shared by Clone2+3+4
- **~0.60**: Mutations shared by Clone3+4
- **~0.40**: Mutations shared by Clone1+2
- **0.15-0.30**: Private mutations of individual clones

### Clonal Matrix

The **clonal_matrix.png** heatmap shows which mutations are present in which clones, ordered by decreasing VAF. This clearly visualizes the evolutionary hierarchy of the tumor.

## Biological Model

### Model Assumptions

1. **Sample purity**: The sum of clonal frequencies represents the tumor fraction
2. **Ploidy**: Diploid assumption (maximum theoretical VAF = 0.5 for heterozygous mutations)
3. **Clonal evolution**: Hierarchical structure (phylogenetic tree)
4. **Technical errors**: Gaussian noise with configurable SD

### Limitations

- Does not simulate Copy Number Alterations (CNA)
- Does not consider explicit contamination with normal cells
- Assumes heterozygous mutations
- Does not include subclonal Loss of Heterozygosity (LOH)

## Use Cases

### Research

- Testing clonal deconvolution algorithms
- Benchmarking variant callers
- Development of tumor phylogenetic analysis methods
- Intratumoral heterogeneity studies

### Education

- Teaching clonal evolution concepts
- Visualizing tumor heterogeneity
- NGS data analysis exercises

### Validation

- Positive control for analysis pipelines
- Sensitivity testing of detection methods
- Validation of bioinformatics tools

## Advanced Examples

### Simulating a Tumor with High Heterogeneity

```r
# 6 subclones with different frequencies
subclone_freqs <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.25)

# Increase number of mutations
n_mut_per_clone <- c(30, 40, 50, 60, 70, 50)

# More complex hierarchical structure
n_mut_shared <- list(
  "1 2 3 4 5 6" = 20,  # Founder (alternative to n_mut_founder)
  "2 3 4 5 6" = 15,
  "3 4 5 6" = 12,
  "4 5 6" = 10,
  "5 6" = 8,
  "1 2" = 5
)
```

### Simulating Low Purity Tumor

```r
# Low frequencies to simulate normal contamination
subclone_freqs <- c(0.05, 0.10, 0.12, 0.13)  # Sum = 0.40 (40% tumor cells)

# Increase noise to simulate greater errors
noise_sd <- 0.03
```

### Export for Downstream Analysis

```r
# After generating the data
data <- read.csv("simulated_mutational_profile.csv")

# Convert to VCF-like format
vcf_like <- data[, c("Chromosome", "Position", "Ref", "Alt", "VAF", "Depth", "Alt_reads")]

# Save
write.table(vcf_like, "mutations.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

## Troubleshooting

### Error: "package 'ggplot2' is not available"

```bash
# Install missing packages
R -e "install.packages(c('ggplot2', 'tidyr'), repos='https://cloud.r-project.org')"
```

### Warning: "removed rows containing missing values"

Normal if some mutations have very low or high VAFs that are truncated to 0.01-0.99.

### Frequencies don't sum to 1

This is intentional! It represents the possibility of contamination with normal cells. If you want 100% purity, ensure that `sum(subclone_freqs) = 1.0`.

## Future Extensions

- [ ] Support for Copy Number Variations (CNV)
- [ ] Loss of Heterozygosity (LOH) simulation
- [ ] Non-hierarchical (reticulate) evolution models
- [ ] Direct export to VCF format
- [ ] Shiny interface for interactive use
- [ ] RNA-seq data simulation with allelic expression

## Contributing

Contributions, bug reports, and feature requests are welcome!

### How to Contribute

1. Fork the repository
2. Create a branch for your feature (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is released under the MIT license. See the `LICENSE` file for details.

## Authors

Created for bioinformatics applications in computational oncology.

## Citation

If you use this tool in your research, please cite:

```
Tumor Clonal Simulator - Tool for simulating tumor mutational profiles
with hierarchical clonal structure
```

## Contact

For questions, suggestions, or collaborations, open an issue on GitHub.

## References

### Recommended Reading

1. McGranahan N, Swanton C. **Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future.** Cell. 2017
2. Dentro SC, Wedge DC, Van Loo P. **Principles of Reconstructing the Subclonal Architecture of Cancers.** Cold Spring Harb Perspect Med. 2017
3. Schwarz RF, et al. **Phylogenetic Quantification of Intra-tumour Heterogeneity.** PLoS Comput Biol. 2014

### Related Tools

- **PyClone**: Statistical analysis of tumor subclones
- **SciClone**: Identification of mutation clusters
- **ABSOLUTE**: Estimation of tumor purity and ploidy
- **THetA**: Inference of copy number and purity

---

**Note**: This is a simulation tool. The generated data is synthetic and for educational/research purposes.
