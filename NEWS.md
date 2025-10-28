# ClonalSim 0.99.0

## New Features

* Initial Bioconductor submission
* Simulate tumor clonal evolution with hierarchical structure
* Realistic biological noise model using Beta distribution
* Technical sequencing noise with negative binomial depth variation
* Binomial sampling for stochastic read counts
* Sequencing error simulation (Illumina-like)
* S4 class for simulation results
* Export to GRanges and VCF formats
* Comprehensive visualization suite:
  - VAF distribution plots
  - Clonal matrix heatmap
  - Fish plot for temporal evolution
* Modular, configurable noise parameters
* Full Roxygen2 documentation
* Comprehensive vignettes for educational use
* Unit tests with >80% coverage

## Implementation Details

* Biological heterogeneity modeled with Beta distribution
* Depth overdispersion with negative binomial distribution
* Realistic read sampling with binomial distribution
* Configurable evolutionary hierarchies
* Support for founder, shared, and private mutations

## Future Plans

* Copy number alteration (CNA) simulation
* Loss of heterozygosity (LOH) support
* Subclonal CNA modeling
* Integration with mutational signature analysis tools
