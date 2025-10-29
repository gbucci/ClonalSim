---
title: 'ClonalSim: An R Package for Configurable Simulation of Tumor Clonal Evolution and Sequencing Noise for Bioinformatics Benchmarking'
tags:
  - R
  - Bioconductor
  - computational biology
  - cancer genomics
  - somatic variant calling
  - simulation
authors:
  - name: Gabriele Bucci
    orcid: https://orcid.org/0000-0001-9838-7204
    affiliation: 1
affiliations:
  - name: IRCCS Ospedale San Raffaele, Milan, Italy
    index: 1
date: 28 October 2025 # Date of final version/Zenodo archive
# This is mandatory: replace with the DOI of the specific version submitted (Zenodo/Figshare)
bibliography: paper.bib 
---

## Summary

The accurate detection and characterization of somatic mutations, particularly those with low Variant Allele Frequency (VAF), is a fundamental challenge in cancer genomics. Algorithms designed to infer tumor subclonal architecture or perform high-sensitivity variant calling must be rigorously benchmarked against realistic datasets. However, existing _in silico_ simulation methods often employ simplistic models for both clonal structure and technical noise, leading to potentially over-optimistic performance evaluations \citep{schwarz2014phylogenetic}.

**ClonalSim** is an R/Bioconductor package designed to address this critical gap by generating synthetic mutational profiles that closely mimic real-world tumor sequencing data. It supports the creation of complex mutational landscapes that include:

1. **Hierarchical Clonal Structure:** Configuration of founder, shared, and private mutations according to user-defined evolutionary trees and cellular prevalence \citep{mcgranahan2017clonal}.
2. **High-Fidelity Noise Modeling:** A sophisticated, multi-staged stochastic model that accurately separates and simulates biological and technical sources of VAF variation.

ClonalSim serves as an essential tool for developing and validating next-generation sequencing (NGS) analysis pipelines, particularly for applications targeting minimal residual disease (MRD) or ultra-low frequency variants.

## Functionality and Architectural Details

ClonalSim’s core strength lies in its **modular and configurable noise architecture** \citep{claudemd_noise}. Unlike approaches that apply a single, uniform noise factor, ClonalSim uses a two-stage process that reflects the physical processes of tumor evolution and sequencing:

1. **Biological Heterogeneity (True VAF variation):** The true VAF for mutations within a given subclone is not fixed but sampled from a **Beta distribution**. This approach models the intrinsic biological dispersion (or "spread") of VAFs more realistically than Gaussian models, controlled by a user-defined *concentration* parameter \citep{newsmddetails}.

2. **Technical Noise (Sequencing and Sampling):**
    * **Sequencing Depth:** Read depth (`Depth`) for each mutation site is drawn from a **Negative Binomial distribution**. This is critical as it simulates the **overdispersion** (higher variance than the mean) commonly observed in real NGS data, which is poorly captured by simpler Poisson models \citep{newsmddetails}.
    * **Allele Counts:** The observed alternate read count (`Alt_reads`) is generated via **Binomial sampling** applied to the True VAF and the simulated Depth. This final step incorporates a configurable sequencing error rate (e.g., $0.1\%$ for Illumina platforms), effectively blurring the line between low-VAF variants and technical background noise \citep{claudemd_noise}.

The package leverages the R language's efficiency and integrates with the Bioconductor ecosystem by employing S4 classes for simulation results. Output formats include data frames, VCF, and GRanges objects, facilitating seamless integration with downstream tools like [PyClone](https://github.com/PyClone-Team/PyClone) and [SciClone](https://github.com/genome/sciclone) for benchmarking \citep{readmerelated}.

ClonalSim also includes a comprehensive visualization suite, featuring VAF density plots, clonal matrix heatmaps, and **Fish plots** (using the `fishplot` package) to dynamically visualize tumor evolution \citep{newsmddetails}. The package is supported by extensive unit tests with greater than 80% coverage to ensure stability and accuracy \citep{newsmddetails}.

## Acknowledgements

(Optional: Add any grants, funding, or specific lab acknowledgements here.)

## References

