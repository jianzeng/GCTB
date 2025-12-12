# GCTB
Genome-wide Complex Trait Bayesian analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-2.5.5-blue.svg)](https://github.com/jianzeng/GCTB)

### Overview

[GCTB](http://cnsgenomics.com/software/gctb) is a versatile command-line software suite for complex trait analysis using genome-wide SNP data. It implements a family of Bayesian mixture models that jointly fit all SNP effects and supports both individual-level genotype and phenotype data as well as GWAS summary statistics. The suite includes:

* **(S)BayesS** – for estimating key genetic architecture parameters such as SNP-based heritability, polygenicity, and the relationship between effect sizes and minor allele frequencies. 
  - Individual-level: [doi: 10.1038/s41588-018-0101-4](https://doi.org/10.1038/s41588-018-0101-4)
  - Summary-level: [doi: 10.1038/s41467-021-21446-3](https://doi.org/10.1038/s41467-021-21446-3)

* **(S)BayesR and (S)BayesRC** – for polygenic score prediction, with (S)BayesRC incorporating functional genomic annotations.
  - (S)BayesR: [doi: 10.1038/s41467-019-12653-0](https://doi.org/10.1038/s41467-019-12653-0)
  - (S)BayesRC: [doi: 10.1038/s41588-024-01704-y](https://doi.org/10.1038/s41588-024-01704-y)

* **GWFM** – for genome-wide fine-mapping of causal variants. [doi:10.1101/2024.07.18.24310667](https://doi.org/10.1101/2024.07.18.24310667)

### Requirements

* **Operating System**: Linux or macOS
* **Compiler**: C++ compiler with C++11 support (e.g., g++, clang++)
* **Dependencies**:
  - [OpenMP](https://www.openmp.org/) (for parallel computing)
  - [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) (header-only library, no compilation needed)
  - [Boost](http://www.boost.org/users/download) (header-only library, no compilation needed)

### Clone the Repository

Clone the GCTB repository using one of the following methods:

**Using HTTPS:**
```bash
git clone https://github.com/jianzeng/GCTB.git
cd GCTB
```

**Using SSH:**
```bash
git clone git@github.com:jianzeng/GCTB.git
cd GCTB
```

**Note**: If you want to work with a specific branch (e.g., `dev` for development version), you can switch to it after cloning:
```bash
git checkout dev
```

### Installation

#### Linux

1. Download [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](http://www.boost.org/users/download)
2. Edit the paths to Eigen3 and Boost in `scr/Makefile`:
   ```makefile
   EIGEN = /path/to/eigen-3.4.0
   BOOST = /path/to/boost_1_81_0
   ```
3. Load OpenMP library (consult your system administrator if needed)
4. Compile:
   ```bash
   cd scr
   make
   ```
5. The compiled binary `gctb` will be in the `scr/` directory

#### macOS

1. Install dependencies using Homebrew:
   ```bash
   brew install libomp
   brew install eigen
   brew install boost
   ```
2. Edit the paths in `scr/Makefile` to point to the Homebrew-installed libraries, or use the system paths
3. Compile:
   ```bash
   cd scr
   make
   ```

**Note**: If you encounter compilation errors, ensure that:
- OpenMP is properly installed and linked
- Eigen3 and Boost paths are correctly specified in the Makefile
- Your compiler supports C++11 standard

### Quick Start

#### Individual-level analysis (BayesS)

```bash
gctb --bfile data/genotype \
     --pheno data/phenotype.phen \
     --bayes S \
     --out results/bayess
```

#### Summary-level analysis (SBayesR)

```bash
gctb --gwas-summary data/summary.ma \
     --ldm-eigen data/ldm \
     --sbayes R \
     --out results/sbayesr
```

### Documentation

* **Tutorial**: See `tutorial/GCTB_tutorial.Rmd` for a comprehensive tutorial with examples
* **Online Documentation**: [https://gctbhub.cloud.edu.au/software/gctb/#Tutorial](https://gctbhub.cloud.edu.au/software/gctb/#Tutorial)
* **Website**: [http://cnsgenomics.com/software/gctb](http://cnsgenomics.com/software/gctb)

### Key Features

* Supports both individual-level and summary-level (GWAS) data
* Multiple LD matrix formats: full, block, sparse, and eigen-decomposed
* Efficient low-rank approximations for large-scale analyses
* Multi-chain MCMC for improved convergence assessment
* Annotation-stratified analysis (SBayesRC)
* Genome-wide fine-mapping (GWFM)

### Citation

If you use GCTB in your research, please cite the relevant publications:

* **BayesS**: Zeng et al. (2018) Nature Genetics. [doi: 10.1038/s41588-018-0101-4](https://doi.org/10.1038/s41588-018-0101-4)
* **SBayesS**: Zeng et al. (2021) Nature Communications. [doi: 10.1038/s41467-021-21446-3](https://doi.org/10.1038/s41467-021-21446-3)
* **SBayesR**: Lloyd-Jones et al. (2019) Nature Communications. [doi: 10.1038/s41467-019-12653-0](https://doi.org/10.1038/s41467-019-12653-0)
* **SBayesRC**: Zeng et al. (2024) Nature Genetics. [doi: 10.1038/s41588-024-01704-y](https://doi.org/10.1038/s41588-024-01704-y)
* **GWFM**: Zeng et al. (2024) bioRxiv. [doi:10.1101/2024.07.18.24310667](https://doi.org/10.1101/2024.07.18.24310667)

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Contact

For inquiries, please contact: **Jian Zeng** <j.zeng@uq.edu.au>

### Version

Current version: **2.5.5** (Last updated: December 12, 2025)
