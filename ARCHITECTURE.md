# GCTB Python/C++ Hybrid Architecture

## Overview: Separation of Concerns

The design philosophy is **Python for interface, C++ for computation**.

```
┌─────────────────────────────────────────────────────────────┐
│                       USER INTERFACE                        │
│                      (Python Layer)                         │
│  - CLI parsing (click)                                      │
│  - Input validation                                         │
│  - File path handling                                       │
│  - Progress bars (tqdm)                                     │
│  - Results visualization (matplotlib)                       │
│  - Pandas DataFrame integration                             │
│  - High-level workflow orchestration                        │
└─────────────────────────────────────────────────────────────┘
                            ↕
              [pybind11 bindings layer]
                            ↕
┌─────────────────────────────────────────────────────────────┐
│                   COMPUTATIONAL CORE                        │
│                      (C++ Layer)                            │
│  - Data structures (SnpInfo, IndInfo, Data)                │
│  - Matrix operations (Eigen)                                │
│  - MCMC sampling engine                                     │
│  - Bayesian model implementations                           │
│  - LD matrix computations                                   │
│  - Statistical calculations                                 │
│  - Parallel processing (OpenMP)                             │
└─────────────────────────────────────────────────────────────┐
```

---

## Detailed Architecture

### 1. Python Layer Components

#### A. CLI Interface (`python/gctb/cli.py`)
**Purpose:** User-facing command-line interface

```python
# Current: MISSING
# Target structure:

import click
from . import _core

@click.group()
def main():
    """GCTB: Genome-wide Complex Trait Bayesian Analysis"""
    pass

@main.command()
@click.option('--bfile', required=True, help='PLINK binary file prefix')
@click.option('--pheno', required=True, help='Phenotype file')
@click.option('--bayes', type=click.Choice(['C', 'R', 'S']), help='Bayes type')
@click.option('--chain-length', default=3000, help='MCMC chain length')
@click.option('--burnin', default=1000, help='Burnin iterations')
@click.option('--out', default='gctb', help='Output prefix')
def bayes(bfile, pheno, bayes, chain_length, burnin, out):
    """Run Bayes analysis on individual-level data"""
    from .workflows import run_bayes_analysis
    run_bayes_analysis(bfile, pheno, bayes, chain_length, burnin, out)

@main.command()
@click.option('--ldm', required=True, help='LD matrix directory')
@click.option('--gwas-summary', required=True, help='GWAS summary file')
@click.option('--sbayes', type=click.Choice(['R', 'S', 'RC']), help='SBayes type')
def sbayes(ldm, gwas_summary, sbayes, **kwargs):
    """Run SBayes analysis on summary statistics"""
    from .workflows import run_sbayes_analysis
    run_sbayes_analysis(ldm, gwas_summary, sbayes, **kwargs)

# etc...
```

**Replaces:** All the string parsing in `options.cpp` (832 lines!)

---

#### B. Configuration Management (`python/gctb/config.py`)
**Purpose:** Type-safe configuration objects

```python
# Current: MISSING
# Target:

from dataclasses import dataclass
from typing import Optional
from pathlib import Path

@dataclass
class BayesConfig:
    """Configuration for Bayes analysis"""
    bed_file: Path
    phenotype_file: Path
    bayes_type: str = 'C'
    chain_length: int = 3000
    burnin: int = 1000
    thin: int = 10
    heritability: float = 0.5
    pi: float = 0.01
    output_prefix: str = 'gctb'
    
    def validate(self):
        """Validate configuration"""
        if not self.bed_file.exists():
            raise FileNotFoundError(f"BED file not found: {self.bed_file}")
        # ... more validation

@dataclass
class SBayesConfig:
    """Configuration for SBayes analysis"""
    ldm_dir: Path
    gwas_summary_file: Path
    sbayes_type: str = 'R'
    # ... etc
```

**Replaces:** Options class initialization in C++

---

#### C. Workflow Orchestration (`python/gctb/workflows.py`)
**Purpose:** High-level analysis workflows

```python
# Current: MISSING
# Target:

from pathlib import Path
import pandas as pd
from tqdm import tqdm
from . import _core
from .config import BayesConfig

def run_bayes_analysis(bfile: str, pheno: str, bayes_type: str, 
                       chain_length: int, burnin: int, output: str):
    """
    High-level wrapper for Bayes analysis.
    
    This is pure Python orchestration calling C++ core.
    """
    
    # 1. Validate inputs (Python)
    config = BayesConfig(
        bed_file=Path(bfile + '.bed'),
        phenotype_file=Path(pheno),
        bayes_type=bayes_type,
        chain_length=chain_length,
        burnin=burnin,
        output_prefix=output
    )
    config.validate()
    
    # 2. Load data (C++ core)
    print(f"Loading data from {bfile}...")
    data = _core.Data()
    data.read_fam_file(str(bfile + '.fam'))
    data.read_bim_file(str(bfile + '.bim'))
    data.read_bed_file(False, str(bfile + '.bed'))
    data.read_phenotype_file(str(pheno), 1)
    
    print(f"Loaded {data.num_snps} SNPs, {data.num_inds} individuals")
    
    # 3. Build model (C++ core)
    print(f"Building Bayes{bayes_type} model...")
    model = _core.build_model(
        data=data,
        bayes_type=bayes_type,
        heritability=config.heritability,
        pi=config.pi
    )
    
    # 4. Run MCMC (C++ core with Python progress bar)
    print("Running MCMC...")
    with tqdm(total=chain_length, desc="MCMC iterations") as pbar:
        results = _core.run_mcmc(
            model=model,
            chain_length=chain_length,
            burnin=burnin,
            callback=lambda i: pbar.update(1)  # Progress callback
        )
    
    # 5. Process results (Python + pandas)
    print("Processing results...")
    results_df = results_to_dataframe(results, data)
    
    # 6. Save results (Python)
    results_df.to_csv(f"{output}.snpRes", sep='\t', index=False)
    print(f"Results saved to {output}.snpRes")
    
    # 7. Visualize (Python)
    plot_results(results_df, output)
    
    return results_df
```

**Replaces:** The if-else logic in `main.cpp` (lines 52-476)

---

#### D. Data Utilities (`python/gctb/data_utils.py`)
**Purpose:** Helper functions for data manipulation

```python
# Current: MISSING
# Target:

import pandas as pd
import numpy as np
from . import _core

def results_to_dataframe(mcmc_results, data) -> pd.DataFrame:
    """Convert C++ MCMC results to pandas DataFrame"""
    
    snps = data.get_incd_snp_info_vec()
    
    df = pd.DataFrame({
        'SNP': [snp.ID for snp in snps],
        'CHR': [snp.chrom for snp in snps],
        'POS': [snp.physPos for snp in snps],
        'A1': [snp.a1 for snp in snps],
        'A2': [snp.a2 for snp in snps],
        'BETA': mcmc_results.posterior_mean,
        'PIP': mcmc_results.pip,
        'VAR_EXPLAINED': mcmc_results.posterior_mean ** 2  # simplified
    })
    
    return df

def load_gwas_summary(filepath: str) -> pd.DataFrame:
    """Load and validate GWAS summary statistics"""
    df = pd.read_csv(filepath, sep='\t')
    
    # Validate required columns
    required = ['SNP', 'A1', 'A2', 'BETA', 'SE', 'N']
    if not all(col in df.columns for col in required):
        raise ValueError(f"Missing required columns. Need: {required}")
    
    return df
```

---

#### E. Visualization (`python/gctb/visualization.py`)
**Purpose:** Result visualization

```python
# Current: MISSING
# Target:

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_manhattan(results_df: pd.DataFrame, output: str):
    """Create Manhattan plot of posterior inclusion probabilities"""
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Color by chromosome
    colors = sns.color_palette("husl", results_df['CHR'].nunique())
    
    for i, (chrom, group) in enumerate(results_df.groupby('CHR')):
        ax.scatter(group.index, -np.log10(1 - group['PIP']), 
                   c=[colors[i]], s=10, alpha=0.6, label=f'Chr{chrom}')
    
    ax.set_xlabel('SNP Index')
    ax.set_ylabel('-log10(1-PIP)')
    ax.set_title('Posterior Inclusion Probability')
    plt.tight_layout()
    plt.savefig(f"{output}_manhattan.png", dpi=300)
    plt.close()

def plot_mcmc_trace(mcmc_samples, output: str):
    """Plot MCMC trace for variance components"""
    # ... implementation
```

---

### 2. C++ Core Components (What Stays in C++)

#### A. Core Classes (Keep as-is, just bind them)

```cpp
// All in scr/*.cpp - KEEP UNCHANGED
// Just expose via pybind11

class SnpInfo { ... };      // SNP metadata
class IndInfo { ... };      // Individual metadata  
class Data { ... };         // Main data container
class Model { ... };        // Bayesian models
class MCMC { ... };         // MCMC engine
class McmcSamples { ... };  // MCMC results storage
```

**These do NOT change** - we just create Python bindings.

---

#### B. Computation Kernels (Pure C++, stay in C++)

```cpp
// Matrix operations - STAY IN C++
void Data::makeLDmatrix(...) { ... }
void Data::buildSparseMME(...) { ... }

// MCMC sampling - STAY IN C++
void MCMC::run(...) { ... }
void Model::sampleUnknowns() { ... }

// Statistical calculations - STAY IN C++
float Stat::ranf() { ... }
float Stat::snorm() { ... }
```

**No Python versions needed** - these are called from Python via bindings.

---

### 3. Binding Layer (`python/bindings/bindings.cpp`)

**Purpose:** Connect Python to C++

```cpp
// Current status: PARTIAL
// What exists:
✅ SnpInfo - fully bound
✅ IndInfo - fully bound
✅ Data - partially bound (read methods only)
✅ Timer - fully bound

// What's MISSING:
❌ Model class
❌ MCMC class  
❌ McmcSamples class
❌ Options class (may not need - replace with Python)
❌ Many Data methods (60+ methods not bound yet)
❌ Callback support (for progress bars)
❌ Exception translation (for better error messages)
```

---

## Current State Analysis

### ✅ What's Complete

#### 1. Infrastructure (100%)
- [x] Build system (CMakeLists.txt)
- [x] Python packaging (setup.py, pyproject.toml)
- [x] Virtual environment
- [x] Test framework
- [x] Git branch structure

#### 2. Basic Bindings (30%)
- [x] SnpInfo class (100%)
- [x] IndInfo class (100%)
- [x] Data class (20% - only 12 of 60+ methods)
- [x] Timer utility (100%)
- [ ] Model class (0%)
- [ ] MCMC class (0%)
- [ ] McmcSamples class (0%)

#### 3. Data I/O (60%)
- [x] FAM file reading ✅
- [x] BIM file reading ✅
- [~] BED file reading (works but throws exception)
- [ ] Phenotype file reading
- [ ] Covariate file reading
- [ ] GWAS summary reading
- [ ] LD matrix reading

#### 4. Python Interface (0%)
- [ ] CLI (cli.py)
- [ ] Config classes (config.py)
- [ ] Workflows (workflows.py)
- [ ] Utilities (data_utils.py)
- [ ] Visualization (visualization.py)

---

## What's Missing: Detailed Breakdown

### Missing Category 1: C++ Bindings

#### A. Data Class - 48 Unbounded Methods

```cpp
// File I/O methods - NEED BINDING
void readCovariateFile(string)
void readRandomCovariateFile(string)
void readGwasSummaryFile(...)
void readLDmatrixInfoFile(string)
void readLDmatrixBinFile(string)
void readGeneticMapFile(string)

// SNP filtering - NEED BINDING
void includeSnp(string)
void excludeSnp(string)
void includeChr(unsigned)
void excludeAmbiguousSNP()
void excludeSNPwithMaf(float, float)

// Matrix operations - NEED BINDING
void makeLDmatrix(...)
void buildSparseMME(...)
void getZPZmat()

// Results output - NEED BINDING
void outputSnpResults(...)
void outputFixedEffects(...)
```

**Strategy:** Bind incrementally as needed for workflows.

---

#### B. Model Class - Completely Unbound

```cpp
// In model.hpp - 0% bound

class Model {
public:
    // Core methods
    void sampleUnknowns();
    void computeResiduals();
    
    // Need to expose:
    Data *data;
    VectorXf snpEffects;
    float varGenotypic;
    float varResidual;
    // ... many more
};

// Multiple model types
class BayesC : public Model { ... };
class BayesR : public Model { ... };
class BayesS : public Model { ... };
// ... etc
```

**Challenge:** Model is a base class with many derived classes.

**Strategy:** 
1. Bind base Model class first
2. Use factory function: `Model* build_model(data, type, params)`
3. Avoid exposing inheritance to Python

---

#### C. MCMC Class - Completely Unbound

```cpp
// In mcmc.hpp - 0% bound

class McmcSamples {
public:
    MatrixXf datMat;         // Dense storage
    SpMat datMatSp;          // Sparse storage
    VectorXf posteriorMean;
    VectorXf pip;
    // Need all these exposed
};

class MCMC {
public:
    vector<McmcSamples*> run(Model &model, 
                              unsigned chainLength,
                              unsigned burnin, 
                              unsigned thin,
                              bool print,
                              unsigned outputFreq,
                              string title,
                              bool writeBinPosterior,
                              bool writeTxtPosterior);
};
```

**Strategy:** 
1. Bind McmcSamples as simple container
2. Add conversion to Python dict/pandas
3. Add progress callback support

---

### Missing Category 2: Python Layer (100% Missing)

#### A. CLI Interface - 0% Complete
- Need: ~200 lines of click decorators
- Replaces: 832 lines of C++ string parsing
- Time: ~4 hours

#### B. Workflow Functions - 0% Complete  
- Need: ~300 lines for main workflows
- Replaces: ~500 lines of main.cpp logic
- Time: ~6 hours

#### C. Utilities - 0% Complete
- Need: ~200 lines of data conversion
- New capability (not in C++ version)
- Time: ~3 hours

#### D. Visualization - 0% Complete
- Need: ~150 lines of plotting code
- New capability
- Time: ~3 hours

---

### Missing Category 3: Error Handling

#### Current Issues

```python
# This throws generic exception:
data.read_bed_file(False, "test.bed")
# RuntimeError: Caught an unknown exception!

# Should throw:
# FileFormatError: BED file magic number mismatch
# Or: IOError: Cannot open test.bed
```

**What's needed:**

```cpp
// In bindings.cpp - ADD THIS

PYBIND11_MODULE(_core, m) {
    // Register exception translations
    py::register_exception<FileFormatError>(m, "FileFormatError");
    py::register_exception<DataMismatchError>(m, "DataMismatchError");
    
    // Wrap functions with try-catch
    m.def("read_bed_wrapper", [](Data& data, bool noscale, string file) {
        try {
            data.readBedFile(noscale, file);
        } catch (const std::exception& e) {
            throw py::value_error(e.what());
        } catch (...) {
            throw py::runtime_error("Unknown error reading BED file");
        }
    });
}
```

---

## Connection Pattern: Python ↔ C++

### Pattern 1: Direct Method Call

```python
# Python code
data = gctb.Data()
data.read_fam_file("test.fam")  # Direct call to C++

# C++ (bindings.cpp)
py::class_<Data>(m, "Data")
    .def("read_fam_file", &Data::readFamFile);

# C++ (data.cpp)
void Data::readFamFile(const string &famFile) {
    // Implementation
}
```

**Data flow:** Python → pybind11 → C++ method → return to Python

---

### Pattern 2: Factory Pattern (for complex objects)

```python
# Python code
model = gctb.build_model(data, bayes_type='R', heritability=0.5)

# C++ (bindings.cpp)
m.def("build_model", [](Data& data, string type, float h2) {
    GCTB gctb_engine(/* config */);
    return gctb_engine.buildModel(data, "", "", type, 0, h2, ...);
});

# C++ (gctb.cpp)
Model* GCTB::buildModel(...) {
    if (type == "R") return new BayesR(...);
    if (type == "C") return new BayesC(...);
    // ...
}
```

**Data flow:** Python → wrapper function → factory → return pointer → Python

---

### Pattern 3: Callback Pattern (for progress)

```python
# Python code
def progress_callback(iteration):
    pbar.update(1)

results = gctb.run_mcmc(model, callback=progress_callback)

# C++ (bindings.cpp)
m.def("run_mcmc", [](Model& model, py::function callback) {
    MCMC mcmc;
    for (int i = 0; i < iterations; i++) {
        // Do MCMC step
        model.sampleUnknowns();
        
        // Call Python callback
        if (i % 10 == 0) {
            py::gil_scoped_acquire acquire;  // Get GIL
            callback(i);
        }
    }
});
```

**Data flow:** C++ → acquire GIL → call Python → release GIL → continue C++

---

### Pattern 4: Data Conversion (Eigen ↔ NumPy)

```python
# Python code
snp_effects = results.posterior_mean  # Gets numpy array

# C++ (bindings.cpp) - pybind11 does this automatically!
py::class_<McmcSamples>(m, "McmcSamples")
    .def_readonly("posterior_mean", &McmcSamples::posteriorMean);
    // pybind11 auto-converts Eigen::VectorXf → numpy.ndarray

# C++ (mcmc.hpp)
VectorXf posteriorMean;  // Eigen vector
```

**Data flow:** Eigen vector → pybind11 auto-convert → NumPy array (zero-copy!)

---

## Detailed Implementation Plan

### Phase 1: Fix Current Issues (2 hours)

#### Task 1.1: Fix BED file exception (30 min)
```cpp
// Add proper exception handling in bindings
```

#### Task 1.2: Add exception translation (30 min)
```cpp
// Register custom exceptions
```

#### Task 1.3: Test all data I/O (1 hour)
```python
# Comprehensive I/O tests
```

---

### Phase 2: Expand Core Bindings (2 days = 16 hours)

#### Task 2.1: Bind Model class (4 hours)
- Base Model class
- Factory function
- Key properties (varGenotypic, varResidual, etc.)
- Test model creation

#### Task 2.2: Bind MCMC class (4 hours)  
- MCMC::run() method
- McmcSamples class
- Conversion to Python dict/DataFrame
- Test simple MCMC run

#### Task 2.3: Expand Data bindings (4 hours)
- Add remaining file I/O methods
- Add SNP filtering methods
- Add result output methods
- Test each method

#### Task 2.4: Add callback support (2 hours)
- Progress callback for MCMC
- Test with tqdm

#### Task 2.5: Complete testing (2 hours)
- Write tests for all bindings
- Fix any issues found

---

### Phase 3: Python Layer (2 days = 16 hours)

#### Task 3.1: CLI interface (4 hours)
```python
# gctb/cli.py - ~200 lines
# - bayes command
# - sbayes command  
# - ldmatrix command
# - predict command
```

#### Task 3.2: Workflow functions (6 hours)
```python
# gctb/workflows.py - ~300 lines
# - run_bayes_analysis()
# - run_sbayes_analysis()
# - make_ld_matrix()
# - predict()
```

#### Task 3.3: Utilities & visualization (4 hours)
```python
# gctb/data_utils.py - ~150 lines
# gctb/visualization.py - ~150 lines
```

#### Task 3.4: Integration testing (2 hours)
- Test complete workflows end-to-end
- Compare outputs with C++ version

---

### Phase 4: Validation & Polish (2 days = 16 hours)

#### Task 4.1: Validation (8 hours)
- Run all workflows on test data
- Compare Python vs C++ outputs bit-by-bit
- Document any differences

#### Task 4.2: Documentation (4 hours)
- Docstrings for all Python functions
- API reference
- Usage examples

#### Task 4.3: Polish (4 hours)
- Error messages
- Logging
- Performance optimization

---

## Total Time Estimate

| Phase | Hours | Days |
|-------|-------|------|
| Phase 1: Fix issues | 2 | 0.25 |
| Phase 2: Core bindings | 16 | 2 |
| Phase 3: Python layer | 16 | 2 |
| Phase 4: Validation | 16 | 2 |
| **Total** | **50** | **6.25** |

**With 8 hours/day → 6-7 working days**
**With 4-6 hours/day → 10-14 working days**

✅ **2-week timeline is achievable!**

---

## Summary: Architecture Principles

### 1. **Python = Interface**
- CLI parsing
- Configuration
- Orchestration
- Visualization
- Error handling for users

### 2. **C++ = Computation**
- All math stays in C++
- All matrix operations stay in C++
- All MCMC stays in C++
- No performance loss

### 3. **Bindings = Bridge**
- Thin wrapper layer
- Auto type conversion (Eigen ↔ NumPy)
- Progress callbacks
- Exception translation

### 4. **Keep C++ Unchanged**
- Minimal changes to existing code
- Just expose what's needed
- Add new code only in Python

---

## Next Steps Discussion

Before we proceed with fixing the BED exception, what would you like to clarify:

1. **Architecture questions?**
2. **Specific patterns you want to discuss?**
3. **Alternative approaches to consider?**
4. **Priority of tasks?**

Or should I proceed with the detailed plan as outlined?

