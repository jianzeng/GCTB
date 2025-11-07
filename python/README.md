# GCTB Python Interface

Python interface for Genome-wide Complex Trait Bayesian Analysis (GCTB) with C++ computational core.

## Quick Start

### Installation

#### Prerequisites

1. **C++ Compiler** with C++14 support (GCC 5+, Clang 3.4+, MSVC 2015+)
2. **CMake** 3.15 or higher
3. **Python** 3.8 or higher
4. **Required libraries:**
   - Eigen3
   - Boost
   - OpenMP
   - pybind11

#### Install dependencies (macOS)

```bash
brew install cmake eigen boost libomp
pip install pybind11
```

#### Install dependencies (Linux)

```bash
# Ubuntu/Debian
sudo apt-get install cmake libeigen3-dev libboost-dev

# CentOS/RHEL
sudo yum install cmake eigen3-devel boost-devel

pip install pybind11
```

#### Build and install GCTB

```bash
# From the repository root
pip install -e .
```

### Basic Usage

```python
import gctb

# Create data object
data = gctb.Data()

# Read PLINK format data
data.read_fam_file("test/data/uk10k_chr1_1mb.fam")
data.read_bim_file("test/data/uk10k_chr1_1mb.bim")
data.read_bed_file(False, "test/data/uk10k_chr1_1mb.bed")

print(f"Loaded {data.num_snps} SNPs and {data.num_inds} individuals")

# Access SNP information
snps = data.get_snp_info_vec()
for i, snp in enumerate(snps[:5]):  # First 5 SNPs
    print(f"SNP {i}: {snp.ID}, chr{snp.chrom}:{snp.physPos}")
```

## Testing

Run tests with pytest:

```bash
cd python
pytest tests/ -v
```

## Development Status

**Current Status:** 🚧 Under active development

### Implemented
- ✅ Basic data structures (SnpInfo, IndInfo, Data)
- ✅ PLINK file I/O
- ✅ Timer utilities

### In Progress
- 🔨 MCMC engine bindings
- 🔨 Model class bindings
- 🔨 CLI interface
- 🔨 Analysis workflows

### Planned
- 📋 Visualization functions
- 📋 Complete CLI parity with C++ version
- 📋 Jupyter notebook examples
- 📋 Comprehensive documentation

## Architecture

```
python/
├── gctb/              # Python package
│   ├── __init__.py    # Package initialization
│   ├── cli.py         # Command-line interface
│   └── _core.so       # C++ extension module (compiled)
├── bindings/          # pybind11 bindings
│   └── bindings.cpp   # C++ to Python bindings
├── tests/             # Test suite
└── examples/          # Example scripts
```

## Contributing

This is a conversion of the existing C++ GCTB codebase to provide a Python interface. The C++ computational core remains unchanged for performance.

## License

MIT License - See LICENSE file for details

## Authors

- Jian Zeng
- Luke Lloyd-Jones
- Zhili Zheng
- Shouye Liu

## References

Original GCTB: http://cnsgenomics.com/software/gctb

