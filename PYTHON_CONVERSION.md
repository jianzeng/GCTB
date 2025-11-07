# GCTB Python Conversion Project

## Overview

This document tracks the conversion of GCTB from pure C++ to a Python/C++ hybrid architecture.

**Goal:** Python interface for usability + C++ core for performance

**Timeline:** 2 weeks (aggressive sprint)

## Architecture

```
┌─────────────────────────────────────┐
│         Python Layer                │
│  - CLI (click)                      │
│  - Workflows                        │
│  - Data validation                  │
│  - Visualization                    │
└─────────────────────────────────────┘
              ↕ pybind11
┌─────────────────────────────────────┐
│         C++ Core (unchanged)        │
│  - Data structures                  │
│  - Matrix operations (Eigen)        │
│  - MCMC engine                      │
│  - Statistical computations         │
└─────────────────────────────────────┘
```

## Project Structure

```
GCTB/
├── scr/                      # Original C++ source (unchanged)
├── python/
│   ├── gctb/                 # Python package
│   │   ├── __init__.py       # Package initialization
│   │   ├── cli.py            # CLI (TODO)
│   │   └── workflows.py      # Analysis workflows (TODO)
│   ├── bindings/
│   │   └── bindings.cpp      # pybind11 bindings
│   ├── tests/                # Python tests
│   └── examples/             # Usage examples (TODO)
├── CMakeLists.txt            # Build system
├── setup.py                  # Python packaging
├── build.sh                  # Build script
└── requirements.txt          # Python dependencies
```

## Progress Tracker

### ✅ Phase 1: Foundation (Day 1) - COMPLETED
- [x] Create Python branch
- [x] Set up project structure
- [x] Create CMakeLists.txt (with macOS SDK configuration)
- [x] Create setup.py
- [x] Create pyproject.toml
- [x] Create initial pybind11 bindings
- [x] Create test framework
- [x] Create build script
- [x] Virtual environment set up
- [x] All dependencies installed
- [x] CMake configuration successful

**Status:** ⚠️ BLOCKED by system issue - Xcode CommandLineTools broken

**Issue:** C++ compiler cannot find standard library headers. This is a system-level problem, not our code.
**Solution:** User needs to run: `sudo rm -rf /Library/Developer/CommandLineTools && xcode-select --install`

**Time spent:** ~3 hours (2 hours on system debugging)
**Once fixed:** Compilation should work immediately, tests can run

### 🔨 Phase 2: Core Bindings (Days 2-7) - IN PROGRESS
- [ ] Test initial build
- [ ] Fix compilation issues
- [ ] Bind Data class (complete)
- [ ] Bind Model class
- [ ] Bind MCMC class
- [ ] Bind Options class
- [ ] Bind utility classes

### 📋 Phase 3: Workflows (Days 8-10) - PLANNED
- [ ] Convert options.cpp → Python CLI
- [ ] Implement Bayes workflow
- [ ] Implement SBayes workflow
- [ ] Implement LDmatrix workflow
- [ ] Add progress bars (tqdm)

### 📋 Phase 4: Testing & Validation (Days 11-14) - PLANNED
- [ ] Run all test cases
- [ ] Validate outputs match C++ version
- [ ] Performance benchmarking
- [ ] Documentation
- [ ] Examples

## Currently Bound Classes

### SnpInfo
- ✅ Basic structure
- ✅ All major fields accessible
- ❓ Methods need testing

### IndInfo
- ✅ Basic structure
- ✅ All major fields accessible
- ❓ Methods need testing

### Data
- ✅ Constructor
- ✅ File I/O methods (fam, bim, bed, pheno)
- ✅ Basic properties
- ⚠️ Many methods still need binding
- ❓ Matrix access patterns need testing

### Timer
- ✅ Full functionality

## Next Steps

### Immediate (Next 2 hours)
1. Run `./build.sh` to test compilation
2. Fix any compilation errors
3. Run basic tests
4. Verify can import gctb in Python

### Tomorrow (Day 2)
1. Load test data successfully
2. Expand Data class bindings
3. Start Model class bindings
4. Begin MCMC class bindings

### This Week (Days 3-7)
1. Complete all core class bindings
2. Start Python CLI development
3. Implement first complete workflow (Bayes)

## Build Instructions

### Quick build:
```bash
./build.sh
```

### Manual build:
```bash
# Set environment variables
export EIGEN3_INCLUDE_DIR=/usr/local/include/eigen3
export BOOST_LIB=/usr/local/include

# Install
pip install -e .
```

### Test:
```bash
python3 -c "import gctb; print(gctb.__version__)"
cd python && pytest tests/ -v
```

## Known Issues

- [ ] None yet (initial setup)

## Notes

- Keeping C++ code completely unchanged for safety
- Using pybind11 for automatic type conversions
- Eigen matrices should convert to numpy automatically
- Need to test GIL release for long-running computations
- MPI support deferred to later phase

## References

- pybind11 docs: https://pybind11.readthedocs.io/
- Eigen integration: https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
- Original GCTB: http://cnsgenomics.com/software/gctb

