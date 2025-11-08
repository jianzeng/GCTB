# GCTB C++ CLI (legacy)

This directory contains the original standalone GCTB command-line application. It is preserved so you (or downstream users) can continue to build and run the native C++ binary even though the Python extension lives elsewhere.

## Contents

- `Makefile` – single-target build of the `gctb` executable.
- `*.cpp/.hpp` – full C++ source for the CLI, data loaders, MCMC engine, LD utilities, etc.
- `deprecated/` – modules that are no longer used by the Python extension (`hsq`, `predict`, MPI helpers) but remain available for the legacy binary.

## Prerequisites

| Dependency | Notes |
|------------|-------|
| C++14 compiler | `clang++` (macOS) or `g++` (Linux) |
| Eigen3 headers | e.g. `/opt/homebrew/include/eigen3` on macOS/Homebrew |
| Boost headers | optional, used for formatting and helpers |
| OpenMP library | `libomp` via Homebrew on macOS |

Set the include paths before building. Example on macOS/Homebrew:

```bash
export EIGEN3_INCLUDE_DIR="/opt/homebrew/include/eigen3"
export BOOST_LIB="/opt/homebrew/include"
```

## Build

```bash
cd scr
make clean
make
```

If the include paths are not exported, you can pass them on the command line:

```bash
make EIGEN_PATH=/opt/homebrew/include/eigen3 BOOST_PATH=/opt/homebrew/include
```

`make` produces the `gctb` binary in the same directory.

## Usage

The CLI behaviour is unchanged from the original project. See the existing documentation or run one of the standard commands, for example:

```bash
./gctb --bfile data --pheno pheno.txt --bayes S --chain-length 1000 --burn-in 100 --out results
```

> **Note:** the legacy binary has no `--help` flag; invoking it without arguments prints the banner and exits.

## Relationship to the Python extension

- The Python bindings use the refactored sources under `../src/`.
- This `scr/` tree is untouched so you can still ship/compare the original C++ package.
- Changes made for the Python extension (e.g. removing `hsq`/`predict` from the bindings) do **not** affect this directory; use it exactly as before.

## Cleaning

To remove build artefacts:

```bash
make clean
```

---

If you need to make adjustments for new compilers or library locations, update the `EIGEN_PATH` and `BOOST_PATH` variables in the Makefile or provide them on the command line.

## Testing with sample data

After building, you can verify the binary with the bundled example dataset:

```bash
# From scr/
./gctb --bfile ../test/data/uk10k_chr1_1mb \
       --pheno ../test/data/test.phen \
       --bayes S --chain-length 1100 --burn-in 100 \
       --seed 12345 --out /tmp/gctb_cpp_example
```

Expected terminal output excerpt (BayesS on 6,717 SNPs, 3,642 individuals):

```
******************************************************************
* GCTB 2.5.2                                                     *
* Genome-wide Complex Trait Bayesian analysis                    *
* Authors: Jian Zeng, Luke Lloyd-Jones, Zhili Zheng, Shouye Liu  *
* MIT License                                                    *
******************************************************************

... progress omitted ...

Posterior statistics from MCMC samples:

              Mean            SD             
        Pi    0.013780        0.002318       
    NnzSnp    91.739998       11.202110      
   SigmaSq    0.835103        0.328335       
         S    -1.116399       0.087933       
    ResVar    112.284843      3.214086       
    GenVar    107.837502      3.376415       
       hsq    0.489883        0.010217       

Computational time: 0:0:4
```

The run completes in a few seconds on modern hardware and writes `.parRes`/`.snpRes` files under `/tmp/gctb_cpp_example*`. Use this to confirm the legacy executable still behaves as expected.
