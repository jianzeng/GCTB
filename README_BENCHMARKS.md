# GCTB Benchmark Suite

Complete toolkit for comparing Python and C++ implementations.

---

## 📦 What's Included

| File | Type | Purpose |
|------|------|---------|
| `benchmark_comparison.py` | Python | Comprehensive automated benchmarks |
| `compare_outputs.sh` | Shell | Quick side-by-side comparison |
| `BENCHMARKING_GUIDE.md` | Docs | Detailed instructions, sample output, OpenMP notes |
| `README_BENCHMARKS.md` | Docs | Quick orientation (this file) |

---

## 🚀 Quick Start

### Option 1: Automated Python Benchmark (Recommended)

```bash
cd /Users/haocheng/Github/GCTB
source venv/bin/activate
python3 benchmark_comparison.py
```

**What it does:**
- Runs unit tests (data loading, model creation, MCMC)
- Runs 3 full comparisons (BayesC short/medium, BayesR)
- Compares speed and accuracy
- Generates detailed report (see `BENCHMARKING_GUIDE.md` for interpretation)

**Time:** 3-5 minutes

---

### Option 2: Quick Shell Comparison

```bash
cd /Users/haocheng/Github/GCTB
./compare_outputs.sh
```

**What it does:**
- Runs one BayesS analysis with both implementations
- Shows execution time comparison
- Points to output files for manual inspection

**Time:** ~40 seconds

---

### Option 3: Manual Comparison

**C++ version:**
```bash
./scr/gctb \
    --bfile test/data/uk10k_chr1_1mb \
    --pheno test/data/test.phen \
    --bayes S \
    --chain-length 1100 \
    --burn-in 100 \
    --out cpp_result
```

**Python version:**
```bash
cd python
python3 -m gctb.cli bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes S \
    --chain-length 1100 \
    --burnin 100 \
    --out py_result
```

**Compare:**
```bash
echo "C++ h²:"; tail -1 cpp_result.parRes | awk '{print $1}'
echo "Python h²:"; tail -1 py_result.parRes | awk '{print $1}'
```

(Additional walkthroughs, screenshots, and example diffs live in `BENCHMARKING_GUIDE.md`.)

---

## 📊 What Gets Benchmarked

### Performance Metrics

1. **Data Loading Speed**
   - Reading PLINK files (.fam, .bim, .bed)
   - Parsing phenotypes
   - Memory allocation

2. **Model Creation Speed**
   - Model initialization
   - Prior setup
   - Memory preparation

3. **MCMC Execution Speed**
   - Per-iteration time
   - Total chain time
   - Sampling efficiency

### Accuracy Metrics

1. **Parameter Estimates**
   - Heritability (h²)
   - Proportion non-zero (π)
   - Genetic variance
   - Residual variance
   - Number of non-zero SNPs

2. **SNP Effects**
   - Effect size distribution
   - Posterior inclusion probabilities

### OpenMP Scaling (details in guide)
- Verifies the Python build links OpenMP correctly
- Documents expected speedups for small/medium/large datasets

---

## 🎯 Expected Results

### Performance

**Speed:** Python ≈ C++ (within 5%)

Both implementations use the same C++ computational core. Python adds minimal overhead (~1-2%).

### Accuracy

**Estimates:** Within MCMC variance (~1-5% difference)

MCMC is stochastic (uses random sampling). Different runs will produce slightly different results, even for the same implementation. Longer chains or fixed seeds reduce variation.

---

## 📋 Validation Checklist

Before publishing results, verify:

### ✅ Prerequisites
- [ ] C++ GCTB compiled (`make` in `scr/`)
- [ ] Python package installed (`pip install -e .`)
- [ ] Test data available (`test/data/uk10k_chr1_1mb.*`)
- [ ] Virtual environment activated

### ✅ Performance
- [ ] Python speed within 10% of C++
- [ ] Data loading < 3 seconds
- [ ] MCMC (1100 iter) < 30 seconds
- [ ] No memory leaks

### ✅ Accuracy
- [ ] Heritability within expected range (0.4-0.6 for test data)
- [ ] Pi (proportion non-zero) reasonable (~0.01)
- [ ] SNP effects have correct distribution
- [ ] Output files formatted correctly

### ✅ Reliability
- [ ] No crashes or errors (aside from known cosmetic warnings)
- [ ] Reproducible results
- [ ] Consistent across runs
- [ ] Handles edge cases

---

## 🐛 Troubleshooting

- **“C++ GCTB not found”** → `cd scr && make`
- **“Python package not installed”** → `pip install -e .`
- **“Python much slower”** → check OpenMP (`echo $OMP_NUM_THREADS`), build flags (`-O3`, `-DNDEBUG`), and virtualenv.
- **“Results very different”** → remember stochastic MCMC; run longer chains or multiple replicates.

More detailed guidance (including full command transcripts and example diff outputs) lives in `BENCHMARKING_GUIDE.md`.

---

## 📈 Performance Tips

Use shorter chains for sanity checks, longer chains for production, and `--seed` for reproducibility. Thinning (`--thin`) helps with very long runs. For custom datasets, adapt the scripts in this directory and follow the reporting templates in the detailed guide.

---

## 📝 Reporting Results

When sharing benchmark outcomes, include:
- System specs (CPU, RAM, OS)
- Data size (SNPs × individuals)
- Model type and chain settings
- Execution times for both implementations
- Parameter comparisons (h², π, GenVar, etc.)

Example snippets are provided near the end of `BENCHMARKING_GUIDE.md`.

---

## 🎉 Summary

- `./compare_outputs.sh` – fastest sanity check (~40 s)
- `python3 benchmark_comparison.py` – automated, thorough (~5 min)
- Manual workflow – customize chain length or datasets as needed

Python and C++ outputs should agree within normal MCMC variance and execute in essentially the same wall-clock time. When in doubt, start with the shell script and dig deeper with the detailed guide.

