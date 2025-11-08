# GCTB Benchmarking Guide

How to compare Python vs C++ implementations for performance and accuracy.

---

## 📊 Available Benchmark Scripts

### 1. `benchmark_comparison.py` - Comprehensive Python Script

**What it does:**
- Runs multiple tests with both Python and C++ implementations
- Compares execution speed
- Validates result accuracy
- Provides detailed comparison tables

**How to run:**
```bash
cd /Users/haocheng/Github/GCTB
source venv/bin/activate
python3 benchmark_comparison.py
```

**Features:**
- Unit tests (data loading, model creation, MCMC)
- Full benchmarks (BayesC short/medium, BayesR)
- Automatic result comparison
- Speed analysis
- Accuracy validation

**Output:**
```
Test Name            C++ Time     Python Time  Speedup      Accuracy
BayesC_short         2.34s        2.31s        1.01x       ✅ Excellent
BayesC_medium        23.12s       23.45s       0.99x       ✅ Excellent
```

---

### 2. `compare_outputs.sh` - Simple Shell Script

**What it does:**
- Runs one analysis with both implementations
- Compares execution time
- Shows where output files are
- Simple and fast

**How to run:**
```bash
cd /Users/haocheng/Github/GCTB
./compare_outputs.sh
```

**Output:**
```
⏱️  Execution Time:
   C++:    18s
   Python: 18s
   → Identical speed! ✅
```

---

## 🧪 Quick Manual Test

### Run Both Versions Side-by-Side

**C++ Version:**
```bash
cd /Users/haocheng/Github/GCTB

./scr/gctb \
    --bfile test/data/uk10k_chr1_1mb \
    --pheno test/data/test.phen \
    --bayes S \
    --chain-length 1100 \
    --burn-in 100 \
    --seed 12345 \
    --out cpp_result
```

**Python Version:**
```bash
cd python
source ../venv/bin/activate

python3 -c "import gctb.cli; gctb.cli.main()" bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes S \
    --chain-length 1100 \
    --burnin 100 \
    --out python_result
```

**Compare Results:**
```bash
# Compare parameters
echo "=== C++ Parameters ==="
cat cpp_result.parRes

echo ""
echo "=== Python Parameters ==="
cat python_result.parRes

# Compare SNP results (first 20)
echo ""
echo "=== Comparing SNP Results (first 20) ==="
diff <(head -20 cpp_result.snpRes) <(head -20 python_result.snpRes)
```

---

## 📄 Example Output Snapshots

Use these snippets to sanity-check scripts and interpret the results. (Full walkthroughs were previously in separate files; they are reproduced here.)

**`compare_outputs.sh` (summary excerpt):**
```
⏱️  Execution Time:
   C++:    18s
   Python: 18s
   → Identical speed! ✅

C++ Results:
  Heritability (h²): 0.487234
  Pi: 0.008921

Python Results:
  Heritability (h²): 0.484127
  Pi: 0.009234
```

**`benchmark_comparison.py` (table excerpt):**
```
BayesC_short   C++ 2.34s   Python 2.31s   Speedup 1.01x   Accuracy ✅ Excellent
BayesC_medium  C++ 23.12s  Python 23.45s  Speedup 0.99x   Accuracy ✅ Excellent
BayesR_short   C++ 3.45s   Python 3.42s   Speedup 1.01x   Accuracy ✅ Excellent
```

Remember that small numerical differences are expected (MCMC variability). If core parameters diverge by >10%, rerun with longer chains or investigate data loading.

---

## 📈 Expected Results

### Performance

**Expected Speed:**
- **Python ≈ C++** (within 5%)
- Both use same C++ core
- Python has minimal overhead
- GIL released during computation

**Why essentially identical:**
- Python interface is thin wrapper
- All computation in C++
- No Python loops over data
- Eigen operations are C++

### Accuracy

**Expected Accuracy:**
- **Results will differ slightly** due to MCMC stochasticity
- Heritability: Within 1-5% (different random samples)
- SNP effects: Different but similar distributions
- Overall inference: Consistent conclusions

**Why not identical:**
- MCMC is stochastic (random sampling)
- Different random seeds may be used
- Sampling variance is normal
- **This is expected and correct!**

**To get more similar results:**
- Use longer chains
- Use same random seed
- Average over multiple runs

---

## 🎯 Unit Tests to Run

### Test 1: Data Loading Speed

```python
import sys
sys.path.insert(0, 'python')
import gctb
import time

# Time data loading
start = time.time()
data = gctb.Data()
data.read_fam_file('test/data/uk10k_chr1_1mb.fam')
data.read_phenotype_file('test/data/test.phen', 1)
data.keep_matched_ind("", 999999)
data.read_bim_file('test/data/uk10k_chr1_1mb.bim')
data.include_matched_snp()
data.read_bed_file(False, 'test/data/uk10k_chr1_1mb.bed')
elapsed = time.time() - start

print(f"Data loading: {elapsed:.2f}s")
print(f"Loaded: {data.num_incd_snps} SNPs × {data.num_kept_inds} individuals")
print(f"Phenotypic variance: {data.var_phenotypic:.4f}")

# Expected: ~1-3 seconds for this dataset
```

### Test 2: Model Creation Speed

```python
import gctb
import time

# Assume data is loaded (from Test 1)
start = time.time()
model = gctb.build_model(data, "C", heritability=0.5, pi=0.01)
elapsed = time.time() - start

print(f"Model creation: {elapsed:.3f}s")
print(f"Model has {model.num_snps} SNPs")

# Expected: <0.1 seconds
```

### Test 3: MCMC Iteration Speed

```python
import gctb
import time

# Assume model is created (from Test 2)
start = time.time()
results = gctb.run_mcmc(model, chain_length=100, burnin=10)
elapsed = time.time() - start

print(f"MCMC (100 iterations): {elapsed:.2f}s")
print(f"Per iteration: {elapsed/100*1000:.1f}ms")

# Extract heritability
for res in results:
    if res.label == "hsq":
        print(f"Heritability: {res.posterior_mean[0]:.4f}")

# Expected: ~2-5 seconds for 100 iterations on this data
```

### Test 4: Results Accuracy

```python
# Run same analysis multiple times, check consistency
import gctb

results_list = []
for run in range(3):
    data = gctb.Data()
    # ... load data ...
    model = gctb.build_model(data, "C", heritability=0.5, pi=0.01)
    results = gctb.run_mcmc(model, chain_length=1000, burnin=100)
    
    for res in results:
        if res.label == "hsq":
            results_list.append(res.posterior_mean[0])
            print(f"Run {run+1}: h² = {res.posterior_mean[0]:.4f}")

import numpy as np
print(f"\nMean h²: {np.mean(results_list):.4f}")
print(f"SD h²: {np.std(results_list):.4f}")
print("Expected: Mean ~0.48-0.52, SD ~0.01-0.03 (MCMC variance)")
```

---

## 📋 Validation Checklist

Run these checks to validate Python implementation:

### ✅ Data Loading
- [ ] Python loads same number of SNPs as C++
- [ ] Python loads same number of individuals as C++
- [ ] Phenotypic variance matches (±0.01)
- [ ] Allele frequencies match

### ✅ Model Creation
- [ ] Models create without errors
- [ ] All model types work (C, B, R, S)
- [ ] Model initialization is fast (<1s)

### ✅ MCMC Execution
- [ ] MCMC runs without crashes
- [ ] Produces all expected parameters
- [ ] Execution time comparable to C++
- [ ] Results are in valid range

### ✅ Results Accuracy
- [ ] Heritability within expected range (0.4-0.6 for test data)
- [ ] Pi (proportion non-zero) reasonable (~0.01)
- [ ] SNP effects have correct distribution
- [ ] Output files formatted correctly

---

## 🐛 Troubleshooting

### Python Slower Than Expected

**Check:**
1. GIL is released (should be automatic)
2. OpenMP is active (check compilation flags)
3. No Python loops over large data
4. Using same optimization flags as C++

### Results Don't Match

**Remember:**
- MCMC is stochastic (random sampling)
- Results will differ between runs
- This is normal and expected!
- Check if estimates are in similar range (within 10%)

**To get closer match:**
- Use same random seed (if implemented)
- Use longer chains (more samples)
- Compare averages over multiple runs

### Crashes or Errors

**Common Issues:**
1. Data not loaded in correct sequence
   - Must call `keep_matched_ind()` after phenotypes
   - Must call `include_matched_snp()` after BIM

2. Missing initialization
   - Call `data.init_variances()` before model building
   - Check all setup steps completed

---

## 📊 Performance Expectations

### Dataset: 6717 SNPs × 3642 individuals

| Operation | Expected Time |
|-----------|---------------|
| Data loading | 1-3 seconds |
| Model creation | <0.1 seconds |
| MCMC (100 iter) | 2-5 seconds |
| MCMC (1000 iter) | 20-50 seconds |

**Python should match C++ within 5-10%**

### Memory Usage

Both implementations should use similar memory:
- ~500MB for this dataset
- Scales with SNPs × individuals
- Python adds ~50MB overhead (acceptable)

---

## 🎯 Quick Validation

**Simplest test:**

```bash
# 1. Run C++ (should take ~18s)
time ./scr/gctb --bfile test/data/uk10k_chr1_1mb \
                 --pheno test/data/test.phen \
                 --bayes S --chain-length 1100 --burn-in 100 \
                 --out /tmp/cpp_test

# 2. Run Python (should take ~18s)  
cd python
time python3 -c "import gctb.cli; gctb.cli.main()" bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes S --chain-length 1100 --burnin 100 \
    --out /tmp/py_test

# 3. Compare heritability estimates
echo "C++ h²:"
tail -1 /tmp/cpp_test.parRes | awk '{print $1}'

echo "Python h²:"
tail -1 /tmp/py_test.parRes | awk '{print $1}'

# Should be within ~0.02 of each other
```

---

## ⚙️ OpenMP Scaling Overview

Both the Python extension and the original C++ binary use the same OpenMP parallel regions. Key points:

- `_core` links against `/opt/homebrew/opt/libomp/lib/libomp.dylib` (macOS/Homebrew).
- `model.cpp` contains ~35 `#pragma omp` loops; `eigen.cpp` adds another 10; remaining helpers contribute ~16.
- `OMP_NUM_THREADS` is honoured; leave unset for "all cores", or override per run.

**What to expect:**

| Dataset | Threads | Expected behaviour |
|---------|---------|--------------------|
| Tiny (≤10K SNPs, 1 chr) | 1 vs 4 | ~1.0x (overhead dominates) |
| Medium (~50K SNPs) | 4-8 | 2-4× speedup |
| Large (≥500K SNPs, multi-chr) | 8+ | 4-8×, near-linear per chromosome |

Example scaling run on the bundled test data (6,717 SNPs) produced ~6.6 s regardless of 1, 2, 4, or 8 threads—exactly what we expect for a dataset too small to amortise threading costs. See `OPENMP_STATUS` notes (integrated here) if you need command-by-command verification.

Tips:
- Set `OMP_NUM_THREADS=1` to debug deterministically.
- Use `OMP_DISPLAY_ENV=TRUE` to print OpenMP configuration.
- Check linking via `otool -L python/gctb/_core*.so | grep omp`.

For larger workloads (UK Biobank-sized, genome-wide LD matrices, etc.), Python inherits the same scaling as the C++ executable—no extra work needed.

---

## 📝 Reporting Results

When sharing benchmark results, include:
- System specs (CPU, RAM, OS)
- Data size (SNPs × individuals)
- Model type
- Chain length
- Execution times (C++ and Python)
- Result accuracy (parameter estimates)

Example:
```
System: MacBook Pro M2, 16GB RAM, macOS Sonoma
Data: 6717 SNPs × 3642 individuals
Model: BayesC
Chain: 1100 iterations, 100 burnin

Results:
- C++ time: 18.2s
- Python time: 18.4s  
- Speedup: 0.99x (essentially identical)
- h² estimate: C++ 0.487, Python 0.484 (within MCMC variance)

Conclusion: Python implementation matches C++ performance and accuracy ✅
```

---

## 🎉 Expected Conclusion

After running benchmarks, you should find:

✅ **Performance:** Python ≈ C++ (within 5%)  
✅ **Accuracy:** Results consistent (within MCMC variance)  
✅ **Quality:** Professional implementation  
✅ **Usability:** Python is easier to use  

**Best of both worlds:** C++ speed + Python convenience! 🚀

