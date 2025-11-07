# Progress Bars and Convergence Diagnostics Guide

**Version 1.1 Features** - Optional progress tracking and MCMC diagnostics

---

## 🎯 Overview

Two new **optional** features have been added to GCTB Python v1.1:

1. **Progress Bars** - Visual feedback during long MCMC runs
2. **Convergence Diagnostics** - Statistical tests to assess MCMC quality

**Both are OFF by default** to maintain backward compatibility.

---

## 📊 Progress Bars

### What They Do

Show estimated progress during MCMC execution:

```
  MCMC Progress:  45%|████████░░░░░░░░| 4500/10000 [03:25<04:10, 20ms/iter]
```

### How to Enable

**CLI:**
```bash
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --progress
```

**Python API:**
Currently only available through CLI.

### Important Notes

⚠️ **The progress bar is an estimate**
- MCMC runs in C++ (faster than Python)
- GIL is released, so Python can't track exact progress
- Progress bar shows estimated completion based on time
- Actual progress may vary slightly

✅ **When to use:**
- Long chains (>5000 iterations)
- Interactive sessions
- Want to monitor progress

❌ **When NOT to use:**
- Batch/script mode
- Short chains (<1000 iterations)
- Redirecting output to files

### Requirements

Requires `tqdm` package:
```bash
pip install tqdm
```

Already included in `requirements.txt`.

---

## 🔬 Convergence Diagnostics

### What They Do

Assess whether your MCMC chain has converged:

```
==================================================
Convergence Diagnostics:
==================================================
Chain length: 1000, Burn-in: 100
Effective samples: 900

Parameter       Mean         ESS   Geweke Z  Status
--------------------------------------------------
hsq             0.4876       145     -1.23    ✓ Good
Pi              0.0089        98      0.87    ✓ Good
GenVar          0.4521       156     -0.92    ✓ Good
ResVar          0.4762       134      1.12    ✓ Good
```

### How to Enable

**CLI:**
```bash
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --diagnostics
```

**Python API:**
```python
import gctb
from gctb.diagnostics import check_convergence, recommend_chain_length

# Run MCMC
results = gctb.run_mcmc(model, chain_length=1000, burnin=100)

# Check convergence
diag = check_convergence(results, chain_length=1000, burnin=100, verbose=True)

# Get recommendations
needs_longer, message, _ = recommend_chain_length(diag)
if needs_longer:
    print(message)
```

### Metrics Explained

#### 1. **Effective Sample Size (ESS)**

Measures independent information in correlated samples.

**Interpretation:**
- **> 400:** Excellent - highly reliable estimates
- **> 100:** Good - sufficient for most applications
- **> 50:** Fair - acceptable but consider longer chain
- **< 50:** Poor - chain too short, run longer

**Why it matters:**
- MCMC samples are correlated (not independent)
- ESS tells you "effective" number of independent samples
- Higher ESS = more reliable parameter estimates

#### 2. **Geweke Z-Score**

Compares early vs late portions of the chain.

**Interpretation:**
- **|Z| < 2:** Good - chain converged
- **|Z| < 3:** Fair - possibly converged
- **|Z| > 3:** Poor - chain not converged, run longer

**Why it matters:**
- Checks if chain reached stationary distribution
- Early and late parts should have similar means
- Large Z indicates chain still "wandering"

### Requirements

Requires `scipy` package:
```bash
pip install scipy
```

Already included in `requirements.txt`.

---

## 🚀 Usage Examples

### Example 1: Quick Analysis (No Extras)

```bash
# Default: No progress bar, no diagnostics (fast)
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 1000 --out results
```

**Best for:**
- Batch jobs
- Scripts
- Quick tests

---

### Example 2: Interactive Analysis (Progress Bar)

```bash
# Enable progress bar to monitor long runs
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 10000 --progress --out results
```

**Best for:**
- Interactive sessions
- Long chains
- Want to see progress

---

### Example 3: Quality Check (Diagnostics)

```bash
# Enable diagnostics to check convergence
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 1000 --diagnostics --out results
```

**Best for:**
- First analysis of new data
- Publication-quality results
- Unsure of chain length

---

### Example 4: Full Features (Both)

```bash
# Enable both progress and diagnostics
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 5000 --progress --diagnostics --out results
```

**Best for:**
- Important analyses
- Want to monitor AND validate
- Interactive quality-controlled runs

---

## 📋 When to Use Each Feature

### Use Progress Bars When:

✅ Running interactively (terminal)  
✅ Chain > 5000 iterations  
✅ Want to see time remaining  
✅ Checking if job is still running  

❌ **Don't use when:**
- Running in scripts/batch
- Output redirected to file
- Chain < 1000 iterations (too fast)

### Use Diagnostics When:

✅ First analysis of new data  
✅ Unsure of appropriate chain length  
✅ Want publication-quality results  
✅ Results seem unexpected  

❌ **Don't use when:**
- Very confident in chain length
- Batch processing many datasets
- Want fastest possible run

---

## 🎓 Understanding Diagnostics Output

### Good Convergence Example

```
Parameter       Mean         ESS   Geweke Z  Status
--------------------------------------------------
hsq             0.4876       245     -0.83    ✓ Good
Pi              0.0089       198      0.67    ✓ Good
```

**Interpretation:**
- ✅ ESS > 100 (sufficient samples)
- ✅ |Z| < 2 (converged)
- ✅ Chain is adequate

**Action:** Proceed with analysis ✓

---

### Fair Convergence Example

```
Parameter       Mean         ESS   Geweke Z  Status
--------------------------------------------------
hsq             0.4912        78      1.45    ○ Fair
Pi              0.0087        65      2.12    ○ Fair
```

**Interpretation:**
- ⚠️ ESS 50-100 (borderline)
- ⚠️ |Z| < 3 (possibly converged)
- ⚠️ Results likely OK but not optimal

**Action:** Consider 2x longer chain for better precision

---

### Poor Convergence Example

```
Parameter       Mean         ESS   Geweke Z  Status
--------------------------------------------------
hsq             0.5234        23      3.87    ⚠ Poor
Pi              0.0102        19      4.23    ⚠ Poor
```

**Interpretation:**
- ❌ ESS < 50 (insufficient)
- ❌ |Z| > 3 (not converged)
- ❌ Chain too short

**Action:** Run 3-5x longer chain! Current results unreliable.

---

## 🔧 Troubleshooting

### "Warning: tqdm not installed"

**Problem:** Progress bar feature unavailable.

**Solution:**
```bash
pip install tqdm
```

---

### "Warning: scipy not installed"

**Problem:** Diagnostics feature unavailable.

**Solution:**
```bash
pip install scipy
```

---

### Progress bar doesn't move smoothly

**This is normal!**
- Progress bar is estimated (C++ runs independently)
- Updates every 0.5 seconds
- Final completion is accurate

---

### Diagnostics say "Poor" but results look reasonable

**Possible reasons:**
1. Chain genuinely too short (run longer)
2. Complex posterior (needs more samples)
3. Thin parameter too high (losing samples)

**Solutions:**
- Increase chain length 2-5x
- Reduce thinning (thin=5 instead of thin=10)
- Run multiple chains and average

---

### ESS values are very low (<20)

**This is serious!**

**Causes:**
- Chain too short
- High autocorrelation in samples
- Model complexity too high for data

**Solutions:**
1. Increase chain length dramatically (10x)
2. Check if model is appropriate for data
3. Consider different prior settings

---

## 💡 Best Practices

### 1. Start with Diagnostics

When analyzing new data:
```bash
# First run: check convergence
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 2000 --diagnostics --out test
```

If diagnostics are poor, increase chain length:
```bash
# Second run: longer chain
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 10000 --diagnostics --out final
```

---

### 2. Use Progress for Long Runs

Once you know the right chain length:
```bash
# Production run with progress
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt \
    --bayes C --chain-length 10000 --progress --out results
```

---

### 3. Batch Mode: No Extras

For batch processing:
```bash
# Fast batch processing (no overhead)
for i in {1..100}; do
    python3 -m gctb.cli bayes --bfile data_$i --pheno pheno_$i.txt \
        --bayes C --chain-length 5000 --out results_$i
done
```

---

## 📊 Performance Impact

### Progress Bar

**Overhead:** ~2-5%
- Minimal impact (runs in Python while C++ computes)
- Slight threading overhead
- Negligible for long chains

### Diagnostics

**Overhead:** <1%
- Runs after MCMC (not during)
- Only processes parameter samples (not SNPs)
- Negligible impact

**Both features have minimal performance cost!**

---

## 🎉 Summary

### Quick Reference

| Feature | Flag | Default | Best For |
|---------|------|---------|----------|
| Progress Bar | `--progress` | OFF | Interactive, long runs |
| Diagnostics | `--diagnostics` | OFF | Quality checks, new data |

### Recommended Usage

**For exploration:**
```bash
--progress --diagnostics
```

**For production:**
```bash
# (no flags, fastest)
```

**For publication:**
```bash
--diagnostics  # validate, but no progress bar
```

---

## 📚 Additional Resources

- **GCTB Documentation:** `README_PYTHON.md`
- **CLI Reference:** `python3 -m gctb.cli --help`
- **API Reference:** `python3 -c "import gctb; help(gctb.diagnostics)"`

---

**Features are optional by default - use them when they help!** ✅

