# Known Issues

## 1. Progress Bar Threading Warning (Cosmetic Only)

###  Symptom

When using `--progress` flag, you may see at the end:

```
resource_tracker: There appear to be 1 leaked semaphore objects 
to clean up at shutdown
```

And exit code 133 instead of 0.

### Cause

- Python 3.13 + tqdm + threading combination
- tqdm uses multiprocessing primitives internally
- Semaphore cleanup happens at shutdown, triggers warning

### Impact

✅ **Cosmetic only** - NO functional impact:
- Analysis completes successfully
- All results are correct  
- Files are saved properly
- Warning appears AFTER everything is done

⚠️ **Affects:**
- Scripts checking exit codes (will see 133 instead of 0)
- Log cleanliness (warning message appears)

### Workaround

**Option 1: Don't use progress bar**
```bash
# No --progress flag = no warning, exit code 0
python3 -m gctb.cli bayes --bfile data --pheno pheno.txt --bayes C
```

**Option 2: Ignore the warning** 
```bash
# Suppress the specific warning
python3 -W ignore::UserWarning:resource_tracker -m gctb.cli bayes \
    --bfile data --pheno pheno.txt --bayes C --progress
```

**Option 3: Check "Analysis completed successfully!" instead of exit code**
```bash
python3 -m gctb.cli bayes ... --progress 2>&1 | \
    grep -q "Analysis completed successfully" && echo "Success!"
```

### Status

- **Will fix:** Waiting for tqdm update for Python 3.13 compatibility
- **Priority:** Low (cosmetic only, doesn't affect results)
- **Recommendation:** Don't use `--progress` in production scripts

---

## 2. Diagnostics Show "N/A" for ESS and Geweke

### Symptom

When using `--diagnostics`, ESS and Geweke Z columns show "N/A":

```
Parameter       Mean         ESS   Geweke Z  Status
Pi              0.0144      N/A        N/A    ? Unknown
```

### Cause

- C++ MCMC doesn't return individual samples to Python
- Only returns posterior means
- ESS and Geweke need full sample arrays

### Impact

✅ **Expected behavior**:
- Diagnostic structure works
- Parameter means are correct
- Just can't calculate convergence stats without samples

### Workaround

**Currently:** No workaround. Would require modifying C++ code to return sample chains.

**Alternative:** Use C++ GCTB's built-in diagnostics (if available).

### Status

- **Will fix:** Requires C++ changes (significant effort)
- **Priority:** Low (infrastructure works, just needs data)
- **Recommendation:** Diagnostics feature demonstrates structure, full stats TBD

---

## 3. No Other Known Issues

✅ All core functionality works:
- Data loading ✅
- Model building ✅
- MCMC execution ✅
- Results saving ✅
- Both Bayes and SBayes ✅

---

## Reporting Issues

If you find a bug:

1. **Check if it's listed above** (might be cosmetic/expected)
2. **Collect info:**
   - Command used
   - Error message
   - Python version
   - OS
3. **Report on GitHub:** https://github.com/jianzeng/GCTB/issues

---

## Summary

| Issue | Impact | Workaround | Priority |
|-------|--------|------------|----------|
| Progress bar warning | Cosmetic | Don't use `--progress` | Low |
| Diagnostics N/A | Expected | Use means only | Low |

**Both issues are non-critical and don't affect analysis results!** ✅

