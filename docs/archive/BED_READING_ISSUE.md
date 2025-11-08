# BED File Reading Issue - TODO

## Status: DEFERRED (Fix Later)

**Decision Date:** Nov 7, 2025
**Decision:** Defer to after Model/MCMC bindings are complete
**Priority:** Medium (needed for individual-level analysis, not for SBayes)

---

## Issue Description

**Symptom:** Segmentation fault when reading PLINK BED files

**Current State:**
```python
data = gctb.Data()
data.read_fam_file('test.fam')           # ✅ Works
data.build_kept_individuals()            # ✅ Works
data.read_bim_file('test.bim')           # ✅ Works
data.include_matched_snp()               # ✅ Works
data.read_bed_file(False, 'test.bed')    # ❌ Segfault

# Output before crash:
# "Reading PLINK BED file from [test.bed] in SNP-major format ..."
# Segmentation fault: 11
```

**What Works:**
- ✅ All error checks pass (file exists, format valid, SNPs/individuals initialized)
- ✅ File opens successfully
- ✅ Header reads correctly
- ✅ Gets past all validation

**What Fails:**
- ❌ Crashes during actual genotype reading
- ❌ Python process terminates (segfault)

---

## Technical Details

**Code Location:** `scr/data.cpp` lines 194-290

**Key Code:**
```cpp
void Data::readBedFile(const bool noscale, const string &bedFile) {
    // Lines 197-198: Checks pass ✓
    if (numIncdSnps == 0) throw ("Error: No SNP is retained");
    if (numKeptInds == 0) throw ("Error: No individual is retained");
    
    // Lines 200-202: Matrix allocation
    Z.resize(numKeptInds, numIncdSnps);      // Genotype matrix
    ZPZdiag.resize(numIncdSnps);             // Diagonal  
    snp2pq.resize(numIncdSnps);              // Allele frequencies
    
    // Lines 205+: File reading
    FILE *in = fopen(bedFile.c_str(), "rb");  // Opens successfully
    // ... somewhere after this, crashes
}
```

**Current Bindings:**
```cpp
// python/bindings/bindings.cpp line 105-115
.def("read_bed_file", [](Data& data, bool noscale, const std::string& bed_file) {
    try {
        data.readBedFile(noscale, bed_file);
    } catch (const std::string& e) {
        throw std::runtime_error(e);
    } catch (const char* e) {
        throw std::runtime_error(e);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Error reading BED file: ") + e.what());
    }
}, py::arg("noscale"), py::arg("bed_file"), "Read PLINK .bed file")
```

---

## Hypotheses (To Investigate Later)

### Hypothesis 1: Missing Matrix Initialization
**Evidence:** 
- C++ code directly does `Z.resize(numKeptInds, numIncdSnps)`
- May need additional Eigen matrix setup

**Test:**
```cpp
// Try explicitly initializing:
data.Z.resize(numKeptInds, numIncdSnps);
data.Z.setZero();  // Explicit initialization
```

### Hypothesis 2: Vector Mismatch
**Evidence:**
- Code iterates over `incdSnpInfoVec` and `keptIndInfoVec`
- These vectors might not be properly sized/initialized

**Test:**
```python
# Check before calling read_bed_file:
print(f"numIncdSnps: {data.num_incd_snps}")
print(f"numKeptInds: {data.num_kept_inds}")
print(f"incdSnpInfoVec size: {len(data.get_incd_snp_info_vec())}")
print(f"keptIndInfoVec size: {len(data.get_kept_ind_info_vec())}")
```

### Hypothesis 3: Missing Prerequisite Call
**Evidence:**
- `gctb.cpp` shows specific sequence before `readBedFile()`
- May need `data.initVariances()` or similar

**Test:**
```python
# Try calling in exact C++ sequence:
data.read_fam_file(...)
data.read_bim_file(...)
data.includeMatchedSnp()  # ✓ We do this
data.initVariances(0.5, 0.0)  # ← Missing? Need to bind this
data.readBedFile(...)
```

### Hypothesis 4: Pointer Access Issues
**Evidence:**
- Code accesses `snpInfo->included`, `indInfo->kept` flags
- If vectors aren't fully populated, could crash

**Test:**
```cpp
// Add safety checks in C++ code or Python wrapper:
for (auto* snp : data.incdSnpInfoVec) {
    if (snp == nullptr) throw("Null SNP pointer");
}
```

---

## How to Debug (When We Come Back)

### Step 1: Add Debug Output
```cpp
// Modify data.cpp temporarily:
void Data::readBedFile(...) {
    cout << "DEBUG: Starting readBedFile" << endl;
    cout << "DEBUG: numIncdSnps = " << numIncdSnps << endl;
    cout << "DEBUG: numKeptInds = " << numKeptInds << endl;
    
    Z.resize(numKeptInds, numIncdSnps);
    cout << "DEBUG: Z.resize done" << endl;
    
    // ... add more debug lines to find crash location
}
```

### Step 2: Check C++ Workflow
```bash
# Look at how original C++ CLI does it:
grep -A 20 "readBedFile" scr/gctb.cpp
# See what setup steps happen before readBedFile
```

### Step 3: Bind Missing Methods
```cpp
// May need to add:
.def("init_variances", &Data::initVariances)
.def("compute_allele_freq", &Data::computeAlleleFreq)
// etc.
```

### Step 4: Test with Smaller Data
```python
# Create tiny test file:
# 10 SNPs x 10 individuals
# Easier to debug with small data
```

---

## Workarounds (Current)

### For SBayes (Summary Statistics):
```python
# No genotypes needed!
data = gctb.Data()
data.read_gwas_summary_file("summary.txt")
data.read_ld_matrix_file("ldm.bin")
# Ready for SBayes analysis
```

### For Testing Without Genotypes:
```python
# Can test most functionality without BED:
data = gctb.Data()
data.read_fam_file("test.fam")
data.read_bim_file("test.bim")
data.build_kept_individuals()
data.include_matched_snp()
# Have metadata, just no genotypes
```

---

## When to Fix

**Trigger Conditions:**
1. ✅ Model class bindings complete
2. ✅ MCMC class bindings complete
3. ✅ SBayes workflow working
4. Starting Bayes (individual-level) workflow

**Or:**
- User specifically needs individual-level analysis
- Someone has insight into the issue
- Quick fix becomes apparent

**Estimated Fix Time:** 2-4 hours of focused debugging

---

## Related Files

- `scr/data.cpp` (lines 194-290) - The actual reading code
- `scr/data.hpp` (line 541) - Function declaration
- `scr/gctb.cpp` (lines 20-40) - How it's called in original code
- `python/bindings/bindings.cpp` (lines 105-115) - Current binding
- `python/tests/test_data_io.py` - Test that currently fails

---

## Notes for Future Me

1. **Don't spend too long on this** - If not obvious after 2 hours, ask for help
2. **Check valgrind** - Memory debugging tool can pinpoint exact issue
3. **Compare with original** - Run original C++ CLI to confirm test data is valid
4. **Ask on Stack Overflow** - Eigen + pybind11 issue might be known
5. **Consider alternative** - Could wrap at higher level (whole workflow, not individual method)

---

## Success Criteria (When Fixed)

```python
# This should work:
import gctb

data = gctb.Data()
data.read_fam_file('test/data/uk10k_chr1_1mb.fam')
data.build_kept_individuals()
data.read_bim_file('test/data/uk10k_chr1_1mb.bim')
data.include_matched_snp()
data.read_bed_file(False, 'test/data/uk10k_chr1_1mb.bed')  # ← Should work

print(f"✅ Loaded {data.num_incd_snps} SNPs x {data.num_kept_inds} individuals")
# Should print: ✅ Loaded 6717 SNPs x 3642 individuals
```

---

**Last Updated:** Nov 7, 2025
**Status:** Documented, deferred, will fix after Model/MCMC

