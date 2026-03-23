# ApproxBayesAPP: Julia vs C++ Implementation Comparison

## Critical Differences Found

### 1. **Data Structure - MAJOR ISSUE**

**Julia Implementation:**
- `betaArray[trait]` has size `my_nsnp * nCategory` (line 285)
- Structure: `[b1_c1, b2_c1, ..., b_nsnp_c1, b1_c2, b2_c2, ..., b_nsnp_c2, ...]`
- Each SNP has **separate beta/delta/alpha values for EACH category** it belongs to
- `markerIndex = (cat - 1) * my_nsnp + true_marker_num` (line 530)

**C++ Implementation:**
- `betaMatrix` has size `nSNPs × 2` (one beta per SNP per trait)
- `deltaMatrix` has size `nSNPs × 2` (one delta per SNP per trait)
- **Only one set of effects per SNP**, not per category

**Impact:** This is a fundamental architectural difference. The Julia code allows SNPs to have different effects in different annotation categories, while the C++ code assumes one effect per SNP.

### 2. **Sampling Loop Structure - CRITICAL**

**Julia Implementation (lines 519-598):**
```julia
for marker = MarkerOrder
    for cat = 1:nCategory  # Loop through ALL categories
        if annoindexm[cat]  # If SNP belongs to this category
            # Sample separate beta/delta/alpha for THIS category
            markerIndex = (cat - 1) * my_nsnp + true_marker_num
            # ... sampling logic ...
        end
    end
end
```

**C++ Implementation (lines 8569-8669):**
```cpp
for (unsigned idx = 0; idx < markerOrder.size(); ++idx) {
    unsigned i = markerOrder[idx];
    // Find which category this SNP belongs to (assume one primary category)
    unsigned cat = 0;
    for (unsigned c = 0; c < annoVec.size(); ++c) {
        if (annoVec[c] != 0.0) {
            cat = c;
            break;  // Only uses FIRST category found!
        }
    }
    // Sample ONCE per SNP, not per category
}
```

**Impact:** C++ only samples once per SNP using the first category found, while Julia samples separately for each category the SNP belongs to.

### 3. **Residual Variance (R) Sampling**

**Julia Implementation (lines 635-651):**
- Samples R per block using `sample_variance_sumstats`
- Has special tuning logic:
  - Checks threshold: `thres = sum(ssq_blk_cat[b][traiti, :]) / totalvarg_blk[b][traiti, traiti]`
  - If `thres > 1.1`, uses sampled R; otherwise sets diagonal to 1.0
  - Tunes covariance: `Rcov = Rcor * sqrt(R_blk[b][1, 1] * R_blk[b][2, 2])`

**C++ Implementation (lines 8429-8477):**
- Samples R per block
- **Missing the threshold check and covariance tuning logic**
- Simply averages R across blocks

**Impact:** R sampling behavior differs, which could affect convergence.

### 4. **Sample Size Handling**

**Julia Implementation:**
- `nInd = my_nGWAS_dict[blk]` is an array `[nInd1, nInd2]` (line 499)
- Uses `nInd[traiti]` and `nInd[traitj]` separately for each trait

**C++ Implementation:**
- `nGWASblocks[blk]` appears to be a single value
- Uses same value for both traits: `float n1 = nGWASblocks[b]; float n2 = nGWASblocks[b];`

**Impact:** If traits have different sample sizes, C++ implementation is incorrect.

### 5. **SSE Computation for Variance Sampling**

**Julia Implementation (lines 751-760):**
```julia
for cat in 1:nCategory
    beta_i = betaArray[traiti][((cat-1)*my_nsnp+1):((cat-1)*my_nsnp+my_nsnp)]
    # Extracts beta values for THIS category from flattened array
    SSE_vec[cat][traiti, traitj] = dot(beta_i, beta_j)
end
```

**C++ Implementation (lines 8705-8719):**
```cpp
for (unsigned c = 0; c < numCategories; ++c) {
    for (unsigned i = 0; i < snpEffects.size; ++i) {
        if (data.annoMat(i, c) != 0.0) {
            // Uses same betaMatrix for all categories
            beta << snpEffects.betaMatrix(i, 0), snpEffects.betaMatrix(i, 1);
            SSE_vec[c](0,0) += beta[0] * beta[0];
        }
    }
}
```

**Impact:** C++ uses the same beta values for all categories, while Julia uses category-specific betas.

### 6. **Genetic Variance Computation**

**Julia Implementation (lines 607-624):**
```julia
for cat in 1:nCategory
    alphaArray_c = [alphaArray[i][(cat-1)*my_nsnp.+SNPIndexb] for i in 1:nTraits]
    # Extracts category-specific alpha values
    what_array_c[traiti][:] = XAb * alphaArray_c[traiti]
    varg_blk_cat[b][traiti, cat] = dot(what_array_c[traiti], what_array_c[traiti])
end
```

**C++ Implementation (lines 8479-8514):**
```cpp
for (unsigned c = 0; c < numCategories; ++c) {
    for (unsigned i = blockStart; i <= blockEnd; ++i) {
        if (annoMatrix(i, c) != 0.0) {
            // Uses same alphaMatrix for all categories
            what1_c += Qi * alphaMatrix(i, 0);
        }
    }
}
```

**Impact:** Same issue - C++ uses one alpha per SNP, Julia uses category-specific alphas.

### 7. **Continuous Annotations (nCon)**

**Julia Implementation:**
- Has special handling for continuous annotations (lines 520-524, 608-612)
- Uses different `xArrayc` and `xpxc` for continuous vs categorical annotations

**C++ Implementation:**
- **No support for continuous annotations**

**Impact:** C++ cannot handle continuous annotations.

### 8. **Pi State Indexing**

**Julia Implementation (lines 580-587):**
```julia
pi_index = 1
for key in keys(BigPi)
    if δ == key
        nLoci_array_vec[cat][pi_index] += 1
    end
    pi_index += 1
end
```
- Uses dictionary keys: `[1.0; 1.0]`, `[1.0; 0.0]`, `[0.0; 1.0]`, `[0.0; 0.0]`

**C++ Implementation (lines 8693):**
```cpp
unsigned idx = (d[0] < 0.5 ? 2 : 0) + (d[1] < 0.5 ? 1 : 0);
```
- States: 0=[1,1], 1=[1,0], 2=[0,1], 3=[0,0]
- **Index mapping appears reversed!** Julia: [1,1]=1, [1,0]=2, [0,1]=3, [0,0]=4
- C++: [1,1]=0, [1,0]=1, [0,1]=2, [0,0]=3

**Impact:** Pi state indexing may be incorrect.

## Additional Issues

### 9. **C12 Formula**

**Julia Implementation (line 545):**
```julia
C12 = Ginv12 + xpxc[marker] * Diagonal(δ[nok]) * Rinv[k, nok]
```
- `nok` is a vector with one element (since nTraits=2)
- `Diagonal(δ[nok])` is a 1×1 matrix
- This simplifies to: `Ginv12 + xpxc[marker] * δ[nok] * Rinv[k, nok]`

**C++ Implementation (line 8625):**
```cpp
float C12 = Ginv12 + xpx * delta_new[nok] * Rinv(k, nok);
```
- Uses `delta_new[nok]` which is correct
- **However**, this uses the NEW delta value, while Julia might be using the current state. Need to verify.

### 10. **Sample Size per Block**

**Julia Implementation:**
- `nInd = my_nGWAS_dict[blk]` is `[nInd1, nInd2]` - array of sample sizes per trait (line 499)
- Uses `nInd[traiti]` and `nInd[traitj]` separately

**C++ Implementation:**
- `nGWASblocks[blk]` is a single float value
- Uses same value for both traits: `float n1 = nGWASblocks[b]; float n2 = nGWASblocks[b];`

**Impact:** If traits have different sample sizes, C++ implementation is incorrect.

## Summary

The C++ implementation **does NOT correctly capture** the Julia implementation. The main issues are:

1. **CRITICAL - Architecture mismatch**: Julia uses category-specific effects per SNP (size `nSNPs × nCategory`), C++ uses one effect per SNP (size `nSNPs × 2`). This is the most fundamental difference.

2. **CRITICAL - Missing category loop**: Julia samples each SNP for EACH category it belongs to, C++ only samples once per SNP using the first category found.

3. **Missing R tuning logic**: C++ doesn't have the threshold check (`thres > 1.1`) and covariance tuning that Julia has.

4. **Sample size handling**: Julia supports different sample sizes per trait per block (`[nInd1, nInd2]`), C++ assumes same sample size.

5. **No continuous annotation support**: C++ doesn't handle continuous annotations (nCon), which Julia fully supports.

6. **Pi indexing**: Need to verify the state mapping is correct (Julia uses dictionary keys, C++ uses array indices).

7. **SSE computation**: C++ uses same beta values for all categories, Julia uses category-specific betas from flattened array.

8. **Genetic variance computation**: Same issue - C++ uses one alpha per SNP, Julia uses category-specific alphas.

**Conclusion:** The C++ implementation appears to be a simplified version that doesn't fully capture the annotation-stratified nature of the model. It needs significant restructuring to match the Julia implementation, particularly:
- Changing data structures to support category-specific effects
- Adding the category loop in SNP effect sampling
- Adding R tuning logic
- Supporting different sample sizes per trait
- Adding continuous annotation support

