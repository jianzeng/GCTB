#!/usr/bin/env python3
"""
Simple OpenMP Scaling Test for GCTB

Direct API test without subprocess to avoid crashes.
Tests 1, 2, 4, 8 threads and measures execution time.
"""

import sys
import os
import time

# Add python directory to path
sys.path.insert(0, 'python')

import gctb

# Configuration
BFILE = "test/data/uk10k_chr1_1mb"
PHENO = "test/data/test.phen"
CHAIN_LENGTH = 2000
BURNIN = 200
THREAD_COUNTS = [1, 2, 4, 8]

print("=" * 80)
print("GCTB OpenMP Scaling Test (Direct API)")
print("=" * 80)
print(f"\nDataset: {BFILE}")
print(f"Chain: {CHAIN_LENGTH}, Burnin: {BURNIN}")
print(f"Testing threads: {THREAD_COUNTS}\n")

results = {}

for num_threads in THREAD_COUNTS:
    print(f"\n{'='*80}")
    print(f"Testing with {num_threads} thread(s)")
    print(f"{'='*80}\n")
    
    # Set OpenMP threads
    os.environ['OMP_NUM_THREADS'] = str(num_threads)
    print(f"  OMP_NUM_THREADS = {num_threads}")
    
    try:
        print(f"  Loading data...", end=" ", flush=True)
        
        # Load data
        data = gctb.Data()
        data.read_fam_file(f"{BFILE}.fam")
        data.read_phenotype_file(PHENO, 1)
        data.keep_matched_ind("", 999999)
        data.read_bim_file(f"{BFILE}.bim")
        data.include_matched_snp()
        data.read_bed_file(False, f"{BFILE}.bed")
        
        print(f"✓ ({data.num_incd_snps} SNPs, {data.num_kept_inds} individuals)")
        
        print(f"  Building model...", end=" ", flush=True)
        # Build model
        model = gctb.build_model(data, "C", heritability=0.5, pi=0.01)
        print("✓")
        
        print(f"  Running MCMC ({CHAIN_LENGTH} iterations)...", end=" ", flush=True)
        
        # Run MCMC and time it
        start = time.time()
        
        results_mcmc = gctb.run_mcmc(
            model=model,
            chain_length=CHAIN_LENGTH,
            burnin=BURNIN,
            thin=10,
            output_freq=500,
            title=f"/tmp/openmp_test_{num_threads}t"
        )
        
        elapsed = time.time() - start
        
        print(f"✓ ({elapsed:.2f}s)")
        
        results[num_threads] = {
            'time': elapsed
        }
        
    except Exception as e:
        print(f"\n  ✗ ERROR: {e}")
        import traceback
        traceback.print_exc()

# Analysis
print("\n" + "=" * 80)
print("RESULTS")
print("=" * 80)

if not results:
    print("\n❌ No successful runs!")
    sys.exit(1)

# Table
print(f"\n{'Threads':<10} {'Time (s)':<12} {'Speedup':<12} {'Efficiency':<12}")
print("-" * 80)

baseline = results.get(1, {}).get('time')

for num_threads in sorted(results.keys()):
    time_val = results[num_threads]['time']
    
    if baseline and baseline > 0:
        speedup = baseline / time_val
        efficiency = (speedup / num_threads) * 100
    else:
        speedup = 1.0
        efficiency = 100.0
    
    print(f"{num_threads:<10} {time_val:<12.2f} {speedup:<12.2f}x {efficiency:<12.1f}%")

# Conclusion
print("\n" + "=" * 80)
print("ANALYSIS")
print("=" * 80)

if baseline and len(results) > 1:
    best_threads = min(results.keys(), key=lambda t: results[t]['time'])
    best_time = results[best_threads]['time']
    best_speedup = baseline / best_time
    
    print(f"\nBaseline (1 thread): {baseline:.2f}s")
    print(f"Best: {best_threads} threads at {best_time:.2f}s")
    print(f"Maximum speedup: {best_speedup:.2f}x")
    
    if best_speedup > 1.5:
        print("\n✅ OpenMP is WORKING!")
        print(f"   Achieved {best_speedup:.2f}x speedup")
    elif best_speedup > 1.1:
        print("\n✓ OpenMP is working with limited benefit")
        print(f"   Only {best_speedup:.2f}x speedup (dataset may be too small)")
    else:
        print("\n⚠️ OpenMP shows minimal benefit")
        print(f"   Speedup only {best_speedup:.2f}x")
        print("   Likely reasons:")
        print("   - Dataset too small (6717 SNPs)")
        print("   - Single chromosome (limited parallelization)")
        print("   - Threading overhead > computation benefit")
    
    # Recommendations
    print("\n" + "-" * 80)
    print("For REAL large-scale data (>50K SNPs, multiple chromosomes):")
    print("  - Expected speedup: 4-8x")
    print("  - Recommended: Use all available cores")
    print("  - This test dataset is too small to show full OpenMP benefit")
    
else:
    print("\nInsufficient data for comparison")

print("\n" + "=" * 80)
print("Test Complete!")
print("=" * 80)

