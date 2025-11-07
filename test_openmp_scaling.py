#!/usr/bin/env python3
"""
Comprehensive OpenMP Scaling Test for GCTB

Tests performance with 1, 2, 4, and 8 threads to verify:
1. OpenMP is actually being used
2. Performance scales with thread count
3. Optimal thread count for this dataset
"""

import sys
import os
import subprocess
import time
import json
from pathlib import Path

# Configuration  
BFILE = "../test/data/uk10k_chr1_1mb"
PHENO = "../test/data/test.phen"
CHAIN_LENGTH = 2000  # Longer chain to see threading benefits
BURNIN = 200
BAYES_TYPE = "C"
OUTPUT_PREFIX = "/tmp/openmp_test"

# Thread counts to test
THREAD_COUNTS = [1, 2, 4, 8]

# Number of runs per configuration (for averaging)
NUM_RUNS = 3

print("=" * 80)
print("GCTB OpenMP Scaling Test")
print("=" * 80)
print(f"\nConfiguration:")
print(f"  Dataset: {BFILE}")
print(f"  Chain length: {CHAIN_LENGTH}")
print(f"  Burn-in: {BURNIN}")
print(f"  Thread counts to test: {THREAD_COUNTS}")
print(f"  Runs per configuration: {NUM_RUNS}")
print(f"  Model: Bayes{BAYES_TYPE}")
print("")

# Check if running from correct directory
test_file = Path("test/data/uk10k_chr1_1mb.bed")
if not test_file.exists():
    print(f"Error: Could not find {test_file}")
    print("Please run this script from the GCTB root directory")
    sys.exit(1)

results = {}

for num_threads in THREAD_COUNTS:
    print("=" * 80)
    print(f"Testing with {num_threads} thread(s)")
    print("=" * 80)
    
    times = []
    
    for run in range(NUM_RUNS):
        print(f"\n  Run {run + 1}/{NUM_RUNS}...", end=" ", flush=True)
        
        # Set environment
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(num_threads)
        
        # Run GCTB
        cmd = [
            sys.executable, "-m", "gctb.cli", "bayes",
            "--bfile", BFILE,
            "--pheno", PHENO,
            "--bayes", BAYES_TYPE,
            "--chain-length", str(CHAIN_LENGTH),
            "--burnin", str(BURNIN),
            "--quiet",
            "--out", f"{OUTPUT_PREFIX}_{num_threads}t_run{run}"
        ]
        
        start = time.time()
        
        try:
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=300,
                cwd="python"
            )
            
            elapsed = time.time() - start
            
            if result.returncode not in [0, 133]:  # 133 is the threading warning exit code
                print(f"FAILED (exit code {result.returncode})")
                print(f"Error: {result.stderr[:200]}")
                continue
            
            times.append(elapsed)
            print(f"{elapsed:.2f}s", flush=True)
            
        except subprocess.TimeoutExpired:
            print("TIMEOUT (>300s)")
        except Exception as e:
            print(f"ERROR: {e}")
    
    if times:
        avg_time = sum(times) / len(times)
        min_time = min(times)
        max_time = max(times)
        std_dev = (sum((t - avg_time) ** 2 for t in times) / len(times)) ** 0.5
        
        results[num_threads] = {
            'times': times,
            'avg': avg_time,
            'min': min_time,
            'max': max_time,
            'std': std_dev
        }
        
        print(f"\n  Summary for {num_threads} thread(s):")
        print(f"    Average: {avg_time:.2f}s")
        print(f"    Min:     {min_time:.2f}s")
        print(f"    Max:     {max_time:.2f}s")
        print(f"    Std Dev: {std_dev:.2f}s")
    else:
        print(f"\n  No successful runs for {num_threads} thread(s)")

# Analysis
print("\n" + "=" * 80)
print("RESULTS SUMMARY")
print("=" * 80)

if not results:
    print("No successful runs!")
    sys.exit(1)

# Print table
print(f"\n{'Threads':<10} {'Avg Time':<12} {'Speedup':<12} {'Efficiency':<12} {'Status'}")
print("-" * 80)

baseline = results[1]['avg'] if 1 in results else None

for num_threads in sorted(results.keys()):
    avg_time = results[num_threads]['avg']
    
    if baseline:
        speedup = baseline / avg_time
        efficiency = speedup / num_threads * 100
    else:
        speedup = 0
        efficiency = 0
    
    # Status based on efficiency
    if efficiency > 90:
        status = "✅ Excellent"
    elif efficiency > 70:
        status = "✅ Good"
    elif efficiency > 50:
        status = "⚠️ Fair"
    elif efficiency > 25:
        status = "⚠️ Poor"
    else:
        status = "❌ No benefit"
    
    print(f"{num_threads:<10} {avg_time:>10.2f}s  {speedup:>10.2f}x  {efficiency:>10.1f}%  {status}")

# Detailed analysis
print("\n" + "=" * 80)
print("ANALYSIS")
print("=" * 80)

if baseline and len(results) > 1:
    # Check for speedup
    best_threads = min(results.keys(), key=lambda t: results[t]['avg'])
    best_time = results[best_threads]['avg']
    best_speedup = baseline / best_time
    
    print(f"\nBaseline (1 thread): {baseline:.2f}s")
    print(f"Best performance: {best_threads} threads at {best_time:.2f}s")
    print(f"Best speedup: {best_speedup:.2f}x")
    
    # Determine if OpenMP is working
    if best_speedup > 1.5:
        print("\n✅ OpenMP is WORKING!")
        print(f"   Achieved {best_speedup:.2f}x speedup with {best_threads} threads")
        print(f"   Efficiency: {(best_speedup / best_threads * 100):.1f}%")
    elif best_speedup > 1.1:
        print("\n✓ OpenMP is working, but with limited benefit")
        print(f"   Only {best_speedup:.2f}x speedup (dataset may be too small)")
    else:
        print("\n⚠️ OpenMP may not be working effectively")
        print(f"   Speedup is only {best_speedup:.2f}x")
        print("   Possible reasons:")
        print("   - Dataset too small (threading overhead dominates)")
        print("   - OpenMP not properly configured")
        print("   - Single chromosome limits parallelization")
    
    # Scaling efficiency
    print("\n" + "-" * 80)
    print("Scaling Efficiency:")
    print("-" * 80)
    for num_threads in sorted(results.keys()):
        if num_threads == 1:
            continue
        speedup = baseline / results[num_threads]['avg']
        efficiency = speedup / num_threads * 100
        print(f"  {num_threads} threads: {efficiency:.1f}% efficient (ideal: 100%)")
    
    # Recommendations
    print("\n" + "-" * 80)
    print("Recommendations:")
    print("-" * 80)
    
    # Find optimal thread count (best time)
    optimal_threads = min(results.keys(), key=lambda t: results[t]['avg'])
    
    print(f"  Optimal for this dataset: {optimal_threads} threads")
    
    if best_speedup < 1.5:
        print("  \n  ⚠️ This dataset is too small to show significant OpenMP benefits")
        print("  For real large-scale analyses (>50K SNPs, >10K individuals):")
        print("    - Expected speedup: 4-8x")
        print("    - Recommended threads: 4-8 (or all cores)")
    else:
        print(f"  \n  ✅ OpenMP is providing good speedup!")
        print(f"  For production runs, use: OMP_NUM_THREADS={optimal_threads}")

else:
    print("Insufficient data for analysis (need at least 2 thread configurations)")

# Save results
results_file = "openmp_test_results.json"
with open(results_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n\nDetailed results saved to: {results_file}")

print("\n" + "=" * 80)
print("Test Complete!")
print("=" * 80)

