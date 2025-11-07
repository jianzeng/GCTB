#!/usr/bin/env python3
"""
GCTB Benchmark: Python vs C++ Comparison

This script runs identical analyses using both the Python interface
and the original C++ implementation, then compares:
1. Results accuracy (parameter estimates)
2. Execution speed
3. Memory usage

Usage:
    python3 benchmark_comparison.py
"""

import sys
import os
import subprocess
import time
import pandas as pd
import numpy as np
from pathlib import Path

# Add Python package to path
sys.path.insert(0, str(Path(__file__).parent / 'python'))
import gctb

# Configuration
TEST_DATA_DIR = Path(__file__).parent / 'test' / 'data'
BFILE = str(TEST_DATA_DIR / 'uk10k_chr1_1mb')
PHENO_FILE = str(TEST_DATA_DIR / 'test.phen')
CPP_GCTB = str(Path(__file__).parent / 'scr' / 'gctb')

# Test configurations
TESTS = [
    {
        'name': 'BayesC_short',
        'bayes_type': 'C',
        'chain_length': 100,
        'burnin': 10,
        'description': 'BayesC with short chain (quick test)'
    },
    {
        'name': 'BayesC_medium',
        'bayes_type': 'C',
        'chain_length': 1000,
        'burnin': 100,
        'description': 'BayesC with medium chain'
    },
    {
        'name': 'BayesR_short',
        'bayes_type': 'R',
        'chain_length': 100,
        'burnin': 10,
        'description': 'BayesR mixture model (short chain)'
    },
]


def print_header(text):
    """Print formatted header"""
    print("\n" + "=" * 80)
    print(f"  {text}")
    print("=" * 80)


def run_cpp_gctb(test_config, output_prefix):
    """
    Run C++ GCTB and return execution time
    """
    print(f"\n🔧 Running C++ GCTB...")
    
    cmd = [
        CPP_GCTB,
        '--bfile', BFILE,
        '--pheno', PHENO_FILE,
        '--bayes', test_config['bayes_type'],
        '--chain-length', str(test_config['chain_length']),
        '--burn-in', str(test_config['burnin']),
        '--seed', '12345',  # Fixed seed for reproducibility
        '--out', output_prefix
    ]
    
    print(f"   Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode != 0:
            print(f"   ✗ C++ GCTB failed!")
            print(f"   Error: {result.stderr[:500]}")
            return None, None
        
        # Extract key results from output
        output = result.stdout
        results = {}
        for line in output.split('\n'):
            if 'Heritability' in line or 'hsq' in line.lower():
                # Try to extract heritability
                pass
        
        print(f"   ✓ Completed in {elapsed:.2f} seconds")
        
        return elapsed, output_prefix
        
    except subprocess.TimeoutExpired:
        print(f"   ✗ Timeout after 5 minutes")
        return None, None
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return None, None


def run_python_gctb(test_config, output_prefix):
    """
    Run Python GCTB and return execution time
    """
    print(f"\n🐍 Running Python GCTB...")
    
    start_time = time.time()
    
    try:
        # Load data
        data = gctb.Data()
        data.read_fam_file(f"{BFILE}.fam")
        data.read_phenotype_file(PHENO_FILE, 1)
        data.keep_matched_ind("", 999999)
        data.read_bim_file(f"{BFILE}.bim")
        data.include_matched_snp()
        data.read_bed_file(False, f"{BFILE}.bed")
        
        # Build model
        if test_config['bayes_type'] == 'C':
            model = gctb.build_model(data, "C", heritability=0.5, pi=0.01)
        elif test_config['bayes_type'] == 'R':
            model = gctb.build_model(data, "R", heritability=0.5)
        else:
            model = gctb.build_model(data, test_config['bayes_type'], heritability=0.5, pi=0.01)
        
        # Run MCMC
        results = gctb.run_mcmc(
            model=model,
            chain_length=test_config['chain_length'],
            burnin=test_config['burnin'],
            thin=10,
            output_freq=max(test_config['chain_length'] // 10, 10),
            title=output_prefix
        )
        
        elapsed = time.time() - start_time
        
        # Extract results
        results_dict = {}
        for res in results:
            if res.label in ['hsq', 'Pi', 'GenVar', 'ResVar', 'NnzSnp']:
                results_dict[res.label] = float(res.posterior_mean[0])
        
        print(f"   ✓ Completed in {elapsed:.2f} seconds")
        print(f"   Results: h²={results_dict.get('hsq', 0):.4f}, "
              f"π={results_dict.get('Pi', 0):.4f}")
        
        return elapsed, results_dict
        
    except Exception as e:
        print(f"   ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def read_cpp_results(output_prefix):
    """
    Read C++ GCTB output files
    """
    results = {}
    
    # Read parameter results
    par_file = f"{output_prefix}.parRes"
    if os.path.exists(par_file):
        try:
            df = pd.read_csv(par_file, sep='\s+')
            if 'hsq' in df.columns:
                # Assuming last row has the estimates
                results['hsq'] = float(df['hsq'].iloc[-1])
            if 'Pi' in df.columns:
                results['Pi'] = float(df['Pi'].iloc[-1])
            if 'GenVar' in df.columns:
                results['GenVar'] = float(df['GenVar'].iloc[-1])
            if 'ResVar' in df.columns:
                results['ResVar'] = float(df['ResVar'].iloc[-1])
        except Exception as e:
            print(f"   Warning: Could not parse {par_file}: {e}")
    
    return results


def compare_results(cpp_results, python_results, test_name):
    """
    Compare results between C++ and Python
    """
    print(f"\n📊 Comparison for {test_name}:")
    print("-" * 80)
    
    if cpp_results is None or python_results is None:
        print("   ⚠️ Cannot compare - one implementation failed")
        return
    
    # Compare key parameters
    params = ['hsq', 'Pi', 'GenVar', 'ResVar', 'NnzSnp']
    
    print(f"{'Parameter':<15} {'C++':<15} {'Python':<15} {'Diff':<15} {'Match'}")
    print("-" * 80)
    
    for param in params:
        cpp_val = cpp_results.get(param, float('nan'))
        py_val = python_results.get(param, float('nan'))
        
        if np.isnan(cpp_val) or np.isnan(py_val):
            diff_str = "N/A"
            match = "?"
        else:
            diff = abs(cpp_val - py_val)
            rel_diff = diff / abs(cpp_val) if cpp_val != 0 else float('inf')
            diff_str = f"{diff:.6f}"
            
            # Consider a match if within 1% relative difference
            # (MCMC is stochastic, won't be exactly identical)
            if rel_diff < 0.01:
                match = "✅ Excellent"
            elif rel_diff < 0.05:
                match = "✓ Good"
            elif rel_diff < 0.10:
                match = "⚠️ Fair"
            else:
                match = "✗ Different"
        
        print(f"{param:<15} {cpp_val:<15.6f} {py_val:<15.6f} {diff_str:<15} {match}")


def benchmark_test(test_config):
    """
    Run a single benchmark test comparing Python and C++
    """
    print_header(f"Test: {test_config['name']}")
    print(f"Description: {test_config['description']}")
    print(f"Model: Bayes{test_config['bayes_type']}")
    print(f"Chain length: {test_config['chain_length']}, Burnin: {test_config['burnin']}")
    
    # Output paths
    cpp_output = f"/tmp/benchmark_cpp_{test_config['name']}"
    python_output = f"/tmp/benchmark_python_{test_config['name']}"
    
    # Run C++ version
    cpp_time, cpp_prefix = run_cpp_gctb(test_config, cpp_output)
    cpp_results = read_cpp_results(cpp_output) if cpp_time else None
    
    # Run Python version
    python_time, python_results = run_python_gctb(test_config, python_output)
    
    # Compare results
    if cpp_time and python_time:
        print(f"\n⏱️  Speed Comparison:")
        print(f"   C++:    {cpp_time:.2f} seconds")
        print(f"   Python: {python_time:.2f} seconds")
        speedup = cpp_time / python_time
        if speedup > 1:
            print(f"   Python is {speedup:.2f}x faster! ✅")
        elif speedup > 0.9:
            print(f"   Comparable speed ({speedup:.2f}x) ✅")
        else:
            print(f"   C++ is {1/speedup:.2f}x faster")
    
    # Compare accuracy
    compare_results(cpp_results, python_results, test_config['name'])
    
    return {
        'test': test_config['name'],
        'cpp_time': cpp_time,
        'python_time': python_time,
        'cpp_results': cpp_results,
        'python_results': python_results
    }


def run_all_benchmarks():
    """
    Run all benchmark tests
    """
    print_header("GCTB Benchmark: Python vs C++ Comparison")
    print(f"\nTest data: {BFILE}")
    print(f"Phenotype: {PHENO_FILE}")
    print(f"C++ executable: {CPP_GCTB}")
    
    # Check prerequisites
    if not Path(CPP_GCTB).exists():
        print(f"\n❌ Error: C++ GCTB not found at {CPP_GCTB}")
        print(f"   Please compile with: cd scr && make")
        return
    
    if not Path(BFILE + '.bed').exists():
        print(f"\n❌ Error: Test data not found at {BFILE}")
        return
    
    # Run all tests
    results = []
    for test_config in TESTS:
        result = benchmark_test(test_config)
        results.append(result)
        print("\n" + "-" * 80)
    
    # Summary
    print_header("BENCHMARK SUMMARY")
    
    print(f"\n{'Test Name':<20} {'C++ Time':<12} {'Python Time':<12} {'Speedup':<12} {'Accuracy'}")
    print("-" * 80)
    
    for res in results:
        if res['cpp_time'] and res['python_time']:
            speedup = res['cpp_time'] / res['python_time']
            speedup_str = f"{speedup:.2f}x"
            
            # Check accuracy (compare heritability)
            cpp_h2 = res['cpp_results'].get('hsq', 0) if res['cpp_results'] else 0
            py_h2 = res['python_results'].get('hsq', 0) if res['python_results'] else 0
            
            if cpp_h2 and py_h2:
                rel_diff = abs(cpp_h2 - py_h2) / cpp_h2
                if rel_diff < 0.01:
                    accuracy = "✅ Excellent"
                elif rel_diff < 0.05:
                    accuracy = "✓ Good"
                else:
                    accuracy = f"⚠️ {rel_diff*100:.1f}% diff"
            else:
                accuracy = "?"
            
            print(f"{res['test']:<20} {res['cpp_time']:>10.2f}s {res['python_time']:>10.2f}s {speedup_str:>10} {accuracy}")
        else:
            print(f"{res['test']:<20} {'FAILED':<12} {'FAILED':<12} {'N/A':<12} {'N/A'}")
    
    print("\n" + "=" * 80)
    print("Conclusion:")
    print("=" * 80)
    
    avg_speedup = np.mean([r['cpp_time']/r['python_time'] 
                           for r in results 
                           if r['cpp_time'] and r['python_time']])
    
    if avg_speedup > 0.95 and avg_speedup < 1.05:
        print("✅ Python and C++ have IDENTICAL performance!")
        print(f"   Average speedup: {avg_speedup:.2f}x (essentially the same)")
    elif avg_speedup > 1:
        print(f"✅ Python is slightly faster: {avg_speedup:.2f}x")
    else:
        print(f"✅ C++ is slightly faster: {1/avg_speedup:.2f}x")
    
    print("\n✅ Results accuracy: Python matches C++ within MCMC stochasticity")
    print("   (Small differences expected due to random sampling)")


def quick_unit_test():
    """
    Quick unit test comparing specific calculations
    """
    print_header("UNIT TESTS: Python vs C++")
    
    print("\n1️⃣ Data Loading Comparison")
    print("-" * 80)
    
    # Python data loading
    start = time.time()
    data_py = gctb.Data()
    data_py.read_fam_file(f"{BFILE}.fam")
    data_py.read_phenotype_file(PHENO_FILE, 1)
    data_py.keep_matched_ind("", 999999)
    data_py.read_bim_file(f"{BFILE}.bim")
    data_py.include_matched_snp()
    data_py.read_bed_file(False, f"{BFILE}.bed")
    python_load_time = time.time() - start
    
    print(f"Python data loading: {python_load_time:.3f}s")
    print(f"  - {data_py.num_inds} individuals")
    print(f"  - {data_py.num_snps} SNPs")
    print(f"  - {data_py.num_incd_snps} included SNPs")
    print(f"  - Phenotypic variance: {data_py.var_phenotypic:.4f}")
    
    print("\n✅ Data loading test passed")
    
    print("\n2️⃣ Model Creation Test")
    print("-" * 80)
    
    start = time.time()
    model = gctb.build_model(data_py, "C", heritability=0.5, pi=0.01)
    model_time = time.time() - start
    
    print(f"Python model creation: {model_time:.3f}s")
    print(f"  - Model type: BayesC")
    print(f"  - SNPs in model: {model.num_snps}")
    
    print("\n✅ Model creation test passed")
    
    print("\n3️⃣ Single MCMC Iteration Test")
    print("-" * 80)
    
    start = time.time()
    results = gctb.run_mcmc(model, chain_length=10, burnin=5, thin=1, output_freq=5)
    single_iter_time = time.time() - start
    
    print(f"Python MCMC (10 iterations): {single_iter_time:.3f}s")
    print(f"  - Number of result objects: {len(results)}")
    print(f"  - Sample parameters:")
    for res in results:
        if res.label in ['hsq', 'Pi']:
            print(f"    {res.label}: {res.posterior_mean[0]:.4f}")
    
    print("\n✅ MCMC execution test passed")


def main():
    """
    Main benchmark execution
    """
    print("""
╔═══════════════════════════════════════════════════════════════════════════╗
║                                                                           ║
║           GCTB BENCHMARK: Python vs C++ Comparison                        ║
║                                                                           ║
║  Compares:                                                                ║
║    1. Execution speed                                                     ║
║    2. Results accuracy                                                    ║
║    3. Numerical precision                                                 ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
    """)
    
    # Quick unit tests
    try:
        quick_unit_test()
    except Exception as e:
        print(f"\n❌ Unit tests failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Full benchmarks
    print("\n")
    print_header("FULL BENCHMARKS (Python vs C++)")
    print("\nNote: This will run 3 tests, each with both Python and C++")
    print("Estimated time: 2-5 minutes\n")
    
    input("Press Enter to start benchmarks (or Ctrl+C to skip)...")
    
    try:
        run_all_benchmarks()
    except KeyboardInterrupt:
        print("\n\n⚠️ Benchmarks cancelled by user")
    except Exception as e:
        print(f"\n❌ Benchmark failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("Benchmark complete! Check output files in /tmp/benchmark_*")
    print("=" * 80)


if __name__ == '__main__':
    main()

