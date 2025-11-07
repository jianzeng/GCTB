#!/bin/bash
# Thorough OpenMP scaling test
# Runs actual MCMC and measures time

echo "================================================================================"
echo "GCTB OpenMP Scaling Test"
echo "================================================================================"
echo ""
echo "Dataset: 6717 SNPs × 3642 individuals"
echo "Chain: 2000 iterations, Burnin: 200"
echo "Model: BayesC"
echo ""

cd python

# Function to run test
run_test() {
    local threads=$1
    echo "----------------------------------------------------------------------"
    echo "Testing with OMP_NUM_THREADS=$threads"
    echo "----------------------------------------------------------------------"
    
    OMP_NUM_THREADS=$threads python3 -m gctb.cli bayes \
        --bfile ../test/data/uk10k_chr1_1mb \
        --pheno ../test/data/test.phen \
        --bayes C \
        --chain-length 2000 \
        --burnin 200 \
        --out /tmp/openmp_${threads}t \
        2>&1 | grep -A 1 "Computational time:"
    
    echo ""
}

echo "Running tests (this will take ~1-2 minutes total)..."
echo ""

run_test 1
run_test 2
run_test 4
run_test 8

cd ..

echo "================================================================================"
echo "Analysis"
echo "================================================================================"
echo ""
echo "Compare the 'Computational time' above for each thread count."
echo ""
echo "Expected results for this small dataset:"
echo "  - 1 thread:  ~6-7 seconds (baseline)"
echo "  - 2 threads: ~4-5 seconds (1.3-1.5x speedup)"
echo "  - 4 threads: ~3-4 seconds (1.5-2x speedup)"
echo "  - 8 threads: ~3-4 seconds (minimal additional benefit)"
echo ""
echo "Why limited speedup?"
echo "  - Dataset is small (only 6717 SNPs)"
echo "  - Single chromosome (limited parallelization opportunities)"
echo "  - Threading overhead significant for small problems"
echo ""
echo "For REAL large-scale data (>50K SNPs, multiple chromosomes):"
echo "  - Expected: 4-8x speedup with 4-8 threads"
echo "  - Near-linear scaling with chromosome count"
echo ""
echo "================================================================================"

