#!/bin/bash
# Run OpenMP scaling tests
# Each thread count in a fresh process

echo "================================================================================"
echo "GCTB OpenMP Scaling Test"
echo "================================================================================"
echo ""
echo "Running MCMC with different thread counts (2000 iterations each)"
echo "Dataset: test/data/uk10k_chr1_1mb (6717 SNPs, 3642 individuals)"
echo ""

cd python

# Test 1 thread
echo "================================================================================"
echo "Test 1: OMP_NUM_THREADS=1"
echo "================================================================================"
time OMP_NUM_THREADS=1 python3 -m gctb.cli bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes C \
    --chain-length 2000 \
    --burnin 200 \
    --quiet \
    --out /tmp/openmp_1t \
    2>&1 | grep -E "(completed|Heritability)"

echo ""

# Test 2 threads
echo "================================================================================"
echo "Test 2: OMP_NUM_THREADS=2"
echo "================================================================================"
time OMP_NUM_THREADS=2 python3 -m gctb.cli bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes C \
    --chain-length 2000 \
    --burnin 200 \
    --quiet \
    --out /tmp/openmp_2t \
    2>&1 | grep -E "(completed|Heritability)"

echo ""

# Test 4 threads
echo "================================================================================"
echo "Test 3: OMP_NUM_THREADS=4"
echo "================================================================================"
time OMP_NUM_THREADS=4 python3 -m gctb.cli bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes C \
    --chain-length 2000 \
    --burnin 200 \
    --quiet \
    --out /tmp/openmp_4t \
    2>&1 | grep -E "(completed|Heritability)"

echo ""

# Test 8 threads
echo "================================================================================"
echo "Test 4: OMP_NUM_THREADS=8"
echo "================================================================================"
time OMP_NUM_THREADS=8 python3 -m gctb.cli bayes \
    --bfile ../test/data/uk10k_chr1_1mb \
    --pheno ../test/data/test.phen \
    --bayes C \
    --chain-length 2000 \
    --burnin 200 \
    --quiet \
    --out /tmp/openmp_8t \
    2>&1 | grep -E "(completed|Heritability)"

cd ..

echo ""
echo "================================================================================"
echo "Summary"
echo "================================================================================"
echo ""
echo "Check the 'real' time for each test above."
echo "If OpenMP is working, you should see speedup with more threads."
echo ""
echo "Expected for this small dataset:"
echo "  - 1.5-2x speedup with 2-4 threads"
echo "  - Limited benefit beyond 4 threads (dataset too small)"
echo ""
echo "For large datasets (>50K SNPs), expect 4-8x speedup!"
echo "================================================================================"

