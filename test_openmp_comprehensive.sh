#!/bin/bash
# Comprehensive OpenMP test with venv

echo "================================================================================"
echo "GCTB OpenMP Scaling Test - Comprehensive"
echo "================================================================================"
echo ""

# Activate venv
source venv/bin/activate

cd python

echo "Testing 1 thread..."
OMP_NUM_THREADS=1 bash -c 'time python3 -m gctb.cli bayes --bfile ../test/data/uk10k_chr1_1mb --pheno ../test/data/test.phen --bayes C --chain-length 2000 --burnin 200 --out /tmp/omp1 > /tmp/omp1.log 2>&1' 2>&1 | grep "real"

echo ""
echo "Testing 2 threads..."
OMP_NUM_THREADS=2 bash -c 'time python3 -m gctb.cli bayes --bfile ../test/data/uk10k_chr1_1mb --pheno ../test/data/test.phen --bayes C --chain-length 2000 --burnin 200 --out /tmp/omp2 > /tmp/omp2.log 2>&1' 2>&1 | grep "real"

echo ""
echo "Testing 4 threads..."
OMP_NUM_THREADS=4 bash -c 'time python3 -m gctb.cli bayes --bfile ../test/data/uk10k_chr1_1mb --pheno ../test/data/test.phen --bayes C --chain-length 2000 --burnin 200 --out /tmp/omp4 > /tmp/omp4.log 2>&1' 2>&1 | grep "real"

echo ""
echo "Testing 8 threads..."
OMP_NUM_THREADS=8 bash -c 'time python3 -m gctb.cli bayes --bfile ../test/data/uk10k_chr1_1mb --pheno ../test/data/test.phen --bayes C --chain-length 2000 --burnin 200 --out /tmp/omp8 > /tmp/omp8.log 2>&1' 2>&1 | grep "real"

echo ""
echo "================================================================================"
echo "Results Summary"
echo "================================================================================"
echo "Compare the 'real' times above."
echo ""
echo "If OpenMP is working:"
echo "  - 2 threads should be ~1.3-1.5x faster than 1 thread"
echo "  - 4 threads should be ~1.5-2x faster than 1 thread"
echo ""
echo "For this small dataset, expect LIMITED speedup due to:"
echo "  - Small size (6717 SNPs)"
echo "  - Single chromosome"
echo "  - Threading overhead"
echo ""
echo "Large datasets (>50K SNPs) will show 4-8x speedup!"
echo "================================================================================"
