#!/bin/bash
# Compare Python vs C++ GCTB outputs
# 
# This script runs the same analysis with both implementations
# and compares the results.

set -e

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║     GCTB Comparison: Python vs C++ Implementation            ║"
echo "╚═══════════════════════════════════════════════════════════════╝"
echo ""

# Configuration
BFILE="test/data/uk10k_chr1_1mb"
PHENO="test/data/test.phen"
CHAIN=1100
BURNIN=100
SEED=12345

# Output directories
CPP_OUT="/tmp/gctb_cpp_test"
PY_OUT="/tmp/gctb_py_test"

echo "Configuration:"
echo "  Data: $BFILE"
echo "  Phenotype: $PHENO"
echo "  Chain length: $CHAIN"
echo "  Burn-in: $BURNIN"
echo "  Seed: $SEED"
echo ""

# ============================================================================
# Test 1: C++ GCTB
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1️⃣  Running C++ GCTB (Original Implementation)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ ! -f "scr/gctb" ]; then
    echo "❌ C++ GCTB not compiled!"
    echo "   Compile with: cd scr && make"
    exit 1
fi

TIME_START=$(date +%s)

./scr/gctb \
    --bfile $BFILE \
    --pheno $PHENO \
    --bayes S \
    --chain-length $CHAIN \
    --burn-in $BURNIN \
    --seed $SEED \
    --out $CPP_OUT \
    > /tmp/cpp_output.log 2>&1

TIME_END=$(date +%s)
CPP_TIME=$((TIME_END - TIME_START))

echo "✅ C++ GCTB completed in ${CPP_TIME}s"
echo ""

# Extract key results
if [ -f "${CPP_OUT}.parRes" ]; then
    echo "C++ Results:"
    head -2 ${CPP_OUT}.parRes | tail -1 | awk '{
        print "  Heritability (h²): " $1
        print "  Pi: " $2  
        print "  Genetic Var: " $5
        print "  Residual Var: " $6
    }'
    echo ""
fi

# ============================================================================
# Test 2: Python GCTB
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2️⃣  Running Python GCTB (New Implementation)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

cd python
source ../venv/bin/activate

TIME_START=$(date +%s)

python3 -c "import gctb.cli; gctb.cli.main()" bayes \
    --bfile ../$BFILE \
    --pheno ../$PHENO \
    --bayes S \
    --chain-length $CHAIN \
    --burnin $BURNIN \
    --out $PY_OUT \
    --quiet \
    > /tmp/python_output.log 2>&1

TIME_END=$(date +%s)
PY_TIME=$((TIME_END - TIME_START))

echo "✅ Python GCTB completed in ${PY_TIME}s"
echo ""

cd ..

# Extract key results
if [ -f "${PY_OUT}.parRes" ]; then
    echo "Python Results:"
    tail -1 ${PY_OUT}.parRes | awk '{
        print "  Heritability (h²): " $1
        print "  Pi: " $2
        print "  Genetic Var: " $5  
        print "  Residual Var: " $6
    }'
    echo ""
fi

# ============================================================================
# Comparison
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3️⃣  Comparison Results"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "⏱️  Execution Time:"
echo "   C++:    ${CPP_TIME}s"
echo "   Python: ${PY_TIME}s"

if [ $PY_TIME -lt $CPP_TIME ]; then
    SPEEDUP=$(echo "scale=2; $CPP_TIME / $PY_TIME" | bc)
    echo "   → Python is ${SPEEDUP}x faster! ✅"
elif [ $PY_TIME -eq $CPP_TIME ]; then
    echo "   → Identical speed! ✅"
else
    SLOWDOWN=$(echo "scale=2; $PY_TIME / $CPP_TIME" | bc)
    echo "   → C++ is ${SLOWDOWN}x faster"
fi

echo ""
echo "📊 Output Files Comparison:"
echo ""

# Compare file sizes
if [ -f "${CPP_OUT}.snpRes" ] && [ -f "${PY_OUT}.snpRes" ]; then
    CPP_SIZE=$(wc -l < ${CPP_OUT}.snpRes)
    PY_SIZE=$(wc -l < ${PY_OUT}.snpRes)
    echo "   SNP Results (.snpRes):"
    echo "     C++:    $CPP_SIZE lines"
    echo "     Python: $PY_SIZE lines"
    
    if [ $CPP_SIZE -eq $PY_SIZE ]; then
        echo "     → Same number of SNPs ✅"
    else
        echo "     → Different number of SNPs ⚠️"
    fi
fi

echo ""
echo "📁 Output Files Location:"
echo "   C++:    ${CPP_OUT}.*"
echo "   Python: ${PY_OUT}.*"
echo ""

echo "To compare parameter estimates in detail:"
echo "  diff ${CPP_OUT}.parRes ${PY_OUT}.parRes"
echo ""
echo "To compare SNP results:"
echo "  diff <(head -20 ${CPP_OUT}.snpRes) <(head -20 ${PY_OUT}.snpRes)"
echo ""

echo "╔═══════════════════════════════════════════════════════════════╗"
echo "║                    BENCHMARK COMPLETE                         ║"
echo "╚═══════════════════════════════════════════════════════════════╝"

