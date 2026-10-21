#!/usr/bin/env bash
# Compare RD with its explicit independent-model alias in both sampling modes.
# Usage: bash scripts/test_rd_default.sh GCTB GWAS Q8_LD_DIR ANNOTATIONS
set -euo pipefail
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 GCTB GWAS Q8_LD_DIR ANNOTATIONS" >&2
    exit 2
fi
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/gctb-rd-default.XXXXXX")
echo "Test outputs: $test_dir"
for chains in 1 2; do
    for mode in gibbs full; do
        sampler_args=(--alg Gibbs)
        if [ "$mode" = full ]; then sampler_args+=(--joint-annot); fi
        for model in RD RD_indep; do
            prefix="$test_dir/${model}_${chains}_${mode}"
            "$1" --sbayes "$model" --gwas-summary "$2" --ldm-eigen-q8 "$3" \
                --annot "$4" --ldm-eigen-cutoff 0.995 --num-chains "$chains" \
                --thread "$chains" --chain-length 40 --burn-in 20 --thin 1 \
                --out-freq 20 --seed 123 "${sampler_args[@]}" --out "$prefix" \
                > "$prefix.stdout" 2>&1
            grep -q 'Using independent annotation inclusion indicators' "$prefix.stdout"
            grep -q 'AnnoPIP_p2' "$prefix.parSetRes"
            grep -q 'AnnoPi_p2' "$prefix.parRes"
            if grep -Eq '^[[:space:]]*AnnoPIP[[:space:]]' "$prefix.parSetRes"; then
                echo "Unexpected shared AnnoPIP output" >&2
                exit 1
            fi
        done
        for suffix in parRes parSetRes snpRes; do
            cmp "$test_dir/RD_${chains}_${mode}.$suffix" \
                "$test_dir/RD_indep_${chains}_${mode}.$suffix"
        done
        echo "RD equals RD_indep: chains=$chains, mode=$mode"
    done
done
