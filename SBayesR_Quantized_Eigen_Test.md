# SBayesR Test: Original vs Quantized Eigen Data

This note tests SBayesR on the tutorial simulated GWAS dataset using:

- Original float eigen LD blocks: `tutorial/data/ldm/block*.eigen.bin`
- Quantized U eigen LD blocks: `.eigen.q8.bin`
- Entropy-compressed quantized U eigen LD blocks: `.eigen.q8e.bin`

The test uses one fixed eigen cutoff so it exercises the eigen reader directly rather than the cutoff-tuning workflow.

## Inputs

Run from the repository root:

```bash
GCTB=./scr/gctb
SUMSTATS=tutorial/data/sim.ma
FLOAT_LDM=tutorial/data/ldm
OUTDIR=tutorial/results/sbayesr_quantized_eigen_test

mkdir -p "${OUTDIR}"
```

Build the executable:

```bash
make -C scr
```

## Create Quantized Eigen Data

Create ordinary 8-bit U-column quantized eigen files:

```bash
"${GCTB}" \
  --quantize-eigen \
  --bits 8 \
  "${FLOAT_LDM}" \
  "${OUTDIR}/ldm_q8"
```

Create entropy-compressed 8-bit U-column quantized eigen files:

```bash
"${GCTB}" \
  --quantize-eigen \
  --bits 8 \
  --entropy \
  "${FLOAT_LDM}" \
  "${OUTDIR}/ldm_q8e"
```

These commands should copy `ldm.info` and `snp.info` and produce files like:

```text
block1563.eigen.q8.bin
block1563.eigen.q8e.bin
```

## Run SBayesR

Use the same seed, MCMC settings, and eigen cutoff for all runs.

```bash
COMMON_ARGS="--gwas-summary ${SUMSTATS} --sbayes R --ldm-eigen-cutoff 0.995 --chain-length 2000 --burn-in 1000 --seed 123"
```

Original float eigen:

```bash
"${GCTB}" ${COMMON_ARGS} \
  --ldm-eigen "${FLOAT_LDM}" \
  --out "${OUTDIR}/sbayesr_float"
```

Quantized q8:

```bash
"${GCTB}" ${COMMON_ARGS} \
  --ldm-eigen-q8 "${OUTDIR}/ldm_q8" \
  --out "${OUTDIR}/sbayesr_q8"
```

Entropy-compressed quantized q8e:

```bash
"${GCTB}" ${COMMON_ARGS} \
  --ldm-eigen-q8e "${OUTDIR}/ldm_q8e" \
  --out "${OUTDIR}/sbayesr_q8e"
```

## Compare Outputs

The q8 and q8e runs should be numerically identical or extremely close because q8e is only lossless zlib compression of the same int8 matrix. The float and q8 runs are expected to be close, not bit-identical, because q8 quantization is lossy.

Run this R snippet from the repository root:

```r
outdir <- "tutorial/results/sbayesr_quantized_eigen_test"

read_snp <- function(prefix) {
  read.table(file.path(outdir, paste0(prefix, ".snpRes")), header = TRUE)
}

read_par <- function(prefix) {
  read.table(file.path(outdir, paste0(prefix, ".parRes")), header = TRUE)
}

float <- read_snp("sbayesr_float")
q8    <- read_snp("sbayesr_q8")
q8e   <- read_snp("sbayesr_q8e")

stopifnot(identical(float$Name, q8$Name))
stopifnot(identical(q8$Name, q8e$Name))

compare_effect <- function(a, b, label) {
  d <- a$A1Effect - b$A1Effect
  data.frame(
    comparison = label,
    correlation = cor(a$A1Effect, b$A1Effect),
    max_abs_diff = max(abs(d)),
    mean_abs_diff = mean(abs(d)),
    rmse = sqrt(mean(d^2))
  )
}

effect_summary <- rbind(
  compare_effect(float, q8, "float_vs_q8"),
  compare_effect(q8, q8e, "q8_vs_q8e")
)

print(effect_summary)

par_float <- read_par("sbayesr_float")
par_q8    <- read_par("sbayesr_q8")
par_q8e   <- read_par("sbayesr_q8e")

print(par_float)
print(par_q8)
print(par_q8e)
```

## Pass Criteria

This is a smoke/regression test rather than a strict statistical equivalence proof.

Expected:

- All three SBayesR runs finish successfully.
- `q8_vs_q8e` has zero or near-zero differences.
- `float_vs_q8` has high effect correlation and small absolute differences.
- Parameter summaries in `.parRes` are broadly consistent between float and q8.

For a faster parser/reader smoke test, reduce `--chain-length` and `--burn-in`, for example:

```bash
--chain-length 20 --burn-in 10
```

For real validation, use longer chains and inspect convergence diagnostics as usual.
