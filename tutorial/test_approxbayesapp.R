#!/usr/bin/env Rscript
# Quick test script for ApproxBayesAPP

# Set paths
data_dir <- "data"
result_dir <- "results"
gctb <- "../scr/gctb"

if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

# Load data
bed_file <- file.path(data_dir, "1000G_eur_chr22")
fam <- read.table(paste0(bed_file, ".fam"), stringsAsFactors = FALSE,
                  col.names = c("FID", "IID", "FATHER", "MOTHER", "SEX", "PHENO"))
bim_all <- read.table(paste0(bed_file, ".bim"), stringsAsFactors = FALSE,
                      col.names = c("CHR", "SNP", "CM", "BP", "ALLELE_A", "ALLELE_B"))

# Simulation parameters
h2_1      <- 0.99   # heritability trait 1
pi_1      <- 0.001  # polygenicity trait 1 (fraction of causal SNPs)
h2_2      <- 0.99   # heritability trait 2
pi_2      <- 0.001  # polygenicity trait 2 (fraction of causal SNPs)

sim_files <- c(file.path(data_dir, "sim.ma"), file.path(data_dir, "sim_trait2.ma"))
if (!all(file.exists(sim_files))) {
  if (!require("BEDMatrix", quietly = TRUE)) install.packages("BEDMatrix")
  library(BEDMatrix)

  bed  <- BEDMatrix(bed_file)
  geno <- as.matrix(bed)
  geno[is.na(geno)] <- 0   # impute missing as 0

  n_ind <- nrow(geno)

  # Filter SNPs with MAF < 0.01
  maf_all <- colMeans(geno) / 2
  maf_all <- pmin(maf_all, 1 - maf_all)
  keep_snp <- maf_all >= 0.01
  geno <- geno[, keep_snp]
  bim  <- bim_all[keep_snp, ]
  n_snp <- ncol(geno)
  af    <- colMeans(geno) / 2

  cat("Loaded:", n_ind, "individuals,", n_snp, "SNPs after filtering\n")

  # Helper: fast vectorised marginal GWAS
  gwas_fast <- function(geno, y) {
    geno_c   <- sweep(geno, 2, colMeans(geno))   # mean-centre
    XtX_diag <- colSums(geno_c^2)
    Xty      <- drop(crossprod(geno_c, y - mean(y)))
    beta     <- Xty / XtX_diag
    yty      <- sum((y - mean(y))^2)
    rss      <- pmax(yty - beta^2 * XtX_diag, 0)
    se       <- sqrt(rss / ((length(y) - 2) * XtX_diag))
    tstat    <- beta / se
    pval     <- 2 * pt(-abs(tstat), df = length(y) - 2)
    list(beta = beta, se = se, pval = pval)
  }

  set.seed(42)

  # ---- Trait 1 ----
  n_causal_1  <- round(pi_1 * n_snp)
  causal_1    <- sample(n_snp, n_causal_1)
  beta_true_1 <- rnorm(n_causal_1) / sqrt(n_causal_1)
  g_1 <- drop(geno[, causal_1] %*% beta_true_1)
  g_1 <- g_1 * sqrt(h2_1 / var(g_1))
  e_1 <- rnorm(n_ind, 0, sqrt(1 - h2_1))
  y_1 <- g_1 + e_1
  cat(sprintf("Trait 1: n_causal = %d (pi = %.2f), h2 = %.2f, actual h2 = %.3f\n",
              n_causal_1, pi_1, h2_1, var(g_1) / var(y_1)))

  gw1 <- gwas_fast(geno, y_1)
  ss1 <- data.frame(SNP = bim$SNP, A1 = bim$ALLELE_B, A2 = bim$ALLELE_A,
                    FREQ = af, BETA = gw1$beta, SE = gw1$se, P = gw1$pval, N = n_ind)
  write.table(ss1, file.path(data_dir, "sim.ma"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("Generated trait 1 summary statistics\n")

  # ---- Trait 2 ----
  n_causal_2  <- round(pi_2 * n_snp)
  causal_2    <- sample(n_snp, n_causal_2)
  beta_true_2 <- rnorm(n_causal_2) / sqrt(n_causal_2)
  g_2 <- drop(geno[, causal_2] %*% beta_true_2)
  g_2 <- g_2 * sqrt(h2_2 / var(g_2))
  e_2 <- rnorm(n_ind, 0, sqrt(1 - h2_2))
  y_2 <- g_2 + e_2
  cat(sprintf("Trait 2: n_causal = %d (pi = %.2f), h2 = %.2f, actual h2 = %.3f\n",
              n_causal_2, pi_2, h2_2, var(g_2) / var(y_2)))

  gw2 <- gwas_fast(geno, y_2)
  ss2 <- data.frame(SNP = bim$SNP, A1 = bim$ALLELE_B, A2 = bim$ALLELE_A,
                    FREQ = af, BETA = gw2$beta, SE = gw2$se, P = gw2$pval, N = n_ind)
  write.table(ss2, file.path(data_dir, "sim_trait2.ma"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("Generated trait 2 summary statistics\n")
} else {
  cat("Using existing summary statistics (delete sim.ma / sim_trait2.ma to re-simulate)\n")
}

# Check if eigen LD matrices exist
ldm_path <- file.path(data_dir, "ldm")
eigen_files <- list.files(ldm_path, pattern = "\\.eigen\\.bin$", full.names = FALSE)
if (length(eigen_files) == 0) {
  cat("Eigen LD matrices not found. Generating...\n")
  cat("Please run: gctb --ldm data/ldm --gwas-summary data/sim.ma --make-ldm-eigen --out data/ldm\n")
  stop("Eigen LD matrices not found")
} else {
  cat("Found", length(eigen_files), "eigen LD matrix files\n")
}

# Run ApproxBayesAPP
# DEBUG: use the same summary stats for both traits to verify symmetric output
debug_same_trait <- FALSE

cat("\nRunning ApproxBayesAPP...\n")
sumstats_path    <- file.path(data_dir, "sim.ma")
sumstats_path_2  <- if (debug_same_trait) sumstats_path else file.path(data_dir, "sim_trait2.ma")
if (debug_same_trait) cat("[DEBUG] Using same summary stats for both traits\n")
bivariate_sumstats <- paste0(sumstats_path, ",", sumstats_path_2)
out_prefix       <- file.path(result_dir, "approxbayesapp")
log_file         <- paste0(out_prefix, ".log")

gctb_args <- c(
  "--gwas-summary", bivariate_sumstats,
  "--ldm-eigen",    ldm_path,
  "--sbayes",       "APP",
  "--chain-length", "200",
  "--burn-in",      "100",
  "--out",          out_prefix
)

# Capture stdout+stderr; tee to console and log file simultaneously
log_lines <- character(0)
con <- pipe(paste(c(gctb, gctb_args, "2>&1"), collapse = " "), open = "r")
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  cat(line, "\n", sep = "")
  log_lines <- c(log_lines, line)
}
close(con)
writeLines(log_lines, log_file)
cat("\nLog saved to:", log_file, "\n")

# Check results
if (file.exists(paste0(out_prefix, ".snpRes"))) {
  cat("\n✓ ApproxBayesAPP completed successfully!\n")

  snp_res  <- read.table(paste0(out_prefix, ".snpRes"),  header = TRUE)
  par_res  <- read.table(paste0(out_prefix, ".parRes"),  header = TRUE)

  cat("\nTop SNPs by T1 PIP:\n")
  snp_t1 <- snp_res[order(-snp_res$PIP_T1), ]
  print(head(snp_t1[snp_t1$PIP_T1 > 0, c("Name","A1Frq_T1","A1Effect_T1","SE_T1","VarExplained_T1","PIP_T1")], 10))

  cat("\nTop SNPs by T2 PIP:\n")
  snp_t2 <- snp_res[order(-snp_res$PIP_T2), ]
  print(head(snp_t2[snp_t2$PIP_T2 > 0, c("Name","A1Frq_T2","A1Effect_T2","SE_T2","VarExplained_T2","PIP_T2")], 10))

  cat("\nParameter results:\n")
  print(par_res)
} else {
  cat("\n✗ ApproxBayesAPP failed - output files not found\n")
}

