#!/usr/bin/env python3
"""
GCTB Command Line Interface

Provides user-friendly commands for Bayesian genomic analysis.
"""

import click
import sys
from pathlib import Path
from . import _core as gctb

__version__ = "1.0.0"


def load_plink_data(bfile: str, pheno: str, mphen: int = 1, 
                    keep_ind_file: str = "", keep_ind_max: int = 999999,
                    verbose: bool = True):
    """
    Load PLINK format data (individual-level).
    
    Parameters
    ----------
    bfile : str
        PLINK binary file prefix
    pheno : str
        Phenotype file path
    mphen : int
        Phenotype column to use (1-indexed)
    keep_ind_file : str
        File with individuals to keep (optional)
    keep_ind_max : int
        Maximum number of individuals to keep
    verbose : bool
        Print progress messages
    """
    data = gctb.Data()
    
    if verbose:
        click.echo("Loading PLINK data...")
    
    data.read_fam_file(f"{bfile}.fam")
    data.read_phenotype_file(pheno, mphen)
    data.keep_matched_ind(keep_ind_file, keep_ind_max)
    data.read_bim_file(f"{bfile}.bim")
    data.include_matched_snp()
    data.read_bed_file(False, f"{bfile}.bed")
    
    if verbose:
        click.echo(f"  ✓ Loaded: {data.num_incd_snps} SNPs × {data.num_kept_inds} individuals")
        click.echo(f"  ✓ Phenotypic variance: {data.var_phenotypic:.4f}")
    
    return data


@click.group()
@click.version_option(version=__version__)
def main():
    """
    GCTB: Genome-wide Complex Trait Bayesian Analysis
    
    Python interface with C++ computational core.
    
    Examples:
    
        # Individual-level analysis
        gctb bayes --bfile data --pheno pheno.txt --bayes C
        
        # Summary statistics analysis  
        gctb sbayes --ldm ldm_dir --gwas summary.txt --sbayes R
    """
    pass


@main.command()
@click.option('--bfile', required=True, type=str,
              help='PLINK binary file prefix (.bed/.bim/.fam)')
@click.option('--pheno', required=True, type=click.Path(exists=True),
              help='Phenotype file')
@click.option('--bayes', type=click.Choice(['C', 'B', 'R', 'S'], case_sensitive=False),
              default='C', help='Bayes model type (default: C)')
@click.option('--chain-length', default=1000, type=int,
              help='MCMC chain length (default: 1000)')
@click.option('--burnin', default=100, type=int,
              help='Burn-in iterations (default: 100)')
@click.option('--thin', default=10, type=int,
              help='Thinning interval (default: 10)')
@click.option('--pi', default=0.01, type=float,
              help='Prior probability of non-zero effect (default: 0.01)')
@click.option('--hsq', default=0.5, type=float,
              help='Heritability (default: 0.5)')
@click.option('--mphen', default=1, type=int,
              help='Phenotype column to use (default: 1)')
@click.option('--out', default='gctb', type=str,
              help='Output prefix (default: gctb)')
@click.option('--verbose/--quiet', default=True,
              help='Print progress messages')
@click.option('--progress/--no-progress', default=False,
              help='Show progress bar during MCMC (default: off)')
@click.option('--diagnostics/--no-diagnostics', default=False,
              help='Run convergence diagnostics after MCMC (default: off)')
def bayes(bfile, pheno, bayes, chain_length, burnin, thin, pi, hsq, mphen, out, verbose, progress, diagnostics):
    """
    Run Bayesian analysis on individual-level data.
    
    Requires PLINK format files (.bed/.bim/.fam) and phenotype file.
    
    Example:
        gctb bayes --bfile test/data/uk10k_chr1_1mb \\
                   --pheno test/data/test.phen \\
                   --bayes C --chain-length 1000 --out results
    """
    try:
        if verbose:
            click.echo("=" * 70)
            click.echo(f"GCTB Bayes{bayes.upper()} Analysis")
            click.echo("=" * 70)
        
        # Load data
        data = load_plink_data(bfile, pheno, mphen, verbose=verbose)
        
        # Build model
        if verbose:
            click.echo(f"\nBuilding Bayes{bayes.upper()} model...")
        
        model = gctb.build_model(data, bayes.upper(), heritability=hsq, pi=pi)
        
        if verbose:
            click.echo(f"  ✓ Model created with {model.num_snps} SNPs")
        
        # Run MCMC
        if verbose:
            click.echo(f"\nRunning MCMC...")
            click.echo(f"  Chain length: {chain_length}")
            click.echo(f"  Burn-in: {burnin}")
            click.echo(f"  Thinning: {thin}")
            click.echo("")
        
        # Progress bar (if enabled)
        if progress:
            try:
                from tqdm import tqdm
                import time
                import threading
                
                # Estimate time per iteration (very rough estimate)
                time_per_iter = 0.02  # 20ms per iteration (will vary by data size)
                total_time = chain_length * time_per_iter
                
                if verbose:
                    click.echo("  Note: Progress bar is an estimate (MCMC runs in C++)")
                
                # Create progress bar
                pbar = tqdm(total=chain_length, desc="  MCMC Progress", 
                           unit="iter", ncols=80)
                
                # Run MCMC in background, update progress bar
                results_container = []
                error_container = []
                
                def run_mcmc_thread():
                    try:
                        res = gctb.run_mcmc(
                            model=model,
                            chain_length=chain_length,
                            burnin=burnin,
                            thin=thin,
                            output_freq=max(chain_length // 10, 100),
                            title=out
                        )
                        results_container.append(res)
                    except Exception as e:
                        error_container.append(e)
                
                thread = threading.Thread(target=run_mcmc_thread)
                thread.start()
                
                # Update progress bar (estimated)
                start_time = time.time()
                while thread.is_alive():
                    elapsed = time.time() - start_time
                    estimated_iters = int(elapsed / time_per_iter)
                    current_pos = min(estimated_iters, chain_length)
                    pbar.n = current_pos
                    pbar.refresh()
                    time.sleep(0.5)
                
                # Ensure completion
                pbar.n = chain_length
                pbar.refresh()
                pbar.close()
                
                thread.join()
                
                if error_container:
                    raise error_container[0]
                
                results = results_container[0]
                
            except ImportError:
                if verbose:
                    click.echo("  Warning: tqdm not installed, progress bar disabled")
                results = gctb.run_mcmc(
                    model=model,
                    chain_length=chain_length,
                    burnin=burnin,
                    thin=thin,
                    output_freq=max(chain_length // 10, 100),
                    title=out
                )
        else:
            # No progress bar
            results = gctb.run_mcmc(
                model=model,
                chain_length=chain_length,
                burnin=burnin,
                thin=thin,
                output_freq=max(chain_length // 10, 100),
                title=out
            )
        
        if verbose:
            click.echo("\n  ✓ MCMC completed!")
        
        # Print summary results
        if verbose:
            click.echo("\n" + "=" * 70)
            click.echo("Results Summary:")
            click.echo("=" * 70)
            
            for res in results:
                if res.label == "hsq":
                    click.echo(f"  Heritability (h²):        {res.posterior_mean[0]:.4f}")
                elif res.label == "Pi":
                    click.echo(f"  Proportion non-zero (π):  {res.posterior_mean[0]:.4f}")
                elif res.label == "NnzSnp":
                    click.echo(f"  Non-zero SNPs:            {res.posterior_mean[0]:.1f}")
                elif res.label == "GenVar":
                    click.echo(f"  Genetic variance:         {res.posterior_mean[0]:.4f}")
                elif res.label == "ResVar":
                    click.echo(f"  Residual variance:        {res.posterior_mean[0]:.4f}")
        
        # Convergence diagnostics (if enabled)
        if diagnostics:
            try:
                from .diagnostics import check_convergence, recommend_chain_length
                diag_results = check_convergence(results, chain_length, burnin, verbose=verbose)
                
                # Recommendations
                needs_longer, message, _ = recommend_chain_length(diag_results)
                if verbose and needs_longer:
                    click.echo(f"\n{message}\n")
            except ImportError:
                if verbose:
                    click.echo("\n  Warning: scipy not installed, diagnostics unavailable")
            except Exception as e:
                if verbose:
                    click.echo(f"\n  Warning: Diagnostics failed: {e}")
        
        # Save results
        if verbose:
            click.echo(f"\nSaving results to {out}.*")
        
        # Save parameter results
        save_parameter_results(results, out)
        
        # Save SNP results
        save_snp_results(results, data, out)
        
        if verbose:
            click.echo(f"  ✓ Results saved to {out}.parRes and {out}.snpRes")
            click.echo("\n" + "=" * 70)
            click.echo("Analysis completed successfully!")
            click.echo("=" * 70)
        
    except Exception as e:
        click.echo(f"\nError: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command()
@click.option('--ldm', required=True, type=str,
              help='LD matrix directory or info file prefix')
@click.option('--gwas-summary', required=True, type=click.Path(exists=True),
              help='GWAS summary statistics file (.ma format)')
@click.option('--sbayes', type=click.Choice(['C', 'R', 'S'], case_sensitive=False),
              default='R', help='SBayes model type (default: R)')
@click.option('--chain-length', default=1000, type=int,
              help='MCMC chain length (default: 1000)')
@click.option('--burnin', default=100, type=int,
              help='Burn-in iterations (default: 100)')
@click.option('--thin', default=10, type=int,
              help='Thinning interval (default: 10)')
@click.option('--hsq', default=0.5, type=float,
              help='Heritability (default: 0.5)')
@click.option('--pi', default=0.01, type=float,
              help='Prior probability (default: 0.01)')
@click.option('--out', default='sbayes', type=str,
              help='Output prefix (default: sbayes)')
@click.option('--verbose/--quiet', default=True)
@click.option('--progress/--no-progress', default=False,
              help='Show progress bar during MCMC (default: off)')
@click.option('--diagnostics/--no-diagnostics', default=False,
              help='Run convergence diagnostics after MCMC (default: off)')
def sbayes(ldm, gwas_summary, sbayes, chain_length, burnin, thin, hsq, pi, out, verbose, progress, diagnostics):
    """
    Run SBayes analysis on GWAS summary statistics.
    
    Requires:
      1. LD matrix (.info and .bin files)
      2. GWAS summary statistics (.ma format)
    
    Example:
        gctb sbayes --ldm ldm_files/ldm \\
                    --gwas-summary summary.ma \\
                    --sbayes R --out sbayes_results
    """
    try:
        if verbose:
            click.echo("=" * 70)
            click.echo(f"GCTB SBayes{sbayes.upper()} Analysis (Summary Statistics)")
            click.echo("=" * 70)
        
        # Load summary statistics data
        data = load_summary_data(ldm, gwas_summary, verbose)
        
        # Build ApproxBayes model
        if verbose:
            click.echo(f"\nBuilding ApproxBayes{sbayes.upper()} model...")
        
        model = gctb.build_model_summary(data, sbayes.upper(), heritability=hsq, pi=pi)
        
        if verbose:
            click.echo(f"  ✓ Model created with {model.num_snps} SNPs")
        
        # Run MCMC
        if verbose:
            click.echo(f"\nRunning MCMC...")
            click.echo(f"  Chain length: {chain_length}")
            click.echo(f"  Burn-in: {burnin}")
            click.echo("")
        
        # Progress bar (if enabled)
        if progress:
            try:
                from tqdm import tqdm
                import time
                import threading
                
                # Estimate time per iteration
                time_per_iter = 0.02  # 20ms per iteration (rough estimate)
                
                if verbose:
                    click.echo("  Note: Progress bar is an estimate (MCMC runs in C++)")
                
                # Create progress bar
                pbar = tqdm(total=chain_length, desc="  MCMC Progress", 
                           unit="iter", ncols=80)
                
                # Run MCMC in background
                results_container = []
                error_container = []
                
                def run_mcmc_thread():
                    try:
                        res = gctb.run_mcmc(
                            model=model,
                            chain_length=chain_length,
                            burnin=burnin,
                            thin=thin,
                            output_freq=max(chain_length // 10, 100),
                            title=out
                        )
                        results_container.append(res)
                    except Exception as e:
                        error_container.append(e)
                
                thread = threading.Thread(target=run_mcmc_thread)
                thread.start()
                
                # Update progress bar (estimated)
                start_time = time.time()
                while thread.is_alive():
                    elapsed = time.time() - start_time
                    estimated_iters = int(elapsed / time_per_iter)
                    current_pos = min(estimated_iters, chain_length)
                    pbar.n = current_pos
                    pbar.refresh()
                    time.sleep(0.5)
                
                # Ensure completion
                pbar.n = chain_length
                pbar.refresh()
                pbar.close()
                
                thread.join()
                
                if error_container:
                    raise error_container[0]
                
                results = results_container[0]
                
            except ImportError:
                if verbose:
                    click.echo("  Warning: tqdm not installed, progress bar disabled")
                results = gctb.run_mcmc(
                    model=model,
                    chain_length=chain_length,
                    burnin=burnin,
                    thin=thin,
                    output_freq=max(chain_length // 10, 100),
                    title=out
                )
        else:
            # No progress bar
            results = gctb.run_mcmc(
                model=model,
                chain_length=chain_length,
                burnin=burnin,
                thin=thin,
                output_freq=max(chain_length // 10, 100),
                title=out
            )
        
        if verbose:
            click.echo("\n  ✓ MCMC completed!")
        
        # Print summary
        if verbose:
            click.echo("\n" + "=" * 70)
            click.echo("Results Summary:")
            click.echo("=" * 70)
            
            for res in results:
                if res.label == "hsq":
                    click.echo(f"  Heritability (h²):        {res.posterior_mean[0]:.4f}")
                elif res.label == "Pi":
                    click.echo(f"  Proportion non-zero (π):  {res.posterior_mean[0]:.4f}")
                elif res.label == "NnzSnp":
                    click.echo(f"  Non-zero SNPs:            {res.posterior_mean[0]:.1f}")
        
        # Convergence diagnostics (if enabled)
        if diagnostics:
            try:
                from .diagnostics import check_convergence, recommend_chain_length
                diag_results = check_convergence(results, chain_length, burnin, verbose=verbose)
                
                # Recommendations
                needs_longer, message, _ = recommend_chain_length(diag_results)
                if verbose and needs_longer:
                    click.echo(f"\n{message}\n")
            except ImportError:
                if verbose:
                    click.echo("\n  Warning: scipy not installed, diagnostics unavailable")
            except Exception as e:
                if verbose:
                    click.echo(f"\n  Warning: Diagnostics failed: {e}")
        
        # Save results
        if verbose:
            click.echo(f"\nSaving results to {out}.*")
        
        save_parameter_results(results, out)
        save_snp_results(results, data, out)
        
        if verbose:
            click.echo(f"  ✓ Results saved")
            click.echo("\n" + "=" * 70)
            click.echo("SBayes analysis completed successfully!")
            click.echo("=" * 70)
        
    except Exception as e:
        click.echo(f"\nError: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def load_summary_data(ldm_prefix: str, gwas_file: str, verbose: bool = True):
    """
    Load GWAS summary statistics and LD matrix.
    
    Parameters
    ----------
    ldm_prefix : str
        LD matrix file prefix (expects .info and .bin files)
    gwas_file : str
        GWAS summary statistics file
    verbose : bool
        Print progress messages
    """
    data = gctb.Data()
    
    # Load LD matrix info (this also loads SNP information)
    if verbose:
        click.echo("Loading LD matrix...")
    
    try:
        data.read_ld_matrix_info_file(f"{ldm_prefix}.info")
    except:
        # Try without extension
        data.read_ld_matrix_info_file(ldm_prefix)
    
    if verbose:
        click.echo(f"  ✓ LD matrix info loaded: {data.num_snps} SNPs")
    
    # Load GWAS summary statistics
    if verbose:
        click.echo("Loading GWAS summary statistics...")
    
    data.read_gwas_summary_file(
        gwas_file=gwas_file,
        af_diff=0.2,         # Default: filter if AF differs by >0.2
        maf_min=0.01,        # Default: minimum MAF 0.01
        maf_max=0.0,         # Default: no max filter
        pvalue_threshold=1.0, # Default: no p-value filter
        impute_n=False,
        remove_outlier_n=False
    )
    
    # Build included SNP list (CRITICAL for summary stats!)
    data.include_matched_snp()
    
    if verbose:
        click.echo(f"  ✓ GWAS summary loaded: {data.num_incd_snps} matched SNPs")
    
    # Load LD matrix binary data
    if verbose:
        click.echo("Loading LD matrix binary data...")
    
    try:
        data.read_ld_matrix_bin_file(f"{ldm_prefix}.bin")
    except:
        data.read_ld_matrix_bin_file(ldm_prefix)
    
    if verbose:
        click.echo(f"  ✓ LD matrix loaded")
    
    # Build sparse mixed model equations
    if verbose:
        click.echo("Building sparse mixed model equations...")
    
    data.build_sparse_mme(sample_overlap=False, noscale=False)
    
    if verbose:
        click.echo(f"  ✓ Sparse MME built")
        click.echo(f"\nReady for analysis: {data.num_incd_snps} SNPs")
    
    return data


def save_parameter_results(results, output_prefix):
    """Save non-SNP parameters to file"""
    with open(f"{output_prefix}.parRes", 'w') as f:
        # Header
        param_names = []
        for res in results:
            if res.label != "SnpEffects":
                param_names.append(res.label)
        
        f.write('\t'.join(param_names) + '\n')
        
        # Values (posterior means)
        values = []
        for res in results:
            if res.label != "SnpEffects":
                values.append(f"{res.posterior_mean[0]:.6f}")
        
        f.write('\t'.join(values) + '\n')


def save_snp_results(results, data, output_prefix):
    """Save SNP results to file"""
    snps = data.get_incd_snp_info_vec()
    
    # Find SNP effects and PIP
    snp_effects = None
    snp_pip = None
    
    for res in results:
        if res.label == "SnpEffects":
            snp_effects = res.posterior_mean
            snp_pip = res.pip
            break
    
    if snp_effects is None:
        return
    
    with open(f"{output_prefix}.snpRes", 'w') as f:
        # Header
        f.write("Index\tName\tChrom\tPosition\tA1\tA2\tFreq\tBeta\tPIP\n")
        
        # Data
        for i, snp in enumerate(snps):
            f.write(f"{i+1}\t{snp.ID}\t{snp.chrom}\t{snp.physPos}\t"
                   f"{snp.a1}\t{snp.a2}\t{snp.af:.4f}\t"
                   f"{snp_effects[i]:.6f}\t{snp_pip[i]:.4f}\n")


if __name__ == '__main__':
    main()

