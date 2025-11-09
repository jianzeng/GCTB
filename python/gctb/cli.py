#!/usr/bin/env python3
"""
GCTB Command Line Interface

Provides user-friendly commands for Bayesian genomic analysis.
"""

import click
import sys
from pathlib import Path
from typing import Optional, Sequence
from . import _core as gctb
from . import workflows

__version__ = "1.2.0"


def load_plink_data(
    bfile: str,
    pheno: Optional[str] = None,
    mphen: int = 1,
    keep_ind_file: Optional[str] = None,
    keep_ind_max: int = 999999,
    covariate_file: Optional[str] = None,
    random_covariate_file: Optional[str] = None,
    residual_diag_file: Optional[str] = None,
    include_snp_file: Optional[str] = None,
    exclude_snp_file: Optional[str] = None,
    include_chr: Optional[int] = None,
    include_block: Optional[int] = None,
    exclude_region_file: Optional[str] = None,
    skeleton_snp_file: Optional[str] = None,
    annotation_file: Optional[str] = None,
    annotation_transpose: bool = False,
    allow_multi_annotation: bool = True,
    continuous_annotation_file: Optional[str] = None,
    flank: int = 0,
    eqtl_file: Optional[str] = None,
    genetic_map_file: Optional[str] = None,
    ld_block_info_file: Optional[str] = None,
    exclude_ambiguous_snp: bool = False,
    exclude_mhc: bool = False,
    noscale: bool = False,
    read_genotypes: bool = True,
    verbose: bool = True,
):
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
    
    fam_path = f"{bfile}.fam"
    data.read_fam_file(fam_path)

    if covariate_file:
        data.read_covariate_file(covariate_file)
    if random_covariate_file:
        data.read_random_covariate_file(random_covariate_file)
    if residual_diag_file:
        data.read_residual_diag_file(residual_diag_file)

    pheno_path: Optional[str] = pheno
    if pheno_path is None:
        fallback = Path(fam_path)
        if fallback.exists():
            pheno_path = str(fallback)

    if keep_ind_file and pheno_path is None:
        raise click.BadParameter("Providing --keep-ind requires a phenotype file.")

    if pheno_path:
        data.read_phenotype_file(pheno_path, mphen)
        data.keep_matched_ind(keep_ind_file or "", keep_ind_max)
    else:
        data.build_kept_individuals()

    data.read_bim_file(f"{bfile}.bim")

    if include_snp_file:
        data.include_snp(include_snp_file)
    if exclude_snp_file:
        data.exclude_snp(exclude_snp_file)
    if include_chr:
        data.include_chr(int(include_chr))
    if include_block:
        data.include_block(int(include_block))
    if exclude_region_file:
        data.exclude_region(exclude_region_file)
    if exclude_mhc:
        data.exclude_mhc()
    if exclude_ambiguous_snp:
        data.exclude_ambiguous_snp()
    if skeleton_snp_file:
        data.include_skeleton_snp(skeleton_snp_file)
    if genetic_map_file:
        data.read_genetic_map_file(genetic_map_file)
    if ld_block_info_file:
        data.read_ld_block_info_file(ld_block_info_file)

    annotation_loaded = False
    if annotation_file:
        data.read_annotation_file(annotation_file, annotation_transpose, allow_multi_annotation)
        annotation_loaded = True
    elif continuous_annotation_file:
        data.read_annotation_file_format2(continuous_annotation_file, flank * 1000, eqtl_file or "")
        annotation_loaded = True

    data.include_matched_snp()

    if annotation_loaded and data.num_annos:
        data.set_annotation_info()

    if read_genotypes:
        data.read_bed_file(noscale, f"{bfile}.bed")

    if verbose:
        click.echo(f"  ✓ Loaded: {data.num_incd_snps} SNPs × {data.num_kept_inds} individuals")
        if data.num_kept_inds:
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
@click.option('--pheno', required=True, type=click.Path(exists=True, dir_okay=False),
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
@click.option('--keep-ind', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File listing individuals to keep (FID IID)')
@click.option('--keep-ind-max', default=999999, type=int,
              help='Maximum individuals to keep when filtering')
@click.option('--covar', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Covariate file')
@click.option('--random-covar', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Random covariate file')
@click.option('--residual', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Residual diagonal file')
@click.option('--include-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of SNPs to include')
@click.option('--exclude-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of SNPs to exclude')
@click.option('--include-chr', default=None, type=int,
              help='Restrict analysis to a chromosome')
@click.option('--include-block', default=None, type=int,
              help='Restrict analysis to an LD block (requires --ld-block-info)')
@click.option('--exclude-region', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of genomic regions to exclude')
@click.option('--skeleton-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Skeleton SNP list for LD operations')
@click.option('--annotation', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Categorical annotation file')
@click.option('--annotation-transpose', is_flag=True, default=False,
              help='Annotation file is transposed (rows=annotations)')
@click.option('--continuous-annotation', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Continuous annotation file (format 2)')
@click.option('--flank', default=0, type=int,
              help='Flanking distance (kb) when reading continuous annotations')
@click.option('--eqtl', default=None, type=click.Path(exists=True, dir_okay=False),
              help='eQTL file for continuous annotations')
@click.option('--genetic-map', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Genetic map file')
@click.option('--ld-block-info', default=None, type=click.Path(exists=True, dir_okay=False),
              help='LD block information file')
@click.option('--exclude-ambiguous', is_flag=True, default=False,
              help='Exclude ambiguous SNPs (A/T, C/G)')
@click.option('--exclude-mhc', is_flag=True, default=False,
              help='Exclude the MHC region')
@click.option('--noscale/--scale', default=False,
              help='Disable genotype scaling when reading BED (default: scale)')
@click.option('--out', default='gctb', type=str,
              help='Output prefix (default: gctb)')
@click.option('--verbose/--quiet', default=True,
              help='Print progress messages')
@click.option('--progress/--no-progress', default=False,
              help='Show progress bar during MCMC (default: off)')
@click.option('--diagnostics/--no-diagnostics', default=False,
              help='Run convergence diagnostics after MCMC (default: off)')
def bayes(
    bfile,
    pheno,
    bayes,
    chain_length,
    burnin,
    thin,
    pi,
    hsq,
    mphen,
    keep_ind,
    keep_ind_max,
    covar,
    random_covar,
    residual,
    include_snp,
    exclude_snp,
    include_chr,
    include_block,
    exclude_region,
    skeleton_snp,
    annotation,
    annotation_transpose,
    continuous_annotation,
    flank,
    eqtl,
    genetic_map,
    ld_block_info,
    exclude_ambiguous,
    exclude_mhc,
    noscale,
    out,
    verbose,
    progress,
    diagnostics,
):
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
        data = load_plink_data(
            bfile=bfile,
            pheno=pheno,
            mphen=mphen,
            keep_ind_file=keep_ind,
            keep_ind_max=keep_ind_max,
            covariate_file=covar,
            random_covariate_file=random_covar,
            residual_diag_file=residual,
            include_snp_file=include_snp,
            exclude_snp_file=exclude_snp,
            include_chr=include_chr,
            include_block=include_block,
            exclude_region_file=exclude_region,
            skeleton_snp_file=skeleton_snp,
            annotation_file=annotation,
            annotation_transpose=annotation_transpose,
            continuous_annotation_file=continuous_annotation,
            flank=flank,
            eqtl_file=eqtl,
            genetic_map_file=genetic_map,
            ld_block_info_file=ld_block_info,
            exclude_ambiguous_snp=exclude_ambiguous,
            exclude_mhc=exclude_mhc,
            noscale=noscale,
            read_genotypes=True,
            verbose=verbose,
        )
        
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
                
                thread = threading.Thread(target=run_mcmc_thread, daemon=True)
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
                
                # Wait for thread to complete
                thread.join(timeout=3600)  # 1 hour timeout
                
                # Ensure completion and cleanup
                pbar.n = chain_length
                pbar.refresh()
                pbar.close()
                
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


@main.group()
def ldmatrix():
    """LD matrix utilities built on the C++ core."""
    pass


@ldmatrix.command("make")
@click.option('--bfile', required=True, type=str,
              help='PLINK binary file prefix (.bed/.bim/.fam)')
@click.option('--pheno', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Phenotype file (defaults to .fam if omitted)')
@click.option('--mphen', default=1, type=int,
              help='Phenotype column to use (default: 1)')
@click.option('--out', required=True, type=str,
              help='Output prefix for LD matrix files')
@click.option('--type', 'ldm_type', default='sparse',
              type=click.Choice(['sparse', 'full', 'band', 'shrunk', 'block'], case_sensitive=False),
              help='LD matrix type to build (default: sparse)')
@click.option('--chisq-threshold', default=10.0, type=float,
              help='Chi-square threshold for pruning (default: 10)')
@click.option('--ld-threshold', default=0.0, type=float,
              help='LD r^2 threshold for pruning (default: 0.0)')
@click.option('--window-width', default=1000000, type=int,
              help='Window width in base pairs (default: 1,000,000)')
@click.option('--snp-range', default='', type=str,
              help='Restrict to SNP index range (e.g. "1-1000")')
@click.option('--write-txt/--no-write-txt', default=False,
              help='Also write LD matrix in text format')
@click.option('--keep-ind', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File listing individuals to keep')
@click.option('--keep-ind-max', default=999999, type=int,
              help='Maximum individuals to keep when filtering')
@click.option('--covar', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Covariate file')
@click.option('--random-covar', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Random covariate file')
@click.option('--residual', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Residual diagonal file')
@click.option('--include-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of SNPs to include')
@click.option('--exclude-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of SNPs to exclude')
@click.option('--include-chr', default=None, type=int,
              help='Restrict analysis to a chromosome')
@click.option('--include-block', default=None, type=int,
              help='Restrict analysis to an LD block (requires --ld-block-info)')
@click.option('--exclude-region', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of genomic regions to exclude')
@click.option('--skeleton-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Skeleton SNP list')
@click.option('--annotation', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Categorical annotation file')
@click.option('--annotation-transpose', is_flag=True, default=False,
              help='Annotation file is transposed (rows=annotations)')
@click.option('--continuous-annotation', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Continuous annotation file (format 2)')
@click.option('--flank', default=0, type=int,
              help='Flanking distance (kb) for continuous annotations')
@click.option('--eqtl', default=None, type=click.Path(exists=True, dir_okay=False),
              help='eQTL file for continuous annotations')
@click.option('--genetic-map', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Genetic map file')
@click.option('--ld-block-info', default=None, type=click.Path(exists=True, dir_okay=False),
              help='LD block information file')
@click.option('--exclude-ambiguous', is_flag=True, default=False,
              help='Exclude ambiguous SNPs (A/T, C/G)')
@click.option('--exclude-mhc', is_flag=True, default=False,
              help='Exclude the MHC region')
@click.option('--noscale/--scale', default=False,
              help='Disable genotype scaling when reading BED (default: scale)')
@click.option('--verbose/--quiet', default=True,
              help='Print progress messages')
def make_ldmatrix(
    bfile,
    pheno,
    mphen,
    out,
    ldm_type,
    chisq_threshold,
    ld_threshold,
    window_width,
    snp_range,
    write_txt,
    keep_ind,
    keep_ind_max,
    covar,
    random_covar,
    residual,
    include_snp,
    exclude_snp,
    include_chr,
    include_block,
    exclude_region,
    skeleton_snp,
    annotation,
    annotation_transpose,
    continuous_annotation,
    flank,
    eqtl,
    genetic_map,
    ld_block_info,
    exclude_ambiguous,
    exclude_mhc,
    noscale,
    verbose,
):
    """Build an LD matrix directly from PLINK genotypes."""
    try:
        if window_width < 0:
            raise click.BadParameter("window-width must be non-negative.")
        if not 0.0 <= ld_threshold <= 1.0:
            raise click.BadParameter("ld-threshold must be between 0 and 1.")
        ldm_type = ldm_type.lower()

        data = load_plink_data(
            bfile=bfile,
            pheno=pheno,
            mphen=mphen,
            keep_ind_file=keep_ind,
            keep_ind_max=keep_ind_max,
            covariate_file=covar,
            random_covariate_file=random_covar,
            residual_diag_file=residual,
            include_snp_file=include_snp,
            exclude_snp_file=exclude_snp,
            include_chr=include_chr,
            include_block=include_block,
            exclude_region_file=exclude_region,
            skeleton_snp_file=skeleton_snp,
            annotation_file=annotation,
            annotation_transpose=annotation_transpose,
            continuous_annotation_file=continuous_annotation,
            flank=flank,
            eqtl_file=eqtl,
            genetic_map_file=genetic_map,
            ld_block_info_file=ld_block_info,
            exclude_ambiguous_snp=exclude_ambiguous,
            exclude_mhc=exclude_mhc,
            noscale=noscale,
            read_genotypes=False,
            verbose=verbose,
        )

        if verbose:
            click.echo(f"\nBuilding {ldm_type.upper()} LD matrix -> {out}")

        bed_file = f"{bfile}.bed"
        data.make_ld_matrix(
            bed_file,
            ldm_type,
            chisq_threshold,
            ld_threshold,
            window_width,
            snp_range or "",
            out,
            write_txt,
        )

        if verbose:
            click.echo("  ✓ LD matrix construction complete")
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
@click.option('--chains', default=1, type=int,
              help='Number of MCMC chains (>=1). Default: 1')
@click.option('--include-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of SNPs to include')
@click.option('--exclude-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of SNPs to exclude')
@click.option('--include-chr', default=None, type=int,
              help='Restrict analysis to a chromosome')
@click.option('--include-block', default=None, type=int,
              help='Restrict analysis to an LD block')
@click.option('--exclude-region', default=None, type=click.Path(exists=True, dir_okay=False),
              help='File of genomic regions to exclude')
@click.option('--skeleton-snp', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Skeleton SNP list for LD operations')
@click.option('--annotation', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Categorical annotation file')
@click.option('--annotation-transpose', is_flag=True, default=False,
              help='Annotation file is transposed (rows=annotations)')
@click.option('--continuous-annotation', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Continuous annotation file (format 2)')
@click.option('--flank', default=0, type=int,
              help='Flanking distance (kb) for continuous annotations')
@click.option('--eqtl', default=None, type=click.Path(exists=True, dir_okay=False),
              help='eQTL file for continuous annotations')
@click.option('--ldscore', default=None, type=click.Path(exists=True, dir_okay=False),
              help='LD score file')
@click.option('--window-file', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Window definition file')
@click.option('--genetic-map', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Genetic map file (for LD post-processing)')
@click.option('--exclude-ambiguous', is_flag=True, default=False,
              help='Exclude ambiguous SNPs (A/T, C/G)')
@click.option('--exclude-mhc', is_flag=True, default=False,
              help='Exclude the MHC region')
@click.option('--multi-ldm/--single-ldm', default=False,
              help='Treat --ldm as multi-LD matrix directory (default: single)')
@click.option('--ldm-txt/--ldm-bin', default=False,
              help='Read LD matrix from .txt instead of .bin (default: bin)')
@click.option('--sample-overlap/--no-sample-overlap', default=False,
              help='Specify whether summary stats include sample overlap')
@click.option('--noscale/--scale', default=False,
              help='Disable scaling when building sparse MME (default: scale)')
@click.option('--rsq-threshold', default=1.0, type=float,
              help='LD r^2 threshold for pruning/binning')
@click.option('--bin-snp/--no-bin-snp', default=False,
              help='Bin SNPs by LD r^2 threshold instead of pruning')
@click.option('--out', default='sbayes', type=str,
              help='Output prefix (default: sbayes)')
@click.option('--verbose/--quiet', default=True)
@click.option('--progress/--no-progress', default=False,
              help='Show progress bar during MCMC (default: off)')
@click.option('--diagnostics/--no-diagnostics', default=False,
              help='Run convergence diagnostics after MCMC (default: off)')
def sbayes(
    ldm,
    gwas_summary,
    sbayes,
    chain_length,
    burnin,
    thin,
    hsq,
    pi,
    chains,
    include_snp,
    exclude_snp,
    include_chr,
    include_block,
    exclude_region,
    skeleton_snp,
    annotation,
    annotation_transpose,
    continuous_annotation,
    flank,
    eqtl,
    ldscore,
    window_file,
    genetic_map,
    exclude_ambiguous,
    exclude_mhc,
    multi_ldm,
    ldm_txt,
    sample_overlap,
    noscale,
    rsq_threshold,
    bin_snp,
    out,
    verbose,
    progress,
    diagnostics,
):
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
        data = load_summary_data(
            ldm_prefix=ldm,
            gwas_file=gwas_summary,
            verbose=verbose,
            include_snp_file=include_snp,
            exclude_snp_file=exclude_snp,
            include_chr=include_chr,
            include_block=include_block,
            exclude_region_file=exclude_region,
            skeleton_snp_file=skeleton_snp,
            annotation_file=annotation,
            annotation_transpose=annotation_transpose,
            continuous_annotation_file=continuous_annotation,
            flank=flank,
            eqtl_file=eqtl,
            ldscore_file=ldscore,
            window_file=window_file,
            genetic_map_file=genetic_map,
            exclude_ambiguous_snp=exclude_ambiguous,
            exclude_mhc=exclude_mhc,
            multi_ldm=multi_ldm,
            read_ldm_txt=ldm_txt,
            sample_overlap=sample_overlap,
            noscale=noscale,
            rsq_threshold=rsq_threshold,
            bin_snp=bin_snp,
            title=out,
        )
        
        if chains > 1:
            if verbose:
                click.echo(f"\nRunning {chains} independent chains...")
                if progress:
                    click.echo("  (Progress bar disabled for multi-chain mode)")
            results_by_chain = workflows.run_multi_chain_sbayes(
                data,
                sbayes_type=sbayes.upper(),
                heritability=hsq,
                pi=pi,
                num_chains=chains,
                chain_length=chain_length,
                burnin=burnin,
                thin=thin,
                output_prefix=out,
                random_start=True,
                verbose=verbose,
            )
            results = results_by_chain[0]
            if verbose:
                click.echo(f"  ✓ Completed {chains} chains (Gelman–Rubin written to {out}.gelman)")
        else:
            # Build ApproxBayes model
            if verbose:
                click.echo(f"\nBuilding ApproxBayes{sbayes.upper()} model...")
            
            model = gctb.build_model_summary(
                data,
                sbayes.upper(),
                heritability=hsq,
                pi=pi,
            )
            
            if verbose:
                click.echo(f"  ✓ Model created with {model.num_snps} SNPs")
                click.echo(f"\nRunning MCMC...")
                click.echo(f"  Chain length: {chain_length}")
                click.echo(f"  Burn-in: {burnin}")
                click.echo("")
            
            def _run_single_chain() -> Sequence[gctb.McmcSamples]:
                return gctb.run_mcmc(
                    model=model,
                    chain_length=chain_length,
                    burnin=burnin,
                    thin=thin,
                    output_freq=max(chain_length // 10, 100),
                    title=out,
                )
            
            if progress:
                try:
                    from tqdm import tqdm
                    import time
                    import threading

                    time_per_iter = 0.02  # rough estimate

                    if verbose:
                        click.echo("  Note: Progress bar is an estimate (MCMC runs in C++)")

                    pbar = tqdm(total=chain_length, desc="  MCMC Progress", unit="iter", ncols=80)
                    results_container: list[Sequence[gctb.McmcSamples]] = []
                    error_container: list[Exception] = []

                    def run_thread():
                        try:
                            results_container.append(_run_single_chain())
                        except Exception as exc:
                            error_container.append(exc)

                    thread = threading.Thread(target=run_thread, daemon=True)
                    thread.start()

                    start_time = time.time()
                    while thread.is_alive():
                        elapsed = time.time() - start_time
                        estimated_iters = int(elapsed / time_per_iter)
                        pbar.n = min(estimated_iters, chain_length)
                        pbar.refresh()
                        time.sleep(0.5)

                    thread.join(timeout=3600)
                    pbar.n = chain_length
                    pbar.refresh()
                    pbar.close()

                    if error_container:
                        raise error_container[0]

                    results = results_container[0]
                except ImportError:
                    if verbose:
                        click.echo("  Warning: tqdm not installed, progress bar disabled")
                    results = _run_single_chain()
            else:
                results = _run_single_chain()
        
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


def load_summary_data(
    ldm_prefix: str,
    gwas_file: str,
    verbose: bool = True,
    include_snp_file: Optional[str] = None,
    exclude_snp_file: Optional[str] = None,
    include_chr: Optional[int] = None,
    include_block: Optional[int] = None,
    exclude_region_file: Optional[str] = None,
    skeleton_snp_file: Optional[str] = None,
    annotation_file: Optional[str] = None,
    annotation_transpose: bool = False,
    allow_multi_annotation: bool = True,
    continuous_annotation_file: Optional[str] = None,
    flank: int = 0,
    eqtl_file: Optional[str] = None,
    ldscore_file: Optional[str] = None,
    window_file: Optional[str] = None,
    genetic_map_file: Optional[str] = None,
    exclude_ambiguous_snp: bool = False,
    exclude_mhc: bool = False,
    multi_ldm: bool = False,
    read_ldm_txt: bool = False,
    sample_overlap: bool = False,
    noscale: bool = False,
    rsq_threshold: float = 1.0,
    bin_snp: bool = False,
    title: Optional[str] = None,
):
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
    
    if verbose:
        click.echo("Loading LD matrix metadata...")
    
    if multi_ldm:
        data.read_multi_ld_matrix_info_file(ldm_prefix)
    else:
        info_path = f"{ldm_prefix}.info"
        try:
            data.read_ld_matrix_info_file(info_path)
        except Exception:
            data.read_ld_matrix_info_file(ldm_prefix)
    
    if verbose:
        click.echo(f"  ✓ LD matrix info loaded: {data.num_snps} SNPs")
    
    if include_snp_file:
        data.include_snp(include_snp_file)
    if exclude_snp_file:
        data.exclude_snp(exclude_snp_file)
    if include_chr:
        data.include_chr(int(include_chr))
    if include_block:
        data.include_block(int(include_block))
    if exclude_region_file:
        data.exclude_region(exclude_region_file)
    if exclude_mhc:
        data.exclude_mhc()
    if exclude_ambiguous_snp:
        data.exclude_ambiguous_snp()
    if skeleton_snp_file:
        data.include_skeleton_snp(skeleton_snp_file)
    if genetic_map_file:
        data.read_genetic_map_file(genetic_map_file)
    if ldscore_file:
        data.read_ldscore_file(ldscore_file)
    if window_file:
        data.read_window_file(window_file)

    annotation_loaded = False
    if annotation_file:
        data.read_annotation_file(annotation_file, annotation_transpose, allow_multi_annotation)
        annotation_loaded = True
    elif continuous_annotation_file:
        data.read_annotation_file_format2(continuous_annotation_file, flank * 1000, eqtl_file or "")
        annotation_loaded = True

    # Load GWAS summary statistics
    if verbose:
        click.echo("Loading GWAS summary statistics...")
    
    data.read_gwas_summary_file(
        gwas_file=gwas_file,
        af_diff=0.2,
        maf_min=0.01,
        maf_max=0.0,
        pvalue_threshold=1.0,
        impute_n=False,
        remove_outlier_n=False,
    )
    
    data.include_matched_snp()
    
    if annotation_loaded and data.num_annos:
        data.set_annotation_info()
    
    if verbose:
        click.echo(f"  ✓ GWAS summary loaded: {data.num_incd_snps} matched SNPs")
    
    if verbose:
        click.echo("Loading LD matrix data...")
    
    if read_ldm_txt:
        txt_path = f"{ldm_prefix}.txt"
        try:
            data.read_ld_matrix_txt_file(txt_path)
        except Exception:
            data.read_ld_matrix_txt_file(ldm_prefix)
    elif multi_ldm:
        data.read_multi_ld_matrix_bin_file(ldm_prefix)
    else:
        bin_path = f"{ldm_prefix}.bin"
        try:
            data.read_ld_matrix_bin_file(bin_path)
        except Exception:
            data.read_ld_matrix_bin_file(ldm_prefix)
    
    if rsq_threshold < 1.0 and not bin_snp:
        data.filter_snp_by_ld_rsq(rsq_threshold)
        data.include_matched_snp()

    if bin_snp:
        data.bin_snp_by_ld_rsq(rsq_threshold, title or Path(gwas_file).stem)

    if window_file:
        data.bin_snp_by_window_id()
    
    if verbose:
        click.echo("  ✓ LD matrix loaded")
        click.echo("Building sparse mixed model equations...")
    
    data.build_sparse_mme(sample_overlap=sample_overlap, noscale=noscale)
    
    if verbose:
        click.echo("  ✓ Sparse MME built")
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

