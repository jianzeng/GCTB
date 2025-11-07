"""
GCTB Convergence Diagnostics

Provides tools to assess MCMC convergence and effective sample size.
"""

import numpy as np
from typing import List, Dict, Tuple


def geweke_test(samples: np.ndarray, first: float = 0.1, last: float = 0.5) -> Tuple[float, float]:
    """
    Geweke convergence diagnostic.
    
    Compares means from early and late portions of the chain.
    Z-score should be < 2 for convergence.
    
    Parameters
    ----------
    samples : np.ndarray
        MCMC samples (after burnin)
    first : float
        Fraction of chain to use for first window (default: 0.1 = first 10%)
    last : float
        Fraction of chain to use for last window (default: 0.5 = last 50%)
    
    Returns
    -------
    z_score : float
        Z-statistic
    p_value : float
        Two-tailed p-value
    """
    n = len(samples)
    
    if n < 20:
        return np.nan, np.nan
    
    # First window
    n_first = int(n * first)
    first_window = samples[:n_first]
    mean_first = np.mean(first_window)
    var_first = np.var(first_window, ddof=1) / n_first
    
    # Last window
    n_last = int(n * last)
    last_window = samples[-n_last:]
    mean_last = np.mean(last_window)
    var_last = np.var(last_window, ddof=1) / n_last
    
    # Z-score
    se = np.sqrt(var_first + var_last)
    if se == 0:
        return 0.0, 1.0
    
    z_score = (mean_first - mean_last) / se
    
    # P-value (two-tailed)
    from scipy import stats
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    
    return z_score, p_value


def effective_sample_size(samples: np.ndarray, max_lag: int = None) -> float:
    """
    Calculate effective sample size (ESS).
    
    ESS measures independent information in correlated samples.
    Higher is better. ESS ≈ n means independent samples.
    
    Parameters
    ----------
    samples : np.ndarray
        MCMC samples
    max_lag : int
        Maximum lag for autocorrelation (default: min(n//2, 1000))
    
    Returns
    -------
    ess : float
        Effective sample size
    """
    n = len(samples)
    
    if n < 10:
        return n
    
    if max_lag is None:
        max_lag = min(n // 2, 1000)
    
    # Compute autocorrelation
    samples_centered = samples - np.mean(samples)
    var = np.var(samples, ddof=1)
    
    if var == 0:
        return n
    
    # Autocorrelation at each lag
    rho = np.zeros(max_lag)
    for lag in range(max_lag):
        if lag >= n - 1:
            break
        rho[lag] = np.mean(samples_centered[:n-lag] * samples_centered[lag:]) / var
    
    # Sum autocorrelations until they become negative or insignificant
    tau = 1.0
    for lag in range(1, max_lag):
        if rho[lag] < 0.05:  # Stop when autocorrelation becomes small
            break
        tau += 2 * rho[lag]
    
    ess = n / tau
    
    return max(1.0, ess)  # ESS should be at least 1


def check_convergence(results, chain_length: int, burnin: int, 
                      verbose: bool = True) -> Dict[str, Dict]:
    """
    Check convergence for all parameters.
    
    Parameters
    ----------
    results : list
        MCMC results from run_mcmc()
    chain_length : int
        Total chain length
    burnin : int
        Burn-in period
    verbose : bool
        Print detailed results
    
    Returns
    -------
    diagnostics : dict
        Dictionary of diagnostics for each parameter
    """
    try:
        from scipy import stats
        has_scipy = True
    except ImportError:
        has_scipy = False
        if verbose:
            print("  Warning: scipy not available, some diagnostics unavailable")
    
    diagnostics = {}
    
    # Effective samples after burnin
    effective_samples = chain_length - burnin
    
    for res in results:
        if res.label == "SnpEffects":
            continue  # Skip SNP effects (too many to diagnose individually)
        
        param_name = res.label
        
        # Get samples (these are already post-burnin)
        try:
            samples = np.array(res.samples) if hasattr(res, 'samples') else None
        except:
            samples = None
        
        diag = {
            'mean': float(res.posterior_mean[0]) if hasattr(res, 'posterior_mean') else np.nan,
            'convergence': 'Unknown',
            'ess': np.nan,
            'geweke_z': np.nan,
            'geweke_p': np.nan
        }
        
        # If we have individual samples, compute diagnostics
        if samples is not None and len(samples) > 10:
            # Effective sample size
            if has_scipy:
                ess = effective_sample_size(samples)
                diag['ess'] = ess
                
                # Geweke test
                z_score, p_value = geweke_test(samples)
                diag['geweke_z'] = z_score
                diag['geweke_p'] = p_value
                
                # Assess convergence
                if not np.isnan(z_score):
                    if abs(z_score) < 2.0 and ess > 100:
                        diag['convergence'] = 'Good'
                    elif abs(z_score) < 3.0 and ess > 50:
                        diag['convergence'] = 'Fair'
                    else:
                        diag['convergence'] = 'Poor'
                else:
                    diag['convergence'] = 'Unknown'
        
        diagnostics[param_name] = diag
    
    # Print summary if verbose
    if verbose:
        print("\n" + "=" * 70)
        print("Convergence Diagnostics:")
        print("=" * 70)
        print(f"Chain length: {chain_length}, Burn-in: {burnin}")
        print(f"Effective samples: {effective_samples}")
        print("")
        
        if has_scipy:
            print(f"{'Parameter':<15} {'Mean':>10} {'ESS':>8} {'Geweke Z':>10} {'Status':<15}")
            print("-" * 70)
            
            for param, diag in diagnostics.items():
                ess_str = f"{diag['ess']:.0f}" if not np.isnan(diag['ess']) else "N/A"
                z_str = f"{diag['geweke_z']:.2f}" if not np.isnan(diag['geweke_z']) else "N/A"
                
                # Status symbol
                if diag['convergence'] == 'Good':
                    status = "✓ Good"
                elif diag['convergence'] == 'Fair':
                    status = "○ Fair"
                elif diag['convergence'] == 'Poor':
                    status = "⚠ Poor"
                else:
                    status = "? Unknown"
                
                print(f"{param:<15} {diag['mean']:>10.4f} {ess_str:>8} {z_str:>10} {status:<15}")
        else:
            print(f"{'Parameter':<15} {'Mean':>10}")
            print("-" * 30)
            for param, diag in diagnostics.items():
                print(f"{param:<15} {diag['mean']:>10.4f}")
        
        print("")
        print("Interpretation:")
        print("  ESS (Effective Sample Size): Higher is better")
        print("    > 400: Excellent")
        print("    > 100: Good")
        print("    > 50:  Fair")
        print("    < 50:  Poor (consider longer chain)")
        print("")
        print("  Geweke Z-score: Should be |Z| < 2 for convergence")
        print("    |Z| < 2: Good (chain converged)")
        print("    |Z| < 3: Fair (possibly converged)")
        print("    |Z| > 3: Poor (chain not converged, run longer)")
        print("=" * 70)
    
    return diagnostics


def recommend_chain_length(diagnostics: Dict[str, Dict]) -> Tuple[bool, str, int]:
    """
    Recommend whether chain should be longer based on diagnostics.
    
    Parameters
    ----------
    diagnostics : dict
        Output from check_convergence()
    
    Returns
    -------
    needs_longer : bool
        Whether chain should be run longer
    message : str
        Recommendation message
    recommended_length : int
        Recommended chain length (0 if current is sufficient)
    """
    # Check if any parameters have poor convergence
    poor_params = []
    fair_params = []
    min_ess = float('inf')
    
    for param, diag in diagnostics.items():
        if diag['convergence'] == 'Poor':
            poor_params.append(param)
        elif diag['convergence'] == 'Fair':
            fair_params.append(param)
        
        if not np.isnan(diag['ess']):
            min_ess = min(min_ess, diag['ess'])
    
    if poor_params:
        message = f"⚠️  Poor convergence detected for: {', '.join(poor_params[:3])}\n"
        message += "    Recommendation: Run 2-5x longer chain"
        return True, message, 0  # Multiply current by 2-5
    
    if fair_params and min_ess < 100:
        message = f"○ Fair convergence for: {', '.join(fair_params[:3])}\n"
        message += "    Recommendation: Consider 2x longer chain for better precision"
        return True, message, 0
    
    if min_ess < 50:
        message = f"⚠️  Low effective sample size (ESS={min_ess:.0f})\n"
        message += "    Recommendation: Run longer chain (2-3x current length)"
        return True, message, 0
    
    message = "✓ Convergence looks good! Chain length is sufficient."
    return False, message, 0

