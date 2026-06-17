"""
Generate threshold curves for dual-polarization coincidence trigger.

Usage:
    python generate_threshold_curves_dualpol.py --n-beams 100 --target-rate 0.1
"""

import numpy as np
from scipy.optimize import curve_fit
import argparse
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import json

def exponential_rate(threshold, A, B):
    """
    Model: Rate(T) = A * exp(-B * T)
    
    Parameters
    ----------
    threshold : float or array
        Normalized power threshold
    A : float
        Pre-exponential factor (rate at T=0)
    B : float
        Decay constant (steepness)
    
    Returns
    -------
    rate : float or array
        Predicted rate in Hz
    """
    return A * np.exp(-B * threshold)


def fit_threshold_curve(thresholds, rates, sigma_rates=None, fit_range=None):
    """
    Fit exponential curve to threshold vs rate data.
    
    Parameters
    ----------
    thresholds : array
        Threshold values
    rates : array
        Measured rates (Hz)
    sigma_rates : array or None
        Uncertainties on rates (for weighted fit)
    fit_range : tuple (thr_min, thr_max) or None
        Threshold range to use for fitting
    
    Returns
    -------
    popt : array [A, B]
        Fitted parameters
    pcov : array (2, 2)
        Covariance matrix
    fit_mask : array (bool)
        Mask indicating which points were used in fit
    chisq_reduced : float
        Reduced chi-squared statistic
    """
    # Select fit region
    if fit_range is not None:
        thr_min, thr_max = fit_range
        fit_mask = (thresholds >= thr_min) & (thresholds <= thr_max) & (rates > 0)
    else:
        fit_mask = rates > 0
    
    if fit_mask.sum() < 3:
        raise ValueError(f"Not enough valid points for fitting: {fit_mask.sum()}")
    
    thr_fit = thresholds[fit_mask]
    rate_fit = rates[fit_mask]
    
    # Weighted fit if uncertainties provided
    if sigma_rates is not None:
        sigma_fit = sigma_rates[fit_mask]
        # Avoid division by zero
        sigma_fit[sigma_fit == 0] = np.min(sigma_fit[sigma_fit > 0])
        popt, pcov = curve_fit(exponential_rate, thr_fit, rate_fit, 
                              sigma=sigma_fit, absolute_sigma=True,
                              p0=[rate_fit[0], 1.0], maxfev=10000)
        
        # Calculate reduced chi-squared
        residuals = rate_fit - exponential_rate(thr_fit, *popt)
        chisq = np.sum((residuals / sigma_fit)**2)
        ndof = len(thr_fit) - 2  # 2 parameters
        chisq_reduced = chisq / ndof if ndof > 0 else np.nan
    else:
        popt, pcov = curve_fit(exponential_rate, thr_fit, rate_fit,
                              p0=[rate_fit[0], 1.0], maxfev=10000)
        chisq_reduced = np.nan
    
    return popt, pcov, fit_mask, chisq_reduced


def compute_coincidence_rate(rate_lhcp, rate_rhcp, window_time_ns=100):
    """
    Compute coincidence trigger rate for dual-pol system.
    
    For independent Poisson processes:
        Rate_coincidence ≈ 2 × Rate_LHCP × Rate_RHCP × τ_window
    
    Parameters
    ----------
    rate_lhcp : float or array
        LHCP trigger rate (Hz)
    rate_rhcp : float or array
        RHCP trigger rate (Hz)
    window_time_ns : float
        Coincidence window (ns)
    
    Returns
    -------
    rate_coincidence : float or array
        Coincidence rate (Hz)
    """
    window_time_s = window_time_ns * 1e-9
    return 2 * rate_lhcp * rate_rhcp * window_time_s


def solve_threshold_for_rate(A, B, target_rate):
    """
    Solve for threshold given target rate: T = -ln(Rate/A) / B
    
    Parameters
    ----------
    A, B : float
        Fitted exponential parameters
    target_rate : float
        Target rate (Hz)
    
    Returns
    -------
    threshold : float
        Required threshold value
    """
    if target_rate <= 0 or target_rate >= A:
        return np.nan
    return -np.log(target_rate / A) / B


def main():
    parser = argparse.ArgumentParser(description="Generate dual-pol threshold curves")
    parser.add_argument('--noise-dir', default='noise', help='Directory with power files')
    parser.add_argument('--window', type=int, default=160, help='Power window (samples)')
    parser.add_argument('--step', type=int, default=40, help='Power step (samples)')
    parser.add_argument('--fs', type=float, default=4e9, help='Sample rate (Hz)')
    parser.add_argument('--thr-min', type=float, default=0.5, help='Min threshold to scan')
    parser.add_argument('--thr-max', type=float, default=10.0, help='Max threshold to scan')
    parser.add_argument('--thr-step', type=float, default=0.05, help='Threshold step')
    parser.add_argument('--n-beams', type=int, default=100, help='Number of active beams (overridden by --beam-data-file)')
    parser.add_argument('--beam-data-file', type=str, default=None,
                       help='Path to beam optimization .npy file (e.g. plots/optimal_beams_sprint6.npy). '
                            'If provided, n_beams is set to the number of selected beams in the file.')
    parser.add_argument('--target-rate', type=float, nargs='+', default=[0.1], 
                       help='Target global trigger rates (Hz)')
    parser.add_argument('--coincidence-window', type=float, default=20.0,
                       help='Coincidence window (ns)')
    parser.add_argument('--fit-thr-min', type=float, default=None,
                       help='Min threshold for fit region (overrides auto-detection)')
    parser.add_argument('--fit-thr-max', type=float, default=None,
                       help='Max threshold for fit region (overrides auto-detection)')
    parser.add_argument('--fit-min-hits', type=float, default=300,
                       help='Min hits for auto fit region')
    parser.add_argument('--fit-max-hits', type=float, default=2e4,
                       help='Max hits for auto fit region')
    parser.add_argument('--output-dir', default='noise', help='Output directory')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    parser.add_argument('--plot-thr-max', type=float, default=5.3, 
                       help='Max threshold to show in plots (default 3.5)')
    parser.add_argument('--suffix', type=str, default='',
                       help='Suffix of power files to load, e.g. "_filt2_750MHz" '
                            '(default: no suffix = first-filter-only files)')
    parser.add_argument('--compare-suffix', type=str, default=None,
                       help='Load a second set of power files with this suffix for sanity-check overlay, '
                            'e.g. "_nofilter" to compare filtered vs unfiltered noise rates')
    
    args = parser.parse_args()
    
    # Override n_beams from beam optimization file if provided
    if args.beam_data_file is not None:
        beam_data = np.load(args.beam_data_file, allow_pickle=True).item()
        args.n_beams = len(beam_data['selected_beams'])
        print(f"Loaded beam count from {args.beam_data_file}: n_beams = {args.n_beams}")
    
    print("=" * 70)
    print("DUAL-POL THRESHOLD CURVE ANALYSIS")
    print("=" * 70)
    print(f"N_beams = {args.n_beams}")
    print(f"Coincidence window: {args.coincidence_window} ns")
    print("")
    
    # Load power files
    noise_path = Path(args.noise_dir)
    suffix = args.suffix
    lhcp_file = noise_path / f'power_lhcp_{args.window}_{args.step}{suffix}.npy'
    rhcp_file = noise_path / f'power_rhcp_{args.window}_{args.step}{suffix}.npy'
    
    if not lhcp_file.exists() or not rhcp_file.exists():
        print(f"Error: Power files not found in {args.noise_dir}/")
        print(f"  Missing: {lhcp_file} or {rhcp_file}")
        print(f"  Run generate_power_dualpol.py first")
        sys.exit(1)
    
    power_lhcp = np.load(lhcp_file)
    power_rhcp = np.load(rhcp_file)
    
    print(f"Loaded power sums:")
    print(f"  LHCP: {len(power_lhcp):,} frames")
    print(f"  RHCP: {len(power_rhcp):,} frames")
    print(f"  Window: {args.window} samples, Step: {args.step} samples")
    print("")
    
    # Calculate observation time
    dt = args.step / args.fs
    obs_time = len(power_lhcp) * dt
    print(f"Observation time: {obs_time:.3f} sec")
    print("")
    
    # Threshold scan
    print(f"Scanning thresholds: {args.thr_min} to {args.thr_max} step {args.thr_step}")
    thresholds = np.arange(args.thr_min, args.thr_max + 0.5*args.thr_step, args.thr_step)
    
    hits_lhcp = np.array([(power_lhcp >= thr).sum() for thr in thresholds])
    hits_rhcp = np.array([(power_rhcp >= thr).sum() for thr in thresholds])
    
    rate_lhcp = hits_lhcp / obs_time
    rate_rhcp = hits_rhcp / obs_time
    
    # Poisson uncertainties
    sigma_lhcp = np.sqrt(hits_lhcp) / obs_time
    sigma_rhcp = np.sqrt(hits_rhcp) / obs_time
    
    # Handle zero counts
    zero_lhcp = hits_lhcp == 0
    zero_rhcp = hits_rhcp == 0
    if np.any(zero_lhcp):
        sigma_lhcp[zero_lhcp] = 1.0 / obs_time
    if np.any(zero_rhcp):
        sigma_rhcp[zero_rhcp] = 1.0 / obs_time
    
    print(f"  LHCP: {hits_lhcp.sum():,} total hits")
    print(f"  RHCP: {hits_rhcp.sum():,} total hits")
    print("")
    
    # Fit exponential curves
    print("Fitting exponential curves: Rate(T) = A * exp(-B * T)")
    
    # Auto-select fit region based on counts (linear region in semilog)
    if args.fit_thr_min is None or args.fit_thr_max is None:
        # Find region with reasonable statistics: 100 < hits < 1e6
        min_hits = args.fit_min_hits
        max_hits = args.fit_max_hits
        good_region = (hits_lhcp > min_hits) & (hits_lhcp < max_hits)
        
        if good_region.sum() > 5:
            fit_thr_min_auto = thresholds[good_region][0]
            fit_thr_max_auto = thresholds[good_region][-1]
            print(f"  Auto-detected linear region: T ∈ [{fit_thr_min_auto:.2f}, {fit_thr_max_auto:.2f}]")
            print(f"    ({min_hits} < counts < {max_hits:.0e})")
            fit_range = (fit_thr_min_auto, fit_thr_max_auto)
        else:
            print(f"  Warning: Could not auto-detect fit region, using all data")
            fit_range = None
    else:
        fit_range = (args.fit_thr_min, args.fit_thr_max)
        print(f"  Manual fit range: T ∈ [{args.fit_thr_min}, {args.fit_thr_max}]")
    
    popt_lhcp, pcov_lhcp, mask_lhcp, chisq_lhcp = fit_threshold_curve(
        thresholds, rate_lhcp, sigma_lhcp, fit_range)
    popt_rhcp, pcov_rhcp, mask_rhcp, chisq_rhcp = fit_threshold_curve(
        thresholds, rate_rhcp, sigma_rhcp, fit_range)
    
    A_lhcp, B_lhcp = popt_lhcp
    A_rhcp, B_rhcp = popt_rhcp
    
    print(f"\n  LHCP: Rate(T) = {A_lhcp:.3e} * exp(-{B_lhcp:.4f} * T)")
    print(f"        Used {mask_lhcp.sum()}/{len(thresholds)} points")
    print(f"        χ²/dof = {chisq_lhcp:.3f}")
    print(f"  RHCP: Rate(T) = {A_rhcp:.3e} * exp(-{B_rhcp:.4f} * T)")
    print(f"        Used {mask_rhcp.sum()}/{len(thresholds)} points")
    print(f"        χ²/dof = {chisq_rhcp:.3f}")
    print("")
    
    # Compute per-beam coincidence rate for each threshold
    rate_coincidence_per_beam = compute_coincidence_rate(
        rate_lhcp, rate_rhcp, args.coincidence_window)
    
    # Global rate with N_beams
    rate_global = args.n_beams * rate_coincidence_per_beam
    
    # Solve for thresholds at target rates
    print(f"Target global rates (N_beams={args.n_beams}):")
    results = []
    
    for target in args.target_rate:
        # Target per-beam coincidence rate
        target_per_beam = target / args.n_beams
        
        # This is tricky: we need to solve for (T_L, T_R) such that:
        # 2 * Rate_L(T_L) * Rate_R(T_R) * τ = target_per_beam
        # For simplicity, assume symmetric: T_L = T_R = T
        # Then: 2 * [A_L * exp(-B_L*T)] * [A_R * exp(-B_R*T)] * τ = target_per_beam
        # 2 * A_L * A_R * τ * exp(-(B_L+B_R)*T) = target_per_beam
        # exp(-(B_L+B_R)*T) = target_per_beam / (2 * A_L * A_R * τ)
        # T = -ln(target_per_beam / (2*A_L*A_R*τ)) / (B_L+B_R)
        
        tau_s = args.coincidence_window * 1e-9
        arg = target_per_beam / (2 * A_lhcp * A_rhcp * tau_s)
        
        if arg > 0 and arg < 1:
            T_sym = -np.log(arg) / (B_lhcp + B_rhcp)
            rate_L_at_T = exponential_rate(T_sym, A_lhcp, B_lhcp)
            rate_R_at_T = exponential_rate(T_sym, A_rhcp, B_rhcp)
            rate_coinc_check = compute_coincidence_rate(rate_L_at_T, rate_R_at_T, 
                                                        args.coincidence_window)
            rate_global_check = args.n_beams * rate_coinc_check
            
            results.append({
                'target_global_rate_hz': target,
                'target_per_beam_rate_hz': target_per_beam,
                'threshold_lhcp': T_sym,
                'threshold_rhcp': T_sym,
                'rate_lhcp_at_threshold': rate_L_at_T,
                'rate_rhcp_at_threshold': rate_R_at_T,
                'coincidence_rate_per_beam': rate_coinc_check,
                'global_rate_check': rate_global_check
            })
            
            print(f"\n  Target: {target} Hz global")
            print(f"    Per-beam: {target_per_beam:.3e} Hz")
            print(f"    Threshold (symmetric): T_L = T_R = {T_sym:.3f}")
            print(f"    Rate_L(T): {rate_L_at_T:.3e} Hz")
            print(f"    Rate_R(T): {rate_R_at_T:.3e} Hz")
            print(f"    Coincidence per beam: {rate_coinc_check:.3e} Hz")
            print(f"    Global rate: {rate_global_check:.3e} Hz ✓")
        else:
            print(f"\n  Target: {target} Hz - WARNING: Cannot solve (out of range)")
            results.append({
                'target_global_rate_hz': target,
                'threshold_lhcp': np.nan,
                'threshold_rhcp': np.nan,
                'error': 'out_of_range'
            })
    
    # Save results
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)
    
    file_suffix = args.suffix  # e.g. '' or '_filt2_750MHz'
    results_file = output_path / f'threshold_analysis_dualpol{file_suffix}.json'
    with open(results_file, 'w') as f:
        json.dump({
            'n_beams': args.n_beams,
            'coincidence_window_ns': args.coincidence_window,
            'observation_time_sec': obs_time,
            'power_file_suffix': file_suffix,
            'fit_parameters': {
                'lhcp': {'A': float(A_lhcp), 'B': float(B_lhcp)},
                'rhcp': {'A': float(A_rhcp), 'B': float(B_rhcp)}
            },
            'thresholds': results
        }, f, indent=2)
    
    print(f"\n{'='*70}")
    print(f"Saved: {results_file}")
    
    # Save threshold scan data
    scan_file = output_path / f'threshold_scan_dualpol{file_suffix}.npz'
    np.savez(scan_file,
             thresholds=thresholds,
             rate_lhcp=rate_lhcp,
             rate_rhcp=rate_rhcp,
             rate_coincidence_per_beam=rate_coincidence_per_beam,
             rate_global=rate_global,
             sigma_lhcp=sigma_lhcp,
             sigma_rhcp=sigma_rhcp,
             fit_lhcp=exponential_rate(thresholds, A_lhcp, B_lhcp),
             fit_rhcp=exponential_rate(thresholds, A_rhcp, B_rhcp))
    print(f"Saved: {scan_file}")
    
    # Plot if requested
    if args.plot:
        plots_path = Path('plots')
        plots_path.mkdir(exist_ok=True)

        # Restrict plotting to region with good statistics
        plot_mask = thresholds <= args.plot_thr_max
        thr_plot = thresholds[plot_mask]
        rate_lhcp_plot = rate_lhcp[plot_mask]
        rate_rhcp_plot = rate_rhcp[plot_mask]
        sigma_lhcp_plot = sigma_lhcp[plot_mask]
        sigma_rhcp_plot = sigma_rhcp[plot_mask]
        rate_coinc_plot = rate_coincidence_per_beam[plot_mask]
        rate_global_plot = rate_global[plot_mask]
        
        # Also mask the fit masks
        mask_lhcp_plot = mask_lhcp[plot_mask]
        mask_rhcp_plot = mask_rhcp[plot_mask]
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Per-polarization rates
        ax = axes[0, 0]
        ax.errorbar(thr_plot, rate_lhcp_plot, yerr=sigma_lhcp_plot, fmt='o', ms=3, 
                   alpha=0.3, color='gray', label='LHCP data')
        ax.errorbar(thr_plot[mask_lhcp_plot], rate_lhcp_plot[mask_lhcp_plot], 
                   yerr=sigma_lhcp_plot[mask_lhcp_plot], fmt='o', ms=4, 
                   alpha=0.8, color='blue', label='Fit region')
        ax.plot(thr_plot, exponential_rate(thr_plot, A_lhcp, B_lhcp),
               'r-', lw=2, label=f'Fit: {A_lhcp:.2e}*exp(-{B_lhcp:.3f}*T)')
        ax.text(0.05, 0.05, f'χ²/dof = {chisq_lhcp:.3f}', 
               transform=ax.transAxes, fontsize=11, verticalalignment='bottom',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_yscale('log')
        ax.set_xlabel('Threshold (normalized power)')
        ax.set_ylabel('Per-beam Rate [Hz]')
        ax.set_title('LHCP Threshold Curve')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        ax = axes[0, 1]
        ax.errorbar(thr_plot, rate_rhcp_plot, yerr=sigma_rhcp_plot, fmt='o', ms=3,
                   alpha=0.3, color='gray', label='RHCP data')
        ax.errorbar(thr_plot[mask_rhcp_plot], rate_rhcp_plot[mask_rhcp_plot],
                   yerr=sigma_rhcp_plot[mask_rhcp_plot], fmt='o', ms=4,
                   alpha=0.8, color='blue', label='Fit region')
        ax.plot(thr_plot, exponential_rate(thr_plot, A_rhcp, B_rhcp),
               'r-', lw=2, label=f'Fit: {A_rhcp:.2e}*exp(-{B_rhcp:.3f}*T)')
        ax.text(0.05, 0.05, f'χ²/dof = {chisq_rhcp:.3f}', 
               transform=ax.transAxes, fontsize=11, verticalalignment='bottom',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.set_yscale('log')
        ax.set_xlabel('Threshold (normalized power)')
        ax.set_ylabel('Per-beam Rate [Hz]')
        ax.set_title('RHCP Threshold Curve')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Coincidence rates
        ax = axes[1, 0]
        ax.plot(thr_plot, rate_coinc_plot, 'b-', lw=2, label='Per-beam coincidence')
        ax.set_yscale('log')
        ax.set_xlabel('Threshold (normalized power)')
        ax.set_ylabel('Coincidence Rate [Hz]')
        ax.set_title(f'Coincidence Rate (τ={args.coincidence_window} ns)')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Global rates
        ax = axes[1, 1]
        ax.plot(thr_plot, rate_global_plot, 'g-', lw=2, label=f'Global (N={args.n_beams} beams)')
        for res in results:
            if 'global_rate_check' in res:
                if res['threshold_lhcp'] <= args.plot_thr_max:
                    ax.axhline(res['target_global_rate_hz'], ls='--', alpha=0.5,
                              label=f"Target: {res['target_global_rate_hz']} Hz")
                    ax.axvline(res['threshold_lhcp'], ls='--', alpha=0.5)
        ax.set_yscale('log')
        ax.set_xlabel('Threshold (normalized power)')
        ax.set_ylabel('Global Trigger Rate [Hz]')
        ax.set_title(f'Global Rate (N_beams={args.n_beams})')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        plot_file = plots_path / f'threshold_curves_dualpol{file_suffix}.png'
        plt.savefig(plot_file, dpi=150)
        print(f"Saved: {plot_file}")
        plt.show()

    # -------------------------------------------------------------------------
    # Optional comparison: overlay a second power file (e.g. unfiltered)
    # -------------------------------------------------------------------------
    if args.compare_suffix is not None:
        sfx = args.compare_suffix
        cmp_lhcp_file = noise_path / f'power_lhcp_{args.window}_{args.step}{sfx}.npy'
        cmp_rhcp_file = noise_path / f'power_rhcp_{args.window}_{args.step}{sfx}.npy'

        if not cmp_lhcp_file.exists() or not cmp_rhcp_file.exists():
            print(f"\nWarning: comparison files not found ({cmp_lhcp_file}), skipping compare plot.")
        else:
            cmp_lhcp = np.load(cmp_lhcp_file)
            cmp_rhcp = np.load(cmp_rhcp_file)
            cmp_obs = len(cmp_lhcp) * dt

            cmp_hits_lhcp = np.array([(cmp_lhcp >= thr).sum() for thr in thresholds])
            cmp_hits_rhcp = np.array([(cmp_rhcp >= thr).sum() for thr in thresholds])
            cmp_rate_lhcp = cmp_hits_lhcp / cmp_obs
            cmp_rate_rhcp = cmp_hits_rhcp / cmp_obs

            print(f"\n--- Comparison ({sfx}) ---")
            print(f"  LHCP: {cmp_hits_lhcp.sum():,} total hits over {cmp_obs:.3f} sec")
            print(f"  RHCP: {cmp_hits_rhcp.sum():,} total hits over {cmp_obs:.3f} sec")

            # Overlay plot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            label_base = 'filtered (default)'
            label_cmp = sfx.strip('_')

            for ax, r_base, r_cmp, pol in [
                (ax1, rate_lhcp, cmp_rate_lhcp, 'LHCP'),
                (ax2, rate_rhcp, cmp_rate_rhcp, 'RHCP')
            ]:
                pm = thresholds <= args.plot_thr_max
                ax.semilogy(thresholds[pm], r_base[pm], 'b-o', ms=3, lw=1.5, label=label_base)
                ax.semilogy(thresholds[pm], r_cmp[pm], 'r--s', ms=3, lw=1.5, label=label_cmp)
                ax.set_xlabel('Threshold (normalized power)')
                ax.set_ylabel('Per-beam Rate [Hz]')
                ax.set_title(f'{pol} noise rate: filter comparison')
                ax.legend()
                ax.grid(True, alpha=0.3)

            plt.tight_layout()
            cmp_plot_file = plots_path / f'threshold_compare{sfx}.png'
            plt.savefig(cmp_plot_file, dpi=150)
            print(f"Saved: {cmp_plot_file}")
            plt.show()

    print("=" * 70)


if __name__ == '__main__':
    main()
