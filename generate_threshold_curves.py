import numpy
from scipy.optimize import curve_fit
import coherent_sum
import tools.CoRaLs_geometry as aso_geometry
import json
import myplot
import matplotlib.pyplot as plt
import sys, os, argparse

def generatePowerSums(noise, window=160, step=40, save=True):
    """
    Build power sums with the SAME normalization used in snr_scan:
      1) Coherently sum antenna voltages
      2) Divide the coherent sum by sqrt(N_ant)
      3) Power frame = (1/window) * sum(v^2) over the window
    Accepts:
      noise: 1D (already coherent & normalized) OR 2D (n_ant, n_samples)
    """
    arr = numpy.asarray(noise)
    if arr.ndim == 2:
        n_ant = arr.shape[0]
        coh = arr.sum(axis=0) / numpy.sqrt(n_ant)
    else:
        n_ant = 1
        coh = arr
    power, _ = coherent_sum.powerSum(coh, window, step)
    if save:
        outpath = f'/home/aknicholas/Python Projects/Simulation/coralstrigger/noise/power_{window}_{step}.npy'
        numpy.save(outpath, power[0] if isinstance(power, tuple) else power)
        print(f"[generatePowerSums] Saved {outpath}  n_ant={n_ant}  norm=sqrt(N_ant)")
    return power


if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Generate threshold curve from a power sums .npy file")
    default_power = os.path.join(os.path.dirname(__file__), 'noise', 'power_160_40_10s.npy')
    parser.add_argument('--power-file', nargs='?', const=None,
                        default=default_power,
                        help='Optional path to power_<window>_<step>.npy. If used without value, uses default.')
    parser.add_argument('--window', type=int, default=160, help='Window (samples) used to make that file (informational)')
    parser.add_argument('--step', type=int, default=40, help='Step (samples) used in power sums')
    parser.add_argument('--fs', type=float, default=4e9, help='Sample rate (Hz)')
    parser.add_argument('--thr-min', type=float, default=0.2)
    parser.add_argument('--thr-max', type=float, default=8.0)
    parser.add_argument('--thr-step', type=float, default=0.05)

    # Original count-based fit region controls (kept)
    parser.add_argument('--fit-min-count', type=float, default=500, help='Min exceedance count to include in auto fit')
    parser.add_argument('--fit-max-count', type=float, default=6e4, help='Max exceedance count to include in auto fit')

    # New: manual / interactive fit controls
    parser.add_argument('--fit-thr-min', type=float, default=None, help='Override: lower threshold for fit region')
    parser.add_argument('--fit-thr-max', type=float, default=None, help='Override: upper threshold for fit region')
    parser.add_argument('--interactive-fit', action='store_true',
                        help='Show preliminary plot with errors then prompt for fit range')

    parser.add_argument('--n-beams', type=int, default=1, help='Number of active beams (for global false rate)')
    parser.add_argument('--target-global-rate', type=float, nargs='*',
                        default=[0.1],
                        help='List of target global false rates (Hz) to solve thresholds for')
    args = parser.parse_args()

    if args.power_file is None:
        args.power_file = default_power
        print(f"No path supplied after --power-file; using default: {args.power_file}")

    if not os.path.exists(args.power_file):
        print(f"Power file not found: {args.power_file}")
        sys.exit(1)

    powersums = numpy.load(args.power_file)
    print(f"Loaded power sums: {powersums.shape[0]} windows from {args.power_file}")
    print(f"Window={args.window} samples  Step={args.step} samples  fs={args.fs/1e9:.3f} GHz")

    # Threshold scan
    thresh_array = numpy.arange(args.thr_min, args.thr_max + 0.5*args.thr_step, args.thr_step)
    hits = numpy.array([(powersums >= thr).sum() for thr in thresh_array])

    total_windows = len(powersums)
    dt = args.step / args.fs
    obs_time = total_windows * dt
    per_beam_rate = hits / obs_time

    # Poisson errors: sigma_count = sqrt(N); sigma_rate = sqrt(N)/obs_time
    sigma_rate = numpy.sqrt(hits) / obs_time
    # For zero counts, assign an upper-limit style error (approx) to avoid zero weight
    zero = hits == 0
    if numpy.any(zero):
        sigma_rate[zero] = 1.0 / obs_time  # conservative

    # Optionally show preliminary plot for interactive range selection
    if args.interactive_fit:
        plt.figure(figsize=(7,5))
        plt.errorbar(thresh_array, per_beam_rate, yerr=sigma_rate, fmt='o', ms=4, capsize=3, label='Per-beam rate')
        plt.yscale('log')
        plt.grid(True, which='both', alpha=0.3)
        plt.xlabel('Normalized Power Threshold')
        plt.ylabel('Per-beam Noise Rate [Hz]')
        plt.title('Select Fit Range (Close window to continue)')
        plt.tight_layout()
        plt.show()
        try:
            user_min = input(f"Enter fit lower threshold (blank to keep {args.fit_thr_min}): ").strip()
            if user_min:
                args.fit_thr_min = float(user_min)
            user_max = input(f"Enter fit upper threshold (blank to keep {args.fit_thr_max}): ").strip()
            if user_max:
                args.fit_thr_max = float(user_max)
        except KeyboardInterrupt:
            print("\nAborted.")
            sys.exit(0)

    # Build fit mask
    if args.fit_thr_min is not None or args.fit_thr_max is not None:
        lo = args.fit_thr_min if args.fit_thr_min is not None else thresh_array.min()
        hi = args.fit_thr_max if args.fit_thr_max is not None else thresh_array.max()
        fit_mask = (thresh_array >= lo) & (thresh_array <= hi) & (hits > 0)
        print(f"Manual fit range: {lo} <= T <= {hi}  (points: {fit_mask.sum()})")
    else:
        fit_mask = (hits >= args.fit_min_count) & (hits <= args.fit_max_count) & (hits > 0)
        print(f"Auto fit mask by counts: {fit_mask.sum()} points "
              f"(counts in [{args.fit_min_count},{args.fit_max_count}])")

    fit_T = thresh_array[fit_mask]
    fit_rates = per_beam_rate[fit_mask]
    fit_sigma = sigma_rate[fit_mask]

    if fit_T.size < 3:
        print("Not enough points in fit region; adjust thresholds or counts.")
        sys.exit(1)

    # Model R(T)=A*exp(bT)
    def rate_model(T, A, b):
        return A * numpy.exp(b * T)

    # Initial guesses via log-linear least squares
    logR = numpy.log(fit_rates)
    w = 1.0 / fit_sigma
    p_lin = numpy.polyfit(fit_T, logR, 1, w=w)
    b0 = p_lin[0]
    A0 = numpy.exp(p_lin[1])

    # curve_fit with weights (sigma=fit_sigma)
    try:
        popt, pcov = curve_fit(rate_model, fit_T, fit_rates,
                               p0=[A0, b0],
                               sigma=fit_sigma,
                               absolute_sigma=True,
                               maxfev=10000)
    except RuntimeError as e:
        print(f"Fit failed: {e}")
        sys.exit(1)

    A_fit, b_fit = popt
    A_err, b_err = numpy.sqrt(numpy.diag(pcov))
    cov_A_b = pcov[0,1]
    corr_A_b = cov_A_b / (A_err * b_err) if A_err > 0 and b_err > 0 else numpy.nan

    # Goodness-of-fit
    model_fit = rate_model(fit_T, *popt)
    chi2 = numpy.sum(((fit_rates - model_fit)/fit_sigma)**2)
    ndf = fit_T.size - 2
    red_chi2 = chi2 / ndf if ndf > 0 else numpy.nan

    print(f"Fit: R(T) = A exp(b T)")
    print(f"  A = {A_fit:.3e} ± {A_err:.3e}")
    print(f"  b = {b_fit:.3e} ± {b_err:.3e}")
    print(f"  Corr(A,b) = {corr_A_b:.3f}")
    print(f"  chi2/ndf = {chi2:.2f}/{ndf} = {red_chi2:.2f}")

    # Threshold solver + uncertainty propagation
    # global_rate = n_beams * A exp(bT)
    # T = (ln(Rg) - ln(A n_beams))/b
    def solve_threshold_and_err(Rg):
        T_val = (numpy.log(Rg) - numpy.log(A_fit * args.n_beams)) / b_fit
        # partial derivatives
        dT_dA = -1.0 / (A_fit * b_fit)
        dT_db = -(numpy.log(Rg) - numpy.log(A_fit * args.n_beams)) / (b_fit**2)
        var_T = (dT_dA**2) * (A_err**2) + (dT_db**2) * (b_err**2) + 2 * dT_dA * dT_db * cov_A_b
        T_err = numpy.sqrt(var_T) if var_T > 0 else 0.0
        return T_val, T_err

    # Full model across scanned thresholds
    model_rate_full = rate_model(thresh_array, A_fit, b_fit)

    os.makedirs('plots', exist_ok=True)

    # Determine dynamic x-limits
    nonzero_idx = numpy.where(hits > 0)[0]
    if nonzero_idx.size:
        x_left = thresh_array[nonzero_idx[0]]
        x_right = thresh_array[nonzero_idx[-1]]
        if x_left == x_right:
            x_left -= args.thr_step
            x_right += args.thr_step
    else:
        x_left = thresh_array[0]; x_right = thresh_array[-1]
    pad = 0.05 * (x_right - x_left) if (x_right - x_left) > 0 else args.thr_step
    x_min_plot = max(thresh_array[0], x_left - pad)
    x_max_plot = min(thresh_array[-1], x_right + pad)
    # Remove entries where per_beam_rate is zero
    nonzero_mask = per_beam_rate > 1
    trimmed_thresh_array = thresh_array[nonzero_mask]
    trimmed_per_beam_rate = per_beam_rate[nonzero_mask]
    print(trimmed_thresh_array)
    print(trimmed_per_beam_rate)
    trimmed_sigma_rate = sigma_rate[nonzero_mask]
    print("Target global rates (Hz):", args.target_global_rate)
    # Plot with error bars + fit
    plt.figure(figsize=(7,5))
    plt.errorbar(trimmed_thresh_array, trimmed_per_beam_rate, yerr=trimmed_sigma_rate,
                 fmt='o', ms=4, capsize=3, label='Per-beam empirical')
    plt.semilogy(thresh_array, model_rate_full, '-', label='Fit')
    plt.semilogy(fit_T, fit_rates, 's', ms=6, label='Fit region')
    plt.xlabel('Normalized Power Threshold (P / (N σ²))')
    plt.ylabel('Per-beam Noise Rate [Hz]')
    plt.grid(True, which='both', alpha=0.3)
    # Plot target global rates as red horizontal lines and show intersections with fit line
    for gr in args.target_global_rate:
        T_needed, T_err = solve_threshold_and_err(gr)
        plt.axhline(gr / args.n_beams, color='red', linestyle='--', alpha=0.7, label=f'Target {gr:g} Hz' if gr == args.target_global_rate[0] else None)
        plt.plot([T_needed], [gr / args.n_beams], 'ro', markersize=7, label=f'Intersection {gr:g} Hz' if gr == args.target_global_rate[0] else None)
        plt.annotate(f'{T_needed:.2f}', xy=(T_needed, gr / args.n_beams), xytext=(5, 0), textcoords='offset points',
                     color='red', fontsize=8, va='center', ha='left', bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))
    plt.legend()
    plt.xlim(x_min_plot, x_max_plot+0.9)
    plt.ylim(1e-4, 1e14)

    txt = (f"A = {A_fit:.3e} ± {A_err:.1e}\n"
           f"b = {b_fit:.3e} ± {b_err:.1e}\n"
           f"corr = {corr_A_b:.2f}\n"
           f"χ²/ndf = {red_chi2:.2f}")
    plt.annotate(txt, xy=(0.02,0.02), xycoords='axes fraction',
                 fontsize=9, ha='left', va='bottom',
                 bbox=dict(boxstyle='round', fc='white', alpha=0.8))
    if args.fit_thr_min is not None or args.fit_thr_max is not None:
        plt.axvspan(
            args.fit_thr_min if args.fit_thr_min is not None else x_min_plot,
            args.fit_thr_max if args.fit_thr_max is not None else x_max_plot,
            color='orange', alpha=0.15, label='Manual fit span'
        )
    plt.tight_layout()
    fit_plot = f'plots/threshold_fit_{args.n_beams}_beams_with_errors.jpg'
    plt.savefig(fit_plot, dpi=150)
    print(f"Saved plot: {fit_plot}")

    # Histogram (log y)
    plt.figure(figsize=(7,5))
    plt.hist(powersums, bins=thresh_array, alpha=0.6)
    plt.yscale('log')
    plt.xlabel('Power')
    plt.ylabel('Counts')
    plt.grid(True, which='both', alpha=0.3)
    plt.xlim(x_min_plot, x_max_plot)
    plt.tight_layout()
    hist_plot = f'plots/threshold_{args.n_beams}_beams.jpg'
    plt.savefig(hist_plot, dpi = 150)
    print(f"Saved histogram: {hist_plot}")
    plt.show()

    # Prepare JSON output
    targets_out = []
    for gr in args.target_global_rate:
        T_needed, T_err = solve_threshold_and_err(gr)
        print(f"Target {gr:g} Hz: threshold {T_needed:.3f} ± {T_err:.3f}")
        targets_out.append({
            "rate": float(gr),
            "threshold": float(T_needed),
            "threshold_unc": float(T_err)
        })

    model_output = {
        "version": 3,
        "window": args.window,
        "step": args.step,
        "fs": args.fs,
        "n_beams": args.n_beams,
        "power_file": args.power_file,
        "fit_method": "weighted_curve_fit_exp",
        "fit_region": {
            "manual_thr_min": args.fit_thr_min,
            "manual_thr_max": args.fit_thr_max,
            "fit_min_count": args.fit_min_count,
            "fit_max_count": args.fit_max_count,
            "n_points": int(fit_T.size)
        },
        "fit_coefficients": {
            "A": float(A_fit),
            "A_unc": float(A_err),
            "b": float(b_fit),
            "b_unc": float(b_err),
            "cov_A_b": float(cov_A_b),
            "corr_A_b": float(corr_A_b),
            "chi2": float(chi2),
            "ndf": int(ndf),
            "chi2_ndf": float(red_chi2)
        },
        "target_global_rates": targets_out,
        "thresh_scan": {
            "thresholds": thresh_array.tolist(),
            "counts": hits.tolist(),
            "rates": per_beam_rate.tolist(),
            "rate_unc": sigma_rate.tolist()
        }
    }

    os.makedirs('plots', exist_ok=True)
    out_json = f'plots/threshold_model_{args.n_beams}_beams.json'
    with open(out_json, 'w') as f:
        json.dump(model_output, f, indent=2)
    print(f"Saved threshold fit JSON with uncertainties: {out_json}")
