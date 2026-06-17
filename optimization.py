"""
Optimize SNR 50% efficiency for various parameters

Current working models require:
Power Windowing
Filtering methods
Wave parameters (phi,psi,theta)
Coincidence window

I'm thinking I'll start with single variable plots and then if I want to get fancy I'll do multivar heatmaps.
I need to find a way to change our filter, I don't fully understand FIR enough to build my own from scratch
"""
import numpy as np
import matplotlib.pyplot as plt
import tools.CoRaLs_geometry as geometry
import tools.constants as constants
import payload_signal as payload
import coherent_sum as trigger
import noise
import tools.filters as filters
import json
import argparse
import sys
import os
import snr_scan_dualpol as snr_scan


def _interpolate_snr50(snr, eff):
    """Return SNR at 50% coincidence efficiency by linear interpolation, or NaN."""
    idx = np.where(eff >= 0.5)[0]
    if len(idx) > 0:
        k = idx[0]
        if k > 0:
            x0, x1 = snr[k - 1], snr[k]
            y0, y1 = eff[k - 1], eff[k]
            return x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else x1
        return snr[0]
    return np.nan


def optimize_windows(phi, theta, psi, thresholds_per_window, coincidence_window,
                     snr_grid, n_trials, window_grid, step, save_filename=None,
                     antennas=None):
    """
    Scan over power-sum window sizes and find the SNR at 50% coincidence
    efficiency for each window.

    The x-axis is window size in ns; the y-axis is the SNR required to reach
    50% efficiency (lower = better).

    Parameters:
    -----------
    phi, theta, psi : float
        Signal arrival angles (deg)
    thresholds_per_window : dict or float
        Per-window thresholds: dict mapping window (samples, int) → threshold
        float, or a single float used for all windows. Using a per-window dict
        (one threshold JSON per window size) is strongly preferred — a fixed
        threshold makes smaller windows look artificially better because the
        noise tail statistics change with window size.
    coincidence_window : float
        Coincidence time window (ns)
    snr_grid : array
        Per-antenna SNR values to test at each window size
    n_trials : int
        Number of trials per SNR point
    window_grid : tuple
        (min_samples, max_samples, samples_per_step) — range of power-sum
        window sizes to scan
    step : int
        Power-sum step size (samples); held fixed across the scan
    save_filename : str, optional
        Base filename for saving results (.npy + .txt)
    antennas : list of int or None
        Physical antenna indices to use (default: all four)

    Returns:
    --------
    dict with keys: windows, windows_ns, snr50, phi, theta, psi,
                    thresholds_per_window, step, snr_grid, n_trials,
                    coincidence_window
    """

    win_min, win_max, win_step = window_grid
    windows = np.arange(win_min, win_max + 1, win_step, dtype=int)

    # Drop windows smaller than the stride — those leave gaps between frames.
    invalid = windows < step
    if invalid.any():
        print(f"WARNING: skipping {invalid.sum()} window value(s) smaller than "
              f"stride ({step} samples): {windows[invalid].tolist()}")
        windows = windows[~invalid]
    if len(windows) == 0:
        raise ValueError(
            f"No valid window sizes remain after removing values < step ({step}). "
            f"Increase --window-min or decrease --step.")

    # ns per ADC sample
    sample_ns = geometry.ritc_sample_step

    print(f"\n=== Window Optimization ===")
    print(f"Direction: phi={phi:.1f}°, theta={theta:.1f}°, psi={psi:.1f}°")
    if isinstance(thresholds_per_window, dict):
        print(f"Thresholds: per-window (calibrated)  |  Coincidence window: {coincidence_window:.1f} ns")
    else:
        print(f"Threshold: {float(thresholds_per_window):.4f}  |  Coincidence window: {coincidence_window:.1f} ns")
    print(f"Window range: {windows[0]}–{windows[-1]} samples "
          f"({windows[0]*sample_ns:.1f}–{windows[-1]*sample_ns:.1f} ns)  "
          f"[stride fixed at {step} samples = {step*sample_ns:.1f} ns]")
    print(f"SNR points: {len(snr_grid)}  |  Trials: {n_trials}")

    snr50_list = []

    for window in windows:
        window_ns = window * sample_ns
        if isinstance(thresholds_per_window, dict):
            thr = thresholds_per_window.get(
                int(window),
                next(iter(thresholds_per_window.values()))
            )
        else:
            thr = float(thresholds_per_window)
        print(f"\n--- Window = {window} samples ({window_ns:.1f} ns), threshold={thr:.4f} ---")

        res = snr_scan.run_snr_scan_dualpol(
            phi=phi,
            theta=theta,
            psi=psi,
            threshold_lhcp=thr,
            threshold_rhcp=thr,
            coincidence_window=coincidence_window,
            snr_grid=snr_grid,
            n_trials=n_trials,
            window=window,
            step=step,
            save_filename=None,
            antennas=antennas,
        )

        # Interpolate SNR at 50% coincidence efficiency
        snr_50 = _interpolate_snr50(res['snr'], res['eff_coinc'])

        snr50_list.append(snr_50)
        print(f"  → SNR50 = {snr_50:.3f}")

    windows_ns = windows * sample_ns

    opt_results = {
        'windows':              windows,
        'windows_ns':           windows_ns,
        'snr50':                np.array(snr50_list),
        'phi':                  phi,
        'theta':                theta,
        'psi':                  psi,
        'thresholds_per_window': thresholds_per_window,
        'step':                 step,
        'snr_grid':             snr_grid,
        'n_trials':             n_trials,
        'coincidence_window':   coincidence_window,
        'antennas':             antennas,
    }

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        np.save(f'plots/{save_filename}.npy', opt_results)

        thr_note = ('per-window' if isinstance(thresholds_per_window, dict)
                    else f'{float(thresholds_per_window):.3f}')
        data = np.column_stack([windows, windows_ns, snr50_list])
        header = (
            f"Window optimization at (phi={phi:.1f}, theta={theta:.1f}, psi={psi:.1f})\n"
            f"Threshold={thr_note}, Power-step={step}, "
            f"Coinc={coincidence_window:.1f} ns, Trials={n_trials}\n"
            "Window(samples)\tWindow(ns)\tSNR50"
        )
        np.savetxt(f'plots/{save_filename}.txt', data, fmt='%.6g', header=header)
        print(f"\nSaved: plots/{save_filename}.npy and .txt")

    return opt_results


def optimize_filter_freq(phi, theta, psi, fc_hz_list, thresholds_per_fc,
                         coincidence_window, snr_grid, n_trials,
                         window, step, save_filename=None, antennas=None):
    """
    Scan over filter cutoff frequencies and find SNR at 50% coincidence efficiency.

    The first entry in fc_hz_list is treated as the baseline (no second filter —
    only the always-on Shannon-Whitaker 1.5 GHz front-end filter is applied).
    Every subsequent entry applies a second lowpass FIR at that cutoff frequency
    after beamforming, lowering the noise floor and potentially improving sensitivity.

    Parameters
    ----------
    phi, theta, psi : float
        Signal arrival angles (deg).
    fc_hz_list : list of float
        Cutoff frequencies in Hz to scan, e.g. [1.5e9, 1.25e9, 1.0e9, 750e6].
        The first value is the no-second-filter baseline.
    thresholds_per_fc : dict or float
        Mapping {fc_hz: threshold} for LHCP and RHCP (assumed symmetric).
        Pass a single float to use the same threshold for all frequencies.
        Thresholds should be calibrated at each (fc, window, step) combination
        via generate_threshold_curves_dualpol.py.
    coincidence_window : float
        Coincidence time window (ns).
    snr_grid : array
        Per-antenna SNR values to test.
    n_trials : int
        Trials per SNR point.
    window, step : int
        Power-sum window and stride (samples).
    save_filename : str or None
        Base filename; saves plots/<save_filename>_filter_opt.npy/.txt.
    antennas : list of int or None
        Physical antenna indices (default: all four).

    Returns
    -------
    dict with keys:
        fc_hz, fc_mhz, snr50, noise_rms_f1, noise_rms_f2,
        phi, theta, psi, window, step, coincidence_window, n_trials
    """
    sample_ns = geometry.ritc_sample_step
    ref_fc = fc_hz_list[0]

    print(f"\n=== Filter Frequency Optimization ===")
    print(f"Direction: phi={phi:.1f}°, theta={theta:.1f}°, psi={psi:.1f}°")
    print(f"Window: {window} samples ({window*sample_ns:.1f} ns)  |  "
          f"Step: {step} samples ({step*sample_ns:.1f} ns)")
    print(f"FC values: {[f/1e6 for f in fc_hz_list]} MHz")

    snr50_list     = []
    noise_rms_f1_list = []
    noise_rms_f2_list = []

    for i, fc in enumerate(fc_hz_list):
        fc_mhz = fc / 1e6
        apply_second = (fc != ref_fc)
        fc_second_arg = fc if apply_second else None

        if isinstance(thresholds_per_fc, dict):
            thr = thresholds_per_fc.get(fc, thresholds_per_fc[ref_fc])
        else:
            thr = float(thresholds_per_fc)

        print(f"\n--- fc = {fc_mhz:.0f} MHz  (second filter: {'ON' if apply_second else 'OFF'})  "
              f"threshold = {thr:.4f} ---")

        res = snr_scan.run_snr_scan_dualpol(
            phi=phi,
            theta=theta,
            psi=psi,
            threshold_lhcp=thr,
            threshold_rhcp=thr,
            coincidence_window=coincidence_window,
            snr_grid=snr_grid,
            n_trials=n_trials,
            window=window,
            step=step,
            save_filename=None,
            antennas=antennas,
            apply_second_filter=apply_second,
            fc_second=fc_second_arg,
        )

        snr_50 = _interpolate_snr50(res['snr'], res['eff_coinc'])
        snr50_list.append(snr_50)
        noise_rms_f1_list.append(res.get('noise_rms_f1'))
        noise_rms_f2_list.append(res.get('noise_rms_f2'))
        print(f"  → SNR50 = {snr_50:.3f}")

    fc_hz_arr  = np.array(fc_hz_list)
    fc_mhz_arr = fc_hz_arr / 1e6

    opt_results = {
        'fc_hz':              fc_hz_arr,
        'fc_mhz':             fc_mhz_arr,
        'snr50':              np.array(snr50_list),
        'noise_rms_f1':       np.array(noise_rms_f1_list, dtype=object),
        'noise_rms_f2':       np.array(noise_rms_f2_list, dtype=object),
        'phi':                phi,
        'theta':              theta,
        'psi':                psi,
        'window':             window,
        'step':               step,
        'coincidence_window': coincidence_window,
        'n_trials':           n_trials,
        'antennas':           antennas,
    }

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        np.save(f'plots/{save_filename}_filter_opt.npy', opt_results)
        data = np.column_stack([fc_mhz_arr, snr50_list])
        header = (
            f"Filter frequency optimization at (phi={phi:.1f}, theta={theta:.1f}, psi={psi:.1f})\n"
            f"Window={window} samp, Step={step} samp, Coinc={coincidence_window:.1f} ns, "
            f"Trials={n_trials}\n"
            "FC(MHz)\tSNR50"
        )
        np.savetxt(f'plots/{save_filename}_filter_opt.txt', data, fmt='%.6g', header=header)
        print(f"\nSaved: plots/{save_filename}_filter_opt.npy and .txt")

    return opt_results


def optimize_step(phi, theta, psi, thresholds_per_step, coincidence_window,
                  snr_grid, n_trials, window, step_grid,
                  save_filename=None, antennas=None):
    """
    Scan over power-sum stride sizes (step) and find SNR at 50% coincidence efficiency.

    The window size is held fixed while the stride (number of samples between
    successive power-sum frames) is varied.  A smaller stride increases the
    number of frames per event and the time resolution of the coincidence check,
    at the cost of higher computational load.

    Parameters
    ----------
    phi, theta, psi : float
        Signal arrival angles (deg).
    thresholds_per_step : dict or float
        Per-step thresholds: dict mapping step (samples, int) -> threshold float,
        or a single float used for all steps. Using a per-step dict is strongly
        preferred because smaller strides produce more frames/sec and require a
        higher threshold to hold the same false rate in Hz.
    coincidence_window : float
        Coincidence time window (ns).
    snr_grid : array
        Per-antenna SNR values to test.
    n_trials : int
        Trials per SNR point.
    window : int
        Power-sum window size (samples); held fixed.
    step_grid : tuple
        (min_samples, max_samples, scan_step) for the stride scan range.
    save_filename : str or None
        Base filename; saves plots/<save_filename>_step_opt.npy/.txt.
    antennas : list of int or None
        Physical antenna indices (default: all four).

    Returns
    -------
    dict with keys:
        steps, steps_ns, snr50, window, window_ns,
        phi, theta, psi, thresholds_per_step, coincidence_window, n_trials
    """
    step_min, step_max, step_scan_step = step_grid
    steps = np.arange(step_min, step_max + 1, step_scan_step, dtype=int)
    sample_ns = geometry.ritc_sample_step
    window_ns = window * sample_ns

    # Drop strides larger than the window — those leave gaps between frames.
    invalid = steps > window
    if invalid.any():
        print(f"WARNING: skipping {invalid.sum()} stride value(s) larger than "
              f"window ({window} samples): {steps[invalid].tolist()}")
        steps = steps[~invalid]
    if len(steps) == 0:
        raise ValueError(
            f"No valid stride sizes remain after removing values > window ({window}). "
            f"Decrease --step-max or increase --window.")

    print(f"\n=== Power-Sum Stride Optimization ===")
    print(f"Direction: phi={phi:.1f}°, theta={theta:.1f}°, psi={psi:.1f}°")
    if isinstance(thresholds_per_step, dict):
        print(f"Thresholds: per-step (calibrated)  |  Coincidence window: {coincidence_window:.1f} ns")
    else:
        print(f"Threshold: {float(thresholds_per_step):.4f}  |  Coincidence window: {coincidence_window:.1f} ns")
    print(f"Window: {window} samples ({window_ns:.1f} ns)  [fixed]")
    print(f"Stride range: {steps[0]}–{steps[-1]} samples "
          f"({steps[0]*sample_ns:.1f}–{steps[-1]*sample_ns:.1f} ns)  "
          f"[window fixed at {window} samples = {window_ns:.1f} ns]")
    print(f"SNR points: {len(snr_grid)}  |  Trials: {n_trials}")

    snr50_list = []

    for step in steps:
        step_ns = step * sample_ns
        if isinstance(thresholds_per_step, dict):
            thr = thresholds_per_step.get(
                int(step),
                next(iter(thresholds_per_step.values()))
            )
        else:
            thr = float(thresholds_per_step)
        print(f"\n--- Stride = {step} samples ({step_ns:.1f} ns), threshold={thr:.4f} ---")

        res = snr_scan.run_snr_scan_dualpol(
            phi=phi,
            theta=theta,
            psi=psi,
            threshold_lhcp=thr,
            threshold_rhcp=thr,
            coincidence_window=coincidence_window,
            snr_grid=snr_grid,
            n_trials=n_trials,
            window=window,
            step=step,
            save_filename=None,
            antennas=antennas,
        )

        snr_50 = _interpolate_snr50(res['snr'], res['eff_coinc'])
        snr50_list.append(snr_50)
        print(f"  → SNR50 = {snr_50:.3f}")

    steps_ns = steps * sample_ns

    opt_results = {
        'steps':               steps,
        'steps_ns':            steps_ns,
        'snr50':               np.array(snr50_list),
        'window':              window,
        'window_ns':           window_ns,
        'phi':                 phi,
        'theta':               theta,
        'psi':                 psi,
        'thresholds_per_step': thresholds_per_step,
        'coincidence_window':  coincidence_window,
        'n_trials':            n_trials,
        'snr_grid':            snr_grid,
        'antennas':            antennas,
    }

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        np.save(f'plots/{save_filename}_step_opt.npy', opt_results)
        thr_note = ('per-step' if isinstance(thresholds_per_step, dict)
                    else f'{float(thresholds_per_step):.3f}')
        data = np.column_stack([steps, steps_ns, snr50_list])
        header = (
            f"Stride optimization at (phi={phi:.1f}, theta={theta:.1f}, psi={psi:.1f})\n"
            f"Threshold={thr_note}, Window={window} samp ({window_ns:.1f} ns), "
            f"Coinc={coincidence_window:.1f} ns, Trials={n_trials}\n"
            "Stride(samples)\tStride(ns)\tSNR50"
        )
        np.savetxt(f'plots/{save_filename}_step_opt.txt', data, fmt='%.6g', header=header)
        print(f"\nSaved: plots/{save_filename}_step_opt.npy and .txt")

    return opt_results


def scan_falserate_vs_window(window_grid, step, phi=0.0, theta=0.0,
                             noise_dir='noise', duration_sec=None,
                             chunk_duration_sec=0.05,
                             threshold_grid=None,
                             coincidence_window=100.0, n_beams=100,
                             save_filename=None):
    """
    Scan false-alarm rate vs. threshold for multiple power-sum window sizes,
    sharing the expensive filter+beamform+circular-pol work across all windows.

    For each chunk of noise:
      1. Apply Shannon-Whitaker filter
      2. Coherent-sum + circular-pol conversion  →  LHCP, RHCP arrays
      3. Apply powerSum at every window size in one pass

    Then for each window compute::

        rate(T) = (power >= T).sum() / obs_time
        rate_coinc = 2 * rate_lhcp * rate_rhcp * tau
        rate_global = n_beams * rate_coinc

    Parameters
    ----------
    window_grid : tuple
        (min_samples, max_samples, samples_per_step)
    step : int
        Power-sum sliding step (samples); held fixed.
    phi, theta : float
        Beamforming direction (deg).
    noise_dir : str
        Directory containing dualpol_noise_hpol/vpol.npy and metadata.
    duration_sec : float or None
        How much noise to process.  None = use all available.
    chunk_duration_sec : float
        Process in chunks of this size (default 50 ms).
    threshold_grid : array or None
        Threshold values to scan.  Defaults to 0.5 … 10.0 step 0.05.
    coincidence_window : float
        Coincidence time window (ns).
    n_beams : int
        Number of active beams for global rate estimate.
    save_filename : str or None
        Base name for .npz output saved into plots/.

    Returns
    -------
    dict with keys:
        windows, threshold_grid, results, phi, theta, step,
        coincidence_window, n_beams
        results[w] has: rate_lhcp, rate_rhcp, rate_coinc_per_beam,
                        rate_global, obs_time, n_frames
    """
    from pathlib import Path
    from payload_signal import getRemappedDelays
    import math

    win_min, win_max, win_step = window_grid
    windows = np.arange(win_min, win_max + 1, win_step, dtype=int)
    max_win = int(windows.max())

    if threshold_grid is None:
        threshold_grid = np.arange(0.5, 10.05, 0.05)
    threshold_grid = np.asarray(threshold_grid)

    # ---- load noise ----
    noise_path = Path(noise_dir)
    meta_file = noise_path / 'dualpol_noise_metadata.npy'
    if meta_file.exists():
        meta = np.load(meta_file, allow_pickle=True).item()
        sample_rate_hz = meta['sample_rate_GHz'] * 1e9
        print(f"Noise metadata: {meta['duration_sec']} sec at {meta['sample_rate_GHz']} GHz")
    else:
        import tools.CoRaLs_geometry as _cg
        sample_rate_hz = _cg.ritc_sample_rate * 1e9
        print(f"No metadata, using default {sample_rate_hz/1e9} GHz")

    h_mmap = np.load(noise_path / 'dualpol_noise_hpol.npy', mmap_mode='r')
    v_mmap = np.load(noise_path / 'dualpol_noise_vpol.npy', mmap_mode='r')
    n_ant = h_mmap.shape[0]

    total_available = h_mmap.shape[1]
    if duration_sec is None:
        total_samples = total_available
    else:
        total_samples = min(int(duration_sec * sample_rate_hz), total_available)

    chunk_samples = int(chunk_duration_sec * sample_rate_hz)
    n_chunks = math.ceil(total_samples / chunk_samples)
    dt_ns = 1.0 / sample_rate_hz * 1e9

    delays = getRemappedDelays(phi, theta, antennas=list(range(n_ant)))

    print(f"\n=== False-Rate Window Scan ===")
    print(f"Direction: phi={phi}°, theta={theta}°")
    print(f"Windows: {win_min}–{win_max} samples, step {win_step}  ({len(windows)} points)")
    print(f"Power-step: {step} samples | Coincidence: {coincidence_window} ns | N_beams: {n_beams}")
    print(f"Noise: {total_samples/sample_rate_hz:.2f} sec in {n_chunks} chunks")

    # accumulators: list of arrays per window
    power_lhcp_acc = {w: [] for w in windows}
    power_rhcp_acc = {w: [] for w in windows}
    # tails for cross-chunk continuity (shared, trimmed to max_win + 2*step)
    tail_lhcp = None
    tail_rhcp = None
    tail_len = max_win + step * 2

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_samples
        end   = min(start + chunk_samples, total_samples)

        h_chunk = np.array(h_mmap[:, start:end]).astype(np.float64)
        v_chunk = np.array(v_mmap[:, start:end]).astype(np.float64)

        # filter
        h_chunk = filters.apply_Shannon_Whitaker_filter(h_chunk, fs=sample_rate_hz)
        v_chunk = filters.apply_Shannon_Whitaker_filter(v_chunk, fs=sample_rate_hz)

        # coherent sum + circular pol
        tb = np.arange(end - start) * dt_ns
        lhcp, rhcp, _ = trigger.coherentSum_dualpol(
            h_chunk, v_chunk, tb, delays,
            downsample=False, output='circular', apply_filter=False,
        )

        # prepend tail for continuity
        if tail_lhcp is not None:
            lhcp = np.concatenate([tail_lhcp, lhcp])
            rhcp = np.concatenate([tail_rhcp, rhcp])

        # apply powerSum at every window size
        for w in windows:
            if len(lhcp) >= w:
                pl, _ = trigger.powerSum(lhcp, window=int(w), step=step)
                pr, _ = trigger.powerSum(rhcp, window=int(w), step=step)
                power_lhcp_acc[w].append(pl)
                power_rhcp_acc[w].append(pr)

        # update shared tail
        if len(lhcp) >= tail_len:
            tail_lhcp = lhcp[-tail_len:]
            tail_rhcp = rhcp[-tail_len:]
        else:
            tail_lhcp = lhcp
            tail_rhcp = rhcp

        print(f"  Chunk {chunk_idx+1}/{n_chunks}: {end-start:,} samples  "
              f"({(end/sample_rate_hz):.2f} s)")

    # ---- concatenate & compute rates ----
    frame_dt_s = step / sample_rate_hz  # time per power frame
    results = {}
    for w in windows:
        if not power_lhcp_acc[w]:
            continue
        pl = np.concatenate(power_lhcp_acc[w])
        pr = np.concatenate(power_rhcp_acc[w])
        obs_time = len(pl) * frame_dt_s

        rate_l = np.array([(pl >= thr).sum() for thr in threshold_grid], dtype=float) / obs_time
        rate_r = np.array([(pr >= thr).sum() for thr in threshold_grid], dtype=float) / obs_time
        tau_s  = coincidence_window * 1e-9
        rate_coinc  = 2.0 * rate_l * rate_r * tau_s
        rate_global = n_beams * rate_coinc

        results[int(w)] = {
            'rate_lhcp':           rate_l,
            'rate_rhcp':           rate_r,
            'rate_coinc_per_beam': rate_coinc,
            'rate_global':         rate_global,
            'obs_time':            obs_time,
            'n_frames':            len(pl),
        }
        print(f"  Window {w:4d} samples ({w*dt_ns:.1f} ns): "
              f"{len(pl):,} frames, {obs_time:.2f} s obs time")

    scan = {
        'windows':           windows,
        'threshold_grid':    threshold_grid,
        'results':           results,
        'phi':               phi,
        'theta':             theta,
        'step':              step,
        'coincidence_window': coincidence_window,
        'n_beams':           n_beams,
        'sample_rate_hz':    sample_rate_hz,
    }

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_falserate.npz'
        np.savez(out, **{k: v for k, v in scan.items() if k != 'results'},
                 **{f'rate_global_{w}': results[w]['rate_global'] for w in results},
                 **{f'rate_lhcp_{w}':  results[w]['rate_lhcp']  for w in results})
        print(f"\nSaved: {out}")

    return scan


def plot_falserate_vs_window(scan, threshold_marker=None, rate_target=None,
                             save_filename=None):
    """
    Semilog plot of global coincidence rate vs. threshold for each window size.

    Parameters
    ----------
    scan : dict
        Output of scan_falserate_vs_window().
    threshold_marker : float or None
        Draw a vertical line at this threshold value.
    rate_target : float or None
        Draw a horizontal line at this global rate (Hz).
    save_filename : str or None
        Save to plots/<save_filename>_falserate.png.

    Returns
    -------
    (fig, ax)
    """
    windows        = scan['windows']
    thr_grid       = scan['threshold_grid']
    results        = scan['results']
    phi            = scan['phi']
    theta          = scan['theta']
    step           = scan['step']
    n_beams        = scan['n_beams']
    sample_ns      = geometry.ritc_sample_step
    step_ns        = step * sample_ns

    cmap   = plt.cm.viridis
    colors = [cmap(i / max(len(windows) - 1, 1)) for i in range(len(windows))]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    ax_global = axes[0]
    ax_single = axes[1]

    for i, w in enumerate(windows):
        if w not in results:
            continue
        r = results[w]
        label = f'{w} samp ({w*sample_ns:.1f} ns)'
        c = colors[i]

        # global coincidence rate
        valid = r['rate_global'] > 0
        if valid.any():
            ax_global.semilogy(thr_grid[valid], r['rate_global'][valid],
                               '-', color=c, linewidth=1.4, label=label)

        # single-channel LHCP rate
        valid_l = r['rate_lhcp'] > 0
        if valid_l.any():
            ax_single.semilogy(thr_grid[valid_l], r['rate_lhcp'][valid_l],
                               '-', color=c, linewidth=1.4, label=label)

    for ax in axes:
        if threshold_marker is not None:
            ax.axvline(threshold_marker, color='red', linestyle='--',
                       linewidth=1.5, label=f'T = {threshold_marker}')
        if rate_target is not None:
            ax.axhline(rate_target, color='orange', linestyle=':',
                       linewidth=1.5, label=f'Target {rate_target} Hz')
        ax.set_xlabel('Threshold', fontsize=12)
        ax.grid(True, which='both', alpha=0.25)

    ax_global.set_ylabel('Global Coincidence Rate (Hz)', fontsize=12)
    ax_global.set_title(
        f'Global Rate vs Threshold (N_beams={n_beams})\n'
        f'φ={phi}°, θ={theta}°  |  step={step} samp ({step_ns:.1f} ns)', fontsize=11)
    ax_global.legend(fontsize=7, ncol=2, loc='upper right')

    ax_single.set_ylabel('Single-Channel Rate (Hz)', fontsize=12)
    ax_single.set_title(
        f'LHCP Per-Channel Rate vs Threshold\n'
        f'φ={phi}°, θ={theta}°  |  step={step} samp ({step_ns:.1f} ns)', fontsize=11)
    ax_single.legend(fontsize=7, ncol=2, loc='upper right')

    plt.tight_layout()

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_falserate.png'
        plt.savefig(out, dpi=150)
        print(f"Saved: {out}")

    plt.show()
    return fig, axes


def plot_window_optimization(opt_results, save_filename=None):
    """
    Plot SNR50 vs power-sum window size from an optimize_windows run.

    Parameters:
    -----------
    opt_results : dict
        Output of optimize_windows()
    save_filename : str, optional
        Base filename; saves <save_filename>_window_opt.png

    Returns:
    --------
    (fig, ax)
    """
    windows_ns = opt_results['windows_ns']
    snr50      = opt_results['snr50']
    phi        = opt_results['phi']
    theta      = opt_results['theta']
    psi        = opt_results['psi']
    _tpw       = opt_results.get('thresholds_per_window', opt_results.get('threshold'))
    thr_label  = ('per-window' if isinstance(_tpw, dict)
                  else f'{float(_tpw):.3f}')
    step       = opt_results['step']
    sample_ns  = geometry.ritc_sample_step
    step_ns    = step * sample_ns

    fig, ax = plt.subplots(figsize=(9, 5))

    # Main curve
    valid = ~np.isnan(snr50)
    ax.plot(windows_ns[valid], snr50[valid], 'o-',
            color='steelblue', linewidth=2, markersize=6, label='SNR @ 50% efficiency')

    # Mark the best (minimum SNR50) window
    if valid.any():
        best_idx = int(np.nanargmin(snr50))
        ax.plot(windows_ns[best_idx], snr50[best_idx], '*',
                color='red', markersize=14, zorder=5,
                label=f'Best: {windows_ns[best_idx]:.1f} ns  (SNR50 = {snr50[best_idx]:.2f})')
        ax.axvline(windows_ns[best_idx], color='red', linestyle=':', alpha=0.4, linewidth=1)

    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3, linewidth=1)

    ax.set_xlabel('Power Sum Window (ns)', fontsize=12)
    ax.set_ylabel('SNR at 50% Coincidence Efficiency', fontsize=12)
    ax.set_title(
        f'Window Optimization: φ={phi:.1f}°, θ={theta:.1f}°, ψ={psi:.1f}°\n'
        f'Threshold = {thr_label}  |  Step = {step} samples ({step_ns:.1f} ns)',
        fontsize=12,
    )
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_window_opt.png'
        plt.savefig(out, dpi=150)
        print(f"Saved: {out}")

    plt.show()
    return fig, ax


def plot_filter_optimization(opt_results, save_filename=None):
    """
    Plot SNR50 vs filter cutoff frequency from an optimize_filter_freq run.

    Parameters
    ----------
    opt_results : dict
        Output of optimize_filter_freq().
    save_filename : str or None
        Base filename; saves plots/<save_filename>_filter_opt.png.

    Returns
    -------
    (fig, ax)
    """
    fc_mhz  = opt_results['fc_mhz']
    snr50   = opt_results['snr50']
    phi     = opt_results['phi']
    theta   = opt_results['theta']
    psi     = opt_results['psi']
    window  = opt_results['window']
    step    = opt_results['step']
    sample_ns = geometry.ritc_sample_step
    window_ns = window * sample_ns
    step_ns   = step   * sample_ns

    fig, ax = plt.subplots(figsize=(9, 5))

    valid = ~np.isnan(snr50)
    ax.plot(fc_mhz[valid], snr50[valid], 'o-',
            color='darkorchid', linewidth=2, markersize=7,
            label='SNR @ 50% efficiency')

    if valid.any():
        best_idx = int(np.nanargmin(snr50))
        ax.plot(fc_mhz[best_idx], snr50[best_idx], '*',
                color='red', markersize=14, zorder=5,
                label=f'Best: {fc_mhz[best_idx]:.0f} MHz  '
                      f'(SNR50 = {snr50[best_idx]:.3f})')
        ax.axvline(fc_mhz[best_idx], color='red', linestyle=':', alpha=0.4, linewidth=1)

    ax.set_xlabel('Filter Cutoff Frequency (MHz)', fontsize=12)
    ax.set_ylabel('SNR at 50% Coincidence Efficiency', fontsize=12)
    ax.set_title(
        f'Filter Frequency Optimization: φ={phi:.1f}°, θ={theta:.1f}°, ψ={psi:.1f}°\n'
        f'Window = {window} samp ({window_ns:.1f} ns)  |  '
        f'Step = {step} samp ({step_ns:.1f} ns)',
        fontsize=12,
    )
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_filter_opt.png'
        plt.savefig(out, dpi=150)
        print(f"Saved: {out}")

    plt.show()
    return fig, ax


def plot_step_optimization(opt_results, save_filename=None):
    """
    Plot SNR50 vs power-sum stride size from an optimize_step run.

    Parameters
    ----------
    opt_results : dict
        Output of optimize_step().
    save_filename : str or None
        Base filename; saves plots/<save_filename>_step_opt.png.

    Returns
    -------
    (fig, ax)
    """
    steps_ns  = opt_results['steps_ns']
    snr50     = opt_results['snr50']
    phi       = opt_results['phi']
    theta     = opt_results['theta']
    psi       = opt_results['psi']
    _tps      = opt_results.get('thresholds_per_step', opt_results.get('threshold'))
    thr_label = ('per-step' if isinstance(_tps, dict) else f'{float(_tps):.3f}')
    window    = opt_results['window']
    window_ns = opt_results['window_ns']

    fig, ax = plt.subplots(figsize=(9, 5))

    valid = ~np.isnan(snr50)
    ax.plot(steps_ns[valid], snr50[valid], 's-',
            color='steelblue', linewidth=2, markersize=7,
            label='SNR @ 50% efficiency')

    if valid.any():
        best_idx = int(np.nanargmin(snr50))
        ax.plot(steps_ns[best_idx], snr50[best_idx], '*',
                color='red', markersize=14, zorder=5,
                label=f'Best: {steps_ns[best_idx]:.1f} ns  '
                      f'(SNR50 = {snr50[best_idx]:.3f})')
        ax.axvline(steps_ns[best_idx], color='red', linestyle=':', alpha=0.4, linewidth=1)

    ax.set_xlabel('Power-Sum Stride (ns)', fontsize=12)
    ax.set_ylabel('SNR at 50% Coincidence Efficiency', fontsize=12)
    ax.set_title(
        f'Stride Optimization: φ={phi:.1f}°, θ={theta:.1f}°, ψ={psi:.1f}°\n'
        f'Threshold = {thr_label}  |  Window = {window} samp ({window_ns:.1f} ns)',
        fontsize=12,
    )
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_step_opt.png'
        plt.savefig(out, dpi=150)
        print(f"Saved: {out}")

    plt.show()
    return fig, ax


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Optimization scans for CoRaLS trigger parameters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Optimization mode
    parser.add_argument('--mode',
                        choices=['windows', 'step', 'filter', 'falserate'],
                        default='windows',
                        help='Optimization axis: windows=scan window size, '
                             'step=scan power-sum stride, '
                             'filter=scan filter cutoff frequency, '
                             'falserate=false-rate vs threshold curves')

    # ---------- Wave direction ----------
    parser.add_argument('--phi',   type=float, default=0.0,
                        help='Azimuth angle (deg)')
    parser.add_argument('--theta', type=float, default=0.0,
                        help='Elevation angle (deg; -90=nadir, 0=horizon, +90=zenith)')
    parser.add_argument('--psi',   type=float, default=45.0,
                        help='Polarization angle (deg; 0=H-pol, 90=V-pol)')

    # ---------- Threshold / coincidence ----------
    parser.add_argument('--threshold', type=float, default=None,
                        help='Normalized power threshold (LHCP and RHCP). '
                             'Overridden by the first entry in --threshold-json if provided.')
    parser.add_argument('--coincidence-window', type=float, default=None,
                        help='Coincidence time window (ns). '
                             'Overridden by coincidence_window_ns in --threshold-json if present.')

    # ---------- Threshold / n_beams from JSON ----------
    parser.add_argument('--threshold-json', type=str, default=None,
                        help='Path to threshold_analysis_dualpol.json. Reads n_beams, '
                             'coincidence_window_ns, and the threshold for --threshold-index. '
                             'CLI --threshold and --coincidence-window take precedence if given.')
    parser.add_argument('--threshold-index', type=int, default=None,
                        help='Index into the "thresholds" list in --threshold-json. '
                             'Omit (or 0) to use the first entry.')
    parser.add_argument('--n-beams', type=int, default=None,
                        help='Number of active beams for global rate estimate. '
                             'Read from --threshold-json if not specified here.')

    # ---------- SNR scan grid ----------
    parser.add_argument('--snr-min',  type=float, default=0.2, help='Minimum SNR')
    parser.add_argument('--snr-max',  type=float, default=6.0, help='Maximum SNR')
    parser.add_argument('--snr-step', type=float, default=0.2, help='SNR step size')
    parser.add_argument('--trials',   type=int,   default=500,
                        help='Trials per SNR point')

    # ---------- Window scan parameters ----------
    parser.add_argument('--window-min',  type=int, default=40,
                        help='Minimum power-sum window (samples) [windows mode]')
    parser.add_argument('--window-max',  type=int, default=320,
                        help='Maximum power-sum window (samples) [windows mode]')
    parser.add_argument('--window-step', type=int, default=40,
                        help='Step size for the window scan (samples) [windows mode]')
    parser.add_argument('--step',        type=int, default=40,
                        help='Power-sum stride (samples); fixed during windows/filter scans')
    parser.add_argument('--window-threshold-jsons', type=str, nargs='+', default=None,
                        help='Per-window threshold JSON files, one per window value in the '
                             'scan grid (in order, after filtering out windows < stride). '
                             'Each JSON is loaded for --threshold-index. '
                             'If omitted, --threshold-json / --threshold is used for all windows.')

    # ---------- Stride (step) scan parameters ----------
    parser.add_argument('--step-min',       type=int, default=10,
                        help='Minimum stride (samples) [step mode]')
    parser.add_argument('--step-max',       type=int, default=160,
                        help='Maximum stride (samples) [step mode]')
    parser.add_argument('--step-scan-step', type=int, default=10,
                        help='Increment for stride scan (samples) [step mode]')
    parser.add_argument('--window',         type=int, default=160,
                        help='Fixed power-sum window (samples) used in step/filter modes')
    parser.add_argument('--step-threshold-jsons', type=str, nargs='+', default=None,
                        help='Per-step threshold JSON files, one per stride value in the '
                             'scan grid (in order, after filtering out strides > window). '
                             'Each JSON is loaded for --threshold-index. '
                             'If omitted, --threshold-json / --threshold is used for all strides.')

    # ---------- Filter frequency scan parameters ----------
    parser.add_argument('--filter-freqs-ghz', type=float, nargs='+',
                        default=[1.5, 1.25, 1.0, 0.75],
                        help='Filter cutoff frequencies in GHz to scan [filter mode]. '
                             'The first value is the baseline (no second filter).')
    parser.add_argument('--filter-threshold-jsons', type=str, nargs='+', default=None,
                        help='Threshold JSON files for each non-baseline filter frequency, '
                             'in the same order as --filter-freqs-ghz[1:]. '
                             'If omitted, the reference threshold is used for all fc values.')

    # ---------- Antenna selection ----------
    parser.add_argument('--antennas', type=int, nargs='+', default=None,
                        help='Physical antenna indices (default: all four)')

    # ---------- False-rate scan parameters ----------
    parser.add_argument('--noise-dir',    type=str,   default='noise',
                        help='Directory containing dualpol noise files')
    parser.add_argument('--duration',     type=float, default=None,
                        help='Seconds of noise to process (default: all available)')
    parser.add_argument('--chunk-duration', type=float, default=0.05,
                        help='Chunk size for noise processing (sec)')
    parser.add_argument('--thr-min',      type=float, default=0.5,
                        help='Minimum threshold to scan [falserate mode]')
    parser.add_argument('--thr-max',      type=float, default=10.0,
                        help='Maximum threshold to scan [falserate mode]')
    parser.add_argument('--thr-step',     type=float, default=0.05,
                        help='Threshold step [falserate mode]')
    parser.add_argument('--threshold-marker', type=float, default=None,
                        help='Mark this threshold on the falserate plot')
    parser.add_argument('--rate-target',  type=float, default=None,
                        help='Draw horizontal line at this global rate (Hz)')

    # ---------- Output ----------
    parser.add_argument('--save',    type=str, default='optimization',
                        help='Base filename for saved results and plots')
    parser.add_argument('--no-plot', action='store_true', help='Skip plotting')

    args = parser.parse_args()

    # ---- Resolve n_beams, coincidence_window, and threshold from JSON ----
    n_beams_resolved           = 63      # CoRaLS default
    coincidence_window_resolved = 100.0  # ns default
    threshold_resolved          = None

    if args.threshold_json is not None:
        with open(args.threshold_json, 'r') as _f:
            _thr_data = json.load(_f)
        n_beams_resolved            = _thr_data.get('n_beams', n_beams_resolved)
        coincidence_window_resolved = _thr_data.get('coincidence_window_ns',
                                                     coincidence_window_resolved)
        _entries = _thr_data['thresholds']
        _idx = args.threshold_index if args.threshold_index is not None else 0
        threshold_resolved = _entries[_idx]['threshold_lhcp']
        print(f"Loaded from {args.threshold_json}: "
              f"n_beams={n_beams_resolved}, "
              f"coinc_window={coincidence_window_resolved:.1f} ns, "
              f"threshold={threshold_resolved:.4f} (entry {_idx})")

    # CLI overrides take precedence
    if args.n_beams is not None:
        n_beams_resolved = args.n_beams
    if args.coincidence_window is not None:
        coincidence_window_resolved = args.coincidence_window
    if args.threshold is not None:
        threshold_resolved = args.threshold
    if threshold_resolved is None:
        threshold_resolved = 3.12  # fallback if no JSON and no --threshold
        print(f"No --threshold or --threshold-json provided; using default threshold="
              f"{threshold_resolved}")

    snr_grid = np.arange(args.snr_min, args.snr_max + 0.5 * args.snr_step, args.snr_step)

    # -------------------------------------------------------------------------
    if args.mode == 'windows':
        # Build per-window threshold dict if calibrated JSONs are provided.
        # Replicate the same window grid + guard logic used inside optimize_windows
        # so indices align with --window-threshold-jsons order.
        _win_all = np.arange(args.window_min, args.window_max + 1, args.window_step, dtype=int)
        _valid_wins = _win_all[_win_all >= args.step]
        if args.window_threshold_jsons:
            if len(args.window_threshold_jsons) != len(_valid_wins):
                parser.error(
                    f'--window-threshold-jsons: expected {len(_valid_wins)} files '
                    f'(one per valid window in {_valid_wins.tolist()}), '
                    f'got {len(args.window_threshold_jsons)}')
            _idx = args.threshold_index if args.threshold_index is not None else 0
            thresholds_per_window = {}
            for win, json_path in zip(_valid_wins, args.window_threshold_jsons):
                with open(json_path, 'r') as _f:
                    _jdata = json.load(_f)
                thresholds_per_window[int(win)] = _jdata['thresholds'][_idx]['threshold_lhcp']
                print(f"  window={win}: threshold={thresholds_per_window[int(win)]:.4f} "
                      f"(from {json_path}, entry {_idx})")
        else:
            thresholds_per_window = threshold_resolved

        opt = optimize_windows(
            phi=args.phi,
            theta=args.theta,
            psi=args.psi,
            thresholds_per_window=thresholds_per_window,
            coincidence_window=coincidence_window_resolved,
            snr_grid=snr_grid,
            n_trials=args.trials,
            window_grid=(args.window_min, args.window_max, args.window_step),
            step=args.step,
            save_filename=args.save,
            antennas=args.antennas,
        )

        if not args.no_plot:
            plot_window_optimization(opt, save_filename=args.save)

        best_idx = int(np.nanargmin(opt['snr50']))
        print(f"\n=== Optimization Complete ===")
        print(f"Best window: {opt['windows'][best_idx]} samples "
              f"({opt['windows_ns'][best_idx]:.1f} ns)  "
              f"→  SNR50 = {opt['snr50'][best_idx]:.3f}")

    # -------------------------------------------------------------------------
    elif args.mode == 'step':
        # Build per-step threshold dict if calibrated JSONs are provided.
        # Replicate the same step grid + guard logic used inside optimize_step
        # so indices align with --step-threshold-jsons order.
        _step_all = np.arange(args.step_min, args.step_max + 1, args.step_scan_step, dtype=int)
        _valid_steps = _step_all[_step_all <= args.window]
        if args.step_threshold_jsons:
            if len(args.step_threshold_jsons) != len(_valid_steps):
                parser.error(
                    f'--step-threshold-jsons: expected {len(_valid_steps)} files '
                    f'(one per valid stride in {_valid_steps.tolist()}), '
                    f'got {len(args.step_threshold_jsons)}')
            _idx = args.threshold_index if args.threshold_index is not None else 0
            thresholds_per_step = {}
            for s, json_path in zip(_valid_steps, args.step_threshold_jsons):
                with open(json_path, 'r') as _f:
                    _jdata = json.load(_f)
                thresholds_per_step[int(s)] = _jdata['thresholds'][_idx]['threshold_lhcp']
                print(f"  step={s}: threshold={thresholds_per_step[int(s)]:.4f} "
                      f"(from {json_path}, entry {_idx})")
        else:
            thresholds_per_step = threshold_resolved

        opt = optimize_step(
            phi=args.phi,
            theta=args.theta,
            psi=args.psi,
            thresholds_per_step=thresholds_per_step,
            coincidence_window=coincidence_window_resolved,
            snr_grid=snr_grid,
            n_trials=args.trials,
            window=args.window,
            step_grid=(args.step_min, args.step_max, args.step_scan_step),
            save_filename=args.save,
            antennas=args.antennas,
        )

        if not args.no_plot:
            plot_step_optimization(opt, save_filename=args.save)

        best_idx = int(np.nanargmin(opt['snr50']))
        print(f"\n=== Optimization Complete ===")
        print(f"Best stride: {opt['steps'][best_idx]} samples "
              f"({opt['steps_ns'][best_idx]:.1f} ns)  "
              f"→  SNR50 = {opt['snr50'][best_idx]:.3f}")

    # -------------------------------------------------------------------------
    elif args.mode == 'filter':
        fc_hz_list = [f * 1e9 for f in args.filter_freqs_ghz]

        # Build per-fc threshold map
        thresholds_per_fc = {fc_hz_list[0]: threshold_resolved}
        if args.filter_threshold_jsons:
            for fc, json_path in zip(fc_hz_list[1:], args.filter_threshold_jsons):
                with open(json_path, 'r') as _fj:
                    _jdata = json.load(_fj)
                _idx = args.threshold_index if args.threshold_index is not None else 0
                thresholds_per_fc[fc] = _jdata['thresholds'][_idx]['threshold_lhcp']
                print(f"  fc={fc/1e6:.0f} MHz: threshold={thresholds_per_fc[fc]:.4f} "
                      f"(from {json_path})")
        else:
            for fc in fc_hz_list[1:]:
                thresholds_per_fc[fc] = threshold_resolved
            if len(fc_hz_list) > 1:
                print("NOTE: No --filter-threshold-jsons provided; using reference threshold "
                      f"({threshold_resolved:.4f}) for all filter frequencies. "
                      "Re-calibrate thresholds per fc with generate_threshold_curves_dualpol.py "
                      "for accurate results.")

        opt = optimize_filter_freq(
            phi=args.phi,
            theta=args.theta,
            psi=args.psi,
            fc_hz_list=fc_hz_list,
            thresholds_per_fc=thresholds_per_fc,
            coincidence_window=coincidence_window_resolved,
            snr_grid=snr_grid,
            n_trials=args.trials,
            window=args.window,
            step=args.step,
            save_filename=args.save,
            antennas=args.antennas,
        )

        if not args.no_plot:
            plot_filter_optimization(opt, save_filename=args.save)

        best_idx = int(np.nanargmin(opt['snr50']))
        print(f"\n=== Optimization Complete ===")
        print(f"Best filter frequency: {opt['fc_mhz'][best_idx]:.0f} MHz  "
              f"→  SNR50 = {opt['snr50'][best_idx]:.3f}")

    # -------------------------------------------------------------------------
    elif args.mode == 'falserate':
        thr_grid = np.arange(args.thr_min, args.thr_max + 0.5 * args.thr_step, args.thr_step)

        scan = scan_falserate_vs_window(
            window_grid=(args.window_min, args.window_max, args.window_step),
            step=args.step,
            phi=args.phi,
            theta=args.theta,
            noise_dir=args.noise_dir,
            duration_sec=args.duration,
            chunk_duration_sec=args.chunk_duration,
            threshold_grid=thr_grid,
            coincidence_window=coincidence_window_resolved,
            n_beams=n_beams_resolved,
            save_filename=args.save,
        )

        if not args.no_plot:
            plot_falserate_vs_window(
                scan,
                threshold_marker=args.threshold_marker,
                rate_target=args.rate_target,
                save_filename=args.save,
            )

