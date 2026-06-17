"""
plot_scope_data.py
------------------
Discover, display, and optionally compare Rohde & Schwarz oscilloscope
reference-curve exports with the CoRaLs antenna impulse response.

Typical usage
-------------
# Plot all waveforms from the default data/ directory:
    python plot_scope_data.py

# Also overlay the impulse response from payload_signal:
    python plot_scope_data.py --impulse impulse/corals_impulse_sci.txt

# Specify a different data directory:
    python plot_scope_data.py --data-dir /path/to/data

# Save the figure instead of (or in addition to) showing it:
    python plot_scope_data.py --save scope_overview
"""

import argparse
import os
import glob
import re
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ---------------------------------------------------------------------------
# Metadata parser
# ---------------------------------------------------------------------------

def parse_metadata(csv_path: str) -> dict:
    """
    Parse the R&S Rohde CPgReferenceCurveAttributes header file.
    Returns a dict with float values where possible.
    """
    meta = {}
    with open(csv_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.rstrip(':').split(':')
            if len(parts) >= 2:
                key = parts[0]
                val = parts[1]
                try:
                    meta[key] = float(val)
                except ValueError:
                    meta[key] = val
    # Convenient derived quantities
    meta['XStart_ns']      = meta.get('XStart', 0.0) * 1e9
    meta['XStop_ns']       = meta.get('XStop',  0.0) * 1e9
    meta['Resolution_ns']  = meta.get('Resolution', 4e-12) * 1e9
    meta['RecordLength']   = int(meta.get('RecordLength', 1000))
    return meta


# ---------------------------------------------------------------------------
# Waveform loader
# ---------------------------------------------------------------------------

def load_waveform(wfm_path: str, meta: dict):
    """
    Load raw voltage samples from a .Wfm.csv file.

    Returns
    -------
    time_ns : ndarray  — time axis in nanoseconds
    voltage : ndarray  — voltage in Volts
    """
    voltage = np.loadtxt(wfm_path)
    n = len(voltage)
    x_start = meta.get('XStart', 0.0)
    x_stop  = meta.get('XStop', 0.0)
    time_s  = np.linspace(x_start, x_stop, n)
    time_ns = time_s * 1e9
    return time_ns, voltage


# ---------------------------------------------------------------------------
# Pair discovery
# ---------------------------------------------------------------------------

def discover_pairs(data_dir: str) -> list[dict]:
    """
    Walk *data_dir* and return one dict per (header, waveform) pair, sorted
    by filename so curves from the same session appear together.

    Each dict has keys: 'stem', 'header', 'wfm', 'channel', 'timestamp'
    """
    header_files = sorted(glob.glob(os.path.join(data_dir, '*.csv')))
    # Exclude *.Wfm.csv headers — those are waveform data, not metadata
    header_files = [p for p in header_files if not p.endswith('.Wfm.csv')]

    pairs = []
    for hdr in header_files:
        stem = hdr[:-4]        # strip .csv
        wfm  = stem + '.Wfm.csv'
        if not os.path.exists(wfm):
            continue
        basename = os.path.basename(stem)
        # Try to parse channel & timestamp from the filename pattern
        # e.g. RefCurve_2026-04-23_1_200237
        m = re.search(r'_(\d+)_(\d+)$', basename)
        ch  = int(m.group(1)) if m else -1
        ts  = m.group(2)      if m else ''
        pairs.append({
            'stem':      basename,
            'header':    hdr,
            'wfm':       wfm,
            'channel':   ch,
            'timestamp': ts,
        })
    return pairs


# ---------------------------------------------------------------------------
# Optional impulse response loader (uses payload_signal infrastructure)
# ---------------------------------------------------------------------------

def load_impulse_response(impulse_path: str, prep=True):
    """
    Load an impulse response file, returning (time_ns, voltage).

    prep=True  (default): load via payload_signal.loadImpulse / prepImpulse,
                          which resamples to 0.25 ns (4 GHz ADC grid).
    prep=False           : load the raw file directly — preserves the
                          original time step, giving full resolution for
                          comparison plots.
    """
    if not os.path.exists(impulse_path):
        print(f"[warn] Impulse file not found: {impulse_path}")
        return None
    if prep:
        try:
            import payload_signal as payload
            imp = payload.loadImpulse(impulse_path)
            imp = payload.prepImpulse(imp)
            return imp.time, imp.voltage
        except Exception:
            pass
    # Raw load (also fallback when prep=True fails)
    try:
        dat = np.loadtxt(impulse_path)
        if dat.ndim == 2 and dat.shape[1] >= 2:
            return dat[:, 0], dat[:, 1]
        return np.arange(len(dat)), dat
    except Exception as e:
        print(f"[warn] Could not load impulse file {impulse_path}: {e}")
        return None


# ---------------------------------------------------------------------------
# Main plotting routine
# ---------------------------------------------------------------------------

def plot_scope_overview(pairs, impulse=None, save_filename=None):
    """
    Plot every discovered waveform in a grid, one panel per file.
    Optionally overlay the impulse response (rescaled to amplitude-match)
    on panels whose time window is compatible (< 20 ns span).
    """
    n = len(pairs)
    if n == 0:
        print("No waveform pairs found.")
        return

    ncols = min(n, 2)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(7 * ncols, 4 * nrows),
                             squeeze=False)

    for idx, pair in enumerate(pairs):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        meta    = parse_metadata(pair['header'])
        t, v    = load_waveform(pair['wfm'], meta)
        span_ns = meta['XStop_ns'] - meta['XStart_ns']

        # ----- main waveform -----
        ax.plot(t, v * 1e3, color='#1f77b4', linewidth=1.2,
                label='Scope waveform')

        # ----- optional impulse overlay -----
        if impulse is not None:
            t_imp, v_imp = impulse
            # Only overlay if the impulse fits within this window
            if t_imp[-1] - t_imp[0] <= span_ns * 1.5:
                # Shift impulse so its peak aligns with the waveform peak,
                # and rescale amplitude to match the waveform peak.
                peak_idx_imp = np.argmax(np.abs(v_imp))
                peak_idx_v   = np.argmax(np.abs(v))
                t_shift      = t[peak_idx_v] - t_imp[peak_idx_imp]
                scale        = (np.max(np.abs(v)) /
                                (np.max(np.abs(v_imp)) + 1e-30))
                ax.plot(t_imp + t_shift, v_imp * scale * 1e3,
                        color='#d62728', linewidth=1.0, linestyle='--',
                        alpha=0.75, label='Impulse response (rescaled)')

        # ----- labels & annotations -----
        trace_type = meta.get('TraceType', '?')
        res_ns     = meta['Resolution_ns']
        n_samp     = meta['RecordLength']

        title = (f"{pair['stem']}\n"
                 f"Ch {pair['channel']}  |  {trace_type}  |  "
                 f"res={res_ns*1e3:.0f} ps  |  "
                 f"span={span_ns:.1f} ns  |  N={n_samp}")
        ax.set_title(title, fontsize=8.5, pad=4)
        ax.set_xlabel('Time (ns)', fontsize=9)
        ax.set_ylabel('Voltage (mV)', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)
        if impulse is not None:
            ax.legend(fontsize=7, loc='upper right')

    # Hide any unused axes
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle('Oscilloscope Reference Curves — CoRaLs Antenna Measurements',
                 fontsize=13, y=1.01)
    plt.tight_layout()

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_scope_overview.png'
        plt.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved: {out}")

    plt.show()


# ---------------------------------------------------------------------------
# Antenna response vs impulse comparison with residuals
# ---------------------------------------------------------------------------

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# IMPULSE TIME OFFSET  (nanoseconds)
# Shift the impulse response left/right on the time axis to align it with
# the antenna-response waveforms.  Positive = shift right (later in time).
# Change this constant or pass --impulse-offset <value> on the command line.
IMPULSE_TIME_OFFSET_NS = 0.0

# PLOT X-AXIS WINDOW  (nanoseconds)
# Set the start and stop of the time axis shown in the ant-resp plots.
# None = use the full scope waveform range.
# Override per-run with --xstart / --xstop on the command line.
PLOT_XSTART_NS = None   # e.g. -10.0
PLOT_XSTOP_NS  = None   # e.g.  80.0
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def _interp_to_grid(t_src, v_src, t_grid):
    """Linearly interpolate v_src onto t_grid, with zeros outside the range."""
    return np.interp(t_grid, t_src, v_src, left=0.0, right=0.0)


def plot_ant_response_vs_impulse(data_dir, impulse_path,
                                  impulse_time_offset_ns=IMPULSE_TIME_OFFSET_NS,
                                  xstart_ns=PLOT_XSTART_NS,
                                  xstop_ns=PLOT_XSTOP_NS,
                                  prep_impulse=True,
                                  save_filename=None):
    """
    For each ant_resp*_waveform.csv in *data_dir*:
      - top panel   : scope waveform + amplitude-normalised impulse overlay
      - bottom panel: residual (scope − impulse_rescaled) in mV

    The impulse is:
      1. Loaded via payload_signal.loadImpulse / prepImpulse (falls back to
         plain loadtxt).
      2. Shifted in time by *impulse_time_offset_ns*.
      3. Amplitude-rescaled to match the scope waveform peak before
         subtraction.

    To adjust the time alignment see IMPULSE_TIME_OFFSET_NS above or pass
    --impulse-offset on the command line.
    """
    # ---- find ant_resp files -----------------------------------------------
    wfm_files = sorted(glob.glob(
        os.path.join(data_dir, 'ant_resp*_waveform.csv')))
    if not wfm_files:
        print(f"[warn] No ant_resp*_waveform.csv files found in {data_dir}")
        return

    # ---- load impulse -------------------------------------------------------
    imp = load_impulse_response(impulse_path, prep=prep_impulse)
    if imp is None:
        print("[error] Could not load impulse response — aborting.")
        return
    t_imp_raw, v_imp_raw = imp

    # Apply the time offset (shift the impulse along the time axis)
    t_imp = t_imp_raw + impulse_time_offset_ns

    n = len(wfm_files)

    # ---- pre-compute impulse FFT once (native resolution) ------------------
    dt_imp = (t_imp_raw[-1] - t_imp_raw[0]) / (len(t_imp_raw) - 1)   # ns
    freq_imp_ghz = np.fft.rfftfreq(len(v_imp_raw), d=dt_imp)          # GHz
    ampl_imp     = np.abs(np.fft.rfft(v_imp_raw))
    eps = 1e-30
    db_imp = 20 * np.log10(ampl_imp / (np.max(ampl_imp) + eps) + eps)

    # ---- figure layout: 3 rows × n cols  [time | residual | FFT] -----------
    fig = plt.figure(figsize=(9 * n, 14))
    outer = gridspec.GridSpec(1, n, figure=fig, hspace=0.05, wspace=0.38)

    for col, wfm_path in enumerate(wfm_files):
        stem        = os.path.basename(wfm_path).replace('_waveform.csv', '')
        hdr_path    = wfm_path.replace('_waveform.csv', '_settings.csv')
        meta        = parse_metadata(hdr_path)
        t_scope, v_scope = load_waveform(wfm_path, meta)
        span_ns   = meta['XStop_ns'] - meta['XStart_ns']
        res_ps    = meta['Resolution_ns'] * 1e3
        trace_type = meta.get('TraceType', '')
        n_samp    = meta['RecordLength']

        # inner grid: [time overlay | residual | FFT]
        inner = gridspec.GridSpecFromSubplotSpec(
            3, 1, subplot_spec=outer[col],
            height_ratios=[3, 1, 2], hspace=0.12)
        ax_top = fig.add_subplot(inner[0])
        ax_res = fig.add_subplot(inner[1], sharex=ax_top)
        ax_fft = fig.add_subplot(inner[2])

        # ---- amplitude-normalise impulse onto scope time grid --------------
        v_imp_interp = _interp_to_grid(t_imp, v_imp_raw, t_scope)
        scope_peak   = np.max(np.abs(v_scope))
        imp_peak     = np.max(np.abs(v_imp_interp))
        scale = scope_peak / (imp_peak + 1e-30)
        v_imp_scaled = v_imp_interp * scale

        residual = (v_scope - v_imp_scaled) * 1e3   # mV

        # ---- top panel: time-domain overlay --------------------------------
        # Build human-readable title from filename + scope metadata
        # e.g. "ant_resp1" → "Antenna Response 1"
        display_name = re.sub(r'ant_resp(\d+)', r'Antenna Response \1', stem)
        display_name = display_name.replace('_', ' ').strip()
        panel_title  = (f"{display_name}\n"
                        f"res={res_ps:.0f} ps  |  span={span_ns:.0f} ns  |  "
                        f"{n_samp} samples"
                        + (f"  |  {trace_type}" if trace_type else ""))

        imp_label = (f'Simulated impulse response\n'
                     f'(rescaled \u00d7{scale*1e3:.2f} mV/V, '
                     f't\u2011shift {impulse_time_offset_ns:+.2f} ns)')
        ax_top.plot(t_scope, v_scope * 1e3,
                    color='#1f77b4', linewidth=1.4, label='Measured antenna response')
        ax_top.plot(t_scope, v_imp_scaled * 1e3,
                    color='#d62728', linewidth=1.2, linestyle='--', label=imp_label)
        ax_top.axhline(0, color='gray', linewidth=0.6, alpha=0.5)
        ax_top.set_ylabel('Voltage (mV)', fontsize=10)
        ax_top.legend(fontsize=8, loc='upper right')
        ax_top.grid(True, alpha=0.3)
        ax_top.set_title(panel_title, fontsize=10, pad=4)
        plt.setp(ax_top.get_xticklabels(), visible=False)

        # ---- middle panel: residual ----------------------------------------
        ax_res.plot(t_scope, residual,
                    color='#2ca02c', linewidth=1.0, label='Scope − Impulse')
        ax_res.axhline(0, color='gray', linewidth=0.6, alpha=0.5)
        ax_res.fill_between(t_scope, residual, alpha=0.15, color='#2ca02c')
        rms = np.sqrt(np.mean(residual**2))
        ax_res.set_ylabel('Residual (mV)', fontsize=9)
        ax_res.set_xlabel('Time (ns)', fontsize=10)
        ax_res.legend(fontsize=8, loc='upper right',
                      title=f'RMS = {rms:.3f} mV', title_fontsize=7)
        ax_res.grid(True, alpha=0.3)

        # ---- apply time-axis window ----------------------------------------
        x0 = xstart_ns if xstart_ns is not None else t_scope[0]
        x1 = xstop_ns  if xstop_ns  is not None else t_scope[-1]
        ax_top.set_xlim(x0, x1)   # ax_res shares the same x axis

        # ---- bottom panel: FFT power spectra (dB, normalised) --------------
        dt_scope = (t_scope[-1] - t_scope[0]) / (len(t_scope) - 1)   # ns
        freq_scope_ghz = np.fft.rfftfreq(len(v_scope), d=dt_scope)    # GHz
        ampl_scope     = np.abs(np.fft.rfft(v_scope))
        db_scope = 20 * np.log10(ampl_scope / (np.max(ampl_scope) + eps) + eps)

        nyq_scope = 1.0 / (2.0 * dt_scope)   # GHz
        xlim_mhz  = min(nyq_scope, 3.0) * 1e3

        ax_fft.plot(freq_scope_ghz * 1e3, db_scope,
                    color='#1f77b4', linewidth=1.2, label='Measured spectrum (scope)')
        ax_fft.plot(freq_imp_ghz * 1e3, db_imp,
                    color='#d62728', linewidth=1.0, linestyle='--',
                    label='Simulated impulse spectrum')
        ax_fft.axhline(-3,  color='gray', linestyle=':',  linewidth=0.8, alpha=0.7)
        ax_fft.axhline(-10, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax_fft.text(xlim_mhz * 0.98, -3  + 1.5, '−3 dB',  fontsize=7,
                    color='gray', ha='right', va='bottom')
        ax_fft.text(xlim_mhz * 0.98, -10 + 1.5, '−10 dB', fontsize=7,
                    color='gray', ha='right', va='bottom')
        ax_fft.set_xlabel('Frequency (MHz)', fontsize=10)
        ax_fft.set_ylabel('Power spectral density (dB, norm. to peak)', fontsize=9)
        ax_fft.set_ylim(-80, 5)
        ax_fft.set_xlim(0, xlim_mhz)
        ax_fft.legend(fontsize=8, loc='upper right')
        ax_fft.grid(True, alpha=0.3)

    # Reduce top gap for single-column figures; give more room for multi-column
    suptitle_top = 0.96 if n == 1 else 0.985
    imp_source = os.path.basename(impulse_path) if impulse_path else 'unknown'
    fig.suptitle(
        f'Antenna Response vs Simulated Impulse — {imp_source}\n'
        f'Time offset: {impulse_time_offset_ns:+.2f} ns  '
        f'| Prep/resample: {"yes" if prep_impulse else "no (native resolution)"}',
        fontsize=12, y=suptitle_top)
    plt.tight_layout(rect=[0, 0, 1, suptitle_top - 0.01])

    if save_filename:
        os.makedirs('plots', exist_ok=True)
        out = f'plots/{save_filename}_ant_vs_impulse.png'
        plt.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved: {out}")

    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot R&S oscilloscope reference-curve exports.')
    parser.add_argument('--data-dir', type=str,
                        default=os.path.join(os.path.dirname(__file__), 'data'),
                        help='Directory containing .csv / .Wfm.csv pairs '
                             '(default: data/ next to this script)')
    parser.add_argument('--impulse', type=str, default=None,
                        help='Path to an impulse-response file to overlay '
                             '(e.g. impulse/corals_impulse_sci.txt)')
    parser.add_argument('--impulse-offset', type=float,
                        default=IMPULSE_TIME_OFFSET_NS,
                        help='Time offset (ns) to shift the impulse response '
                             'on the x-axis (default: IMPULSE_TIME_OFFSET_NS '
                             'constant in this file)')
    parser.add_argument('--xstart', type=float, default=PLOT_XSTART_NS,
                        help='Start of the plotted time window (ns). '
                             'Default: full scope range (PLOT_XSTART_NS in file)')
    parser.add_argument('--xstop', type=float, default=PLOT_XSTOP_NS,
                        help='End of the plotted time window (ns). '
                             'Default: full scope range (PLOT_XSTOP_NS in file)')
    parser.add_argument('--no-prep', action='store_true',
                        help='Skip payload_signal.prepImpulse and load the '
                             'impulse file at its native resolution instead '
                             'of resampling to 0.25 ns. Use this to compare '
                             'at full waveform resolution.')
    parser.add_argument('--ant-resp', action='store_true',
                        help='Plot ant_resp files vs impulse with residuals '
                             'instead of the full overview grid')
    parser.add_argument('--save', type=str, default=None,
                        help='Base filename for the saved figure')
    args = parser.parse_args()

    if args.ant_resp:
        if args.impulse is None:
            parser.error('--ant-resp requires --impulse')
        plot_ant_response_vs_impulse(
            data_dir=args.data_dir,
            impulse_path=args.impulse,
            impulse_time_offset_ns=args.impulse_offset,
            xstart_ns=args.xstart,
            xstop_ns=args.xstop,
            prep_impulse=not args.no_prep,
            save_filename=args.save,
        )
    else:
        pairs = discover_pairs(args.data_dir)
        print(f"Found {len(pairs)} waveform pair(s) in {args.data_dir}:")
        for p in pairs:
            meta = parse_metadata(p['header'])
            print(f"  {p['stem']}  ch={p['channel']}  "
                  f"span={meta['XStop_ns']-meta['XStart_ns']:.1f} ns  "
                  f"res={meta['Resolution_ns']*1e3:.0f} ps  "
                  f"type={meta.get('TraceType','?')}")
        imp = load_impulse_response(args.impulse, prep=not args.no_prep) if args.impulse else None
        plot_scope_overview(pairs, impulse=imp, save_filename=args.save)
