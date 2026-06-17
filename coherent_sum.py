import numpy
import matplotlib.pyplot as plt
import tools.CoRaLs_geometry as aso_geometry
import tools.constants as constants
import matplotlib.gridspec as gridspec
import tools.waveform as waveform
import tools.filters as filters
import payload_signal as payload
import math
import noise
import os

def coherentSum(waveforms, timebase, delays, downsample=False, channel_mask=None):
    '''
    Coherently sum waveforms with geometric delays.
    
    Args:
        waveforms: (n_antennas, n_samples) array 
        timebase: time array (ns)
        delays: (n_antennas,) delay array (ns) - 1D array
        downsample: if True, downsample to RITC sampling rate
        channel_mask: optional boolean mask for which antennas to include
    
    Returns:
        coh_sum: 1D coherent sum
        timebase: time array
    '''
    # Validate input dimensions
    if waveforms.ndim != 2:
        raise ValueError(f"Expected 2D waveforms [antennas, samples], got shape {waveforms.shape}")
    
    n_ant, n_samp = waveforms.shape
    
    # Default to using all antennas
    if channel_mask is None:
        channel_mask = numpy.ones(n_ant, dtype=bool)
    
    # Decimate factor for downsampling to RITC rate
    decimate_factor = int(aso_geometry.ritc_sample_step / (timebase[1] - timebase[0]))
    dt = timebase[1] - timebase[0]
    
    if downsample:
        # Downsample to RITC sampling rate
        output_length = int(n_samp / decimate_factor)
        coh_sum = numpy.zeros(output_length)
        
        for ant_idx in range(n_ant):
            if channel_mask[ant_idx]:
                _wave = waveforms[ant_idx][::decimate_factor]
                _delay = -int(numpy.round(delays[ant_idx] / (dt * decimate_factor)))
                coh_sum += numpy.roll(_wave[:len(coh_sum)], _delay)
        
        timebase = timebase[::decimate_factor][:output_length]
    
    else:
        # Use input sampling rate
        coh_sum = numpy.zeros(n_samp)
        
        for ant_idx in range(n_ant):
            if channel_mask[ant_idx]:
                _wave = waveforms[ant_idx]
                _delay = -int(numpy.round(delays[ant_idx] / dt))
                coh_sum += numpy.roll(_wave, _delay)
    
    return coh_sum, timebase


def coherentSum_dualpol(waveforms_h, waveforms_v, timebase, delays,
                        downsample=True, channel_mask=None, output='circular',
                        apply_filter=True, digitize_first=True,
                        apply_second_filter=False, fc_second=750e6,
                        return_intermediates=False):
    '''
    Coherent sum for dual polarization with circular basis conversion.

    Hardware signal chain:
    1. Downsample/digitize to 4 GHz ADC rate (if downsample=True)
    2. Apply Shannon-Whitaker digital lowpass filter at 1.5 GHz (if apply_filter=True)
    2b. Optional second lowpass filter at fc_second (if apply_second_filter=True)
    3. Apply geometric delays and coherent sum each polarization
    4. FFT to frequency domain for circular polarization conversion
    5. Apply π/2 phase shift to V-pol (multiply by j)
    6. Compute LHCP = (H + j*V)/√2 and RHCP = (H - j*V)/√2
    7. IFFT back to time domain

    Args:
        waveforms_h: (n_antennas, n_samples) H-pol array
        waveforms_v: (n_antennas, n_samples) V-pol array
        timebase: time array (ns)
        delays: (n_antennas,) delay array (ns) - same for both pols
        downsample: if True, decimate to RITC sampling rate (4 GHz)
        channel_mask: optional boolean mask for which antennas to include
        apply_filter: if True, apply first Shannon-Whitaker FIR filter (1.5 GHz)
        digitize_first: if True (default), decimate before filtering (hardware order)
        apply_second_filter: if True, apply second lowpass FIR filter at fc_second
        fc_second: cutoff frequency (Hz) for the second filter (default: 750 MHz)
        output: 'separate' | 'circular'
            'separate': return (coh_h, coh_v, timebase)
            'circular': return (lhcp, rhcp, timebase)
        return_intermediates: if True, append an intermediates dict to the return tuple

    Returns (without return_intermediates):
        - 'separate': (coh_h, coh_v, timebase)
        - 'circular': (lhcp, rhcp, timebase)
    Returns (with return_intermediates=True):
        same as above but with an extra intermediates dict appended;
        intermediates always contains 'coh_h_first_filtered'/'coh_v_first_filtered'
        (after filter 1) and 'coh_h_filtered'/'coh_v_filtered' (after final filter).
    '''

    # HARDWARE ORDER: Digitize → Filter → Sum → Circular

    # Step 1: Digitize/downsample to 4 GHz
    dt_in = timebase[1] - timebase[0]
    decimate_factor = int(aso_geometry.ritc_sample_step / dt_in) if downsample else 1
    fs_adc = 1.0 / (dt_in * decimate_factor) * 1e9  # Hz

    if downsample:
        waveforms_h_dig = waveforms_h[:, ::decimate_factor]
        waveforms_v_dig = waveforms_v[:, ::decimate_factor]
        timebase = timebase[::decimate_factor]
    else:
        waveforms_h_dig = waveforms_h
        waveforms_v_dig = waveforms_v

    # Step 2a: Apply first digital lowpass filter (1.5 GHz Shannon-Whitaker)
    if apply_filter:
        waveforms_h_filt1 = filters.apply_Shannon_Whitaker_filter(waveforms_h_dig)
        waveforms_v_filt1 = filters.apply_Shannon_Whitaker_filter(waveforms_v_dig)
    else:
        waveforms_h_filt1 = waveforms_h_dig
        waveforms_v_filt1 = waveforms_v_dig

    # Step 2b: Optional second digital lowpass filter (default 750 MHz)
    if apply_second_filter:
        waveforms_h_filt2 = filters.apply_Shannon_Whitaker_filter(waveforms_h_filt1, fc=fc_second)
        waveforms_v_filt2 = filters.apply_Shannon_Whitaker_filter(waveforms_v_filt1, fc=fc_second)
    else:
        waveforms_h_filt2 = waveforms_h_filt1
        waveforms_v_filt2 = waveforms_v_filt1

    # Step 3: Coherent sum each polarization using final filtered waveforms
    coh_h, tb = coherentSum(waveforms_h_filt2, timebase, delays, downsample=False, channel_mask=channel_mask)
    coh_v, _  = coherentSum(waveforms_v_filt2, timebase, delays, downsample=False, channel_mask=channel_mask)

    # Capture post-digitize and per-filter sums for diagnostics
    coh_h_digitized, _ = coherentSum(waveforms_h_dig,   timebase, delays, downsample=False, channel_mask=channel_mask)
    coh_v_digitized, _ = coherentSum(waveforms_v_dig,   timebase, delays, downsample=False, channel_mask=channel_mask)
    coh_h_first_filt,_ = coherentSum(waveforms_h_filt1, timebase, delays, downsample=False, channel_mask=channel_mask)
    coh_v_first_filt,_ = coherentSum(waveforms_v_filt1, timebase, delays, downsample=False, channel_mask=channel_mask)

    if output == 'separate':
        result = (coh_h, coh_v, tb)
        if return_intermediates:
            inter = {
                'coh_h_digitized':      coh_h_digitized,
                'coh_v_digitized':      coh_v_digitized,
                'coh_h_first_filtered': coh_h_first_filt,
                'coh_v_first_filtered': coh_v_first_filt,
                'coh_h_filtered':       coh_h,   # final (may equal first if no second filter)
                'coh_v_filtered':       coh_v,
                'apply_second_filter':  apply_second_filter,
                'fc_second':            fc_second,
                'decimate_factor':      decimate_factor,
                'fs_adc':               fs_adc,
            }
            return result + (inter,)
        return result

    elif output == 'circular':
        # Step 4-6: FFT → phase shift V → circular decomposition
        H_fft = numpy.fft.rfft(coh_h)
        V_fft = numpy.fft.rfft(coh_v)
        freqs  = numpy.fft.rfftfreq(len(coh_h), d=(tb[1] - tb[0]) * 1e-9)  # Hz

        # Apply π/2 phase shift to V-pol: multiply by j = exp(j*π/2)
        V_shifted = 1j * V_fft

        # Combine into circular basis
        LHCP_fft = (H_fft + V_shifted) / numpy.sqrt(2)
        RHCP_fft = (H_fft - V_shifted) / numpy.sqrt(2)

        # Step 7: IFFT back to time domain
        lhcp = numpy.fft.irfft(LHCP_fft, n=len(coh_h))
        rhcp = numpy.fft.irfft(RHCP_fft, n=len(coh_h))

        result = (lhcp, rhcp, tb)
        if return_intermediates:
            inter = {
                'coh_h_digitized':      coh_h_digitized,
                'coh_v_digitized':      coh_v_digitized,
                'coh_h_first_filtered': coh_h_first_filt,
                'coh_v_first_filtered': coh_v_first_filt,
                'coh_h_filtered':       coh_h,   # final (may equal first if no second filter)
                'coh_v_filtered':       coh_v,
                'H_fft':                H_fft,
                'V_fft':                V_fft,
                'V_shifted':            V_shifted,
                'LHCP_fft':             LHCP_fft,
                'RHCP_fft':             RHCP_fft,
                'freqs':                freqs,
                'apply_second_filter':  apply_second_filter,
                'fc_second':            fc_second,
                'decimate_factor':      decimate_factor,
                'fs_adc':               fs_adc,
            }
            return result + (inter,)
        return result

    else:
        raise ValueError(f"Unknown output mode: {output}. Use 'separate' or 'circular'.")

def check_coincidence(triggers_lhcp, triggers_rhcp, min_overlap_frames=2, 
                      max_gap_frames=1):
    '''
    Check for time-domain coincidence between LHCP and RHCP triggers.
    
    Hardware implementation: Both LHCP and RHCP must exceed their thresholds
    with time-domain overlap to generate a valid trigger.
    
    Args:
        triggers_lhcp: boolean array of LHCP above-threshold frames
        triggers_rhcp: boolean array of RHCP above-threshold frames
        min_overlap_frames: minimum number of overlapping frames for trigger
        max_gap_frames: maximum gap allowed within a trigger window (default: 1)
    
    Returns:
        trigger_indices: list of frame indices where coincidence occurs (trigger time)
        trigger_windows: list of (start, stop) tuples for each trigger
    '''
    from scipy.ndimage import label
    
    # Find frames where both LHCP and RHCP exceed threshold
    coincidence = triggers_lhcp & triggers_rhcp
    
    # Label connected components (handles gaps up to max_gap_frames)
    labeled, n_triggers = label(coincidence)
    
    trigger_windows = []
    trigger_indices = []
    
    for trigger_id in range(1, n_triggers + 1):
        trigger_frames = numpy.where(labeled == trigger_id)[0]
        
        # Check minimum overlap duration
        if len(trigger_frames) >= min_overlap_frames:
            start = trigger_frames[0]
            stop = trigger_frames[-1]
            
            # Check for acceptable gaps (no large interruptions)
            gaps = numpy.diff(trigger_frames)
            if numpy.all(gaps <= max_gap_frames + 1):
                trigger_windows.append((start, stop))
                trigger_indices.append(start)  # Trigger on first coincidence frame
    
    return trigger_indices, trigger_windows

def compute_coincidence_rate(rate_lhcp, rate_rhcp, window_time_ns=100):
    '''
    Estimate coincidence trigger rate from individual polarization rates.
    
    For independent Poisson processes:
        Rate_coincidence ≈ 2 × Rate_LHCP × Rate_RHCP × τ_window
    
    Args:
        rate_lhcp: LHCP trigger rate (Hz)
        rate_rhcp: RHCP trigger rate (Hz)
        window_time_ns: coincidence window duration (ns)
    
    Returns:
        coincidence_rate: expected coincidence rate (Hz)
    '''
    window_time_s = window_time_ns * 1e-9
    return 2 * rate_lhcp * rate_rhcp * window_time_s


def GimmeInfo(waveforms):
    '''
    Print waveform array information. Updated for 2D arrays [antennas, samples].
    '''
    if waveforms.ndim == 2:
        print("Waveforms shape: {} (antennas, samples)".format(waveforms.shape))
        print("Number of antennas: {}".format(waveforms.shape[0]))
        print("Samples per antenna: {}".format(waveforms.shape[1]))
    elif waveforms.ndim == 3:
        print("WARNING: 3D waveforms detected!")
        print("Shape: {} (sectors, rings, samples)".format(waveforms.shape))
    else:
        print("Waveforms shape: {}".format(waveforms.shape))

def powerSum(coh_sum, window=32, step=16):
    '''
    Sliding-window power: sum(|v|^2) over each window, stepped by "step".
    Returns (power/window, num_frames).
    '''
    coh_sum = numpy.asarray(coh_sum)
    if coh_sum.ndim != 1:
        raise ValueError(f"powerSum expects 1D array; got shape {coh_sum.shape}")
    if window <= 0 or step <= 0:
        raise ValueError(f"window and step must be positive; got window={window}, step={step}")
    if len(coh_sum) < window:
        raise ValueError(f"Input length {len(coh_sum)} smaller than window {window}")

    # number of frames
    num_frames = int((len(coh_sum) - window) // step) + 1
    if num_frames <= 0:
        raise ValueError("Computed num_frames <= 0; check window/step/length")

    # Power: handle complex (|v|^2)
    coh_power = numpy.abs(coh_sum)**2

    # Strided windows
    s0 = coh_power.strides[0]
    coh_sum_windowed = numpy.lib.stride_tricks.as_strided(
        coh_power,
        shape=(num_frames, window),
        strides=(s0*step, s0)
    )

    power = numpy.sum(coh_sum_windowed, axis=1) / window
    return power.astype(numpy.float64), num_frames


if __name__ == '__main__':
    """
    Test script using impulse data.
    """
    print("=" * 70)
    print("Beginning coherent_sum.py test")
    print("=" * 70)
    
    # =========================================================================
    # ADJUSTABLE PARAMETERS
    # =========================================================================
    
    # Pointing direction (degrees)
    phi_deg = 0.0       # Azimuth
    theta_deg = -30.0   # Elevation (nadir = -90)
    
    # SNR for signal injection
    snr = 5.0
    
    # Polarization angle (degrees)
    psi_deg = 45.0      # 0°=H-pol, 45°=equal H/V, 90°=V-pol
    
    # Antenna configuration (None = all antennas)
    antennas = None
    
    # Power sum window parameters
    window_samples = 160  # Standard RITC window
    step_samples = 40     # Standard RITC step
    
    # Coincidence trigger parameters
    min_overlap_frames = 2
    threshold_percentile = 50  # Use 50th percentile of power as threshold
    
    # Display options
    plot_results = True
    
    # =========================================================================
    # LOAD REAL IMPULSE DATA
    # =========================================================================
    print(f"\nLoading impulse data...")
    impulse = payload.loadImpulse('impulse/corals_impulse_sci.txt')
    impulse = payload.prepImpulse(impulse, highpass_cutoff=0.15, lowpass_cutoff=2.0)
    print(f"  Impulse loaded: {len(impulse.voltage)} samples")
    print(f"  Time step: {impulse.dt} ns")
    print(f"  Duration: {len(impulse.voltage) * impulse.dt:.1f} ns")
    
    # Load beam patterns
    print(f"\nLoading beam patterns...")
    # Sprint 2: Load separate H-pol and V-pol beam patterns
    eplane_h = payload.beamPattern(plot=False, which_plane='E', which_pol='H')
    hplane_h = payload.beamPattern(plot=False, which_plane='H', which_pol='H')
    eplane_v = payload.beamPattern(plot=False, which_plane='E', which_pol='V')
    hplane_v = payload.beamPattern(plot=False, which_plane='H', which_pol='V')
    print(f"  Loaded H-pol beam patterns (E-plane, H-plane)")
    print(f"  Loaded V-pol beam patterns (E-plane, H-plane)")
    
    # =========================================================================
    # GENERATE DUAL-POL WAVEFORMS (Sprint 2)
    # =========================================================================
    print(f"\n" + "=" * 70)
    print(f"Generating dual-pol waveforms: phi={phi_deg}°, theta={theta_deg}°, SNR={snr}, psi={psi_deg}°")
    print("=" * 70)
    
    waveforms, timebase, multipliers = payload.getPayloadWaveforms_dualpol(
        phi=phi_deg,
        el=theta_deg,
        impulse=impulse,
        beam_patterns_h=(eplane_h, hplane_h),
        beam_patterns_v=(eplane_v, hplane_v),
        antennas=antennas,
        snr=snr,
        noise=None,
        psi=psi_deg,
        plot=False
    )
    
    # Split into H and V polarizations
    waveforms_h = waveforms[:4]
    waveforms_v = waveforms[4:]
    
    delays = payload.getRemappedDelays(phi=phi_deg, el=theta_deg, antennas=antennas)
    
    print(f"  H-pol waveforms shape: {waveforms_h.shape}")
    print(f"  V-pol waveforms shape: {waveforms_v.shape}")
    print(f"  Timebase: {len(timebase)} samples, dt={timebase[1]-timebase[0]:.4f} ns")
    print(f"  Delays: {delays} ns")
    
    # =========================================================================
    # TEST: Single-pol coherent sum (using V-pol as reference)
    # =========================================================================
    print(f"\n" + "=" * 70)
    print(f"TEST: coherentSum() with V-pol waveforms")
    print("=" * 70)
    
    coh_sum, tb = coherentSum(waveforms_v, timebase, delays, downsample=False)
    
    print(f"✓ Coherent sum computed")
    print(f"  Output shape: {coh_sum.shape}")
    print(f"  Peak amplitude: {numpy.max(numpy.abs(coh_sum)):.3e}")
    print(f"  RMS: {numpy.sqrt(numpy.mean(coh_sum**2)):.3e}")
    
    # Compute power sum
    power, n_frames = powerSum(coh_sum, window=window_samples, step=step_samples)
    print(f"  Power frames: {n_frames}")
    print(f"  Max power: {numpy.max(power):.3e}")
    print(f"  Mean power: {numpy.mean(power):.3e}")
    
    # =========================================================================
    # TEST: Dual-pol with H/V decomposition 
    # =========================================================================
    print(f"\n" + "=" * 70)
    print(f"TEST: Dual-pol with real H/V decomposition")
    print("=" * 70)
    
    # Use H and V waveforms
    print(f"  Using H-pol and V-pol waveforms from getPayloadWaveforms_dualpol()")
    
    # Compute LHCP and RHCP using hardware-correct signal chain
    lhcp, rhcp, tb_dualpol = coherentSum_dualpol(
        waveforms_h, waveforms_v, timebase, delays, 
        downsample=True, apply_filter=True,
        output='circular'
    )
    
    print(f"✓ Circular polarization computed")
    print(f"  LHCP shape: {lhcp.shape}")
    print(f"  RHCP shape: {rhcp.shape}")
    print(f"  LHCP peak: {numpy.max(numpy.abs(lhcp)):.3e}")
    print(f"  RHCP peak: {numpy.max(numpy.abs(rhcp)):.3e}")
    
    # Power sums for each polarization
    power_lhcp, _ = powerSum(numpy.abs(lhcp), window=window_samples, step=step_samples)
    power_rhcp, _ = powerSum(numpy.abs(rhcp), window=window_samples, step=step_samples)
    
    # Set thresholds
    threshold_lhcp = numpy.percentile(power_lhcp, threshold_percentile)
    threshold_rhcp = numpy.percentile(power_rhcp, threshold_percentile)
    
    print(f"  LHCP threshold ({threshold_percentile}%ile): {threshold_lhcp:.3e}")
    print(f"  RHCP threshold ({threshold_percentile}%ile): {threshold_rhcp:.3e}")
    
    # Detect triggers
    triggers_lhcp = power_lhcp > threshold_lhcp
    triggers_rhcp = power_rhcp > threshold_rhcp
    
    print(f"  LHCP triggers: {numpy.sum(triggers_lhcp)} frames")
    print(f"  RHCP triggers: {numpy.sum(triggers_rhcp)} frames")
    
    # Check coincidence
    trigger_idx, trigger_win = check_coincidence(
        triggers_lhcp, triggers_rhcp, min_overlap_frames=min_overlap_frames
    )
    
    print(f"  Coincidence events: {len(trigger_idx)}")
    if len(trigger_idx) > 0:
        for i, (idx, (start, stop)) in enumerate(zip(trigger_idx, trigger_win)):
            print(f"    Event {i+1}: frame {idx}, window ({start}, {stop})")
    
    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    if plot_results:
        print(f"\n" + "=" * 70)
        print("Generating plots...")
        print("=" * 70)
        
        fig, axes = plt.subplots(4, 1, figsize=(14, 12))
        
        # Plot 1: Single antenna waveform vs coherent sum
        axes[0].plot(timebase, waveforms_v[0], 'gray', alpha=0.5, linewidth=0.5, label='Antenna 0 V-pol')
        axes[0].plot(tb, coh_sum, 'b-', linewidth=1.5, label='Coherent Sum')
        axes[0].set_xlabel('Time (ns)')
        axes[0].set_ylabel('Amplitude')
        axes[0].set_title(f'Real Impulse: Antenna 0 V-pol vs Coherent Sum (phi={phi_deg}°, theta={theta_deg}°)')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Plot 2: Power vs time (single-pol V-pol)
        power_time = numpy.arange(n_frames) * step_samples * (timebase[1] - timebase[0])
        axes[1].plot(power_time, power, 'k-', linewidth=2, label='V-pol Power')
        axes[1].set_xlabel('Time (ns)')
        axes[1].set_ylabel('Power')
        axes[1].set_title(f'Single-Pol Power Sum (window={window_samples}, step={step_samples})')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        axes[1].set_yscale('log')
        
        # Plot 3: LHCP and RHCP power
        # Compute time axis for dual-pol (uses digitized timebase)
        n_frames_dualpol = len(power_lhcp)
        power_time_dualpol = numpy.arange(n_frames_dualpol) * step_samples * (tb_dualpol[1] - tb_dualpol[0])
        axes[2].plot(power_time_dualpol, power_lhcp, 'r-', linewidth=2, alpha=0.7, label='LHCP Power')
        axes[2].plot(power_time_dualpol, power_rhcp, 'b-', linewidth=2, alpha=0.7, label='RHCP Power')
        axes[2].axhline(threshold_lhcp, color='r', linestyle='--', alpha=0.5, label=f'LHCP Threshold')
        axes[2].axhline(threshold_rhcp, color='b', linestyle='--', alpha=0.5, label=f'RHCP Threshold')
        axes[2].set_xlabel('Time (ns)')
        axes[2].set_ylabel('Power')
        axes[2].set_title('Dual-Pol Power (Circular Basis)')
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)
        axes[2].set_yscale('log')
        
        # Plot 4: Coincidence triggers
        axes[3].fill_between(power_time_dualpol, 0, triggers_lhcp.astype(int),
                             color='red', alpha=0.3, label='LHCP Triggers')
        axes[3].fill_between(power_time_dualpol, 0, triggers_rhcp.astype(int),
                             color='blue', alpha=0.3, label='RHCP Triggers')
        axes[3].fill_between(power_time_dualpol, 0, (triggers_lhcp & triggers_rhcp).astype(int),
                             color='purple', alpha=0.7, label='Coincidence')
        for idx in trigger_idx:
            axes[3].axvline(idx * step_samples * (tb_dualpol[1] - tb_dualpol[0]), 
                           color='green', linestyle='--', linewidth=2, alpha=0.8)
        axes[3].set_xlabel('Time (ns)')
        axes[3].set_ylabel('Trigger State')
        axes[3].set_title(f'Coincidence Triggers ({len(trigger_idx)} events, min_overlap={min_overlap_frames})')
        axes[3].set_ylim(-0.1, 1.2)
        axes[3].legend()
        axes[3].grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_file = 'plots/coherent_sum_real_impulse.png'
        plt.savefig(output_file, dpi=150)
        print(f" Plots saved to: {output_file}")
        plt.show()
