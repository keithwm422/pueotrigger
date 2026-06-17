"""
Generate power sums for dual-polarization noise with Shannon-Whitaker filtering.

- Loads pre-generated dual-pol noise files (H-pol and V-pol)
- Applies Shannon-Whitaker digital lowpass filter to each antenna
- Performs coherent beamforming in H/V basis
- Converts to circular polarization (LHCP/RHCP)
- Computes power sums for both LHCP and RHCP
- Saves power distributions for threshold curve fitting

Usage:
    python generate_power_dualpol.py --duration 0.5 --window 160 --step 40
"""

import numpy as np
import os
import argparse
from pathlib import Path
import tools.CoRaLs_geometry as corals_geometry
import tools.filters as filters
import coherent_sum

def load_dualpol_noise(noise_dir='noise', duration_sec=None, start_sample=0):
    """
    Load pre-generated dual-pol noise files.
    
    Parameters
    ----------
    noise_dir : str
        Directory containing noise files
    duration_sec : float or None
        Duration to load in seconds. If None, load all available.
    start_sample : int
        Starting sample index (for loading subsets)
    
    Returns
    -------
    noise_h : ndarray (4, n_samples)
        H-pol noise for 4 antennas
    noise_v : ndarray (4, n_samples)
        V-pol noise for 4 antennas
    sample_rate_hz : float
        Sampling rate in Hz
    """
    noise_path = Path(noise_dir)
    h_file = noise_path / 'dualpol_noise_hpol.npy'
    v_file = noise_path / 'dualpol_noise_vpol.npy'
    
    if not h_file.exists() or not v_file.exists():
        raise FileNotFoundError(f"Noise files not found in {noise_dir}/. Run generate_dualpol_noise_chunked.py first.")
    
    # Load metadata to get sample rate
    metadata_file = noise_path / 'dualpol_noise_metadata.npy'
    if metadata_file.exists():
        metadata = np.load(metadata_file, allow_pickle=True).item()
        sample_rate_hz = metadata['sample_rate_GHz'] * 1e9
        print(f"Loaded metadata: {metadata['duration_sec']} sec, {metadata['sample_rate_GHz']} GHz")
    else:
        sample_rate_hz = corals_geometry.ritc_sample_rate * 1e9
        print(f"Warning: metadata not found, using default sample rate {sample_rate_hz/1e9} GHz")
    
    # Load noise arrays (memory-mapped for efficiency)
    noise_h_full = np.load(h_file, mmap_mode='r')
    noise_v_full = np.load(v_file, mmap_mode='r')
    
    print(f"Noise files: {noise_h_full.shape} (antennas, samples)")
    
    # Determine sample range
    if duration_sec is None:
        end_sample = noise_h_full.shape[1]
    else:
        n_samples = int(duration_sec * sample_rate_hz)
        end_sample = min(start_sample + n_samples, noise_h_full.shape[1])
    
    # Load subset into RAM
    noise_h = np.array(noise_h_full[:, start_sample:end_sample])
    noise_v = np.array(noise_v_full[:, start_sample:end_sample])
    
    actual_duration = noise_h.shape[1] / sample_rate_hz
    print(f"Loaded {actual_duration:.3f} sec ({noise_h.shape[1]:,} samples) starting at sample {start_sample}")
    
    return noise_h, noise_v, sample_rate_hz


def generate_power_dualpol_chunked(noise_dir='noise', duration_sec=0.5, 
                                   chunk_duration_sec=0.05,
                                   sample_rate_hz=None, window=160, step=40,
                                   phi_deg=0, theta_deg=-30, apply_filter=True,
                                   apply_second_filter=False, fc_second=750e6):
    """
    Generate power sums for dual-pol noise with beamforming - CHUNKED PROCESSING.
    
    Processes noise in chunks to avoid memory overflow.
    
    Parameters
    ----------
    noise_dir : str
        Directory containing noise files
    duration_sec : float
        Total duration to process
    chunk_duration_sec : float
        Process in chunks of this duration (default 0.05 = 50 ms)
    sample_rate_hz : float or None
        Sampling rate (if None, read from metadata)
    window : int
        Power sum window size in samples
    step : int
        Power sum step size in samples
    phi_deg : float
        Azimuth angle for beamforming (degrees)
    theta_deg : float
        Elevation angle for beamforming (degrees)
    apply_filter : bool
        Apply Shannon-Whitaker digital lowpass filter (1.5 GHz)
    apply_second_filter : bool
        Apply second lowpass filter after Shannon-Whitaker (default: False)
    fc_second : float
        Cutoff frequency of second filter in Hz (default: 750 MHz)
    
    Returns
    -------
    power_lhcp : ndarray
        Power frames for LHCP
    power_rhcp : ndarray
        Power frames for RHCP
    n_frames : int
        Number of power frames
    """
    from payload_signal import getRemappedDelays
    from pathlib import Path
    
    # Load metadata
    noise_path = Path(noise_dir)
    metadata_file = noise_path / 'dualpol_noise_metadata.npy'
    if metadata_file.exists():
        metadata = np.load(metadata_file, allow_pickle=True).item()
        if sample_rate_hz is None:
            sample_rate_hz = metadata['sample_rate_GHz'] * 1e9
        print(f"Loaded metadata: {metadata['duration_sec']} sec, {metadata['sample_rate_GHz']} GHz")
    else:
        if sample_rate_hz is None:
            sample_rate_hz = corals_geometry.ritc_sample_rate * 1e9
        print(f"Warning: metadata not found, using default sample rate {sample_rate_hz/1e9} GHz")
    
    # Open noise files as memory-mapped (don't load into RAM yet)
    h_file = noise_path / 'dualpol_noise_hpol.npy'
    v_file = noise_path / 'dualpol_noise_vpol.npy'
    
    if not h_file.exists() or not v_file.exists():
        raise FileNotFoundError(f"Noise files not found in {noise_dir}/")
    
    noise_h_mmap = np.load(h_file, mmap_mode='r')
    noise_v_mmap = np.load(v_file, mmap_mode='r')
    
    total_samples_available = noise_h_mmap.shape[1]
    total_samples = int(duration_sec * sample_rate_hz)
    total_samples = min(total_samples, total_samples_available)
    
    chunk_samples = int(chunk_duration_sec * sample_rate_hz)
    n_chunks = int(np.ceil(total_samples / chunk_samples))
    
    print(f"\nChunked processing:")
    print(f"  Total duration: {duration_sec} sec ({total_samples:,} samples)")
    print(f"  Chunk duration: {chunk_duration_sec} sec ({chunk_samples:,} samples)")
    print(f"  Number of chunks: {n_chunks}")
    
    print(f"\nBeamforming configuration:")
    print(f"  Direction: phi={phi_deg}°, theta={theta_deg}°")
    print(f"  Filter 1 (1.5 GHz Shannon-Whitaker): {'enabled' if apply_filter else 'disabled'}")
    print(f"  Filter 2 ({fc_second/1e6:.0f} MHz lowpass): {'enabled' if apply_second_filter else 'disabled'}")
    
    # Get geometric delays for beamforming direction
    n_ant = noise_h_mmap.shape[0]
    ant_indices = list(range(n_ant))
    delays = getRemappedDelays(phi_deg, theta_deg, antennas=ant_indices)
    print(f"  Delays: {delays} ns")
    
    # Process in chunks
    power_lhcp_chunks = []
    power_rhcp_chunks = []
    tail_lhcp = None
    tail_rhcp = None
    
    for chunk_idx in range(n_chunks):
        start_sample = chunk_idx * chunk_samples
        end_sample = min(start_sample + chunk_samples, total_samples)
        actual_chunk_size = end_sample - start_sample
        
        progress = (chunk_idx + 1) / n_chunks * 100
        print(f"\nChunk {chunk_idx+1}/{n_chunks} ({progress:.1f}%): samples {start_sample:,} to {end_sample:,}")
        
        # Load chunk into RAM
        noise_h_chunk = np.array(noise_h_mmap[:, start_sample:end_sample])
        noise_v_chunk = np.array(noise_v_mmap[:, start_sample:end_sample])
        
        # Apply Shannon-Whitaker filter (1.5 GHz)
        if apply_filter:
            noise_h_chunk = filters.apply_Shannon_Whitaker_filter(noise_h_chunk, fs=sample_rate_hz)
            noise_v_chunk = filters.apply_Shannon_Whitaker_filter(noise_v_chunk, fs=sample_rate_hz)
        
        # Apply second lowpass filter (e.g. 750 MHz)
        if apply_second_filter:
            noise_h_chunk = filters.apply_Shannon_Whitaker_filter(noise_h_chunk, fs=sample_rate_hz, fc=fc_second)
            noise_v_chunk = filters.apply_Shannon_Whitaker_filter(noise_v_chunk, fs=sample_rate_hz, fc=fc_second)
        
        # Create timebase for this chunk
        dt = 1.0 / sample_rate_hz * 1e9  # ns
        timebase_chunk = np.arange(actual_chunk_size) * dt
        
        # Coherent sum with circular conversion
        lhcp_chunk, rhcp_chunk, _ = coherent_sum.coherentSum_dualpol(
            noise_h_chunk, noise_v_chunk, timebase_chunk, delays,
            downsample=False, output='circular', apply_filter=False  # Already filtered
        )
        
        # Prepend tail from previous chunk for continuity
        if tail_lhcp is not None:
            lhcp_chunk = np.concatenate([tail_lhcp, lhcp_chunk])
            rhcp_chunk = np.concatenate([tail_rhcp, rhcp_chunk])
        
        # Compute power sums for this chunk
        if len(lhcp_chunk) >= window:
            power_lhcp_chunk, _ = coherent_sum.powerSum(lhcp_chunk, window=window, step=step)
            power_rhcp_chunk, _ = coherent_sum.powerSum(rhcp_chunk, window=window, step=step)
            
            power_lhcp_chunks.append(power_lhcp_chunk)
            power_rhcp_chunks.append(power_rhcp_chunk)
            
            # Keep tail for next chunk (overlap to ensure continuity)
            keep = window + step * 2
            tail_lhcp = lhcp_chunk[-keep:] if keep < len(lhcp_chunk) else lhcp_chunk
            tail_rhcp = rhcp_chunk[-keep:] if keep < len(rhcp_chunk) else rhcp_chunk
            
            print(f"  Power frames this chunk: {len(power_lhcp_chunk):,}")
        else:
            # Not enough samples yet, accumulate
            tail_lhcp = lhcp_chunk
            tail_rhcp = rhcp_chunk
    
    # Concatenate all chunks
    print(f"\nConcatenating {len(power_lhcp_chunks)} chunks...")
    power_lhcp = np.concatenate(power_lhcp_chunks)
    power_rhcp = np.concatenate(power_rhcp_chunks)
    n_frames = len(power_lhcp)
    
    print(f"  Total power frames: {n_frames:,} ({n_frames * step / sample_rate_hz:.3f} sec equivalent)")
    
    return power_lhcp, power_rhcp, n_frames


def main():
    parser = argparse.ArgumentParser(description="Generate dual-pol power sums for threshold analysis")
    parser.add_argument('--noise-dir', default='noise', help='Directory with noise files')
    parser.add_argument('--duration', type=float, default=None, help='Duration to process (sec), default: use all available')
    parser.add_argument('--chunk-duration', type=float, default=0.05, help='Chunk size (sec)')
    parser.add_argument('--window', type=int, default=160, help='Power window size (samples)')
    parser.add_argument('--step', type=int, default=40, help='Power step size (samples)')
    parser.add_argument('--phi', type=float, default=0.0, help='Azimuth angle (deg)')
    parser.add_argument('--theta', type=float, default=-30.0, help='Elevation angle (deg)')
    parser.add_argument('--no-filter', action='store_true', help='Disable Shannon-Whitaker filter')
    parser.add_argument('--apply-second-filter', action='store_true', help='Apply second 750 MHz lowpass filter after Shannon-Whitaker')
    parser.add_argument('--fc-second', type=float, default=750e6, help='Second filter cutoff frequency in Hz (default: 750e6)')
    parser.add_argument('--suffix', type=str, default='', help='Suffix appended to output filenames, e.g. "_nofilter" to distinguish runs')
    parser.add_argument('--output-dir', default='noise', help='Output directory')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("DUAL-POL POWER GENERATION FOR THRESHOLD ANALYSIS (CHUNKED)")
    print("=" * 70)
    print(f"Threshold curves with placeholder N_beams=100")
    print("")
    
    # Get available noise duration if not specified
    if args.duration is None:
        from pathlib import Path
        metadata_file = Path(args.noise_dir) / 'dualpol_noise_metadata.npy'
        if metadata_file.exists():
            metadata = np.load(metadata_file, allow_pickle=True).item()
            args.duration = metadata['duration_sec']
            print(f"Using all available noise: {args.duration:.3f} sec")
        else:
            args.duration = 0.5  # fallback
            print(f"Warning: Metadata not found, using {args.duration} sec")
    
    # Generate power sums in chunks
    power_lhcp, power_rhcp, n_frames = generate_power_dualpol_chunked(
        noise_dir=args.noise_dir,
        duration_sec=args.duration,
        chunk_duration_sec=args.chunk_duration,
        window=args.window, step=args.step,
        phi_deg=args.phi, theta_deg=args.theta,
        apply_filter=not args.no_filter,
        apply_second_filter=args.apply_second_filter,
        fc_second=args.fc_second,
    )
    
    # Save results
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)
    
    if args.suffix:
        suffix = args.suffix
    elif args.no_filter:
        suffix = '_nofilter'
    elif args.apply_second_filter:
        suffix = f'_filt2_{int(args.fc_second/1e6)}MHz'
    else:
        suffix = ''
    lhcp_file = output_path / f'power_lhcp_{args.window}_{args.step}{suffix}.npy'
    rhcp_file = output_path / f'power_rhcp_{args.window}_{args.step}{suffix}.npy'
    
    sample_rate_hz = corals_geometry.ritc_sample_rate * 1e9
    
    np.save(lhcp_file, power_lhcp.astype(np.float32))
    np.save(rhcp_file, power_rhcp.astype(np.float32))
    
    print(f"\n{'=' * 70}")
    print("SAVED:")
    print(f"  ✓ LHCP power: {lhcp_file} ({power_lhcp.nbytes/1e6:.1f} MB)")
    print(f"  ✓ RHCP power: {rhcp_file} ({power_rhcp.nbytes/1e6:.1f} MB)")
    
    # Print statistics
    print(f"\n{'=' * 70}")
    print("STATISTICS:")
    print(f"  LHCP: mean={power_lhcp.mean():.4f}, std={power_lhcp.std():.4f}")
    print(f"  RHCP: mean={power_rhcp.mean():.4f}, std={power_rhcp.std():.4f}")
    print(f"  Frames: {n_frames:,}")
    print(f"  Time equivalent: {n_frames * args.step / sample_rate_hz:.3f} sec")
    print("=" * 70)


if __name__ == '__main__':
    main()
