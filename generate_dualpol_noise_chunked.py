"""
Generate dual-polarization noise files for 8 channels (4 antennas × 2 pols).

Usage:
    python generate_dualpol_noise_chunked.py --duration 0.5 --chunk-duration 0.01
    
    --duration: Total duration in seconds (default 0.5 sec = 2 GB)
    --chunk-duration: Size of each chunk in seconds (default 0.01 sec = 10 ms)
"""

import numpy as np
import argparse
import os
from pathlib import Path
import noise
import tools.CoRaLs_geometry as corals_geometry

def generate_dualpol_noise_chunked(duration_sec=0.5, output_dir='noise', vrms=1.0, 
                                   chunk_duration_sec=0.01):
    """
    Generate independent white thermal noise for all 8 dual-pol channels using chunked streaming.
    
    Noise is generated as white (flat spectrum). 
    
    Parameters
    ----------
    duration_sec : float
        Total duration of noise in seconds. Default is 0.5 second.
    output_dir : str
        Directory to save noise files. Default is 'noise/'.
    vrms : float
        RMS voltage for noise. Default is 1.0.
    chunk_duration_sec : float
        Duration of each memory chunk in seconds. Default is 0.01 (10 ms).
        Smaller chunks use less memory but take longer.

    Noise arrays are saved directly to disk to avoid memory overflow.
    """
    
    # Calculate required samples
    sample_rate_GHz = corals_geometry.ritc_sample_rate  # 4 GHz
    sample_step_ns = corals_geometry.ritc_sample_step   # 0.25 ns
    total_samples = int(duration_sec * 1e9 / sample_step_ns)
    chunk_samples = int(chunk_duration_sec * 1e9 / sample_step_ns)
    n_chunks = int(np.ceil(total_samples / chunk_samples))
    
    # Calculate chunk memory requirements
    fbins_chunk = int(2**np.ceil(np.log2(chunk_samples)))
    chunk_mem_GB = (8 * fbins_chunk * 16) / (1024**3)
    total_file_size_GB = (total_samples * 8 * 8) / (1024**3)  # 8 channels, float64 = 8 bytes
    
    print(f"Chunk FFT bins: {fbins_chunk:,}")
    print(f"Memory per chunk: ~{chunk_mem_GB:.2f} GB")
    print(f"Final file size: ~{total_file_size_GB:.2f} GB")
    print("")
    
    # Check if chunk size is reasonable
    if chunk_mem_GB > 10:
        print(f"ERROR: Chunk size ({chunk_mem_GB:.1f} GB) is too large!")
        print(f"Use smaller --chunk-duration (currently {chunk_duration_sec})")
        print(f"Recommended: --chunk-duration 0.001 (1 ms) or smaller")
        return None
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Initialize output files for streaming
    h_file = output_path / 'dualpol_noise_hpol.npy'
    v_file = output_path / 'dualpol_noise_vpol.npy'
    time_file = output_path / 'dualpol_noise_time.npy'
    
    # Create memory-mapped arrays for direct disk writing
    print("Creating memory-mapped output files...")
    noise_h_mmap = np.lib.format.open_memmap(h_file, mode='w+', dtype=np.float64, shape=(4, total_samples))
    noise_v_mmap = np.lib.format.open_memmap(v_file, mode='w+', dtype=np.float64, shape=(4, total_samples))
    
    # Initialize noise generator (once for all chunks)
    print("Initializing thermal noise generator...")
    thermal_noise = noise.ThermalNoise(v_rms=vrms,
        fbins=fbins_chunk,
        time_domain_sampling_rate=sample_step_ns
    )
    
    # Generate noise in chunks
    print("")
    print(f"Generating {n_chunks} chunks (streaming to disk)...")
    all_correlations = []
    
    for chunk_idx in range(n_chunks):
        # Calculate sample range for this chunk
        start_sample = chunk_idx * chunk_samples
        end_sample = min(start_sample + chunk_samples, total_samples)
        actual_chunk_size = end_sample - start_sample
        
        # Progress update
        progress = (chunk_idx + 1) / n_chunks * 100
        elapsed_time = (chunk_idx + 1) * chunk_duration_sec
        print(f"  Chunk {chunk_idx+1}/{n_chunks} ({progress:.1f}%): t={elapsed_time:.3f}s, samples {start_sample:,}-{end_sample:,}")
        
        # Generate 8 independent noise traces for this chunk
        passband, time_chunk, voltage = thermal_noise.makeNoiseWaveform(ntraces=8)
        noise_voltage = voltage.real[:, :actual_chunk_size]
        
        # Write directly to memory-mapped arrays (streams to disk)
        noise_h_mmap[0:4, start_sample:end_sample] = noise_voltage[0:4, :]
        noise_v_mmap[0:4, start_sample:end_sample] = noise_voltage[4:8, :]
        
        # Verify independence for first chunk only
        if chunk_idx == 0:
            print("  Verifying channel independence (first chunk):")
            for i in range(4):
                corr_hv = np.corrcoef(noise_voltage[i], noise_voltage[i+4])[0, 1]
                all_correlations.append(abs(corr_hv))
                print(f"    Ant {i}: H-V correlation = {corr_hv:.6f}")
    
    # Flush to disk
    print("")
    print("Flushing buffers to disk...")
    del noise_h_mmap
    del noise_v_mmap
    
    # Generate time array
    print("Generating time array...")
    time = np.arange(total_samples) * sample_step_ns
    np.save(time_file, time)
    
    # Calculate statistics
    max_corr = max(all_correlations) if all_correlations else 0.0
    print(f"Maximum H-V correlation: {max_corr:.6f}")
    
    if max_corr > 0.01:
        print(f" WARNING: Correlation {max_corr:.6f} exceeds threshold (0.01)")
    
    # Get file sizes
    h_size_GB = os.path.getsize(h_file) / (1024**3)
    v_size_GB = os.path.getsize(v_file) / (1024**3)
    
    print("")
    print("Files saved:")
    print(f" H-pol: {h_file} ({h_size_GB:.3f} GB)")
    print(f" V-pol: {v_file} ({v_size_GB:.3f} GB)")
    print(f" Time:  {time_file}")
    
    # Save metadata
    metadata = {
        'duration_sec': duration_sec,
        'total_samples': total_samples,
        'sample_rate_GHz': sample_rate_GHz,
        'sample_step_ns': sample_step_ns,
        'vrms': vrms,
        'n_antennas': 4,
        'shape_h': (4, total_samples),
        'shape_v': (4, total_samples),
        'max_correlation': float(max_corr),
        'chunk_duration_sec': chunk_duration_sec,
        'n_chunks': n_chunks
    }
    metadata_file = output_path / 'dualpol_noise_metadata.npy'
    np.save(metadata_file, metadata)
    print(f"  ✓ Metadata: {metadata_file}")
    
    print("")
    print("=" * 70)
    print("NOISE GENERATION COMPLETE")
    print("=" * 70)
    print(f"Total size: {(h_size_GB + v_size_GB):.3f} GB")
    print(f"Duration: {duration_sec} sec ({duration_sec*1000:.1f} ms)")
    print(f"Load with:")
    print(f"  noise_h = np.load('{h_file}')")
    print(f"  noise_v = np.load('{v_file}')")
    print(f"  time = np.load('{time_file}')")
    
    return None  # Data is on disk, not returned


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate dual-pol noise for 8 channels (chunked)')
    parser.add_argument('--duration', type=float, default=0.5,
                        help='Total duration in seconds (default: 0.5 sec = ~2 GB)')
    parser.add_argument('--chunk-duration', type=float, default=0.01,
                        help='Chunk duration in seconds (default: 0.01 = 10 ms, ~400 MB RAM)')
    parser.add_argument('--output-dir', type=str, default='noise',
                        help='Output directory (default: noise/)')
    parser.add_argument('--vrms', type=float, default=1.0,
                        help='RMS voltage (default: 1.0)')
    parser.add_argument('--fmin', type=float, default=0.1,
                        help='Minimum frequency in GHz (default: 0.1)')
    parser.add_argument('--fmax', type=float, default=2.0,
                        help='Maximum frequency in GHz (default: 2.0)')
    
    args = parser.parse_args()
    
    generate_dualpol_noise_chunked(
        duration_sec=args.duration,
        chunk_duration_sec=args.chunk_duration,
        output_dir=args.output_dir,
        vrms=args.vrms,
    )
