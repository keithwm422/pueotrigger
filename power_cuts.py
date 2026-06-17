import numpy
import matplotlib.pyplot as plt
import tools.CoRaLs_geometry as aso_geometry
import payload_signal as payload
import coherent_sum as sum
import os
import sys
import argparse



# ============================================================================
# SPRINT 6: BEAM OPTIMIZATION FOR DUAL-POL
# ============================================================================

def generate_beam_candidates(theta_min=-45, theta_max=45, phi_min=-45, phi_max=45,
                            resolution=5, radius=None):
    """
    Generate systematic grid of candidate beam directions.

    Boresight is at (phi=0, theta=0). Coordinates follow CoRaLs_geometry.py:
      phi   = azimuth from +x axis (degrees)
      theta = elevation from horizon (degrees; 0=horizon, -90=nadir)

    Parameters:
    -----------
    theta_min, theta_max : float
        Elevation range in degrees (default ±45° around boresight)
    phi_min, phi_max : float
        Azimuth range in degrees (default ±45° around boresight)
    resolution : float
        Angular spacing between candidate beam centers (degrees)
    radius : float or None
        If set, only keep candidates within this angular radius of boresight (0, 0).
        Useful for a circular FOV constraint instead of a rectangular grid.

    Returns:
    --------
    candidates : list of (phi, theta) tuples
    """
    phi_vals   = numpy.arange(phi_min,   phi_max   + 0.5*resolution, resolution)
    theta_vals = numpy.arange(theta_min, theta_max + 0.5*resolution, resolution)

    candidates = []
    for phi in phi_vals:
        for theta in theta_vals:
            if radius is not None:
                if numpy.hypot(phi, theta) > radius:
                    continue
            candidates.append((float(phi), float(theta)))

    radius_str = f", radius≤{radius}°" if radius is not None else ""
    print(f"Generated {len(candidates)} candidate beams{radius_str}:")
    print(f"  Theta: [{theta_min}°, {theta_max}°] in {resolution}° steps")
    print(f"  Phi:   [{phi_min}°, {phi_max}°] in {resolution}° steps")
    if radius is not None:
        print(f"  Circular filter: radius={radius}°")

    return candidates


def compute_beam_pattern_dualpol(phi_c, th_c, span_phi=45, span_theta=45, res_deg=5,
                                 window=160, step=40, impulse=None, 
                                 beam_patterns_h=None, beam_patterns_v=None,
                                 trigger_sectors_phi=None, psi=45.0,
                                 ringmask=(1,1,1,1), output_polarization='total',
                                 return_all_pols=False,
                                 apply_filter=True, apply_second_filter=None):
    """
    
    Uses actual waveforms and coherent_sum.py
    
    Parameters:
    -----------
    phi_c, th_c : float
        Center direction for beam steering
    span_phi, span_theta : float
        Angular span to scan around each beam (degrees)
    res_deg : float
        Angular resolution (degrees)
    window, step : int
        Power sum window parameters
    impulse : Waveform
        Impulse response
    beam_patterns_h, beam_patterns_v : tuples
        (eplane, hplane) for each polarization
    trigger_sectors_phi : list
        Antenna indices
    psi : float
        Polarization angle (degrees)
    ringmask : tuple
        Antenna mask
    output_polarization : str
        'lhcp', 'rhcp', or 'total' (default: 'total')
    return_all_pols : bool
        If True, compute and return all three polarizations (LHCP, RHCP, Total)
        regardless of output_polarization setting
    apply_filter : bool
        Apply 1.5 GHz Shannon-Whitaker FIR before beamforming (default: True)
    apply_second_filter : bool or None
        Apply second 750 MHz FIR after the first.  None means use
        payload.APPLY_SECOND_FILTER (the module-level toggle).

    Returns:
    --------
    dict with PHI, TH, heatmap, heatmap_db, delays_used
    If return_all_pols=True, also includes heatmap_lhcp, heatmap_rhcp, heatmap_total
    and corresponding _db versions
    """
    if trigger_sectors_phi is None:
        trigger_sectors_phi = [0, 1, 2, 3]
    
    ringmask = list(ringmask)
    # Define scan range — symmetric around beam center, clamped to physical limits
    phi_min = phi_c - span_phi
    phi_max = phi_c + span_phi
    th_min  = max(-90, th_c - span_theta)
    th_max  = min( 90, th_c + span_theta)
    
    phiscan = numpy.arange(phi_min, phi_max + 1e-9, res_deg)
    thscan = numpy.arange(th_min, th_max + 1e-9, res_deg)
    PHI, TH = numpy.meshgrid(phiscan, thscan)
    
    # Initialize heatmaps
    if return_all_pols:
        heatmap_lhcp = numpy.zeros_like(PHI, dtype=float)
        heatmap_rhcp = numpy.zeros_like(PHI, dtype=float)
        heatmap_total = numpy.zeros_like(PHI, dtype=float)
    else:
        heatmap = numpy.zeros_like(PHI, dtype=float)
    
    # Steering delays for beam center
    delays_center = payload.getRemappedDelays(phi_c, th_c, trigger_sectors_phi)
    delays_center_q = numpy.round(delays_center / aso_geometry.ritc_sample_step) * aso_geometry.ritc_sample_step
    
    if apply_second_filter is None:
        apply_second_filter = payload.APPLY_SECOND_FILTER

    print(f"Computing dual-pol beam pattern at (φ={phi_c}°, θ={th_c}°)")
    filt_label = ('none' if not apply_filter
                  else ('1.5G+0.75G' if apply_second_filter else '1.5G only'))
    print(f"  Filter: {filt_label}")
    print(f"  Scan range: φ=[{phi_min}°, {phi_max}°], θ=[{th_min}°, {th_max}°]")
    print(f"  Resolution: {res_deg}° ({len(thscan)}×{len(phiscan)} = {len(thscan)*len(phiscan)} points)")
    
    for i, th in enumerate(thscan):
        print(f"  Row {i+1}/{len(thscan)}: θ={th:.1f}°...", end='', flush=True)
        for j, ph in enumerate(phiscan):
            # Generate dual-pol waveforms for this sky direction (returns 8-channel array)
            waveforms, timebase, _ = payload.getPayloadWaveforms_dualpol(
                ph, th, impulse, 
                beam_patterns_h, beam_patterns_v,
                antennas=trigger_sectors_phi,
                snr=1e10,  # No noise
                psi=psi
            )
            
            # Split into H-pol (channels 0-3) and V-pol (channels 4-7)
            waveforms_h = waveforms[:4]
            waveforms_v = waveforms[4:]
            
            # Apply coherent beamforming (converts to LHCP/RHCP)
            # Use hardware-correct order: digitize → filter → 2nd filter → sum → circular
            lhcp_wf, rhcp_wf, tb = sum.coherentSum_dualpol(
                waveforms_h, waveforms_v, timebase, delays_center_q,
                downsample=False, channel_mask=ringmask, output='circular',
                apply_filter=apply_filter, digitize_first=True,
                apply_second_filter=apply_second_filter,
                fc_second=payload.FC_SECOND_FILTER
            )
            
            # Compute sliding-window power
            if return_all_pols:
                # Compute all three polarizations
                power_lhcp, _ = sum.powerSum(lhcp_wf, window=window, step=step)
                power_rhcp, _ = sum.powerSum(rhcp_wf, window=window, step=step)
                heatmap_lhcp[i, j] = numpy.max(power_lhcp)
                heatmap_rhcp[i, j] = numpy.max(power_rhcp)
                heatmap_total[i, j] = numpy.max(power_lhcp + power_rhcp)
            else:
                # Compute requested polarization only
                if output_polarization == 'lhcp':
                    power, _ = sum.powerSum(lhcp_wf, window=window, step=step)
                elif output_polarization == 'rhcp':
                    power, _ = sum.powerSum(rhcp_wf, window=window, step=step)
                else:  # 'total'
                    power_lhcp, _ = sum.powerSum(lhcp_wf, window=window, step=step)
                    power_rhcp, _ = sum.powerSum(rhcp_wf, window=window, step=step)
                    power = power_lhcp + power_rhcp
                
                heatmap[i, j] = numpy.max(power)
        
        if return_all_pols:
            print(f" done (LHCP max={numpy.max(heatmap_lhcp[i, :]):.2e}, RHCP max={numpy.max(heatmap_rhcp[i, :]):.2e})")
        else:
            print(f" done (max={numpy.max(heatmap[i, :]):.2e})")
    
    # Convert to dB
    eps = 1e-12
    if return_all_pols:
        heatmap_lhcp_db = 10.0 * numpy.log10(numpy.maximum(heatmap_lhcp, eps) / numpy.max(heatmap_total))
        heatmap_rhcp_db = 10.0 * numpy.log10(numpy.maximum(heatmap_rhcp, eps) / numpy.max(heatmap_total))
        heatmap_total_db = 10.0 * numpy.log10(numpy.maximum(heatmap_total, eps) / numpy.max(heatmap_total))
        
        return {
            'PHI': PHI, 'TH': TH,
            'heatmap_lhcp': heatmap_lhcp,
            'heatmap_rhcp': heatmap_rhcp,
            'heatmap_total': heatmap_total,
            'heatmap_lhcp_db': heatmap_lhcp_db,
            'heatmap_rhcp_db': heatmap_rhcp_db,
            'heatmap_total_db': heatmap_total_db,
            'delays_used': delays_center_q,
            'beam_phi': phi_c,
            'beam_theta': th_c,
            'psi': psi,
            'apply_filter': apply_filter,
            'apply_second_filter': apply_second_filter,
        }
    else:
        heatmap_db = 10.0 * numpy.log10(numpy.maximum(heatmap, eps) / numpy.max(heatmap))

        return {
            'PHI': PHI, 'TH': TH,
            'heatmap': heatmap,
            'heatmap_db': heatmap_db,
            'delays_used': delays_center_q,
            'beam_phi': phi_c,
            'beam_theta': th_c,
            'psi': psi,
            'output_pol': output_polarization,
            'apply_filter': apply_filter,
            'apply_second_filter': apply_second_filter,
        }


def compute_coverage_mask(beam_result, threshold_db=-3.0):
    """
    Determine which sky directions are covered by this beam.
    
    Parameters:
    -----------
    beam_result : dict
        Output from compute_beam_pattern_dualpol()
    threshold_db : float
        Coverage threshold in dB (default: -3 dB)
    
    Returns:
    --------
    coverage_mask : boolean array
        True where beam exceeds threshold
    """
    coverage_mask = beam_result['heatmap_db'] >= threshold_db
    coverage_fraction = numpy.sum(coverage_mask) / coverage_mask.size
    
    print(f"  Beam at (φ={beam_result['beam_phi']}°, θ={beam_result['beam_theta']}°): "
          f"{coverage_fraction*100:.1f}% coverage at {threshold_db} dB")
    
    return coverage_mask


def compute_all_coverage(beam_candidates, span_phi=45, span_theta=45, resolution=5,
                        window=160, step=40, threshold_db=-3.0, psi=45.0,
                        phi_range=(-90, 90), theta_range=(-90, 90),
                        save_progress=True):
    """
    Compute coverage masks for all candidate beams using full dual-pol processing.
    
    This is computationally expensive - uses actual waveforms and beamforming.
    All beams are interpolated onto a common sky grid for consistent coverage analysis.
    
    Parameters:
    -----------
    beam_candidates : list of (phi, theta)
        Candidate beam directions
    span_phi, span_theta : float
        How far to scan around each beam (degrees)
    resolution : float
        Angular resolution (degrees)
    window, step : int
        Power sum parameters
    threshold_db : float
        Coverage threshold (dB)
    psi : float
        Polarization angle (degrees)
    phi_range, theta_range : tuple
        (min, max) for common sky grid
    save_progress : bool
        Save after each beam (default: True)
    
    Returns:
    --------
    beam_results : list of dicts
        Beam pattern data for each candidate
    coverage_masks : list of boolean arrays
        Coverage mask for each candidate (on common grid)
    sky_grid : tuple (PHI, TH)
        Common sky grid for all beams
    """
    from scipy.interpolate import RegularGridInterpolator
    
    print("="*70)
    print("COMPUTING BEAM COVERAGE (DUAL-POL TIME-DOMAIN)")
    print("="*70)
    print(f"Number of beams: {len(beam_candidates)}")
    print(f"Scan span: ±{span_phi}° phi, ±{span_theta}° theta")
    print(f"Resolution: {resolution}°")
    print(f"Threshold: {threshold_db} dB")
    print(f"Polarization: ψ={psi}°")
    print(f"Common sky grid: φ={phi_range}, θ={theta_range}")
    print("="*70)
    
    # Create common sky grid
    phi_common = numpy.arange(phi_range[0], phi_range[1] + 1e-9, resolution)
    theta_common = numpy.arange(theta_range[0], theta_range[1] + 1e-9, resolution)
    PHI_common, TH_common = numpy.meshgrid(phi_common, theta_common)
    sky_grid = (PHI_common, TH_common)
    
    # Load beam patterns and impulse once
    eplane_h = payload.beamPattern(plot=False, which_plane='E', which_pol='H')
    hplane_h = payload.beamPattern(plot=False, which_plane='H', which_pol='H')
    eplane_v = payload.beamPattern(plot=False, which_plane='E', which_pol='V')
    hplane_v = payload.beamPattern(plot=False, which_plane='H', which_pol='V')
    beam_patterns_h = (eplane_h, hplane_h)
    beam_patterns_v = (eplane_v, hplane_v)
    
    impulse = payload.loadImpulse('impulse/corals_impulse_sci.txt')
    impulse = payload.prepImpulse(impulse)
    
    trigger_sectors_phi = [0, 1, 2, 3]
    ringmask = (1, 1, 1, 1)
    
    beam_results = []
    coverage_masks = []
    
    for idx, (phi_c, th_c) in enumerate(beam_candidates):
        print(f"\n[{idx+1}/{len(beam_candidates)}] Computing beam at (φ={phi_c}°, θ={th_c}°)")
        
        result = compute_beam_pattern_dualpol(
            phi_c, th_c, span_phi, span_theta, resolution,
            window, step, impulse, beam_patterns_h, beam_patterns_v,
            trigger_sectors_phi, psi=psi, ringmask=ringmask
        )
        
        # Interpolate onto common grid
        phi_local = result['PHI'][0, :]  # First row
        theta_local = result['TH'][:, 0]  # First column
        heatmap_db_local = result['heatmap_db']
        
        # Create interpolator
        interp = RegularGridInterpolator(
            (theta_local, phi_local), heatmap_db_local,
            method='linear', bounds_error=False, fill_value=-numpy.inf
        )
        
        # Interpolate to common grid
        points = numpy.stack([TH_common.ravel(), PHI_common.ravel()], axis=-1)
        heatmap_db_common = interp(points).reshape(PHI_common.shape)
        
        # Compute coverage mask on common grid
        mask = heatmap_db_common >= threshold_db
        coverage_fraction = numpy.sum(mask) / mask.size
        
        print(f"  Beam at (φ={phi_c}°, θ={th_c}°): "
              f"{coverage_fraction*100:.1f}% coverage at {threshold_db} dB")
        
        beam_results.append(result)
        coverage_masks.append(mask)
        
        # Save progress
        if save_progress:
            os.makedirs('plots', exist_ok=True)
            save_data = {
                'beam_results': beam_results,
                'coverage_masks': coverage_masks,
                'sky_grid': sky_grid,
                'candidates': beam_candidates[:idx+1],
                'threshold_db': threshold_db,
                'psi': psi
            }
            numpy.save('plots/beam_coverage_progress.npy', save_data)
    
    print("\n" + "="*70)
    print("COVERAGE COMPUTATION COMPLETE")
    print("="*70)
    
    return beam_results, coverage_masks, sky_grid


def select_optimal_beams(coverage_masks, beam_candidates, target_coverage=0.95):
    """
    Select minimum set of beams to achieve target sky coverage.
    
    Uses greedy algorithm: repeatedly select beam that covers most uncovered sky.
    
    Parameters:
    -----------
    coverage_masks : list of boolean arrays
        Coverage mask for each candidate beam
    beam_candidates : list of (phi, theta)
        Candidate beam directions
    target_coverage : float
        Desired fraction of sky coverage (default: 0.95 = 95%)
    
    Returns:
    --------
    selected_indices : list of int
        Indices of selected beams
    selected_beams : list of (phi, theta)
        Directions of selected beams
    coverage_history : list of float
        Coverage fraction after each beam added
    """
    print("\n" + "="*70)
    print("BEAM SELECTION OPTIMIZATION")
    print("="*70)
    print(f"Target coverage: {target_coverage*100:.1f}%")
    print(f"Total candidates: {len(beam_candidates)}")
    
    n_beams = len(coverage_masks)
    total_pixels = coverage_masks[0].size
    
    # Initialize: no coverage
    total_coverage = numpy.zeros_like(coverage_masks[0], dtype=bool)
    remaining_beams = set(range(n_beams))
    selected_indices = []
    coverage_history = []
    
    while True:
        current_coverage_fraction = numpy.sum(total_coverage) / total_pixels
        coverage_history.append(current_coverage_fraction)
        
        print(f"\nIteration {len(selected_indices)+1}:")
        print(f"  Current coverage: {current_coverage_fraction*100:.2f}%")
        
        if current_coverage_fraction >= target_coverage:
            print(f"  ✓ Target coverage achieved!")
            break
        
        if len(remaining_beams) == 0:
            print(f"  ⚠ No more beams available (coverage: {current_coverage_fraction*100:.2f}%)")
            break
        
        # Find beam that covers most new sky
        best_beam = None
        best_new_coverage = 0
        
        for beam_idx in remaining_beams:
            new_coverage = numpy.sum(coverage_masks[beam_idx] & ~total_coverage)
            if new_coverage > best_new_coverage:
                best_new_coverage = new_coverage
                best_beam = beam_idx
        
        # Check if any beam adds coverage
        if best_beam is None:
            print(f"  ⚠ No remaining beams add coverage (stuck at {current_coverage_fraction*100:.2f}%)")
            break
        
        # Add best beam
        selected_indices.append(best_beam)
        total_coverage |= coverage_masks[best_beam]
        remaining_beams.remove(best_beam)
        
        phi, theta = beam_candidates[best_beam]
        new_fraction = best_new_coverage / total_pixels
        print(f"  Selected beam {best_beam}: (φ={phi}°, θ={theta}°)")
        print(f"  New coverage: {new_fraction*100:.2f}% ({best_new_coverage} pixels)")
        print(f"  Total beams: {len(selected_indices)}")
    
    selected_beams = [beam_candidates[i] for i in selected_indices]
    
    print("\n" + "="*70)
    print("OPTIMIZATION COMPLETE")
    print("="*70)
    print(f"Selected {len(selected_beams)} beams for {coverage_history[-1]*100:.2f}% coverage")
    print("\nSelected beam directions:")
    for i, (phi, theta) in enumerate(selected_beams):
        print(f"  {i+1}. φ={phi:6.1f}°, θ={theta:6.1f}°")
    
    return selected_indices, selected_beams, coverage_history


def plot_coverage_analysis(beam_results, coverage_masks, selected_indices, 
                           sky_grid, coverage_history, threshold_db=-3.0,
                           filename='beam_coverage_optimization.png'):
    """
    Visualize beam coverage optimization results.
    
    Creates multi-panel plot showing:
    - Individual beam patterns
    - Combined coverage map
    - Coverage vs number of beams
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    
    PHI, TH = sky_grid
    n_selected = len(selected_indices)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 1, hspace=0.35)

    # --- Build total coverage ---
    total_coverage = numpy.zeros_like(coverage_masks[0], dtype=int)
    for idx in selected_indices:
        total_coverage += coverage_masks[idx].astype(int)

    # Panel 1: Binary coverage map (any beam ≥ threshold)
    ax_binary = fig.add_subplot(gs[0])
    binary = total_coverage > 0
    ax_binary.pcolormesh(PHI, TH, binary.astype(float), cmap='Blues', vmin=0, vmax=1,
                         shading='auto')
    # Overlay beam centers
    for idx in selected_indices:
        phi_b = beam_results[idx]['beam_phi']
        th_b  = beam_results[idx]['beam_theta']
        ax_binary.plot(phi_b, th_b, 'r*', markersize=8, markeredgecolor='white', markeredgewidth=0.8)
    ax_binary.set_xlabel('φ (azimuth) [deg]', fontsize=12)
    ax_binary.set_ylabel('θ (elevation) [deg]', fontsize=12)
    covered_frac = numpy.sum(binary) / binary.size
    ax_binary.set_title(f'Binary Coverage: {n_selected} Beams — {covered_frac*100:.1f}% of sky grid ≥ {threshold_db} dB',
                        fontsize=13, fontweight='bold')
    ax_binary.grid(True, alpha=0.3)

    # Panel 2: Overlap count map
    ax_combined = fig.add_subplot(gs[1])
    masked_cov = numpy.ma.masked_equal(total_coverage, 0)
    im = ax_combined.pcolormesh(PHI, TH, masked_cov, cmap='viridis',
                                vmin=1, vmax=numpy.max(total_coverage),
                                shading='auto')
    plt.colorbar(im, ax=ax_combined, label='Number of beams covering')
    for idx in selected_indices:
        phi_b = beam_results[idx]['beam_phi']
        th_b  = beam_results[idx]['beam_theta']
        ax_combined.plot(phi_b, th_b, 'r*', markersize=8, markeredgecolor='white', markeredgewidth=0.8)
    ax_combined.set_xlabel('φ (azimuth) [deg]', fontsize=12)
    ax_combined.set_ylabel('θ (elevation) [deg]', fontsize=12)
    ax_combined.set_title(f'Overlap Count (uncovered = white) | {threshold_db} dB threshold',
                          fontsize=13, fontweight='bold')
    ax_combined.grid(True, alpha=0.3)

    # Panel 3: Coverage vs number of beams
    ax_history = fig.add_subplot(gs[2])
    ax_history.plot(range(1, len(coverage_history)+1), 
                   numpy.array(coverage_history) * 100, 'o-', linewidth=2)
    ax_history.axhline(95, color='red', linestyle='--', label='95% target')
    ax_history.set_xlabel('Number of Beams', fontsize=12)
    ax_history.set_ylabel('Sky Coverage [%]', fontsize=12)
    ax_history.set_title('Coverage vs Beams', fontsize=13, fontweight='bold')
    ax_history.grid(True, alpha=0.3)
    ax_history.legend()
    ax_history.set_ylim(0, 105)
    
    plt.suptitle('Dual-Pol Beam Coverage Optimization (Time-Domain)', 
                fontsize=15, fontweight='bold')
    
    os.makedirs('plots', exist_ok=True)
    filepath = os.path.join('plots', filename)
    plt.savefig(filepath, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {filepath}")
    plt.close()
    
    return fig


def plot_beam_pattern_dualpol(result, output_path='plots/beam_pattern_dualpol.png',
                              contour_db=-3.0):
    """
    Plot traditional beam pattern heatmaps for dual-pol (LHCP/RHCP/Total).
    
    Creates a 3-panel figure showing oval-shaped beam profiles like the old
    coherent_sum.py single-pol visualizations.
    
    Parameters:
    -----------
    result : dict
        Output from compute_beam_pattern_dualpol() with return_all_pols=True
        Must contain: PHI, TH, heatmap_lhcp_db, heatmap_rhcp_db, heatmap_total_db
    output_path : str
        Where to save the figure
    contour_db : float
        Contour level in dB (default: -3.0 for half-power beamwidth)
    """
    import matplotlib.pyplot as plt
    
    PHI = result['PHI']
    TH = result['TH']
    lhcp_db = result['heatmap_lhcp_db']
    rhcp_db = result['heatmap_rhcp_db']
    total_db = result['heatmap_total_db']
    
    phi_c = result['beam_phi']
    th_c = result['beam_theta']
    psi = result['psi']
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Common colormap settings
    vmin, vmax = -30, 0
    cmap = 'viridis'
    
    # Panel 1: LHCP
    ax = axes[0]
    pcm = ax.pcolormesh(PHI, TH, lhcp_db, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    cs = ax.contour(PHI, TH, lhcp_db, levels=[contour_db], colors='red', linewidths=2)
    ax.plot(phi_c, th_c, 'r*', markersize=20, markeredgecolor='white', markeredgewidth=2)
    ax.set_xlabel('Azimuth φ (deg)', fontsize=12)
    ax.set_ylabel('Elevation θ (deg)', fontsize=12)
    ax.set_title(f'LHCP Beam Pattern\nψ={psi}°', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.colorbar(pcm, ax=ax, label='Power (dB)')
    
    # Panel 2: RHCP
    ax = axes[1]
    pcm = ax.pcolormesh(PHI, TH, rhcp_db, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    cs = ax.contour(PHI, TH, rhcp_db, levels=[contour_db], colors='red', linewidths=2)
    ax.plot(phi_c, th_c, 'r*', markersize=20, markeredgecolor='white', markeredgewidth=2)
    ax.set_xlabel('Azimuth φ (deg)', fontsize=12)
    ax.set_ylabel('Elevation θ (deg)', fontsize=12)
    ax.set_title(f'RHCP Beam Pattern\nψ={psi}°', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.colorbar(pcm, ax=ax, label='Power (dB)')
    
    # Panel 3: Total
    ax = axes[2]
    pcm = ax.pcolormesh(PHI, TH, total_db, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    cs = ax.contour(PHI, TH, total_db, levels=[contour_db], colors='red', linewidths=2)
    ax.plot(phi_c, th_c, 'r*', markersize=20, markeredgecolor='white', markeredgewidth=2)
    ax.set_xlabel('Azimuth φ (deg)', fontsize=12)
    ax.set_ylabel('Elevation θ (deg)', fontsize=12)
    ax.set_title(f'Total (LHCP+RHCP) Beam Pattern\nψ={psi}°', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.colorbar(pcm, ax=ax, label='Power (dB)')
    
    # Overall title
    fig.suptitle(f'Dual-Pol Beam Pattern: φ={phi_c}°, θ={th_c}° | Contour: {contour_db} dB',
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nBeam pattern visualization saved to {output_path}")


if __name__=='__main__':

    # ---------------- Argument Parsing ----------------
    parser = argparse.ArgumentParser(description="Dual-pol beam pattern simulation and optimization")
    # --- Single-beam visualization (dual-pol by default) ---
    parser.add_argument('--single', action='store_true',
                        help='Single-beam mode: visualize beam pattern with LHCP/RHCP/Total heatmaps')
    parser.add_argument('--phi', type=float, default=0.0, help='Beam center phi (deg)')
    parser.add_argument('--theta', type=float, default=0.0, help='Beam center theta (deg)')
    parser.add_argument('--span-phi', type=float, default=15.0, help='+/- scan span in phi (deg) [--single and --optimize-beams]')
    parser.add_argument('--span-theta', type=float, default=15.0, help='+/- scan span in theta (deg) [--single and --optimize-beams]')
    parser.add_argument('--resolution', type=float, default=1.0, help='Angular resolution for beam pattern scan (deg)')
    parser.add_argument('--window', type=int, default=160, help='PowerSum window (samples)')
    parser.add_argument('--pstep', type=int, default=40, help='PowerSum step (samples)')
    # --- Dual-pol beam optimization ---
    parser.add_argument('--optimize-beams', action='store_true',
                        help='Run dual-pol greedy beam coverage optimization')
    # Candidate beam grid — centered on boresight (phi=0, theta=0)
    parser.add_argument('--theta-min', type=float, default=-45,
                        help='Minimum elevation for candidates (deg, default: -45)')
    parser.add_argument('--theta-max', type=float, default=45,
                        help='Maximum elevation for candidates (deg, default: 45)')
    parser.add_argument('--phi-min', type=float, default=-45,
                        help='Minimum azimuth for candidates (deg, default: -45)')
    parser.add_argument('--phi-max', type=float, default=45,
                        help='Maximum azimuth for candidates (deg, default: 45)')
    parser.add_argument('--radius', type=float, default=None,
                        help='Circular FOV radius from boresight (0,0) to filter candidates (deg)')
    parser.add_argument('--beam-resolution', type=float, default=5,
                        help='Angular spacing between candidate beam centers (deg, default: 5)')
    parser.add_argument('--threshold-db', type=float, default=-3.0,
                        help='Coverage threshold in dB (default: -3.0)')
    parser.add_argument('--target-coverage', type=float, default=0.95,
                        help='Target sky coverage fraction (default: 0.95 = 95%%)')
    parser.add_argument('--psi', type=float, default=45.0,
                        help='Polarization angle in degrees (0=H-pol, 45=balanced, 90=V-pol, default: 45)')
    parser.add_argument('--output-pol', choices=['lhcp', 'rhcp', 'total'], default='total',
                        help='Output polarization for optimization (default: total)')

    args = parser.parse_args()
    
    # ================ SINGLE-BEAM VISUALIZATION MODE ================
    if args.single:
        phi_c = float(args.phi); th_c = float(args.theta)
        dphi = float(args.span_phi); dth = float(args.span_theta)
        angular_res = float(args.resolution)
        window = int(args.window); step = int(args.pstep)

        # Dual-pol mode: LHCP/RHCP visualization
        print("\n" + "="*70)
        print("SINGLE-BEAM DUAL-POL VISUALIZATION")
        print("="*70)
        print(f"Beam center: φ={phi_c}°, θ={th_c}°")
        print(f"Scan span: φ=±{dphi}°, θ=±{dth}°")
        print(f"Resolution: {angular_res}°")
        print(f"Polarization angle: ψ={args.psi}°")
        print(f"Window/Step: {window}/{step} samples")
        print(f"Threshold: {args.threshold_db} dB")
        print("="*70 + "\n")
        
        # Load impulse and beam patterns
        print("Loading impulse and beam patterns...")
        impulse = payload.loadImpulse('impulse/corals_impulse_sci.txt')
        impulse = payload.prepImpulse(impulse, highpass_cutoff=0.15, lowpass_cutoff=2)
        
        eplane_h = payload.beamPattern(plot=False, which_plane='E', which_pol='H')
        hplane_h = payload.beamPattern(plot=False, which_plane='H', which_pol='H')
        eplane_v = payload.beamPattern(plot=False, which_plane='E', which_pol='V')
        hplane_v = payload.beamPattern(plot=False, which_plane='H', which_pol='V')
        beam_patterns_h = (eplane_h, hplane_h)
        beam_patterns_v = (eplane_v, hplane_v)
        print("  Loaded successfully.")
        
        # Common kwargs shared by both scan calls
        _scan_kw = dict(
            phi_c=phi_c, th_c=th_c,
            span_phi=dphi, span_theta=dth, res_deg=angular_res,
            window=window, step=step,
            impulse=impulse,
            beam_patterns_h=beam_patterns_h,
            beam_patterns_v=beam_patterns_v,
            trigger_sectors_phi=[0, 1, 2, 3],
            psi=args.psi,
            ringmask=(1, 1, 1, 1),
            return_all_pols=True,
        )

        # --- Unfiltered scan ---
        print(f"\nComputing UNFILTERED dual-pol beam pattern...")
        result_unfilt = compute_beam_pattern_dualpol(apply_filter=False, **_scan_kw)

        # --- Filtered scan ---
        print(f"\nComputing FILTERED dual-pol beam pattern...")
        result_filt = compute_beam_pattern_dualpol(
            apply_filter=True,
            apply_second_filter=payload.APPLY_SECOND_FILTER,
            **_scan_kw
        )

        os.makedirs('plots', exist_ok=True)
        base = f'phi{phi_c:.0f}_theta{th_c:.0f}_psi{args.psi:.0f}'

        # --- Individual heatmap figures (one per filter state) ---
        plot_beam_pattern_dualpol(
            result_unfilt,
            output_path=os.path.join('plots', f'beam_pattern_dualpol_unfilt_{base}.png'),
            contour_db=args.threshold_db
        )
        plot_beam_pattern_dualpol(
            result_filt,
            output_path=os.path.join('plots', f'beam_pattern_dualpol_filt_{base}.png'),
            contour_db=args.threshold_db
        )

        # --- 2×3 comparison heatmap figure ---
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        filt_label = ('1.5G+0.75G' if payload.APPLY_SECOND_FILTER else '1.5G only')
        row_labels = ['Unfiltered', f'Filtered ({filt_label})']
        pol_keys = ['heatmap_lhcp_db', 'heatmap_rhcp_db', 'heatmap_total_db']
        pol_titles = ['LHCP', 'RHCP', 'Total']
        vmin_cmp, vmax_cmp = -30, 0

        fig_cmp, axes_cmp = plt.subplots(2, 3, figsize=(18, 10),
                                          sharex=True, sharey=True)
        for row_idx, (res, rlabel) in enumerate([(result_unfilt, row_labels[0]),
                                                  (result_filt,   row_labels[1])]):
            for col_idx, (pkey, ptitle) in enumerate(zip(pol_keys, pol_titles)):
                ax = axes_cmp[row_idx, col_idx]
                pcm = ax.pcolormesh(res['PHI'], res['TH'], res[pkey],
                                    cmap='viridis', vmin=vmin_cmp, vmax=vmax_cmp,
                                    shading='auto')
                ax.contour(res['PHI'], res['TH'], res[pkey],
                           levels=[args.threshold_db], colors='red', linewidths=1.5)
                ax.plot(phi_c, th_c, 'r*', markersize=14,
                        markeredgecolor='white', markeredgewidth=1.5)
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.3)
                if row_idx == 1:
                    ax.set_xlabel('Azimuth φ (deg)', fontsize=11)
                if col_idx == 0:
                    ax.set_ylabel(f'{rlabel}\nElevation θ (deg)', fontsize=11)
                title = ptitle if row_idx > 0 else ptitle
                ax.set_title(f'{ptitle}', fontsize=12, fontweight='bold')
                plt.colorbar(pcm, ax=ax, label='Power (dB)')

        fig_cmp.suptitle(
            f'Filter Comparison — φ={phi_c}°, θ={th_c}°, ψ={args.psi}°'
            f' | Threshold: {args.threshold_db} dB',
            fontsize=14, fontweight='bold'
        )
        plt.tight_layout()
        cmp_path = os.path.join('plots', f'beam_comparison_{base}.png')
        fig_cmp.savefig(cmp_path, dpi=150, bbox_inches='tight')
        plt.close(fig_cmp)
        print(f"  Comparison heatmap saved → {cmp_path}")

        # --- 1D cross-section beamwidth figure ---
        PHI  = result_filt['PHI']
        TH   = result_filt['TH']
        phi_axis  = PHI[0, :]   # 1-D phi values
        th_axis   = TH[:, 0]    # 1-D theta values

        # Nearest row/column to beam centre
        i_th  = int(numpy.argmin(numpy.abs(th_axis  - th_c)))
        j_phi = int(numpy.argmin(numpy.abs(phi_axis - phi_c)))

        def _halfpower_width(axis, profile_db):
            """Return -3 dB full-width.  Returns NaN if no crossing found."""
            idx_peak = numpy.argmax(profile_db)
            half = profile_db[idx_peak] - 3.0
            # left crossing
            left_idxs  = numpy.where(profile_db[:idx_peak] < half)[0]
            right_idxs = numpy.where(profile_db[idx_peak:] < half)[0]
            if len(left_idxs) == 0 or len(right_idxs) == 0:
                return float('nan')
            left_val  = numpy.interp(half,
                                     [profile_db[left_idxs[-1]], profile_db[left_idxs[-1]+1]],
                                     [axis[left_idxs[-1]],        axis[left_idxs[-1]+1]])
            right_val = numpy.interp(half,
                                     [profile_db[idx_peak + right_idxs[0] - 1],
                                      profile_db[idx_peak + right_idxs[0]]],
                                     [axis[idx_peak + right_idxs[0] - 1],
                                      axis[idx_peak + right_idxs[0]]])
            return abs(right_val - left_val)

        fig_bw, axes_bw = plt.subplots(1, 2, figsize=(14, 5))
        for res, lbl, ls, clr in [(result_unfilt, 'Unfiltered', '--', 'steelblue'),
                                   (result_filt,   filt_label,   '-',  'darkorange')]:
            total_db = res['heatmap_total_db']
            phi_slice = total_db[i_th, :]
            th_slice  = total_db[:, j_phi]

            bw_phi = _halfpower_width(phi_axis, phi_slice)
            bw_th  = _halfpower_width(th_axis,  th_slice)

            lbl_phi = f'{lbl}  (BW={bw_phi:.2f}°)' if numpy.isfinite(bw_phi) else f'{lbl}  (BW=N/A)'
            lbl_th  = f'{lbl}  (BW={bw_th:.2f}°)'  if numpy.isfinite(bw_th)  else f'{lbl}  (BW=N/A)'

            axes_bw[0].plot(phi_axis, phi_slice, ls=ls, color=clr, lw=2, label=lbl_phi)
            axes_bw[1].plot(th_axis,  th_slice,  ls=ls, color=clr, lw=2, label=lbl_th)

            print(f"  [{lbl}]  φ-beamwidth (at θ={th_axis[i_th]:.1f}°): "
                  f"{bw_phi:.2f}°  |  θ-beamwidth (at φ={phi_axis[j_phi]:.1f}°): {bw_th:.2f}°")

        for ax, xlabel, fixed_label in [
                (axes_bw[0], 'Azimuth φ (deg)',   f'θ = {th_axis[i_th]:.1f}°'),
                (axes_bw[1], 'Elevation θ (deg)', f'φ = {phi_axis[j_phi]:.1f}°')]:
            ax.axhline(args.threshold_db, color='gray', lw=1, ls=':', label=f'{args.threshold_db} dB')
            ax.set_xlabel(xlabel, fontsize=12)
            ax.set_ylabel('Total Power (dB)', fontsize=12)
            ax.set_title(f'1-D Slice  ({fixed_label})', fontsize=12)
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.3)

        fig_bw.suptitle(
            f'Beamwidth Cross-Sections — φ={phi_c}°, θ={th_c}°, ψ={args.psi}°',
            fontsize=13, fontweight='bold'
        )
        plt.tight_layout()
        bw_path = os.path.join('plots', f'beamwidth_{base}.png')
        fig_bw.savefig(bw_path, dpi=150, bbox_inches='tight')
        plt.close(fig_bw)
        print(f"  Beamwidth figure saved → {bw_path}")

        print("\n" + "="*70)
        print("VISUALIZATION COMPLETE")
        print("="*70)
        sys.exit(0)

    # ---------------- Beam Optimization Mode ----------------
    if args.optimize_beams:
        print("\n" + "="*70)
        print("DUAL-POL GREEDY BEAM OPTIMIZATION")
        print("="*70)

        # Generate candidate beams
        candidates = generate_beam_candidates(
            theta_min=args.theta_min,
            theta_max=args.theta_max,
            phi_min=args.phi_min,
            phi_max=args.phi_max,
            resolution=args.beam_resolution,
            radius=args.radius
        )
        
        # Compute coverage for all candidates
        beam_results, coverage_masks, sky_grid = compute_all_coverage(
            candidates,
            span_phi=args.span_phi,
            span_theta=args.span_theta,
            resolution=args.resolution,
            window=args.window,
            step=args.pstep,
            threshold_db=args.threshold_db,
            psi=args.psi,
            phi_range=(args.phi_min, args.phi_max),
            theta_range=(args.theta_min, args.theta_max),
            save_progress=True
        )
        
        # Select optimal beam set
        selected_indices, selected_beams, coverage_history = select_optimal_beams(
            coverage_masks,
            candidates,
            target_coverage=args.target_coverage
        )
        
        # Plot results
        plot_coverage_analysis(
            beam_results,
            coverage_masks,
            selected_indices,
            sky_grid,
            coverage_history,
            threshold_db=args.threshold_db,
            filename='beam_coverage_optimization.png'
        )
        
        # Save optimal beam set
        os.makedirs('plots', exist_ok=True)
        results_file = 'plots/optimal_beams_sprint6.npy'
        numpy.save(results_file, {
            'selected_beams': selected_beams,
            'selected_indices': selected_indices,
            'coverage_history': coverage_history,
            'candidates': candidates,
            'beam_results': beam_results,
            'coverage_masks': coverage_masks,
            'sky_grid': sky_grid,
            'threshold_db': args.threshold_db,
            'psi': args.psi,
            'output_pol': args.output_pol,
            'target_coverage': args.target_coverage
        })
        print(f"\nSaved optimal beam set to: {results_file}")
        
        sys.exit(0)
