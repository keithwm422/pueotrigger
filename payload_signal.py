import numpy
import tools.waveform as waveform
import tools.delays as delays
import tools.CoRaLs_geometry as aso_geometry
from scipy import interpolate
from scipy.signal import lfilter, butter, cheby1
import matplotlib.pyplot as plt
import coherent_sum

# ─────────────────────────────────────────────────────────────────────────────
# FILTER CONFIGURATION TOGGLE
#   APPLY_SECOND_FILTER : bool
#       True  → apply a second lowpass FIR at FC_SECOND_FILTER after the
#               primary 1.5 GHz Shannon-Whitaker filter.
#       False → use only the primary 1.5 GHz filter (legacy behaviour).
#   FC_SECOND_FILTER : float (Hz)
#       Cutoff frequency for the second filter.  Per Filter_plan.md: 750 MHz.
# ─────────────────────────────────────────────────────────────────────────────
APPLY_SECOND_FILTER = True
FC_SECOND_FILTER    = 750e6   # Hz

def loadImpulse(filename='impulse/corals_impulse_sci.txt'):

    dat=numpy.loadtxt(filename)
    impulse=waveform.Waveform(dat[:,1], time=dat[:,0])
    print(impulse.gimmeInfo(filename))
    return impulse

def prepImpulse(impulse, upsample=10, filter=False, highpass_cutoff=0.28, lowpass_cutoff=1.10 ):
    '''
    Prepare impulse: zeropad and normalize.
    
    NOTE: Analog filtering removed - digital Shannon-Whitaker lowpass filter
    will be applied after beamforming (matching hardware implementation).
    
    Args:
        impulse: Waveform object
        upsample: unused (kept for backward compatibility)
        filter: unused (kept for backward compatibility) 
        highpass_cutoff: unused (kept for backward compatibility)
        lowpass_cutoff: unused (kept for backward compatibility)
    '''
    impulse.zeropad(4096)
    impulse.fft()
    impulse.time = impulse.time-impulse.time[0]
    
    # Set Vpp = 1 (normalize)
    impulse.voltage = impulse.voltage / (numpy.max(impulse.voltage) - numpy.min(impulse.voltage))

    #start = numpy.argmax(impulse.voltage)-3200
    #impulse.takeWindow([start, start+10000])
    #impulse.time = impulse.time-impulse.time[0]

    
    return impulse

# beamPattern interpolation from datafiles or function. For corals right now this is function not data.
# we should break this into 2 sections for E_plane and H_plane of antennas and then we can do beam pattern for them based on the plane and both phi and el for each plane...
# we use h_plane and e_plane as the same resp, just h_plane is rotated by 90 for resp with respect to el I guess...

# This function returns a E_plane beamPattern tuple that has both a el and phi resp (2D spline)
# don't really need a spline since its all a function currently, but eventually need a spline for real resp data included...
# This function returns a H_plane beamPattern tuple that has both a el and phi resp

def beamPattern(plot=False, which_plane='E', which_pol='V', selected_az=None, selected_el=None):
    '''
    CoRaLS proxy beam pattern
    '''
    num_az = 361
    az = numpy.linspace(-180, 180, num_az)
    num_el = 181
    el = numpy.linspace(-90, 90, num_el)
    resp = 1 * numpy.ones((num_el, num_az))
    el_i = 0
    while el_i < len(el):
        az_i = 0
        while az_i < len(az):
            resp[el_i, az_i] = 10 * numpy.log10(
                resp[el_i, az_i] * (numpy.cos(numpy.radians(el[el_i])) ** (6)) * numpy.cos(numpy.radians((az[az_i]) / (2.0))) ** (10)
            )
            az_i += 1
        el_i += 1

    if which_plane == 'E':
        if which_pol == 'V':
            # V-pol: E-plane is vertical (uses elevation angle)
            plane_interp = interpolate.interp1d(el, resp[:, int(num_az / 2)], kind='cubic')
        elif which_pol == 'H':
            # H-pol: E-plane is horizontal (uses azimuth angle)
            plane_interp = interpolate.interp1d(az, resp[int(num_el / 2), :], kind='cubic')
    elif which_plane == 'E1':
        plane_interp = interpolate.interp1d(el, resp[:, int(num_az / 2)], kind='cubic')
    elif which_plane == 'H':
        if which_pol == 'V':
            # V-pol: H-plane is horizontal (uses azimuth angle)
            plane_interp = interpolate.interp1d(az, resp[int(num_el / 2), :], kind='cubic')
        elif which_pol == 'H':
            # H-pol: H-plane is vertical (uses elevation angle)
            plane_interp = interpolate.interp1d(el, resp[:, int(num_az / 2)], kind='cubic')
    elif which_plane == 'H1':
        plane_interp = interpolate.interp1d(el, resp[:, int(num_az / 2)], kind='cubic')

    if plot:
        az_i = 0
        while az_i < num_az:
            plt.plot(el, resp[:, az_i])
            if az_i == int(num_az / 2):
                plt.plot(el, resp[:, az_i], '--', label=str(which_plane))
            az_i += 10

        # Highlight the selected azimuth and elevation if provided
        if selected_az is not None:
            plt.axvline(x=selected_az, color='red', linestyle=':', label=f'Selected az={selected_az}°')
        if selected_el is not None:
            plt.axvline(x=selected_el, color='blue', linestyle=':', label=f'Selected el={selected_el}°')

        plt.legend(loc='upper left')
        plt.grid(True)
        plt.xlabel('off-boresight angle [deg.]')
        plt.ylabel('amplitude [dB]')
        plt.xlim([-50, 50])
        plt.ylim([-10, 1])
        # Add annotation box
        if selected_az is not None or selected_el is not None:
            txt = ""
            if selected_az is not None:
                txt += f"Azimuth: {selected_az}°\n"
            if selected_el is not None:
                txt += f"Elevation: {selected_el}°"
            plt.gca().text(0.05, 0.95, txt, transform=plt.gca().transAxes, fontsize=10,
                           verticalalignment='top', bbox=dict(boxstyle="round", fc="w", alpha=0.7))
        plt.show()
    return plane_interp
def dBtoVoltsAtten(db_value):
    '''
    does what it says
    '''
    atten_fraction = 10**(db_value/20)
    return atten_fraction

def getRemappedDelays(phi, el, antennas=None):
    '''
    converts delays from delays.getAllDelays(), from a dict to a numpy array
    antennas: list of antenna indices (e.g., [0,1,2,3])
    '''
    if antennas is None:
        antennas = [0, 1, 2, 3]
    
    delay = delays.getAllDelays([phi], [el], antennas=antennas)
    delays_remapped = numpy.zeros(len(antennas))

    for ant_idx, ant_name in enumerate([f'Ant{a+1}' for a in antennas]):
        delays_remapped[ant_idx] = delay[0]['delays'][ant_name]

    return delays_remapped


def getPayloadWaveforms(phi, el, impulse, beam_pattern, antennas=None, snr=1, noise=None, plot=False, downsample=False):

    if antennas is None:
        antennas = [0, 1, 2, 3]
    
    delay = delays.getAllDelays([phi], [el], antennas=antennas)
    
    # Construct trigger_waves: shape is [num_antennas, num_samples]
    trigger_waves = numpy.zeros((len(antennas), len(impulse.voltage)))
    
    print("phi = {:.2f}, theta = {:.2f}".format(round(delay[0]["phi"], 2), round(delay[0]["theta"], 2)))
    
    # multiplier stores beam pattern attenuation for each antenna
    multiplier=[]
    
    # Generate waveforms for each antenna
    for ant_idx, ant_name in enumerate([f'Ant{a+1}' for a in antennas]):
        # Get delay for this antenna
        ant_delay = delay[0]['delays'][ant_name]
        delay_samples = int(numpy.round(ant_delay / impulse.dt))
        
        # Calculate off-boresight angles relative to antenna's pointing direction
        phi_interp = phi - aso_geometry.phi_ant[antennas[ant_idx]]
        # Wrap phi_interp to [-180, 180]
        if phi_interp > 180 and phi_interp < 360:
            phi_interp -= 360
        elif phi_interp > 360:
            phi_interp -= 540
        elif phi_interp < -180 and phi_interp > -360:
            phi_interp += 360
        elif phi_interp < -360:
            phi_interp += 540
        
        el_interp = el - aso_geometry.theta_ant[antennas[ant_idx]]
        # Wrap el_interp to [-90, 90]
        if el_interp > 90 and el_interp < 180:
            el_interp -= 180
        elif el_interp > 180:
            el_interp -= 270
        elif el_interp < -90 and el_interp > -180:
            el_interp += 180
        elif el_interp < -180:
            el_interp += 270
        
        # Apply delayed impulse with beam pattern attenuation
        # beam_pattern[0] = E-plane response (elevation dependent)
        # beam_pattern[1] = H-plane response (azimuth dependent)
        # Currently both are multiplied together into a single waveform
        trigger_waves[ant_idx] = \
            numpy.roll(impulse.voltage * 2 * snr, delay_samples) * \
            dBtoVoltsAtten(beam_pattern[1](phi_interp)) * \
            dBtoVoltsAtten(beam_pattern[0](el_interp))
        
        # Store combined beam pattern attenuation
        multiplier.append(dBtoVoltsAtten(beam_pattern[1](phi_interp)) * \
            dBtoVoltsAtten(beam_pattern[0](el - aso_geometry.theta_ant[antennas[ant_idx]])))
        
        # Add noise if provided
        if noise is not None:
            print("noise shape " + str(noise.shape))
            trigger_waves[ant_idx] += noise[ant_idx]

    '''
    #numpyfied waveform generation. Factor of 2 multipler since impulse Vpp putatively normalized to = 1.0
    trigger_waves2 = numpy.roll(numpy.tile(impulse.voltage,len(trigger_sectors)*len(ring_map)).reshape((len(trigger_sectors)*len(ring_map), len(impulse.voltage))) 
                                * 2 * snr, (numpy.round(delay2.flatten() / impulse.dt)).astype(numpy.int)).reshape((len(trigger_sectors), len(ring_map), len(impulse.voltage)))  
    '''
    #print(phi, el)
    #print ('impulse time length: ', len(impulse.time))
    # Downsample if requested
    if downsample:
        #print ('impulse time length: ', len(impulse.time))
        trigger_waves, new_time = downsamplePayload(impulse.time, trigger_waves)
        #impulse.time = new_time  # for any later code that reads impulse.time
        new_dt = new_time[1] - new_time[0]
        #print(len(trigger_waves[0][0]))
        # Debug print of delays in original vs down‑sampled indices
        '''
        print("Delay sample counts (original vs downsampled):")
        for ant_key, raw_delay in delay[0]['delays'].items():
            orig_samples = int(round(raw_delay / impulse.dt))
            down_samples = int(round(raw_delay / new_dt))
            print(f"  Ant {ant_key}: orig {orig_samples}, down {down_samples}")
        '''
        # Replace the impulse time with downsampled time
        timebase = new_time
    else:
        timebase = impulse.time
    # Make sure we use the (possibly downsampled) time vector

    #n_samples = trigger_waves.shape[2]  # should now equal len(timebase)
    #print("timebase len:", len(timebase), "waveform len:", n_samples)
    
    if plot:
        # Plot waveforms for each antenna
        n_ant = len(antennas)
        fig, axes = plt.subplots(n_ant, 1, figsize=(10, 3*n_ant), sharex=True)
        if n_ant == 1:
            axes = [axes]
        
        for ant_idx in range(n_ant):
            axes[ant_idx].plot(timebase, trigger_waves[ant_idx], '''b-''', lw=0.8)
            axes[ant_idx].set_title(f"Antenna {antennas[ant_idx]}")
            axes[ant_idx].set_ylabel("Voltage [V]")
            axes[ant_idx].grid(True, alpha=0.3)
        
        axes[-1].set_xlabel("Time [ns]")
        fig.suptitle(f"Antenna Waveforms: φ = {phi}°, θ = {el}°", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        import os
        os.makedirs('''plots''', exist_ok=True)
        plt.savefig("plots/antenna_waveforms.png", dpi=150)
        plt.show()



    return trigger_waves, timebase, multiplier

def getPayloadWaveforms_dualpol(phi, el, impulse, beam_patterns_h, beam_patterns_v, 
                                antennas=None, snr=1, noise=None, 
                                psi=45.0, plot=False):
    # Default 4 physical antennas if not specified
    if antennas is None:
        antennas = [0, 1, 2, 3]
    
    n_antennas = len(antennas)
    n_samples = len(impulse.voltage)
    
    # Compute polarization projection coefficients
    # For linear polarization at angle psi from H-pol axis:
    psi_rad = numpy.deg2rad(psi)
    pol_h = numpy.cos(psi_rad)  # Projection onto H-pol channel
    pol_v = numpy.sin(psi_rad)  # Projection onto V-pol channel
    
    # Threshold very small coefficients to exactly zero (avoid floating point precision issues)
    if numpy.abs(pol_h) < 1e-10:
        pol_h = 0.0
    if numpy.abs(pol_v) < 1e-10:
        pol_v = 0.0
    
    # Initialize output array for all 8 channels
    waveforms = numpy.zeros((8, n_samples))
    multipliers = numpy.zeros(8)
    
    # Get geometric delays (same for both polarizations at each location)
    delay = delays.getAllDelays([phi], [el], antennas=antennas)
    
    print("phi = {:.2f}, theta = {:.2f}".format(round(delay[0]["phi"], 2), round(delay[0]["theta"], 2)))
    
    # Generate waveforms for each physical location
    for ant_idx, ant_num in enumerate(antennas):
        # Get delay for this antenna
        ant_name = f'Ant{ant_num+1}'
        ant_delay = delay[0]['delays'][ant_name]
        delay_samples = int(numpy.round(ant_delay / impulse.dt))
        
        # All 8 channels share the same physical boresight pointing direction.
        # phi_ant/theta_ant are all zero and are intentionally NOT used here — the
        # H vs V pol distinction is captured entirely by the choice of beam pattern
        # function: beam_patterns_h and beam_patterns_v have their E/H planes swapped,
        # which correctly models the 90° rotation of the V-pol element around boresight.
        # Geometric delay differences (from x/y/z positions) are handled by delays.py.
        phi_interp   = numpy.clip(phi, -180, 180)
        theta_interp = numpy.clip(el, -90, 90)
        
        # H-polarization waveform (channels 0-3, aligned at phi=0°)
        # H-pol: E-plane spans azimuth (narrow), H-plane spans elevation (wide)
        # NOTE: impulse is normalized to Vpp=1. snr parameter is kept for legacy callers
        # but snr_scan_dualpol.py sets snr=1 and applies its own amplitude scaling externally.
        waveforms[ant_num] = numpy.roll(impulse.voltage * snr, delay_samples) * \
            pol_h * \
            dBtoVoltsAtten(beam_patterns_h[0](phi_interp)) * \
            dBtoVoltsAtten(beam_patterns_h[1](theta_interp))
        
        # Store H-pol beam attenuation
        multipliers[ant_num] = pol_h * \
            dBtoVoltsAtten(beam_patterns_h[0](phi_interp)) * \
            dBtoVoltsAtten(beam_patterns_h[1](theta_interp))
        
        # V-polarization waveform (channels 4-7, aligned at phi=90°)
        # V-pol: E-plane spans elevation (narrow), H-plane spans azimuth (wide)
        # Pattern indices are swapped because antenna is rotated 90°
        waveforms[ant_num + 4] = numpy.roll(impulse.voltage * snr, delay_samples) * \
            pol_v * \
            dBtoVoltsAtten(beam_patterns_v[0](theta_interp)) * \
            dBtoVoltsAtten(beam_patterns_v[1](phi_interp))
        
        # Store V-pol beam attenuation
        multipliers[ant_num + 4] = pol_v * \
            dBtoVoltsAtten(beam_patterns_v[0](theta_interp)) * \
            dBtoVoltsAtten(beam_patterns_v[1](phi_interp))
        
        # Add noise if provided
        if noise is not None:
            waveforms[ant_num] += noise[ant_idx]


    timebase = impulse.time
    
    if plot:
        # Create side-by-side plots for H-pol and V-pol
        fig, axes = plt.subplots(n_antennas, 2, figsize=(12, 3*n_antennas))
        if n_antennas == 1:
            axes = axes.reshape(1, -1)
        
        for ant_idx in range(n_antennas):
            ant_num = antennas[ant_idx]
            # H-pol plot (channels 0-3)
            axes[ant_idx, 0].plot(timebase, waveforms[ant_num], 'b-', lw=0.8)
            axes[ant_idx, 0].set_title(f"Ant {ant_num} H-pol (Ch{ant_num})")
            axes[ant_idx, 0].set_ylabel("Voltage [V]")
            axes[ant_idx, 0].grid(True, alpha=0.3)
            
            # V-pol plot (channels 4-7)
            axes[ant_idx, 1].plot(timebase, waveforms[ant_num + 4], 'r-', lw=0.8)
            axes[ant_idx, 1].set_title(f"Ant {ant_num} V-pol (Ch{ant_num+4})")
            axes[ant_idx, 1].grid(True, alpha=0.3)
        
        axes[-1, 0].set_xlabel("Time [ns]")
        axes[-1, 1].set_xlabel("Time [ns]")
        fig.suptitle(f"Dual-Pol Waveforms: φ = {phi}°, θ = {el}°, ψ = {psi}°", fontsize=16)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        import os
        os.makedirs('plots', exist_ok=True)
        plt.savefig("plots/dualpol_waveforms.png", dpi=150)
        plt.show()
    
    return waveforms, timebase, multipliers

def downsamplePayload(time, trigger_waves):
    '''
    Downsample waveforms from upsampled rate to RITC sampling rate (250 ps)
    
    Args:
        time: time array (ns)
        trigger_waves: (num_antennas, num_samples) array
    
    Returns:
        downsampled trigger_waves and time arrays
    '''
    decimate_factor = int(aso_geometry.ritc_sample_step/((time[1]-time[0])))
    #print('decimate factor in downsamplePayload = ' ,decimate_factor)
    #print('time length pre decimation ', len(time))
    trigger_waves = trigger_waves[:,::decimate_factor]  # Updated for 2D array [antennas, samples]
    time = time[::decimate_factor]
    #print('time length post decimation ', len(time))
    #print('trigger_waves length post decimation ', len(trigger_waves[0]))
    return trigger_waves, time

if __name__ =='__main__':
    import matplotlib.pyplot as plt
    import sys
    import argparse
    sys.path.insert(0, '../')
    import coherent_sum as csum
    import tools.filters as filters

    def _save_axes_individually(fig, base_name, plots_dir='plots'):
        """Save each titled axes in fig as its own PNG, stripping step-number prefixes."""
        import re, os
        from matplotlib.transforms import Bbox

        os.makedirs(plots_dir, exist_ok=True)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        dpi = fig.get_dpi()
        pad_px = 14

        for ax in fig.get_axes():
            title = ax.get_title()
            if not title:
                continue
            # Strip leading step prefix: "1. ", "2a. ", "10d. ", etc.
            clean_title = re.sub(r'^\d+[a-zA-Z]?\.\s+', '', title)
            # Filename: use first line only, sanitize to alphanumeric/underscore
            fname = re.sub(r'[^a-zA-Z0-9]+', '_',
                           clean_title.split('\n')[0]).strip('_').lower()
            out_path = os.path.join(plots_dir, f'{base_name}_{fname}.png')
            bbox_disp = ax.get_tightbbox(renderer)
            if bbox_disp is None:
                continue
            bbox_in = Bbox([
                [(bbox_disp.x0 - pad_px) / dpi, (bbox_disp.y0 - pad_px) / dpi],
                [(bbox_disp.x1 + pad_px) / dpi, (bbox_disp.y1 + pad_px) / dpi],
            ])
            ax.set_title(clean_title, fontweight='bold')
            fig.savefig(out_path, dpi=150, bbox_inches=bbox_in)
            ax.set_title(title, fontweight='bold')
            print(f'  Saved: {out_path}')

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Dual-pol signal processing visualization')
    parser.add_argument('--phi', type=float, default=20, help='Azimuth angle (deg)')
    parser.add_argument('--theta', type=float, default=10, help='Elevation angle (deg)')
    parser.add_argument('--psi', type=float, default=45, help='Polarization angle (deg): 0=H, 90=V')
    parser.add_argument('--plots', type=str, default='all', 
                       help='Comma-separated plot rows to show (1-6) or "all". E.g., "1,3,5" or "4,5,6"')
    parser.add_argument('--no-save', action='store_true', help='Do not save plot to file')
    parser.add_argument('--second-filter', dest='second_filter', action='store_true', default=None,
                       help='Enable second 750 MHz lowpass filter (overrides module APPLY_SECOND_FILTER)')
    parser.add_argument('--no-second-filter', dest='second_filter', action='store_false',
                       help='Disable second filter (overrides module APPLY_SECOND_FILTER)')
    parser.add_argument('--snr', type=float, default=None,
                       help='Add independent white noise to each antenna at this SNR '
                            '(Vpp/2\u03c3 per channel). Omit for noise-free chain.')
    parser.add_argument('--save-subplots', action='store_true',
                       help='Save every subplot panel as an individual PNG '
                            '(strips step-number prefixes from titles, for LaTeX subfigure use).')
    args = parser.parse_args()

    # Resolve toggle: CLI flag beats module constant
    apply_second_filter = APPLY_SECOND_FILTER if args.second_filter is None else args.second_filter
    
    # Parse plot selection
    if args.plots == 'all':
        show_rows = [1, 2, 3, 4, 5, 6]
    else:
        show_rows = sorted([int(x.strip()) for x in args.plots.split(',')])
        show_rows = [r for r in show_rows if 1 <= r <= 6]
    
    phi = args.phi
    el = args.theta
    psi = args.psi
    
    print("\n" + "="*70)
    print("DUAL-POL SIGNAL PROCESSING CHAIN")
    print("="*70)
    print(f"Sky direction: phi={phi}°, theta={el}°")
    print(f"Polarization: psi={psi}° (0=H-pol, 90=V-pol)")
    print(f"Showing plots: rows {show_rows}")
    print(f"Second filter (750 MHz): {'ENABLED' if apply_second_filter else 'DISABLED'}")
    print("="*70 + "\n")
    
    # ========== STEP 1: RAW IMPULSE ==========
    print("STEP 1: Load raw impulse response")
    impulse = loadImpulse()
    impulse = prepImpulse(impulse)
    print(f"  Impulse: {len(impulse.voltage)} samples, dt={impulse.dt:.6f} ns")
    print(f"  Vpp = {numpy.max(impulse.voltage) - numpy.min(impulse.voltage):.3f} (normalized)")
    
    # ========== STEP 2: ANTENNA WAVEFORMS (TIME DOMAIN -  ) ==========
    print("\nSTEP 2: Generate antenna waveforms with beam patterns and delays")
    eplane_h = beamPattern(plot=False, which_plane='E', which_pol='H')
    hplane_h = beamPattern(plot=False, which_plane='H', which_pol='H')
    eplane_v = beamPattern(plot=False, which_plane='E', which_pol='V')
    hplane_v = beamPattern(plot=False, which_plane='H', which_pol='V')
    beam_patterns_h = (eplane_h, hplane_h)
    beam_patterns_v = (eplane_v, hplane_v)
    
    waveforms, timebase, multipliers = getPayloadWaveforms_dualpol(
        phi, el, impulse, beam_patterns_h, beam_patterns_v, 
        antennas=[0, 1, 2, 3], snr=5, noise=None, 
        psi=psi, plot=False
    )
    print(f"  Generated 8 channels (4 H-pol + 4 V-pol)")
    print(f"  All waveforms are   voltages (physical measurements)")
    print(f"  Beam multipliers (pol × beam pattern):")
    print(f"    H-pol: {multipliers[0]:.4f}, {multipliers[1]:.4f}, {multipliers[2]:.4f}, {multipliers[3]:.4f}")
    print(f"    V-pol: {multipliers[4]:.4f}, {multipliers[5]:.4f}, {multipliers[6]:.4f}, {multipliers[7]:.4f}")

    # ========== NOISE OVERLAY (optional) ==========
    _noise_label = 'no noise'
    if args.snr is not None:
        print(f"\nNOISE OVERLAY: SNR = {args.snr:.2f}  (Vpp / 2\u03c3 per antenna)")
        _rng = numpy.random.default_rng()
        _vpp_vals = []
        for _ch in range(8):
            _vpp = numpy.max(waveforms[_ch]) - numpy.min(waveforms[_ch])
            if _vpp < 1e-10:
                # zero-signal channel (e.g. H-pol when psi=90°) — skip noise
                continue
            _vrms = _vpp / (2.0 * args.snr)
            waveforms[_ch] = waveforms[_ch] + _rng.normal(0.0, _vrms, waveforms.shape[1])
            _vpp_vals.append(_vpp)
            print(f"  Ch{_ch}: Vpp={_vpp:.4f}  vrms_noise={_vrms:.4f}")
        _noise_label = f'SNR={args.snr:.1f}'
        print(f"  White noise added; gets bandlimited by downstream filter.")

    # ========== STEP 3: COHERENT SUM (TIME DOMAIN -  ) ==========
    print("\nSTEP 3: Coherent sum with geometric delays")
    delays_ns = getRemappedDelays(phi, el, antennas=[0, 1, 2, 3])
    delays_q = delays.quantize_delays_ns(delays_ns)
    print(f"  Delays: {delays_ns}")
    print(f"  Quantized: {delays_q}")
    
    # Get pre-digitization coherent sums for visualization (STEP 3)
    coh_h, coh_v, tb = csum.coherentSum_dualpol(
        waveforms[:4], waveforms[4:], timebase, delays_q,
        downsample=False, channel_mask=[1,1,1,1], output='separate', apply_filter=False
    )
    print(f"✓ Step 3: Coherent sum (pre-digitization)")
    print(f"  H-sum: max={numpy.max(numpy.abs(coh_h)):.3e} V")
    print(f"  V-sum: max={numpy.max(numpy.abs(coh_v)):.3e} V")
    print(f"  Sampling rate: {1/(tb[1]-tb[0])*1e-9/1e9:.3f} GHz")
    
    # ========== STEPS 4-10: DIGITIZE, FILTER, AND CIRCULAR DECOMPOSITION ==========
    # Use coherentSum_dualpol() with hardware-correct signal chain and return intermediates
    print("\n✓ Steps 4-10: Hardware signal chain (digitize → filter → sum → circular)")
    lhcp, rhcp, tb_digitized, inter = csum.coherentSum_dualpol(
        waveforms[:4], waveforms[4:], timebase, delays_q,
        downsample=True, channel_mask=[1,1,1,1], output='circular',
        apply_filter=True, apply_second_filter=apply_second_filter,
        fc_second=FC_SECOND_FILTER, return_intermediates=True
    )
    
    # Extract intermediate values for visualization
    coh_h_digitized     = inter['coh_h_digitized']
    coh_v_digitized     = inter['coh_v_digitized']
    coh_h_first_filtered = inter['coh_h_first_filtered']
    coh_v_first_filtered = inter['coh_v_first_filtered']
    coh_h_filtered      = inter['coh_h_filtered']   # final (= first if no second filter)
    coh_v_filtered      = inter['coh_v_filtered']
    H_fft = inter['H_fft']
    V_fft = inter['V_fft']
    V_shifted = inter['V_shifted']
    LHCP_fft = inter['LHCP_fft']
    RHCP_fft = inter['RHCP_fft']
    freqs = inter['freqs']
    decimate_factor = inter['decimate_factor']
    fs_adc = inter['fs_adc']
    dt_adc = (tb_digitized[1] - tb_digitized[0]) * 1e-9  # ns to s
    
    # Get filter coefficients for visualization
    filter_coeffs = filters.get_Shannon_Whitaker_coeffs(fs=4e9)
    nfft_filter = 4096
    filter_fft_plot = numpy.fft.rfft(filter_coeffs, n=nfft_filter)
    freqs_filter = numpy.fft.rfftfreq(nfft_filter, d=dt_adc)
    filter_response_db = 20*numpy.log10(numpy.abs(filter_fft_plot) / numpy.max(numpy.abs(filter_fft_plot)))

    # Second filter response (750 MHz)
    filter2_coeffs = filters.get_Shannon_Whitaker_coeffs(fs=4e9, fc=FC_SECOND_FILTER)
    filter2_fft_plot = numpy.fft.rfft(filter2_coeffs, n=nfft_filter)
    filter2_response_db = 20*numpy.log10(numpy.abs(filter2_fft_plot) / numpy.max(numpy.abs(filter2_fft_plot)))
    # Combined (cascaded) response
    filter_combined_db = filter_response_db + filter2_response_db
    
    # Compute power spectra
    power_lhcp_freq = numpy.abs(LHCP_fft)**2
    power_rhcp_freq = numpy.abs(RHCP_fft)**2
    power_total_freq = power_lhcp_freq + power_rhcp_freq
    
    # Print summary
    print(f"  Decimation: {decimate_factor}×")
    print(f"  ADC rate: {fs_adc/1e9:.3f} GHz")
    print(f"  Digitized samples: {len(coh_h_digitized)}")
    print(f"  H-pol max: {numpy.max(numpy.abs(coh_h_digitized)):.3e} V")
    print(f"  V-pol max: {numpy.max(numpy.abs(coh_v_digitized)):.3e} V")
    print(f"  LHCP power: {numpy.sum(power_lhcp_freq):.3e}")
    print(f"  RHCP power: {numpy.sum(power_rhcp_freq):.3e}")
    print(f"  Pol ratio: {numpy.sum(power_lhcp_freq)/numpy.sum(power_rhcp_freq):.3f}")
    
    # ========== VISUALIZATION (all intermediates now available) ==========
    # Note: lhcp and rhcp are already time-domain signals from irfft
    print("\n✓ Step 10: IFFT already computed by coherentSum_dualpol()")
    print(f"  lhcp and rhcp are   time-domain signals")
    print(f"  These are reconstructed waveforms from circular spectra")
    print(f"  Power = signal^2")
    
    # ========================================================================
    # CREATE COMPREHENSIVE DIAGNOSTIC PLOT
    # ========================================================================
    
    n_rows = len(show_rows)
    fig = plt.figure(figsize=(32, 4*n_rows))
    gs = fig.add_gridspec(n_rows, 4, hspace=0.5, wspace=0.4)
    
    # Map row number (1-6) to subplot row index (0 to n_rows-1)
    row_map = {row_num: idx for idx, row_num in enumerate(show_rows)}
    
    # Row 1: Time-domain antenna signals
    if 1 in show_rows:
        r = row_map[1]
        ax1 = fig.add_subplot(gs[r, 0])
        ax1.plot(timebase, impulse.voltage, 'k-', linewidth=1.5)
        ax1.set_title('1. Raw Impulse', fontweight='bold')
        ax1.set_xlabel('Time [ns]')
        ax1.set_ylabel('Voltage [V]')
        ax1.grid(True, alpha=0.3)
        
        ax2 = fig.add_subplot(gs[r, 1])
    for i in range(4):
        ax2.plot(timebase, waveforms[i], label=f'Ch{i}', alpha=0.6)
    ax2.set_title('2a. H-pol Antennas', fontweight='bold')
    ax2.set_xlabel('Time [ns]')
    ax2.set_ylabel('Voltage [V]')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    ax3 = fig.add_subplot(gs[0, 2])
    for i in range(4):
        ax3.plot(timebase, waveforms[i+4], label=f'Ch{i+4}', alpha=0.6)
    ax3.set_title('2b. V-pol Antennas', fontweight='bold')
    ax3.set_xlabel('Time [ns]')
    ax3.set_ylabel('Voltage [V]')
    ax3.legend(loc='upper right', fontsize=8)
    ax3.grid(True, alpha=0.3)
    
    ax4 = fig.add_subplot(gs[0, 3])
    ax4.plot(tb, coh_h.real, 'b-', linewidth=2, label='H')
    ax4.plot(tb, coh_v.real, 'r-', linewidth=2, label='V')
    ax4.set_title('3. Coherent Sum', fontweight='bold')
    ax4.set_xlabel('Time [ns]')
    ax4.set_ylabel('Voltage [V]')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    
    # Row 2: Digitization and filtering
    if 2 in show_rows:
        r = row_map[2]
        ax5 = fig.add_subplot(gs[r, 0])
        ax5.plot(tb_digitized, coh_h_digitized, 'b-', linewidth=1.5, alpha=0.7, label='H digitized')
        ax5.plot(tb_digitized, coh_v_digitized, 'r-', linewidth=1.5, alpha=0.7, label='V digitized')
        ax5.set_title('4. Digitized @ 4 GHz', fontweight='bold')
        ax5.set_xlabel('Time [ns]')
        ax5.set_ylabel('Voltage [V]')
        ax5.legend(fontsize=8)
        ax5.grid(True, alpha=0.3)
        ax5.text(0.02, 0.98, f'ADC: {fs_adc/1e9:.1f} GHz\nDecimate: 1/{decimate_factor}', 
                 transform=ax5.transAxes, va='top', fontsize=8,
                 bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
        
        ax7 = fig.add_subplot(gs[r, 1])
        ax7.plot(freqs_filter/1e9, filter_response_db, 'b-', linewidth=2, label='Filter 1 (1.5 GHz)')
        if apply_second_filter:
            ax7.plot(freqs_filter/1e9, filter2_response_db, 'g-', linewidth=2, label=f'Filter 2 ({FC_SECOND_FILTER/1e9:.2f} GHz)')
            ax7.plot(freqs_filter/1e9, filter_combined_db, 'k--', linewidth=2, label='Cascaded')
        ax7.axhline(-3, color='red', linestyle='--', alpha=0.7, label='-3 dB')
        ax7.set_title('5a. Filter Frequency Response', fontweight='bold')
        ax7.set_xlabel('Frequency [GHz]')
        ax7.set_ylabel('Gain [dB]')
        ax7.set_xlim([0, 2])
        ax7.set_ylim([-60, 5])
        ax7.legend(fontsize=7)
        ax7.grid(True, alpha=0.3)
        
        ax8 = fig.add_subplot(gs[r, 2])
        H_max_ref = numpy.max(numpy.abs(numpy.fft.rfft(coh_h_digitized))) + 1e-20

        H_unfilt_db      = 20*numpy.log10((numpy.abs(numpy.fft.rfft(coh_h_digitized)) + 1e-20)/H_max_ref)
        H_filt1_db       = 20*numpy.log10((numpy.abs(numpy.fft.rfft(coh_h_first_filtered)) + 1e-20)/H_max_ref)
        H_filt_final_db  = 20*numpy.log10((numpy.abs(H_fft) + 1e-20)/H_max_ref)

        ax8.plot(freqs/1e9, H_unfilt_db,     'b:', linewidth=1,   label='H raw',         alpha=0.4)
        ax8.plot(freqs/1e9, H_filt1_db,      'b-', linewidth=1.5, label='H filt1 (1.5G)', alpha=0.8)
        if apply_second_filter:
            ax8.plot(freqs/1e9, H_filt_final_db, 'b-', linewidth=2.5, label='H filt2 (0.75G)', alpha=1.0)
        ax8.plot(freqs_filter/1e9, filter_response_db, 'gray', linewidth=1, linestyle='--', label='Filter 1', alpha=0.6)
        if apply_second_filter:
            ax8.plot(freqs_filter/1e9, filter_combined_db, 'k--', linewidth=1.5, label='Cascaded', alpha=0.7)
        ax8.set_title('5b. H-pol FFT: Filter Stages', fontweight='bold')
        ax8.set_xlabel('Frequency [GHz]')
        ax8.set_ylabel('Magnitude [dB]')
        ax8.set_xlim([0, 2])
        ax8.set_ylim([-60, 5])
        ax8.legend(fontsize=7, loc='upper right')
        ax8.grid(True, alpha=0.3)
        
        ax6 = fig.add_subplot(gs[r, 3])
        ax6.plot(tb_digitized, coh_h_first_filtered, 'b-',  linewidth=1.5, label='H filt1 (1.5 GHz)', alpha=0.7)
        ax6.plot(tb_digitized, coh_v_first_filtered, 'r-',  linewidth=1.5, label='V filt1 (1.5 GHz)', alpha=0.7)
        if apply_second_filter:
            ax6.plot(tb_digitized, coh_h_filtered, 'b-', linewidth=2.5, label='H filt2 (0.75 GHz)', alpha=1.0)
            ax6.plot(tb_digitized, coh_v_filtered, 'r-', linewidth=2.5, label='V filt2 (0.75 GHz)', alpha=1.0)
            ax6.set_title('6. Once- & Twice-Filtered', fontweight='bold')
        else:
            ax6.set_title('6. Filtered (1.5 GHz FIR)', fontweight='bold')
        ax6.set_xlabel('Time [ns]')
        ax6.set_ylabel('Voltage [V]')
        ax6.legend(fontsize=7)
        ax6.grid(True, alpha=0.3)
    
    # Row 3: FFT spectra for H and V polarizations
    if 3 in show_rows:
        r = row_map[3]
        ax_hfft = fig.add_subplot(gs[r, 0])
        ax_hfft.plot(freqs/1e9, 20*numpy.log10(numpy.abs(H_fft)/numpy.max(numpy.abs(H_fft))), 'b-', linewidth=2)
        ax_hfft.set_title('6a. H-pol FFT Spectrum', fontweight='bold')
        ax_hfft.set_xlabel('Frequency [GHz]')
        ax_hfft.set_ylabel('Magnitude [dB]')
        ax_hfft.set_xlim([0, 2])
        ax_hfft.set_ylim([-60, 5])
        ax_hfft.grid(True, alpha=0.3)
        ax_hfft.text(0.5, 0.9, 'FFT(H filtered)', 
                     transform=ax_hfft.transAxes, ha='center', fontsize=9,
                     bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        
        ax_vfft = fig.add_subplot(gs[r, 1])
        if numpy.max(numpy.abs(V_fft)) > 1e-10:  # Only plot if V has signal
            ax_vfft.plot(freqs/1e9, 20*numpy.log10(numpy.abs(V_fft)/numpy.max(numpy.abs(V_fft))), 'r-', linewidth=2)
        else:
            ax_vfft.axhline(0, color='gray', linestyle='--', alpha=0.5)
            ax_vfft.text(0.5, 0.5, 'No V-pol signal\n(psi=0°)', 
                         transform=ax_vfft.transAxes, ha='center', va='center',
                         fontsize=12, bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
        ax_vfft.set_title('6b. V-pol FFT Spectrum', fontweight='bold')
        ax_vfft.set_xlabel('Frequency [GHz]')
        ax_vfft.set_ylabel('Magnitude [dB]')
        ax_vfft.set_xlim([0, 2])
        ax_vfft.set_ylim([-60, 5])
        ax_vfft.grid(True, alpha=0.3)
        if numpy.max(numpy.abs(V_fft)) > 1e-10:
            ax_vfft.text(0.5, 0.9, 'FFT(V filtered)', 
                         transform=ax_vfft.transAxes, ha='center', fontsize=9,
                         bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
    
    # Row 4: Phase shift and circular decomposition
    if 4 in show_rows:
        r = row_map[4]
        ax9 = fig.add_subplot(gs[r, 0])
        # Phase of V before and after shift
        phase_V = numpy.angle(V_fft)
        phase_V_shifted = numpy.angle(V_shifted)
        ax9.plot(freqs/1e9, phase_V, 'r-', alpha=0.5, label='V filtered')
        ax9.plot(freqs/1e9, phase_V_shifted, 'r-', linewidth=2, label='V × j (π/2 shift)')
        ax9.set_title('7. Phase Shift (V-pol)', fontweight='bold')
        ax9.set_xlabel('Frequency [GHz]')
        ax9.set_ylabel('Phase [rad]')
        ax9.set_xlim([0, 2])
        ax9.legend(fontsize=8)
        ax9.grid(True, alpha=0.3)
        ax9.text(0.5, 0.5, 'Multiply by j\n= +90° phase', 
                 transform=ax9.transAxes, ha='center', va='center', fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='orange', alpha=0.5))
        
        ax10 = fig.add_subplot(gs[r, 1])
        ax10.plot(freqs/1e9, 20*numpy.log10(numpy.abs(LHCP_fft)/numpy.max(numpy.abs(LHCP_fft))), 'b-', linewidth=2)
        ax10.set_title('8a. LHCP Spectrum', fontweight='bold')
        ax10.set_xlabel('Frequency [GHz]')
        ax10.set_ylabel('Magnitude [dB]')
        ax10.set_xlim([0, 2])
        ax10.grid(True, alpha=0.3)
        ax10.text(0.5, 0.9, 'LHCP = (H + j×V)/√2', 
                 transform=ax10.transAxes, ha='center', fontsize=9,
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        
        ax11 = fig.add_subplot(gs[r, 2])
        ax11.plot(freqs/1e9, 20*numpy.log10(numpy.abs(RHCP_fft)/numpy.max(numpy.abs(RHCP_fft))), 'r-', linewidth=2)
        ax11.set_title('8b. RHCP Spectrum', fontweight='bold')
        ax11.set_xlabel('Frequency [GHz]')
        ax11.set_ylabel('Magnitude [dB]')
        ax11.set_xlim([0, 2])
        ax11.grid(True, alpha=0.3)
        ax11.text(0.5, 0.9, 'RHCP = (H - j×V)/√2', 
                  transform=ax11.transAxes, ha='center', fontsize=9,
                  bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
        
        ax12 = fig.add_subplot(gs[r, 3])
        ax12.semilogy(freqs/1e9, power_lhcp_freq, 'b-', linewidth=2, label='LHCP')
        ax12.semilogy(freqs/1e9, power_rhcp_freq, 'r-', linewidth=2, label='RHCP')
        ax12.set_title('9. Power Spectrum', fontweight='bold')
        ax12.set_xlabel('Frequency [GHz]')
        ax12.set_ylabel('Power [V²]')
        ax12.set_xlim([0, 2])
        ax12.legend()
        ax12.grid(True, alpha=0.3)
        ax12.text(0.5, 0.9, f'Ratio: {numpy.sum(power_lhcp_freq)/numpy.sum(power_rhcp_freq):.3f}', 
                  transform=ax12.transAxes, ha='center', fontsize=9,
                  bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    
    # ---- Power-sum + coincidence logic (used in Row 5 panel 10d) ----
    _ps_window, _ps_step = 160, 40                           # hardware parameters (samples)
    pw_lhcp_sum, _n_ps_fr = csum.powerSum(lhcp, window=_ps_window, step=_ps_step)
    pw_rhcp_sum, _         = csum.powerSum(rhcp, window=_ps_window, step=_ps_step)
    # Frame time axis: start of each sliding window
    _dt_ns   = (tb_digitized[1] - tb_digitized[0])           # ns per ADC sample
    _step_ns = _ps_step   * _dt_ns                           # 10 ns between frames
    _win_ns  = _ps_window * _dt_ns                           # 40 ns window width
    ps_frame_t = tb_digitized[0] + numpy.arange(_n_ps_fr) * _step_ns
    # Load threshold from JSON (same powerSum units, no n_ant normalization)
    import json as _json, os as _os
    _thr_lhcp = _thr_rhcp = None
    _thr_label = 'threshold'
    _json_path = 'noise/threshold_analysis_dualpol.json'
    if _os.path.exists(_json_path):
        with open(_json_path) as _jf:
            _jd = _json.load(_jf)
        for _jentry in _jd.get('thresholds', []):
            if abs(_jentry.get('target_global_rate_hz', -1) - 0.1) < 0.01:
                _thr_lhcp = _jentry['threshold_lhcp']
                _thr_rhcp = _jentry['threshold_rhcp']
                _thr_label = (f"T={_thr_lhcp:.2f} "
                              f"({_jd.get('n_beams',63)} beams, 0.1 Hz)")
                break
    if _thr_lhcp is None:
        _thr_lhcp = _thr_rhcp = 0.5 * max(pw_lhcp_sum.max(), pw_rhcp_sum.max())
        _thr_label = '50% peak (no JSON)'
    _coinc_l = pw_lhcp_sum >= _thr_lhcp
    _coinc_r = pw_rhcp_sum >= _thr_rhcp
    _coinc   = _coinc_l & _coinc_r

    # Row 5: Time domain circular signals (reconstructed from circular spectra)
    if 5 in show_rows:
        r = row_map[5]
        ax13 = fig.add_subplot(gs[r, 0])
        ax13.plot(tb_digitized, lhcp, 'b-', linewidth=1.5)
        ax13.set_title('10a. LHCP Time-Domain', fontweight='bold')
        ax13.set_xlabel('Time [ns]')
        ax13.set_ylabel('Voltage [V]')
        ax13.grid(True, alpha=0.3)
        #ax13.text(0.02, 0.98, 'Reconstructed from\nLHCP spectrum', 
        #          transform=ax13.transAxes, va='top', fontsize=8,
        #          bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
        
        ax14 = fig.add_subplot(gs[r, 1])
        ax14.plot(tb_digitized, rhcp, 'r-', linewidth=1.5)
        ax14.set_title('10b. RHCP Time-Domain', fontweight='bold')
        ax14.set_xlabel('Time [ns]')
        ax14.set_ylabel('Voltage [V]')
        ax14.grid(True, alpha=0.3)
        #ax14.text(0.02, 0.98, 'Reconstructed from\nRHCP spectrum', 
        #          transform=ax14.transAxes, va='top', fontsize=8,
        #          bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))
        
        ax15 = fig.add_subplot(gs[r, 2])
        ax15.plot(tb_digitized, lhcp, 'b-', linewidth=2, label='LHCP')
        ax15.plot(tb_digitized, rhcp, 'r-', linewidth=2, label='RHCP')
        ax15.set_title('10c. Overlaid Signals', fontweight='bold')
        ax15.set_xlabel('Time [ns]')
        ax15.set_ylabel('Voltage [V]')
        ax15.legend()
        ax15.grid(True, alpha=0.3)
    
        ax16 = fig.add_subplot(gs[r, 3])
        # Instantaneous power as light background traces
        ax16.plot(tb_digitized, lhcp**2, color='steelblue', linewidth=0.8,
                  alpha=0.25, label='_nolegend_')
        ax16.plot(tb_digitized, rhcp**2, color='firebrick',  linewidth=0.8,
                  alpha=0.25, label='_nolegend_')
        # Sliding-window power sum (step plot, one value per 10 ns step)
        # Append a closing point so the last frame gets its full visible width
        _ps_t_closed = numpy.append(ps_frame_t, ps_frame_t[-1] + _step_ns)
        _pw_l_closed = numpy.append(pw_lhcp_sum, pw_lhcp_sum[-1])
        _pw_r_closed = numpy.append(pw_rhcp_sum, pw_rhcp_sum[-1])
        ax16.step(_ps_t_closed, _pw_l_closed, where='post',
                  color='steelblue', linewidth=2.0, label='LHCP power sum')
        ax16.step(_ps_t_closed, _pw_r_closed, where='post',
                  color='firebrick',  linewidth=2.0, label='RHCP power sum')
        # Threshold line from noise analysis
        ax16.axhline(_thr_lhcp, color='k', linestyle='--', linewidth=1.5,
                     label=_thr_label)
        # Shade each 10 ns step interval by trigger state
        #   green  = both LHCP and RHCP above threshold (coincidence)
        #   blue   = LHCP only above threshold
        #   red    = RHCP only above threshold
        for _i, (_t, _cl, _cr, _cc) in enumerate(
                zip(ps_frame_t, _coinc_l, _coinc_r, _coinc)):
            _t1, _t2 = float(_t), float(_t) + _step_ns
            if _cc:
                ax16.axvspan(_t1, _t2, alpha=0.30, color='limegreen', zorder=-1)
            elif _cl:
                ax16.axvspan(_t1, _t2, alpha=0.20, color='steelblue', zorder=-1)
            elif _cr:
                ax16.axvspan(_t1, _t2, alpha=0.20, color='firebrick',  zorder=-1)
        # Annotation box
        ax16.text(0.02, 0.98,
                  f'Coinc frames : {int(_coinc.sum())}/{_n_ps_fr}\n'
                  f'LHCP trig    : {int(_coinc_l.sum())}\n'
                  f'RHCP trig    : {int(_coinc_r.sum())}\n'
                  f'Green=coinc  Blue/Red=single-pol',
                  transform=ax16.transAxes, va='top', fontsize=7,
                  bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
        ax16.set_title('10d. Power Sum & Coincidence', fontweight='bold')
        ax16.set_xlabel('Time [ns]')
        ax16.set_ylabel('Power [V²/sample]')
        ax16.legend(fontsize=7, loc='upper right')
        ax16.grid(True, alpha=0.3)
        # xlim spans exactly the frames present - no clipping, count and shading are consistent
        _t_end = ps_frame_t[-1] + _win_ns   # end of last window
        ax16.set_xlim([tb_digitized[0], _t_end])
    
    # Row 6: Summary and comparisons
    if 6 in show_rows:
        r = row_map[6]
        ax17 = fig.add_subplot(gs[r, 0:2])
        # Processing chain flow diagram
        ax17.axis('off')
        flow_text = """
        SIGNAL PROCESSING CHAIN SUMMARY:
        
        TIME DOMAIN ( ):
        1. Impulse → 2. Antenna waveforms → 3. Coherent sum
        ↓
        FREQUENCY DOMAIN (COMPLEX):
        4. FFT → 5. Filter → 6. Filtered spectrum → 7. Phase shift V by π/2
        ↓
        8. Circular decomposition: LHCP = (H + j×V)/√2, RHCP = (H - j×V)/√2
        ↓
        9. Power spectrum = |LHCP|², |RHCP|²
        ↓
        10. IFFT (optional) → Complex time-domain → Power = |signal|²
        """
        ax17.text(0.05, 0.95, flow_text, transform=ax17.transAxes, 
                  va='top', ha='left', fontsize=10, family='monospace',
                  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        ax18 = fig.add_subplot(gs[r, 2:4])
        # Summary statistics
        summary_text = f"""
        RESULTS SUMMARY:
        
        Input: φ={phi}°, θ={el}°, ψ={psi}° 
        
        Coherent Sum (  voltages):
          H-pol: {numpy.max(numpy.abs(coh_h.real)):.3e} V
          V-pol: {numpy.max(numpy.abs(coh_v.real)):.3e} V
        
        Circular Power (integrated over spectrum):
          LHCP: {numpy.sum(power_lhcp_freq):.3e} V²
          RHCP: {numpy.sum(power_rhcp_freq):.3e} V²
          Ratio: {numpy.sum(power_lhcp_freq)/numpy.sum(power_rhcp_freq):.3f}
        
        Physical Interpretation:
          • ψ=0°: Pure H-pol → LHCP = RHCP (equal power)
          • ψ=45°: Linear ±45° → Circular (LHCP or RHCP dominant)
          • ψ=90°: Pure V-pol → LHCP = RHCP (equal power)
        """
        ax18.text(0.05, 0.95, summary_text, transform=ax18.transAxes,
                  va='top', ha='left', fontsize=11, family='monospace',
                  bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        ax18.axis('off')
    
    _filt_label = '1.5G+750M' if apply_second_filter else '1.5G only'
    plt.suptitle(f'Complete Dual-Pol Signal Processing Chain: φ={phi}°, θ={el}°, ψ={psi}°  '
                 f'[filter: {_filt_label}  |  {_noise_label}]',
                 fontsize=16, fontweight='bold', y=0.995)
    
    if not args.no_save:
        plt.savefig('plots/dualpol_processing_chain.png', dpi=300, bbox_inches='tight')
        print(f"\n{'='*70}")
        print(f"Saved comprehensive diagnostic plot to:")
        print(f"  plots/dualpol_processing_chain.png")
        print(f"{'='*70}\n")

    if args.save_subplots:
        print("Saving processing chain subplots individually...")
        _save_axes_individually(fig, 'processing_chain')
        print("")

    # =========================================================================
    # AUXILIARY NOISE FIGURE — amplitude before and after filter stages
    # =========================================================================
    # Generate a short (4096 sample) white-noise reference at the ADC rate.
    # Only the first 256 samples are shown in the time-domain panels so the
    # plot stays readable.  The full 4096-point FFT gives clean spectral view.
    numpy.random.seed(42)
    n_noise       = 4096
    n_noise_show  = 256          # samples shown in time-domain panel
    noise_raw     = numpy.random.normal(0.0, 1.0, n_noise)
    t_noise_ns    = numpy.arange(n_noise) * (1.0 / fs_adc) * 1e9  # ns

    noise_filt1   = filters.apply_Shannon_Whitaker_filter(noise_raw)
    rms_raw       = numpy.sqrt(numpy.mean(noise_raw  **2))
    rms_filt1     = numpy.sqrt(numpy.mean(noise_filt1**2))

    if apply_second_filter:
        noise_filt2 = filters.apply_Shannon_Whitaker_filter(noise_filt1, fc=FC_SECOND_FILTER)
        rms_filt2   = numpy.sqrt(numpy.mean(noise_filt2**2))
        n_noise_cols = 3
    else:
        n_noise_cols = 2

    freqs_noise = numpy.fft.rfftfreq(n_noise, d=1.0/fs_adc) / 1e9  # GHz

    def psd_db(sig):
        S = numpy.abs(numpy.fft.rfft(sig))**2
        return 10*numpy.log10(S / numpy.max(S) + 1e-14)

    fig_noise, axes_noise = plt.subplots(2, n_noise_cols, figsize=(5*n_noise_cols, 7))
    fig_noise.suptitle('Noise Amplitude: Filter Stage Comparison\n'
                        f'(white noise, {fs_adc/1e9:.1f} GHz ADC, {n_noise} samples, '
                        f'showing first {n_noise_show} in time domain)',
                        fontsize=13, fontweight='bold')

    # ── Row 0: time-domain waveforms (short window) ──────────────────────────
    slice_t = slice(0, n_noise_show)
    t_short  = t_noise_ns[slice_t]

    axes_noise[0, 0].plot(t_short, noise_raw[slice_t], 'k-', linewidth=0.7)
    axes_noise[0, 0].set_title(f'Raw white noise\nRMS = {rms_raw:.3f}', fontweight='bold')
    axes_noise[0, 0].set_xlabel('Time [ns]'); axes_noise[0, 0].set_ylabel('Amplitude')
    axes_noise[0, 0].grid(True, alpha=0.3)
    axes_noise[0, 0].axhline( rms_raw, color='red',  linestyle='--', linewidth=1, label='+1σ')
    axes_noise[0, 0].axhline(-rms_raw, color='red',  linestyle='--', linewidth=1, label='-1σ')
    axes_noise[0, 0].legend(fontsize=8)

    axes_noise[0, 1].plot(t_short, noise_filt1[slice_t], 'b-', linewidth=0.7)
    axes_noise[0, 1].set_title(f'After Filter 1 (1.5 GHz)\nRMS = {rms_filt1:.3f}  '
                                f'({100*rms_filt1/rms_raw:.1f}% of raw)', fontweight='bold')
    axes_noise[0, 1].set_xlabel('Time [ns]'); axes_noise[0, 1].set_ylabel('Amplitude')
    axes_noise[0, 1].grid(True, alpha=0.3)
    axes_noise[0, 1].axhline( rms_filt1, color='red', linestyle='--', linewidth=1)
    axes_noise[0, 1].axhline(-rms_filt1, color='red', linestyle='--', linewidth=1)

    if apply_second_filter:
        axes_noise[0, 2].plot(t_short, noise_filt2[slice_t], 'g-', linewidth=0.7)
        axes_noise[0, 2].set_title(f'After Filter 2 ({FC_SECOND_FILTER/1e9:.3f} GHz)\n'
                                    f'RMS = {rms_filt2:.3f}  ({100*rms_filt2/rms_raw:.1f}% of raw)',
                                    fontweight='bold')
        axes_noise[0, 2].set_xlabel('Time [ns]'); axes_noise[0, 2].set_ylabel('Amplitude')
        axes_noise[0, 2].grid(True, alpha=0.3)
        axes_noise[0, 2].axhline( rms_filt2, color='red', linestyle='--', linewidth=1)
        axes_noise[0, 2].axhline(-rms_filt2, color='red', linestyle='--', linewidth=1)

    # ── Row 1: power spectral density ────────────────────────────────────────
    axes_noise[1, 0].plot(freqs_noise, psd_db(noise_raw), 'k-', linewidth=0.8)
    axes_noise[1, 0].set_title('PSD — Raw', fontweight='bold')
    axes_noise[1, 0].set_xlabel('Frequency [GHz]'); axes_noise[1, 0].set_ylabel('Power [dB, norm.]')
    axes_noise[1, 0].set_xlim([0, 2]); axes_noise[1, 0].set_ylim([-40, 5])
    axes_noise[1, 0].grid(True, alpha=0.3)

    axes_noise[1, 1].plot(freqs_noise, psd_db(noise_filt1), 'b-', linewidth=0.8)
    axes_noise[1, 1].plot(freqs_filter/1e9, filter_response_db, 'r--', linewidth=1.5,
                           label='Filter 1 resp.', alpha=0.7)
    axes_noise[1, 1].set_title('PSD — After Filter 1', fontweight='bold')
    axes_noise[1, 1].set_xlabel('Frequency [GHz]'); axes_noise[1, 1].set_ylabel('Power [dB, norm.]')
    axes_noise[1, 1].set_xlim([0, 2]); axes_noise[1, 1].set_ylim([-40, 5])
    axes_noise[1, 1].legend(fontsize=8); axes_noise[1, 1].grid(True, alpha=0.3)

    if apply_second_filter:
        axes_noise[1, 2].plot(freqs_noise, psd_db(noise_filt2), 'g-', linewidth=0.8)
        axes_noise[1, 2].plot(freqs_filter/1e9, filter_combined_db, 'r--', linewidth=1.5,
                               label='Cascaded resp.', alpha=0.7)
        axes_noise[1, 2].set_title('PSD — After Filter 2', fontweight='bold')
        axes_noise[1, 2].set_xlabel('Frequency [GHz]'); axes_noise[1, 2].set_ylabel('Power [dB, norm.]')
        axes_noise[1, 2].set_xlim([0, 2]); axes_noise[1, 2].set_ylim([-40, 5])
        axes_noise[1, 2].legend(fontsize=8); axes_noise[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()
    if not args.no_save:
        plt.savefig('plots/noise_filter_comparison.png', dpi=150, bbox_inches='tight')
        print(f"Saved noise comparison plot to: plots/noise_filter_comparison.png")

    if args.save_subplots:
        print("Saving noise filter subplots individually...")
        _save_axes_individually(fig_noise, 'noise_filter')
        print("")

    # Print RMS reduction summary
    print(f"\nNoise RMS (white input, {fs_adc/1e9:.1f} GHz):")
    print(f"  Raw         : {rms_raw:.4f}  (100%)")
    print(f"  After filt1 : {rms_filt1:.4f}  ({100*rms_filt1/rms_raw:.1f}%)")
    if apply_second_filter:
        print(f"  After filt2 : {rms_filt2:.4f}  ({100*rms_filt2/rms_raw:.1f}%)")

    plt.show()
