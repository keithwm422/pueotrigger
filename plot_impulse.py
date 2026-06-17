import numpy
import myplot
import matplotlib.pyplot as plt
import tools.constants as constants
import tools.CoRaLs_geometry as aso_geometry
import payload_signal as payload
import noise
import math


directory='impulse/'
files = []
files.append(directory+'corals_impulse_sci.txt')
#files.append(directory+'corals_impulse_downsampled.txt') #Patrick's impulse respones
#files.append(directory+'corals_impulse_downsampled1.txt') #Patrick's impulse respones
#files.append(directory+'corals_impulse_sci.txt') 
#files.append(directory+'triggerA3.txt') #A3 trig
#files.append(directory+'mcm_trigger_path_pulse.txt') #ejo mcm A4 tests
#files.append(directory+'mcm_input_pulse.txt') #ejo mcm A4 tests
#dat = numpy.loadtxt('impulse/corals_impulse_sci.txt')
#print(dat.shape)
#print(dat[:5])
#dat = numpy.loadtxt('impulse/corals_impulse.txt')
#print(dat.shape)
#print(dat[:5])
label=['Corals IR','Corals IR 0.05 downsample', 'A3 trig. IR', 'A4 trigtest','A4 trigtest, input ']

for i,f in enumerate(files):
    impulse = payload.loadImpulse(f)
    impulse = payload.prepImpulse(impulse)

    print("Voltage min/max after prepImpulse:", numpy.min(impulse.voltage), numpy.max(impulse.voltage))
    print("First 10 voltage samples:", impulse.voltage[:10])

    print('timestep of input pulse [ns]:', impulse.dt)

    plt.figure(1)
    plt.plot(impulse.time, impulse.voltage, label=label[i])
    plt.xlabel('Time [ns]')
    plt.ylabel('amplitude [arb V]')
    
    # --- FFT plot ---
    impulse.fft()  # Ensure FFT is up to date
    abs_fft = numpy.abs(impulse.ampl)
    eps = 1e-20

    y_vals = 20 * numpy.log10((abs_fft + eps) / (numpy.max(abs_fft) + eps)) + 3
    # Find the maximum frequency where y > -80 dB
    valid_indices = numpy.where(y_vals > -80)[0]
    if valid_indices.size > 0:
        max_freq = impulse.freq[valid_indices[-1]]
    else:
        max_freq = 5  # fallback if all values are below -80 dB

    plt.figure(2)
    plt.plot(impulse.freq, y_vals, label=label[i])
    plt.grid(True)
    plt.xlim([0, max_freq])
    plt.ylim([-80, 5])
    plt.xlabel('Freq [GHz]')
    plt.ylabel('amplitude [dB]')

#noise profile 2:
thermal_noise = noise.ThermalNoise(0.26, 0.95, filter_order=(10,10), v_rms=1.0, 
                                   fbins=2**12, 
                                   time_domain_sampling_rate=0.01)

plt.plot(thermal_noise.frequencies, 20*numpy.log10(thermal_noise.amplitudes/numpy.max(thermal_noise.amplitudes)), '--', c='red', label='sim noise')

plt.legend()

plt.figure(1)
plt.legend()
plt.show()
