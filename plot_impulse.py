import numpy
import myplot
import matplotlib.pyplot as plt
import tools.constants as constants
import tools.aso_geometry as aso_geometry
import payload_signal as payload
import noise
import math


directory='impulse/'
files = []
files.append(directory+'triggerTF_02TH.txt') #A4 trig, used in simulation studies
files.append(directory+'triggerA3.txt') #A3 trig
#files.append(directory+'mcm_trigger_path_pulse.txt') #ejo mcm A4 tests
#files.append(directory+'mcm_input_pulse.txt') #ejo mcm A4 tests

label=['A4 trig. IR', 'A3 trig. IR', 'A4 trigtest','A4 trigtest, input ']

for i,f in enumerate(files):
    impulse = payload.loadImpulse(f)
    impulse = payload.prepImpulse(impulse)

    print('timestep of input pulse [ns]:', impulse.dt)

    plt.figure(1)
    plt.plot(impulse.time, impulse.voltage, label=label[i])
    plt.xlabel('Time [ns]')
    plt.ylabel('amplitude [arb V]')
    
    plt.figure(2)
    abs_fft = numpy.abs(impulse.ampl)
    plt.plot(impulse.freq, 20*numpy.log10(abs_fft/numpy.max(abs_fft))+3, label=label[i])
    plt.grid(True)
    plt.xlim([0,2])
    plt.ylim([-25,5])
    plt.xlabel('Freq [GHz]')
    plt.ylabel('amplitude [dB]')


#noise profile 2:
thermal_noise = noise.ThermalNoise(0.26, 1.05, filter_order=(10,10), v_rms=1.0, 
                                   fbins=2**12, 
                                   time_domain_sampling_rate=0.01)

plt.plot(thermal_noise.frequencies, 20*numpy.log10(thermal_noise.amplitudes/numpy.max(thermal_noise.amplitudes)), '--', c='red', label='sim noise')

plt.legend()

plt.figure(1)
plt.legend()
plt.show()
