import numpy as np
import noise
import tools.CoRaLs_geometry as corals_geometry

vrms = 1.0
#v2: noise profile '2'
thermal_noise = noise.ThermalNoise(0.26, .95, filter_order=(10,10), v_rms=1.0, 
                                   fbins=2**26, 
                                   time_domain_sampling_rate=corals_geometry.ritc_sample_step)
#####
#v1: noise profile '1'
#thermal_noise = noise.ThermalNoise(0.26, 0.9, filter_order=(10,6), v_rms=1.0, 
#                                   fbins=2**26, 
#                                   time_domain_sampling_rate=aso_geometry.ritc_sample_step)

noise = thermal_noise.makeNoiseWaveform(ntraces=10)
np.save('simulated_noise.npy', noise[2])