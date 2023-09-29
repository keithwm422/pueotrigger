import numpy
import myplot
import matplotlib.pyplot as plt
import tools.aso_geometry as aso_geometry
import tools.constants as constants
import payload_signal as payload
import coherent_sum as trigger
import noise

#pick phi, theta (incoming wave angle)
#phi = 15
#theta = -10
phi=numpy.arange(0,51,10)
theta=numpy.arange(-40,10,4)

## angle scan not yet implemented
#theta_range = numpy.arange(-40, 10, 4)
#phi_range = numpy.arange(0, 51,10) 

if __name__=='__main__':
    import sys

    if len(sys.argv) == 2:
        save_filename = sys.argv[1]
    else:
        save_filename = 'snr_scan_data'
        
    eplane, hplane = payload.beamPattern()

    impulse = payload.loadImpulse('impulse/triggerTF_02TH.txt') # this is the freq db response of trigger. Can we load in the new PUEO antenna?
    impulse = payload.prepImpulse(impulse)

    ### pick phi sectors to include in trigger
    #phi_sectors = [2,3,4]
    phi_sectors = [1,2,3,4]
    #phi_sectors = [2]

    ### pick rings to include in trigger
    ringmask=[1,1,1,1]

    ### pick power sum window
    window = 32# number of samples in a window

    step=window/2

    ## number of events to throw per step
    num_of_events_per_snr_step = 200
    ## snr steps to scan
    snr_scan = numpy.arange(0.2, 4.0, 0.1)
    data_to_save = []
    data_to_save.append(snr_scan)
    ## thresholds picked from fits provided by generate_threshold_curves.py
    ## these are given as a list to match the number of noise profiles tested
    threshold = [4.2,4.0]
    #threshold = [6.8,6.3]
    #delays,waveforms,timebase,mulitplier will be new array, with phi size on one axis and theta on the other?
    print(phi.shape[0])
    print(type(theta.shape))
    delays=numpy.empty((phi.shape[0],theta.shape[0])) # declare shape based on theta and phi which are tuples
    # access like delays[iphi,itheta] for example delays[1,12] 
    waveforms=numpy.empty((phi.shape[0],theta.shape[0]))
    timebase=numpy.empty((phi.shape[0],theta.shape[0]))
    multiplier=numpy.empty((phi.shape[0],theta.shape[0]))
    # we can loop over all incoming wave angles.
    phi_i=0
    theta_i=0
    while phi_i<len(phi):
        theta_i=0
        while theta_i<len(theta):
            delays[phi_i,theta_i]=payload.getRemappedDelays(phi[phi_i], theta[theta_i], phi_sectors) # generate the delays per phi sector, and for each incoming wave angle phi[phi_i] and theta[theta_i]
            waveforms[phi_i,theta_i], timebase[phi_i,theta_i], multiplier[phi_i,theta_i] = payload.getPayloadWaveforms(phi[phi_i], theta[theta_i], phi_sectors, impulse, (eplane, hplane), plot=False, downsample=True) # find the waveform at the payload for each
            theta_i+=1
        phi_i+=1
    
    #delays = payload.getRemappedDelays(phi, theta, phi_sectors)
    # all of these also will need to get bigger to array sizes given by phi and theta sizes?
    #waveforms, timebase, multiplier = payload.getPayloadWaveforms(phi, theta, phi_sectors, impulse, (eplane, hplane), plot=False, downsample=True)

    # only need one thermal noise generatd for all angles (but version might be different, newer versions too?)
    # All waveforms have same waveforms.shape so just use first one
    ## noise profile 1
    #thermal_noise_v1 = noise.ThermalNoise(0.26, 0.9, filter_order=(10,6), v_rms=1.0, 
    #                                   fbins=waveforms[0,0].shape[2], 
    #                                   time_domain_sampling_rate=aso_geometry.ritc_sample_step)

    ## noise profile 2
    thermal_noise_v2 = noise.ThermalNoise(0.26, 1.05, filter_order=(10,10), v_rms=1.0, 
                                       fbins=waveforms[0,0].shape[2], 
                                       time_domain_sampling_rate=aso_geometry.ritc_sample_step)

    noise_list=[] # list which now has entries that are different sizes...
    #noise_list.append(thermal_noise_v1.makeNoiseWaveform(ntraces=num_of_events_per_snr_step*waveforms.shape[0]*waveforms.shape[1]))
    noise_list.append(thermal_noise_v2.makeNoiseWaveform(ntraces=num_of_events_per_snr_step*waveforms[0,0].shape[0]*waveforms[0,0].shape[1]))

    # first pick phi_i and theta_i for testing...
    myphi=15
    mytheta=-10
    phi_index=(numpy.abs(phi-myphi)).argmin
    theta_index=(numpy.abs(theta-mytheta)).argmin
    print('your phi: ', phi[phi_index])
    print('your theta: ', theta[theta_index])
    for j in range(len(noise_list)):
        hits=[]
        for snr in snr_scan:
            _hits = 0
            for i in range(num_of_events_per_snr_step):
                event_noise = numpy.reshape(numpy.real(noise_list[j][2][i*waveforms.shape[0]*waveforms.shape[1]:i*waveforms.shape[0]*waveforms.shape[1]+\
                                                                         waveforms.shape[0]*waveforms.shape[1]]), (waveforms.shape[0], waveforms.shape[1], waveforms.shape[2]))
                
                coh_sum, timebase_coh_sum  = trigger.coherentSum(waveforms*snr+event_noise, timebase, delays, ringmask=ringmask)
                #print(coh_sum, timebase_coh_sum)
                power,_ = trigger.powerSum(coh_sum/numpy.sqrt(waveforms.shape[0]*numpy.sum(ringmask)), window=window, step=step)

                if numpy.max(power) > threshold[j]:
                    _hits = _hits+1

            print (snr, float(_hits)/num_of_events_per_snr_step)
            hits.append(float(_hits)/num_of_events_per_snr_step)

        data_to_save.append(numpy.array(hits))
        data_to_save.append(numpy.array(hits)*(numpy.max(multiplier)))
                              
        plt.plot(snr_scan, hits)
        plt.plot(snr_scan*(numpy.max(multiplier)), hits)

    plt.ylim([-.1,1.1])
    plt.grid(True)


    numpy.savetxt(save_filename+'.txt', numpy.array(data_to_save))
    
    plt.show()
