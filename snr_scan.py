import numpy
import myplot
import matplotlib.pyplot as plt
import tools.aso_geometry as aso_geometry
import tools.constants as constants
import payload_signal as payload
import coherent_sum as trigger
import noise
import copy
import tools.waveform as waveform

#pick phi, theta (incoming wave angle)
#phi = 15
#theta = -10
phi=numpy.arange(-5,51,10)
theta=numpy.arange(-40,10,5)

## angle scan not yet implemented
#theta_range = numpy.arange(-40, 10, 4)
#phi_range = numpy.arange(0, 51,10) 
def gimmeWaveformIndex(myphi,mytheta,phi,theta):
    #myphi=15
    #mytheta=-10
    phi_index=numpy.argmin(numpy.abs(phi-myphi))
    theta_index=numpy.argmin(numpy.abs(theta-mytheta))
    #phi_i=0 and theta_i =0 is waveforms_list[0], phi_i=0 and theta_i=1 would be waveforms_list[1], phi_i=1 and theta_i=0 would be waveforms_list[len(theta)] so (phi_index(len_theta)) + theta_index
    #print(phi_index)
    #print('your phi: ', phi[phi_index])
    #print('phi index is: ', phi_index)
    #print('your theta: ', theta[theta_index])
    #print('theta index is: ', theta_index)
    #print('waveforms list length is: ',len(waveforms_list))
    waveforms_index=(phi_index*len(theta)) + theta_index
    return waveforms_index
    #print('waveform index: ', waveforms_index)
    #print(phi)
    #print(theta)
#def gimme_phi_theta(waveforms_i):
#    
def gimmeHits(noise_list, snr_scan, num_of_events_per_snr_step, waveforms_list,timebase_list,delays_list,ringmask,window,step,threshold,multiplier_list,waveforms_index,save_filename):
    waveforms=waveforms_list[waveforms_index]
    data_to_save_i=[]
    for j in range(len(noise_list)):
        hits=[]
        for snr in snr_scan:
            _hits = 0
            for i in range(num_of_events_per_snr_step):
                event_noise = numpy.reshape(numpy.real(noise_list[j][2][i*waveforms.shape[0]*waveforms.shape[1]:i*waveforms.shape[0]*waveforms.shape[1]+\
                                                                         waveforms.shape[0]*waveforms.shape[1]]), (waveforms.shape[0], waveforms.shape[1], waveforms.shape[2]))
                
                coh_sum, timebase_coh_sum  = trigger.coherentSum(waveforms_list[waveforms_index]*snr+event_noise, timebase_list[waveforms_index], delays_list[waveforms_index], ringmask=ringmask)
                #print(coh_sum, timebase_coh_sum)
                power,_ = trigger.powerSum(coh_sum/numpy.sqrt(waveforms.shape[0]*numpy.sum(ringmask)), window=window, step=step)
                if numpy.max(power) > threshold[j]:
                    _hits = _hits+1
            print (snr, float(_hits)/num_of_events_per_snr_step)
            hits.append(float(_hits)/num_of_events_per_snr_step)
        data_to_save_i.append(numpy.array(hits))
        data_to_save_i.append(numpy.array(hits)*(numpy.max(multiplier_list[waveforms_index])))
    numpy.savetxt(save_filename+str(waveforms_index)+'.txt', numpy.array(data_to_save_i))
    return hits
    
if __name__=='__main__':
    import sys

    if len(sys.argv) == 2:
        save_filename = sys.argv[1]
    else:
        save_filename = 'snr_scan_data'
    # testing sample rate in waveform init
    N_samples=4096
    freq_i=100*10**6 # 100 MHz
    omega=2*numpy.pi*freq_i
    #now omega*t is the argument in the sine wave, so we can get real_times from this
    period=1/freq_i
    print(period)
    # how many periods do you measure for, a random number maybe?
    rng = numpy.random.default_rng()
    num_periods=1+(rng.random()*3)
    print(period*num_periods)

    # need sampling rate now...
    sampling_period= 3*10**(-9)# lets do 500MHz or 50Mega samples / sec
    times=numpy.linspace(0,N_samples*sampling_period,N_samples)
    print(times)
    volty=numpy.sin(omega*times)
    real_times=numpy.linspace(0,N_samples*sampling_period,int(100*N_samples))
    print(real_times[-1])
    real_volty=numpy.sin(omega*real_times)
    reaL_siney=waveform.Waveform(real_volty,real_times) # what is the sampling rate now? if we had 4096 samples in 2 periods of a 100MHz sine wave then 

    #real_times=numpy.linspace(0,period,N_samples)
    #siney=waveform.Waveform(volty, sampling_rate=(sampling_period)) # what is the sampling rate now? if we had 4096 samples in 2 periods of a 100MHz sine wave then 
    siney=waveform.Waveform(volty, times) # what is the sampling rate now? if we had 4096 samples in 2 periods of a 100MHz sine wave then 

    #plt.scatter(siney.time,siney.voltage, label='waveform')
    plt.scatter(reaL_siney.time,reaL_siney.voltage, label='real wave')
    plt.scatter(times,volty, label='sampled wave')

    #plt.scatter(siney.time,siney.voltage, label='sampled wave')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.xlim([0,period])
    plt.show()
    quit()
    eplane, hplane = payload.beamPattern(plot=False)
    #first we have the antennas Eplane and Hplane decibels versus angle which can be graphed here
    print(type(eplane))

    impulse_d = payload.loadImpulse('impulse/triggerTF_02TH.txt',plot=True) # this is the freq db response of trigger. Can we load in the new PUEO antenna?
    impulse_X = copy.deepcopy(impulse_d)
    impulse_d = payload.prepImpulse(impulse_d,plot=True)
    #plt.scatter(impulse_X.time,impulse_X.voltage, label='loaded file')
    plt.scatter(impulse_d.time,impulse_d.voltage, label='waveform')
#       plt.plot(interp_angle, eplane_vpol_interp(interp_angle), label='vpolinterp Eplane')
#       plt.plot(interp_angle, hplane_vpol_interp(interp_angle), label='vpolinterp Hplane')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.show()

    quit()
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
    #delays=numpy.empty((phi.shape[0],theta.shape[0])) # declare shape based on theta and phi which are tuples
    # access like delays[iphi,itheta] for example delays[1,12]

    #waveforms=numpy.empty((phi.shape[0],theta.shape[0]))
    #timebase=numpy.empty((phi.shape[0],theta.shape[0]))
    #multiplier=numpy.empty((phi.shape[0],theta.shape[0]))
    impulse_list=[]
    delays_list=[]
    waveforms_list=[]
    timebase_list=[]
    multiplier_list=[]
    waveforms_all=[]
    # we can loop over all incoming wave angles.
    phi_i=0
    theta_i=0
    while phi_i<len(phi):
        theta_i=0
        while theta_i<len(theta):
            delays=payload.getRemappedDelays(phi[phi_i], theta[theta_i], phi_sectors) # generate the delays per phi sector, and for each incoming wave angle phi[phi_i] and theta[theta_i]
            #impulse should be copied over for every incoming wave, but getPayloadWaveforms overwrites it, so recopy to default...
            impulse=copy.deepcopy(impulse_d)
            waveforms, timebase, multiplier = payload.getPayloadWaveforms(phi[phi_i], theta[theta_i], phi_sectors, impulse, (eplane, hplane), plot=False, downsample=True) # find the waveform at the payload for each
            # must add to list for now
            print('waveforms shape: ', waveforms.shape[0], ' ', waveforms.shape[1], ' ', waveforms.shape[2])
            impulse_list.append(impulse)
            delays_list.append(delays)
            waveforms_list.append(waveforms)
            timebase_list.append(timebase)
            multiplier_list.append(multiplier)
            waveforms_all.append(gimmeWaveformIndex(phi[phi_i],theta[theta_i],phi,theta))
            theta_i+=1
        phi_i+=1
    # first pick phi_i and theta_i for testing...
    #waveforms_index2=gimmeWaveformIndex(15,-10,phi,theta)
    #waveforms_index1=gimmeWaveformIndex(25,-10,phi,theta)
    waveforms=waveforms_list[0]

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
                                       fbins=waveforms.shape[2], 
                                       time_domain_sampling_rate=aso_geometry.ritc_sample_step)

    noise_list=[] # list which now has entries that are different sizes...
    #noise_list.append(thermal_noise_v1.makeNoiseWaveform(ntraces=num_of_events_per_snr_step*waveforms.shape[0]*waveforms.shape[1]))
    noise_list.append(thermal_noise_v2.makeNoiseWaveform(ntraces=num_of_events_per_snr_step*waveforms.shape[0]*waveforms.shape[1]))
    hits=[]
    for wavey in waveforms_all: 
        hits.append(gimmeHits(noise_list, snr_scan, num_of_events_per_snr_step, waveforms_list,timebase_list,delays_list,ringmask,window,step,threshold,multiplier_list,wavey,save_filename))
    #hits1=gimmeHits(noise_list, snr_scan, num_of_events_per_snr_step, waveforms_list,timebase_list,delays_list,ringmask,window,step,threshold,multiplier_list,waveforms_index1,save_filename)
    #hits2=gimmeHits(noise_list, snr_scan, num_of_events_per_snr_step, waveforms_list,timebase_list,delays_list,ringmask,window,step,threshold,multiplier_list,waveforms_index2,save_filename)
    
    plt.plot(snr_scan, hits[0],label="phi25")
    plt.plot(snr_scan*(numpy.max(multiplier_list[0])), hits[0],label="phi25")
    plt.plot(snr_scan, hits[1],label="phi15")
    plt.plot(snr_scan*(numpy.max(multiplier_list[1])), hits[1],label="phi15")
    plt.legend()
    plt.ylim([-.1,1.1])
    plt.grid(True)    
    plt.show()
