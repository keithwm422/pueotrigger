import numpy
import tools.waveform as waveform
import tools.delays as delays
import tools.CoRaLs_geometry as aso_geometry
from scipy import interpolate
from scipy.signal import lfilter, butter, cheby1
import matplotlib.pyplot as plt

#def h_plane_vpol_beam_pattern():
    
ring_map = {
    'B'  : 0,
    'T'  : 1,}
ring_map_inv = {v: k for k, v in ring_map.copy().items()}
    
def loadImpulse(filename='impulse/triggerTF_02TH.txt'):

    dat=numpy.loadtxt(filename)
    impulse=waveform.Waveform(dat[:,1], time=dat[:,0])
    return impulse

def prepImpulse(impulse, upsample=10, filter=True, highpass_cutoff=0.28, lowpass_cutoff=1.10 ):
    '''
    upsample, center, and filter impulse
    highpass_cutoff [GHz]
    '''
    impulse.zeropad(4096)
    #impulse.takeWindow([300, 1024+200])
    impulse.fft()
    impulse.upsampleFreqDomain(upsample)
    impulse.time = impulse.time-impulse.time[0]
    
    if filter:
        #highpass
        filtercoeff = cheby1(4, rp=0.5, Wn=highpass_cutoff/impulse.freq[-1], btype='highpass')
        impulse = waveform.Waveform(lfilter(filtercoeff[0], filtercoeff[1], impulse.voltage), 
                                    time=impulse.time)
        impulse.fft()
        '''
        #lowpass
        filtercoeff = cheby1(4, rp=0.5, Wn=lowpass_cutoff/impulse.freq[-1], btype='lowpass')
        impulse = waveform.Waveform(lfilter(filtercoeff[0], filtercoeff[1], impulse.voltage), 
                                    time=impulse.time)
        impulse.fft()
        '''
    #set Vpp = 1
    impulse.voltage = impulse.voltage / (numpy.max(impulse.voltage) - numpy.min(impulse.voltage))

    start = numpy.argmax(impulse.voltage)-3200
    impulse.takeWindow([start, start+10000])
    impulse.time = impulse.time-impulse.time[0]

    
    return impulse

#beamPattern interpolation from datafiles or function. For corals right now this is function not data.
# we should break this into 2 sections for E_plane and H_plane of antennas and then we can do beam pattern for them based on the plane and both phi and el for each plane...
# we use h_plane and e_plane as the same resp, just h_plane is rotated by 90 for resp with respect to el I guess...

# This function returns a E_plane beamPattern tuple that has both a el and phi resp (2D spline)
# don't really need a spline since its all a function currently, but eventually need a spline for real resp data included...
# This function returns a H_plane beamPattern tuple that has both a el and phi resp

def beamPattern(plot=False,which_plane='E', which_pol='V'):
    '''
    CoRaLS proxy beam pattern
    '''
    num_az=361
    az=numpy.linspace(-180,180,num_az)
    #print(len(az))
    #print(az)
    num_el=181
    el=numpy.linspace(-90,90,num_el)
    #el = -90:90;
    #fc = [150*10**6,300*10**6,450*10**6,600*10**6]; # units in MHz
    #fc2 = 150e6:50e6:850e6;
    #fc2=np.arange(150*10**6,850*10**6,50*10**6) # step size of 50MHz
    #c = physconst('Lightspeed');
    #c_speed=299792458 #m/s
    resp = 1*numpy.ones((num_el,num_az))
    el_i = 0 
    az_i = 0
    while el_i<len(el):
        az_i = 0
        while az_i<len(az):
            resp[el_i,az_i] = 10*numpy.log10(resp[el_i,az_i]*(numpy.cos(numpy.radians(el[el_i]))**(6))*numpy.cos(numpy.radians((az[az_i])/(2.0)))**(10))
            az_i+=1
        el_i+=1
    '''
    angle=
    angle = [-50,-40,-30,-20,-10,0,10,20,30,40,50]
    eplane_vpol = [-7.5, -5, -3, -1.5, -.5, 0, -.5, -1.5, -3, -5, -7.5]
    hplane_vpol = [-16, -11, -6, -3, -1, 0, -1, -3, -6, -11, -16]
    interp_angle = numpy.arange(-50,51,1)    
    eplane_vpol_interp = interpolate.interp1d(angle, eplane_vpol, kind='cubic')
    hplane_vpol_interp = interpolate.interp1d(angle, hplane_vpol, kind='cubic')
    '''
    if which_plane=='E':
        if which_pol=='V':
            # eplane and vpol response at a phi?
            plane_interp=interpolate.interp1d(el, resp[:,int(num_az/2)], kind='cubic')
        elif which_pol=='H':
            # eplane and hpol response at a el?
            # now this is obviously wrong but I am not sure what to put for the hpol response of the eplane of the antenna...
            # so for now stick with vpol responses...
            plane_interp=interpolate.interp1d(az, resp[int(num_el/2),:], kind='cubic')
    elif which_plane=='E1':
        plane_interp=interpolate.interp1d(el, resp[:,int(num_az/2)], kind='cubic')
    elif which_plane=='H':
        if which_pol=='V':
            # hplane and vpol response at a el seems to be what was happening for the antenna responses
            plane_interp=interpolate.interp1d(az, resp[int(num_el/2),:], kind='cubic')
        elif which_pol=='H':
            # hplane and hpol response at a az?
            # now this is obviously wrong but I am not sure what to put for the hpol response of the eplane of the antenna...
            # so for now stick with vpol responses...
            plane_interp=interpolate.interp1d(el, resp[:,int(num_az/2)], kind='cubic')
    elif which_plane=='H1':
        plane_interp = interpolate.interp1d(el, resp[:,int(num_az/2)], kind='cubic')
    #we can pick a phi I suppose for this plot. but I am not sure if this is phi or theta. It may have been symmetric?
    if plot:
        #plt.plot(interp_angle, eplane_vpol_interp(interp_angle), label='vpol Eplane')
        #plt.plot(interp_angle, hplane_vpol_interp(interp_angle), label='vpol Hplane')
        az_i=0
        while az_i < num_az:
            #plt.plot(el, resp[:,int(num_az/2)])
            plt.plot(el, resp[:,az_i])
            if az_i==int(num_az/2):
                plt.plot(el, resp[:,az_i],'--',label=str(which_plane))
            az_i+=10
        plt.legend(loc='upper left')
        plt.grid(True)
        plt.xlabel('off-boresight angle [deg.]')
        plt.ylabel('amplitude [dB]')
        plt.xlim([-50,50])
        plt.ylim([-10,1])

        plt.show()
    return plane_interp
def dBtoVoltsAtten(db_value):
    '''
    does what it says
    '''
    atten_fraction = 10**(db_value/20)
    return atten_fraction

def getRemappedDelays(phi, el, trigger_sectors):
    '''
    converts delays from delays.getAllDelays(), from a dict to a numpy array
    '''
    delay = delays.getAllDelays([phi], [el], phi_sectors=trigger_sectors) #gets delays at all antennas
    delays_remapped = numpy.zeros((len(trigger_sectors), 4))

    for i in delay[0]['delays']:

        delays_remapped[i[0]-numpy.min(trigger_sectors),ring_map[i[1]]] = delay[0]['delays'][i]

    return delays_remapped


def getPayloadDelays(phi,el,trigger_sectors):
    delay = delays.getAllDelays([phi], [el], phi_sectors=trigger_sectors) #gets delays at all antennas     
    print(type(delay))
    print(delay)
    return delay_all

def getSinglePayloadWaveform(phi, el, trigger_sectors, impulse, beam_pattern, snr=1, noise=None, plot=False, downsample=False):
    delay = delays.getAllDelays([phi], [el], phi_sectors=trigger_sectors) #gets delays at all antennas     
    print(type(delay))
    print(delay)
    #construct trigger_waves to keep the wfs that passed the trigger
    trigger_waves=numpy.zeros((len(trigger_sectors), 4, len(impulse.voltage)))
    print("trigger waves.shape is: {}".format(trigger_waves.shape))

    # the trigger waves isn't going to have the right shape now, because the code below expects each number index to be either a top or bottom, whereas corals doesn't have that kind of geometry

    #I dont know what multiplier does here yet
    multiplier=[]
    
    #not yet optimized for speed
    for delay_i in delay[0]['delays']:
        # calculate the phi and el and correct it to be within the interpolation range [(-90,90) or (-180,180)]
        phi_interp=phi-aso_geometry.phi_ant[i[0]-1]
        if phi_interp < -90:
            delta_phi=phi_interp+90
            phi_interp=-90-(phi_interp+90)
        elif phi_interp>90:
            phi_interp=-180

        # this line below: trigger_waves elements are calculated based on beam_pattern which is a tuple of eplane,hplane interp functions. beam_pattern[1] is hplane and [0] is eplane
        # hplane is evaluated with phi positions of antennas, and eplane is evaluated at el of antennas. The arguments to these are actually the "off-boresight" angle
        # so the argument is an angle that is the difference of incoming wave and the antenna's angle tilt in phi and theta 
        # the numpy.roll function will move the elements forward or backward along the axis, effectively forcing the delay to be taken into account for the impulse.
        # roll the impulse.voltage which is scaled by 2*snr with the delay: int(numpy.round(delay[0]['delays'][i] / impulse.dt))
        # we also multiply that impulse.voltage that was rolled by the beam_pattern, b/c the antenna direction and the angle of the impulse wave should be taken into acct
        # 
        #trigger waves is calculated with 2 voltage fractions multiplied in, one from eplane and one from hplane. Trigger wave must be on both?
        delay_i[0]-numpy.min(trigger_sectors)
        trigger_waves[i[0]-numpy.min(trigger_sectors),ring_map[i[1]]] = \
            numpy.roll(impulse.voltage * 2 * snr, int(numpy.round(delay[0]['delays'][i] / impulse.dt)))
        
        #there is a multiplier here that is hplane[dphi(phi)] and eplane[del(el)].
        #multiplier.append(dBtoVoltsAtten(beam_pattern[1](phi-aso_geometry.phi_ant[i[0]-1])) * \
        #    dBtoVoltsAtten(beam_pattern[0](el -aso_geometry.theta_ant[0])))
        

        #add in the noise if its included.
        if noise is not None:
            print("noise shape " + str(noise.shape))
            print("length of i in payload signal is: {}".format(len(i)))
            print("length of ring map in payload signal is: {}".format(len(ring_map)))
            trigger_waves[i[0]-numpy.min(trigger_sectors),ring_map[i[1]]] += \
                                noise[(i[0]-numpy.min(trigger_sectors))*len(ring_map) + ring_map[i[1]]]

    '''
    #numpyfied waveform generation. Factor of 2 multipler since impulse Vpp putatively normalized to = 1.0
    trigger_waves2 = numpy.roll(numpy.tile(impulse.voltage,len(trigger_sectors)*len(ring_map)).reshape((len(trigger_sectors)*len(ring_map), len(impulse.voltage))) 
                                * 2 * snr, (numpy.round(delay2.flatten() / impulse.dt)).astype(numpy.int)).reshape((len(trigger_sectors), len(ring_map), len(impulse.voltage)))  
    '''
    if downsample:
        trigger_waves, impulse.time = downsamplePayload(impulse.time, trigger_waves)
                       
    if plot:
        fig, ax = plt.subplots(len(ring_map), len(trigger_sectors))
        for i in range(len(trigger_sectors)):
            for j in range(len(ring_map)):
                if j != 0:
                    ax[len(ring_map)-j-1,i].set_xticklabels([])
                if i != 0:
                    ax[len(ring_map)-j-1,i].set_yticklabels([])

                ax[len(ring_map)-j-1,i].plot(impulse.time, trigger_waves[i,j], label=str(i)+ring_map_inv[j], c='black', lw=1, alpha=0.7)
                #ax[len(ring_map)-j-1,i].plot(impulse.time, trigger_waves2[i,j],  c='black', lw=1, alpha=0.7)

                ax[len(ring_map)-j-1,i].legend(loc='upper right')
                #ax[len(ring_map)-j-1,i].set_ylim([-snr-1,snr+1])

        plt.suptitle('phi = '+str(phi)+'deg.  theta = '+str(el)+'deg.', fontsize=20)
        #plt.tight_layout()
        plt.show()
    
    return trigger_waves, impulse.time, multiplier
def getPayloadWaveforms(phi, el, trigger_sectors, impulse, beam_pattern, snr=1, noise=None, plot=False, downsample=False):
    '''
    return waveforms for a single phi, el
    waveforms is 3d (N,M,P) array, where N=number trigger_sectors, M=number of rings, P=size of impulse

    typically, impulse is still upsampled here in order to align waveforms to a good 
    approximation of the true phi, el
    '''
    delay = delays.getAllDelays([phi], [el], phi_sectors=trigger_sectors) #gets delays at all antennas
    #delay2 = getRemappedDelays([phi], [el], trigger_sectors) #gets delays at all antennas
    print(type(delay))
    print(delay)
    #construct trigger_waves to keep the wfs that passed the trigger
    trigger_waves=numpy.zeros((len(trigger_sectors), 4, len(impulse.voltage)))
    print("trigger waves.shape is: {}".format(trigger_waves.shape))

    #I dont know what multiplier does here yet
    multiplier=[]
    
    #not yet optimized for speed
    for i in delay[0]['delays']:
        print("i[1] for delay is: {}".format(i[1]))
        print("ring map i[1] is: {}".format(ring_map[i[1]]))
        print("i[0] is {}".format(i[0]))
        print("numpy.min(trigger sectors) is: {}".format(numpy.min(trigger_sectors)))
        #,ring_map[i[1]]
        print("argument to noise: {}".format((i[0]-numpy.min(trigger_sectors))*len(ring_map) + ring_map[i[1]]))

        # calculate the phi and el and correct it to be within the interpolation range [(-90,90) or (-180,180)]
        phi_interp=phi-aso_geometry.phi_ant[i[0]-1]
        if phi_interp < -90:
            delta_phi=phi_interp+90
            phi_interp=-90-(phi_interp+90)
        elif phi_interp>90:
            phi_interp=-180

        # this line below: trigger_waves elements are calculated based on beam_pattern which is a tuple of eplane,hplane interp functions. beam_pattern[1] is hplane and [0] is eplane
        # hplane is evaluated with phi positions of antennas, and eplane is evaluated at el of antennas. The arguments to these are actually the "off-boresight" angle
        # so the argument is an angle that is the difference of incoming wave and the antenna's angle tilt in phi and theta 
        # the numpy.roll function will move the elements forward or backward along the axis, effectively forcing the delay to be taken into account for the impulse.
        # roll the impulse.voltage which is scaled by 2*snr with the delay: int(numpy.round(delay[0]['delays'][i] / impulse.dt))
        # we also multiply that impulse.voltage that was rolled by the beam_pattern, b/c the antenna direction and the angle of the impulse wave should be taken into acct
        # 
        #trigger waves is calculated with 2 voltage fractions multiplied in, one from eplane and one from hplane. Trigger wave must be on both?
            
        print("roll argument as integer for this trigger-wave-delay: {}".format(int(numpy.round(delay[0]['delays'][i] / impulse.dt))))
        trigger_waves[i[0]-numpy.min(trigger_sectors),ring_map[i[1]]] = \
            numpy.roll(impulse.voltage * 2 * snr, int(numpy.round(delay[0]['delays'][i] / impulse.dt))) * \
            dBtoVoltsAtten(beam_pattern[1](phi-aso_geometry.phi_ant[i[0]-1])) * \
            dBtoVoltsAtten(beam_pattern[0](el -aso_geometry.theta_ant[0]))
        
        #there is a multiplier here that is hplane[dphi(phi)] and eplane[del(el)].
        multiplier.append(dBtoVoltsAtten(beam_pattern[1](phi-aso_geometry.phi_ant[i[0]-1])) * \
            dBtoVoltsAtten(beam_pattern[0](el -aso_geometry.theta_ant[0])))
        

        #add in the noise if its included.
        if noise is not None:
            print("noise shape " + str(noise.shape))
            print("length of i in payload signal is: {}".format(len(i)))
            print("length of ring map in payload signal is: {}".format(len(ring_map)))
            trigger_waves[i[0]-numpy.min(trigger_sectors),ring_map[i[1]]] += \
                                noise[(i[0]-numpy.min(trigger_sectors))*len(ring_map) + ring_map[i[1]]]

    '''
    #numpyfied waveform generation. Factor of 2 multipler since impulse Vpp putatively normalized to = 1.0
    trigger_waves2 = numpy.roll(numpy.tile(impulse.voltage,len(trigger_sectors)*len(ring_map)).reshape((len(trigger_sectors)*len(ring_map), len(impulse.voltage))) 
                                * 2 * snr, (numpy.round(delay2.flatten() / impulse.dt)).astype(numpy.int)).reshape((len(trigger_sectors), len(ring_map), len(impulse.voltage)))  
    '''
    if downsample:
        trigger_waves, impulse.time = downsamplePayload(impulse.time, trigger_waves)
                                        
    if plot:
        fig, ax = plt.subplots(len(ring_map), len(trigger_sectors))
        for i in range(len(trigger_sectors)):
            for j in range(len(ring_map)):
                if j != 0:
                    ax[len(ring_map)-j-1,i].set_xticklabels([])
                if i != 0:
                    ax[len(ring_map)-j-1,i].set_yticklabels([])

                ax[len(ring_map)-j-1,i].plot(impulse.time, trigger_waves[i,j], label=str(i)+ring_map_inv[j], c='black', lw=1, alpha=0.7)
                #ax[len(ring_map)-j-1,i].plot(impulse.time, trigger_waves2[i,j],  c='black', lw=1, alpha=0.7)

                ax[len(ring_map)-j-1,i].legend(loc='upper right')
                #ax[len(ring_map)-j-1,i].set_ylim([-snr-1,snr+1])

        plt.suptitle('phi = '+str(phi)+'deg.  theta = '+str(el)+'deg.', fontsize=20)
        #plt.tight_layout()
        plt.show()
    
    return trigger_waves, impulse.time, multiplier


def downsamplePayload(time, trigger_waves):

    decimate_factor = int(aso_geometry.ritc_sample_step/((time[1]-time[0])))

    trigger_waves = trigger_waves[:,:,::decimate_factor]
    time = time[::decimate_factor]

    return trigger_waves, time                                  

def gimmePlotsOriginal(impulse):
            impulse.fft()
            plt.figure(1)
            plt.plot(impulse.time_zeroed, impulse.voltage, '-.', ms=2)
            plt.title("time domain?") 
            plt.figure(2)
            plt.plot(impulse.freq, numpy.abs(impulse.ampl), '-.', ms=2)
            impulse.upsampleFreqDomain(2)
            plt.figure(1)
            plt.plot(impulse.time_zeroed, impulse.voltage, '--', ms=1)
            #plt.xlim([0,60])
            plt.figure(2)
            plt.plot(impulse.freq, numpy.abs(impulse.ampl), '--', ms=1)   
            #plt.xlim([0,3])
            plt.title("frequency?") 
            #plt.figure(4)
            #plt.plot(impulse.time, numpy.correlate(impulse.voltage, impulse.voltage, "same"))
            plt.show()


def gimmePlots(impulse,impulse2):
            impulse.fft()
            impulse2.fft()
            plt.figure(1)
            plt.plot(impulse.time_zeroed, impulse.voltage, '-.', ms=2)
            plt.title("time domain?") 
            plt.figure(2)
            plt.plot(impulse.freq, numpy.abs(impulse.ampl), '-.', ms=2)
            #impulse.upsampleFreqDomain(2)
            plt.figure(1)
            plt.plot(impulse2.time, impulse2.voltage, '--', ms=1)
            #plt.xlim([0,60])
            plt.figure(2)
            plt.plot(impulse2.freq, numpy.abs(impulse2.ampl), '--', ms=1)   
            #plt.xlim([0,3])
            plt.title("frequency?") 
            #plt.figure(4)
            #plt.plot(impulse.time, numpy.correlate(impulse.voltage, impulse.voltage, "same"))
            plt.show()

def gimmePlotsImpulse(impulse):
            impulse.fft()
            plt.figure(1)
            plt.plot(impulse.time_zeroed, impulse.voltage, '-.', ms=2)
            plt.title("time domain")
            plt.figure(2)
            plt.plot(impulse.freq, numpy.abs(impulse.ampl), '-.', ms=2)
            plt.title("frequency") 
            plt.xlim([0,2])
            #plt.figure(4)
            #plt.plot(impulse.time, numpy.correlate(impulse.voltage, impulse.voltage, "same"))
            plt.show()

def gimmePlots2Impulses(impulse,impulse2):
            impulse.fft()
            impulse2.fft()
            fig, ax1 = plt.subplots()
            ax1.plot(impulse.time_zeroed, impulse.voltage, '-.', ms=2,color='red')
            fig.suptitle("time domain?") 
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            ax2.plot(impulse2.time, impulse2.voltage, '--', ms=1,color='blue')
            #ax1.set_ylim([-0.5,0.5])
            #ax2.set_ylim([-100,100])
            figB, ax1B = plt.subplots()
            ax1B.plot(impulse.freq, numpy.abs(impulse.ampl), '-.', ms=2,color='red')
            ax2B = ax1B.twinx()  # instantiate a second axes that shares the same x-axis
            ax2B.plot(impulse2.freq, numpy.abs(impulse2.ampl), '--', ms=1,color='blue')   
            #plt.xlim([0,3])
            figB.suptitle("frequency?") 
            #plt.figure(4)
            #plt.plot(impulse.time, numpy.correlate(impulse.voltage, impulse.voltage, "same"))
            plt.show()


def loadPlaneWave(freq_i=100*10**6):
           # testing sample rate in waveform init
    N_samples=4096
    #freq_i=100*10**6 # 100 MHz
    omega=2*numpy.pi*freq_i
    #now omega*t is the argument in the sine wave, so we can get real_times from this
    period=1/freq_i
    print("period: {}".format(period))
    # how many periods do you measure for, a random number maybe?
    rng = numpy.random.default_rng()
    num_periods=1+(rng.random()*3)
    print("plotting {} seconds of wave ".format(period*num_periods))

    # need sampling rate now...
    sampling_period= 3*10**(-10)# lets do 500MHz or 50Mega samples / sec
    times=numpy.linspace(0,N_samples*sampling_period,N_samples)
    #print(times)
    volty=numpy.sin(omega*times)
    real_times=numpy.linspace(0,N_samples*sampling_period,int(100*N_samples))
    #print(real_times[-1])
    real_volty=numpy.sin(omega*real_times)
    reaL_siney=waveform.Waveform(real_volty,real_times*10**9) # what is the sampling rate now? if we had 4096 samples in 2 periods of a 100MHz sine wave then 

    #real_times=numpy.linspace(0,period,N_samples)
    #siney=waveform.Waveform(volty, sampling_rate=(sampling_period)) # what is the sampling rate now? if we had 4096 samples in 2 periods of a 100MHz sine wave then 
    siney=waveform.Waveform(volty, times*10**9) # what is the sampling rate now? if we had 4096 samples in 2 periods of a 100MHz sine wave then 

    #plt.scatter(siney.time,siney.voltage, label='waveform')
    plt.scatter(reaL_siney.time,reaL_siney.voltage, label='real wave')
    plt.scatter(siney.time,siney.voltage, label='sampled wave')

    #plt.scatter(siney.time,siney.voltage, label='sampled wave')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.xlabel('Time [ns]')
    plt.ylabel('Voltage')
    #plt.xlim([0,period*10**9])
    plt.title("waveform for sine wave")
    plt.show()
    return siney

if __name__=="__main__":
    
    import noise

    # for corals, the beamPatterns get more complex, so we need to break apart to multiple calls
    eplane = beamPattern(plot=True,which_plane='E',which_pol='V')
    hplane = beamPattern(plot=True,which_plane='H',which_pol='V')
    #eplane = beamPattern(which_plane='E',which_pol='V')
    #hplane = beamPattern(which_plane='H',which_pol='V')

    #upsampling of my plane wave... checking the original impulse...looks like impulse after upsample should be compared to prepImpulse...
    #should make some plots of the upsampling and prepimpulse before and after? well gimmePlots shouldn't upsample anyways....
    #maybe we remove the upsampling from gimmeplots, right now I am upsampling like way too many times by calling gimmePlots...

    #load a sinusiodal plane wave?
    #impulse=loadPlaneWave(500*10**6) # needs Hz right now

    # on 03/06/24 investigating differences between IR originally put in "TriggerTF_02TH.txt" b/c payload signals look weird with anything but it....

    #load new corals lpda impulse
    #impulse = loadImpulse('impulse/coralsLPDA_impResponse.txt')
    #gimmePlotsImpulse(impulse)
    #impulse2 = loadImpulse('impulse/coralsLPDA_impResponse.txt')
    #gimmePlotsOriginal(impulse)
    #impulse = prepImpulse(impulse)
    #gimmePlots2Impulses(impulse,impulse2)

    # or load impulse event
    impulse = loadImpulse('impulse/triggerTF_02TH.txt')
    gimmePlotsImpulse(impulse)

    #impulse2 = loadImpulse('impulse/triggerTF_02TH.txt')
    impulse = prepImpulse(impulse)
    #impulse.gimmeInfo()
    #gimmePlots(impulse,impulse2)
    #gimmePlotsImpulse(impulse,impulse2)
    
    thermal_noise = noise.ThermalNoise(0.28, .95, filter_order=(10,10), v_rms=1.0, 
                                       fbins=len(impulse.voltage), 
                                       time_domain_sampling_rate=impulse.dt)

    noise = thermal_noise.makeNoiseWaveform(ntraces=aso_geometry.num_antennas)

    '''
    plt.figure()
    plt.plot(impulse.time, impulse.voltage*8+numpy.real(noise[2])[0]) #snr=4

    plt.figure()
    plt.plot(thermal_noise.frequencies, 20*numpy.log10(thermal_noise.amplitudes/numpy.max(thermal_noise.amplitudes)))
    plt.plot(impulse.freq, 20*numpy.log10(numpy.abs(impulse.ampl) / numpy.max(numpy.abs(impulse.ampl))))

    plt.figure()
    sample=50
    #plt.hist(np.real(mynoise[2][:, sample]), bins=np.arange(-6, 6.1, 0.1), normed=True, histtype='step')
    plt.hist(numpy.real(noise[2].flatten()), bins=numpy.arange(-6, 6.1, 0.1), normed=True, histtype='step')
    x = numpy.linspace(-5., 5., 1.e4)
    g = numpy.exp(-0.5 * x**2 * 1**-2) / (numpy.sqrt(2. * numpy.pi))
    plt.plot(x, g, label='Normal Distribution')
    #plt.yscale('log')
    plt.ylim([1.e-2, 1.e0])
    plt.xlabel('V / $\sigma$', size=18)
    plt.ylabel('PDF', size=18)
    plt.tight_layout()
    plt.show()
    '''
    
    #plot a boresight SNR of 5
    #waveform is incoming at phi and el, so maybe better if beam pattern is 2D now?
    #this worked in some way
    #getPayloadWaveforms(22.5, -25, [1,2,3,4], impulse, (eplane, hplane), snr=5, noise=numpy.real(noise[2]), plot=True)
    #try more channels?
    #phi sectors
    #plot needs to make phi sectors use Corals geometry class, not use just ring_map where it assumes there is a top and bottom antenna at every phi...
    trigger_sectors_phi=[1,2,3,4]
    #trigger_sectors_phi=[1,2,3,4,6,7,8]
    print(len(trigger_sectors_phi))

    #getSinglePayloadWaveform(22.5, -25, trigger_sectors_phi, impulse, (eplane, hplane), snr=5, plot=True)

    #getPayloadWaveforms(22.5, -25, trigger_sectors_phi, impulse, (eplane, hplane), snr=5, noise=numpy.real(noise[2]), plot=True)
    #try without noise and just plane wave as impulse
    getPayloadWaveforms(22.5, -25, trigger_sectors_phi, impulse, (eplane, hplane), snr=5, plot=True)
    #getPayloadWaveforms(22.5, -25, trigger_sectors_phi, impulse, (eplane, hplane), snr=5, noise=numpy.real(noise[2]), plot=True)

    #I think I need to include just the plane wave (i.e. no noise added and maybe not even impulse/ impulse response) to see how the delay and all that works for a "trigger_wave" in the getPayload function
    #the trigger_wave here will nto really be anything about a trigger but a check to see if the delays and all that are working correctly...
      
    #getPayloadWaveforms(22.5, -25, [1,2,3,4], impulse, (eplane, hplane), snr=5,  plot=True)


