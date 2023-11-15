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

def beamPattern(plot=False):
    '''
    anita vpol beam battern
    '''
    angle = [-50,-40,-30,-20,-10,0,10,20,30,40,50]
    eplane_vpol = [-7.5, -5, -3, -1.5, -.5, 0, -.5, -1.5, -3, -5, -7.5]
    hplane_vpol = [-16, -11, -6, -3, -1, 0, -1, -3, -6, -11, -16]
    
    interp_angle = numpy.arange(-50,51,1)
    
    eplane_vpol_interp = interpolate.interp1d(angle, eplane_vpol, kind='cubic')
    hplane_vpol_interp = interpolate.interp1d(angle, hplane_vpol, kind='cubic')
    
    if plot:
        plt.plot(interp_angle, eplane_vpol_interp(interp_angle), label='vpol Eplane')
        plt.plot(interp_angle, hplane_vpol_interp(interp_angle), label='vpol Hplane')
        plt.legend(loc='upper left')
        plt.grid(True)
        plt.xlabel('off-boresight angle [deg.]')
        plt.ylabel('amplitude [dB]')

        plt.show()

    return eplane_vpol_interp, hplane_vpol_interp

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
    trigger_waves=numpy.zeros((len(trigger_sectors), 4, len(impulse.voltage)))

    multiplier=[]
    
    #not yet optimized for speed
    for i in delay[0]['delays']:
        print(i[1])
        print(ring_map[i[1]])
        print(i[0]-numpy.min(trigger_sectors),ring_map[i[1]])
        print(trigger_waves.shape)
        print(noise.shape)
        trigger_waves[i[0]-numpy.min(trigger_sectors),ring_map[i[1]]] = \
            numpy.roll(impulse.voltage * 2 * snr, int(numpy.round(delay[0]['delays'][i] / impulse.dt))) * \
            dBtoVoltsAtten(beam_pattern[1](phi-aso_geometry.phi_ant[i[0]-1])) * \
            dBtoVoltsAtten(beam_pattern[0](el -aso_geometry.theta_ant[0]))

        multiplier.append(dBtoVoltsAtten(beam_pattern[1](phi-aso_geometry.phi_ant[i[0]-1])) * \
            dBtoVoltsAtten(beam_pattern[0](el -aso_geometry.theta_ant[0])))
        
        if noise is not None:
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
                ax[len(ring_map)-j-1,i].set_ylim([-snr-1,snr+1])

        plt.suptitle('phi = '+str(phi)+'deg.  theta = '+str(el)+'deg.', fontsize=20)
        #plt.tight_layout()
        plt.show()
    
    return trigger_waves, impulse.time, multiplier


def downsamplePayload(time, trigger_waves):

    decimate_factor = int(aso_geometry.ritc_sample_step/((time[1]-time[0])))

    trigger_waves = trigger_waves[:,:,::decimate_factor]
    time = time[::decimate_factor]

    return trigger_waves, time                                  
    


if __name__=="__main__":
    
    import noise

    eplane, hplane = beamPattern(plot=True)

    impulse = loadImpulse('impulse/triggerTF_02TH.txt')

    impulse = prepImpulse(impulse)
    
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
    getPayloadWaveforms(22.5, -25, [1,2,3,4], impulse, (eplane, hplane), snr=5, noise=numpy.real(noise[2]), plot=True)
    #getPayloadWaveforms(22.5, -25, [1,2,3,4], impulse, (eplane, hplane), snr=5,  plot=True)


