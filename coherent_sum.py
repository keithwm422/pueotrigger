import numpy
import myplot
import matplotlib.pyplot as plt
import tools.aso_geometry as aso_geometry
import tools.constants as constants
import payload_signal as payload
import math

def coherentSum(waveforms, timebase, delays, downsample=False, ringmask=[1,1,1,1]):
    '''
    mask = ring mask [bb, b, m, t]
    '''
    
    #decimates to specified sample rate
    decimate_factor = int(aso_geometry.ritc_sample_step/((timebase[1]-timebase[0])))

    #downsample to ritc sampling, if input is upsampled:
    if downsample == True:
        coh_sum = numpy.zeros(int(len(waveforms[0,0]) / decimate_factor)) # choice to cast as int, so 333.333 becomes 333
        
        for i in range(waveforms.shape[0]):
            for j in range(waveforms.shape[1]):
                if ringmask[j]:
                    _wave = waveforms[i,j][::decimate_factor]
                    _delay = -int(numpy.round(delays[i,j]/aso_geometry.ritc_sample_step))
                    coh_sum = coh_sum + numpy.roll(_wave[:len(coh_sum)], _delay)

        timebase = timebase[::decimate_factor]
        timebase = timebase[:len(coh_sum)]

    #use input sampling:
    else:
        coh_sum = numpy.zeros(len(waveforms[0,0]))
        for i in range(waveforms.shape[0]):
            for j in range(waveforms.shape[1]):
                if ringmask[j]:
                    _wave = waveforms[i,j]
                    _delay = -int(numpy.round(delays[i,j]/(timebase[1]-timebase[0])))
                    coh_sum = coh_sum + numpy.roll(_wave, _delay)

        coh_sum = coh_sum[::decimate_factor]
        timebase = timebase[::decimate_factor]

    return coh_sum, timebase

def powerSum(coh_sum, window=32, step=16):
    '''
    calculate power summed over a length defined by 'window', overlapping at intervals defined by 'step'
    '''
    num_frames = int(math.floor((len(coh_sum)-window) / step))
   # print(num_frames)
    #print(window)
    coh_sum_squared = (coh_sum * coh_sum).astype(numpy.int64)
   # print(coh_sum_squared)
    #print(coh_sum_squared.strides[0]*step)
   #print(coh_sum_squared.strides[0])
    coh_sum_windowed = numpy.lib.stride_tricks.as_strided(coh_sum_squared, (num_frames, window),
                                                          (int(coh_sum_squared.strides[0]*step), coh_sum_squared.strides[0]))
    power = numpy.sum(coh_sum_windowed, axis=1)

    return power.astype(numpy.float64)/window, num_frames

                
if __name__=='__main__':

    phi = 22.5
    theta = -10

    eplane, hplane = payload.beamPattern(plot=False)
    
    impulse = payload.loadImpulse('impulse/triggerTF_02TH.txt')
    impulse = payload.prepImpulse(impulse)

    print ('timestep of input pulse [ns]:', impulse.dt)
    
    delays = payload.getRemappedDelays(phi, theta, [1,2,3,4])
    waveforms, timebase,_ = payload.getPayloadWaveforms(phi, theta, [1,2,3,4], impulse, (eplane, hplane))

    coh_sum, timebase_coh_sum  = coherentSum(waveforms, timebase, delays, True)
    #coh_sum2, timebase_coh_sum2  = coherentSum(waveforms, timebase, delays, False)

    #print len(coh_sum2), len(timebase_coh_sum2)

    fig, ax = plt.subplots(2,1, sharex=True)
    #ax[0].plot(timebase_coh_sum, coh_sum,  ':o', c='black',label='ideal', ms=2)
    #ax[0].plot(timebase_coh_sum2, coh_sum2,  '-o', c='red',label='ritc sample step', alpha=0.6, ms=4)
    ax[0].plot(timebase_coh_sum, coh_sum,  ':o', c='black', ms=2)
    ax[0].legend()
    ax[0].set_ylabel('16-ant coherent sum')

    power, frames = powerSum(coh_sum)
    #power2, frames2 = powerSum(coh_sum2)
    
    #ax[1].plot(timebase_coh_sum[:frames*16:16]+timebase_coh_sum[8], power, '--', c='black', label='ideal')
    #ax[1].plot(timebase_coh_sum2[:frames*16:16]+timebase_coh_sum2[8], power2, '-', c='red', label='ritc sample step', alpha=0.6)
    ax[1].plot(timebase_coh_sum[:frames*16:16]+timebase_coh_sum[8], power, '--', c='black')
    ax[1].legend()
    ax[1].set_xlabel('Time [ns]')
    ax[1].set_ylabel('Power [arb]')
    ax[1].set_ylim([-1, 31])

    #scan phi, 45 degrees:
    phi_scan_power_true = []
    phi_scan_power_ritc = []

    phiscan_start = 0
    elscan_start = -30
    
    for phiscan in range(phiscan_start, 45, 1):
        waveforms, timebase,_ = payload.getPayloadWaveforms(phiscan, theta, [1,2,3,4], impulse, (eplane, hplane))
        coh_sum, timebase_coh_sum  = coherentSum(waveforms, timebase, delays, False, ringmask = [1,1,1,1])
        coh_sum2, timebase_coh_sum2  = coherentSum(waveforms, timebase, delays, False, ringmask = [1,1,1,0])
        power, frames = powerSum(coh_sum)
        power2, frames2  = powerSum(coh_sum2)
        phi_scan_power_true.append(numpy.max(power))
        phi_scan_power_ritc.append(numpy.max(power2))

    el_scan_power_true = []
    el_scan_power_ritc = []

    for elscan in range(elscan_start,10,1):
        waveforms, timebase,_ = payload.getPayloadWaveforms(phi, elscan, [1,2,3,4], impulse, (eplane, hplane))
        coh_sum, timebase_coh_sum = coherentSum(waveforms, timebase, delays, False, ringmask = [1,1,1,1])
        coh_sum2, timebase_coh_sum2  = coherentSum(waveforms, timebase, delays, False, ringmask = [1,1,1,0])
        power, frames = powerSum(coh_sum)
        power2, frames2 = powerSum(coh_sum2)

        el_scan_power_true.append(numpy.max(power))
        el_scan_power_ritc.append(numpy.max(power2))    
    

    fig, ax = plt.subplots(2,1)
    ax[0].plot(numpy.array(range(45),dtype=float)-phiscan_start-22.5, 10*numpy.log10(phi_scan_power_true/numpy.max(phi_scan_power_true)),'--', c='black', label='4-ring')
    ax[0].plot(numpy.array(range(45),dtype=float)-phiscan_start-22.5, 10*numpy.log10(phi_scan_power_ritc/numpy.max(phi_scan_power_ritc)), '-', c='red', label='3-ring, no top', alpha=0.6)
    ax[0].set_ylabel('Power (normalized) [dB]')
    ax[0].set_xlabel('Phi [degrees]')
    ax[0].set_ylim([-10, 1])
    ax[0].legend(loc='upper left')
    ax[0].grid(True)

    ax[1].plot(numpy.array(range(40),dtype=float)+elscan_start+10, 10*numpy.log10(el_scan_power_true/numpy.max(el_scan_power_true)),'--', c='black', label='4-ring')
    ax[1].plot(numpy.array(range(40),dtype=float)+elscan_start+10, 10*numpy.log10(el_scan_power_ritc/numpy.max(el_scan_power_ritc)), '-', c='red', label='3-ring, no top', alpha=0.6)
    ax[1].set_ylabel('Power (normalized) [dB]')
    ax[1].set_xlabel('theta [degrees]')
    ax[1].set_ylim([-10, 1])
    ax[1].legend(loc='upper left')  
    ax[1].grid(True)
    plt.show()
                                              
 
