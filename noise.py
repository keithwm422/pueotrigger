import numpy as np
import scipy.stats

class ThermalNoise:
    '''
    units: ns, GHz
    '''
    ############################################################
    def __init__(self, fmin, fmax, fbins=2048, v_rms=1.0, normalize=True,
                 filter_order=(1,1),      #(hipass order, lowpass order)
                 time_domain_sampling_rate=0.1):

        self.vrms = v_rms
        self.n  = fbins
        self.nyq_freq = 1. / (2. * time_domain_sampling_rate)      
        if fmax > self.nyq_freq:
            print ('class ThermalNoise::invalid fmax given input time_domain_sampling_rate')
            return

        self.df= 2. * self.nyq_freq / fbins
        f = np.arange(0, self.nyq_freq + self.df, step=self.df, dtype=float)

        self.frequencies = np.hstack((f, -f[1:len(f)-1][::-1]))
        self.amplitudes = np.zeros(fbins, dtype=float)

        if filter_order == (0,0):
             self.amplitudes[int(np.floor(fbins/2*(fmin/fnyq))):int(np.floor(fbins/2*(fmax/fnyq)))]= 1.0        
        else:
            ##Butterworth type roll-off
            hipass_rolloff = 2 * filter_order[0]
            lopass_rolloff = 2 * filter_order[1]

            for i in range(1,fbins//2+1):
                self.amplitudes[i] = np.sqrt(1./(1+pow((fmin/f[i]),hipass_rolloff)) * \
                                             1./(1+pow((f[i]/fmax),lopass_rolloff)))
        if normalize:
            ##normalize to Vrms
            positive_definite_frequencies = np.ceil(len(self.frequencies)/2)
            n = 0
            for i in range(int(positive_definite_frequencies)):
                n += pow(self.amplitudes[i],2)

            self.amplitudes *= (len(self.amplitudes)+1) * self.vrms * np.sqrt(2) / np.sqrt(n)

    ############################################################
    def makeNoiseWaveform(self, ntraces=1, use_rayleigh=True):
        '''
        make simulated noise waveform in the time-domain
        '''
        passband = np.zeros((ntraces, self.n), dtype=complex)
        ramplitude = np.ones((ntraces, self.n))
        if use_rayleigh:
            ramplitude = scipy.stats.rayleigh.rvs(size=(ntraces, self.n)) / np.sqrt(2)
            
        passband = np.tile(self.amplitudes, ntraces).reshape(ntraces, self.n) * ramplitude *\
                   np.exp(1j * np.random.uniform(0., 2. * np.pi, (ntraces,self.n)))
        noise_waveforms_voltage = np.fft.ifft(passband) 
        noise_waveforms_time    = np.arange(0., self.n/( 2. * self.nyq_freq), 1/( 2. * self.nyq_freq))

        return passband, noise_waveforms_time, noise_waveforms_voltage

    def makeNotNoiseWaveform(self, ntraces=1, use_rayleigh=True):
        '''
        make simulated non-noise waveform in the time-domain to see trends in the simulation
        '''
        passband = np.zeros((ntraces, self.n), dtype=complex)
        ramplitude = np.ones((ntraces, self.n))
        if use_rayleigh:
            ramplitude = scipy.stats.rayleigh.rvs(size=(ntraces, self.n)) / np.sqrt(2)
            
        passband = np.tile(self.amplitudes, ntraces).reshape(ntraces, self.n) * ramplitude *\
                   np.exp(1j * np.random.uniform(0., 2. * np.pi, (ntraces,self.n)))
        noise_waveforms_voltage = np.fft.ifft(passband)
        print("noise waveforms voltage type is: {}".format(type(noise_waveforms_voltage)))
        print("noise waveforms voltage shape is: {}".format(noise_waveforms_voltage.shape))
        #print("noise waveforms voltage {}".format(noise_waveforms_voltage))
        print("noise waveforms voltage last 20 are {}".format(noise_waveforms_voltage[-20:]))
        noise_waveforms_voltage[:,-200:]=(noise_waveforms_voltage[:,-200:]+10)*1.1
        noise_waveforms_time    = np.arange(0., self.n/( 2. * self.nyq_freq), 1/( 2. * self.nyq_freq))

        return passband, noise_waveforms_time, noise_waveforms_voltage


#################################################################
if __name__ == '__main__':
    import myplot
    import matplotlib.pyplot as plt
    import sys
    import payload_signal

    vrms = 1.0
    thermal_noise = ThermalNoise(0.3, 1.2, filter_order=(4,6), v_rms=vrms, fbins=4096, time_domain_sampling_rate=0.2)
    mynoise = thermal_noise.makeNoiseWaveform(ntraces=1)
    #np.save('simulated_noise_2pow26_0p2ns.npy', mynoise[2])
    #sys.exit()

    
    
    '''
    plt.figure()
    plt.plot(thermal_noise.frequencies, 20*np.log10(thermal_noise.amplitudes/np.max(thermal_noise.amplitudes)))
    plt.xlim([0, 2])
    plt.ylim([-20, 2])
    plt.ylabel('Amplitude [dB]', size=16)
    plt.xlabel('Frequency [GHz]', size=16)
    '''

    plt.figure()
    a=np.where(thermal_noise.frequencies > 0.4)[0]
    #plt.hist(np.real(mynoise[0][:, a])/np.cos(np.angle(mynoise[0][:, a])), 50, normed=True, histtype='step')
    plt.hist(np.abs(mynoise[0][:, a[0]])/(np.mean(abs(mynoise[0][:,a[0]]))*np.sqrt(2 / np.pi)), 50, density=True, histtype='step', label='from freq. domain')
    plt.hist(np.abs(np.fft.ifft(mynoise[2])[:,a[0]])/(np.mean(np.abs(np.fft.ifft(mynoise[2])[:,a[0]]))*np.sqrt(2 / np.pi)), 50, density=True, histtype='step', label='inverse FFT from time-domain')

    #print np.mean(abs(mynoise[0][:,a[0]]))*np.sqrt(2 / np.pi)
    #print np.mean(np.abs(np.fft.ifft(mynoise[2])[:,a[0]]))*np.sqrt(2 / np.pi)
    x = np.linspace(0., 5., int(1.e4))
    r = scipy.stats.rayleigh.pdf(x)
    plt.plot(x, r, label='Rayleigh Distribution')
    plt.ylabel('PDF', size=16)
    plt.xlabel('Normed amplitude of bin, {:.1f} MHz < f < {:.1f} MHz'.format(1.e3*thermal_noise.frequencies[a[0]], 1.e3*thermal_noise.frequencies[a[1]]), size=16)
    plt.legend()
    plt.tight_layout()
    
    plt.figure()
    plt.plot(mynoise[1], np.real(mynoise[2][0]), '-', ms=2, color='black')
    plt.ylabel('V / $\sigma$')
    plt.xlabel('Time [ns]')

    print('making voltage histogram...')
    plt.figure()
    sample=50
    #plt.hist(np.real(mynoise[2][:, sample]), bins=np.arange(-6, 6.1, 0.1), normed=True, histtype='step')
    plt.hist(np.real(mynoise[2].flatten()), bins=np.arange(-6, 6.1, 0.1), density=True, histtype='step')
    x = np.linspace(-5., 5., int(1.e4))
    g = np.exp(-0.5 * x**2 * vrms**-2) / (np.sqrt(2. * np.pi) * vrms)
    plt.plot(x, g, label='Normal Distribution')
    plt.yscale('log')
    plt.ylim([1.e-5, 1.e0])
    plt.xlabel('V / $\sigma$', size=18)
    plt.ylabel('PDF', size=18)
    plt.tight_layout()
            
    plt.show()
    
    
    
    
        
        
