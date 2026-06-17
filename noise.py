import numpy as np
import scipy.stats

class ThermalNoise:
    '''
    units: ns, GHz
    '''
    ############################################################
    def __init__(self, fbins=2048, v_rms=1.0, normalize=True,     
                 time_domain_sampling_rate=0.1):
        '''
        Generate white thermal noise in frequency domain.
        Noise is generated as white noise.
        '''
        self.vrms = v_rms
        self.n  = fbins
        self.nyq_freq = 1. / (2. * time_domain_sampling_rate)      

        self.df= 2. * self.nyq_freq / fbins
        f = np.arange(0, self.nyq_freq + self.df, step=self.df, dtype=float)

        self.frequencies = np.hstack((f, -f[1:len(f)-1][::-1]))
        self.amplitudes = np.zeros(fbins, dtype=float)
        self.amplitudes[:fbins//2+1] = 1.0
        
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