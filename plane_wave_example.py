import numpy
import tools.waveform as waveform
import tools.delays as delays
import tools.CoRaLs_geometry as aso_geometry
from scipy import interpolate
from scipy.signal import lfilter, butter, cheby1
import matplotlib.pyplot as plt

def loadImpulse(filename='impulse/triggerTF_02TH.txt'):

    dat=numpy.loadtxt(filename)
    impulse=waveform.Waveform(dat[:,1], time=dat[:,0])
    return impulse

def gimmePlots(impulse):
            impulse.fft()
            plt.figure(1)
            plt.plot(impulse.time_zeroed, impulse.voltage, 's', ms=2)
            plt.figure(2)
            plt.plot(impulse.freq, numpy.abs(impulse.ampl), 's', ms=2)
            impulse.upsampleFreqDomain(2)
            plt.figure(1)
            plt.plot(impulse.time, impulse.voltage, 'o', ms=1)
            plt.figure(2)
            plt.plot(impulse.freq, numpy.abs(impulse.ampl), 'o', ms=1)    
            #plt.figure(4)
            #plt.plot(impulse.time, numpy.correlate(impulse.voltage, impulse.voltage, "same"))
            plt.show()

if __name__=="__main__":
    
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
    sampling_period= 3*10**(-10)# lets do 500MHz or 50Mega samples / sec
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
    plt.scatter(reaL_siney.time*10**9,reaL_siney.voltage, label='real wave')
    plt.scatter(times*10**9,volty, label='sampled wave')

    #plt.scatter(siney.time,siney.voltage, label='sampled wave')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.xlabel('Time [ns]')
    plt.ylabel('Voltage')
    plt.xlim([0,period*10**9])
    plt.show()
    #quit()
    #impulse = loadImpulse()
    # or load impulse event
    impulse = loadImpulse('impulse/triggerTF_02TH.txt')
    
    gimmePlots(impulse)
    #impulse = prepImpulse(impulse)