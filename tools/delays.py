#UWB planewaves->geometric delays
#EJO 11/2016
#   ......cleaned up 12/2017

import myplot  #specific for running on UC midway cluster
import matplotlib.pyplot as plt
import numpy as np
from . constants import *
from . import CoRaLs_geometry as anita

def delay(phi, theta):
    '''
    theta: elevation angle 
    '''
    x_planewave = np.cos(np.radians(theta)) * np.cos(np.radians(phi))
    y_planewave = np.cos(np.radians(theta)) * np.sin(np.radians(phi))
    z_planewave = np.sin(np.radians(theta))
    delays =  -1.0 * ( (anita.x_ant * x_planewave) \
                       + (anita.y_ant * y_planewave) \
                       + (anita.z_ant * z_planewave) ) / c_light
    #print(delays)
    return delays - np.min(delays)

def scanDelays(phi, theta):
    '''
    phi, theta:  numpy arrays of *equal* length
    '''
    t_delays=[]
    for i in range(len(phi)):
        t_delays.append(delay(phi[i], theta[i]))

    return np.array(t_delays)

def getDelays(phi, theta, phi_sectors=range(1,anita.num_phi_sectors), verbose=True):
    '''
    specify phi and theta values (scalars)
    and list of phi sectors of interest (default is all 16)
    print relative delays to terminal if verbose=True
    '''
    t_delay = delay(phi, theta)
    useful = {}
    useful_delays=[]
    for i in range(len(t_delay)):
        for j in phi_sectors:
            if anita.phisector[i] == j:
                useful[str(anita.phisector[i]) + anita.loc[i]] = t_delay[i]
                useful_delays.append(t_delay[i])

    useful_keys= sorted(useful.keys())
    #dump the info in an organized fashion:
    print ('plane wave direction: phi =', phi, 'deg; theta =', theta, 'deg')
    for i in useful_keys:
        #useful[i]-max(useful_delays)
        if verbose:
            #print i, '{0:.1f}'.format(useful[i]-max(useful_delays)), 'ns'
            print (i, '{0:.3f}'.format(useful[i]), 'ns')
    return useful

def getAllDelays(phi, theta, phi_sectors=range(1,anita.num_phi_sectors)):
    '''
    basically same as getDelays, but phi and theta are now numpy arrays,
    and a dictionary of delays is created for each (phi, theta) combination
    '''
    t_delays = scanDelays(phi, theta)
    #print(type(t_delays))
    #print(t_delays.shape)
    data_dict={}
    delays_only=[]
    for k in range(len(t_delays)):
        data_dict[k]={}
        data_dict[k]['phi']=phi[k]
        data_dict[k]['theta']=theta[k]
        data_dict[k]['delays'] = {}

        for i in range(len(t_delays[0])):
            for j in phi_sectors:
                if anita.phisector[i] == j:
                    data_dict[k]['delays'][anita.phisector[i], anita.loc[i]] = t_delays[k][i]
                    #data_dict[k]['delays'][str(anita.phisector[i]) + anita.loc[i]] = t_delays[k][i]

    #[eventually use json to dump to file]
    return data_dict

def plotDelayDictEvent(delay_dict, event):
    '''
    event is an integer
    '''
    if event in delay_dict:
        print ('------------')
        print ('wave theta:', delay_dict[event]['theta'])
        print ('wave phi:  ', delay_dict[event]['phi'])

        phi_sectors=delay_dict[event]['delays'].keys()
        plt.figure()
        for i in phi_sectors:
            if i[1]=='B':
                plt.plot(delay_dict[event]['delays'][i[0],i[1]],int(i[0]), 'o', ms=4, color='blue')
            #elif i[1]=='M':
            #    plt.plot(delay_dict[event]['delays'][i[0],i[1]],int(i[0]), 'o', ms=4, color='green')
            elif i[1]=='T':
                plt.plot(delay_dict[event]['delays'][i[0],i[1]],int(i[0]), 'o', ms=4, color='red')
        plt.grid()
        str_to_plot="Theta: " +str(delay_dict[event]['theta'])+" , Phi: "+str(delay_dict[event]['phi'])
        #print(str_to_plot)
        #plt.text(delay_dict[event]['theta'])
        plt.xlabel('pulse arrival time [ns]')
        plt.ylabel('anita phi sector no.')
        plt.title(str_to_plot)
        plt.tight_layout()
        
    else:
        print ('event specified is not in dataset')

def makeDelayElevationPlot(phi=0.0, phi_sector=1, plot=True):
    '''
    generate relative antenna-pair delays vs incoming plane wave elevation angle 
    for specified phi-sector and plane-wave azimuthal direction
    '''
    thetas = np.arange(-65, 41, 2)
    phis = np.ones(len(thetas)) * phi
    phis2 = np.ones(len(thetas)) * (phi+30.0)
    phis3 = np.ones(len(thetas)) * (phi+330.0)

    t_delays = scanDelays(phis, thetas)
    t_delays2 = scanDelays(phis2, thetas)
    t_delays3 = scanDelays(phis3, thetas)
    phi_sector = int(phi_sector)

    if phi_sector < 1 or phi_sector > anita.num_phi_sectors:
        print ('no such phi-sector')
        return

    if plot:
        fig=plt.figure()
        #plt.plot(thetas, t_delays[:,phi_sector+antennas_per_ring-1] - t_delays[:,phi_sector+antennas_per_ring*2-1], label='mid-bot')
        #plt.plot(thetas, t_delays[:,phi_sector-1] - t_delays[:,phi_sector+antennas_per_ring*2-1], label='top-bot')
        #plt.plot(thetas, t_delays[:,phi_sector-1] - t_delays[:,phi_sector], label='top-bottom')
        plt.plot(thetas, t_delays[:,phi_sector-1] - t_delays[:,phi_sector], label='phi='+str(phi))
        plt.plot(thetas, t_delays2[:,phi_sector-1] - t_delays2[:,phi_sector], label='phi='+str(phi+30.0))
        plt.plot(thetas, t_delays3[:,phi_sector-1] - t_delays3[:,phi_sector], label='phi='+str(phi-30.0))
        plt.legend()
        plt.xlabel('Elevation Angle [deg]')
        plt.ylabel('Delay [ns]')
        plt.title('Delay times: Top-Bottom')
        plt.grid()
        plt.tight_layout()
        #plt.show()
        return fig
    else:
        return np.array((thetas, t_delays[:,phi_sector+antennas_per_ring-1] - t_delays[:,phi_sector+antennas_per_ring*2-1], t_delays[:,phi_sector-1] - t_delays[:,phi_sector+antennas_per_ring*2-1],
                         t_delays[:,phi_sector-1] - t_delays[:,phi_sector+antennas_per_ring-1]))
    
if __name__=='__main__':

    #example usage:
    print(f'anita: {anita.loc!s}')
    print(f'anita: {anita.phisector!s}')
    ## getDelays function:
    phi = 90.0 #11.25
    theta = -50
    #print(delay(phi,theta))
    phi_sectors_of_interest = [1,2,3] #range(1,anita.num_phi_sectors) this will be a top, bottom, and top again..
    getDelays(phi, theta, phi_sectors_of_interest, verbose=True)

    
    ## getAllDelays function:
    phi = np.array([0.0,0.0])
    theta = np.array([-30,0.0])
    phi_sectors_of_interest = [1,2,3,4,5,6,7,8] #range(1,8)
    data_dict=getAllDelays(phi, theta, phi_sectors_of_interest)

    ## read DelayDict (read in output of getAllDelays)
    plotDelayDictEvent(data_dict, 0)
    
    ##make del-el plot
    makeDelayElevationPlot()

    plt.show()
    
