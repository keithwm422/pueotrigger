import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c  # m/s
import myplot
from tools import CoRaLs_geometry as corals

c_light = c * 1e-9  # m/ns

def delay(phi, theta):
    '''
    theta: elevation angle 
    '''
    x_planewave = np.cos(np.radians(theta)) * np.cos(np.radians(phi))
    y_planewave = np.cos(np.radians(theta)) * np.sin(np.radians(phi))
    z_planewave = np.sin(np.radians(theta))
    delays =  -1.0 * ( (corals.x_ant * x_planewave) \
                       + (corals.y_ant * y_planewave) \
                       + (corals.z_ant * z_planewave) ) / c_light
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

def getDelays(phi, theta, antennas=None, verbose=True):
    '''
    specify phi and theta values (scalars)
    and list of antennas of interest (default is all)
    print relative delays to terminal if verbose=True
    '''
    if antennas is None:
        antennas = range(corals.num_antennas)
    
    t_delay = delay(phi, theta)
    useful = {}
    useful_delays=[]
    
    for i in antennas:
        ant_name = f'Ant{i+1}'
        useful[ant_name] = t_delay[i]
        useful_delays.append(t_delay[i])

    useful_keys = sorted(useful.keys())
    # dump the info in an organized fashion:
    print('plane wave direction: phi =', phi, 'deg; theta =', theta, 'deg')
    for i in useful_keys:
        if verbose:
            print(i, '{0:.3f}'.format(useful[i]), 'ns')
    return useful

def getAllDelays(phi, theta, antennas=None):
    '''
    basically same as getDelays, but phi and theta are now numpy arrays,
    and a dictionary of delays is created for each (phi, theta) combination
    '''
    if antennas is None:
        antennas = range(corals.num_antennas)
    
    t_delays = scanDelays(phi, theta)
    data_dict = {}
    
    for k in range(len(t_delays)):
        data_dict[k] = {}
        data_dict[k]['phi'] = phi[k]
        data_dict[k]['theta'] = theta[k]
        data_dict[k]['delays'] = {}

        for i in antennas:
            ant_name = f'Ant{i+1}'
            data_dict[k]['delays'][ant_name] = t_delays[k][i]

    # [eventually use json to dump to file]
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
        plt.ylabel('CoRaLS phi sector no.')
        plt.title(str_to_plot)
        plt.tight_layout()
        
    else:
        print ('event specified is not in dataset')
def quantize_delays_ns(delays_ns, decimation_ps=250):
    """
    Quantize delays to hardware step (default 250 ps) and return quantized delays (ns).
    """
    decimation_ns = decimation_ps * 1e-3
    return np.round(np.asarray(delays_ns) / decimation_ns) * decimation_ns

def recommended_angle_steps(phi, theta, antennas=None, sample_dt_ns=0.25, decimation_ps=250,
                            search_max_deg=5.0, search_step_deg=0.01):
    """
    Find the minimum delta-phi and delta-theta (in deg) that cause at least 1-sample
    change in any selected antenna delay after quantization.

    Returns: (min_phi_step_deg, min_theta_step_deg)
    """
    if antennas is None:
        antennas = range(corals.num_antennas)
    
    # Reference quantized delays in samples
    base = delay(phi, theta)
    base_sel = [base[i] for i in antennas]
    base_q = quantize_delays_ns(base_sel, decimation_ps=decimation_ps)
    base_samp = np.asarray(base_q) / sample_dt_ns

    # Scan phi
    min_phi = None
    for dphi in np.arange(search_step_deg, search_max_deg + 1e-9, search_step_deg):
        t = delay(phi + dphi, theta)
        t_sel = [t[i] for i in antennas]
        t_q = quantize_delays_ns(t_sel, decimation_ps=decimation_ps) / sample_dt_ns
        if np.any(np.abs(t_q - base_samp) >= 1):
            min_phi = dphi
            break

    # Scan theta
    min_theta = None
    for dth in np.arange(search_step_deg, search_max_deg + 1e-9, search_step_deg):
        t = delay(phi, theta + dth)
        t_sel = [t[i] for i in antennas]
        t_q = quantize_delays_ns(t_sel, decimation_ps=decimation_ps) / sample_dt_ns
        if np.any(np.abs(t_q - base_samp) >= 1):
            min_theta = dth
            break

    return min_phi, min_theta

if __name__=='__main__':

    print(f'Number of antennas: {corals.num_antennas}')
    print(f'Antenna positions (x,y,z):')
    for i in range(corals.num_antennas):
        print(f'  Ant{i+1}: ({corals.xpos[i]:.2f}, {corals.ypos[i]:.2f}, {corals.zpos[i]:.2f}) m')
    
    phi = 100
    theta = 10
    antennas_of_interest = list(range(corals.num_antennas))  # all channels
    getDelays(phi, theta, antennas_of_interest, verbose=True)

    sample_dt_ns = corals.ritc_sample_step
    # recommended_angle_steps only needs the 4 physical positions (ch0-3; ch4-7 are co-located)
    min_phi_step, min_theta_step = recommended_angle_steps(
        phi, theta, [0, 1, 2, 3], sample_dt_ns=sample_dt_ns, decimation_ps=250)
    print(f"Minimum phi step for 1 sample change: {min_phi_step} deg")
    print(f"Minimum theta step for 1 sample change: {min_theta_step} deg")