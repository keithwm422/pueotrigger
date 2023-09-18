import numpy
import myplot
import matplotlib.pyplot as plt
import tools.aso_geometry as aso_geometry
import tools.constants as constants
import payload_signal as payload


if __name__=='__main__':

    phi = 22.5
    theta = -20
    
    delays = payload.getRemappedDelays(phi, theta, [1,2,3,4])

    print (delays)
    print (delays / aso_geometry.ritc_sample_step)

    color = ['blue', 'green', 'red', 'black']
    label = ['bb', 'b', 'm', 't']

    plt.figure(1, figsize = (14,6))
    for i in range(4):
        for j in range(4):

            plt.figure(1)
            if i < 1:
                plt.plot(delays[i,j], i+j*.05, 'o', ms=7, color=color[j], label=label[j])
            else:
                plt.plot(delays[i,j], i+j*.05, 'o', ms=7, color=color[j])


    for i in range(30):
        plt.figure(1)
        plt.plot([aso_geometry.ritc_sample_step*i, aso_geometry.ritc_sample_step*i],
                 [-2,5], '--', c='black', alpha=0.5)


    plt.legend(numpoints=1, loc='lower center')
    plt.ylim([-.8,3.8])
    plt.xlim([-.5,10.5])
    plt.ylabel('Phi Sector')
    plt.xlabel('Delay [ns]')
    plt.title('wave phi = '+str(phi)+'deg.   theta = '+str(theta)+'deg.')

    plt.show()
