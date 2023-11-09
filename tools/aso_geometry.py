import numpy as np

num_phi_sectors = 24
num_skirt_rings = 3
num_top_rings = 1
num_antennas = num_phi_sectors * (num_top_rings + num_skirt_rings)

#antenna names / array organization
phisector=[]
loc=[]
for j in ['T', 'M', 'B', 'BB']:
    for i in range(num_phi_sectors):
        phisector.append(i+1)
        loc.append(j)

#azimuthal direction of antennas, degrees
phi_ant  = np.tile(np.arange(0., 360., 360./num_phi_sectors), 4)

#radial position of antennas, meters
radius_top_even = 1.036
radius_top_odd = 0.926
radius_skirt = 2.234

r_ant = np.concatenate(( np.tile(np.array([radius_top_even, radius_top_odd]), num_top_rings * int(num_phi_sectors/2)),
                         np.tile(np.array([radius_skirt]), num_skirt_rings * num_phi_sectors) ))

#vertical position of antennas, meters
z_bb = -5.09
z_b = -4.361
z_m = -3.632
z_top_odd = -0.63 
z_top_even = -0.002 

z_ant = np.concatenate(( np.tile(np.array([z_top_even, z_top_odd]), num_top_rings * int(num_phi_sectors/2)),
                         np.tile(np.array([z_m]), num_phi_sectors),
                         np.tile(np.array([z_b]), num_phi_sectors),
                         np.tile(np.array([z_bb]), num_phi_sectors) ))


x_ant = r_ant * np.cos(np.radians(phi_ant)) # Antenna x positions (m)
y_ant = r_ant * np.sin(np.radians(phi_ant)) # Antenna y positions (m)

#antenna tilt angle: theta, degrees
theta_ant=-10.*np.ones(len(z_ant))
#xpos = [0,5.6,7.5,5.6,0,0,5.6,7.5,5.6,0]
#ypos = [7.5,5.6,0,-5.6,-7.5,7.5,5.6,0,-5.6,-7.5]
#zpos = [-1.3,-8,-1.3,-8,-1.3,-1.3,-8,-1.3,-8,-1.3]
# this is the whole array
xpos = [ 7.5,    0, -7.5,    0,  5.6, 5.6, -5.6, -5.6];
ypos = [   0,  7.5,    0, -7.5, -5.6, 5.6,  5.6, -5.6];
zpos = [-1.3, -1.3, -1.3, -1.3,   -8,  -8,   -8,   -8];

#ritc_sampling
ritc_sample_rate = 2.6 #GHz
ritc_sample_step = 1./ritc_sample_rate #ns

def drawPayload():
    import myplot    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xpos, ypos, zpos, marker='s', color='gray', alpha=.7, s=120)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    #plt.zlabel('z [m]')

    plt.figure(figsize=(6,8))
    plt.plot(xpos, zpos, 's', color='gray', ms=30, alpha=.3)
    plt.xlabel(' x [m]')
    plt.ylabel(' z [m]')
    #plt.ylim([-8, 1])
    plt.show()

if __name__=='__main__':
    drawPayload()
