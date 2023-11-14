import numpy as np

num_phi_sectors = 8
num_skirt_rings = 0
num_top_rings = 0
num_antennas = num_phi_sectors * (num_top_rings + num_skirt_rings)

#antenna names / array organization
phisector=[]
loc=[]
for j in ['T', 'B']:
    for i in range(num_phi_sectors):
        phisector.append(i+1)
        loc.append(j)

#azimuthal direction of antennas, degrees
phi_ant  = np.tile(np.arange(0., 360., 360./num_phi_sectors), 4)

# this is the whole array
xpos = [ 7.5,    0, -7.5,    0,  5.6, 5.6, -5.6, -5.6];
ypos = [   0,  7.5,    0, -7.5, -5.6, 5.6,  5.6, -5.6];
zpos = [-1.3, -1.3, -1.3, -1.3,   -8,  -8,   -8,   -8];

z_ant = np.array(zpos)
x_ant = np.array(xpos) # Antenna x positions (m)
y_ant = np.array(ypos) # Antenna y positions (m)
r_list = []
i=0
while i<len(x_ant):
    r_list.append(np.sqrt(x_ant[i]**2+y_ant[i]**2))#in x y plane is just the sqrt x**2+y**2
    i+=1
r_ant=np.array(r_list)
#antenna tilt angle: theta, degrees
theta_ant=-10.*np.ones(len(z_ant))
#xpos = [   0, 5.6,  7.5,  5.6,    0,    0, 5.6,  7.5,  5.6,    0]
#ypos = [ 7.5, 5.6,    0, -5.6, -7.5,  7.5, 5.6,    0, -5.6, -7.5]
#zpos = [-1.3,  -8, -1.3,   -8, -1.3, -1.3,  -8, -1.3,   -8, -1.3]
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
