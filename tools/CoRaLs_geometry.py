import numpy as np

num_phi_sectors = 8
num_skirt_rings = 0 # nothing like this yet in Corals sim
num_top_rings = 0 # probably should be 1 top ring and 1 bottom ring, but the entire idea of "phi" sectors in ANITA is that antennas from multiple rings will have same phi for triggering purposes. 
#num_antennas = num_phi_sectors * (num_top_rings + num_skirt_rings)
num_antennas=8 # for corals, set to 8 antennas

#antenna names / array organization
# Choice of convention, start with a "top" ring antenna as phi sector (but i starts at 0 in for loop below so it looks wrong but it aint)
phisector=[]
loc=[]
for i in range(num_phi_sectors):
    if i % 2 == 0:
        loc.append('T')
    else: 
        loc.append('B')
    phisector.append(i+1)



#azimuthal direction of antennas, degrees
#phi_ant  = np.tile(np.arange(0., 360., 360./num_phi_sectors), 4)

# this is the whole array
# but the existing code needs these positions to be in the order of the declared phi-sectors locations above. Otherwise, user-interfaces with phisectors calls are completely wrong.
# original arrays
#xpos = [ 7.5,    0, -7.5,    0,  5.6, 5.6, -5.6, -5.6]
#ypos = [   0,  7.5,    0, -7.5, -5.6, 5.6,  5.6, -5.6]
#zpos = [-1.3, -1.3, -1.3, -1.3,   -8,  -8,   -8,   -8]
xpos      = [ 7.5,  5.6,    0, -5.6, -7.5, -5.6,    0,  5.6]
ypos      = [   0,  5.6,  7.5,  5.6,    0, -5.6, -7.5, -5.6]
zpos      = [-1.3,   -8, -1.3,   -8, -1.3,   -8, -1.3,   -8] 
#real
phi_tilt   = [   0,   45,   90,  135,  180, -135,  -90,  -45]
theta_tilt = [ -30,  -60,  -30,  -60,  -30,  -60,  -30,  -60]
#make them all the same
#phi_tilt   = [   45,   45,   45,  45,  45, 45,  45,  45]
#theta_tilt = [ -30,  -30,  -30,  -30,  -30,  -30,  -30,  -30]
z_ant = np.array(zpos)
x_ant = np.array(xpos) # Antenna x positions (m)
y_ant = np.array(ypos) # Antenna y positions (m)
r_list = []
i=0
while i<len(x_ant):
    r_list.append(np.sqrt(x_ant[i]**2+y_ant[i]**2+z_ant[i]**2))#in x y plane is just the sqrt x**2+y**2
    i+=1
r_ant=np.array(r_list)
#antenna tilt angle: theta, degrees
#theta_ant=-10.*np.ones(len(z_ant))
# tilt angles from peter's matlab scripts
phi_ant=np.array(phi_tilt)
theta_ant=np.array(theta_tilt)

#normal_az = [0  90 180 -90 -45 45 135 -135];
#normal_el = [-30 -30 -30 -30 -60 -60 -60 -60];


#xpos = [   0, 5.6,  7.5,  5.6,    0,    0, 5.6,  7.5,  5.6,    0]
#ypos = [ 7.5, 5.6,    0, -5.6, -7.5,  7.5, 5.6,    0, -5.6, -7.5]
#zpos = [-1.3,  -8, -1.3,   -8, -1.3, -1.3,  -8, -1.3,   -8, -1.3]
center_x=np.sum(x_ant)/len(x_ant)
center_y=np.sum(y_ant)/len(y_ant)
center_z=np.sum(z_ant)/len(z_ant)
#ritc_sampling
ritc_sample_rate = 2.6 #GHz
ritc_sample_step = 1./ritc_sample_rate #ns

def drawPayload(incoming_wave=False, phi=0, theta=0):
    import myplot    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xpos, ypos, zpos, marker='v', color='gray', alpha=.7, s=120)
    ax.text(center_x, center_y, center_z, '*', fontsize=12)
    l=0
    while l < len(loc):
        ax.text(xpos[l], ypos[l], zpos[l], loc[l], fontsize=12)
        ax.text(xpos[l]+0.5, ypos[l]+0.5, zpos[l]+0.5, str(l+1), fontsize=12)
        l+=1
    if (incoming_wave):
            r=r_ant[0]
            x_planewave = r* np.cos(np.radians(theta)) * np.cos(np.radians(phi))
            y_planewave = r* np.cos(np.radians(theta)) * np.sin(np.radians(phi))
            #z_planewave = (-0.5* r) +r* np.sin(np.radians(theta))
            z_planewave = r* np.sin(np.radians(theta))
            if(theta==0):
                z_planewave=0
            ax.quiver( x_planewave, y_planewave, z_planewave, # <-- starting point of vector
                90-phi, phi, theta, #  directions of vector in degrees
                color = 'red', alpha = .8, lw = 0.3)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    #plt.zlabel('z [m]')

    plt.figure(figsize=(6,8))
    plt.plot(xpos, zpos, 'v', color='gray', ms=30, alpha=.3)
    l=0
    while l < len(loc):
        plt.text(xpos[l], zpos[l], loc[l], fontsize=12)
        plt.text(xpos[l]+0.5, zpos[l]+0.5, str(l+1), fontsize=12)
        l+=1
    plt.xlabel(' x [m]')
    plt.ylabel(' z [m]')
    #plt.ylim([-8, 1])
    plt.show()

def drawPayloadWavefront(incoming_wave=False, phi=0, theta=0):
    import myplot    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xpos, ypos, zpos, marker='v', color='gray', alpha=.7, s=120)
    ax.text(center_x, center_y, center_z, '*', fontsize=12)
    l=0
    while l < len(loc):
        ax.text(xpos[l], ypos[l], zpos[l], loc[l], fontsize=12)
        ax.text(xpos[l]+0.5, ypos[l]+0.5, zpos[l]+0.5, str(l+1), fontsize=12)
        l+=1
        r=1
        x_planewave = r* np.cos(np.radians(theta)) * np.cos(np.radians(phi))
        y_planewave = r* np.cos(np.radians(theta)) * np.sin(np.radians(phi))
        #z_planewave = (-0.5* r) +r* np.sin(np.radians(theta))
        z_planewave = r* np.sin(np.radians(theta))
        if(theta==0):
            z_planewave=0
        u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
        x = np.abs(1-r_ant[0])*x_planewave + np.cos(u) * np.sin(v)
        y = np.abs(1-r_ant[0])*y_planewave + np.sin(u) * np.sin(v)
        z = np.abs(1-r_ant[0])*z_planewave + np.cos(v)
            

    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    #plt.zlabel('z [m]')

    plt.figure(figsize=(6,8))
    plt.plot(xpos, zpos, 'v', color='gray', ms=30, alpha=.3)
    l=0
    while l < len(loc):
        plt.text(xpos[l], zpos[l], loc[l], fontsize=12)
        plt.text(xpos[l]+0.5, zpos[l]+0.5, str(l+1), fontsize=12)
        l+=1
    plt.xlabel(' x [m]')
    plt.ylabel(' z [m]')
    ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)
    print("(x,y,z)=({},{},{})".format(x_planewave,y_planewave,z_planewave))
    ax.plot([0,r_ant[0]*x_planewave],[0,r_ant[0]*y_planewave],zs=[0,r_ant[0]*z_planewave])
    #plt.ylim([-8, 1])
    plt.show()

if __name__=='__main__':
    #drawPayload(incoming_wave=True, phi=30, theta=-80)
    drawPayloadWavefront(incoming_wave=True, phi=30, theta=-80)
