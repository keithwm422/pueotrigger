import numpy as np

# 8 channels: 4 physical antennas × 2 pols (H,V)
# Channels 0-3: H-pol, Channels 4-7: V-pol
num_antennas = 8

ritc_sample_rate = 4  # GHz
ritc_sample_step = 1/ritc_sample_rate  # ns

# Antenna positions in x-z plane (y=0)
span = 2
d = span/(2*np.sqrt(2))  # d = 1/√2 ≈ 0.707 m
xpos = [0, 0,  0,  0,  0, 0,  0,  0]
ypos = [d, d, -d, -d,  d, d, -d, -d]
zpos = [d, -d, -d, d,  d, -d, -d, d]

phi_tilt   = [0, 0, 0, 0,   # H-pol boresight pointing (phi)
              0, 0, 0, 0]   # V-pol
theta_tilt = [0, 0, 0, 0,   # H-pol boresight pointing (theta)
              0, 0, 0, 0]   # V-pol

x_ant = np.array(xpos)
y_ant = np.array(ypos)
z_ant = np.array(zpos)
r_ant = np.sqrt(x_ant**2 + y_ant**2 + z_ant**2)
phi_ant = np.array(phi_tilt)
theta_ant = np.array(theta_tilt)

center_x = np.mean(x_ant[:4])  # Use first 4 positions for center
center_y = np.mean(y_ant[:4])
center_z = np.mean(z_ant[:4])

def drawPayload(incoming_wave=False, phi=0, theta=0):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(14, 7))
    
    # 3D view
    ax = fig.add_subplot(121, projection='3d')
    
    ax.scatter(xpos[:4], ypos[:4], zpos[:4], marker='^', color='blue', s=300, 
               alpha=0.8, edgecolors='black', linewidth=2, label='H-pol')

    ax.scatter(xpos[4:], ypos[4:], zpos[4:], marker='s', color='red', s=200, 
               alpha=0.6, edgecolors='black', linewidth=2, label='V-pol')
    
    arrow_length = 0.5
    for i in range(4):
        theta_rad = np.radians(theta_ant[i])
        phi_rad = np.radians(phi_ant[i])
        dx = np.cos(theta_rad) * np.cos(phi_rad)
        dy = np.cos(theta_rad) * np.sin(phi_rad)
        dz = np.sin(theta_rad)
        ax.quiver(xpos[i], ypos[i], zpos[i], dx*arrow_length, dy*arrow_length, dz*arrow_length,
                 color='green', arrow_length_ratio=0.3, linewidth=2)
        ax.text(xpos[i], ypos[i], zpos[i]+0.3, f'Ch{i}(H)\nCh{i+4}(V)', 
                fontsize=10, fontweight='bold', ha='center')
    
    ax.scatter([center_x], [center_y], [center_z], marker='*', color='red', s=400, label='Center')
    
    if incoming_wave:
        r = 2.5
        x_planewave = r * np.cos(np.radians(theta)) * np.cos(np.radians(phi))
        y_planewave = r * np.cos(np.radians(theta)) * np.sin(np.radians(phi))
        z_planewave = r * np.sin(np.radians(theta))
        
        ax.quiver(0, 0, 1, x_planewave, y_planewave, z_planewave-1,
                 color='red', arrow_length_ratio=0.2, linewidth=3, alpha=0.7, label=f'Wave φ={phi}°, θ={theta}°')
    
    ax.set_xlabel('X [m]', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y [m]', fontsize=12, fontweight='bold')
    ax.set_zlabel('Z [m]', fontsize=12, fontweight='bold')
    ax.set_title('3D Antenna Array', fontsize=14, fontweight='bold')
    
    max_range = 2.5
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-1, max_range])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2D side view (x-z plane)
    ax2 = fig.add_subplot(122)
    ax2.scatter(xpos, zpos, marker='^', color='blue', s=300, alpha=0.8, edgecolors='black', linewidth=2)
    
    for i in range(num_antennas):
        theta_rad = np.radians(theta_ant[i])
        phi_rad = np.radians(phi_ant[i])
        dx = np.cos(theta_rad) * np.cos(phi_rad)
        dz = np.sin(theta_rad)
        ax2.arrow(xpos[i], zpos[i], dx*arrow_length, dz*arrow_length,
                 head_width=0.2, head_length=0.1, fc='green', ec='green', linewidth=2)
        ax2.text(xpos[i], zpos[i]+0.15, f'A{i+1}', fontsize=12, fontweight='bold', ha='center')
    
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=1)
    ax2.axvline(x=0, color='k', linestyle='--', alpha=0.3, linewidth=1)
    ax2.set_xlabel('X [m]', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Z [m]', fontsize=12, fontweight='bold')
    ax2.set_title('Side View (X-Z Plane)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    ax2.set_ylim([-1, 0.5])
    
    plt.tight_layout()
    plt.show()



def drawWavefrontPlanes(incoming_wave=False, phi=0, theta=0, show_labels=True):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    # Calculate wavefront direction
    r = 1
    x_planewave = r * np.cos(np.radians(theta)) * np.cos(np.radians(phi))
    y_planewave = r * np.cos(np.radians(theta)) * np.sin(np.radians(phi))
    z_planewave = r * np.sin(np.radians(theta))
    
    # Create figure with 2x2 layout
    fig = plt.figure(figsize=(14, 12))
    fig.suptitle(f'CoRaLS Antenna Array Geometry - Orthographic Projection\nWave: φ={phi}°, θ={theta}°', 
                 fontsize=16, fontweight='bold')
    
    # Top view (x-y plane) - upper left
    ax1 = fig.add_subplot(2, 2, 1)
    for idx in range(num_antennas):
        ax1.scatter(xpos[idx], ypos[idx], marker='^', color='blue', s=300, 
                   alpha=0.8, edgecolors='black', linewidth=2)
        if show_labels:
            ax1.text(xpos[idx], ypos[idx]+0.25, f'A{idx+1}', fontsize=12, 
                    fontweight='bold', ha='center')
        # Draw pointing arrow (x-y projection)
        theta_rad = np.radians(theta_ant[idx])
        phi_rad   = np.radians(phi_ant[idx])
        dx = np.cos(theta_rad) * np.cos(phi_rad)
        dy = np.cos(theta_rad) * np.sin(phi_rad)
        proj_mag = np.sqrt(dx**2 + dy**2)
        if proj_mag > 1e-6:
            ax1.arrow(xpos[idx], ypos[idx], dx*0.5, dy*0.5,
                     head_width=0.18, head_length=0.1, fc='green', ec='green',
                     linewidth=2, length_includes_head=True)
    
    # Draw square connecting antennas
    square_x = [xpos[0], xpos[1], xpos[2], xpos[3], xpos[0]]
    square_y = [ypos[0], ypos[1], ypos[2], ypos[3], ypos[0]]
    ax1.plot(square_x, square_y, 'k--', alpha=0.3, linewidth=1)
    
    # Draw incoming wave direction (top view projection)
    if incoming_wave:
        scale = 1.5
        ax1.arrow(0, 0, scale*x_planewave, scale*y_planewave, 
                 head_width=0.25, head_length=0.15, fc='red', ec='red', 
                 alpha=0.7, linewidth=2.5, label='Wave direction')
    
    ax1.scatter([0], [0], marker='*', color='red', s=300, zorder=10)
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.2, linewidth=0.5)
    ax1.axvline(x=0, color='k', linestyle='--', alpha=0.2, linewidth=0.5)
    ax1.set_xlabel('X [m]', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Y [m]', fontsize=12, fontweight='bold')
    ax1.set_title('Top View (X-Y Plane)', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    ax1.legend(loc='upper right')
    ax1.set_xlim([-2.5, 2.5])
    ax1.set_ylim([-2.5, 2.5])
    
    # Front view (x-z plane) - lower left
    ax2 = fig.add_subplot(2, 2, 3)
    # Group antennas by x position to handle overlaps
    x_groups = {}
    for idx in range(num_antennas):
        x_key = round(xpos[idx], 2)
        if x_key not in x_groups:
            x_groups[x_key] = []
        x_groups[x_key].append(idx)
    
    for idx in range(num_antennas):
        ax2.scatter(xpos[idx], zpos[idx], marker='^', color='blue', s=300,
                   alpha=0.8, edgecolors='black', linewidth=2)
        # Draw antenna pointing arrows (x-z projection)
        theta_rad = np.radians(theta_ant[idx])
        phi_rad = np.radians(phi_ant[idx])
        dx = np.cos(theta_rad) * np.cos(phi_rad)
        dz = np.sin(theta_rad)
        ax2.arrow(xpos[idx], zpos[idx], dx*0.5, dz*0.5,
                 head_width=0.2, head_length=0.08, fc='green', ec='green', linewidth=2)
        
        # Offset labels for antennas at same x position
        if show_labels:
            x_key = round(xpos[idx], 2)
            if len(x_groups[x_key]) > 1:
                offset_idx = x_groups[x_key].index(idx)
                x_offset = -0.3 + offset_idx * 0.6  # Spread labels left and right
                ax2.text(xpos[idx] + x_offset, zpos[idx]+0.12, f'A{idx+1}', fontsize=12, 
                        fontweight='bold', ha='center')
            else:
                ax2.text(xpos[idx], zpos[idx]+0.12, f'A{idx+1}', fontsize=12, 
                        fontweight='bold', ha='center')
    
    # Draw incoming wave direction (X-Z projection)
    if incoming_wave:
        # Normalize the X-Z projection to fit in plot
        xz_mag = np.sqrt(x_planewave**2 + z_planewave**2)
        if xz_mag > 0:
            scale = 0.8  # Scale to fit nicely in the view
            arrow_x = scale * (x_planewave / xz_mag)
            arrow_z = scale * (z_planewave / xz_mag)
            ax2.arrow(0, 0.3, arrow_x, arrow_z, 
                     head_width=0.2, head_length=0.12, fc='red', ec='red', 
                     alpha=0.7, linewidth=2.5, label='Wave direction')
            ax2.legend(loc='upper right', fontsize=10)
    
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.4, linewidth=1)
    ax2.axvline(x=0, color='k', linestyle='--', alpha=0.2, linewidth=0.5)
    ax2.set_xlabel('X [m]', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Z [m]', fontsize=12, fontweight='bold')
    ax2.set_title('FRONT VIEW (X-Z Plane)', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_xlim([-2.5, 2.5])
    ax2.set_ylim([-2.5, 2.5])
    
    # Side view (y-z plane) - upper right  
    ax3 = fig.add_subplot(2, 2, 2)
    # Group antennas by y position to handle overlaps
    y_groups = {}
    for idx in range(num_antennas):
        y_key = round(ypos[idx], 2)
        if y_key not in y_groups:
            y_groups[y_key] = []
        y_groups[y_key].append(idx)
    
    for idx in range(num_antennas):
        ax3.scatter(ypos[idx], zpos[idx], marker='^', color='blue', s=300,
                   alpha=0.8, edgecolors='black', linewidth=2)
        # Draw antenna pointing arrows (y-z projection)
        theta_rad = np.radians(theta_ant[idx])
        phi_rad = np.radians(phi_ant[idx])
        dy = np.cos(theta_rad) * np.sin(phi_rad)
        dz = np.sin(theta_rad)
        proj_mag = np.sqrt(dy**2 + dz**2)
        if proj_mag > 1e-6:
            ax3.arrow(ypos[idx], zpos[idx], dy*0.5, dz*0.5,
                     head_width=0.2, head_length=0.08, fc='green', ec='green', linewidth=2)
        else:
            # Pointing direction is perpendicular to this plane (pure +X when phi=0, theta=0)
            ax3.annotate('→X', xy=(ypos[idx], zpos[idx]), fontsize=9, color='green',
                        ha='center', va='center', fontweight='bold')
        
        # Offset labels for antennas at same y position
        if show_labels:
            y_key = round(ypos[idx], 2)
            if len(y_groups[y_key]) > 1:
                offset_idx = y_groups[y_key].index(idx)
                y_offset = -0.3 + offset_idx * 0.6  # Spread labels left and right
                ax3.text(ypos[idx] + y_offset, zpos[idx]+0.12, f'A{idx+1}', fontsize=12, 
                        fontweight='bold', ha='center')
            else:
                ax3.text(ypos[idx], zpos[idx]+0.12, f'A{idx+1}', fontsize=12, 
                        fontweight='bold', ha='center')
    
    # Draw incoming wave direction (Y-Z projection)
    if incoming_wave:
        # Normalize the Y-Z projection to fit in plot
        yz_mag = np.sqrt(y_planewave**2 + z_planewave**2)
        if yz_mag > 0:
            scale = 0.8  # Scale to fit nicely in the view
            arrow_y = scale * (y_planewave / yz_mag)
            arrow_z = scale * (z_planewave / yz_mag)
            ax3.arrow(0, 0.3, arrow_y, arrow_z, 
                     head_width=0.2, head_length=0.12, fc='red', ec='red', 
                     alpha=0.7, linewidth=2.5, label='Wave direction')
            ax3.legend(loc='upper right', fontsize=10)
    
    ax3.axhline(y=0, color='k', linestyle='-', alpha=0.4, linewidth=1)
    ax3.axvline(x=0, color='k', linestyle='--', alpha=0.2, linewidth=0.5)
    ax3.set_xlabel('Y [m]', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Z [m]', fontsize=12, fontweight='bold')
    ax3.set_title('SIDE VIEW (Y-Z Plane)', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_aspect('equal', adjustable='box')
    ax3.set_xlim([-2.5, 2.5])
    ax3.set_ylim([-2.5, 2.5])
    
    plt.tight_layout()
    plt.savefig('plots/geometry.png', dpi=150, bbox_inches='tight')
    print(f"\nWavefront direction (unit vector): ({x_planewave:.3f}, {y_planewave:.3f}, {z_planewave:.3f})")
    print(f"Plot saved to: plots/geometry.png")
    plt.show()

if __name__=='__main__':
    drawWavefrontPlanes(incoming_wave=True, phi=30, theta=20, show_labels= False)
    drawPayload(incoming_wave=False, phi=90, theta=0)
