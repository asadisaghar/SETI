import time
import numpy as np
from math import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.rc("font", size=14, family='serif', weight='normal')
#plt.rc("axes", labelsize=12, titlesize=14)
#plt.rc("xtick", labelsize=12)
#plt.rc("ytick", labelsize=12)
#plt.rc("legend", fontsize=12)
#from scipy.optimize import curve_fit
import random
import tools 
start_time = time.time()
import cProfile

#        CONSTANTS     ##
# ===================== #
G = 4.302e-6            # kpc.M_solar^-1.(km/s)^2 
c_pcyr = 0.306594845    # pc/yr
yr = pi*1.e7            # s
c_pcs = c_pcyr/yr       # pc/s
c = c_pcs
pc = 3.0857e13          # km
Myr = pi*1.e7*1.e6      # s
k = 1000                # kilo!
# ===================== #

#     HANDLES      #
# ================ #
N_disk = int(1.e4)  # Number of particles in disk
#dt = np.logspace(0, 2, num=1)*Myr
dt = 1*Myr
dt_const = np.logspace(-6, 1, num=1)*Myr
VC = np.logspace(-1, 0, num=1)*c
t = 0.  # time (s)
t_Myr = t/Myr  # time (Myr)
#dt = dt_const # time step (s)
dt_Myr = dt/Myr  # time step(Myr)
dist = VC * (dt-dt_const)    # pc
t_f = 1.e1*Myr

# Galactic parameters  #
# =====================#
# DISK
N_thin = int(0.90*N_disk)
N_thick = int(N_disk-N_thin)
h_z_thin = 3.e2  # scale height of the thin disk(pc)(MW:300-400pc)
h_z_thick = 5.e2  # scale height of the thick disk(pc)(MW:1-1.5kpc)
h_rho_disk = 5.e3  # scale length of the thin disk(pc)(MW:3kpc)
I_0_disk = 20.e9  # disk central intensity
alpha_disk = h_rho_disk  # disk scale length(pc)
n_s_disk = 1  # disk Sersic index
#M_disk = 6.e10  # disk mass(M_solar)
#VH = 300. #(km/s)
VH = 300./pc #(pc/s)
aH = 5.e3 #(pc)   

# BULGE
N_bulge = int(0.33000*N_disk)
R_bulge = 2.7e3  # bulge radius(pc)(MW:a few kpc)
I_0_bulge = 5.e9  # bulge central intensity
alpha_bulge = R_bulge/3.  # bulge scale length(pc)
n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
mean_bulge = 200 # (My arbitrary value!!)(km/s) 
sigma_bulge = 130 # Velocity dispercion(km/s)
#M_bulge = 2.e10  # bulge mass(M_solar)

#GALAXY
alpha_gal = alpha_disk
N_gal = N_disk + N_bulge
#M_DM = 3.e12  # dark halo mass(M_solar)
#M_gal = M_disk + M_bulge + M_DM
t_f = 1.e2*Myr  # Final time(s), i.e. 1Gyr
R_opt = 2.*h_rho_disk # (kpc)
V_opt = 300. # (km/s)
n_s_gal = 4

#    Initialisation (t=0)    #
# ========================== #
# POS:Thin disk (Cylindrical)
rho_thin= tools.init_pos(N_thin, 0., h_rho_disk, 'exp')
# disk position initializing (r)
phi_thin = tools.init_pos(N_thin, 0., 2*pi, 'uni')
# disk position initializing (phi)
z_thin = tools.init_z(N_thin, 0., h_z_thin, 'exp')
# thin disk position initializing (z)
# POS:Thick disk (Cylindrical)
rho_thick = tools.init_pos(N_thick, 0, h_rho_disk, 'exp')
# disk position initializing (r)
phi_thick = tools.init_pos(N_thick, 0., 2*pi, 'uni')
# disk position initializing (phi)
z_thick = tools.init_z(N_thick, 0., h_z_thick, 'exp')
# thin disk position initializing (z)
# POS:Disk (Cylindrical)
rho_disk = np.append(rho_thin, rho_thick)
phi_disk = np.append(phi_thin, phi_thick)
z_disk = np.append(z_thin, z_thick)
# POS:Bulge (Spherical)
r_bulge = tools.init_pos(N_bulge, 0., R_bulge, 'uni')
phi_bulge = tools.init_pos(N_bulge, 0., pi, 'uni')
theta_bulge = tools.init_pos(N_bulge, 0., 2.*pi, 'uni')
# POS:Bulge (Cylindrical)
rho_bulge = r_bulge*np.sin(phi_bulge)
z_bulge = r_bulge*np.cos(phi_bulge)
phi_bulge = theta_bulge
# VEL:Disk (Rotation curve analytic relation)
V_disk = tools.v_r_curve(rho_disk, VH, aH)
# VEL:Bulge (Gaussian dist. with given mean and dispersion)
V_bulge = tools.init_norm(mean_bulge, sigma_bulge, N_bulge)
# I(I-band):Disk
I_disk = tools.init_sersic(I_0_disk, alpha_disk, n_s_disk, rho_disk)
# L(I-band):Bulge (Gaussian dist. with given total L)
lgI_I0 = np.power(rho_bulge, 0.25)
I_bulge = I_0_bulge*np.exp(-lgI_I0)
# CS:Galaxy (Colonisation Status)
CS_gal = tools.CS_random(N_gal)
# POS:Galaxy
rho_gal = np.append(rho_disk, rho_bulge)
phi_gal = np.append(phi_disk, phi_bulge)
z_gal = np.append(z_disk, z_bulge)
V_gal = np.append(V_disk, V_bulge)
# This is just to keep the galaxy array consistent, NOT for SB measurement
I_gal = np.append(I_disk, I_bulge)

# I introduce this to be able to plot the total intensity of both bulge and disk
# which is additional in the center (within the bulge)
I_tot = I_bulge.copy()
I_tot.resize(np.shape(I_disk), refcheck=False)
I_tot += I_disk
Li = np.sum(I_gal)

# Galaxy multi-D array
galaxy = np.zeros((6, N_gal))
galaxy[0,:] = rho_gal
galaxy[1,:] = phi_gal
galaxy[2,:] = z_gal
galaxy[3,:] = V_gal
galaxy[4,:] = I_gal
galaxy[5,:] = CS_gal

galaxy_cart2 = np.zeros((2, N_gal))
galaxy_cart2[0,:] = galaxy[0,:]*np.cos(galaxy[1,:])
galaxy_cart2[1,:] = galaxy[0,:]*np.sin(galaxy[1,:])

#fig=plt.figure()
#xy=fig.add_subplot(121)
#xy.scatter(galaxy_cart2[0,:], galaxy_cart2[1,:], 'o',
#        c=np.log10(galaxy[4,:]), cmap=plt.cm.jet_r)
#xz=fig.add_subplot(122)
#xz.scatter(galaxy_cart2[0,:], galaxy[2,:], 'o',
#        c=np.log10(galaxy[4,:]), cmap=plt.cm.jet_r)
#plt.show()


def update():
    global t
    with open('logfile.txt', 'a') as logfile:
    #    UPDATING (t=t+dt)    #
    # ======================= #
        count = 1
        #cnt_p = 0
        col_tot = 0.
        log = np.zeros((3, int(t_f/dt)+1))
        for i in xrange(int(t_f/dt)+1):
        #    print '\n-----------------\nt = %f Myr\n-----------------'%(t/Myr)
            t = i*dt
            log[0, i] = t
        # Colonize the galaxy!
        #    print 'Colonizing the galaxy!'
            galaxy[4,:], galaxy[5,:], colonized, count = tools.col_single(galaxy, dist, count)
        #    ind = np.where(galaxy[5,:]==1)[0]
        #    galaxy[4,:], galaxy[5,:], colonized, count = tools.col_inf(galaxy, dist, count, ind)
            col_tot += colonized
            log[1, i] = col_tot
            log[2, i] = np.sum(galaxy[4,:])
            i += 1
        # Rotate the galaxy!
        #    print ('Rotating the galaxy!')
            galaxy[1,0:N_disk] += np.arcsin(galaxy[3,0:N_disk]*dt/galaxy[0,0:N_disk])
            galaxy_cart2[0,:] = galaxy[0,:]*np.cos(galaxy[1,:])
            galaxy_cart2[1,:] = galaxy[0,:]*np.sin(galaxy[1,:])
            # logfile.write('%f\t%f\t%f\n' %(log[0,:], log[1,:], log[2,:]))

    logfile.close()
    #       TIMING          #
    # ===================== #
    end_time=time.time()
    total_time=end_time-start_time
    #print('\n========================\nIt took %f s!\n========================'%(total_time))

cProfile.run("update()", "stats")
