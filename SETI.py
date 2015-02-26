import random
import cProfile
from math import *
import numpy as np
import tools
# ===================== #
#        CONSTANTS      #
# ===================== #
pc = 1.
sec = 1.
M_solar = 1.

km = 1./3.086e13  # pc
year = 3.154e7*sec  # s
Myr = 1.e6*year  # s
kg = 1./1.1988435e30  # M_solar
cSpeed = 3.e5*km/sec  # pc/s
# ===================== #
#        HANDLES        #
# ===================== #
N_disk = int(1.e4)  # number of particles in disk (same order as N_gal)
dt = np.logspace(-3, 1, num=1)*Myr  # galactic rotation time step
dt_const = np.logspace(-6, 1, num=1)*Myr  # construction time delay
VC = np.logspace(-0.5, 0, num=1)*cSpeed  # probe velocity
t_f = 1.e4*Myr  # time to stop
SingleProbe = True
InfiniteProbe = not(SingleProbe)
coveringFraction = 1.0
RandomStart = False
# Change below only if RandomStart = False
start_r = 8.e2  # radial distance from the galactic center in pc
r_err = 10.
name = 'bulge_cf1_single'
# ===================== #
#  Galactic parameters  #
# ===================== #
# DISK
N_thin = int(0.90*N_disk)
N_thick = int(N_disk-N_thin)
h_z_thin = 3.e2  # scale height of the thin disk (pc)
h_z_thick = 5.e2  # scale height of the thick disk (pc)
h_rho_disk = 5.e3  # scale length of the thin disk (pc)
I_0_disk = 20.e9  # disk central intensity
alpha_disk = h_rho_disk  # disk scale length (pc)
n_s_disk = 1  # disk Sersic index
M_disk = 6.e10  # disk mass

# BULGE
#N_bulge = int(0.33000*N_disk)
## Maybe this is too small after all...
N_bulge = N_disk
R_bulge = 2.7e3  # bulge radius
I_0_bulge = 5.e9  # bulge central intensity
alpha_bulge = R_bulge/3.  # bulge scale length
n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
mean_bulge = 200.*km/sec # (My arbitrary value!!)
sigma_bulge = 130.*km/sec # Velocity dispercion
M_bulge = 2.e10  # bulge mass

#GALAXY
alpha_gal = alpha_disk
N_gal = N_disk + N_bulge
M_DM = 2.e12  # dark halo mass(M_solar)
M_gal = M_disk + M_bulge + M_DM
L2Lstar = np.power((M_gal/2.e12*M_solar), 2)
R_opt = 2.*h_rho_disk # Optical radius
V_opt = 300.*km/sec # V(R_opt)
n_s_gal = 4
# ========================== #
#    Initialisation (t=0)    #
# ========================== #
# POS:Thin disk (Cylindrical)
rho_thin= tools.init_pos(N_thin, 0., h_rho_disk, 'exp')
phi_thin = tools.init_pos(N_thin, 0., 2*pi, 'uni')
z_thin = tools.init_z(N_thin, 0., h_z_thin, 'exp')

# POS:Thick disk (Cylindrical)
rho_thick = tools.init_pos(N_thick, 0, h_rho_disk, 'exp')
phi_thick = tools.init_pos(N_thick, 0., 2*pi, 'uni')
z_thick = tools.init_z(N_thick, 0., h_z_thick, 'exp')

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
V_disk = tools.v_rotational(rho_disk, V_opt, R_opt, L2Lstar)

# VEL:Bulge (Gaussian dist. with given mean and dispersion)
V_bulge = tools.init_norm(mean_bulge, sigma_bulge, N_bulge)

# I(I-band):Disk
I_disk = tools.init_sersic(I_0_disk, alpha_disk, n_s_disk, rho_disk)

# L(I-band):Bulge (Gaussian dist. with given total L)
lgI_I0 = np.power(rho_bulge, 0.25)
I_bulge = I_0_bulge*np.exp(-lgI_I0)

# POS:Galaxy
rho_gal = np.append(rho_disk, rho_bulge)
phi_gal = np.append(phi_disk, phi_bulge)
z_gal = np.append(z_disk, z_bulge)
V_gal = np.append(V_disk, V_bulge)

# This is just to keep the galaxy array consistent, NOT for SB measurement
I_gal = np.append(I_disk, I_bulge)
Li = np.sum(I_gal)

# Galaxy 6-D array
galaxy = np.zeros((6, N_gal))
galaxy[0,:] = rho_gal
galaxy[1,:] = phi_gal
galaxy[2,:] = z_gal
galaxy[3,:] = V_gal
galaxy[4,:] = I_gal
# CS:Galaxy (Colonisation Status)
if RandomStart:
    CS_gal = tools.CS_random(N_gal)
else:
    CS_gal, r_colonizer = tools.CS_manual(N_gal, galaxy, start_r, r_err)
galaxy[5,:] = CS_gal

def update():
    global t
    logfile = open('%s.txt'%(name), 'w')
# ============== #
#    UPDATING    #
# ============== #
    count = 1
    i = 0
    col_tot = 0.
    count_tot = 0.
    while col_tot/N_gal < 0.50:
        t = i*dt
    # Colonize the galaxy!
        dist = VC * (dt-dt_const)
        ## how about the case where dt_const is larger?
        if SingleProbe:
            logg = np.zeros(((int(t_f/dt)+1), 4))
            logg[i, 0] = t/Myr
            galaxy[4,:], galaxy[5,:], colonized, count, ind_dmin = tools.col_single(galaxy, dist, count, coveringFraction)
            col_tot += colonized
            logg[i, 1] = np.sum(galaxy[4,:]/Li)
            logg[i, 2] = np.sum(abs(galaxy[5,:])/N_gal)
            logg[i, 3] = galaxy[0, ind_dmin]
            logfile.write('%e\t%e\t%e\t%f\n' %(logg[i,0], logg[i,1], logg[i,2], logg[i,3]))
        elif InfiniteProbe:
            logg = np.zeros(((int(t_f/dt)+1), 3))
            logg[i, 0] = t/Myr
            ind = np.where(CS_gal==1)[0]
            galaxy[4,:], galaxy[5,:], colonized, count = tools.col_inf(galaxy, dist, count, ind)
            col_tot += colonized
            logg[i, 1] = np.sum(galaxy[4,:]/Li)
            logg[i, 2] = np.sum(abs(galaxy[5,:])/N_gal)
            logfile.write('%e\t%e\t%e\t%f\n' %(logg[i,0], logg[i,1], logg[i,2], logg[i,3]))
    # Rotate the galaxy!
        galaxy[1,0:N_disk] += galaxy[3,0:N_disk]*dt/galaxy[0,0:N_disk]
        galaxy[1,N_disk:-1] += galaxy[3,N_disk:-1]*dt/galaxy[0,N_disk:-1]
        print '\tt = %.2f Myr\n-------%d = %.2f N_gal------'%(t/Myr, col_tot, col_tot/N_gal)

        i += 1

    logfile.close()

cProfile.run("update()", "stats")

# ==================== #
#        PLOTTING      #
# ==================== #
if SingleProbe:
    tools.singleplot(name, N_gal, Li, r_colonizer, VC, dt_const)
elif InfiniteProbe:
    tools.infplot(name, N_gal, Li, r_colonizer, VC, dt_const)

