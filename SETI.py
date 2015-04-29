import random
import cProfile
from math import *
import numpy as np
import tools
import datetime
import contextlib, time
@contextlib.contextmanager
def timer(msg="XXXXXX"):
    start = datetime.datetime.now()
    yield
    end = datetime.datetime.now()
    print msg, end - start
# ===================== #
#        CONSTANTS      #
# ===================== #
pc = 1.
sec = 1.
M_solar = 1.

meter = 1./3.086e16 #km
km = 1./3.086e13  # pc
year = 3.154e7*sec  # s
Myr = 1.e6*year  # s
kg = 1./1.1988435e30  # M_solar
cSpeed = 3.e5*km/sec  # pc/s
Gconst = 6.67e-11 #N.m^2.kg^-2 = kg^1.m^1.s^-2.m^2.kg^-2 = kg^-1.m^3.s^-2 
Gconst = Gconst*meter**3/(kg*sec**2)
# ===================== #
#        HANDLES        #
# ===================== #
N_disk = int(1.e5)  # number of particles in disk (same order as N_gal)
#dt = np.logspace(-1, 1, num=1)*Myr  # galactic rotation time step
dt = 5.*Myr
#dt_const = np.logspace(-12, 1, num=1)*Myr  # construction time delay
dt_const = 1e-12*Myr
#VC = np.logspace(2, 0, num=1)*cSpeed  # probe velocity
VC = 1e-3*cSpeed
t = 0
t_f = 1.e3*Myr  # time to stop
SingleProbe = False
InfiniteProbe = not(SingleProbe)
coveringFraction = 1.0
RandomStart = False
# Change below only if RandomStart = False
start_r = 8e3  # radial distance from the galactic center in pc
r_err = 100.
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
#FIXME values
mean_rho_disk = 50.*km/sec
sigma_rho_disk = 40.*km/sec
mean_z_disk = 10.*km/sec
sigma_z_disk = 3.*km/sec

# BULGE
N_bulge = int(0.33*N_disk)
#N_bulge = 0
R_bulge = 2.7e3  # bulge radius
I_0_bulge = 5.e9  # bulge central intensity
alpha_bulge = R_bulge/3.  # bulge scale length
n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
M_bulge = 9.3e6  # bulge mass #FIXME
#M_bulge = 2.e10  # bulge mass
mean_bulge = 100.*km/sec # (My arbitrary value!!) #FIXME value
sigma_bulge = 60.*km/sec # Velocity dispercion

# Halo (NOT included!)
N_halo = 0
#R_halo = 1.e8
#I_0_halo = 5.e6  # bulge central intensity
#alpha_halo = R_halo/3.  # bulge scale length
#n_s_halo = 5  # Bulge Sersic index(in range: 1.5-10)
#mean_halo = 200.*km/sec # (My arbitrary value!!)
#sigma_halo = 130.*km/sec # Velocity dispercion
#M_halo = 2.e10  # bulge mass

#GALAXY
alpha_gal = alpha_disk
N_gal = N_disk + N_bulge + N_halo
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
r_bulge = tools.init_pos(N_bulge, 0., R_bulge, 'uni')  #galaxy[0]
theta_bulge = tools.init_pos(N_bulge, 0., pi, 'uni')  #galaxy[1]
phi_bulge = tools.init_pos(N_bulge, 0., 2.*pi, 'uni')  #galaxy[2]

# POS:Bulge (Cylindrical)
#rho_bulge = r_bulge*np.sin(theta_bulge)
#z_bulge = r_bulge*np.cos(theta_bulge)
#phi_bulge = phi_bulge

# VEL:Disk (Rotation curve analytic relation)
Vrho_disk = np.zeros(N_disk)
Vphi_disk = np.zeros(N_disk)
Vz_disk = np.zeros(N_disk)

# VEL:Bulge (Gaussian dist. with given mean and dispersion)
Vr_bulge = np.zeros(N_bulge)
Vtheta_bulge = np.zeros(N_bulge)
Vphi_bulge = np.zeros(N_bulge)

# I(I-band):Disk
I_disk = tools.init_sersic(I_0_disk, alpha_disk, n_s_disk, rho_disk)

# L(I-band):Bulge (Gaussian dist. with given total L)
lgI_I0 = np.power(r_bulge, 0.25)
I_bulge = I_0_bulge*np.exp(-lgI_I0)

# POS:Galaxy
P1_gal = np.append(r_bulge, rho_disk)
P2_gal = np.append(theta_bulge, phi_disk)
P3_gal = np.append(phi_bulge, z_disk)

V1_gal = np.append(Vr_bulge, Vrho_disk)
V2_gal = np.append(Vtheta_bulge, Vphi_disk)
V3_gal = np.append(Vphi_bulge, Vz_disk)

# This is just to keep the galaxy array consistent, NOT for SB measurement
I_gal = np.append(I_bulge, I_disk)
Li = np.sum(I_gal)

# Galaxy 6-D array
galaxy = np.zeros((8, N_gal))
galaxy[0] = P1_gal
galaxy[1] = P2_gal
galaxy[2] = P3_gal
galaxy[3] = V2_gal
galaxy[4] = I_gal
# CS:Galaxy (Colonization Status)
if RandomStart:
    CS_gal = tools.CS_random(N_gal)
else:
    CS_gal, r_colonizer = tools.CS_manual(N_gal, galaxy, start_r, r_err)
galaxy[5] = CS_gal
galaxy[6] = V1_gal
galaxy[7] = V3_gal

#print "Initial state"
#print "Writing to file..."
#filename = "galaxy"
#np.save(filename, galaxy)

#def update():
#global t
# ============== #
#    UPDATING    #
# ============== #
#count = 1
i = 0
#col_tot = 0.
#count_tot = 0.
colonized_fraction = np.sum(galaxy[5])/len(galaxy[5])

while colonized_fraction < 0.7:
    print t
    t = i*dt
    # Colonize the galaxy!
    dist = VC * (dt-dt_const)
#    with timer("===========COLONIZING!==========="):
    ## how about the case where dt_const is larger?
    if SingleProbe:
        galaxy[4], galaxy[5], colonized, count, ind_dmin = tools.col_single(galaxy, galaxy_cart, dist, count, coveringFraction)
    elif InfiniteProbe:
        galaxy = tools.col_inf(galaxy, dist, coveringFraction, N_bulge=N_bulge)
#            print i
        # Evaluate bulge velocities (Spherical)
    sign = np.round(np.random.uniform(0,1,N_bulge))*2.-1
    galaxy[6,:N_bulge] = sign*np.random.normal(mean_bulge, 2.*sigma_bulge, N_bulge)
    sign = np.round(np.random.uniform(0,1,N_bulge))*2.-1
    galaxy[3,:N_bulge] = sign*np.random.normal(mean_bulge, 2.*sigma_bulge, N_bulge)
    sign = np.round(np.random.uniform(0,1,N_bulge))*2.-1
    galaxy[7,:N_bulge] = sign*np.random.normal(mean_bulge, 2.*sigma_bulge, N_bulge)
    # Rotate the bulge
    galaxy[0,:N_bulge] += galaxy[6,:N_bulge]*dt # v_r = dr/dt
    galaxy[1,:N_bulge] += galaxy[0,:N_bulge]*galaxy[3,:N_bulge]/dt # v_theta = r*dtheta/dt
    galaxy[2,:N_bulge] += galaxy[0,:N_bulge]*np.sin(galaxy[2,:N_bulge])*galaxy[7,:N_bulge]/dt # v_phi = r*sin(theta)*dphi/dt
    # Evaluate disk velocities (Cylindrical)
    sign = np.round(np.random.uniform(0,1,N_disk))*2.-1
    galaxy[6,N_bulge:] = sign*np.random.normal(mean_rho_disk, 2.*sigma_rho_disk, N_disk)
    galaxy[3,N_bulge:] = tools.v_rotational(rho_disk, V_opt, R_opt, L2Lstar)
    sign = np.round(np.random.uniform(0,1,N_disk))*2.-1
    galaxy[7,N_bulge:] = sign*np.random.normal(mean_z_disk, 2.*sigma_z_disk, N_disk)
    # Rotate the disk
    galaxy[0,N_bulge:] += galaxy[6,N_bulge:]*dt # v_rho = drho_dt
    galaxy[1,N_bulge:] += galaxy[0,N_bulge:]*galaxy[3,N_bulge:]/dt # v_phi = rho*dphi/dt
    galaxy[2,N_bulge:] += galaxy[7,N_bulge:]*dt # v_z = dz/dt
            
    colonized_fraction = np.sum(galaxy[5])/len(galaxy[5])
    i += 1

#    if np.round(colonized_fraction)==0.2 or np.round(colonized_fraction)==0.5:
    print "%.1f Myr \t %.2f colonized"%(t/Myr, colonized_fraction*100.)
    if (t/Myr)%50 == 0:
        print "%.2f colonized"%(colonized_fraction)
        print "Writing to file..."
        filename = "galaxy_%.0f"%(t/Myr)
        np.save(filename, galaxy)


print "%.2f colonized"%(colonized_fraction)
print "Writing to file..."
filename = "galaxy_%.0f"%(t/Myr)
np.save(filename, galaxy)
#cProfile.run("update()", "stats")
