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
# Multiply everything by relevant units below to translate them into pc, s, and M_solar
km2m = 1.e3
m2km = 1./km2m
kpc2pc = 1.e3
pc2kpc = 1./kpc2pc
pc2m = 3.08567758e16
m2pc = 1./pc2m
pc2km = pc2m*m2km
km2pc = km2m*m2pc
Msolar2kg = 1.9891e30
kg2Msolar = 1./Msolar2kg
yr2sec = pi*1e7
sec2yr = 1./yr2sec
Myr2yr = 1.e6
yr2Myr = 1./Myr2yr
Myr2sec = Myr2yr*yr2sec
sec2Myr = 1./Myr2sec

cSpeed = 3.e5 #[km/s]
# ===================== #
#        HANDLES        #
# ===================== #
N_disk = int(5.e4)  # number of particles in disk (same order as N_gal)
#dt = np.logspace(-1, 1, num=1)*Myr  # galactic rotation time step
dt = 0.01*Myr2sec #[s] #do NOT use anything longer than 0.1 Myr, or disk rotation is messed up!!
#dt_const = np.logspace(-12, 1, num=1)*Myr  # construction time delay
dt_const = 1e-12*Myr2sec #[s]
#VC = np.logspace(2, 0, num=1)*cSpeed  # probe velocity
VC = 1e-3*cSpeed #[km/s]
t = 1
t_f = 1.5e3*Myr2sec  # time to stop #[s]
SingleProbe = False
probe = "inf"
InfiniteProbe = not(SingleProbe)
coveringFraction = 1.0
RandomStart = False
# Change below only if RandomStart = False
start_r = 1e3  # radial distance from the galactic center in #[pc]
loc = 1
r_err = 100. #[pc]
# ===================== #
#  Galactic parameters  #
# ===================== #
# DISK
N_thin = int(0.95*N_disk)
N_thick = int(N_disk-N_thin)
h_z_thin = 3.e2  # scale height of the thin disk #[pc]
h_z_thick = 1.e3  # scale height of the thick disk #[pc]
h_rho_disk = 3.e3  # scale length of the thin disk #[pc]
I_0_disk = 20.e9  # disk central intensity
alpha_disk = h_rho_disk  # disk scale length #[pc]
n_s_disk = 1  # disk Sersic index
M_disk = 6.e10  # disk mass #[M_solar]
#FIXME values
mean_rho_disk = 30. #[km/s]
sigma_rho_disk = 30.  #[km/s]
mean_z_disk = 15.  #[km/s]
sigma_z_disk = 15.  #[km/s]

# BULGE
N_bulge = int(0.33*N_disk)
#N_bulge = 0
R_bulge = 1.5e3  # bulge radius #[pc]
I_0_bulge = 5.e9  # bulge central intensity
alpha_bulge = R_bulge/3.  # bulge scale length #[pc]
n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
#M_bulge = 9.3e6  # bulge mass #FIXME
M_bulge = 2.e10  # bulge mass #[M_solar]
mean_bulge = 100. # (My arbitrary value!!) #FIXME value  #[km/s]
sigma_bulge = 100. # Velocity dispercion  #[km/s]
rho_bulge = M_bulge/(4./3.*pi*R_bulge**3)

# Halo
N_halo = int(0.01*N_bulge)
#N_halo = N_bulge
R_halo = 30.e3
I_0_halo = 20.e7  # Halo central intensity
alpha_halo = R_halo/3.  # bulge scale length
n_s_halo = 5  # Bulge Sersic index(in range: 1.5-10)
mean_halo = 200. # (My arbitrary value!!)
sigma_halo = 200. # Velocity dispercion
M_halo = 2.e9  # bulge mass
rho_halo = M_halo/(4./3.*pi*R_halo**3)

#GALAXY
alpha_gal = alpha_disk
N_gal = N_disk + N_bulge + N_halo
M_DM = 2.e11  # dark halo mass(M_solar) #[M_solar]
M_gal = M_disk + M_bulge + M_DM #[M_solar]
L2Lstar = np.power((M_gal/2.e12*M_solar), 2)
R_opt = 2.*h_rho_disk # Optical radius #[pc]
V_opt = 200. # V(R_opt)  #[km/s]
n_s_gal = 4
# ========================== #
#    Initialisation (t=0)    #
# ========================== # 
# POS:Thin disk (Cylindrical)
rho_thin= tools.init_pos(N_thin, 0., h_rho_disk, 'exp') #[pc]
phi_thin = tools.init_pos(N_thin, 0., 2*pi, 'uni')  #[pc]
z_thin = tools.init_z(N_thin, 0., h_z_thin, 'exp')  #[pc]

# POS:Thick disk (Cylindrical)
rho_thick = tools.init_pos(N_thick, 0., h_rho_disk, 'exp')  #[pc]
phi_thick = tools.init_pos(N_thick, 0., 2*pi, 'uni')  #[pc]
z_thick = tools.init_z(N_thick, 0., h_z_thick, 'exp')  #[pc]

# POS:Bulge (Spherical)
r_bulge = abs(np.random.normal(0, R_bulge, N_bulge))  #galaxy[0]  #[pc]
theta_bulge = tools.init_pos(N_bulge, 0., np.nextafter(pi,4), 'uni')  #galaxy[1]  #[pc]
phi_bulge = tools.init_pos(N_bulge, 0., 2.*pi, 'uni')  #galaxy[2]  #[pc]

# POS:Disk (Cylindrical)
rho_disk = np.append(rho_thin, rho_thick)  #[pc]
phi_disk = np.append(phi_thin, phi_thick)  #[pc]
z_disk = np.append(z_thin, z_thick)  #[pc]

# POS:Halo (Spherical)
r_halo = abs(np.random.normal(0, R_halo, N_halo))  #galaxy[0]  #[pc]
theta_halo = tools.init_pos(N_halo, 0., np.nextafter(pi,4), 'uni')  #galaxy[1]  #[pc]
phi_halo = tools.init_pos(N_halo, 0., 2.*pi, 'uni')  #galaxy[2]  #[pc]

# VEL:Bulge (Gaussian dist. with given mean and dispersion)
Vr_bulge = np.zeros(N_bulge) #[km/s]
Vtheta_bulge = np.zeros(N_bulge) #[km/s]
Vphi_bulge = np.zeros(N_bulge) #[km/s]

# VEL:Disk (Rotation curve analytic relation)
Vrho_disk = np.zeros(N_disk) #[km/s]
Vphi_disk = np.zeros(N_disk)  #[km/s]
Vz_disk = np.zeros(N_disk)  #[km/s]

# VEL:Halo (Gaussian dist. with given mean and dispersion)
Vr_halo = np.zeros(N_halo) #[km/s]
Vtheta_halo = np.zeros(N_halo) #[km/s]
Vphi_halo = np.zeros(N_halo) #[km/s]

# L(I-band):Bulge (Gaussian dist. with given total L)
lgI_I0 = np.power(r_bulge, 0.25)
I_bulge = I_0_bulge*np.exp(-lgI_I0)

# I(I-band):Disk
I_disk = tools.init_sersic(I_0_disk, alpha_disk, n_s_disk, rho_disk)

# L(I-band):Halo (Gaussian dist. with given total L)
lgI_I0 = np.power(r_halo, 0.25)
I_halo = I_0_halo*np.exp(-lgI_I0)

# POS:Galaxy
P1_gal = np.append(r_bulge, rho_disk) #[pc, pc]
P1_gal = np.append(P1_gal, r_halo)
P2_gal = np.append(theta_bulge, phi_disk) #[,]
P2_gal = np.append(P2_gal, theta_halo)
P3_gal = np.append(phi_bulge, z_disk) #[,pc]
P3_gal = np.append(P3_gal, phi_halo)

V1_gal = np.append(Vr_bulge, Vrho_disk) #[km/s, km/s]
V1_gal = np.append(V1_gal, Vr_halo)
V2_gal = np.append(Vtheta_bulge, Vphi_disk) #[km/s, km/s]
V2_gal = np.append(V2_gal, Vtheta_halo)
V3_gal = np.append(Vphi_bulge, Vz_disk) #[km/s, km/s]
V3_gal = np.append(V3_gal, Vphi_halo)

# This is just to keep the galaxy array consistent, NOT for SB measurement
I_gal = np.append(I_bulge, I_disk)
I_gal = np.append(I_gal, I_halo)

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
count = 1
i = 0
#col_tot = 0.
#count_tot = 0.
colonized_fraction = abs(np.sum(galaxy[5]))/len(galaxy[5])

while t < t_f:
#    print t
    t = i*dt #[sec]
    # Colonize the galaxy!
    dist = VC * (dt-dt_const)*km2pc #[pc]
#    with timer("===========COLONIZING!==========="):
    ## how about the case where dt_const is larger?
    if abs(colonized_fraction-0.50) >= 0.05:
        if SingleProbe:
            galaxy, count = tools.col_sing(galaxy, dist, count, coveringFraction, N_bulge=N_bulge, N_disk=N_disk)
    #        galaxy[4], galaxy[5], colonized, count, ind_dmin = tools.col_single(galaxy, galaxy_cart, dist, count, coveringFraction)
        elif InfiniteProbe:
            print 'pre-count: %d'%(count)
            galaxy, count = tools.col_inf2(galaxy, dist, count, coveringFraction, N_bulge=N_bulge, N_disk=N_disk)
            print 'post-count: %d'%(count)

    # Evaluate bulge velocities (Spherical)
    galaxy[3,:N_bulge] = tools.v_rotational_unisphere(galaxy[0,:N_bulge], rho_bulge) #[km/s]
    sign = np.round(np.random.uniform(0,1,N_bulge))*2.-1
    galaxy[7,:N_bulge] = sign*np.random.normal(0.0, sigma_bulge, N_bulge) #[km/s]
    # Rotate the bulge
    galaxy[1,:N_bulge] += galaxy[3,:N_bulge]*km2pc*dt/galaxy[0,:N_bulge] # v_theta = r*dtheta/dt  #[pc]
    galaxy[2,:N_bulge] += galaxy[7,:N_bulge]*km2pc*dt/galaxy[0,:N_bulge]*np.sin(galaxy[2,:N_bulge]) # v_phi = r*sin(theta)*dphi/dt #[pc]
    # Oscillate bulge particles (around original r positions only)
    galaxy[0,:N_bulge] = tools.r_bulge_oscillation(galaxy, t, mean_bulge, N_bulge, R_bulge)   
     
    # Evaluate disk velocities (Cylindrical)
    galaxy[3,N_bulge:N_disk+N_bulge] = tools.v_rotational_disk(galaxy[0,N_bulge:N_disk+N_bulge]*pc2km, V_opt, R_opt*pc2km, L2Lstar) #[km/s]
    sign = np.round(np.random.uniform(0,1,N_disk))*2.-1
    galaxy[3,N_bulge:N_bulge+N_disk] += sign*np.random.normal(0.0, sigma_rho_disk, N_disk) #[km/s]
    # Rotate the disk
    galaxy[1,N_bulge:N_disk+N_bulge] += galaxy[3,N_bulge:N_disk+N_bulge]*km2pc*dt/galaxy[0,N_bulge:N_disk+N_bulge] # v_phi = rho*dphi/dt #[pc]
    # Oscillate disk particles (around original rho and z positions only)
    galaxy[0,N_bulge:N_disk+N_bulge] = rho_disk + tools.r_oscillation(galaxy, t, mean_rho_disk, N_bulge, N_disk, amp=h_z_thick)
    galaxy[2,N_bulge:N_disk+N_bulge] = tools.z_oscillation(galaxy, t, mean_z_disk, N_bulge, N_disk, amp=z_disk)    

    # Evaluate halo velocities (Spherical)
    galaxy[3,N_bulge+N_disk:] = tools.v_rotational_unisphere(galaxy[0,N_bulge+N_disk:], rho_halo) #[km/s]
    sign = np.round(np.random.uniform(0,1,N_halo))*2.-1
    galaxy[7,N_bulge+N_disk:] = sign*np.random.normal(0.0, sigma_halo, N_halo) #[km/s]
    # Rotate the halo
    galaxy[1,N_bulge+N_disk:] += galaxy[3,N_bulge+N_disk:]*km2pc*dt/galaxy[0,N_bulge+N_disk:] # v_theta = r*dtheta/dt  #[pc]
    galaxy[2,N_bulge+N_disk:] += galaxy[7,N_bulge+N_disk:]*km2pc*dt/galaxy[0,N_bulge+N_disk:]*np.sin(galaxy[2,N_bulge+N_disk:]) # v_phi = r*sin(theta)*dphi/dt #[pc]
    # Oscillate halo particles (around original r positions only)
    galaxy[0,N_bulge+N_disk:] = tools.r_halo_oscillation(galaxy, t, mean_halo, N_bulge, N_disk, N_halo, R_halo)
            
    colonized_fraction = abs(np.sum(galaxy[5]))/len(galaxy[5])
    i += 1

#    if np.round(colonized_fraction)==0.2 or np.round(colonized_fraction)==0.5:
    print "%.2f Myr \t %.2f colonized"%(t*sec2Myr, colonized_fraction*100.)
    if abs(colonized_fraction-0.50) >= 0.05:
        if int(t*sec2Myr)%0.01 == 0:
            print "%.2f colonized"%(colonized_fraction)
            print "Writing to file..."
            filename = "/home/saas9842/PhD/tmp/%s%d/galaxy_%.2d"%(probe, loc, int(t*sec2Myr*100))
            np.save(filename, galaxy)

    else:
        if int(t*sec2Myr)%100 == 0:
            print "%.2f colonized"%(colonized_fraction)
            print "Writing to file..."
            filename = "/home/saas9842/PhD/tmp/%s%d/galaxy_%.2d"%(probe, loc, int(t*sec2Myr*100))
            np.save(filename, galaxy)

print "%.2f colonized"%(colonized_fraction)
print "Writing to file..."
filename = "/home/saas9842/PhD/tmp/%s%d/galaxy_%.2d"%(probe, loc, int(t*sec2Myr*100))
np.save(filename, galaxy)
#cProfile.run("update()", "stats")
