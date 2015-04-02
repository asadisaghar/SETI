#FIXME
# initialize the galaxy
# let the gal rotate for a Myr or so, dump the whole galaxy array into a file to plot
# turn on the colonization and at each [few] time step[s] dump the location of colonized/colonizer(s) into a txt file
# overplot the colonization patter for whatever time step you want from the dumped measurements 

# ranges of galactic parameters after initialization, why are they not in expected ranges?
# enforce orbital velocities for Vphi_bulge/halo
# realistic Vrho/z_*
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
N_disk = int(1.e6)  # number of particles in disk (same order as N_gal)
dt = np.logspace(-2, 1, num=1)*Myr  # galactic rotation time step
dt_const = np.logspace(-6, 1, num=1)*Myr  # construction time delay
VC = np.logspace(-1, 0, num=1)*cSpeed  # probe velocity
t = 0
t_f = 1.e4*Myr  # time to stop
SingleProbe = False
InfiniteProbe = not(SingleProbe)
coveringFraction = 1.0
RandomStart = False
# Change below only if RandomStart = False
start_r = 1.5e3  # radial distance from the galactic center in pc
r_err = 100.
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
#FIXME values
mean_rho_disk = 200.*km/sec
sigma_rho_disk = 130.*km/sec
mean_z_disk = 20.*km/sec
sigma_z_disk = 13.*km/sec

# BULGE
#N_bulge = int(0.33000*N_disk)
## Maybe this is too small after all...
N_bulge = int(0.01*N_disk)
R_bulge = 2.7e3  # bulge radius
I_0_bulge = 5.e9  # bulge central intensity
alpha_bulge = R_bulge/3.  # bulge scale length
n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
M_bulge = 2.e10  # bulge mass
mean_bulge = 200.*km/sec # (My arbitrary value!!) #FIXME value
sigma_bulge = 130.*km/sec # Velocity dispercion

# Halo (NOT included!)
N_halo = 0
R_halo = 1.e8
I_0_halo = 5.e6  # bulge central intensity
alpha_halo = R_halo/3.  # bulge scale length
n_s_halo = 5  # Bulge Sersic index(in range: 1.5-10)
mean_halo = 200.*km/sec # (My arbitrary value!!)
sigma_halo = 130.*km/sec # Velocity dispercion
M_halo = 2.e10  # bulge mass

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
r_bulge = tools.init_pos(N_bulge, 0., R_bulge, 'uni')
phi_bulge = tools.init_pos(N_bulge, 0., pi, 'uni')
theta_bulge = tools.init_pos(N_bulge, 0., 2.*pi, 'uni')

# POS:Bulge (Cylindrical)
rho_bulge = r_bulge*np.sin(phi_bulge)
z_bulge = r_bulge*np.cos(phi_bulge)
phi_bulge = theta_bulge

# POS:Halo (Spherical)
r_halo = tools.init_pos(N_halo, 0., R_halo, 'uni')
phi_halo = tools.init_pos(N_halo, 0., pi, 'uni')
theta_halo = tools.init_pos(N_halo, 0., 2.*pi, 'uni')

# POS:Halo (Cylindrical)
rho_halo = r_halo*np.sin(phi_halo)
z_halo = r_halo*np.cos(phi_halo)
phi_halo = theta_halo

# VEL:Disk (Rotation curve analytic relation)
Vrho_disk = np.random.normal(mean_rho_disk, 2.*sigma_rho_disk, N_disk)
Vphi_disk = tools.v_rotational(rho_disk, V_opt, R_opt, L2Lstar)
Vz_disk = np.random.normal(mean_z_disk, 2.*sigma_z_disk, N_disk)

# VEL:Bulge (Gaussian dist. with given mean and dispersion)
Vrho_bulge = np.random.normal(mean_bulge, 2.*sigma_bulge, N_bulge)
Vphi_bulge = np.random.normal(mean_bulge, 2.*sigma_bulge, N_bulge)
omega_bulge = np.sqrt(Gconst*M_bulge/(R_bulge)**3)  # [1/s]
Vphi_bulge = r_bulge*km*omega_bulge  # [km/s] 
Vz_bulge = np.random.normal(mean_bulge, 2.*sigma_bulge, N_bulge)

# VEL:Halo (Gaussian dist. with given mean and dispersion)
Vrho_halo = tools.init_norm(mean_halo, sigma_halo, N_halo)
Vphi_halo = tools.init_norm(mean_halo, sigma_halo, N_halo)
Vz_halo = tools.init_norm(mean_halo, sigma_halo, N_halo)

# I(I-band):Disk
I_disk = tools.init_sersic(I_0_disk, alpha_disk, n_s_disk, rho_disk)

# L(I-band):Bulge (Gaussian dist. with given total L)
lgI_I0 = np.power(rho_bulge, 0.25)
I_bulge = I_0_bulge*np.exp(-lgI_I0)

# L(I-band):Halo (Gaussian dist. with given total L)
lgI_I0 = np.power(rho_halo, 0.25)
I_halo = I_0_halo*np.exp(-lgI_I0)

# POS:Galaxy
rho_gal = np.append(rho_disk, rho_bulge)
rho_gal = np.append(rho_gal, rho_halo)
phi_gal = np.append(phi_disk, phi_bulge)
phi_gal = np.append(phi_gal, phi_halo)
z_gal = np.append(z_disk, z_bulge)
z_gal = np.append(z_gal, z_halo)
Vrho_gal = np.append(Vrho_disk, Vrho_bulge)
Vrho_gal = np.append(Vrho_gal, Vrho_halo)
Vphi_gal = np.append(Vphi_disk, Vphi_bulge)
Vphi_gal = np.append(Vphi_gal, Vphi_halo)
Vz_gal = np.append(Vz_disk, Vz_bulge)
Vz_gal = np.append(Vz_gal, Vz_halo)

#FIXME WTF??
# This is just to keep the galaxy array consistent, NOT for SB measurement
I_gal = np.append(I_disk, I_bulge)
I_gal = np.append(I_gal, I_halo)
Li = np.sum(I_gal)

# Galaxy 6-D array
galaxy = np.zeros((8, N_gal))
galaxy[0] = rho_gal
galaxy[1] = phi_gal
galaxy[2] = z_gal
galaxy[3] = Vphi_gal
galaxy[4] = I_gal
# CS:Galaxy (Colonisation Status)
if RandomStart:
    CS_gal = tools.CS_random(N_gal)
else:
    CS_gal, r_colonizer = tools.CS_manual(N_gal, galaxy, start_r, r_err)
galaxy[5] = CS_gal
galaxy[6] = Vrho_gal
galaxy[7] = Vz_gal

#galaxy_inds = {"rho": 0, "phi": 1, "z", 2, "Vphi": 3, "I": 4, "CS": 5, "Vrho": 6, "Vz": 7}


## Initial plotting
x_gal=rho_gal*np.cos(phi_gal)
y_gal=rho_gal*np.sin(phi_gal)
print "Initial"
tools.plot_cont_galaxy(t, x_gal, y_gal, z_gal, I_gal)

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

    while col_tot/N_gal < 0.45:
        t = i*dt
    # Colonize the galaxy!
        dist = VC * (dt-dt_const)
        ## how about the case where dt_const is larger?
        if SingleProbe:
            report = np.zeros(((int(t_f/dt)+1), 4))
            report[i, 0] = t/Myr
            galaxy[4], galaxy[5], colonized, count, ind_dmin = tools.col_single(galaxy, dist, count, coveringFraction)
            col_tot += colonized
            report[i, 1] = np.sum(galaxy[4,:]/Li)
            report[i, 2] = np.sum(abs(galaxy[5,:])/N_gal)
            report[i, 3] = galaxy[0, ind_dmin]
            logfile.write('%e\t%e\t%e\t%f\n' %(report[i,0], report[i,1], report[i,2], report[i,3]))
        elif InfiniteProbe:
            report = np.zeros(((int(t_f/dt)+1), 4))
            report[i, 0] = t/Myr
            ind = np.where(CS_gal==1)[0]
            galaxy[4], galaxy[5], colonized, count = tools.col_inf(galaxy, dist, count, ind)
            col_tot += colonized
            report[i, 1] = np.sum(galaxy[4,:]/Li)
            report[i, 2] = np.sum(abs(galaxy[5,:])/N_gal)
            logfile.write('%e\t%e\t%e\t%f\n' %(report[i,0], report[i,1], report[i,2], report[i,3]))
    # Rotate the galaxy!
        galaxy[0] += galaxy[6]*dt
        galaxy[1] += galaxy[3]*dt/galaxy[0,:]
        galaxy[2] += galaxy[7]*dt

        print '\tt = %.2f Myr\n-------%d = %.2f N_gal------'%(t/Myr, col_tot, col_tot/N_gal)
        if i==500:
            galaxy_cart = np.zeros((3,(N_bulge+N_disk+N_halo)))
            galaxy_cart[0]=galaxy[0]*np.cos(galaxy[1])
            galaxy_cart[1]=galaxy[0]*np.sin(galaxy[1])
            galaxy_cart[2]=galaxy[2]
            print "step%d"%(i)
            tools.plot_cont_galaxy(t, galaxy_cart[0], galaxy_cart[1], galaxy_cart[2], galaxy[4])

        i += 1
#galaxy_inds = {"rho": 0, "phi": 1, "z", 2, "Vphi": 3, "I": 4, "CS": 5, "Vrho": 6, "Vz": 7}

    logfile.close()

cProfile.run("update()", "stats")

# ==================== #
#        PLOTTING      #
# ==================== #
#if SingleProbe:
#    tools.singleplot(name, N_gal, Li, r_colonizer, VC, dt_const)
#elif InfiniteProbe:
#    tools.infplot(name, N_gal, Li, r_colonizer, VC, dt_const)
galaxy_cart = np.zeros((3,(N_bulge+N_disk+N_halo)))
galaxy_cart[0]=galaxy[0]*np.cos(galaxy[1])
galaxy_cart[1]=galaxy[0]*np.sin(galaxy[1])
galaxy_cart[2]=galaxy[2]
print "Final"
tools.plot_cont_galaxy(t, galaxy_cart[0], galaxy_cart[1], galaxy_cart[2], galaxy[4])
#tools.plot_cont_galaxy(galaxy_cart)
