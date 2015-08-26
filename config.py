from math import *
import numpy as np
import astropy.units as u
import astropy.constants as cst

# ===================== #
#        HANDLES        #
# ===================== #
N_disk = int(2.e4)
# time step to dump the galaxy array before the maximum colonization
dt_log =  1.*u.Myr.to(u.second)
#time step to dump the galaxy array after the maximum colonization
dt_log_stall = 20.*dt_log
dt_rotation = 0.2*u.Myr.to(u.second)
dt_colonization = 0.2*u.Myr.to(u.second)
dt_construction = 1e-12*u.Myr.to(u.second)
v_probe = 1e-4*cst.c.to(u.km/u.second)
t_initial = 0.*u.second
t_final = 15.e2*u.Myr.to(u.second)

## What does the galaxy do
galaxy_rotation = True
bulge_rotation = True
disk_rotation = True
disk_oscillation_r = True
disk_oscillation_z = True
stellar_halo_rotation = True

## What does the civilization do
colonization = True
single_probe = False
infinite_probe = not single_probe
max_colonization_fraction = 0.75
covering_fraction_colonizer = 1.0
random_start = True
if not random_start:
# Change below only if random_start = False
    start_r_colonizer = 8.e3*u.parsec  # radial distance from the galactic center
    start_r_error_colonizer = 10.*u.parsec

# ===================== #
#  Galactic parameters  #
# ===================== #
# DISK
N_thin_disk = int(0.95*N_disk)
N_thick_disk = int(N_disk-N_thin_disk)
h_z_thin_disk = 3.e2*u.parsec  # scale height of the thin disk
h_z_thick_disk = 1.e3*u.parsec # scale height of the thick disk
h_rho_disk = 3.e3*u.parsec     # scale length of the thin disk
I_0_disk = 20.e9*u.parsec      # disk central intensity
alpha_disk = h_rho_disk        # disk scale length
n_s_disk = 1                   # disk Sersic index
M_disk = 6.e10*u.M_sun         # disk mass
#FIXME values
mean_vel_rho_disk = 30.*u.km/u.second
sigma_vel_rho_disk = 30.*u.km/u.second
mean_vel_z_disk = 15.*u.km/u.second
sigma_vel_z_disk = 15.*u.km/u.second

# BULGE
N_bulge = int(0.33*N_disk)
R_bulge = 1.5e3*u.parsec             # bulge radius
I_0_bulge = 5.e9*u.L_sun             # bulge central intensity
alpha_bulge = R_bulge/3.*u.parsec    # bulge scale length
n_s_bulge = 5                        # Bulge Sersic index(in range: 1.5-10)
M_bulge = 2.e10*u.M_sun              # bulge mass
mean_vel_bulge = 100.*u.km/u.second  # (My arbitrary value!!) #FIXME
sigma_vel_bulge = 10.                # Velocity dispercion
density_bulge = M_bulge/(4./3.*pi*R_bulge**3)

# [Stellar] Halo
N_stellar_halo = int(0.01*N_disk)
R_stellar_halo = 30.e3*u.parsec
I_0_stellar_halo = 20.e7*u.L_sun       # Halo central intensity
alpha_stellar_halo = R_stellar_halo/3. # bulge scale length
n_s_stellar_halo = 5                   # Bulge Sersic index(in range: 1.5-10)
mean_vel_stellar_halo = 200.*u.km/u.second # (My arbitrary value!!)
sigma_vel_stellar_halo = 10.*u.km/u.second # Velocity dispercion
M_stellar_halo = 2.e9*u.M_sun              # Halo mass
#M_halo += M_bulge + M_disk
density_stellar_halo = M_stellar_halo/(4./3.*pi*R_stellar_halo**3)

# [Dark matter] Halo
M_dark_halo = 2.e11*u.M_sun

# Galaxy as a whole
alpha_galaxy = alpha_disk
N_galaxy = N_disk + N_bulge + N_stellar_halo
M_galaxy = M_disk + M_bulge + M_stellar_halo + M_dark_halo
L2Lstar_galaxy = np.power((M_galaxy/2.e12*cst.M_sun), 2) #L/L_stellar
R_optical_galaxy = 2.*h_rho_disk*u.parsec                # Optical radius
V_optical_galaxy = 200.*u.km/u.second                    # V(R_opt)
n_s_galaxy = 4
