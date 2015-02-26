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
#  Galactic parameters  #
# ===================== #
# DISK
disk_I0 = 20.e9  # disk central intensity
disk_mass = 6.e10  # disk mass
disk_h_rho = 

thin_partNo = int(0.90*particleNo_disk)
thin_h_z = 3.e2  # scale height of the thin disk (pc)
thin_h_rho = 5.e3  # scale length of the thin disk (pc)
thin_alpha = disk_h_rho  # disk scale length (pc)

thick_partNo = int(particleNo_disk-particleNo_thin)
thick_h_z = 5.e2  # scale height of the thick disk (pc)
thick_h_rho = thin_h_rho

thick_alpha = thin_alpha
thin_n_s = 1  # disk Sersic index
thick_n_s = thin_n_s
thin_mass = disk_mass*0.9
thick_mass = disk_mass - thin_mass

# BULGE
#N_bulge = int(0.33000*N_disk)
## Maybe this is too small after all...
particleNo_bulge = particleNo_disk
radius_bulge = 2.7e3  # bulge radius
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

