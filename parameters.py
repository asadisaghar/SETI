#0
class settings(object):
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

#1
class galaxy_type(settings):
    pass
#1.1
class milky_way(galaxy_type):
    # ===================== #
    #  Galactic parameters  #
    # ===================== #
    # DISK
#    particleNo_disk = 0.9*particleNo_total
    disk_I0 = 20.e9  # disk central intensity
    disk_mass = 6.e10  # disk mass
    disk_h_rho = 

#    particleNo_thin_disk = int(0.90*particleNo_disk)
    thin_h_z = 3.e2  # scale height of the thin disk (pc)
    thin_h_rho = 5.e3  # scale length of the thin disk (pc)
    thin_alpha = disk_h_rho  # disk scale length (pc)

#    particleNo_thick_disk = int(particleNo_disk-particleNo_thin)
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
#    particleNo_bulge = particleNo_disk
    radius_bulge = 2.7e3  # bulge radius
    I_0_bulge = 5e9  # bulge central intensity
    alpha_bulge = R_bulge/3.  # bulge scale length
    n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
    mean_bulge = 200.*km/sec # (My arbitrary value!!)
    sigma_bulge = 130.*km/sec # Velocity dispercion
    mass_bulge = 2e10  # bulge mass

    # HALO
#    particleNo_dark_halo = 0.1*particleNo_total # !!!Check the value!!!
#    particleNo_stellar_halo = 0.01*particleNo_total # !!!Check the value!!!
    mass_dark_halo = 2e12  # dark halo mass(M_solar)
    mass_stellar_halo = 1e11 # stellar halo mass(M_solar) !!!Check the value!!!
    radius_dark_halo = 50e3 # (pc) !!!Check the value!!!
    radius_stellar_halo = 15e3 # (pc) !!!Check the value!!!

    #GALAXY
    alpha_gal = alpha_disk
    N_gal = N_disk + N_bulge
    M_gal = M_disk + M_bulge + M_DM
    L2Lstar = np.power((M_gal/2.e12*M_solar), 2)
    R_opt = 2.*h_rho_disk # Optical radius
    V_opt = 300.*km/sec # V(R_opt)
    n_s_gal = 4
    
#2
class civil_location(settings):
    pass

#2.1
class init_in_bulge(civil_location):
    galaxy.alien.colonizer.r_init = np.random.random()*galaxy.mysettings.galaxy_settings.radius_bulge

#2.2
class init_in_bulde(civil_location):
    galaxy.alien.colonizer.r_init = np.random.random()*galaxy.mysettings.galaxy_settings.disk_h_rho

#2.3
class init_in_sun(civil_location):
    galaxy.alien.colonizer.r_init = 8e3 # (pc)

#2.4
class init_in_halo(civil_location):
    galaxy.alien.colonizer.r_init = np.random.random()*galaxy.mysettings.galaxy_settings.radius_stellar_halo

#3
class civil_strategy(settings):
    pass

#3.1
class p2p(civil_strategy):
    pass
    
#3.2
class sphwave(civil_strategy):
    pass

#4
class calc_mode(settings):
    pass
    
#4.1
class univ_rc(calc_mode):
    pass
    
#4.2
class grav_interaction(calc_mode):
    pass

