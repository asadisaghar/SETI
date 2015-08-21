import config
from math import *
import random
RS = 200
random.seed(RS)
import numpy as np
import astropy.units as u
import astropy.constants as cst
import random
import matplotlib.pyplot as plt
plt.rc("font", size=20, family='serif', weight='normal')
plt.rc("axes", labelsize=16, titlesize=20)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)


def init_pos(N, low, high_scale, dist):
    random.seed(RS)
    if dist == 'exp':
        init = np.random.exponential(high_scale, N)
    elif dist == 'uni':
        init = np.random.uniform(low, high_scale, N)

    return init


def init_z(N, low, high_scale, dist):
    random.seed(RS)
    if dist == 'exp':
        init = np.random.exponential(high_scale, N)
    elif dist == 'uni':
        init = np.random.uniform(low, high_scale, N)

    sign = np.random.random_integers(low=0, high=N-1, size=N/2)
    for i in sign:
        init[i] *= -1

    return init

def init_sersic(I_0, alpha, n_s, r):
#    L = I_0*np.exp((-r/alpha)**(1./n_s))
    L = I_0
    return L


def initialize_bulge():
    # POS:Bulge (Spherical)
    global pos_r_bulge
    pos_r_bulge = abs(np.random.normal(0, config.R_bulge, config.N_bulge))
    global pos_theta_bulge 
    pos_theta_bulge = init_pos(config.N_bulge, 0., np.nextafter(pi,4), 'uni')
    global pos_phi_bulge
    pos_phi_bulge = init_pos(config.N_bulge, 0., 2.*pi, 'uni')
    
    # VEL:Bulge (Gaussian dist. with given mean and dispersion)
    global vel_r_bulge
    vel_r_bulge = np.zeros(config.N_bulge) #FIXME
    global vel_theta_bulge
    vel_theta_bulge = np.ones(config.N_bulge)*config.mean_vel_bulge
    global vel_phi_bulge
    vel_phi_bulge = np.ones(config.N_bulge)*config.mean_vel_bulge/5.
    
    # L(I-band):Bulge (Gaussian dist. with given total L)
    lgI_I0 = np.power(pos_r_bulge, 0.25)
    global I_bulge
    I_bulge = config.I_0_bulge*np.exp(-lgI_I0)
    
    
def initialize_disk():
    # POS:Thin disk (Cylindrical)
    rho_thin = init_pos(config.N_thin_disk, 0., config.h_rho_disk, 'exp')
    phi_thin = init_pos(config.N_thin_disk, 0., 2*pi, 'uni')
    z_thin = init_z(config.N_thin_disk, 0., config.h_z_thin_disk, 'exp')

    # POS:Thick disk (Cylindrical)
    rho_thick = init_pos(config.N_thick_disk, 0., config.h_rho_disk, 'exp')
    phi_thick = init_pos(config.N_thick_disk, 0., 2*pi, 'uni')
    z_thick = init_z(config.N_thick_disk, 0., config.h_z_thick_disk, 'exp')

    # POS:Disk (Cylindrical)
    global pos_rho_disk
    pos_rho_disk = np.append(rho_thin, rho_thick)
    global pos_phi_disk
    pos_phi_disk = np.append(phi_thin, phi_thick)
    global pos_z_disk
    pos_z_disk = np.append(z_thin, z_thick)

    # VEL:Disk (Rotation curve analytic relation)
    global vel_rho_disk
    vel_rho_disk = np.zeros(config.N_disk)
    global vel_phi_disk
    vel_phi_disk = np.zeros(config.N_disk)
    global vel_z_disk
    vel_z_disk = np.zeros(config.N_disk)

    # I(I-band):Disk
    global I_disk
    I_disk = init_sersic(config.I_0_disk, config.alpha_disk, config.n_s_disk, pos_rho_disk)


def initialize_stellar_halo():
    # POS:Halo (Spherical)
    global pos_r_stellar_halo
    pos_r_stellar_halo = abs(np.random.normal(config.R_stellar_halo, config.R_stellar_halo, config.N_stellar_halo))
    global pos_theta_stellar_halo
    pos_theta_stellar_halo = init_pos(config.N_stellar_halo, 0., np.nextafter(pi,4), 'uni')
    global pos_phi_stellar_halo
    pos_phi_stellar_halo = init_pos(config.N_stellar_halo, 0., 2.*pi, 'uni')
    
    # VEL:Halo (Gaussian dist. with given mean and dispersion)
    global vel_r_stellar_halo
    vel_r_stellar_halo = np.zeros(config.N_stellar_halo)
    global vel_theta_stellar_halo
    vel_theta_stellar_halo = np.ones(config.N_stellar_halo)*config.mean_vel_stellar_halo
    global vel_phi_stellar_halo
    vel_phi_stellar_halo = np.ones(config.N_stellar_halo)*config.mean_vel_stellar_halo/5.
        
    # L(I-band):Halo (Gaussian dist. with given total L)
    lgI_I0 = np.power(pos_r_stellar_halo, 0.25)
    global I_stellar_halo
    I_stellar_halo = config.I_0_stellar_halo*np.exp(-lgI_I0)


#FIXME
#This function does NOT work because:
def initialize_users_colonizer(N_gal, galaxy, start_r, r_err):
    CS = np.zeros(N_gal)
    random.seed(RS)
    while True:
        try:
            colonizer = np.where((abs(galaxy["pos1"]-start_r)<=r_err))[0]
#Can only apply 'subtract' function to dimensionless quantities when other argument is not a quantity (unless the latter is all zero/infinity/nan)
            N_potential = np.size(colonizer)
            RNDind = np.random.random_integers(0, N_potential-1)
            CS[colonizer[RNDind]] = 1
            r_colonizer = galaxy["pos2", colonizer[RNDind]]
            galaxy["col_state"][colonizer[RNDind]] = True
            break
        except ValueError:
            print 'No particle found! Try a larger value for r_err = %f!'%(r_err)
            r_err = float(raw_input('Enter flexibility in r in pc! '))


def initialize_random_colonizer(N_gal, galaxy):
    CS = np.zeros(N_gal)
    random.seed(RS)
    colonizer = int(np.random.uniform(0, N_gal-1))
    CS[colonizer] = 1
    galaxy["col_state"][colonizer] = True

def initialize_galaxy(galaxy):
    initialize_bulge()
    initialize_disk()
    initialize_stellar_halo()

    # POS:Galaxy
    pos1_gal = np.append(pos_r_bulge, pos_rho_disk)
    pos1_gal = np.append(pos1_gal, pos_r_stellar_halo)
    pos2_gal = np.append(pos_theta_bulge, pos_phi_disk)
    pos2_gal = np.append(pos2_gal, pos_theta_stellar_halo)
    pos3_gal = np.append(pos_phi_bulge, pos_z_disk)
    pos3_gal = np.append(pos3_gal, pos_phi_stellar_halo)
    
    vel1_gal = np.append(vel_r_bulge, vel_rho_disk)
    vel1_gal = np.append(vel1_gal, vel_r_stellar_halo)
    vel2_gal = np.append(vel_theta_bulge, vel_phi_disk)
    vel2_gal = np.append(vel2_gal,vel_theta_stellar_halo)
    vel3_gal = np.append(vel_phi_bulge, vel_z_disk)
    vel3_gal = np.append(vel3_gal, vel_phi_stellar_halo)
    
    # This is just to keep the galaxy array consistent, NOT for SB measurement
    I_gal = np.append(I_bulge, I_disk)
    I_gal = np.append(I_gal, I_stellar_halo)
    
    L_initial = np.sum(I_gal)

    galaxy["pos1"] = pos1_gal
    galaxy["pos2"] = pos2_gal
    galaxy["pos3"] = pos3_gal
    galaxy["vel1"] = vel1_gal
    galaxy["vel2"] = vel2_gal
    galaxy["vel3"] = vel3_gal
#    galaxy["I"] = I_gal #FIXME
    if config.random_start:
        initialize_random_colonizer(config.N_galaxy, galaxy)
    else:
        initialize_users_colonizer(config.N_galaxy, galaxy, 
                                   config.start_r_colonizer, 
                                   config.start_r_error_colonizer)

