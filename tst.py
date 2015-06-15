import random
import cProfile
from math import *
import numpy as np
import tools
import datetime

# ===================== #
#        FUNCTIONS      #
# ===================== #
def v_rotational_unisphere(r, R, mass):
        r *= pc2km
        R *= pc2km
        density = 3.*mass/(4.*pi*(R)**3)
	V_rot = np.sqrt(4./3.*pi*Gconst*density)*r
        print np.max(V_rot)
	return V_rot

def r_bulge_oscillation(galaxy, t, v, N_bulge, amp):
    r = galaxy[0,:N_bulge]
    omega = 2.*pi/(25.*pi*1e7*1e6)
    phase = init_pos(N_bulge, 0., 2.*pi, 'uni')
    galaxy[0,:N_bulge] = amp*np.cos(omega*t+phase)
    return galaxy[0,:N_bulge]

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
RS = 11
random.seed(RS)

# BULGE
N_bulge = 1e2
R_bulge = 1.5e3  # bulge radius #[pc]
I_0_bulge = 5.e9  # bulge central intensity
alpha_bulge = R_bulge/3.  # bulge scale length #[pc]
n_s_bulge = 5  # Bulge Sersic index(in range: 1.5-10)
M_bulge = 2.e10  # bulge mass #[M_solar]
mean_bulge = 200. # (My arbitrary value!!) #FIXME value  #[km/s]
rho_bulge = M_bulge/(4./3.*pi*R_bulge**3)

# POS:Bulge (Spherical)
r = abs(np.random.normal(0, R_bulge, N_bulge))  #galaxy[0]  #[pc]
theta = tools.init_pos(N_bulge, 0., np.nextafter(pi,4), 'uni')  #galaxy[1]  #[pc]
phi = tools.init_pos(N_bulge, 0., 2.*pi, 'uni')  #galaxy[2]  #[pc]

# VEL:Bulge (Gaussian dist. with given mean and dispersion)
Vr = np.zeros(N_bulge) #[km/s]
Vtheta = np.zeros(N_bulge) #[km/s]
Vphi = np.zeros(N_bulge) #[km/s]

galaxy = np.zeros((6, N_bulge))
galaxy[0] = r
galaxy[1] = theta
galaxy[2] = phi
galaxy[3] = Vr
galaxy[4] = Vtheta
galaxy[5] = Vphi


# # Evaluate bulge velocities (Spherical)
# galaxy[3,:N_bulge] = tools.v_rotational_unisphere(galaxy[0,:N_bulge]*pc2km, R_bulge*pc2km, M_bulge) #[km/s]
# sign = np.round(np.random.uniform(0,1,N_bulge))*2.-1
# galaxy[7,:N_bulge] = sign*np.random.normal(0.0, sigma_bulge, N_bulge) #[km/s]
# # Rotate the bulge
# galaxy[1,:N_bulge] += galaxy[3,:N_bulge]*km2pc*dt_r/galaxy[0,:N_bulge] # v_theta = r*dtheta/dt  #[pc]
# galaxy[2,:N_bulge] += galaxy[7,:N_bulge]*km2pc*dt_r/galaxy[0,:N_bulge]*np.sin(galaxy[2,:N_bulge]) # v_phi = r*sin(theta)*dphi/dt #[pc]
# # Oscillate bulge particles (around original r positions only)
# galaxy[0,:N_bulge] = tools.r_bulge_oscillation(galaxy, t, mean_bulge, N_bulge, R_bulge)   

