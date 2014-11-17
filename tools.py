from math import *
import numpy as np
import scipy as sp
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.rc("font", size=16, family='serif', weight='normal')
#plt.rc("axes", labelsize=16, titlesize=16)
#plt.rc("xtick", labelsize=16)
#plt.rc("ytick", labelsize=16)
#plt.rc("legend", fontsize=16)
import sys

#          CONSTANTS       ##
# ========================= #
yr = pi*1.e7            # s
Myr = yr*1.e6      # s
c_pcyr = 0.306594845    # pc/yr
c_pcs = c_pcyr/yr       # pc/s
c = c_pcs
pc = 3.0857e13          # km
kilo = 1000                # kilo!
G = 4.302e-3/(pc**2)            # pc.M_solar^-1.(pc/s)^2 
np.random.seed(10)
almost_black = '#262626'    
# ========================= #

# Initialization (Uniform or Exponential), no signs. rho, and phi of the disk, and all three(r, phi, theta) of the bulge
def    init_pos(N, low, high_scale, dist):
    random.seed(10)
    if dist == 'exp':
        init = np.random.exponential(high_scale, N)
    elif dist == 'uni':
        init = np.random.uniform(low, high_scale, N)

    return init

# Initialization (Uniform or Exponential), with random signs. z component of the disk
def    init_z(N, low, high_scale, dist):
    random.seed(10)
    if dist == 'exp':
        init = np.random.exponential(high_scale, N)
    elif dist == 'uni':
        init = np.random.uniform(low, high_scale, N)

    sign = np.random.random_integers(low=0, high=N-1, size=N/2)
    for i in sign:
        init[i] *= -1
                    
    return init

# Initialization (Normal)
def    init_norm(loc, scale, size):
    random.seed(10)
    init = np.random.normal(loc, scale, size)
    return init

# Calculating velocity as a function of r to force the rotation curve (A)
def v_rotational(r, V_opt, R_opt, L2Lstar):
    V2_opt = np.power(V_opt, 2)
    x = r/R_opt
    a = 1.5*np.power(L2Lstar, 1.5)
    beta = 0.72+0.44*np.log10(L2Lstar)
    # Disk component
    V2_disk = V2_opt*beta*1.97*np.power(x, 1.22)/np.power((x**2+0.78**2), 1.43)
    # DM component
    V2_DM = V2_opt*(1. - beta)*(1. + a**2)*x**2/(x**2 + a**2)

    V2_rot = V2_disk + V2_DM
    return np.sqrt(V2_rot)


#Initializing disk luminosity distribution
def init_sersic(I_0, alpha, n_s, r):
    L = I_0*np.exp((-r/alpha)**(1./n_s))
    return L

#Where in the galaxy do you think the first colonizing civilization(s) arise? 
def CS_random(N_gal):
    CS = np.zeros(N_gal)
    random.seed(10)
    colonizer = int(np.random.uniform(0, N_gal-1))
    CS[colonizer] = 1
    return CS

# Spherical colonization, infinite probes
# As each step, the *closest* site within the sphere of r=dist is colonized
def col_inf(galaxy, dist, count, ind):
    dist *= count
    col_dist = 0
    colonized = 0
    r_col = galaxy[0,ind]
    phi_col = galaxy[1,ind]
    z_col = galaxy[2,ind]
#    print dist
#    print r_col
    if dist<=r_col:
        dphi = np.arcsin(dist/r_col)
        reachable = np.where((abs(galaxy[0,:]-r_col)<=dist) &
                             (abs(galaxy[1,:]-phi_col)<=dphi) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    else:
#            sys.exit('OUCH! Send slower probes or check on them at shorter time steps!')
        reachable = np.where((galaxy[0,:]<=dist) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    N_reachable = np.size(reachable[0])
#    print N_reachable
#    print dist
#    print col_dist
    if N_reachable>0:
        while col_dist < dist/2.:
#            galaxy = update_colonization_inf(galaxy, dist, ind,
#                                             reachable)
            galaxy[5, ind] = -1.
            galaxy[5, reachable] = 1.
            galaxy[4, reachable] = 0.
            colonized += N_reachable
            col_dist += dist
            print '%d new sites!'%(colonized)
    else:
        print 'Dead end! %d-%f'%(count, dist)
        count +=1

    return galaxy[4,:], galaxy[5,:], colonized, count

#def update_colonization_inf(galaxy, dist, ind, reachable):
#        r_col = galaxy[0,ind]
#        phi_col = galaxy[1,ind]
#        z_col = galaxy[2,ind]
#        galaxy[5,ind] = -1.
#        galaxy[5,reachable] = 1.
#        galaxy[4, reachable] = 0
#        return galaxy

# Particle-toparticle colonization, single probe
# As each step, the *closest* site within the sphere of r=dist is colonized,
# ONLY by the sites which are colonized during the FIRST previous step
def col_single(galaxy, dist, count):
    dist *= count
    col_dist = 0
    colonized = 0
#    print 'distance = %e pc'%(dist)
    ind, reachable = calculate_reachable(galaxy, dist, ind=0)
    N_reachable = np.size(reachable[0])
#    print '%d Reachable sites exist! %d-%f'%(N_reachable, count, dist)
    if N_reachable>0:
        while col_dist < dist/2.:
            galaxy, dmin = update_colonization(galaxy, dist, ind,
                                               reachable)
            
            colonized += 1
#            print '%d newly colonized'%(colonized)
#            print '%f pc to go!'%(dist-col_dist)
            col_dist += dmin
            count = 1
    else:
        print 'Dead end! %d-%f'%(count, dist)
        count +=1

    return galaxy[4,:], galaxy[5,:], colonized, count


def calculate_reachable(galaxy, dist, ind):
    # spot the colonizer!
#    ind = np.where(galaxy[5,:]==1)[0]
#    print 'Mission started at # %d: '%(ind[0])
    r_col = galaxy[0,ind]
    phi_col = galaxy[1,ind]
    z_col = galaxy[2,ind]
#    x_col = r_col*cos(phi_col)
#    y_col = r_col*sin(phi_col)
    # mark the reachable sphere, i.e. dr = dist, dz = dist, dph = arcsin(dist/r_col)
    if dist<=r_col:
        dphi = np.arcsin(dist/r_col)
        reachable = np.where((abs(galaxy[0,:]-r_col)<=dist) &
                             (abs(galaxy[1,:]-phi_col)<=dphi) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    else:
#            sys.exit('OUCH! Send slower probes or check on them at shorter time steps!')
        reachable = np.where((galaxy[0,:]<=dist) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    return ind, reachable

def update_colonization(galaxy, dist, ind, reachable):
        r_col = galaxy[0,ind]
        phi_col = galaxy[1,ind]
        z_col = galaxy[2,ind]
        x_col = r_col*cos(phi_col)
        y_col = r_col*sin(phi_col)
        galaxy_cart2 = np.zeros((2, np.size(galaxy)))
        galaxy_cart2[0,reachable] = galaxy[0,reachable]*np.cos(galaxy[1,reachable])
        galaxy_cart2[1,reachable] = galaxy[0,reachable]*np.sin(galaxy[1,reachable])
        d2 = np.zeros((1, np.size(galaxy)))
        d2[0,reachable] = np.power(galaxy_cart2[0,reachable]-x_col,2)
        d2[0,reachable] += np.power(galaxy_cart2[1,reachable]-y_col,2)
        d2[0,reachable] += np.power(galaxy[2,reachable]-z_col,2)
#        print 'size(d2):%d'%(np.size(d2))
        d2min = d2[d2 !=0].min()
#        print d2min
        dmin = np.sqrt(d2min)
        ind_dmin = np.where(d2==d2min)[1]
        galaxy[5,ind] = -1.
        galaxy[5,ind_dmin] = 1.
        galaxy[4, ind_dmin] = 0
#        d2min = np.min(d2[np.nonzero(d2)])
#        print'Mission accomplished at # %d: '%(ind_dmin[0])
#        print '(d_min, r_reachable) = (%f, %f)'%(dmin, dist)
        return galaxy, dmin


