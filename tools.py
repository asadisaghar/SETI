import random
import cProfile
from math import *
import numpy as np
import scipy as sp
import sys
import matplotlib.pyplot as plt
plt.rc("font", size=20, family='serif', weight='normal')
plt.rc("axes", labelsize=16, titlesize=20)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)
RS = 123124
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
colors = ['#d7191c', '#fdae61', '#abd9e9','#2c7bb6']
# ===================== #
# Initialization (Uniform or Exponential), no signs. rho, and phi of the disk, and all three(r, phi, theta) of the bulge
def    init_pos(N, low, high_scale, dist):
    random.seed(RS)
    if dist == 'exp':
        init = np.random.exponential(high_scale, N)
    elif dist == 'uni':
        init = np.random.uniform(low, high_scale, N)

    return init

# Initialization (Uniform or Exponential), with random signs. z component of the disk
def    init_z(N, low, high_scale, dist):
    random.seed(RS)
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
    random.seed(RS)
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
def CS_manual(N_gal, galaxy, start_r, r_err):
    CS = np.zeros(N_gal)
    random.seed(121323)
    while True:
        try:
            colonizer = np.where((abs(galaxy[0,:]-start_r)<=r_err))[0]
            N_potential = np.size(colonizer)
            RNDind = np.random.random_integers(0, N_potential-1)
            CS[colonizer[RNDind]] = 1
            r_colonizer = galaxy[0, colonizer[RNDind]]
            break
        except ValueError:
            print 'No particle found! Try a larger value for r_err = %f!'%(r_err)
            r_err = float(raw_input('Enter flexibility in r in pc! '))
    return CS, r_colonizer

def CS_random(N_gal):
    CS = np.zeros(N_gal)
    random.seed(RS)
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
    if dist<=r_col:
        dphi = np.arcsin(dist/r_col)
        reachable = np.where((abs(galaxy[0,:]-r_col)<=dist) &
                             (abs(galaxy[1,:]-phi_col)<=dphi) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    else:
        reachable = np.where((galaxy[0,:]<=dist) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    N_reachable = np.size(reachable[0])
    if N_reachable>0:
        while col_dist < dist/2.:
            galaxy[5, ind] = -1.
            galaxy[5, reachable] = 1.
            galaxy[4, reachable] = 0.
            colonized += N_reachable
            col_dist += dist

    else:
        count +=1

    return galaxy[4,:], galaxy[5,:], colonized, count

# Particle-toparticle colonization, single probe
# As each step, the *closest* site within the sphere of r=dist is colonized,
# ONLY by the sites which are colonized during the FIRST previous step
def col_single(galaxy, dist, count, coveringFraction):
    dist *= count
    col_dist = 0
    colonized = 0
    ind_dmin = 0
    ind, reachable = calculate_reachable(galaxy, dist, ind=0)
    N_reachable = np.size(reachable[0])
    if N_reachable>0:
        while col_dist < dist/2.:
            galaxy, dmin , ind_dmin = update_colonization(galaxy, dist, ind,
                                               reachable, coveringFraction)
            
            colonized += 1
            col_dist += dmin
            count = 1
    else:
        count +=1

    return galaxy[4,:], galaxy[5,:], colonized, count, ind_dmin


def calculate_reachable(galaxy, dist, ind):
    # spot the colonizer!
    r_col = galaxy[0,ind]
    phi_col = galaxy[1,ind]
    z_col = galaxy[2,ind]
    # mark the reachable space, i.e. dr = dist, dz = dist, dph = arcsin(dist/r_col)
    if dist<=r_col:
        dphi = np.arcsin(dist/r_col)
        reachable = np.where((abs(galaxy[0,:]-r_col)<=dist) &
                             (abs(galaxy[1,:]-phi_col)<=dphi) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    else:
        reachable = np.where((galaxy[0,:]<=dist) &
                             (abs(galaxy[2,:]-z_col)<=dist) &
                             (galaxy[5,:]==0))
    return ind, reachable

def update_colonization(galaxy, dist, ind, reachable, coveringFraction):
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
        d2min = d2[d2 !=0].min()
        # Add try!!
#        while True:
#            try:
#                d2min = d2[d2 !=0].min()
#                dmin = np.sqrt(d2min) 
#                break
#            except ValueError:
#                print 'No particle found! Try a larger value for r_err = %f!'%(r_err)
#                r_err = float(raw_input('Enter flexibility in r in pc! '))

        dmin = np.sqrt(d2min) 
        ind_dmin = np.where(d2==d2min)[1]
        galaxy[5,ind] = -1.
        galaxy[5,ind_dmin] = 1.
        galaxy[4, ind_dmin] *= (1. - coveringFraction) 
        return galaxy, dmin, ind_dmin

# ==================== #
#        PLOTTING      #
# ==================== #
def infplot(name, N_gal, Li, r_colonizer, VC, dt_const):
    fig=plt.figure()
    log = np.loadtxt('%s.txt'%(name))
    ax=fig.add_subplot(111)
    Ns=ax.plot(log[:,0], log[:,2], '-', c=colors[0], label='$N/N_t$', linewidth=3)
    Nl=ax.plot(log[:,0], log[:,1], '-', c=colors[1], label='$L/L_t$', linewidth=3)
#    Ns=ax.scatter(log[:,0], log[:,2], marker='o', facecolors='none', edgecolors=colors[0], label='$N/N_t$', linewidth=2)
#    Nl=ax.scatter(log[:,0], log[:,1], marker='o', facecolors='none', edgecolors=colors[1], label='$L/L_t$', linewidth=2)
    plt.xlim([-log[-1,0]/10.,log[-1,0]*1.1])
    plt.ylim([-0.05,1.05])
    plt.xlabel('t (Myr)', fontsize=20)
    #ylabel('$N/N_t$', fontsize=20)
#    plt.legend(loc=6, frameon=False, fontsize=20)

    plt.title('R$_i$ = %.2f kpc'%(r_colonizer*1.e-3), fontsize=20)
#    plt.text(0.5, 0.5, s='$V_{probe} = %.2f c$\n $dt_{const} = %.1f$ $yr$\n'%(VC/cSpeed, dt_const/year), fontsize=22)
    plt.show()
    plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/%s.png'%(name))

def singleplot(name, N_gal, Li, r_colonizer, VC, dt_const):
    fig=plt.figure()
    log = np.loadtxt('%s.txt'%(name))
    ax=fig.add_subplot(111)
#    Ns=ax.plot(log[:,0], log[:,2], '-', c=colors[0], label='$N/N_t$', linewidth=3)
#    Nl=ax.plot(log[:,0], log[:,1], '-', c=colors[1], label='$L/L_t$', linewidth=3)
    Ns=ax.scatter(log[:,0], log[:,2], marker='o', facecolors='none', edgecolors=colors[0], label='$N/N_t$', linewidth=2)
    Nl=ax.scatter(log[:,0], log[:,1], marker='o', facecolors='none', edgecolors=colors[1], label='$L/L_t$', linewidth=2)
    plt.xlim([-log[-1,0]/10.,log[-1,0]*1.1])
    plt.ylim([-0.05,1.05])
    plt.xlabel('t (Myr)', fontsize=20)
    #ylabel('$N/N_t$', fontsize=20)
#    plt.legend(loc=6, frameon=False, fontsize=20)

    plt.title('R$_i$ = %.2f kpc'%(r_colonizer*1.e-3), fontsize=20)
#    plt.text(0.5, 0.5, s='$V_{probe} = %.2f c$\n $dt_{const} = %.1e$ $Myr$\n'%(VC/cSpeed, dt_const/Myr), fontsize=22)
    plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/%s_50.png'%(name))
    figHIST=plt.figure()
    axHIST=figHIST.add_subplot(111)
    HIST=axHIST.hist(log[:,3])
    plt.show()
