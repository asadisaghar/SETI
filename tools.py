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
RS = 1
# ===================== #
#        CONSTANTS      #
# ===================== #
pc = 1.
sec = 1.
M_solar = 1.
# Multiply everything by relevant units below to translate them into pc, s, and M_solar
meter = 1.e-3 #km
km = 1./3.086e13  # pc
year = 3.154e7*sec  # s
Myr = 1.e6*year  # s
kg = 1./1.1988435e30  # M_solar
cSpeed = 3.e5*km/sec  # pc/s
colors = ['#d7191c', '#fdae61', '#abd9e9','#2c7bb6']
Gconst = 6.67e-11 #N.m^2.kg^-2 = kg^1.m^1.s^-2.m^2.kg^-2 = kg^-1.m^3.s^-2 
Gconst = Gconst*meter**3/(kg*sec**2)

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
random.seed(RS)

import datetime
import contextlib, time
@contextlib.contextmanager
def timer(msg="XXXXXX"):
    start = datetime.datetime.now()
    yield
    end = datetime.datetime.now()
    print msg, end - start

# Initialization (Uniform or Exponential), no signs. rho, and phi of the disk, and all three(r, phi, theta) of the bulge
def init_pos(N, low, high_scale, dist):
    random.seed(RS)
    if dist == 'exp':
        init = np.random.exponential(high_scale, N)
    elif dist == 'uni':
        init = np.random.uniform(low, high_scale, N)

    return init

# Initialization (Uniform or Exponential), with random signs. z component of the disk
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

def delta_vir_calculator(z_l):
    omega_0 = -1.0
    omega_a = 0.0
    h = 0.7                   # H_0 / 100 [Mpc s km-1]
    H0 = h*100./kpc           # [s-1]
    omega_m = 0.27
    omega_lambda = 0.73
    z_l = 0.5
    omega_m_z = omega_m
    x = omega_m_z -1.
    delta_vir = (18.*pi**2 + 82.*x - 39.*x**2)/omega_m_z
    return delta_vir

#logM/M_solar to sigma(km/s) translation
def mass2sigma(mass, delta_vir, z_l):
    M_solar = 1.98892e30      # solar mass in [kg]
    H = H0*np.sqrt((omega_m*(1.+z_l)**3.+(1.-omega_m)*(1.+z_l)**(3.*(1.+omega_0+omega_a))*np.exp(-3.*omega_a*(z_l/(1.+z_l)))))    # Hubble parameter as a function of redshift (from 10.1103/PhysRevD.83.084045)
    rho_c = 3./(8.*pi*G)*H**2.    # critical mass density of the universe at redshift z_l [kg km-3]
    print rho_c
    Mkg = mass*M_solar
    R_vir = ((3./(4.*pi*delta_vir*rho_c))*mass*M_solar)**(1./3.)        # the radius within which density is 200 times the 
    sigma = np.sqrt(3.*G*Mkg/(2.*R_vir))                      # Velocity dispersion in [km s-1]     
    return sigma

# Initialization (Normal)
def init_norm(loc, scale, size):
    random.seed(RS)
    init = np.random.normal(loc, scale, size)
    return init

# Calculating velocity as a function of r to force the rotation curve (A)
def v_rotational_disk(r, V_opt, R_opt, L2Lstar):
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

def v_rotational_unisphere(r, R, mass):
#        r *= pc2km
#        R *= pc2km
        density = 3.*mass/(4.*pi*(R)**3)
	V_rot = np.sqrt(4./3.*pi*Gconst*density)*r
        print np.max(V_rot)
	return V_rot

# Calculating the oscillatory motion of disk particles [along z direction]
def z_oscillation(galaxy, t, v, N_bulge, N_disk, amp):
    z = galaxy[2,N_bulge:N_bulge+N_disk]
    omega = 2.*pi/(25.*pi*1e7*1e6)
    phase = init_pos(N_disk, 0., 2.*pi, 'uni')
    galaxy[2,N_bulge:N_bulge+N_disk] = amp*np.cos(omega*t+phase)
    return galaxy[2,N_bulge:N_bulge+N_disk]

def r_oscillation(galaxy, t, v, N_bulge, N_disk, amp):
    r = galaxy[0,N_bulge:N_bulge+N_disk]
    omega = 2.*pi/(25.*pi*1e7*1e6)
    phase = init_pos(N_disk, 0., 2.*pi, 'uni')
    galaxy[0,N_bulge:N_bulge+N_disk] = amp*np.cos(omega*t+phase)
    return galaxy[0,N_bulge:N_bulge+N_disk]

def r_bulge_oscillation(galaxy, t, v, N_bulge, amp):
    r = galaxy[0,:N_bulge]
    omega = 2.*pi/(25.*pi*1e7*1e6)
    phase = init_pos(N_bulge, 0., 2.*pi, 'uni')
    galaxy[0,:N_bulge] = amp*np.cos(omega*t+phase)
    return galaxy[0,:N_bulge]

def r_halo_oscillation(galaxy, t, v, N_bulge, N_disk, N_halo, amp):
    r = galaxy[0,N_bulge+N_disk:]
    omega = 2.*pi/(300.*pi*1e7*1e6)
    phase = init_pos(N_halo, 0., 2.*pi, 'uni')
    galaxy[0,N_bulge+N_disk:] = amp*np.cos(omega*t+phase)
    return galaxy[0,N_bulge+N_disk:]

#Initializing disk luminosity distribution
def init_sersic(I_0, alpha, n_s, r):
    L = I_0*np.exp((-r/alpha)**(1./n_s))
    return L

#Where in the galaxy do you think the first colonizing civilization(s) arise?
def CS_manual(N_gal, galaxy, start_r, r_err):
    CS = np.zeros(N_gal)
    random.seed(RS)
    while True:
        try:
            colonizer = np.where((abs(galaxy[0]-start_r)<=r_err))[0]
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

# Particle-to-particle colonization, single probe
# As each step, the *closest* site within the sphere of r=dist is colonized,
# ONLY by the sites which are colonized during the FIRST previous step
def col_sing(galaxy, dist_o, count, coveringFraction, N_bulge, N_disk):
    x_gal = np.zeros_like(galaxy[0])
    y_gal = np.zeros_like(galaxy[1])
    z_gal = np.zeros_like(galaxy[2])
    cs_gal = galaxy[5]
    I_gal = galaxy[4]
    
    # Bulge (Spherical to Cartesian)
    r_gal = galaxy[0,:N_bulge]
    theta_gal = galaxy[1,:N_bulge]
    phi_gal_sph = galaxy[2,:N_bulge]

    x_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[:N_bulge] = r_gal*np.cos(theta_gal)
    
    # Disk (Cylindrical to Cartesian)
    rho_gal = galaxy[0,N_bulge:N_bulge+N_disk]
    phi_gal_cyl = galaxy[1,N_bulge:N_bulge+N_disk]
    z_gal_cyl = galaxy[2,N_bulge:N_bulge+N_disk]
    
    x_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.cos(phi_gal_cyl)
    y_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.sin(phi_gal_cyl)
    z_gal[N_bulge:N_bulge+N_disk] = z_gal_cyl

    # Halo (Spherical to Cartesian)
    r_gal = galaxy[0,N_bulge+N_disk:]
    theta_gal = galaxy[1,N_bulge+N_disk:]
    phi_gal_sph = galaxy[2,N_bulge+N_disk:]

    x_gal[N_bulge+N_disk:] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[N_bulge+N_disk:] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[N_bulge+N_disk:] = r_gal*np.cos(theta_gal)

    # Spot the colonizer
    inds  = np.where(cs_gal==1)[0]
    captured = 0
    for ind in inds:
        dist = dist_o*(count-galaxy[5,ind])
        x_col = x_gal[ind]              
        y_col = y_gal[ind]
        z_col = z_gal[ind]

        # Spot potential colonies (particles that have never been colonized)
        pots = np.where(cs_gal==0)[0] #index
        x_pots = x_gal[pots]
        y_pots = y_gal[pots]
        z_pots = z_gal[pots]    
        d = np.sqrt((x_col-x_pots)**2+(y_col-y_pots)**2+(z_col-z_pots)**2) #distance
        d_min = min(d)
        if d_min <= dist:
            ind_new = pots[np.where(d==d_min)[0][0]]
            galaxy[5,ind_new] = count
            galaxy[4,ind_new] *= (1.-coveringFraction)
            captured += 1

    return galaxy, captured

# Spherical colonization, infinite probes
# As each step, all site within the sphere of r=dist is colonized
def col_inf2(galaxy, dist_o, count, coveringFraction, N_bulge, N_disk):
    x_gal = np.zeros_like(galaxy[0])
    y_gal = np.zeros_like(galaxy[1])
    z_gal = np.zeros_like(galaxy[2])
    cs_gal = galaxy[5]
    I_gal = galaxy[4]
    
    # Bulge (Spherical to Cartesian)
    r_gal = galaxy[0,:N_bulge]
    theta_gal = galaxy[1,:N_bulge]
    phi_gal_sph = galaxy[2,:N_bulge]

    x_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[:N_bulge] = r_gal*np.cos(theta_gal)
    
    # Disk (Cylindrical to Cartesian)
    rho_gal = galaxy[0,N_bulge:N_bulge+N_disk]
    phi_gal_cyl = galaxy[1,N_bulge:N_bulge+N_disk]
    z_gal_cyl = galaxy[2,N_bulge:N_bulge+N_disk]
    
    x_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.cos(phi_gal_cyl)
    y_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.sin(phi_gal_cyl)
    z_gal[N_bulge:N_bulge+N_disk] = z_gal_cyl

    # Halo (Spherical to Cartesian)
    r_gal = galaxy[0,N_bulge+N_disk:]
    theta_gal = galaxy[1,N_bulge+N_disk:]
    phi_gal_sph = galaxy[2,N_bulge+N_disk:]

    x_gal[N_bulge+N_disk:] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[N_bulge+N_disk:] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[N_bulge+N_disk:] = r_gal*np.cos(theta_gal)

    # Spot the colonizer
    inds  = np.where(cs_gal != 0)[0]
    captured = 0
    for ind in inds:
        dist = dist_o*(count-galaxy[5,ind])
        x_col = x_gal[ind]
        y_col = y_gal[ind]
        z_col = z_gal[ind]
        # Spot potential colonies (particles that have never been colonized)
        pots = np.where(cs_gal==0)[0] #index
        x_pots = x_gal[pots]
        y_pots = y_gal[pots]
        z_pots = z_gal[pots]    
        d = np.sqrt((x_col-x_pots)**2+(y_col-y_pots)**2+(z_col-z_pots)**2) #distance
        cols = np.where(d<=dist)[0]
        if len(cols)>0:
            indcs = pots[cols]
            galaxy[5,indcs] = count
            galaxy[4,indcs] *= (1.-coveringFraction)
            captured += len(cols)

    return galaxy, captured

def calculate_reachable(galaxy, dist, ind):
    # spot the colonizer!
    r_col = galaxy[0,ind]
    phi_col = galaxy[1,ind]
    z_col = galaxy[2,ind]
    # mark the reachable space, i.e. dr = dist, dz = dist, dph = arcsin(dist/r_col)
    if dist<=r_col:
        dphi = np.arcsin(dist/r_col)
        reachable = np.where((abs(galaxy[0]-r_col)<=dist) &
                             (abs(galaxy[1]-phi_col)<=dphi) &
                             (abs(galaxy[2]-z_col)<=dist) &
                             (galaxy[5]==0))
    else:
        reachable = np.where((galaxy[0]<=dist) &
                             (abs(galaxy[2]-z_col)<=dist) &
                             (galaxy[5]==0))
    return ind, reachable

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
#    plt.show()
#    plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/%s.png'%(name))

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
#    plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/%s_50.png'%(name))
    figHIST=plt.figure()
    axHIST=figHIST.add_subplot(111)
    HIST=axHIST.hist(log[:,3])
#    plt.show()


def plot_part_galaxy(filename, N_bulge, N_disk, N_halo, mode, txt="", ttxt=""):
    tmp = filename.split('.npy')
    t = int(tmp[0].split('_')[1])
    colonized_fraction = int(tmp[0].split('_')[2])
    colPercent = int(colonized_fraction)
    timeMyr = round(t/1000.,0)
    galaxy = np.load('%s'%(filename))
    
    x_gal = np.zeros_like(galaxy[0])
    y_gal = np.zeros_like(galaxy[1])
    z_gal = np.zeros_like(galaxy[2])
    cs_gal = galaxy[5]
    cont = galaxy[4]
    
    # Bulge (Spherical to Cartesian)
    r_gal = galaxy[0,:N_bulge]
    theta_gal = galaxy[1,:N_bulge]
    phi_gal_sph = galaxy[2,:N_bulge]

    x_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[:N_bulge] = r_gal*np.cos(theta_gal)
    
    # Disk (Cylindrical to Cartesian)
    rho_gal = galaxy[0,N_bulge:N_bulge+N_disk]
    phi_gal_cyl = galaxy[1,N_bulge:N_bulge+N_disk]
    z_gal_cyl = galaxy[2,N_bulge:N_bulge+N_disk]
    
    x_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.cos(phi_gal_cyl)
    y_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.sin(phi_gal_cyl)
    z_gal[N_bulge:N_bulge+N_disk] = z_gal_cyl

    # Halo (Spherical to Cartesian)
    r_gal = galaxy[0,-N_halo:]
    theta_gal = galaxy[1,-N_halo:]
    phi_gal_sph = galaxy[2,-N_halo:]

    x_gal[-N_halo:] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[-N_halo:] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[-N_halo:] = r_gal*np.cos(theta_gal)
    
    # Spot the colonizer
    inds  = np.where(cs_gal!=0)[0]

    x_col = x_gal[inds]              
    y_col = y_gal[inds]
    z_col = z_gal[inds]

#    colonized = np.where(cs_gal!=0)[0]
#    colonized_fraction = len(colonized)/len(cs_gal)*1.

    fig = plt.figure(figsize=(15, 8))
    # Face-on
    axfo = plt.subplot(121)
    axfo.set_aspect('equal')
    cmap = plt.cm.spectral_r
#    cmap.set_bad('w', 0.)
    if mode == 'k':
        fo = axfo.scatter(x_gal/1e3, y_gal/1e3, marker='o', c=(cont), edgecolor='None', alpha=0.1, cmap=cmap, s=30)
        focol = axfo.scatter(x_col/1e3, y_col/1e3, marker='o', c='k', alpha=0.3, s=5)
    elif mode == 'w':
        fo = axfo.scatter(x_gal/1e3, y_gal/1e3, marker='o', c=(cont), edgecolor='None', alpha=0.5, cmap=cmap, s=30)
        focol = axfo.scatter(x_col/1e3, y_col/1e3, marker='o', c='w', edgecolor='None', alpha=0.3, s=5)
    elif mode == 'h':
#        fo = axfo.scatter(x_gal/1e3, y_gal/1e3, marker='o', c=(cont), edgecolor='None', alpha=0.5, cmap=cmap, s=30)
#        focol = axfo.scatter(x_col/1e3, y_col/1e3, marker='o', c='w', edgecolor='None', alpha=0.3, s=5)
        ftst = axfo.scatter(x_gal[N_bulge+1]/1e3, y_gal[N_bulge+1]/1e3, marker='o', c='k', edgecolor='None', alpha=1.0, s=80)
        ftst = axfo.scatter(x_gal[16]/1e3, y_gal[16]/1e3, marker='o', c='b', edgecolor='None', alpha=1.0, s=80)
        ftst = axfo.scatter(x_gal[-1]/1e3, y_gal[-1]/1e3, marker='o', c='r', edgecolor='None', alpha=1.0, s=80)
    
    plt.xlabel(r'X (kpc)')
    plt.ylabel(r'Y (kpc)')
    plt.xlim([-3e1, 3e1])
    plt.ylim([-3e1, 3e1])
#    fig.figtext(-2.7e1, 2.4e1, "t = %4.f Myr"%(timeMyr), fontsize=30)
#    plt.text(-2.7e1, 2.4e1, "%s"%(txt), fontsize=30)

#    axfo.set_title("t = %d Myr"%(timeMyr), fontsize=30)
#    print ("time = %s Myr"%(t/100.))
#    plt.title("time = %s Myr"%(t))
#    cb = plt.colorbar(pad=0.2,
#                      orientation='horizontal')
#    cb.set_label(r'$\mathrm{log(L/L_\odot)}$')
    #Edge-on
    axeo = plt.subplot(122)
    axeo.set_aspect('equal')
    cmap = plt.cm.spectral_r
#    cmap.set_bad('w', 0.)
    if mode == 'k':
        eo = axeo.scatter(x_gal/1e3, z_gal/1e3, marker='o', c=(cont), edgecolor='None', alpha=0.1, cmap=cmap, s=30)
        eocol = axeo.scatter(x_col/1e3, z_col/1e3, marker='o', c='k', alpha=0.3, s=5)
    elif mode == 'w':
        eo = axeo.scatter(x_gal/1e3, z_gal/1e3, marker='o', c=(cont), edgecolor='None', alpha=0.5, cmap=cmap, s=30)
        eocol = axeo.scatter(x_col/1e3, z_col/1e3, marker='o', c='w', edgecolor='None', alpha=0.3, s=5)
    elif mode == 'h':
#        eo = axeo.scatter(x_gal/1e3, z_gal/1e3, marker='o', c=(cont), edgecolor='None', alpha=0.5, cmap=cmap, s=30)
#        eocol = axeo.scatter(x_col/1e3, z_col/1e3, marker='o', c='w', edgecolor='None', alpha=0.3, s=5)
        ftst = axeo.scatter(x_gal[N_bulge+1]/1e3, z_gal[N_bulge+1]/1e3, marker='o', c='k', edgecolor='None', alpha=1.0, s=80)
        ftst = axeo.scatter(x_gal[16]/1e3, z_gal[16]/1e3, marker='o', c='b', edgecolor='None', alpha=1.0, s=80)
        ftst = axeo.scatter(x_gal[-1]/1e3, z_gal[-1]/1e3, marker='o', c='r', edgecolor='None', alpha=1.0, s=80)
    
    plt.xlabel(r'X (kpc)')
    plt.ylabel(r'Z (kpc)')
    plt.xlim([-3e1, 3e1])
    plt.ylim([-3e1, 3e1])
#    cb = plt.colorbar(pad=0.2,
#                      orientation='horizontal')
#    cb.set_label(r'$\mathrm{log(L/L_\odot)}$')

#    axeo.set_title("%d"r'%'" colonized"%(colPercent), fontsize=30)

    plt.text(-9.7e1, 2.4e1, "%4.f"%(timeMyr), fontsize=30, horizontalalignment="left", verticalalignment="center")
    plt.text(-9.0e1, 2.4e1, "   Myr", fontsize=30, horizontalalignment="left", verticalalignment="center")

    plt.text(-2.7e1, 2.4e1, "%3.f"%(colPercent), fontsize=30, horizontalalignment="left", verticalalignment="center")
    plt.text(-1.7e1, 2.4e1, "% colonized", fontsize=30, horizontalalignment="left", verticalalignment="center")

#    plt.figtext(0.15, 0.75, "%4.d Myr"%(timeMyr), fontsize=30)
#    plt.figtext(0.6, 0.75, "%3.d%% colonized"%(colPercent), fontsize=30)
#    axfo.set_title("%4.d Myr"%(timeMyr), fontsize=30)
#    axeo.set_title("%3.d%% colonized"%(colPercent), fontsize=30)
#    plt.text(-1.0e1, 2.4e1, "%s"%(ttxt), fontsize=30)
    plt.text(-8.e1, 3.6e1, "Probe Velocity: 0.01% of light speed", fontsize=30)
    print ("Colonized : %s"%(colPercent))
    print ("time : %s"%(timeMyr))
#    plt.savefig("%s.svg"%(filename))
    plt.savefig("%s.png"%(filename))
#    plt.show()

    return galaxy


#def plot_cont_galaxy(t, x_gal, y_gal, z_gal, cont, R0, I_gal, col_frac=0, bin_no=100):
def plot_cont_galaxy(filename, N_bulge, N_disk, bin_no=100): #If you have enough particle resolution, choose 1000
    tmp = filename.split('.npy')
    t = tmp[0].split('galaxy_')[1]
    galaxy = np.load('%s'%(filename))
    
    x_gal = np.zeros_like(galaxy[0])
    y_gal = np.zeros_like(galaxy[1])
    z_gal = np.zeros_like(galaxy[2])
    cs_gal = galaxy[5]
    cont = galaxy[4]/np.sum(galaxy[4])
    
    # Bulge (Spherical to Cartesian)
    r_gal = galaxy[0,:N_bulge]
    theta_gal = galaxy[1,:N_bulge]
    phi_gal_sph = galaxy[2,:N_bulge]

    x_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[:N_bulge] = r_gal*np.cos(theta_gal)
    
    # Disk (Cylindrical to Cartesian)
    rho_gal = galaxy[0,N_bulge:N_bulge+N_disk]
    phi_gal_cyl = galaxy[1,N_bulge:N_bulge+N_disk]
    z_gal_cyl = galaxy[2,N_bulge:N_bulge+N_disk]
    
    x_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.cos(phi_gal_cyl)
    y_gal[N_bulge:N_bulge+N_disk] = rho_gal*np.sin(phi_gal_cyl)
    z_gal[N_bulge:N_bulge+N_disk] = z_gal_cyl

    # Halo (Spherical to Cartesian)
    r_gal = galaxy[0,N_bulge+N_disk:]
    theta_gal = galaxy[1,N_bulge+N_disk:]
    phi_gal_sph = galaxy[2,N_bulge+N_disk:]

    x_gal[N_bulge+N_disk:] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[N_bulge+N_disk:] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[N_bulge+N_disk:] = r_gal*np.cos(theta_gal)
    
    # Spot the colonizer
    inds  = np.where(cs_gal!=0)[0]

    x_col = x_gal[inds]              
    y_col = y_gal[inds]
    z_col = z_gal[inds]

    colonized_fraction = abs(np.sum(cs_gal)/len(cs_gal))

    from astroML.plotting import setup_text_plots
    setup_text_plots(fontsize=16, usetex=False)
    from astroML.stats import binned_statistic_2d

    fig = plt.figure(figsize=(20, 10))
    # Face-on
    axfo = plt.subplot(121)
    cmap = plt.cm.spectral_r
    cmap.set_bad('w', 0.)
    N, xedges, yedges = binned_statistic_2d(x_gal*1e-3, y_gal*1e-3, cont, 'sum', bins=bin_no)
    plt.imshow((N.T), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               aspect='equal', interpolation='nearest', cmap=cmap)
    plt.xlabel(r'X (kpc)')
    plt.ylabel(r'Y (kpc)')
    plt.xlim([-1e1, 1e1])
    plt.ylim([-1e1, 1e1])
    cb = plt.colorbar(pad=0.2,
                      orientation='horizontal')
    cb.set_label(r'$\mathrm{log(L/L_{total})}$')
    clim_min = min(cont)
    clim_max = max(cont) 
    plt.title("time = %s Myr"%(t))
#    plt.clim(clim_min, clim_max)

    # Edge-on
    axeo = plt.subplot(122)
    cmap = plt.cm.spectral_r
    cmap.set_bad('w', 0.)
    N, xedges, zedges=binned_statistic_2d(x_gal*1e-3, z_gal*1e-3, cont, 'sum', bins=bin_no)
    plt.imshow((N.T), origin='lower',
               extent=[xedges[0], xedges[-1], zedges[0], zedges[-1]],
               aspect='equal', interpolation='nearest', cmap=cmap)

    plt.xlabel(r'X (kpc)')
    plt.ylabel(r'Z (kpc)')
    plt.xlim([-1e1, 1e1])
    plt.ylim([-1e1, 1e1])
    cb = plt.colorbar(pad=0.2,
                      orientation='horizontal')
    cb.set_label(r'$\mathrm{log(L/L_{total})}$')
    clim_min = min(cont)
    clim_max = max(cont) 
#    print ("Colonized fraction = %.2f"%(colonized_fraction))
#    plt.clim(clim_min, clim_max)
#    plt.title("Colonized fraction = %.2f"%(abs(colonized_fraction)))
    plt.savefig("%s.svg"%(filename))
    plt.savefig("%s.png"%(filename))
#    plt.show()
    return galaxy
        

def plot_profile(filename, N_bulge, N_disk):
    tmp = filename.split('.npy')
    t = float(tmp[0].split('galaxy_')[1])
    galaxy = np.load('%s'%(filename))

    radius = galaxy[0]/1e3
    luminosity = galaxy[4]
    col = galaxy[5]
    N_gal = len(col)

    fig = plt.figure(figsize=(10,10))

#    axb = fig.add_subplot(231)
#    axb.plot((radius[:N_bulge]), (luminosity[:N_bulge]), '.y')
#    plt.title("Bulge")

    axd = fig.add_subplot(222)
    axd.plot((radius[N_bulge:N_bulge+N_disk]), np.log10(abs(luminosity[N_bulge:N_bulge+N_disk])), '.c')
    plt.xlim([0, 40])
    plt.ylim([0, 11])
    plt.title("Disk")

#    axh = fig.add_subplot(233)
#    axh.plot((radius[N_bulge+N_disk:]), (luminosity[N_bulge+N_disk:]), '.g')
#    plt.title("Halo")

    axt = fig.add_subplot(223)
    axt.hist((luminosity/1e10), bins=10, color='g', normed=True)
    plt.title("Luminosity")
    plt.xlim([0, 2])
    plt.ylim([0, 2])

    axt = fig.add_subplot(224)
    axt.hist(((luminosity/1e10)), bins=10, color='m', cumulative=True, normed=True)
    plt.title("Cumulative Luminosity")
    plt.xlim([0, 2])
    plt.ylim([0, 2])

    axcol = fig.add_subplot(221)
    axcol.hist(col, bins=2, color='y', normed=True)
    plt.title("Colonization")
    plt.ylim([0, 2])
#    axcolp = fig.add_subplot(236)
#    axcolp.plot(radius, col, '.m')
#    plt.title("Col. Profile")
    
    colonized_fraction = abs(np.sum(col)/N_gal)
    print ("Colonized fraction = %.2f"%(colonized_fraction))

    plt.suptitle("time = %.2f Myr"%(t/1000.))
    plt.savefig("%s.svg"%(filename))
    plt.savefig("%s.png"%(filename))
#    plt.show()                                                                                                                    
    return galaxy

# ==================== #
#        TESTING       #
# ==================== #
# Model functions	
def expo(x, a, b, c):
    return a * np.exp(-b * x) + c

def gauss(x, A, mu, sigma):
    val = A*np.exp(-(x-mu)**2/(2.*sigma**2))
    return val

def particle_dist_tst(xdata, ydata, param1, param2, param3):
    import scipy.optimize as spo
    y = func(xdata, param1, param2, param3)
    popt, pcov = spo.curve_fit(func, xdata, ydata)
    return popt, pcov
    
def number_density_fit(data, bins, model, p1=0, p2=0, p3=0, Pdensity=True, fit=True):
    import scipy.optimize as spo
# Histogramming data
    hist, bin_edges = np.histogram(data, bins, density=Pdensity)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
    if fit:
        print('Fitting in process...')
# The initial guess for function parameters
        p0 = [p1, p2, p3]
        if model == 'gauss':
            coeff, var_matrix = spo.curve_fit(gauss, bin_centers, hist, p0=p0)
# Get the fitted curve
            hist_fit = gauss(bin_centers, coeff[0], coeff[1], coeff[2])
        elif model == 'expo':
            coeff, var_matrix = spo.curve_fit(expo, bin_centers, hist, p0=p0)
            hist_fit = expo(bin_centers, coeff[0], coeff[1], coeff[2])

        plt.plot(bin_centers, hist, 'om', label='Test data')
        plt.plot(bin_centers, hist_fit, '-k', linewidth=2, label='Fitted data')
        
    if fit:
        return (coeff[0], coeff[1], coeff[2])
        

