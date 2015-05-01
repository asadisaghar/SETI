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

meter = 1.e-3 #km
km = 1./3.086e13  # pc
year = 3.154e7*sec  # s
Myr = 1.e6*year  # s
kg = 1./1.1988435e30  # M_solar
cSpeed = 3.e5*km/sec  # pc/s
colors = ['#d7191c', '#fdae61', '#abd9e9','#2c7bb6']
Gconst = 6.67e-11 #N.m^2.kg^-2 = kg^1.m^1.s^-2.m^2.kg^-2 = kg^-1.m^3.s^-2 
Gconst = Gconst*meter**3/(kg*sec**2)
# ===================== #
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

# Calculating radial velocity od the disk as a function of r from arXiv:1104.3236v4 (eq. 12) 
def v_radial(r, R_c):
    return 1. - np.exp(-r/R_c)

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

# Spherical colonization, infinite probes
# As each step, all site within the sphere of r=dist is colonized
def col_inf(galaxy, dist, coveringFraction, N_bulge):
    x_gal = np.zeros_like(galaxy[0])
    y_gal = np.zeros_like(galaxy[1])
    z_gal = np.zeros_like(galaxy[2])

    # Bulge (Spherical to Cartesian)
    r_gal = galaxy[0,:N_bulge]
    theta_gal = galaxy[1,:N_bulge]
    phi_gal_sph = galaxy[2,:N_bulge]

    x_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.cos(phi_gal_sph)
    y_gal[:N_bulge] = r_gal*np.sin(theta_gal)*np.sin(phi_gal_sph)
    z_gal[:N_bulge] = r_gal*np.cos(theta_gal)
    
    # Disk (Cylindrical to Cartesian)
    rho_gal = galaxy[0,N_bulge:]
    phi_gal_cyl = galaxy[1,N_bulge:]
    z_gal_cyl = galaxy[2,N_bulge:]
    
    x_gal[N_bulge:] = rho_gal*np.cos(phi_gal_cyl)
    y_gal[N_bulge:] = rho_gal*np.sin(phi_gal_cyl)
    z_gal[N_bulge:] = z_gal_cyl

    # Spot the colonizer
    ind  = np.where(galaxy[5]==1)[0][0]
    z_col = z_gal[ind]
    y_col = y_gal[ind]
    x_col = x_gal[ind]              

    d = np.sqrt((x_col-x_gal)**2+(y_col-y_gal)**2+(z_col-z_gal)**2)
    d_reachable = np.where(d<=dist)[0]
    if len(d_reachable) > 10:
        indcs = d_reachable[np.random.random_integers(0, len(d_reachable)-1, 10)]
    else:
        indcs = d_reachable
#    print d_reachable
#    if galaxy[5,d_reachable].any()==0:
    galaxy[5,indcs]=1

    galaxy[4,indcs] *= (1.-coveringFraction)
    galaxy[5,ind] = -1

#    distances = d[d_reachable]

    print "%d new sites are captured!"%(len(indcs))
#    else:
#        print "all reachable site are already taken"

    return galaxy

# Particle-to-particle colonization, single probe
# As each step, the *closest* site within the sphere of r=dist is colonized,
# ONLY by the sites which are colonized during the FIRST previous step
def col_single(galaxy, galaxy_cart, dist, count, coveringFraction):
    dist *= count
    print "Our current reachable diatance: %.2f pc"%(dist)
    col_dist = 0
    colonized = 0
    ind_dmin = 0
    ind, reachable = calculate_reachable(galaxy, dist, ind=0)
    N_reachable = np.size(reachable[0])
    if N_reachable>0:
        print "%d potential sites!"%(N_reachable)
        print "Searching for the closest site..."
        ##inside update_colonization##
        r_col = galaxy[0,ind]
        phi_col = galaxy[1,ind]
        z_col = galaxy[2,ind]
        x_col = r_col*cos(phi_col)
        y_col = r_col*sin(phi_col)
        d2 = np.power(galaxy_cart[0,reachable]-x_col,2)
        d2 += np.power(galaxy_cart[1,reachable]-y_col,2)
        d2 += np.power(galaxy[2,reachable]-z_col,2)
        d2min = d2[d2 !=0].min()
        dmin = np.sqrt(d2min) 
        ind_dmin = np.where(d2==d2min)[1]
        print "site no. %d will soon be ours!"%(ind_dmin)
        galaxy[5,ind] = -1.
        galaxy[5,ind_dmin] = 1.
        galaxy[4,ind_dmin] *= (1. - coveringFraction)
##
#            galaxy, dmin , ind_dmin = update_colonization(galaxy, dist, ind,
#                                               reachable, coveringFraction)
        colonized += 1
        col_dist += dmin
#        print "We have gone so far: %.2f pc"%(col_dist)
        count = 1
    else:
#        print "No potential sites! Try again in the next time step!!"
        count +=1

    return galaxy[4], galaxy[5], colonized, count, ind_dmin


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


def plot_part_galaxy(filename):
    tmp = filename.split('.npy')
    t = tmp[0].split('_')[1]
    galaxy = np.load('%s'%(filename))
    x_gal = galaxy[0]*np.cos(galaxy[1])
    y_gal = galaxy[0]*np.sin(galaxy[1])
    z_gal = galaxy[2]
    
    cont = galaxy[4]
    cs = galaxy[5]
    colonized_fraction = np.sum(galaxy[5])/len(galaxy[5])

    x_col = galaxy[0,np.where(galaxy[5]==1)[0]]*np.cos(galaxy[1,np.where(galaxy[5]==1)[0]])
    y_col = galaxy[0,np.where(galaxy[5]==1)[0]]*np.sin(galaxy[1,np.where(galaxy[5]==1)[0]])
    z_col = galaxy[2,np.where(galaxy[5]==1)[0]]

    fig = plt.figure(figsize=(20, 10))
    # Face-on
    axfo = plt.subplot(121)
    cmap = plt.cm.jet
    cmap.set_bad('k', 1.)
#    fo = axfo.scatter(x_gal, y_gal, marker='o', c=np.log10(cont), alpha=0.3, cmap=cmap)
    fo = axfo.scatter(x_col, y_col, marker='o', c='k', alpha=0.3)
    plt.xlabel(r'X (pc)')
    plt.ylabel(r'Y (pc)')
    plt.xlim([-5e4, 5e4])
    plt.ylim([-5e4, 5e4])
    plt.title("time = %s Myr"%(t))
#    cb = plt.colorbar(pad=0.2,
#                      orientation='horizontal')
#    cb.set_label(r'$\mathrm{log(L/L_\odot)}$')
    #Edge-on
    axfo = plt.subplot(122)
    cmap = plt.cm.jet
    cmap.set_bad('k', 1.)
#    eo = axfo.scatter(x_gal, z_gal, marker='o', c=np.log10(cont), alpha=0.3, cmap=cmap)
    eo = axfo.scatter(x_col, z_col, marker='o', c='k', alpha=0.3)
    plt.xlabel(r'X (pc)')
    plt.ylabel(r'Z (pc)')
    plt.xlim([-5e4, 5e4])
    plt.ylim([-5e4, 5e4])
#    cb = plt.colorbar(pad=0.2,
#                      orientation='horizontal')
#    cb.set_label(r'$\mathrm{log(L/L_\odot)}$')
    plt.title("Colonized fraction = %.2f"%(colonized_fraction))
    plt.savefig("%s.png"%(filename))
#    plt.show()

    return galaxy


#def plot_cont_galaxy(t, x_gal, y_gal, z_gal, cont, R0, I_gal, col_frac=0, bin_no=100):
def plot_cont_galaxy(filename, bin_no=100): #If you have enough particle resolution, choose 1000
    tmp = filename.split('.npy')
    t = tmp[0].split('_')[1]
    galaxy = np.load('%s'%(filename))
    x_gal = galaxy[0]*np.cos(galaxy[1])
    y_gal = galaxy[0]*np.sin(galaxy[1])
    z_gal = galaxy[2]
    cont = galaxy[4]
    colonized_fraction = np.sum(galaxy[5])/len(galaxy[5])

    from astroML.plotting import setup_text_plots
    setup_text_plots(fontsize=16, usetex=False)
    from astroML.stats import binned_statistic_2d

    fig = plt.figure(figsize=(20, 10))
    # Face-on
    axfo = plt.subplot(121)
    cmap = plt.cm.spectral_r
    cmap.set_bad('w', 1.)
    N, xedges, yedges = binned_statistic_2d(x_gal, y_gal, cont, 'mean', bins=bin_no)
    plt.imshow(np.log10(N.T), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               aspect='equal', interpolation='nearest', cmap=cmap)
    plt.xlabel(r'X (pc)')
    plt.ylabel(r'Y (pc)')
    plt.xlim([-5e4, 5e4])
    plt.ylim([-5e4, 5e4])
    cb = plt.colorbar(pad=0.2,
                      orientation='horizontal')
    cb.set_label(r'$\mathrm{log(L/L_\odot)}$')
    clim_min = np.min(np.log10(cont))
    clim_max = np.max(np.log10(cont)) 
    plt.title("time = %s Myr"%(t))
#    plt.clim(clim_min, clim_max)

    # Edge-on
    axeo = plt.subplot(122)
    cmap = plt.cm.spectral_r
    cmap.set_bad('w', 1.)
    N, xedges, zedges=binned_statistic_2d(x_gal, z_gal, cont, 'mean', bins=bin_no)
    plt.imshow(np.log10(N.T), origin='lower',
               extent=[xedges[0], xedges[-1], zedges[0], zedges[-1]],
               aspect='equal', interpolation='nearest', cmap=cmap)

    plt.xlabel(r'X (pc)')
    plt.ylabel(r'Z (pc)')
    plt.xlim([-1e4, 1e4])
    plt.ylim([-1e4, 1e4])
    cb = plt.colorbar(pad=0.2,
                      orientation='horizontal')
    cb.set_label(r'$\mathrm{log(L/L_\odot)}$')
    clim_min = np.min(np.log10(cont))
    clim_max = np.max(np.log10(cont)) 
#    plt.clim(clim_min, clim_max)
    plt.title("Colonized fraction = %.2f"%(colonized_fraction))
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
        
