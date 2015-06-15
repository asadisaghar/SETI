import datetime
import contextlib, time
import numpy as np
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
dt = 0.1*Myr
#dt_const = np.logspace(-12, 1, num=1)*Myr  # construction time delay
dt_const = 1e-13*Myr
#VC = np.logspace(2, 0, num=1)*cSpeed  # probe velocity
VC = 1e-1*cSpeed
dist = VC * (dt-dt_const)

@contextlib.contextmanager
def timer(msg="XXXXXX"):
    start = datetime.datetime.now()
    yield
    end = datetime.datetime.now()
    print msg, end - start

galaxy = np.load('galaxy.npy')
ind  = np.where(galaxy[5]==1)[0][0]
r_gal = galaxy[0]
phi_gal = galaxy[1]

x_gal = galaxy[0]*np.cos(galaxy[1])
y_gal = galaxy[0]*np.sin(galaxy[1])
z_gal = galaxy[2]

r_col = galaxy[0,ind]
phi_col = galaxy[1,ind]
z_col = galaxy[2,ind]
y_col = y_gal[ind]
x_col = x_gal[ind]

dr = dist # particles with r in [r-dr, r+dr] are reachable
dphi = np.arcsin(dist/r_col) # particles with phi in [phi-dphi, phi+dphi] are reachable

with timer("----------np.where------------"):

    r_reach = np.where((abs(r_gal-r_col)<=dr))[0]
    print len(r_reach)
    phi_reach = np.where((abs(phi_gal-phi_col)<=dphi))[0]
    print len(phi_reach)
    z_reach = np.where((abs(z_gal-z_col)<=dr))[0]
    print len(z_reach)
    CS_reach = np.where(galaxy[5]==0)[0]
    print len(CS_reach)

with timer("----------np.where&------------"):
    reachable = np.where(
                         (abs(galaxy[0]-r_col)<=dist) &
                         (abs(galaxy[1]-phi_col)<=dphi) &
                         (abs(galaxy[2]-z_col)<=dist) &
                         (galaxy[5]==0)
                         )[0]
#    print reachable
    print len(reachable)
   
    
with timer("------------distance-----------"):
    d = np.sqrt((x_col-x_gal)**2+(y_col-y_gal)**2+(z_col-z_gal)**2)
    d_reachable = np.where(d<=dist)[0]
    distances = d[d_reachable]
#    print distances
    print len(d_reachable)
