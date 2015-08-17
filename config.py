from math import *
import numpy as np
import astropy.units as u
import astropy.constants as cst
import random
RS = 200
random.seed(RS)

import cProfile
import datetime
import contextlib, time
@contextlib.contextmanager
def timer(msg="XXXXXX"):
    start = datetime.datetime.now()
    yield
    end = datetime.datetime.now()
    print msg, end - start
# ===================== #
#        HANDLES        #
# ===================== #
N_disk = int(2.e4)
# time step to dump the array into a file before the maximum colonization is achieved  (enter the number in Myr)
dt_log =  1.*u.Myr.to(u.second)
#time step to dump the array into a file after the maximum colonization is achieved  (enter the number in Myr)
dt_log_stall = 20.*dt_log
dt_rotation = 0.2*u.Myr.to(u.second)
dt_colonization = 0.2*u.Myr.to(u.second)
dt_construction = 1e-12*u.Myr.to(u.second)
v_probe = 1e-4*cst.c.to(u.km/u.second)
t_initial = 0*u.second
t_final = 15.e2*u.Myr.to(u.second)
 
colonization = True
disk_rotation = True
disk_oscillation_r = True
disk_oscillation_z = True

SingleProbe = False
InfiniteProbe = not(SingleProbe)

probe = "animatedgal"
if probe == 'inf' or probe == 'sing':
    col_f = 0.75
elif probe == 'sinf' or probe == 'ssing':
    col_f = 0.5
else:
    col_f = 0.75

coveringFraction = 1.0
RandomStart = False
if not RandomStart:
# Change below only if RandomStart = False
    start_r = 8.e3*u.parsec  # radial distance from the galactic center
#    loc = 2
    start_r_error = 10.*u.parsec
