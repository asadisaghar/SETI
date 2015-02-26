import math
import numpy as np
import random
import colgal.galaxy as galaxy
import colgal.disks as disk
import colgal.bulge as bulge
import colgal.halos as halo
import colgal.aliens as alien
import colgal.parameters as pars
# DISK

class Galaxy(object):
    def __init__(self, particle_number, total_mass, disk_args={}):
        self.particle_no = particle_number
        self.total_mass = total_mass
        self.components = []
        self.thin_disk = disk.Disk(self, pars.particleNo_thin, pars.h_z_thin, pars.h_rho_thin, pars.mass_thin, pars.I_0_thin, pars.alpha_disk, pars.n_s_disk)
