import math
import numpy as np
import random
import colgal.galaxy as galaxy
import colgal.disks as disk
import colgal.bulge as bulge
import colgal.halos as halo
import colgal.aliens as alien
import colgal.parameters.galaxy_type as pars        
# DISK

class Galaxy(object):
    def __init__(self, particle_number, civil_location, civil_strategy, vessel_speed, galaxy_type="MW", mode="univ_RC"):
        mygalaxy = getattr(parameters.settings, galaxy_type)
        mycivil_location = getattr(parameters.settings, civil_location)
        mycivil_strategy = getattr(parameters.settings, civil_strategy)
        mycalc_mode = getattr(parameters.settings, mode)
        self.particle_no = particle_number
        # set up galaxy components
        self.components = []
        self.thin_disk = disk.Disk(self, mygalaxy.particleNo_thin, mygalaxy.h_z_thin, mygalaxy.h_rho_thin, mygalaxy.mass_thin, mygalaxy.I_0_thin, mygalaxy.alpha_disk, mygalaxy.n_s_disk)
        self.thick_disk = disk.Disk(self, mygalaxy.particleNo_thick, mygalaxy.h_z_thick, mygalaxy.h_rho_thick, mygalaxy.mass_thick, mygalaxy.I_0_thik, mygalaxy.alpha_disk, mygalaxy.n_s_disk)
        self.bulge = bulge.Bulge(self, mygalaxy.particleNo_bulge, mygalaxy.I_0_bulge, mygalaxy.radius_bulge, mygalaxy.mass_bulge)
        self.stellar_halo = halos.Stellar(self. mygalaxy.particleNo_stellar_halo, mygalaxy.radius_stellar_halo, mygalaxy.mass_stellar_halo, mygalaxy.luminosity_stellar_halo)
        self.dark_halo = halos.Dark(self. mygalaxy.particleNo_dark_halo, mygalaxy.radius_dark_halo, mygalaxy.mass_dark_halo, mygalaxy.luminosity_dark_halo)
