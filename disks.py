import math
import numpy as np
import random
RS = 123124
random.seed(RS)
import colgal.parameters as pars

class Disk(object):
    def __init__(self, particle_number, r_inner=0., scale_height, scale_length, mass, I_0, sersic_alpha, sersic_index):
        self.particle_no = particle_number
        self.r_inner = r_inner
        self.scale_height = scale_height
        self.mass = mass
        self.I_0 = I_0
        self.alpha_sersic = sersic_alpha
        self.n_sersic = sersic_index
# Initialize an exponential disk in a CYLENDRICAL coordinate system
        self.rho = np.random.exponential(self.scale_length, self.particle_no)
        self.phi = np.random.uniform(0., 2*pi, self.particle_no)
        self.height = np.random.exponential(self.scale_height, self.particle_no)
        self.luminosity = self.I_0*np.exp((-rho/self.alpha_sersic)**(1./self.n_sersic))
