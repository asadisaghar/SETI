import math
import numpy as np
import random
import colgal.parameters as pars

class Bulge(object):
    def __init__(self, particle_number, luminosity, effective_radius, mass):
        self.particle_no = particle_number
        self.luminosity = luminosity
        self.r_eff = effective_radius
        self.mass = mass
