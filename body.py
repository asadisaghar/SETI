import config
import tools_and_functions
from math import *
import numpy as np
import cProfile
import datetime
import contextlib, time
@contextlib.contextmanager
def timer(msg="XXXXXX"):
    start = datetime.datetime.now()
    yield
    end = datetime.datetime.now()
    print msg, end - start

galaxy = np.zeros(config.N_galaxy,
                  dtype={"names":["pos1", "pos2", "pos3", "vel1", "vel2", "vel3", "I", "col_state"],
                         "formats":["f4", "f4", "f4", "f4", "f4", "f4", "f4", "bool"]})
tools_and_functions.initialize_galaxy(galaxy)
tools_and_functions.initialize_colonizing(galaxy)
