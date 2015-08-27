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
                  dtype={"names":["pos1", "pos2", "pos3", "vel1", "vel2", "vel3", "I", "col_state", "x", "y", "z"],
                         "formats":["f4", "f4", "f4", "f4", "f4", "f4", "f4", "bool", "f4", "f4", "f4"]})
tools_and_functions.initialize_galaxy(galaxy)
tools_and_functions.initialize_colonization(galaxy)
simulation_time = config.t_initial
steps = 0
colonized_fraction = 0
while simulation_time < config.t_final:
    if config.galaxy_motion and simulation_time%config.dt_rotation==0:
        tools_and_functions.update_galaxy(galaxy, simulation_time)
    if config.galaxy_colonization and simulation_time%config.dt_colonization==0:
        captured = tools_and_functions.colonize_galaxy(galaxy, steps)
        colonized_fraction = tools_and_functions.measure_colonized_fraction(galaxy)
    if colonized_fraction <= config.max_colonization_fraction:
        if simulation_time%config.dt_log==0:
            tools_and_functions.record_galaxy(galaxy, simulation_time)
    elif colonized_fraction > config.max_colonized_fraction:
        if simulation_time%config.dt_log_stall==0:
            tools_and_functions.record_galaxy(galaxy, simulation_time)
    steps += 1
    
