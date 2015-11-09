.. image:: http://img.shields.io/badge/arXiv-1508.02406-orange.svg?style=flat
        :target: http://arxiv.org/abs/1508.02406
# SETI
- see https://github.com/koffol/SETI/issues/1 for best suggested parameter combination
- tools.py and vis.py need to be in the same dir as SETI.py (the main body)
- run SETI.py with desired parameters. This dumps the multi-D "galaxy" array into a file at certain frequencies set by either dt_log or colonization_fraction
- vis.py is to visualize the array using different modes of plot_part_gal function

# SETI.py output file formats:
"strategy+starting point"_"kyr from the start"_"colonized percentage"
strategy: inf :spherical wave
          sing: particle-to-particle
starting point: 1: R0 = 0
                2: R0 = 8kpc
Kyr from the start: in static colonization cases, these are chosen to be 
                    ~500Myr and ~1000Myr after the ~50% colonized fraction reached

# modes of colonization:
    Spherical wave (tools.col_inf2):
        active colonization front
        non-active colonization front (only one colonizer)
    Particle-to-particle (tools.col_sing):
        active colonization front (a particle turns into an active colonizer as soon as captured and stays active) 
        non-active colonization front (the actuve colonizing quality is passed on to the newest captured particle)


# SETI
- see https://github.com/koffol/SETI/issues/1 for best suggested parameter combination
- tools.py and vis.py need to be in the same dir as SETI.py (the main body)
- run SETI.py with desired parameters. This dumps the multi-D "galaxy" array into a file at certain frequencies set by either dt_log or colonization_fraction
- vis.py is to visualize the array using different modes of plot_part_gal function

# SETI.py output file formats:
"strategy+starting point"_"kyr from the start"_"colonized percentage"
strategy: inf :spherical wave
          sing: particle-to-particle
starting point: 1: R0 = 0
                2: R0 = 8kpc
Kyr from the start: in static colonization cases, these are chosen to be 
                    ~500Myr and ~1000Myr after the ~50% colonized fraction reached

# modes of colonization:
    Spherical wave (tools.col_inf2):
        active colonization front
        non-active colonization front (only one colonizer)
    Particle-to-particle (tools.col_sing):
        active colonization front (a particle turns into an active colonizer as soon as captured and stays active) 
        non-active colonization front (the actuve colonizing quality is passed on to the newest captured particle)
