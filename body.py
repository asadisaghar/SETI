import config


# Galaxy 11-D array
galaxy = np.zeros(N_gal, dtype={'names':['pos1', 'pos2', 'pos3', 'vel1', 'vel2', 'vel3', 'I', 'CS'],
                                'formats':['f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'bool']})
# CS:Galaxy (Colonization Status)
if RandomStart:
    CS_gal = tools.CS_random(N_gal)
else:
    CS_gal, r_colonizer = tools.CS_manual(N_gal, galaxy, start_r, r_err)



