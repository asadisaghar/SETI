import os
#here = os.system('pwd')
#os.chdir('/home/saas9842/mn/PhD/numerical_project')
import tools
#os.chdir(here)
probe = "fastsing"
loc = 2
fileformat = "png"

# in_path+"inf2/npy/galaxy_00_00.npy"
# in_path+"inf2/npy/galaxy_1780_13.npy"
# in_path+"/npy/galaxy_1840_53.npy"
# in_path+"/npy/galaxy_2040_73.npy"
# in_path+"/npy/galaxy_100000_75.npy"
# in_path+"/npy/galaxy_400000_75.npy"
# in_path+"/npy/galaxy_800000_75.npy"

#in_path+"/npy/galaxy_1199980_55.npy"

in_path = "/home/saas9842/PhD/tmp/%s%d"%(probe, loc)
out_path = "/home/saas9842/PhD/tmp/%s%d"%(probe, loc)
#in_path = "/home/saas9842/PhD/numerical_project/%s%d/npy"%(probe, loc)
#out_path = "/home/saas9842/PhD/numerical_project/%s%d"%(probe, loc)
#out_path = "/home/saas9842/PhD/numerical_project/re_referee1/low_res2/low_col2/%s%d"%(probe, loc)

#in_path = "/home/saas9842/PhD/numerical_project/re_referee1/low_res2/low_col2/mn_%s%d"%(probe, loc)
#out_path = "/home/saas9842/PhD/numerical_project/re_referee1/low_res2/low_col2/mn_%s%d"%(probe, loc)

out_path_k = out_path+"/"+fileformat+"_k"
#out_path_w = out_path+"/"+fileformat+"_w"
#out_path_h = out_path+"/"+fileformat+"_h"

#out_path = "/home/saas9842/PhD/numerical_project/re_referee1/low_res2/low_col/%s%d"%(probe, loc)
os.system('ls %s/galaxy*.npy > %s/npy.lst'%(in_path, in_path))
galist = open("%s/npy.lst"%(in_path), "r")
#os.system('mkdir %s/png_c'%(out_path))
os.system('mkdir %s/npy'%(out_path))
#os.system('mkdir %s'%(out_path_w))
os.system('mkdir %s'%(out_path_k))
#os.system('mkdir %s'%(out_path_h))
#os.system('mkdir %s/png_p'%(out_path))

N_disk = int(2.e4)  # number of particles in disk (same order as N_gal)
N_bulge = int(0.33*N_disk)
N_halo = int(0.01*N_disk)
for columns in ( raw.strip().split() for raw in galist ):
    print columns[-1]
#    tools.plot_cont_galaxy(columns[-1], N_bulge, N_disk)
#    os.system('mv %s/*.png %s/png_c'%(in_path, out_path))
#    tools.plot_part_galaxy(columns[-1], N_bulge, N_disk, N_halo, mode='w', txt="tst")
#    os.system('mv %s/*.%s %s'%(in_path, fileformat, out_path_w))
    tools.plot_part_galaxy(columns[-1], N_bulge, N_disk, N_halo, mode='k', txt="tst")
    os.system('mv %s/*.%s %s'%(in_path, fileformat, out_path_k))
#    tools.plot_part_galaxy(columns[-1], N_bulge, N_disk, N_halo, mode='h', txt="a)")
#    os.system('mv %s/*.%s %s'%(in_path, fileformat, out_path_h))
#    tools.plot_profile(columns[-1], N_bulge, N_disk)
#    os.system('mv %s/*.png %s/png_p'%(in_path, out_path))
#    os.system('mv %s %s/npy'%(columns[-1], out_path))

