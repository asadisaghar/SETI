import os
import tools

# in_path = "/home/saas9842/PhD/tmp"
# proj = "slowgalfast2"
# a=in_path+"/%s/galaxy_00_00.npy"%(proj)
# b=in_path+"/%s/galaxy_2199_09.npy"%(proj)
# c=in_path+"/%s/galaxy_2999_53.npy"%(proj)
# d=in_path+"/%s/galaxy_252999_53.npy"%(proj)
# e=in_path+"/%s/galaxy_503000_53.npy"%(proj)
# f=in_path+"/%s/galaxy_1030200_53.npy"%(proj)

in_path = "/home/saas9842/PhD/tmp"
proj = "slowgalmid2"
a=in_path+"/%s/galaxy_4600_00.npy"%(proj)
b=in_path+"/%s/galaxy_19000_13.npy"%(proj)
c=in_path+"/%s/galaxy_23800_50.npy"%(proj)
d=in_path+"/%s/galaxy_254200_50.npy"%(proj)
e=in_path+"/%s/galaxy_523000_50.npy"%(proj)
f=in_path+"/%s/galaxy_1022199_50.npy"%(proj)

# in_path = "/home/saas9842/PhD/tmp"
# proj = "slowgalVslow2"
# a=in_path+"/%s/galaxy_99800_00.npy"%(proj)
# b=in_path+"/%s/galaxy_149800_01.npy"%(proj)
# c=in_path+"/%s/galaxy_199799_05.npy"%(proj)
# d=in_path+"/%s/galaxy_249799_18.npy"%(proj)
# e=in_path+"/%s/galaxy_299799_50.npy"%(proj)
# f=in_path+"/%s/galaxy_899800_50.npy"%(proj)
# g=in_path+"/%s/galaxy_2199800_50.npy"%(proj)

# in_path = "/home/saas9842/PhD/tmp"
# proj = "slowgalslow2"
# a=in_path+"/%s/galaxy_87999_01.npy"%(proj)
# b=in_path+"/%s/galaxy_117400_04.npy"%(proj)
# c=in_path+"/%s/galaxy_127200_06.npy"%(proj)
# d=in_path+"/%s/galaxy_156600_51.npy"%(proj)
# e=in_path+"/%s/galaxy_401599_51.npy"%(proj)
# f=in_path+"/%s/galaxy_656400_51.npy"%(proj)
# g=in_path+"/%s/galaxy_1195400_51.npy"%(proj)

# in_path = "/home/saas9842/PhD/tmp"
# proj = "fastinf2"
# a=in_path+"/%s/galaxy_00_00.npy"%(proj)
# b=in_path+"/%s/galaxy_2199_09.npy"%(proj)
# c=in_path+"/%s/galaxy_2999_54.npy"%(proj)
# d=in_path+"/%s/galaxy_252999_54.npy"%(proj)
# e=in_path+"/%s/galaxy_503000_54.npy"%(proj)
# f=in_path+"/%s/galaxy_1030200_54.npy"%(proj)

# in_path = "/home/saas9842/PhD/tmp"
# proj = "slowinf2"
# a=in_path+"/%s/galaxy_87999_01.npy"%(proj)
# b=in_path+"/%s/galaxy_117400_04.npy"%(proj)
# c=in_path+"/%s/galaxy_127200_07.npy"%(proj)
# d=in_path+"/%s/galaxy_156600_48.npy"%(proj)
# e=in_path+"/%s/galaxy_401599_48.npy"%(proj)
# f=in_path+"/%s/galaxy_656400_48.npy"%(proj)
# g=in_path+"/%s/galaxy_1195400_48.npy"%(proj)

# in_path = "/home/saas9842/PhD/tmp"
# proj = "midinf2"
# a=in_path+"/%s/galaxy_00_00.npy"%(proj)
# b=in_path+"/%s/galaxy_19000_12.npy"%(proj)
# c=in_path+"/%s/galaxy_23800_52.npy"%(proj)
# d=in_path+"/%s/galaxy_273400_52.npy"%(proj)
# e=in_path+"/%s/galaxy_523000_52.npy"%(proj)
# f=in_path+"/%s/galaxy_1022199_52.npy"%(proj)

N_disk = int(2.e4)
N_bulge = int(0.33*N_disk)
N_halo = int(0.01*N_disk)

tools.plot_part_galaxy(a, N_bulge, N_disk, N_halo, mode='k', txt="a)")
tools.plot_part_galaxy(b, N_bulge, N_disk, N_halo, mode='k', txt="b)")
tools.plot_part_galaxy(c, N_bulge, N_disk, N_halo, mode='k', txt="c)")
tools.plot_part_galaxy(d, N_bulge, N_disk, N_halo, mode='k', txt="d)")
tools.plot_part_galaxy(e, N_bulge, N_disk, N_halo, mode='k', txt="e)")
tools.plot_part_galaxy(f, N_bulge, N_disk, N_halo, mode='k', txt="f)")
tools.plot_part_galaxy(g, N_bulge, N_disk, N_halo, mode='k', txt="g)")
