import tools
import matplotlib.pyplot as plt

N_disk = int(5.e2)  # number of particles in disk (same order as N_gal)
N_bulge = int(0.33*N_disk)

gi=tools.plot_part_galaxy('galaxy_00.npy')
gf=tools.plot_part_galaxy('galaxy_95.npy')

fig=plt.figure()
ax1=fig.add_subplot(231)
ax1.hist(gi[0,:N_bulge], 100, label='initial')
ax1.hist(gf[0,:N_bulge], 100, alpha=0.2, label='final')
ax1.set_title('Bulge r')

ax2=fig.add_subplot(232)
ax2.hist(gi[0,N_bulge:N_bulge+N_disk], 100)
ax2.hist(gf[0,N_bulge:N_bulge+N_disk], 100, alpha=0.2)
ax2.set_title('Disk r')

ax3=fig.add_subplot(233)
ax3.hist(gi[2,N_bulge:N_bulge+N_disk], 100)
ax3.hist(gf[2,N_bulge:N_bulge+N_disk], 100, alpha=0.2)
ax3.set_title('Disk z')

ax4=fig.add_subplot(234)
ax4.plot(gi[0,N_bulge:N_bulge+N_disk], gi[1,N_bulge:N_bulge+N_disk], 'oy')
ax4.plot(gf[0,N_bulge:N_bulge+N_disk], gf[1,N_bulge:N_bulge+N_disk], 'oc', alpha=0.2)
ax4.set_title('Disk r-phi')

ax5=fig.add_subplot(235)
ax5.plot(gi[2,N_bulge:N_bulge+N_disk], gi[1,N_bulge:N_bulge+N_disk], 'oy')
ax5.plot(gf[2,N_bulge:N_bulge+N_disk], gf[1,N_bulge:N_bulge+N_disk], 'oc', alpha=0.2)
ax5.set_title('Disk z_phi')

ax6=fig.add_subplot(236)
ax6.plot(gi[0,N_bulge:N_bulge+N_disk], gi[2,N_bulge:N_bulge+N_disk], 'oy')
ax6.plot(gf[0,N_bulge:N_bulge+N_disk], gf[2,N_bulge:N_bulge+N_disk], 'oc', alpha=0.2)
ax6.set_title('Disk r-z')

#plt.legend()
plt.show()
