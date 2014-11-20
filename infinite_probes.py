from math import *
import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", size=20, family='serif', weight='normal')
plt.rc("axes", labelsize=16, titlesize=20)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)

#          CONSTANTS       ##
# ========================= #
yr = pi*1.e7            # s
Myr = yr*1.e6      # s
c_pcyr = 0.306594845    # pc/yr
c_pcs = c_pcyr/yr       # pc/s
c = c_pcs
pc = 3.0857e13          # km
kilo = 1000                # kilo!
G = 4.302e-3/(pc**2)            # pc.M_solar^-1.(pc/s)^2 
colors = ['#d7191c', '#fdae61', '#abd9e9','#2c7bb6']

fig=plt.figure()
N_disk = 1e6
N_bulge = int(0.33000*N_disk)
N_gal = N_disk + N_bulge
log = np.loadtxt('bulge_cf1_inf.txt')
ax=fig.add_subplot(111)


#Ns=ax.plot(log[:,0]/Myr, log[:,3]/(N_gal), c=colors[0], label='$N/N_t$', linewidth=3)
#Nl=ax.plot(log[:,0]/Myr, 1. - log[:,3]/(N_gal), c=colors[1], label='$L/L_t$', linewidth=3)

Ns=ax.scatter(log[:,0]/Myr, log[:,3]/N_gal, marker='o', facecolors='none', edgecolors=colors[0], label='$N/N_t$', linewidth=3)
Nl=ax.scatter(log[:,0]/Myr, log[:,2]/Li, marker='o', facecolors='none', edgecolors=colors[1], label='$L/L_t$', linewidth=3)
#Nl=ax.scatter(log[:,0]/Myr, 1. - log[:,3]/(N_gal), marker='o', facecolors='none', edgecolors=colors[1], label='$L/L_t$', linewidth=3)

#Nl=ax.semilogy(log[:,0]/Myr, log[:,3]/(N_gal), '-m', label='$L/L_t$', linewidth=3)
plt.xlim([-50,2100])
plt.ylim([-0.05,1.05])
plt.xlabel('t (Myr)', fontsize=20)
#ylabel('$N/N_t$', fontsize=20)
plt.legend(loc=6, frameon=False, fontsize=20)

plt.title('Infinite probes\n', fontsize=20)
plt.text(x=1055, y=0.365, s='$V_{probe} = 0.01 c$\n $dt_{const} = 1$ $yr$\n', fontsize=22)

plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/bulge_cf1_inf.png')
