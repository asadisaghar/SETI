#        CONSTANTS     ##
# ===================== #
yr = pi*1.e7            # s
Myr = yr*1.e6      # s
c_pcyr = 0.306594845    # pc/yr
c_pcs = c_pcyr/yr       # pc/s
c = c_pcs
pc = 3.0857e13          # km
kilo = 1000                # kilo!
G = 4.302e-3/(pc**2)            # pc.M_solar^-1.(pc/s)^2 
# ===================== #

fig=plt.figure()
N_disk = 1e6
N_bulge = int(0.33000*N_disk)
N_gal = N_disk + N_bulge
log = np.loadtxt('inf_4.txt')
ax=fig.add_subplot(111)
#Ls=ax.semilogy(log[:,0]/Myr, log[:,2], 'm', label='$L/L_t$', linewidth=3)
#Ns=ax.semilogy(log[:,0]/Myr, log[:,1]/(N_gal), '-c', label='$N/N_t$', linewidth=3)

#Ns=ax.scatter(log[:,0]/Myr, log[:,3]/(N_gal), marker='o', facecolors='none', edgecolors='c', label='$N/N_t$', linewidth=3)
#Nl=ax.scatter(log[:,0]/Myr, 1. - log[:,3]/(N_gal), marker='o', facecolors='none', edgecolors='m', label='$L/L_t$', linewidth=3)

Ns=ax.plot(log[:,0]/Myr, log[:,3]/(N_gal), '-c', label='$N/N_t$', linewidth=3)
Nl=ax.plot(log[:,0]/Myr, 1. - log[:,3]/(N_gal), '-m', label='$L/L_t$', linewidth=3)


#Nl=ax.semilogy(log[:,0]/Myr, log[:,3]/(N_gal), '-m', label='$L/L_t$', linewidth=3)
#plt.xlim([-0.005,350])
plt.ylim([-0.05,1.05])
xlabel('t (Myr)', fontsize=20)
#ylabel('$N/N_t$', fontsize=20)
plt.legend(loc=10, frameon=False, fontsize=20)

plt.title('Infinite probes\n', fontsize=20)
plt.text(x=1755, y=0.365, s='$V_{probe} = 0.01 c$\n $dt_{const} = 1$ $yr$\n', fontsize=22)

plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/inf_4.png')


