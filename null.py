import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", size=20, family='serif', weight='normal')
plt.rc("axes", labelsize=16, titlesize=20)
plt.rc("xtick", labelsize=16)
plt.rc("ytick", labelsize=16)
plt.rc("legend", fontsize=16)

ns = [1, 137, 1359]
    # Number of independent measurements, i.e. individual galaxies w/o KIII
refs = ['Single galaxy', 'Annis 99', 'Zackrisson+14']
#colors = ['#fdcc8a', '#fc8d59', '#d7301f']
colors = ['#d7191c', '#fdae61', '#abd9e9','#2c7bb6']
i = 0
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
for n in ns:
    T = 5.5e9*n
#    T = 1.e10*n
        # the time we integgrate over, i.e. the amount of time we think our
        # Universe had the potential to harbor a KIII, but it hasen't!!
    #P = [0.01, 0.1, 0.25, 0.5, 0.7]
    P = np.linspace(0.01, 0.9, 100)
        # the probability of null detection being due
        # to random chance, or in our case,
        # indistinguishable KIII features from physical
    r = -1/T*np.log(P)  # the upper limit on the arising rate of KIIIs in the Univ.
    t = (1/r)*1.e-6
        # A lower limit on the timescale for the arise of a KIII
    lns = ax.semilogy(P, t, label='%s (n=%d)'%(refs[i],n), c=colors[i], linewidth=3)
#    plt.errorbar(P, np.log10(t), yerr=0.1, xerr=None, fmt='-', capsize=3,
#                barsabove=False, lolims=False, uplims=True)
    i += 1

ax.set_xlabel('$P_{null}$')
ax.set_ylabel('$t_{rise}$(Myr)')
plt.legend(loc=2, frameon=False)
ax.tick_params(size=5, which='major', axis='y')
plt.show()
plt.savefig('/home/saas9842/Dropbox/SETI_report/Figs4/null.jpg')
