import numpy as np
import math as mh
import matplotlib.pyplot as plt
import pys2let as ps

#Latex font settings
from matplotlib import rc
from matplotlib import rcParams
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rcParams.update({'font.size': 16})

##Input
L = 3600
J_min = 0
Bs = np.array([60,2,1.3,1.2])
l_transitions = np.array([61,513,2017])

pltscal = 1.25

#Construct valid hybrid tiling
hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits, J, L_bounds = ps.construct_hybrid_tiling(L,J_min,Bs,l_transitions)
res = ps.verify_tiling(L, hybrid_scal_l, hybrid_wav_l.T.ravel(), hybrid_scal_bandlimit, hybrid_wav_bandlimits)
if res == False:
    raise ValueError('\nAn invalid wavelet tiling has been chosen.\n')
else:
    print '\nA valid wavelet tiling has been chosen.\n'

'''fig, ax = plt.subplots(1,1, figsize=(12,4))
ax.plot(scal_tiles,label=r"$j = \mathrm{Scal.}$")
for j in range(jmax+1):
    if j < 6:
        ax.plot(wav_tiles[:,j],'-',label=r"$j = %i$" % j)
        #elif j < 13:
        #ax.plot(wav_tiles[:,j],':',label=r"$j = %i$" % j) #Use dotted lines after colors run out
    elif j < 20:
        ax.plot(wav_tiles[:,j],'--',label=r"$j = %i$" % j) #Use dashes after dots run out
    else:
        ax.plot(wav_tiles[:,j],'-.',label=r"$j = %i$" % j) #Use dash-dot lines after dashes run out
plt.axvline(x=ellmax,color='k',ls='--') #Vertical line at l_max
#plt.axvline(x=3400,color='k',ls='-.') #Vertical line at l = 3400
ax.set_xlim([0,1.25*ellmax])
ax.set_ylim([0,1.2])
ax.set_xlabel(r'$\mathrm{Multipole}$')
#ax.set_title(r'Harmonic response of wavelet kernels: ($\ell_\mathrm{max} = %i$' % L + r',  $\lambda = %.2f$' % B + r',  $j_\mathrm{min} = %i$' % J_min+ r',  $j_\mathrm{max} = %i$)' % J)
ax.legend(prop={'size':16})

plt.savefig("/Users/keir/Documents/s2let_ilc_planck/hybrid_wavelets_harmonic.pdf",bbox_inches='tight')'''

fig, axs = plt.subplots(Bs.size+2, 1, figsize=(8, 14))
axs = axs.ravel()
J_mins = np.zeros((Bs.size,), dtype=np.int32)
J_mins[0] = J_min
for k in range(Bs.size):
    scal_l, wav_l = ps.axisym_wav_l(Bs[k], L, J_mins[k])
    Jt = ps.pys2let_j_max(Bs[k], L, 0)
    axs[k].plot(scal_l, lw=2, ls='dashed')
    for j in range(0,Jt+1-J_mins[k]):
        #if j < 7:
        axs[k].plot(wav_l[:,j], lw=2, ls='solid')
        #else:
        #axs[k].plot(wav_l[:,j], lw=2, ls='-.')
    for kk in range(Bs.size):
        axs[k].axvline(L_bounds[kk], c='k', ls='-.')
    axs[k].axvline(x=L,color='k',ls='-.') #Vertical line at l_max
    axs[k].set_xscale('log')
    axs[k].set_ylim([0, 1.2])
    axs[k].set_xlim([1, pltscal*L])
    axs[k].set_title('Tiling %i'%(k+1)+' of %i'%Bs.size+', with B=%.1f'%Bs[k]+', Jmin=0, and defined at %i'%L_bounds[k]+'-el-%i'%L_bounds[k+1])
for k in [Bs.size, Bs.size+1]:
    axs[k].plot(hybrid_scal_l, lw=2, ls=':',label=r"$j = \mathrm{Scal.}$")
    for j in range(0,J+1):
        if j < 7:
            axs[k].plot(hybrid_wav_l[:,j], lw=2, ls='solid',label=r"$j = %i$" % j)
        else:
            axs[k].plot(hybrid_wav_l[:,j], lw=2, ls='--',label=r"$j = %i$" % j)
    for kk in range(Bs.size):
        axs[k].axvline(L_bounds[kk], c='k', ls='-.')
    axs[k].axvline(x=L,color='k',ls='-.') #Vertical line at l_max
    axs[k].set_ylim([0, 1.2])
    axs[k].set_xlim([1, pltscal*L])
    if k == Bs.size:
        axs[k].set_title('Hybrid wavelet tiling (log scale)')
        axs[k].set_xscale('log')
        #axs[k].legend(prop={'size':16},loc='upper left')
    else:
        axs[k].set_xlim([0, pltscal*L])
        axs[k].set_xlabel(r'$\mathrm{Multipole}$')
        axs[k].set_title('Hybrid wavelet tiling (linear scale)')
        axs[k].legend(prop={'size':10})
plt.tight_layout()

plt.savefig("/Users/keir/Documents/s2let_ilc_planck/hybrid_wavelets_harmonic_Ctest.pdf",bbox_inches='tight')

plt.show()


