import numpy as np
import healpy as hp
import math as mh
import copy as cp
import matplotlib.pyplot as plt
import distinct_colours as dc
import pys2let as ps

# Linear colormap (rainbow) - len = 256
cols = [(0,0,0)]
for x in np.linspace(0,1, 254):
    rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
    gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
    bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
    cols.append((rcol, gcol, bcol))
cols.append((1,1,1))

step = len(cols) / 12 #13
cols = cols[1:-1:step]

#Latex font settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=20.0)

#Increase line thicknesses
plt.rc('axes', linewidth=2.0)
plt.rc('xtick.major', width=2.0)
plt.rc('xtick.minor', width=2.0)
plt.rc('ytick.major', width=2.0)
plt.rc('ytick.minor', width=2.0)
plt.rc('lines', linewidth=2.0)

##Input
L = 3600 #3600
J_min = 6 #0
Bs = np.array([2,1.3,1.2])
l_transitions = np.array([513,2017])

pltscal = 1.3 #1.35 for j=Scal.

#Construct valid hybrid tiling
hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits, J, L_bounds = ps.construct_hybrid_tiling(L,J_min,Bs,l_transitions)
res = ps.verify_tiling(L, hybrid_scal_l, hybrid_wav_l.T.ravel(), hybrid_scal_bandlimit, hybrid_wav_bandlimits)
if res == False:
    raise ValueError('\nAn invalid wavelet tiling has been chosen.\n')
else:
    print '\nA valid wavelet tiling has been chosen.\n'

#cols = dc.get_distinct(7)

fig, ax = plt.subplots(1) #,1, figsize=(12,4))
ax.plot(hybrid_scal_l,color=cols[0],label=r"$\mathrm{Scal.}$")
for j in range(J+1):
    if j < 6:
        ax.plot(hybrid_wav_l[:,j],color=cols[j+1],label=r"$j = %i$" % j)
        #elif j < 13:
        #ax.plot(wav_tiles[:,j],':',label=r"$j = %i$" % j) #Use dotted lines after colors run out
    elif j < 20:
        ax.plot(hybrid_wav_l[:,j],color=cols[j+1],label=r"$j = %i$" % j) #Use dashes after dots run out
    else:
        ax.plot(wav_tiles[:,j],'-.',label=r"$j = %i$" % j) #Use dash-dot lines after dashes run out
ax.axvline(x=3400,ls='--',color='k') #Vertical line at l = 3400
ax.axvline(x=L,ls='--',color='k') #Vertical line at l_max
ax.set_xlim([0,pltscal*L])
ax.set_ylim([0,1.1])
ax.set_xlabel(r'Multipole $\ell$', labelpad = -1) #, size = 30.) #Increased size for talk
#ax.set_title(r'Harmonic response of wavelet kernels: ($\ell_\mathrm{max} = %i$' % L + r',  $\lambda = %.2f$' % B + r',  $j_\mathrm{min} = %i$' % J_min+ r',  $j_\mathrm{max} = %i$)' % J)
ax.legend(prop={'size':17.5},frameon=False) #16.5 for paper #16. for talk

fig.subplots_adjust(left=0.05,right=0.99) #,bottom = 0.125) #Bottom for talk

plt.savefig("/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/hybrid_wavelets_mnras5.pdf") #,bbox_inches='tight')

'''fig, axs = plt.subplots(Bs.size+2, 1, figsize=(8, 14))
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

plt.savefig("/Users/keir/Documents/s2let_ilc_planck/hybrid_wavelets_harmonic_Ctest.pdf",bbox_inches='tight')'''

#plt.show()


