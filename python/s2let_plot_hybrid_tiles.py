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
L = 3999
Bs = [2.0736,1.728,1.44,1.2]
J_min = int(mh.ceil(mh.log(100,B) - 1.)) #Chooses j_min s.t. not less than ell = 100
print "j_min = ", J_min
#J_min = 11

J = ps.pys2let_j_max(B, L, J_min)
scal_ang, wav_ang = ps.axisym_wav_l(B, L, J_min)

#Admissibility
#admiss = (4.*mh.pi)/(2.*np.arange(len(scal_ang))+1.) * (np.square(scal_ang.real) + np.sum(np.square(wav_ang.real),axis=-1)) #Not working!?
admiss = np.square(np.absolute(scal_ang)) + np.sum(np.square(np.absolute(wav_ang)),axis=-1) #Now works

fig, ax = plt.subplots(1,1, figsize=(12,4))
ax.plot(scal_ang,label=r"$j = \mathrm{Scal.}$")
for j in range(J-J_min+1):
    jay = j+J_min
    if j < 6:
        ax.plot(wav_ang[:,j],'-',label=r"$j = %i$" % jay)
        #elif j < 13:
        #ax.plot(wav_ang[:,j],':',label=r"$j = %i$" % jay) #Use dotted lines after colors run out
    elif j < 20:
        ax.plot(wav_ang[:,j],'--',label=r"$j = %i$" % jay) #Use dashes after dots run out
    else:
        ax.plot(wav_ang[:,j],'-.',label=r"$j = %i$" % jay) #Use dash-dot lines after dashes run out
plt.axvline(x=L,color='k',ls='--') #Vertical line at l_max
#plt.axvline(x=3400,color='k',ls='-.') #Vertical line at l = 3400
ax.set_xlim([0,1.25*L])
ax.set_ylim([0,1.2])
ax.set_xlabel(r'$\mathrm{Multipole}$')
#ax.set_title(r'Harmonic response of wavelet kernels: ($\ell_\mathrm{max} = %i$' % L + r',  $\lambda = %.2f$' % B + r',  $j_\mathrm{min} = %i$' % J_min+ r',  $j_\mathrm{max} = %i$)' % J)
ax.legend(prop={'size':16})

plt.savefig("/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/wavelets_harmonic.pdf",bbox_inches='tight')
plt.show()