import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt

#Latex font settings
from matplotlib import rc
from matplotlib import rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 13})

##Changing input for S2LET ILC
scal_maps = np.load('/Users/keir/Documents/s2let_ilc_planck/hybrid_data/planck2015_2_cmb_map_1_scal_300_hybridA_7_1.npy')
wav_fits_root = '/Users/keir/Documents/s2let_ilc_planck/hybrid_data/planck2015_2_cmb_map_1_wav_300_hybridA_7_1_A'
ellmax = 300
scal_lmax = 128
jmin = 7
lamdas = np.array([2,1.4])
l_transitions = np.array([129])
ndir = 1
scal_tiles, wav_tiles, scal_bandlims, wav_bandlims, jmax, l_bounds = construct_hybrid_tiling(ellmax,jmin,lamdas,l_transitions)

# Home made plotting routine! inputs : function f (1D array of MW signal), bandlimit L, plot axis ax, and title
def myplot(f, L, ax, plt, title=''):
    thetas, phis = mw_sampling(L) # Compute the theta and phi values of the MW equiangular grid.
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi)) # Convert the input MW 1D array into 2D array with rows for theta and columns for phi. As simple as that!
    cax = ax.imshow(arr.astype(float),interpolation='nearest',vmin=-1.*LIM,vmax=LIM) # Don't forget to convert to float, otherwise imshow will complain about complex numbers. For polalisation we will have to look into other operations (real part, imag part, angle).
    
    # This is to undersample the grid of x/yticklabels.
    if L > 10:
        step_phi = int(L/3)
        step_theta = int(L/3)
    else:
        step_phi = 3
        step_theta = 3
    
    selec = np.arange(0, nphi, step_phi) # subset of phi tick labels
    ax.set_xticks(selec)
    ax.set_xticklabels([r'$%.1f$' % x for x in phis[selec]])
    ax.set_xlabel(r'$\phi$')
    selec = np.arange(0, ntheta, step_theta) # subset of theta tick labels
    ax.set_yticks(selec)
    ax.set_yticklabels([r'$%.1f$' % x for x in thetas[selec]])
    ax.set_ylabel(r'$\theta$')
    ax.set_title(title)

# Plot scaling function map
fig, ax = plt.subplots(1,1)
LIM = 0.0002
myplot(scal_maps, scal_lmax, ax, r'$\textrm{Scaling function map}$')

# Create giant array figure
fig, axs = plt.subplots(jmax+1, ndir, figsize=(4*ndir, 3*jmax))
axs = axs.ravel()
# Loop through scales j and directions n
offset = 0
LIM = 0.00002
for j in range(jmax+1):
    for n in range(0, ndir):
        #Loading single wavelet maps
        wav_fits = wav_fits_root + str(j) + '_n' + str(n+1) + '.npy'
        wav_maps = np.load(wav_fits)

        bandlim = wav_bandlims[j]
        nelem = bandlim*(2.*bandlim-1.)
        myplot(wav_maps, bandlim, axs[j],plt,title='Scale '+str(j)+'/'+str(jmax)+', direction '+str(n+1)+'/'+str(ndir))
        del wav_maps
        #fig.colorbar(axs[(j-jmin)*ndir+n],ticks=[min(wav_maps[offset:offset+nelem]),max(wav_maps[offset:offset+nelem])])
        offset += nelem

# Pretty adjustment
fig.subplots_adjust(hspace=0.4, wspace=0.5)
#fig.savefig('test_directional_python_wrappers_2.png')
plt.show()


