import healpy as hp
from pys2let import *
import math
import matplotlib.pyplot as plt

##Changing input for S2LET ILC
#scal_maps = np.load('../s2let_ilc/deconv_data/s2let_ilc_dir_para_gauss_wmap_deconv_smoothw_extrapolated_9yr_scal_1024_2_6_3.npy') #WMAP
scal_maps = np.load('/Users/keir/Documents/s2let_ilc_planck/deconv_data/planck_deconv_353_scal_3400_2_6_2.npy') #Planck
#scal_maps = scal_planck*1000 - scal_maps #Residual
wav_maps = np.load('/Users/keir/Documents/s2let_ilc_planck/deconv_data/planck_deconv_353_wav_3400_2_6_2.npy')
#wav_sing_map = np.load('/Users/keir/Documents/s2let_ilc_planck/deconv_data/s2let_ilc_dir_hypatia_planck_deconv_wav_3400_2_6_2_j7_n1.npy')
#R_maps = np.load('/Users/keir/Documents/s2let_ilc/deconv_data/s2let_ilc_dir_para_gauss_wmap_deconv_nosource_smoothw_extrapolated_9yr_wav_1024_2_6_2_j6_n2_Rsingwavscal4.npy')
#R_maps2 = np.load('/Users/keir/Documents/s2let_ilc/deconv_data/s2let_ilc_dir_para_gauss_wmap_deconv_nosource_smoothw_extrapolated_9yr_wav_1024_2_6_2_j6_n1_Rsmoothfatsingwavscal4.npy')
#R_maps_resid = R_maps - R_maps2
ellmax = 3400
scal_lmax = 64
wavparam = 2
ndir = 2
upsample = 0
jmin = 6
jmax = pys2let_j_max(wavparam,ellmax,jmin)
LIM = 0.0002

# Home made plotting routine! inputs : function f (1D array of MW signal), bandlimit L, plot axis ax, and title
def myplot(f, L, ax, plt, title=''):
    thetas, phis = mw_sampling(L) # Compute the theta and phi values of the MW equiangular grid.
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi)) # Convert the input MW 1D array into 2D array with rows for theta and columns for phi. As simple as that!
    cax = ax.imshow(arr.astype(float),interpolation='none',vmin=-1.*LIM,vmax=LIM) # Don't forget to convert to float, otherwise imshow will complain about complex numbers. For polalisation we will have to look into other operations (real part, imag part, angle).
    
    # This is to undersample the grid of x/yticklabels.
    if L > 10:
        step_phi = int(L/4)
        step_theta = int(L/4)
    else:
        step_phi = 3
        step_theta = 3
    
    selec = np.arange(0, nphi, step_phi) # subset of phi tick labels
    ax.set_xticks(selec)
    ax.set_xticklabels(['%.1f' % x for x in phis[selec]])
    ax.set_xlabel(r'$\phi$')
    selec = np.arange(0, ntheta, step_theta) # subset of theta tick labels
    ax.set_yticks(selec)
    ax.set_yticklabels(['%.1f' % x for x in thetas[selec]])
    ax.set_ylabel(r'$\theta$')
    ax.set_title(title)

    #WOULD BE GOOD TO HAVE COLORBARS
    #cbar = plt.colorbar(cax,ticks=[min(f),max(f)],orientation='horizontal')

# Plot scaling function map
fig, ax = plt.subplots(1,1)
myplot(scal_maps, scal_lmax, ax, plt, title='Scaling function map') #NEEDS CORRECTING
'''fig, ax = plt.subplots(1,1)
myplot(f_mw_rec, L, ax, 'Input map converted to MW (reconstructed)')
fig.savefig('test_directional_python_wrappers_1.png')'''


# Plot R/wav_sing map
'''fig, ax = plt.subplots(1,1)
offset, bandlimit, nelem, nelem_wav = wav_ind(7, 1, wavparam, ellmax, ndir, jmin, upsample)
myplot(wav_sing_map, bandlimit, ax, plt, title='j7, n1 map') #NEEDS CORRECTING'''

'''fig, axs = plt.subplots(len(R_maps), len(R_maps[0]), figsize=(4*5, 3*5))
axs = axs.ravel()
for i in xrange(len(R_maps)):
    for j in xrange(len(R_maps[0])):
        offset, bandlimit, nelem, nelem_wav = wav_ind(6, 2, wavparam, ellmax, ndir, jmin, upsample)
        myplot(R_maps2[i,j], bandlimit*2, axs[i*len(R_maps)+j], plt, title='R map [Wav smooth] (j=6,n=1)_' + str(i) + '_' + str(j)) #
        #fig.colorbar(axs[(j-jmin)*ndir+n],ticks=[min(wav_maps[offset:offset+nelem]),max(wav_maps[offset:offset+nelem])])'''
# Pretty adjustment
#fig.subplots_adjust(hspace=0.4, wspace=0.5)
#fig.savefig('test_directional_python_wrappers_2.png')


# Create giant array figure
fig, axs = plt.subplots(jmax-jmin+1, ndir, figsize=(4*ndir, 3*(jmax-jmin)))
axs = axs.ravel()
# Loop through scales j and directions n
for j in range(jmin, jmax+1):
    for n in range(0, ndir):
        # Retreive the boundaries and positions of the right wavelet scale in the giant f_wav array!
        offset, bandlimit, nelem, nelem_wav = wav_ind(j, n, wavparam, ellmax, ndir, jmin, upsample)
        # The right wavelet map corresponding to (j,n) will be f_wav[offset:offset+nelem].
        # It has a band-limit bandlimit
        # nelem_wav is the total number of elements in the j-th scale (i.e., sum of all directions). It's a safety check to verify that we are not forgetting any directions.
        print 'plot id', (j-jmin)*ndir+n, 'j=', j, 'n=', n, 'bounds=',offset, 'to', offset+nelem, 'Total elems:', nelem_wav
        print bandlimit
        # Make the plot!
        myplot(wav_maps[offset:offset+nelem], bandlimit, axs[(j-jmin)*ndir+n],plt,title='Scale '+str(j)+'/'+str(jmax)+', direction '+str(n+1)+'/'+str(ndir))
        #fig.colorbar(axs[(j-jmin)*ndir+n],ticks=[min(wav_maps[offset:offset+nelem]),max(wav_maps[offset:offset+nelem])])

# Pretty adjustment
fig.subplots_adjust(hspace=0.4, wspace=0.5)
#fig.savefig('test_directional_python_wrappers_2.png')
plt.show()
