import numpy as np
import healpy as hp
import math as mh
import copy as cp
import matplotlib.pyplot as plt
import distinct_colours as dc
import pys2let as ps

#Latex font settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=17.0)

#Increase line thicknesses
plt.rc('axes', linewidth=2.0)
plt.rc('xtick.major', width=2.0)
plt.rc('xtick.minor', width=2.0)
plt.rc('ytick.major', width=2.0)
plt.rc('ytick.minor', width=2.0)
plt.rc('lines', linewidth=1.5)

##Input
L = 128 #255 #500 #3600
J_min = 3 #0
#Bs = np.array([14,2]) #60,2]) #,1.3,1.2])
B = 2
#l_transitions = np.array([15]) #61]) #,513,2017])
J = ps.pys2let_j_max(B,L,J_min)
print J

N = 1 #3
nside = 64 #256
spin = 0
upsample = 1

#Construct valid hybrid tiling
'''hybrid_scal_l, hybrid_wav_l, hybrid_scal_bandlimit, hybrid_wav_bandlimits, J, L_bounds = ps.construct_hybrid_tiling(L,J_min,Bs,l_transitions)
res = ps.verify_tiling(L, hybrid_scal_l, hybrid_wav_l.T.ravel(), hybrid_scal_bandlimit, hybrid_wav_bandlimits)
if res == False:
    raise ValueError('\nAn invalid wavelet tiling has been chosen.\n')
else:
    print '\nA valid wavelet tiling has been chosen.\n'
'''

#Plot array
fname = '/Users/keir/Documents/planck2015_2_cmb_realisations/planck2015_2_cmb_map_1.fits'
#fname = '/Users/keir/Documents/s2let_ilc_planck/COM_CompMap_dust-commander_0256_R2.00.fits'

f_ini = hp.read_map(fname) # Initial map
f_lm = hp.map2alm(f_ini, lmax=L-1) # Its alms
f = hp.alm2map(f_lm, nside=nside, lmax=L-1) # Band limited version

# Convert to MW sampling from spherical harmonics
f_mw = ps.alm2map_mw(ps.lm_hp2lm(f_lm, L), L, spin)

print 'Running analysis_lm2wav'
#f_wav, f_scal = ps.analysis_lm2wav_manualtiling(f_lm, L, N, spin, hybrid_scal_l, hybrid_wav_l.T.ravel(), hybrid_scal_bandlimit, hybrid_wav_bandlimits)
f_wav, f_scal = ps.analysis_lm2wav(f_lm, B, L, J_min, N, spin, upsample)

def mollweide_grid(thetas, phis):
   MAX_ITERATIONS = 1000
   TOL = 1e-10
   thetas = np.pi/2 - thetas
   phis = phis - np.pi
   t = thetas
   for it in xrange(MAX_ITERATIONS):
      dt = (t + np.sin(t) - np.pi*np.sin(thetas)) / (1 + np.cos(t))
      t = t - dt
      if np.max(np.abs(dt)) < TOL:
         break
   t = t/2
   x = 2 * np.sqrt(2) / np.pi * phis * np.cos(t)
   y = np.sqrt(2) * np.sin(t)
   return x, y

# Home made plotting routine! inputs : function f (1D array of MW signal), bandlimit L, plot axis ax, and title
def myplot(f, L, ax, title=''):
    thetas, phis = ps.mw_sampling(L)
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi))
    ax.imshow(arr.astype(float))
    '''if L > 10:
        step_phi = int(L/4)
        step_theta = int(L/4)
    else:
        step_phi = 3
        step_theta = 3'''
    #selec = np.arange(0, nphi, step_phi) # subset of phi tick labels
    #selec = np.arange(0,6,1)
    selec = np.arange(0,7,1) * (nphi / max(phis))
    ax.set_xticks(selec.astype(int))
    ax.set_xticklabels(['%.1f' % x for x in np.arange(0,7,1)]) #phis[selec]])
    ax.set_xlabel(r'$\phi$')
    #selec = np.arange(0, ntheta, step_theta) # subset of theta tick labels
    #selec = np.arange(0,6,1)
    selec = np.arange(0,4,1) * (ntheta / max(thetas))
    ax.set_yticks(selec.astype(int))
    ax.set_yticklabels(['%.1f' % x for x in np.arange(0,4,1)]) #thetas[selec]])
    ax.set_ylabel(r'$\theta$')
    ax.set_title(title)

def myplot_moll(f, L, ax, vmin=None, vmax=None, title='', cmap=None):
    thetas, phis = ps.mw_sampling(L)
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi))
    theta_grid, phi_grid = np.meshgrid(thetas, phis)
    x_grid, y_grid = mollweide_grid(theta_grid, phi_grid)
    if vmin is None and vmax is None :
        cc = np.max(np.abs(arr))
        vmin = -cc
        vmax = cc
    ax.pcolormesh(x_grid, y_grid, arr.T, rasterized=True, vmin=vmin, vmax=vmax) #, cmap=cmap)
    ax.axis('off')

    ax.set_title(title, y=0.975)

    #return im

# Plot equiangular map
'''fig, ax = plt.subplots(1) #,1)
myplot_moll(f_mw.real, L, ax, title=r'Input CMB map') #map converted to MW')
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/cmb_orig_mollweide.pdf',dpi=1200)'''
#plt.show()

#LIM = 300
'''hp.mollview(f_ini,unit=r'$\mu\mathrm{K}$',title=r'Input CMB map',max=1000,norm='log') #,min=-1.*LIM,max=LIM) #*1e6
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/dust_orig_hp.pdf')'''

# Create giant array figure
plt.rc('font', family='serif',size=14.0)
fig, axs = plt.subplots(J-J_min,N,sharex=True,sharey=True,figsize=(8./3.,6.4)) #figsize=(8,6.4) #J+1, N) #, figsize=(4*N, 3*(J)))
axs = axs.ravel()
# Loop through scales j and directions n
offset = 0
for j in xrange(0,4): #J-J_min+1): #0, J+1):
    for n in xrange(0, N):
        bandlimit = L #hybrid_wav_bandlimits[j]
        nelem = bandlimit*(2*bandlimit-1)
        print 'plot id', (j)*N+n, 'j=', j, 'n=', n, 'bounds=',offset, 'to', offset+nelem
        # Make the plot!
        myplot_moll(f_wav[offset:offset+nelem].real, bandlimit, axs[(j)*N+n],
            title=r'Scale %i, direction %i' % (j+1, n+1))
        offset += nelem

# Pretty adjustment
#fig.subplots_adjust(hspace=0.4, wspace=0.5)
fig.subplots_adjust(wspace=0, hspace=0.2, left=0.01, right=0.99, bottom=0.01, top=0.95)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/cmb_decomp_mollweide_n1.pdf',dpi=400)
#plt.show()


