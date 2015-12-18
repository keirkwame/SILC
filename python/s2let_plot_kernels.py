execfile('/Users/keir/Documents/s2let_ilc/s2let_ilc_code/python/CBcm.py')

import sys
sys.path.append('/Users/bl/Dropbox/Voids/flaglet-env/src/flagletvoidfinder/flagletvoidfinder')

import matplotlib.ticker
import matplotlib.cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from copy import copy
import numpy as np
import scipy.interpolate
from pys2let import *
#from pyflaglet import *
#from flagtools import *

#Latex font settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=14.0)

#Increase line thicknesses
plt.rc('axes', linewidth=2.0)
plt.rc('xtick.major', width=2.0)
plt.rc('xtick.minor', width=2.0)
plt.rc('ytick.major', width=2.0)
plt.rc('ytick.minor', width=2.0)
plt.rc('lines', linewidth=1.5)

# Deviation around zero colormap (blue--red)
cols = []
for x in np.linspace(0,1, 256):
	rcol = 0.237 - 2.13*x + 26.92*x**2 - 65.5*x**3 + 63.5*x**4 - 22.36*x**5
	gcol = ((0.572 + 1.524*x - 1.811*x**2)/(1 - 0.291*x + 0.1574*x**2))**2
	bcol = 1/(1.579 - 4.03*x + 12.92*x**2 - 31.4*x**3 + 48.6*x**4 - 23.36*x**5)
	cols.append((rcol, gcol, bcol))

cm_plusmin = matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_plusmin", cols)

def mollweide_grid(thetas, phis):
    MAX_ITERATIONS = 1000;
    TOL = 1e-10;
    thetas = np.pi/2 - thetas;
    phis = phis - np.pi
    t = thetas
    for it in range(MAX_ITERATIONS):
       dt = (t + np.sin(t) - np.pi*np.sin(thetas)) / (1 + np.cos(t));
       t = t - dt;
       if(np.max(np.abs(dt)) < TOL):
          break
    t = t/2
    x = 2 * np.sqrt(2) / np.pi * phis * np.cos(t)
    y = np.sqrt(2) * np.sin(t)
    return x, y

def myplot(f, L, ax, vmin=None, vmax=None, title='', cmap=None):
    thetas, phis = mw_sampling(L) 
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi)).astype(float)
    theta_grid, phi_grid = np.meshgrid(thetas, phis)
    if vmin is None and vmax is None :
        cc = np.max(np.abs(arr))
        vmin = -cc
        vmax = cc
    im = ax.pcolormesh(theta_grid, phi_grid, arr.T, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.axis('off')
    return im

'''def myplot_moll(f, L, ax, vmin=None, vmax=None, title='', cmap=None):
    thetas, phis = mw_sampling(L) 
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi)).astype(float)
    theta_grid, phi_grid = np.meshgrid(thetas, phis)
    x_grid, y_grid = mollweide_grid(theta_grid, phi_grid)
    if vmin is None and vmax is None :
        cc = np.max(np.abs(arr))
        vmin = -cc
        vmax = cc
    im = ax.pcolormesh(x_grid, y_grid, arr.T, cmap=cmap, rasterized=True) #vmin=vmin, vmax=vmax
    ax.axis('off')
    
    ax.set_title(title, y=0.975)
    
    return im'''

def myplot_moll(f, L, ax, vmin=None, vmax=None, title='', cmap=None):
    thetas, phis = mw_sampling(L)
    ntheta = len(thetas)
    nphi = len(phis)
    arr = f.reshape((ntheta, nphi)).astype(float)
    theta_grid, phi_grid = np.meshgrid(thetas, phis)
    x_grid, y_grid = mollweide_grid(theta_grid, phi_grid)
    if vmin is None and vmax is None :
        cc = np.max(np.abs(arr))
        vmin = -cc
        vmax = cc
    ax.pcolormesh(x_grid, y_grid, arr.T, rasterized=True, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.axis('off')

    ax.set_title(title, y=0.975)

def wignerfct2(el, em, en, alpha, beta, gamma, dl_beta):
    return dl_beta[el, em+L-1, en+L-1] * np.exp(-1j*em*alpha) * np.exp(-1j*en*gamma)

def rotate(flm, L, alpha, beta, gamma, dl_beta):
    flm_rec = copy(flm)
    for el in range(L):
        fac2 = [flm[el*el+el+en] 
                for en in range(-el, el+1)]
        for em in range(-el, el+1):
            fac1 = [wignerfct2(el, em, en, alpha, beta, gamma, dl_beta) 
                    for en in range(-el, el+1)]
            flm_rec[el*el+el+em] = np.dot(fac1, fac2)
    return flm_rec

##Input
L = 128
B_l = 2
N = 1 #3
spin = 0
J_min_l = 3
J_l = pys2let_j_max(B_l, L, J_min_l)

# Get directional tiling
scal_l, wav_lm = wavelet_tiling(B_l, L, N, J_min_l, spin)

# Need to rotate some
alpha, beta, gamma = (np.pi, np.pi/2, 0)
dl_beta = ssht_dl_beta_risbo(beta, L)

# Now plot
#cmap = plt.get_cmap(cm_plusmin)
#cmap = CB2cm['BurP']
cmap = 'coolwarm'
fig, axs = plt.subplots(J_l-J_min_l,N, sharex=True, sharey=True, figsize=(8./3.,6.4)) #figsize=(8,6.4)
#cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
axs = axs.ravel()
'''j = 0
n = 0
gamma = -1.*np.pi*(2./3.)
f = alm2map_mw(wav_lm[:,j].ravel(), L, spin).reshape((L,2*L-1))
wav_lm_rot = rotate(wav_lm[:,j], L, alpha, beta, gamma, dl_beta)
f_rot = alm2map_mw(wav_lm_rot, L, spin).reshape((L,2*L-1))'''
for j in xrange(0,4): #J_l-J_min_l+1):
    #n = 0
    f = alm2map_mw(wav_lm[:,j].ravel(), L, spin).reshape((L,2*L-1))
    '''wav_lm_rot = rotate(wav_lm[:,j], L, alpha, beta, gamma, dl_beta)
    f_rot = alm2map_mw(wav_lm_rot, L, spin).reshape((L,2*L-1))'''
    for n in xrange(0, N):
        print j, n
        #f = alm2map_mw(wav_lm[:,j].ravel(), L, spin).reshape((L,2*L-1))
        
        gamma = -1.*n*np.pi*(2./3.)
        
        wav_lm_rot = rotate(wav_lm[:,j], L, alpha, beta, gamma, dl_beta)
        f_rot = alm2map_mw(wav_lm_rot, L, spin).reshape((L,2*L-1))
        '''myplot(f.real, L, axs[5*j+0], cmap=cmap, title='Real part of unrotated kernel')
        myplot(f.imag, L, axs[5*j+1], cmap=cmap, title='Imag part of rotated kernel')
        myplot(f_rot.real, L, axs[5*j+2], cmap=cmap, title='Real part of unrotated kernel')
        myplot(f_rot.imag, L, axs[5*j+3], cmap=cmap, title='Imag part of rotated kernel')'''
        myplot_moll(f_rot.real, L, axs[j*N+n], title=r'Scale %i, direction %i' % (j+1, n+1), cmap=cmap)

fig.subplots_adjust(wspace=0, hspace=0.2, left=0.01, right=0.99, bottom=0.01, top=0.95)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/wavelets_spatial_mollweide_n1.pdf',dpi=400)


