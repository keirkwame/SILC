import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
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
L = 139
N = 1
nrows = 1

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

    #ax.set_title(title, y=0.975)

# Plot equiangular map
'''fig, ax = plt.subplots(1) #,1)
myplot_moll(f_mw.real, L, ax, title=r'Input CMB map') #map converted to MW')
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/cmb_orig_mollweide.pdf',dpi=1200)'''
#plt.show()

'''LIM = 150
hp.mollzoom(f_ini[0]*1e6,unit=r'$\mu\mathrm{K}$',title=r'Input CMB Stokes $Q$ map',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_pol/fg_353_q.pdf')
hp.mollzoom(f_ini[1]*1e6,unit=r'$\mu\mathrm{K}$',title=r'Input CMB Stokes $U$ map',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_pol/fg_353_u.pdf')'''

# Create giant array figure
plt.rc('font', family='serif',size=14.0)
fig, axs = plt.subplots(nrows,3,sharex=True,sharey=True,figsize=(8.,6.4*(nrows/4.)))
axs = axs.ravel()
# Loop through scales j and directions n
for j in xrange(0,3):
    for n in xrange(0, N):
        f_wav = np.load('/Users/keir/Documents/spin_silc/wavelet_maps/ffp8_cmb_scl_353_wav_70_hybridTestF_3_1_j' + str(j) + '_n' + str(n+1) + '.npy')
        
        # Make the plot!
        myplot_moll(-1.*f_wav.imag, L, axs[(j)*N+n],title=r'Scale %i, direction %i' % (j+1, n+1))

# Pretty adjustment
fig.subplots_adjust(wspace=0, hspace=0.2, left=0.01, right=0.99, bottom=0.01, top=0.95)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_pol/cmb_decomp_scalar_B_axisym.png') #pdf',dpi=400)


