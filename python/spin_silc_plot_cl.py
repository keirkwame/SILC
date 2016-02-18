import numpy as np
import healpy as hp
import math as mh
import matplotlib.pyplot as plt
import distinct_colours as dc

if __name__ == "__main__":
    #Latex font settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=21.0) #21 for power spectra

    #Increase line thicknesses
    plt.rc('axes', linewidth=2.0)
    plt.rc('xtick.major', width=2.0)
    plt.rc('xtick.minor', width=2.0)
    plt.rc('ytick.major', width=2.0)
    plt.rc('ytick.minor', width=2.0)
    plt.rc('lines', linewidth=1.5)

    ncls = 4
    lmin = 2
    lmax = 701
    xmax = 750
    bin_len = 10
    resid_idx = 3
    
    fname = '/Users/keir/Documents/spin_silc/plots/EEcls_PSfsky.pdf'

    dlconv = (np.arange(lmin,lmax+1) * (np.arange(lmin,lmax+1) + 1.)) / (2. * mh.pi)
    
    #Beams
    planck_beam = hp.gauss_beam(mh.radians(10./60.),lmax=lmax)[lmin:]
    filterbeam = np.concatenate((np.zeros(18),0.5*(1 - np.cos((mh.pi*(np.arange(20,41)-20))/20)),np.ones(lmax-40))) #High-pass filter
    #pixrecip = np.concatenate((np.ones(2),np.reciprocal(hp.pixwin(1024,pol=True)[1][2:lmax+1])))[2:] #P pixwin #Not defined for l < 2

    cls = [None]*ncls
    cls[0] = hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noN_planck_pol_deconv_917_hybridTestD_6_1_recon_EEcls_PSfsky.fits')[lmin:lmax+1]*1.e12
    cls[1] = hp.read_cl('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/nilc_PSfsky_lmax2300_EEcls.fits')[lmin:lmax+1]*1.e12
    cls[2] = hp.read_cl('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/smica_PSfsky_lmax2300_EEcls.fits')[lmin:lmax+1]*1.e12
    cls[3] = (np.loadtxt('/Users/keir/Software/CAMB/planck_lmax2350_totCls.dat',usecols=(2,))[:lmax-1] * planck_beam * planck_beam * filterbeam * filterbeam) / dlconv #(l,TT,EE,BB,TE)
    
    lins = ['-','-','-','--']
    dis_cols = dc.get_distinct(2)
    cols = ['k',dis_cols[0],dis_cols[1],'k']
    labs = [r'SILC $(N = 1)$ [PS mask]',r'NILC [PS mask]',r'SMICA [PS mask]',r'Theory']

    dls = cls * dlconv[None,:]
    dls_bin = np.mean(np.reshape(dls,(ncls,-1,bin_len)),axis=-1)
    ell_bin = np.mean(np.reshape(np.arange(lmin,lmax+1),(-1,bin_len)),axis=-1)
    
    f, (ax0,ax1) = plt.subplots(2,sharex=True,figsize=(8,8))
    for i in xrange(ncls):
        ax0.plot(ell_bin,dls_bin[i],ls=lins[i],color=cols[i],label=labs[i])
    ax0.set_xlim([0,xmax])
    ax0.set_ylim([0,40])
    ax0.set_ylabel(r'$D_{\ell}$ $[{\mu\mathrm{K}}^2]$') #, labelpad = 3)
    ax0.legend(prop={'size':21},frameon=False,loc='upper left')

    for i in xrange(ncls):
        ax1.plot(ell_bin,dls_bin[i]-dls_bin[resid_idx],ls=lins[i],color=cols[i],label=labs[i])
    ax1.set_xlim([0,xmax])
    ax1.set_ylim([0,18])
    ax1.set_xlabel(r'Multipole $\ell$', labelpad = 10)
    ax1.set_ylabel(r'$\Delta D_{\ell}$ $[{\mu\mathrm{K}}^2]$') #, labelpad = -4)

    yticks = ax1.yaxis.get_major_ticks()
    yticks[0].label1.set_visible(False)
    yticks[-1].label1.set_visible(False)

    f.subplots_adjust(hspace=0,right=0.99)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    plt.savefig(fname)
