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

    ncls = 2
    lmin = 2
    lmax = 2041
    xmax = 2100
    bin_len = 10
    resid_idx = 1
    
    fname = '/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_pol/cl_ffp8_N1_2.pdf'

    dlconv = (np.arange(lmin,lmax+1) * (np.arange(lmin,lmax+1) + 1.)) / (2. * mh.pi)
    
    #Beams
    planck_beam = hp.gauss_beam(mh.radians(10./60.),lmax=lmax)[lmin:]
    filterbeam = np.concatenate((np.zeros(18),0.5*(1 - np.cos((mh.pi*(np.arange(20,41)-20))/20)),np.ones(lmax-40))) #High-pass filter
    #pixrecip = np.concatenate((np.ones(2),np.reciprocal(hp.pixwin(1024,pol=True)[1][2:lmax+1])))[2:] #P pixwin #Not defined for l < 2
    #beam353 = np.load('/Users/keir/Documents/spin_silc/beams/planck_bl_353_pr2.npy')[lmin:lmax+1]
    #newbeam = hp.gauss_beam(mh.radians(10./60.),lmax=lmax)[lmin:] / hp.gauss_beam(mh.radians(5./60.),lmax=lmax)[lmin:]

    cls = {}
    dls = {}
    dls_bin = {}
    cross_corr = ['EE','BB','EB']
    camb_idx = [2,3,3]
    k=0

    for i in cross_corr:
        cls[i] = [None]*ncls
        #cls[i][0] = hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_planck_pol_diffusePS_deconv_2253_hybridE_6_1_recon_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)
        #cls[i][1] = hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_planck_pol_diffusePS_deconv_2253_hybridE_6_5_recon_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)

        cls[i][0] = hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noFilter_ffp8_pol_diffusePS_deconv_2253_hybridE_6_1_recon_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)
        #cls[i][1] = hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noFilter_ffp8_pol_diffusePS_deconv_2253_hybridE_6_5_recon_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)
        #cls[i][2] = hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noFilter_ffp8_pol_diffusePS_deconv_2253_hybridE_6_10_recon_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)
        #cls[i][2] = np.concatenate((hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_ffp8_pol_diffusePS_deconv_917_hybridTestD_6_10_recon_' + i + 'cls.fits')[lmin:700],hp.read_cl('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_ffp8_pol_diffusePS_deconv_2253_hybridE_6_1_recon_' + i + 'cls.fits')[700:lmax+1]))*1.e12 / (planck_beam * planck_beam) #WAITING FOR FULL RESOLUTION RESULTS

        #cls[i][1] = hp.read_cl('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/nilc_lmax2300_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)
        #cls[i][2] = hp.read_cl('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/smica_lmax2300_' + i + 'cls.fits')[lmin:lmax+1]*1.e12 / (planck_beam * planck_beam)

        #cls[i][3] = (np.loadtxt('/Users/keir/Software/CAMB/planck_lmax2350_totCls.dat',usecols=(camb_idx[k],))[:lmax-1] * filterbeam * filterbeam) / dlconv #(l,TT,EE,BB,TE)
        cls[i][1] = hp.read_cl('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/ffp8_cmb_scl_000_full_lmax2300_' + i + 'cls_fullSky.fits')[lmin:lmax+1]*1.e12 #* filterbeam * filterbeam
        k=+1
        
        dls[i] = cls[i] * dlconv[None,:]
        dls_bin[i] = np.mean(np.reshape(dls[i],(ncls,-1,bin_len)),axis=-1)

    ell_bin = np.mean(np.reshape(np.arange(lmin,lmax+1),(-1,bin_len)),axis=-1)

    lins = ['-','--'] #['--','-','-'] #['-','-','-','--']
    dis_cols = dc.get_distinct(2)
    cols = ['k','k'] #['k',dis_cols[0],dis_cols[1],'k']
    labs = [r'SILC $(N = 1)$ [FFP8]',r'Input'] #[r'$N = 1$',r'$N = 5$',r'$N = 10$'] #[r'SILC $(N = 1)$',r'NILC',r'SMICA',r'Theory']

    f, (ax0,ax1,ax2) = plt.subplots(3,sharex=True,figsize=(8,12))
    #f, (ax0,ax1) = plt.subplots(2,sharex=True,figsize=(8,8)) #Directionality
    for i in xrange(ncls):
        ax0.plot(ell_bin,dls_bin['EE'][i] - dls_bin['EE'][resid_idx],ls=lins[i],color=dis_cols[1],label=labs[i],linewidth=0.75)
        ax0.plot(ell_bin,dls_bin['EE'][i],ls=lins[i],color=cols[i],label=labs[i])
    ax0.set_ylabel(r'$D_\ell^{EE}$ $[{\mu\mathrm{K}}^2]$',labelpad = 1) #labelpad = 1 for axisymmetric, -2 for directional
    ax0.set_xlim([0,xmax])
    ax0.set_ylim([1.e-2,1.e4])
    #ax0.set_ylim([-3.25,0.25]) #Directionality
    #ax0.set_xscale('log')
    ax0.set_yscale('log')

    for i in xrange(ncls):
        ax1.plot(ell_bin,dls_bin['BB'][i] - dls_bin['BB'][resid_idx],ls=lins[i],color=dis_cols[1],linewidth=0.75) #FFP
        ax1.plot(ell_bin,dls_bin['BB'][i],ls=lins[i],color=cols[i],label=labs[i])
    ax1.set_ylabel(r'$D_\ell^{BB}$ $[{\mu\mathrm{K}}^2]$',labelpad = 1) #labelpad = 1 for axisymmetric, -2 for directional
    ax1.legend(frameon=False,loc='center right')
    ax1.set_xlim([0,xmax])
    ax1.set_ylim([1.e-2,1.e4-1.e3])
    #ax1.set_ylim([-3.25,0.25]) #Directionality
    ax1.set_yscale('log')
    #ax1.set_xlabel(r'Multipole $\ell$', labelpad = 10) #Directionality

    for i in xrange(ncls):
        ax2.plot(ell_bin,dls_bin['EB'][i],ls=lins[i],color=cols[i],label=labs[i])
    ax2.set_xlabel(r'Multipole $\ell$', labelpad = 10)
    ax2.set_ylabel(r'$D_\ell^{EB}$ $[{\mu\mathrm{K}}^2]$',labelpad = 1)
    ax2.set_xlim([0,xmax])

    #Directionality
    '''yticks0 = ax0.yaxis.get_major_ticks()
    yticks0[-1].label1.set_visible(False)
    yticks0[0].label1.set_visible(False)'''

    yticks = ax1.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)
    #yticks[0].label1.set_visible(False) #Directionality

    yticks2 = ax2.yaxis.get_major_ticks()
    yticks2[0].label1.set_visible(False)
    yticks2[-1].label1.set_visible(False)

    f.subplots_adjust(hspace=0,right=0.99)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    plt.savefig(fname)
