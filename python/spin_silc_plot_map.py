import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #Latex font settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=35.0) #17 for maps #40 for gnomview
    
    #Increase line thicknesses
    plt.rc('axes', linewidth=2.0)
    plt.rc('xtick.major', width=2.0)
    plt.rc('xtick.minor', width=2.0)
    plt.rc('ytick.major', width=2.0)
    plt.rc('ytick.minor', width=2.0)
    plt.rc('lines', linewidth=1.5)

    LIM = [15,15]
    fname_root = '/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_pol/big_'
    cmb = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_planck_pol_diffusePS_deconv_2253_hybridE_6_1_recon_EBmaps.fits',field=(1,2))
    lab = r'Spin-SILC'
    fnam = 'silcN1'
    
    '''mask = hp.read_map('/Users/keir/Documents/spin_silc/masks/planck_pol_PSmask_nside1024_05thresh.fits') #PS mask
    cmb[0][np.where(mask==False)] = None
    cmb[1][np.where(mask==False)] = None'''

    pol_labs = [r' [$E$]',r' [$B$]'] #[r' $E$ map',r' $B$ map']
    pol_fnams = ['_Emap.pdf','_Bmap.pdf']
    
    gnom_labs = [r'$E$',r'$B$']

    for k in xrange(2):
        title = lab + pol_labs[k]
        fname = fname_root + fnam + pol_fnams[k]
        hp.mollview(cmb[k]*1.e6,unit=r'$\mu\mathrm{K}$',title=title,min=-1.*LIM[k],max=LIM[k]) #,xsize=700) #,margins=[0,0.1,1,0.9])
        plt.savefig(fname)
        
        #Gnomonic view
        '''hp.gnomview(cmb[k]*1.e6,title=gnom_labs[k]+r' zoom',notext=True,cbar=False)
        plt.savefig(fname[:-4] + '_zoom.pdf')'''
