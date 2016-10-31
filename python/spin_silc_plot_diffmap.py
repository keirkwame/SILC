import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #Latex font settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=17.0) #17 for maps

    #Increase line thicknesses
    plt.rc('axes', linewidth=2.0)
    plt.rc('xtick.major', width=2.0)
    plt.rc('xtick.minor', width=2.0)
    plt.rc('ytick.major', width=2.0)
    plt.rc('ytick.minor', width=2.0)
    plt.rc('lines', linewidth=1.5)

    nmaps = 3
    LIM = 1
    
    fname_root = '/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_pol/'

    maps = [None]*nmaps
    maps[0] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/smica_QUmaps_dg.fits',field=(1,2)) #SMICA
    maps[1] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/nilc_QUmaps_dg.fits',field=(1,2)) #NILC
    maps[2] = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_planck_pol_diffusePS_deconv_2253_hybridE_6_1_recon_QUmaps_dg.fits',field=(1,2)) #(Q,U) #SILC

    '''maps[0] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/ffp8_cmb_scl_000_full_QUmaps_dg.fits',field=(1,2)) #Input
    maps[1] = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noFilter_ffp8_pol_diffusePS_deconv_917_hybridTestD_6_1_recon_QUmaps_dg.fits',field=(1,2)) #N=1'''
    #maps[2] = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_ffp8_pol_diffusePS_deconv_1193_hybridE_6_5_recon_QUmaps_dg.fits',field=(1,2)) #N=5
    
    #mask = hp.read_map('/Users/keir/Documents/spin_silc/masks/UPB77_PR2_nside128_05thresh.fits') #UPB77
    mask = np.ones_like(maps[0][0]) #No mask
    for i in xrange(nmaps):
        for j in xrange(2): #Q,U
            maps[i][j][np.where(mask==False)] = None
    
    labs = [r'SMICA',r'NILC',r'SILC $(N = 1)$'] #[r'input [FFP8]',r'SILC $(N = 1)$']
    fnams = ['smica','nilc','silcN1'] #['input','silcN1']

    pol_labs = [r' [$Q$]',r' [$U$]']
    pol_fnams = ['_QdiffMap.pdf','_UdiffMap.pdf']
    for i in xrange(1,nmaps): #1,2
        for j in xrange(i): #0;0,1
            for k in xrange(2):
                title = labs[i] + r' - ' + labs[j] + pol_labs[k]
                fname = fname_root + fnams[i] + '_' + fnams[j] + pol_fnams[k]
                diffmap = maps[i][k]*1.e6-maps[j][k]*1.e6 #uK
                print fnams[i], '-', fnams[j], pol_fnams[k], 'mean =', np.mean(diffmap[np.where(mask==True)]), 'std =', np.std(diffmap[np.where(mask==True)])
                hp.mollview(diffmap,unit=r'$\mu\mathrm{K}$',title=title,min=-1.*LIM,max=LIM) #1-0;2-0;2-1
                plt.savefig(fname)
