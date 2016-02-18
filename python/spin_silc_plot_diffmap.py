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
    
    fname_root = '/Users/keir/Documents/spin_silc/plots/diffusePS_fwhm50_'

    maps = [None]*nmaps
    maps[0] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/smica_QUmaps_dg.fits',field=(1,2)) #SMICA
    maps[1] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/compsepMaps/nilc_QUmaps_dg.fits',field=(1,2)) #NILC
    maps[2] = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_fwhm50_planck_pol_diffusePS_deconv_917_hybridTestD_6_1_recon_QUmaps_dg.fits',field=(1,2)) #(Q,U) #SILC
    
    labs = [r'SMICA',r'NILC',r'SILC $(N = 1)$']
    fnams = ['smica','nilc','silcN1']

    pol_labs = [r' ($Q$)',r' ($U$)']
    pol_fnams = ['_Qdiffmap.pdf','_Udiffmap.pdf']
    for i in xrange(1,nmaps): #1,2
        for j in xrange(i): #0;0,1
            for k in xrange(2):
                title = labs[i] + r' - ' + labs[j] + pol_labs[k]
                fname = fname_root + fnams[i] + '_' + fnams[j] + pol_fnams[k]
                hp.mollview(maps[i][k]*1.e6-maps[j][k]*1.e6,unit=r'$\mu\mathrm{K}$',title=title,min=-1.*LIM,max=LIM) #1-0;2-0;2-1
                plt.savefig(fname)
