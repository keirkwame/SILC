import numpy as np
import healpy as hp
import multiprocessing as mg
import pys2let as ps

def analworker(i):
    print "This is analysis worker starting for map", i[1]+1, "/", nmaps
    alms = hp.map2alm(i[0],lmax=ellmax-1)
    hp.almxfl(alms,pixrecip,inplace=True) #Correcting for pixwin
    
    wav_maps,scal_maps = ps.analysis_lm2wav(alms,wavparam,ellmax,jmin,ndir,spin,upsample)
    del alms
    np.save(scal_outfits[i[1]],scal_maps)
    del scal_maps
    np.save(wav_outfits[i[1]],wav_maps)
    del wav_maps
    
    return 0

if __name__ == "__main__":
    ##Input
    nprocess = 8
    nmaps = 9 #No. maps (WMAP = 5) (Planck = 9)
    ellmax = 3999 #S2LET parameters - actually band-limits to 1 less
    wavparam = 2
    ndir = 1 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals
    upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
    jmin = 6
    jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

    fitsdir = '/home/keir/s2let_ilc_data/' #'/Users/keir/Documents/s2let_ilc_planck/deconv_data/'
    fitsroot = 'planck_deconv_tapered_' #'ffp6_combined_mc_0000_deconv_' #'planck_deconv_lmax3400_' #'simu_dirty_beam_wmap_9yr_' #'wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    fitscode = ['30','44','70','100','143','217','353','545','857'] #['K','Ka','Q','V','W']
    fitsend = '_pr2.fits' #'.fits'
    fits = [None]*nmaps
    for i in xrange(len(fits)):
        fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

    outdir = fitsdir
    outroot = fitsroot #'ffp6_combined_mc_0000_deconv_' #'planck_deconv_' #'simu_dirty_beam_wmap_9yr_' #'wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    outcode = fitscode #['k','ka','q','v','w']
    scal_outfits = [None]*nmaps
    wav_outfits = [None]*nmaps
    for i in xrange(len(scal_outfits)):
        scal_outfits[i] = outdir + outroot + outcode[i] + '_scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_outfits[i] = outdir + outroot + outcode[i] + '_wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'

    #Load CMB maps
    mapsextra = [None]*(nmaps-1)
    for i in xrange(0,len(mapsextra)):
        mapsextra[i] = (hp.read_map(fits[i+1]),i+1)

    #Calculate band-limited alms and analyse
    print "\nBand-limiting input maps and analysing"
    pixrecip = np.reciprocal(hp.pixwin(hp.get_nside(mapsextra[0][0]))[:ellmax]) #pixwin
    pool = mg.Pool(nprocess)
    anal_output = pool.map(analworker,mapsextra)
    pool.close()
    pool.join()
    del anal_output

    #Calculate wavelet and scaling function maps for each channel
    #Serial version for single map
    '''print "Calculating scaling function and wavelet maps for input map"
    pixrecip = np.reciprocal(hp.pixwin(hp.get_nside(mapsextra[0]))[:ellmax]) #pixwin
    anal_output = analworker(mapsextra)
    del anal_output'''
