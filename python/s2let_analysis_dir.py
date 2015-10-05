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
    
    #Splitting up output wavelet maps
    for j in xrange(jmin,jmax+1):
        for n in xrange(0,ndir):
            offset,scale_lmax,nelem,nelem_wav = ps.wav_ind(j,n,wavparam,ellmax,ndir,jmin,upsample)
            wav_outfits = wav_outfits_root[i[1]] + '_j' + str(j) + '_n' + str(n+1) + '.npy'
            np.save(wav_outfits,wav_maps[offset:offset+nelem])

    del wav_maps
    
    return 0

if __name__ == "__main__":
    ##Input
    nprocess = 9
    nmaps = 9 #No. maps (WMAP = 5) (Planck = 9)
    ellmax = 2500 #S2LET parameters - actually band-limits to 1 less
    wavparam = 1.2
    wavparam_str = '1dot2'
    ndir = 1 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals
    upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
    jmin = 25
    jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

    fitsdir = '/home/keir/s2let_ilc_data/ffp6_pla_data/' #'/Users/keir/Documents/s2let_ilc_planck/deconv_data/'
    fitsroot = 'ffp6_pla_deconv_' #'planck_deconv_tapered_thresh_minusgaussps2_' #'ffp6_fiducial_withPS_tapered_' #'planck_deconv_tapered_noPS_' #'ffp6_fiducial_noPS_tapered_' #'ffp6_combined_mc_0000_deconv_' #'planck_deconv_lmax3400_' #'simu_dirty_beam_wmap_9yr_' #'wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    fitscode = ['30','44','70','100','143','217','353','545','857'] #['K','Ka','Q','V','W']
    fitsend = '.fits' #'_pr2.fits'
    fits = [None]*nmaps
    for i in xrange(len(fits)):
        fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

    outdir = fitsdir
    outroot = fitsroot #'planck_deconv_tapered_pr1_noPS_' #'ffp6_combined_mc_0000_deconv_' #'planck_deconv_' #'simu_dirty_beam_wmap_9yr_' #'wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    outcode = fitscode #['k','ka','q','v','w']
    scal_outfits = [None]*nmaps
    wav_outfits_root = [None]*nmaps
    for i in xrange(len(scal_outfits)):
        scal_outfits[i] = outdir + outroot + outcode[i] + '_scal_' + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_outfits_root[i] = outdir + outroot + outcode[i] + '_wav_' + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir)

    #Load CMB maps
    mapsextra = [None]*nmaps
    #j = 0
    for i in xrange(0,len(mapsextra)):
        mapsextra[i] = (hp.read_map(fits[i]),i)
        #j+=1

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
    i = 0
    mapsextra = (hp.read_map(fits[i]),i)
    pixrecip = np.reciprocal(hp.pixwin(hp.get_nside(mapsextra[0]))[:ellmax]) #pixwin
    anal_output = analworker(mapsextra)
    del anal_output'''
