import numpy as np
import healpy as hp
import multiprocessing as mg
import pys2let as ps

def analworker(i):
    print "This is analysis worker starting for map", i[1]+1, "/", nmaps
    alms = hp.map2alm(i[0],lmax=ellmax-1)
    hp.almxfl(alms,pixrecip,inplace=True) #Correcting for pixwin
    
    #TESTING
    '''cls = hp.alm2cl(alms)
    map_bandlim = hp.alm2map(alms,nside=256,pixwin=True)
    hp.write_map('/Users/keir/Documents/s2let_ilc_planck/planck2015_2_cmb_map_1_bandlim300.fits',map_bandlim)'''

    wav_maps, scal_maps = ps.analysis_lm2wav_manualtiling(alms,ellmax,ndir,spin,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    del alms
    np.save(scal_outfits[i[1]],scal_maps)
    del scal_maps
    
    #Splitting up output wavelet maps
    offset = 0
    for j in xrange(jmax+1):
        for n in xrange(ndir):
            bandlim = wav_bandlims[j]
            nelem = bandlim*(2.*bandlim-1.)
            wav_outfits = wav_outfits_root[i[1]] + '_' + wavparam_code + str(j) + '_n' + str(n+1) + '.npy'
            np.save(wav_outfits,wav_maps[offset:offset+nelem])
            offset += nelem
    del wav_maps

    return 0 #wav_maps,scal_maps #,map_bandlim,alms,cls

if __name__ == "__main__":
    ##Input
    #Set directory structure
    comp = 0
    
    if comp == 0: #Keir's iMac
        nprocess = 3
        fitsdir = '/Users/keir/Documents/s2let_ilc_planck/diffuse_data/'
        #fitsdir = '/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/'
    elif comp == 1: #Hypatia
        nprocess = 9
        fitsdir = '/home/keir/s2let_ilc_data/hybrid_data/'
    
    nmaps = 9 #No. maps (Planck = 9)
    ellmax = 300
    jmin = 6
    lamdas = np.array([2,1.9]) #60,2
    wavparam_code = 'C'
    l_transitions = np.array([256]) #61,2017
    ndir = 5 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals

    fitsroot = 'planck_diffuse_deconv_tapered_thresh_lmax1300_'
    #fitsroot = 'ffp8_diffuse_deconv_tapered_thresh_lmax3600_'
    fitscode = ['30','44','70','100','143','217','353','545','857']
    fitsend = '_pr2.fits'
    fits = [None]*nmaps
    for i in xrange(len(fits)):
        fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

    #outdir = fitsdir
    outdir = '/Users/keir/Documents/s2let_ilc_planck/scal_data/'
    outroot = fitsroot
    outcode = fitscode
    scal_outfits = [None]*nmaps
    wav_outfits_root = [None]*nmaps
    for i in xrange(len(scal_outfits)):
        scal_outfits[i] = outdir + outroot + outcode[i] + '_scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_outfits_root[i] = outdir + outroot + outcode[i] + '_wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    #Construct valid hybrid tiling
    scal_tiles, wav_tiles, scal_bandlims, wav_bandlims, jmax, l_bounds = ps.construct_hybrid_tiling(ellmax,jmin,lamdas,l_transitions)
    res = ps.verify_tiling(ellmax,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    if res == False:
        raise ValueError('\nAn invalid wavelet tiling has been chosen.\n')
    else:
        print '\nA valid wavelet tiling has been chosen.\n'

    #TESTING analysis on CMB map
    '''fits = '/Users/keir/Documents/planck2015_2_cmb_realisations/planck2015_2_cmb_map_1.fits'
    scal_outfits = ['/Users/keir/Documents/s2let_ilc_planck/planck2015_2_cmb_map_1_scal_300_hybridC_6_5.npy']
    wav_outfits_root = ['/Users/keir/Documents/s2let_ilc_planck/planck2015_2_cmb_map_1_wav_300_hybridC_6_5']'''

    #Load CMB maps
    mapsextra = [None]*nmaps
    for i in xrange(len(mapsextra)):
        mapsextra[i] = (hp.read_map(fits[i]),i)
    #mapsextra = [(hp.read_map(fits),0)]

    #Calculate band-limited alms and analyse
    print "\nBand-limiting input maps and analysing"
    pixrecip = np.reciprocal(hp.pixwin(hp.get_nside(mapsextra[0][0]))[:ellmax]) #pixwin
    pool = mg.Pool(nprocess)
    anal_output = pool.map(analworker,mapsextra)
    pool.close()
    pool.join()
    #del anal_output


