import numpy as np
import healpy as hp
import multiprocessing as mg
import pys2let as ps

def analworker(i):
    print "This is analysis worker starting for map", i+1, "/", nmaps
    QU_maps = hp.read_map(fits[i],field=(1,2)) #(Q,U)
    pixrecip = np.concatenate((np.ones(2),np.reciprocal(hp.pixwin(hp.get_nside(QU_maps[0]),pol=True)[1][2:smoothing_lmax]))) #P pixwin #Not defined for l < 2
    pm_alms = hp.map2alm_spin(QU_maps,spin,lmax=smoothing_lmax-1)
    del QU_maps
    hp.almxfl(pm_alms[0],pixrecip,inplace=True) #Correcting for pixwin
    hp.almxfl(pm_alms[1],pixrecip,inplace=True) #Correcting for pixwin
    
    #Reorder to S2LET alms
    pm_alms[0] = ps.lm_hp2lm(pm_alms[0],smoothing_lmax)
    pm_alms[1] = ps.lm_hp2lm(pm_alms[1],smoothing_lmax)
    
    P_alms = -1.*pm_alms[0] -1.j*pm_alms[1] #CHECK THIS IS CORRECT!
    del pm_alms
    
    wav_maps, scal_maps = ps.analysis_lm2wav_manualtiling(P_alms,smoothing_lmax,ndir,spin,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    del P_alms
    np.save(scal_outfits[i],scal_maps)
    del scal_maps
    
    #Splitting up output wavelet maps
    offset = 0
    for j in xrange(jmax+1):
        for n in xrange(ndir):
            bandlim = wav_bandlims[j]
            nelem = bandlim*(2.*bandlim-1.)
            wav_outfits = wav_outfits_root[i] + '_j' + str(j) + '_n' + str(n+1) + '.npy'
            np.save(wav_outfits,wav_maps[offset:offset+nelem])
            offset += nelem
    del wav_maps

    return 0

if __name__ == "__main__":
    ##Input
    #Set directory structure
    comp = 0
    
    if comp == 0: #Keir's iMac
        nprocess = 1
        fitsdir = '/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/'
        outdir = '/Users/keir/Documents/spin_silc/wavelet_maps/'
    elif comp == 1: #Hypatia
        nprocess = 7
        fitsdir = ''
        outdir = ''

    nmaps = 1
    ellmax = 70 #300
    jmin = 3
    lamdas = np.array([2,2.02]) #60,2
    l_transitions = np.array([65]) #61,2017
    wavparam_code = 'TestF'
    ndir = 1 #No. directions for each wavelet scale
    spin = 2

    #Construct valid hybrid tiling
    scal_tiles, wav_tiles, scal_bandlims, wav_bandlims, jmax, l_bounds = ps.construct_hybrid_tiling(ellmax,jmin,lamdas,l_transitions)

    #Remove truncated wavelets
    ntrunc = 1 #1 or 2
    ellmax = wav_bandlims[-1-1*ntrunc] #Lower total bandlimit to match smallest wavelet
    print "New bandlimit =", ellmax
    scal_tiles = scal_tiles[:ellmax]
    wav_tiles = wav_tiles[:ellmax,:-1*ntrunc]
    wav_bandlims = wav_bandlims[:-1*ntrunc]
    jmax = jmax - ntrunc
    print jmax

    #Zeropad wavelets to allow double resolution
    smoothing_lmax = 2*(ellmax-1) + 1
    print "l_max = %i, smoothing l_max = %i" %(ellmax,smoothing_lmax)
    scal_tiles = np.concatenate((scal_tiles,np.zeros(smoothing_lmax-ellmax)))
    print len(scal_tiles)
    wav_tiles = np.concatenate((wav_tiles,np.zeros((smoothing_lmax-ellmax,len(wav_tiles[0])))))
    print len(wav_tiles), len(wav_tiles[0])
    scal_bandlims = smoothing_lmax #2*(scal_bandlims-1) + 1
    print scal_bandlims
    wav_bandlims[:] = smoothing_lmax #2*(wav_bandlims-1) + 1
    print wav_bandlims

    res = ps.verify_tiling(smoothing_lmax,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    if res == False:
        print '\nCheck that above admissibility condition is consistent with zero-padding of wavelets\n'
        pass
    else:
        print '\nA valid wavelet tiling has been chosen.\n'

    fitsroot = 'ffp8_cmb_scl_'
    fitscode = ['353'] #,'044','070','100','143','217','353']
    fitsend = '_full_map.fits'
    fits = [None]*nmaps
    for i in xrange(len(fits)):
        fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

    outroot = fitsroot
    outcode = fitscode
    scal_outfits = [None]*nmaps
    wav_outfits_root = [None]*nmaps
    for i in xrange(len(scal_outfits)):
        scal_outfits[i] = outdir + outroot + outcode[i] + '_scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_outfits_root[i] = outdir + outroot + outcode[i] + '_wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    #Calculate band-limited alms and analyse
    print "Band-limiting input maps and analysing"
    pool = mg.Pool(nprocess)
    anal_output = pool.map(analworker,np.arange(nmaps))
    pool.close()
    pool.join()
