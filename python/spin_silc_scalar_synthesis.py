import numpy as np
import healpy as hp
import math as mh
import pys2let as ps

if __name__ == "__main__":
    ##Input
    #Set directory structure
    comp = 0
    
    if comp == 0: #Keir's iMac
        fitsdir = '/Users/keir/Documents/spin_silc/ilc_maps/'
        outdir = '/Users/keir/Documents/spin_silc/recon_maps/'
    elif comp == 1: #Hypatia
        fitsdir = ''
        outdir = ''

    outnside = 512
    ellmax = 917
    jmin = 6
    lamdas = np.array([2,1.3]) #60,2])
    l_transitions = np.array([513]) #61])
    wavparam_code = 'TestD'
    ndir = 1 #No. directions for each wavelet scale
    spin = 0 #Synthesise with scalar wavelets
    
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
    scal_bandlims = 2*(scal_bandlims-1) + 1
    print scal_bandlims
    wav_bandlims = 2*(wav_bandlims-1) + 1
    print wav_bandlims

    res = ps.verify_tiling(smoothing_lmax,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    if res == False:
        print '\nCheck that above admissibility condition is consistent with zero-padding of wavelets\n'
        pass
    else:
        print '\nA valid wavelet tiling has been chosen.\n'

    fitsroot = 'spin_silc_fwhm10_planck_pol_diffuse_deconv_'
    scal_fits = fitsdir + fitsroot + 'scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_fits_root = fitsdir + fitsroot + 'wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    outroot = fitsroot
    outfits_root = outdir + outroot + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '_recon'

    #Load scaling function map
    scal_map = np.load(scal_fits)

    print "Loading input wavelet maps"
    for j in xrange(jmax+1): #Loading sliced wavelet maps
        for n in xrange(ndir):
            wav_fits = wav_fits_root + '_j' + str(j) + '_n' + str(n+1) + '.npy'
            if j == 0 and n == 0:
                wav_map = np.load(wav_fits)
                continue
            wav_map_part = np.load(wav_fits)
            wav_map = np.concatenate((wav_map,wav_map_part))
            del wav_map_part

    #Synthesise final maps
    print "Synthesising final a_lm"
    e_alms_mw = ps.synthesis_wav2lm_manualtiling(-1.*wav_map.real+0.j,-1.*scal_map.real+0.j,smoothing_lmax,ndir,spin,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims) #Including l=0,1
    b_alms_mw = ps.synthesis_wav2lm_manualtiling(-1.*wav_map.imag+0.j,-1.*scal_map.imag+0.j,smoothing_lmax,ndir,spin,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims) #Add 0 imag part to make complex arrays

    ell = np.arange(2,smoothing_lmax)
    newbeam = hp.gauss_beam(mh.radians(10./60.),lmax=smoothing_lmax-1) / hp.gauss_beam(mh.radians(5./60.),lmax=smoothing_lmax-1)
    filterbeam = np.concatenate((np.zeros(20),0.5*(1 - np.cos((mh.pi*(np.arange(20,41)-20))/20)),np.ones(smoothing_lmax-41))) #High-pass filter
    E_alms_hp = hp.almxfl(ps.lm2lm_hp(e_alms_mw,smoothing_lmax),newbeam*filterbeam) #Reorder and expand to HPX standard
    B_alms_hp = hp.almxfl(ps.lm2lm_hp(b_alms_mw,smoothing_lmax),newbeam*filterbeam) #Same beams as NILC

    print "Calculating final C_l"
    final_cls = hp.alm2cl((E_alms_hp,B_alms_hp)) #(EE,BB,EB)
    clscode = ['EE','BB','EB']
    for i in xrange(len(final_cls)):
        cl_outfits = outfits_root + '_' + clscode[i] + 'cls.fits'
        hp.write_cl(cl_outfits,final_cls[i])

    print "Calculating final maps"
    final_maps = hp.alm2map_spin((E_alms_hp,B_alms_hp),outnside,2,lmax=smoothing_lmax-1) #(Q,U)
    T_map = [np.zeros_like(final_maps[0]),]
    map_outfits = outfits_root + '_QUmaps.fits'
    hp.write_map(map_outfits,T_map+final_maps)

    #Downgraded map
    print "Downgrading maps"
    QU_alms = hp.map2alm(final_maps,lmax=256,pol=False)
    bigbeam = hp.gauss_beam(mh.radians(80./60.),lmax=256) / (hp.gauss_beam(mh.radians(10./60.),lmax=256) * np.concatenate((np.ones(2),hp.pixwin(hp.get_nside(final_maps[0]),pol=True)[1][2:257])))
    hp.almxfl(QU_alms[0],bigbeam,inplace=True)
    hp.almxfl(QU_alms[1],bigbeam,inplace=True)
    dg_maps = hp.alm2map(QU_alms,nside=128,pixwin=True,pol=False)
    T_map_dg = [np.zeros_like(dg_maps[0]),]
    map_outfits = outfits_root + '_QUmaps_dg.fits'
    hp.write_map(map_outfits,T_map_dg+dg_maps)

    #Masked spectra
    print "Calculating masked C_l"
    mask_name = 'NILCconfsky'
    mask = hp.read_map('/Users/keir/Documents/spin_silc/masks/nilc_pol_conmask_nside512_05thresh.fits') #0 where holes
    #mask = hp.read_map('/Users/keir/Documents/spin_silc/masks/UP78_PR2_nside512_05thresh.fits') #0 where holes
    f_sky = np.sum(mask) / len(mask)
    print "f_sky =", f_sky
    masked_maps = [final_maps[0]*mask,final_maps[1]*mask] #(Q*mask,U*mask)
    masked_cls = hp.anafast(T_map+masked_maps,lmax=ellmax-1) #(TT,EE,BB,TE,EB,TB)
    pixrecip = np.concatenate((np.ones(2),np.reciprocal(hp.pixwin(hp.get_nside(masked_maps[0]),pol=True)[1][2:ellmax]))) #P pixwin #Not defined for l < 2
    clsidx = [1,2,4]
    for i in xrange(len(clscode)):
        cl_outfits = outfits_root + '_' + clscode[i] + 'cls_' + mask_name + '.fits'
        hp.write_cl(cl_outfits,masked_cls[clsidx[i]] * pixrecip * pixrecip / f_sky)

    #Some utilities for quick plotting
    
