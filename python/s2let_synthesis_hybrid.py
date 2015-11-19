import numpy as np
import scipy as sp
import healpy as hp
import math as mh
import multiprocessing as mg
import pys2let as ps

def reducealms(almsextra): #(alms,scale_lmax,scale_lmin)
    for em in xrange(almsextra[2]+1): #Up to and including scale_lmin
        startindex = em*almsextra[1] - .5*em*(em-1)
        almsextra[startindex:startindex+almsextra[2]+1-em] = 0.+0.j #Zero-ing alm's below and including scale_lmin
    return almsextra

def varworker(wav_map_indices): #(j,n)
    wav_fits = wav_fits_root + '_j' + str(wav_map_indices[0]) + '_n' + str(wav_map_indices[1]+1) + '.npy'
    map = np.load(wav_fits,mmap_mode='r')
    scale_lmax = wavparam**(wav_map_indices[0]+1) #lambda^(j+1)
    if scale_lmax > ellmax:
        scale_lmax = ellmax
    scale_lmin = wavparam**(wav_map_indices[0]-1) #lambda^(j-1)
    alms = ps.map2alm_mw(map,scale_lmax,spin)
    del map
    alms = reducealms((alms,scale_lmax,scale_lmin))
    map_reduced = ps.alm2map_mw(alms,scale_lmax,spin)
    del alms
    map_squared = np.square(map_reduced.real)
    del map_reduced
    return map_squared

def variance(scal_map):
    #Scaling function calculation
    scal_map_var = np.square(scal_map.real) + 0.j

    #Wavelet calculation
    wav_map_indices = [None]*(jmax+1-jmin)*ndir
    i = 0
    for j in xrange(jmin,jmax+1):
        for n in xrange(0,ndir):
            wav_map_indices[i] = (j,n)
            i+=1
    pool = mg.Pool(nprocess)
    wav_map_var_list = pool.map(varworker,wav_map_indices)
    pool.close()
    pool.join()
    wav_map_var = wav_map_var_list[0] + 0.j
    for i in xrange(1,len(wav_map_var_list)):
        wav_map_var = np.concatenate((wav_map_var,wav_map_var_list[i] + 0.j))
    del wav_map_var_list

    return scal_map_var, wav_map_var

if __name__ == "__main__":
    ##Input
    #Set directory structure
    comp = 0
    
    if comp == 0: #Keir's iMac
        fitsdir = '/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/'
    elif comp == 1: #Hypatia
        fitsdir = '/home/keir/s2let_ilc_data/hybrid_data/'

    nprocess = 1
    outnside = 256
    ellmax = 300
    jmin = 0
    lamdas = np.array([60,2])
    wavparam_code = 'C'
    l_transitions = np.array([61])
    ndir = 5 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals

    fitsroot = 's2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax1300_'
    #fitsroot = 'planck2015_2_cmb_map_1_' #TESTING
    scal_fits = fitsdir + fitsroot + 'scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_fits_root = fitsdir + fitsroot + 'wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)
    
    outfits = fitsdir + fitsroot + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '_recon.fits' #Output is HPX map
    outclfits = fitsdir + fitsroot + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '_recon_cls.fits'
    varfits = fitsdir + fitsroot + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '_variance.fits' #Variance map

    #Construct valid hybrid tiling
    scal_tiles, wav_tiles, scal_bandlims, wav_bandlims, jmax, l_bounds = ps.construct_hybrid_tiling(ellmax,jmin,lamdas,l_transitions)
    res = ps.verify_tiling(ellmax,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    if res == False:
        raise ValueError('\nAn invalid wavelet tiling has been chosen.\n')
    else:
        print '\nA valid wavelet tiling has been chosen.\n'

    #Load scaling function map
    scal_map = np.load(scal_fits)

    print "Loading input wavelet maps"
    for j in xrange(jmax+1): #Loading sliced wavelet maps
        for n in xrange(ndir):
            wav_fits = wav_fits_root + '_' + str(wavparam_code) + str(j) + '_n' + str(n+1) + '.npy'
            if j == 0 and n == 0:
                wav_map = np.load(wav_fits)
                continue
            wav_map_part = np.load(wav_fits)
            wav_map = np.concatenate((wav_map,wav_map_part))
            del wav_map_part

    #VARIANCE calculation
    '''scal_map_var,wav_map_var = variance(scal_map)
    print "Synthesising variance alm's"
    var_alms = ps.synthesis_wav2lm(wav_map_var,scal_map_var,wavparam,ellmax,jmin,ndir,spin,upsample)
    print "Calculating variance map"
    var_map = hp.alm2map(var_alms,nside=outnside,pixwin=True)
    hp.write_map(varfits,var_map)'''

    #Synthesise final map
    print "Synthesising final alm's"
    final_alms = ps.synthesis_wav2lm_manualtiling(wav_map,scal_map,ellmax,ndir,spin,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    print "Calculating final map"
    final_map = hp.alm2map(final_alms,nside=outnside,pixwin=True)
    hp.write_map(outfits,final_map)
    print "Calculating final cl's"
    final_cls = hp.alm2cl(final_alms)
    hp.write_cl(outclfits,final_cls)
    ell = np.arange(len(final_cls))
    invtwopi = 1./(2.*mh.pi)

    #Binning final power spectrum for plotting
    '''binlen = 8 #Should be factor of len(final_cls)
    final_cls_binned = np.mean(np.reshape(final_cls,(-1,binlen)),axis=-1)
    ell_binned = np.mean(np.reshape(ell,(-1,binlen)),axis=-1)'''


