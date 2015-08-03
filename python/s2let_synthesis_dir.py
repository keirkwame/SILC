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
    nprocess = 1
    outnside = 2048
    ellmax = 2500 #S2LET parameters - actually band-limits to 1 less
    wavparam = 2
    wavparam_str = '2'
    ndir = 1 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals
    upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
    jmin = 6
    jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

    fitsdir = '/home/keir/s2let_ilc_data/ffp6_data_withPS/' #'/Users/keir/Documents/s2let_ilc_planck/deconv_data/'
    fitsroot = 's2let_ilc_dir_hypatia_ffp6_fiducial_withPS_tapered_' #'s2let_ilc_dir_hypatia_planck_deconv_tapered_minusgaussps_' #'s2let_ilc_dir_hypatia_memeff_planck_deconv_tapered_pr1_noPS_' #'s2let_ilc_dir_hypatia_memeff_ffp6_fiducial_noPS_tapered_' #'s2let_ilc_dir_hypatia_ffp6_combined_mc_0000_deconv_' #'s2let_ilc_dir_para_gauss_planck_deconv_' #'s2let_ilc_dir_para_gauss_simu_dirty_beam_wmap_9yr_' #'s2let_ilc_dir_para_gauss_wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    scal_fits = fitsdir + fitsroot + 'scal_' + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_fits_root = fitsdir + fitsroot + 'wav_' + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir)

    print "Loading input wavelet maps"
    for j in xrange(jmin,jmax+1): #Loading sliced wavelet maps
        for n in xrange(0,ndir):
            wav_fits = wav_fits_root + '_j' + str(j) + '_n' + str(n+1) + '.npy'
            
            #Override some input maps
            '''if j < 7:
                wav_fits = fitsdir + 's2let_ilc_dir_para_gauss_wmap_deconv_nosource_smoothw_extrapolated_9yr_' + 'wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '_j' + str(j) + '_n' + str(n+1) + '_testoptimise.npy'
            '''
            
            if j == jmin and n == 0:
                wav_map = np.load(wav_fits)
                continue
            wav_map_part = np.load(wav_fits)
            wav_map = np.concatenate((wav_map,wav_map_part))
            del wav_map_part

    outfits = fitsdir + fitsroot + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir) + '_recon.fits' #Output is HPX map
    outclfits = fitsdir + fitsroot + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir) + '_recon_cls.fits'
    varfits = fitsdir + fitsroot + str(ellmax) + '_' + wavparam_str + '_' + str(jmin) + '_' + str(ndir) + '_variance.fits' #Variance map

    #Override some input maps
    #scal_fits = 'deconv_data/s2let_ilc_dir_para_gauss_wmap_deconv_smoothw_extrapolated_9yr_scal_1024_2_6_3.npy'

    #Load scaling function map
    scal_map = np.load(scal_fits)

    #VARIANCE calculation
    '''scal_map_var,wav_map_var = variance(scal_map)
    print "Synthesising variance alm's"
    var_alms = ps.synthesis_wav2lm(wav_map_var,scal_map_var,wavparam,ellmax,jmin,ndir,spin,upsample)
    print "Calculating variance map"
    var_map = hp.alm2map(var_alms,nside=outnside,pixwin=True)
    hp.write_map(varfits,var_map)'''

    #Synthesise final map
    print "Synthesising final alm's"
    final_alms = ps.synthesis_wav2lm(wav_map,scal_map,wavparam,ellmax,jmin,ndir,spin,upsample)
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
