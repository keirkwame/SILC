import numpy as np
import scipy as sp
import healpy as hp
import math as mh
import multiprocessing as mg
import pys2let as ps

##Input
nprocess = 4
outnside = 512
ellmax = 1024 #S2LET parameters - actually band-limits to 1 less
wavparam = 2
ndir = 6 #No. directions for each wavelet scale
spin = 0 #0 for temp, 1 for spin signals
upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
jmin = 6
jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

fitsdir = '/Users/keir/Documents/s2let_ilc/simu_data/'
fitsroot = 's2let_ilc_dir_para_gauss_simu_dirty_beam_wmap_9yr_' #'s2let_ilc_dir_para_gauss_wmap_deconv_nosource_smoothw_extrapolated_9yr_' #'s2let_ilc_dir_para_gauss_planck_deconv_'
fitsend = '.fits'
scal_fits = fitsdir + fitsroot + 'scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
wav_fits_root = fitsdir + fitsroot + 'wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir)

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

outfits = fitsdir + fitsroot + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '_recon.fits' #Output is HPX map
outclfits = fitsdir + fitsroot + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '_recon_cls.fits'

#Override some input maps
#scal_fits = 'deconv_data/s2let_ilc_dir_para_gauss_wmap_deconv_smoothw_extrapolated_9yr_scal_1024_2_6_3.npy'

#Load wavelet and scaling function maps
scal_map = np.load(scal_fits)

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
binlen = 8 #Should be factor of len(final_cls)
final_cls_binned = np.mean(np.reshape(final_cls,(-1,binlen)),axis=-1)
ell_binned = np.mean(np.reshape(ell,(-1,binlen)),axis=-1)
