import numpy as np
import healpy as hp
import multiprocessing as mg
import pys2let as ps

def almworker(i):
    print "This is (map2alm) worker starting for another map"
    alms_worker = hp.map2alm(i,lmax=ellmax-1)
    hp.almxfl(alms_worker,pixrecip,inplace=True) #Correcting for pixwin
    return alms_worker

##Input
nprocess = 4
nmaps = 5 #No. maps (WMAP = 5) (Planck = 9)
ellmax = 1024 #S2LET parameters - actually band-limits to 1 less
wavparam = 2
ndir = 1 #No. directions for each wavelet scale
spin = 0 #0 for temp, 1 for spin signals
upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
jmin = 6
jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

fitsdir = 'deconv_data/'
fitsroot = 'wmap_deconv_smoothw_extrapolated_9yr_' #'planck_deconv_lmax2048_'
fitscode = ['K','Ka','Q','V','W'] #['30','44','70','100','143','217','353','545','857']
fitsend = '.fits' #'_pr2.fits'
fits = [None]*nmaps
for i in xrange(len(fits)):
    fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

outdir = fitsdir
outroot = 'wmap_deconv_smoothw_extrapolated_9yr_' #'planck_deconv_'
outcode = ['k','ka','q','v','w'] #fitscode
scal_outfits = [None]*nmaps
wav_outfits = [None]*nmaps
for i in xrange(len(scal_outfits)):
    scal_outfits[i] = outdir + outroot + outcode[i] + '_scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_outfits[i] = outdir + outroot + outcode[i] + '_wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'

#Load CMB maps
maps = hp.read_map(fits[0])
for i in xrange(1,len(fits)):
    maps = np.vstack((maps,hp.read_map(fits[i])))

#Calculate band-limited alms
print "\nBand-limiting input maps"
pixrecip = np.reciprocal(hp.pixwin(hp.get_nside(maps))[0:ellmax]) #pixwin
'''alms = np.array(hp.map2alm(maps,lmax=ellmax-1,pol=False)) #PARALLELISE & correct for pixellisation? #USE_WEIGHTS
    for i in xrange(len(alms)):
    hp.almxfl(alms[i],pixrecip,inplace=True)'''
pool = mg.Pool(nprocess)
alms = np.array(pool.map(almworker,maps))

#Calculate wavelet and scaling function maps for each channel
scal_maps = [None]*len(alms)
wav_maps = [None]*len(alms)
for i in xrange(len(alms)): #PARALELLISE
    print "Calculating scaling function and wavelet maps for input map", i
    wav_maps[i],scal_maps[i] = ps.analysis_lm2wav(alms[i],wavparam,ellmax,jmin,ndir,spin,upsample)
    np.save(scal_outfits[i],scal_maps[i])
    np.save(wav_outfits[i],wav_maps[i])
scal_maps = np.array(scal_maps)
wav_maps = np.array(wav_maps)

