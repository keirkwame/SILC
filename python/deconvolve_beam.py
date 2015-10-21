import numpy as np
import matplotlib as mp
import healpy as hp
import math as mh
import multiprocessing as mg
import copy as cp

def almworker(i):
    print "This is (map2alm & alm2map) worker starting for another map"
    alms_worker = hp.map2alm(i[0],lmax=ellmax-1)
    hp.almxfl(alms_worker,i[1],inplace=True) #Correcting for pixwin & smoothing
    map_worker = hp.alm2map(alms_worker,nside_out,pixwin=True)
    del alms_worker
    return map_worker

def smoothfunc(x,borig): #x is range of ell's
    xmiddle = .5*(x[0]+x[-1])
    alpha = 3.5/(x[-1] - xmiddle) #numer=5 for WMAP #numer=3 for Planck
    y = (borig * (np.exp(alpha*(x[0] - xmiddle)) + 1.))/(np.exp(alpha*(x - xmiddle)) + 1.)
    return y

if __name__ == "__main__":
    nprocess = 9
    nmaps = 9 #No. maps (WMAP = 5) (Planck = 9)
    nda = 9 #No. differencing assemblies (WMAP = 10) (Planck = [effectively] 9)
    ellmax = 3600 #Max. is 750 due to WMAP beams
    nside_out = 2048

    #Frequency channel map FITS files
    fitsdir = '/home/keir/s2let_ilc_data/maps/' #'/Users/keir/Documents/s2let_ilc_planck/maps/PR2/frequencyMaps/'
    fitsprefix = ['LFI','LFI','LFI','HFI','HFI','HFI','HFI','HFI','HFI'] #['030/LFI','044/LFI','070/LFI','100/HFI','143/HFI','217/HFI','353/HFI','545/HFI','857/HFI']
    fitsroot = '_SkyMap_'
    fitscode = ['030_1024_R2.01','044_1024_R2.01','070_2048_R2.01','100_2048_R2.00','143_2048_R2.00','217_2048_R2.00','353_2048_R2.00','545_2048_R2.00','857_2048_R2.00'] #['030-ffp6_1024','044-ffp6_1024','070-ffp6_1024','100-ffp6_2048','143-ffp6_2048','217-ffp6_2048','353-ffp6_2048','545-ffp6_2048','857-ffp6_2048']
    fitsend = '_full.fits'
    fits = [None]*nmaps
    for i in xrange(len(fits)):
        fits[i] = fitsdir + fitsprefix[i] + fitsroot + fitscode[i] + fitsend

    #WMAP beam transfer function TXT files
    beamdir = '/home/keir/s2let_ilc_data/beams/' #'/Users/keir/Documents/s2let_ilc_planck/beams/'
    beamroot = 'planck_bl_'
    beamcode = ['30','44','70','100','143','217','353','545','857']
    beamend = '_pr2.npy'
    txt = [None]*nda
    for i in xrange(len(txt)):
        txt[i] = beamdir + beamroot + beamcode[i] + beamend

    #Output map FITS files
    outdir = '/home/keir/s2let_ilc_data/hybrid_data/' #'/Users/keir/Documents/s2let_ilc_planck/hybrid_data/'
    outroot = 'planck_deconv_tapered_thresh_lmax3600_'
    outcode = beamcode
    outend = '_pr2.fits'
    outfits = [None]*nmaps
    for i in xrange(len(outfits)):
        outfits[i] = outdir + outroot + outcode[i] + outend

    #Pixel noise coefficients
    '''beam_weights_q = [1./2.245,1./2.131]
    beam_weights_v = [1./3.314,1./2.949]
    beam_weights_w = [1./5.899,1./6.562,1./6.941,1./6.773] #Inverse to sigma_0 parameters
    fwhms = [0.88,0.66,0.51,0.35,0.22] #deg - square root of beam solid angle'''

    #Load maps
    maps = [None]*nmaps
    for i in xrange(len(fits)):
        maps[i] = hp.read_map(fits[i])

    #Unit conversions for 545 & 857 GHz - not necessary for Fiducial/MC simulations/minusgaussps
    maps[-2] = maps[-2] / 58.0356
    maps[-1] = maps[-1] / 2.2681

    #Load beam transfer functions
    rawbeam = [None]*nda
    for i in xrange(nda):
        rawbeam[i] = np.load(txt[i]) #np.loadtxt(txt[i],usecols=[1])
        rawbeam[i] = rawbeam[i][:ellmax] #Truncate beam transfer functions to lmax

    #Calculating beam transfer functions for each channel - specific to WMAP
    '''combbeam = [None]*(nmaps)
    combbeam[0] = rawbeam[0] #K
    combbeam[1] = rawbeam[1] #Ka
    combbeam[2] = np.average(np.array([rawbeam[2],rawbeam[3]]),axis=0,weights=beam_weights_q) #Weighted Q
    combbeam[3] = np.average(np.array([rawbeam[4],rawbeam[5]]),axis=0,weights=beam_weights_v) #Weighted V
    combbeam[4] = np.average(np.array([rawbeam[6],rawbeam[7],rawbeam[8],rawbeam[9]]),axis=0,weights=beam_weights_w) #Weighted W'''

    #Extrapolate K & Ka (for l_max = 1024) comb-beams by Gaussian approximation
    '''kgauss = hp.gauss_beam(mh.radians(fwhms[0]),lmax=ellmax)
    kgaussnorm = kgauss / kgauss[1] #Normalise to b_1 = 1
    kagauss = hp.gauss_beam(mh.radians(fwhms[1]),lmax=ellmax)
    kagaussnorm = kagauss / kagauss[1] #Normalise to b_1 = 1
    kgaussextrap = kgauss[751:] * (combbeam[0][-1]/kgauss[750]) #Gaussian shifted up to join real K-beam
    kagaussextrap = kagauss[851:] * (combbeam[1][-1]/kagauss[850]) #Gaussian shifted to join real Ka-beam
    combbeam[0] = np.concatenate((combbeam[0],kgaussextrap))
    combbeam[1] = np.concatenate((combbeam[1],kagaussextrap))
    combbeam = np.array(combbeam)'''
    
    #Extrapolate 30, 44 & 70 GHz (for l_max >= 2048) raw-beams by Gaussian approximation
    fwhms = [32.33/60.,27.01/60.,13.25/60.] #[32.24/60.,27.00/60.,13.25/60.] #PR2/PR1
    '''gauss30 = hp.gauss_beam(mh.radians(fwhms[0]),lmax=ellmax-1)
    gauss30norm = gauss30 / gauss30[1] #Normalise to b_1 = 1
    gauss44 = hp.gauss_beam(mh.radians(fwhms[1]),lmax=ellmax-1)
    gauss44norm = gauss44 / gauss44[1] #Normalise to b_1 = 1
    gauss70 = hp.gauss_beam(mh.radians(fwhms[2]),lmax=ellmax-1)
    gauss70norm = gauss70 / gauss70[1] #Normalise to b_1 = 1
    gauss30extrap = gauss30norm[2049:] * (rawbeam[0][-1]/gauss30norm[2048]) #Gauss shift to join 30-beam
    gauss44extrap = gauss44norm[2049:] * (rawbeam[1][-1]/gauss44norm[2048]) #Gauss shift to join 44-beam
    gauss70extrap = gauss70norm[2049:] * (rawbeam[2][-1]/gauss70norm[2048]) #Gauss shift to join 70-beam
    rawbeam[0] = np.concatenate((rawbeam[0],gauss30extrap))
    rawbeam[1] = np.concatenate((rawbeam[1],gauss44extrap))
    rawbeam[2] = np.concatenate((rawbeam[2],gauss70extrap))'''
    #Try replacing LFI beams with Gaussian approximations
    rawbeam[0] = hp.gauss_beam(mh.radians(fwhms[0]),lmax=ellmax-1)
    rawbeam[1] = hp.gauss_beam(mh.radians(fwhms[1]),lmax=ellmax-1)
    rawbeam[2] = hp.gauss_beam(mh.radians(fwhms[2]),lmax=ellmax-1)
    rawbeam = np.array(rawbeam)

    #Smoothing W-beam for 're-convolving'
    '''smoothbeam = cp.deepcopy(combbeam[4])
    smoothbegin = 875
    smoothend = 1024
    smoothbeam[smoothbegin:smoothend+1] = smoothfunc(np.arange(smoothbegin,smoothend+1),smoothbeam[smoothbegin])'''
    #For Planck, effective beam is 5' FWHM gaussian
    smoothbeam = hp.gauss_beam(mh.radians(5./60.),lmax=ellmax-1)
    smoothbegin = 3400
    smoothend = ellmax-1
    smoothbeam[smoothbegin:smoothend+1] = smoothfunc(np.arange(smoothbegin,smoothend+1),smoothbeam[smoothbegin])
    smoothbeam_orig = hp.gauss_beam(mh.radians(5./60.),lmax=ellmax-1)

    #Take reciprocal of beam transfer functions
    #pixrecip = np.reciprocal(hp.pixwin(hp.get_nside(maps))[:ellmax+1]) #TESTING pixwin
    pixrecip1024 = np.reciprocal(hp.pixwin(1024)[:ellmax]) #30 & 44 GHz at Nside=1024
    pixrecip2048 = np.reciprocal(hp.pixwin(2048)[:ellmax]) #70 & HFI @ Nside=2048 #Not 70 for MC
    recipcombbeam = np.reciprocal(rawbeam) #combbeam)

    #Thresholding of deconvolved beams to 1/b_l <= 1000
    recipcombbeam[recipcombbeam > 1000] = 1000

    deconbeam = [None]*len(recipcombbeam)
    for i in xrange(2): #len(recipcombbeam)):
        deconbeam[i] = recipcombbeam[i] * smoothbeam * pixrecip1024 #Deconvolving to tapered 5' gaussian beam
    for i in xrange(2,len(recipcombbeam)):
        deconbeam[i] = recipcombbeam[i] * smoothbeam * pixrecip2048 #Deconvolving to tapered 5' gaussian beam
    deconbeam = np.array(deconbeam)

    #Parallelised SHT to and from harmonic space
    pool = mg.Pool(nprocess)
    mapsextra = [(maps[i],deconbeam[i]) for i in xrange(len(maps))]
    del maps
    deconmaps = pool.map(almworker,mapsextra)
    pool.close()
    pool.join()

    #Write new maps to FITS files
    for i in xrange(len(deconmaps)):
        hp.write_map(outfits[i],deconmaps[i])


