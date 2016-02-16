import numpy as np
import healpy as hp
import math as mh
import multiprocessing as mg

def almworker(i):
    print "This is (map2alm & alm2map) worker starting for another map"
    QU_maps = hp.read_map(fits[i],field=(1,2)) #(Q,U)
    print "Calculating a_lm"
    alms_worker = (hp.map2alm(QU_maps[0],lmax=ellmax-1),hp.map2alm(QU_maps[1],lmax=ellmax-1)) #(Q,U)
    del QU_maps
    hp.almxfl(alms_worker[0],deconbeam[i],inplace=True) #Correcting for pixwin & smoothing #Q
    hp.almxfl(alms_worker[1],deconbeam[i],inplace=True) #Correcting for pixwin & smoothing #U
    print "Calculating maps"
    map_worker = (hp.alm2map(alms_worker[0],nside_out,pixwin=True),hp.alm2map(alms_worker[1],nside_out,pixwin=True)) #(Q,U)
    del alms_worker
    T_map = (np.zeros_like(map_worker[0]),)
    hp.write_map(outfits[i],T_map+map_worker)
    del T_map,map_worker
    return 0

def smoothfunc(x,borig): #x is range of ell's
    xmiddle = .5*(x[0]+x[-1])
    alpha = 3.5/(x[-1] - xmiddle)
    y = (borig * (np.exp(alpha*(x[0] - xmiddle)) + 1.))/(np.exp(alpha*(x - xmiddle)) + 1.)
    return y

if __name__ == "__main__":
    nprocess = 4
    nmaps = 7 #No. files
    ellmax = 2300
    nside_out = 2048

    #Frequency channel map FITS files
    fitsdir = '/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/' #frequencyMaps/'
    fitsprefix = ['LFI','LFI','LFI','HFI','HFI','HFI','HFI']
    fitsroot = 'planck_pol_diffuse_' #'_SkyMap_'
    fitscode = ['30','44','70','100','143','217','353'] #['030_1024_R2.01','044_1024_R2.01','070_2048_R2.01','100_2048_R2.02','143_2048_R2.02','217_2048_R2.02','353_2048_R2.02']
    fitsend = '_pr2_nside' #'_full.fits'
    fitssuffix = ['1024.fits','1024.fits','2048.fits','2048.fits','2048.fits','2048.fits','2048.fits']
    fits = [None]*nmaps
    for i in xrange(len(fits)):
        fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend + fitssuffix[i]

    #Beam transfer function Numpy files
    beamdir = '/Users/keir/Documents/spin_silc/beams/'
    beamroot = 'planck_bl_'
    beamcode = ['100','143','217','353']
    beamend = '_pr2.npy'
    txt = [None]*4
    for i in xrange(len(txt)):
        txt[i] = beamdir + beamroot + beamcode[i] + beamend

    #Output map FITS files
    outdir = '/Users/keir/Documents/spin_silc/maps/PR2/deconvolvedMaps/'
    outroot = 'planck_pol_diffuse_deconv_'
    outcode = ['30','44','70','100','143','217','353']
    outend = '_pr2.fits'
    outfits = [None]*nmaps
    for i in xrange(len(outfits)):
        outfits[i] = outdir + outroot + outcode[i] + outend

    #Load beam transfer functions
    rawbeam = [None]*nmaps
    fwhms = np.radians([32.33/60.,27.01/60.,13.25/60.]) #PR2 beams
    rawbeam[0] = hp.gauss_beam(fwhms[0],lmax=ellmax-1)
    rawbeam[1] = hp.gauss_beam(fwhms[1],lmax=ellmax-1)
    rawbeam[2] = hp.gauss_beam(fwhms[2],lmax=ellmax-1)
    for i in xrange(3,nmaps):
        rawbeam[i] = np.load(txt[i-3])[:ellmax]

    #For Planck, effective beam is 5' FWHM gaussian
    smoothbeam = hp.gauss_beam(mh.radians(5./60.),lmax=ellmax-1)
    '''smoothbegin = 2100
    smoothend = ellmax-1
    smoothbeam[smoothbegin:smoothend+1] = smoothfunc(np.arange(smoothbegin,smoothend+1),smoothbeam[smoothbegin])
    np.save('/Users/keir/Documents/spin_silc/beams/spin_silc_planck_bl_pr2.npy',smoothbeam)'''

    #Take reciprocal of beam transfer functions
    pixrecip1024 = np.reciprocal(hp.pixwin(1024)[:ellmax])
    pixrecip2048 = np.reciprocal(hp.pixwin(2048)[:ellmax])
    recipcombbeam = np.reciprocal(rawbeam)

    #Thresholding of deconvolved beams to 1/b_l <= 1000
    recipcombbeam[recipcombbeam > 1000.] = 1000.

    deconbeam = [None]*len(recipcombbeam)
    for i in xrange(2): #30 & 44 GHz
        deconbeam[i] = recipcombbeam[i] * smoothbeam * pixrecip1024 #Deconvolving to tapered 5' gaussian beam
    for i in xrange(2,len(recipcombbeam)): #70 to 353 GHz
        deconbeam[i] = recipcombbeam[i] * smoothbeam * pixrecip2048 #Deconvolving to tapered 5' gaussian beam

    #Parallelised SHT to and from harmonic space
    pool = mg.Pool(nprocess)
    deconmaps = pool.map(almworker,np.arange(nmaps))
    pool.close()
    pool.join()
