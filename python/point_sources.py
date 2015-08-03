import numpy as np
import healpy as hp
import math as mh
import copy as cp
import astropy.io.fits as af

def gauss_source(chanmap,peaktemp,chansig,centre_ang,samppixs):
    sampangs = hp.pix2ang(hp.get_nside(chanmap),samppixs) #theta,phi
    thetadiff2 = (sampangs[0] - centre_ang[0])**2
    phidiff2 = (sampangs[1] - centre_ang[1])**2
    gaussvals = peaktemp * np.exp((thetadiff2 + phidiff2) / (-2. * (chansig**2)))
    chanmap[samppixs] = chanmap[samppixs] - gaussvals
    return chanmap

nmaps = 9 #No. maps (WMAP = 5) (Planck = 9)

#Point source catalogue TXT file
#psc = 'foreground_data/wmap_ptsrc_catalog_9yr_v5p1.txt'
psc = [None]*nmaps
psc[0] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_030_R1.30.fits'
psc[1] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_044_R1.30.fits'
psc[2] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_070_R1.30.fits'
psc[3] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_100_R1.20.fits'
psc[4] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_143_R1.20.fits'
psc[5] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_217_R1.20.fits'
psc[6] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_353_R1.20.fits'
psc[7] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_545_R1.20.fits'
psc[8] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_857_R1.20.fits'

#Frequency channel map FITS files
fitsdir = '/Users/keir/Documents/s2let_ilc_planck/maps/PR1/' #'deconv_data/'
fitsroot = 'planck_filled_' #'wmap_band_deconv_imap_r9_9yr_'
fitscode = ['30','44','70','100','143','217','353','545','857'] #['K','Ka','Q','V','W']
fitsend = '_pr1.fits' #'_v5.fits'
fits = [None]*nmaps
for i in xrange(len(fits)):
    fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

#Output map FITS files
outdir = fitsdir
outroot = 'planck_filled_minusgaussps_' #'wmap_deconv_nosource_9yr_'
outcode = fitscode
outend = '_pr1.fits' #'.fits'
outfits = [None]*nmaps
residfits = [None]*nmaps
for i in xrange(len(outfits)):
    outfits[i] = outdir + outroot + outcode[i] + outend
    residfits[i] = outdir + outroot + outcode[i] + '_resid' + outend

#Flux conversion factors and FWHM's
fluxfacs = np.array([1./23.5099, 1./55.7349 ,1./129.1869, 1./244.0960, 1./371.7327, 1./483.6874, 1./287.4517, 1./58.0356, 1./2.2681]) #K_CMB / (MJy/sr) #np.array([262.7,211.9,219.6,210.1,179.2]) #uK Jy^(-1)
facs = (4. * mh.log(2.) * 1.e-9 * 10800 * 10800 * fluxfacs) / (mh.pi**3)
fwhms = [32.38,27.10,13.30,9.88,7.18,4.87,4.65,4.72,4.39] #arcmin (Planck) #[0.88,0.66,0.51,0.35,0.22] #deg - sqrt(solid_angle) #Planck - each source also has own eff. Gauss. FWHM

#Load maps
maps = [None]*nmaps
resid = [None]*nmaps
for i in xrange(nmaps):
    maps[i] = hp.read_map(fits[i])
#Unit conversions for 545 & 857 GHz - not necessary for Fiducial/MC simulations
maps[-2] = maps[-2] / 58.0356
maps[-1] = maps[-1] / 2.2681
origmaps = cp.deepcopy(maps)

samppixsall = np.array([0])

for j in xrange(nmaps): #len(sigmas)): #Loop over maps
    #Load catalogue
    full_catalog = af.open(psc[j])[1].data
    catalog =  np.array([full_catalog['GLON'],full_catalog['GLAT'],full_catalog['DETFLUX'],full_catalog['GAUFLUX'],full_catalog['GAU_FWHM_EFF'],full_catalog['EXTENDED']]).T #GLON,GLAT,mJy,mJy,arcmin,flag #np.loadtxt(psc,usecols=(2,3,4,5,6,7,8)) #lon,lat,K,Ka,Q,V,W

    coords = catalog[:,:2] #lon,lat
    coords[:,1] = 90. - coords[:,1] #phi,theta
    coords = np.radians(coords) #Convert to radians
    catalog[:,3][np.isnan(catalog[:,3])] = 0. #Set sentinel values (WMAP = -9.9, Planck = nan) to 0
    catalog[:,3][catalog[:,3] < 0.] = 0. #Throw out spurious negative fluxes
    catalog[:,4][np.isnan(catalog[:,4])] = 1.
    temps = (catalog[:,2] * facs[j]) / np.square(fwhms[j]) #catalog[:,3]) #Using DETFLUX
    temps[catalog[:,-1] == 1] = (catalog[:,3][catalog[:,-1] == 1] * facs[j]) / np.square(catalog[:,4][catalog[:,-1] == 1]) #Using GAUFLUX if EXTENDED
    #catalog[:,2:] * fluxfacs[:,None].T * 0.001 #mK
    sigmas = np.radians(catalog[:,4] / 60.) / (2.*mh.sqrt(2.*mh.log(2.))) #np.radians(fwhms) / (2.*mh.sqrt(2.*mh.log(2.))) #sigma in radians
    #sigmas = np.zeros_like(temps)
    sigmas[catalog[:,-1] == 0] = np.radians(fwhms[j] / 60.) / (2.*mh.sqrt(2.*mh.log(2.))) #Use beamFWHM if not EXTENDED

    for i in xrange(len(coords)): #Loop over point sources
        print 'Subtracting Gaussian profile for point source', i+1, '/', len(coords), 'map', j
        samppixs = hp.query_disc(hp.get_nside(maps[j]),hp.ang2vec(coords[i,1],coords[i,0]),5.*sigmas[i])
        maps[j] = gauss_source(maps[j],temps[i],sigmas[i],(coords[i,1],coords[i,0]),samppixs)
        #samppixsall = np.concatenate((samppixsall,samppixs)) #Slows program down
        
    #Save point source residuals
    resid[j] = origmaps[j] - maps[j]

'''testmap = cp.deepcopy(maps[j])
testmap[:] = 0.
testmap[samppixsall[1:]] = 1.'''

#Write new maps to FITS files
for i in xrange(len(maps)):
    hp.write_map(outfits[i],maps[i])
    hp.write_map(residfits[i],resid[i])

