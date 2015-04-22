import numpy as np
import healpy as hp
import math as mh
import copy as cp

def gauss_source(chanmap,peaktemp,chansig,centre_ang,samppixs):
    sampangs = hp.pix2ang(hp.get_nside(chanmap),samppixs) #theta,phi
    thetadiff2 = (sampangs[0] - centre_ang[0])**2
    phidiff2 = (sampangs[1] - centre_ang[1])**2
    gaussvals = peaktemp * np.exp((thetadiff2 + phidiff2) / (-2. * (chansig**2)))
    chanmap[samppixs] = chanmap[samppixs] - gaussvals
    return chanmap

nmaps = 5 #No. maps (WMAP = 5)

#Point source catalogue TXT file
psc = 'foreground_data/wmap_ptsrc_catalog_9yr_v5p1.txt'

#Frequency channel map FITS files
fitsdir = 'deconv_data/'
fitsroot = 'wmap_band_deconv_imap_r9_9yr_'
fitscode = ['K','Ka','Q','V','W']
fitsend = '_v5.fits'
fits = [None]*nmaps
for i in xrange(len(fits)):
    fits[i] = fitsdir + fitsroot + fitscode[i] + fitsend

#Output map FITS files
outdir = fitsdir
outroot = 'wmap_deconv_nosource_9yr_'
outcode = fitscode
outend = '.fits'
outfits = [None]*nmaps
residfits = [None]*nmaps
for i in xrange(len(outfits)):
    outfits[i] = outdir + outroot + outcode[i] + outend
    residfits[i] = outdir + outroot + outcode[i] + '_resid' + outend

#Flux conversion factors and FWHM's
fluxfacs = np.array([262.7,211.9,219.6,210.1,179.2]) #uK Jy^(-1)
fwhms = [0.88,0.66,0.51,0.35,0.22] #deg - sqrt(solid_angle)

#Load maps
maps = hp.read_map(fits[0])
for i in xrange(1,len(fits)):
    maps = np.vstack((maps,hp.read_map(fits[i])))
origmaps = cp.deepcopy(maps)

#Load catalogue
catalog = np.loadtxt(psc,usecols=(2,3,4,5,6,7,8)) #lon,lat,K,Ka,Q,V,W


coords = catalog[:,0:2] #lon,lat
coords[:,1] = 90. - coords[:,1] #phi,theta
coords = np.radians(coords) #Convert to radians
catalog[:,2:][catalog[:,2:] == -9.9] = 0. #Set sentinel values to 0
temps = catalog[:,2:] * fluxfacs[:,None].T * 0.001 #mK
sigmas = np.radians(fwhms) / (2.*mh.sqrt(2.*mh.log(2.))) #sigma in radians

samppixsall = np.array([0])
for i in xrange(len(coords)): #Loop over point sources
    for j in xrange(len(sigmas)): #Loop over maps
        print 'Subtracting Gaussian profile for point source', i, 'map', j
        samppixs = hp.query_disc(hp.get_nside(maps),hp.ang2vec(coords[i,1],coords[i,0]),5.*sigmas[j])
        maps[j] = gauss_source(maps[j],temps[i,j],sigmas[j],(coords[i,1],coords[i,0]),samppixs)
        #samppixsall = np.concatenate((samppixsall,samppixs)) #Slows program down

#Save point source residuals
resid = origmaps - maps

#Write new maps to FITS files
for i in xrange(len(maps)):
    hp.write_map(outfits[i],maps[i])
    hp.write_map(residfits[i],resid[i])

