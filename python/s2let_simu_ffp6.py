import numpy as np
import healpy as hp
import math as mh

##Input
fitsdir = '/home/keir/global/project/projectdirs/cmb/data/planck2013/ffp6/fiducial/components/' #'/Users/keir/Documents/s2let_ilc_planck/components/'
fitsroot = 'ffp6_'
fitscode = ['030','044','070','100','143','217','353','545','857']
fitsend = '_nominal_map.fits'

outdir = '/home/keir/s2let_ilc_data/ffp6_data_withPS/' #'/Users/keir/Documents/s2let_ilc_planck/ffp6_data/'
outroot = 'ffp6_fiducial_withPS_'
outcode = fitscode
outend = '.fits'

lfi = np.zeros((3*4,hp.nside2npix(1024))) #Pre-allocate arrays [list of 12 maps]
hfi = np.zeros((6*4,hp.nside2npix(2048))) #List of 24 maps

#Loading maps
for i in xrange(3):
    cmb_fits = fitsdir + fitsroot + 'cmb_lensed_' + fitscode[i] + fitsend
    lfi[i,:] = hp.read_map(cmb_fits)
    fg_fits = fitsdir + fitsroot + 'foreground_' + fitscode[i] + fitsend
    lfi[i+3,:] = hp.read_map(fg_fits)
    ps_fits = fitsdir + fitsroot + 'ps_' + fitscode[i] + fitsend
    lfi[i+6,:] = hp.read_map(ps_fits)
    noise_fits = fitsdir + fitsroot + 'noise_' + fitscode[i] + fitsend
    lfi[i+9,:] = hp.read_map(noise_fits)

for i in xrange(6):
    cmb_fits = fitsdir + fitsroot + 'cmb_lensed_' + fitscode[i+3] + fitsend
    hfi[i,:] = hp.read_map(cmb_fits)
    fg_fits = fitsdir + fitsroot + 'foreground_' + fitscode[i+3] + fitsend
    hfi[i+6,:] = hp.read_map(fg_fits)
    ps_fits = fitsdir + fitsroot + 'ps_' + fitscode[i+3] + fitsend
    hfi[i+12,:] = hp.read_map(ps_fits)
    noise_fits = fitsdir + fitsroot + 'noise_' + fitscode[i+3] + fitsend
    hfi[i+18,:] = hp.read_map(noise_fits)

#Filling in bad pixels - LFI
badval = -1.63750E+30
print "Filling in bad pixels in LFI\n"
badpixs = np.where(lfi == badval) #Output is tuple of 2 arrays
badvecs = hp.pix2vec(1024,badpixs[1]) #Output is tuple of 3 arrays
lfi_masked = np.ma.masked_values(lfi,badval) #Masking badvals
totnum = len(badvecs[0])
for i in xrange(len(badvecs[0])): #Loop over bad pixels
    print "Filling in bad pixel no.", i+1, "/", totnum
    baddisc = hp.query_disc(1024,(badvecs[0][i],badvecs[1][i],badvecs[2][i]),mh.radians(1.)) #array
    newval = np.ma.average(lfi_masked[badpixs[0][i],baddisc]) #Avg. over disc excluding badvals
    lfi[badpixs[0][i],badpixs[1][i]] = newval #Update original array with new value
#HFI
print "\nFilling in bad pixels in HFI\n"
badpixs2 = np.where(hfi == badval) #Output is tuple of 2 arrays
badvecs2 = hp.pix2vec(2048,badpixs2[1]) #Output is tuple of 3 arrays
hfi_masked = np.ma.masked_values(hfi,badval) #Masking badvals
totnum2 = len(badvecs2[0])
for i in xrange(len(badvecs2[0])): #Loop over bad pixels
    print "Filling in bad pixel no.", i+1, "/", totnum2
    baddisc = hp.query_disc(2048,(badvecs2[0][i],badvecs2[1][i],badvecs2[2][i]),mh.radians(1.))#arr
    newval = np.ma.average(hfi_masked[badpixs2[0][i],baddisc]) #Avg. over disc excluding badval
    hfi[badpixs2[0][i],badpixs2[1][i]] = newval #Update original array with new value

lfi_combined = lfi[:3,:] + lfi[3:6,:] + lfi[9:,:] #- lfi[6:9,:] #CMB + FG [astrophysical components & PS] + Noise - PS
hfi_combined = hfi[:6,:] + hfi[6:12,:] + hfi[18:,:] #- hfi[12:18,:]

#Saving maps
print "\nSaving output maps"
for i in xrange(3):
    outfits = outdir + outroot + outcode[i] + outend
    hp.write_map(outfits,lfi_combined[i,:])

for i in xrange(6):
    outfits = outdir + outroot + outcode[i+3] + outend
    hp.write_map(outfits,hfi_combined[i,:])


