import numpy as np
import healpy as hp

##Input
fitsdir = '/Users/keir/Documents/s2let_ilc_planck/components/'
fitsroot = 'ffp6_'
fitscode = ['030','044','070','100','143','217','353','545','857']
fitsend = '_nominal_map.fits'

outdir = '/Users/keir/Documents/s2let_ilc_planck/ffp6_data/'
outroot = 'ffp6_fiducial_noPS_'
outcode = fitscode
outend = '.fits'

cmb_lfi = np.zeros((3,hp.nside2npix(1024))) #Pre-allocate arrays
fg_lfi = np.zeros((3,hp.nside2npix(1024)))
ps_lfi = np.zeros((3,hp.nside2npix(1024)))
noise_lfi = np.zeros((3,hp.nside2npix(1024)))
cmb_hfi = np.zeros((6,hp.nside2npix(2048)))
fg_hfi = np.zeros((6,hp.nside2npix(2048)))
ps_hfi = np.zeros((6,hp.nside2npix(2048)))
noise_hfi = np.zeros((6,hp.nside2npix(2048)))

#Loading maps
for i in xrange(3):
    cmb_fits = fitsdir + fitsroot + 'cmb_lensed_' + fitscode[i] + fitsend
    cmb_lfi[i,:] = hp.read_map(cmb_fits)
    fg_fits = fitsdir + fitsroot + 'foreground_' + fitscode[i] + fitsend
    fg_lfi[i,:] = hp.read_map(cmb_fits)
    ps_fits = fitsdir + fitsroot + 'ps_' + fitscode[i] + fitsend
    ps_lfi[i,:] = hp.read_map(cmb_fits)
    noise_fits = fitsdir + fitsroot + 'noise_' + fitscode[i] + fitsend
    noise_lfi[i,:] = hp.read_map(cmb_fits)

for i in xrange(6):
    cmb_fits = fitsdir + fitsroot + 'cmb_lensed_' + fitscode[i+3] + fitsend
    cmb_hfi[i,:] = hp.read_map(cmb_fits)
    fg_fits = fitsdir + fitsroot + 'foreground_' + fitscode[i+3] + fitsend
    fg_hfi[i,:] = hp.read_map(cmb_fits)
    ps_fits = fitsdir + fitsroot + 'ps_' + fitscode[i+3] + fitsend
    ps_hfi[i,:] = hp.read_map(cmb_fits)
    noise_fits = fitsdir + fitsroot + 'noise_' + fitscode[i+3] + fitsend
    noise_hfi[i,:] = hp.read_map(cmb_fits)

lfi = cmb_lfi + fg_lfi + noise_lfi - ps_lfi #CMB + FG [astrophysical components & PS] + Noise - PS
hfi = cmb_hfi + fg_hfi + noise_hfi - ps_hfi

#Saving maps
for i in xrange(3):
    outfits = outdir + outroot + outcode[i] + outend
    hp.write_map(outfits,lfi[i,:])

for i in xrange(6):
    outfits = outdir + outroot + outcode[i+3] + outend
    hp.write_map(outfits,hfi[i,:])


