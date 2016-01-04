import numpy as np
import healpy as hp
import math as mh
import copy as cp

if __name__ == "__main__":
    niter = 20000
    old_map = hp.read_map('/Users/kwame/Documents/s2let_ilc_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_1_recon_inpaint800.fits')
    mask = np.load('/Users/kwame/Documents/s2let_ilc_data/nilc_pr1_builtmask_holes_ring_800.npy') #0 where holes
    new_map = cp.deepcopy(old_map)
    #nside = hp.get_nside(mask)
    #hole_pixs = np.where(mask == 0)[0]
    nside = 2048
    hole_pixs = mask[0]
    '''theta = mh.pi / 2.
    phi = 0.
    rad_arcmin = 21
    hole_pixs = hp.query_disc(nside,hp.ang2vec(theta,phi),np.radians(rad_arcmin/60.))'''
    hole_neighbours = hp.get_all_neighbours(nside,hole_pixs) #8 x N_pix_holes

    new_map[hole_pixs] = np.mean(old_map) #Start with mean map value
    for i in xrange(niter):
        print "Iteration no.", i+1, "/", niter
        new_map[hole_pixs] = np.mean(new_map[hole_neighbours],axis=0)

    resid_map = new_map - old_map

    hp.write_map('/Users/kwame/Documents/s2let_ilc_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_1_recon_inpaint800GD20000.fits',new_map)