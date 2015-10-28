import numpy as np
import healpy as hp
import math as mh
import copy as cp

if __name__ == "__main__":
    niter = 2000
    old_map = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/hybrid_data/s2let_ilc_planck_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon.fits')
    mask = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits') #0 where holes
    new_map = cp.deepcopy(old_map)
    nside = hp.get_nside(mask)
    hole_pixs = np.where(mask == 0)[0]
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

    hp.write_map('/Users/keir/Documents/s2let_ilc_planck/hybrid_data/s2let_ilc_planck_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_diffuse.fits',new_map)