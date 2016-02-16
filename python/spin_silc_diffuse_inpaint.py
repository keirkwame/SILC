import numpy as np
import healpy as hp
import copy as cp

if __name__ == "__main__":
    nside = 2048
    niter = 2000
    
    old_map = np.zeros((8,hp.nside2npix(nside))) #HFI all at once
    old_map[0],old_map[1] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_100_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    old_map[2],old_map[3] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_143_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    old_map[4],old_map[5] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_217_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    old_map[6],old_map[7] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_353_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    #mask = hp.reorder(hp.read_map('/Users/keir/Documents/spin_silc/masks/nilc_pol_conmask.fits'),n2r=True) #0 where holes
    mask = hp.read_map('/Users/keir/Documents/spin_silc/masks/nilc_pol_conmask_nside2048.fits') #0 where holes
    
    new_map = cp.deepcopy(old_map)
    hole_pixs = np.where(mask == 0)[0]
    hole_neighbours = hp.get_all_neighbours(nside,hole_pixs) #8 x N_pix_holes

    new_map[:,hole_pixs] = np.mean(old_map,axis=-1)[:,None] #Start with mean map value
    for i in xrange(niter):
        print "Iteration no.", i+1, "/", niter
        new_map[:,hole_pixs] = np.mean(new_map[:,hole_neighbours],axis=1) #Average over neighbouring pixs

    resid_map = new_map - old_map
    T_map = np.zeros_like(new_map[0])

    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_100_pr2_nside2048.fits',(T_map,new_map[0],new_map[1]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_143_pr2_nside2048.fits',(T_map,new_map[2],new_map[3]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_217_pr2_nside2048.fits',(T_map,new_map[4],new_map[5]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_353_pr2_nside2048.fits',(T_map,new_map[6],new_map[7]))