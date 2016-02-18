import numpy as np
import healpy as hp
import copy as cp

if __name__ == "__main__":
    nside = 512
    niter = 2000
    
    old_map = np.zeros((2,hp.nside2npix(nside))) #Final maps
    old_map[0],old_map[1] = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_fwhm50_planck_pol_diffusePS_deconv_917_hybridTestD_6_1_recon_QUmaps.fits',field=(1,2)) #(Q,U)

    '''old_map = np.zeros((4,hp.nside2npix(nside))) #N_side=1024 maps all at once
    old_map[0],old_map[1] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/LFI_SkyMap_030_1024_R2.01_full.fits',field=(1,2)) #(Q,U)
    old_map[2],old_map[3] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/LFI_SkyMap_044_1024_R2.01_full.fits',field=(1,2)) #(Q,U)'''
    
    '''old_map = np.zeros((10,hp.nside2npix(nside))) #N_side=2048 maps all at once
    old_map[0],old_map[1] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/LFI_SkyMap_070_2048_R2.01_full.fits',field=(1,2)) #(Q,U)
    old_map[2],old_map[3] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_100_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    old_map[4],old_map[5] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_143_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    old_map[6],old_map[7] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_217_2048_R2.02_full.fits',field=(1,2)) #(Q,U)
    old_map[8],old_map[9] = hp.read_map('/Users/keir/Documents/spin_silc/maps/PR2/frequencyMaps/HFI_SkyMap_353_2048_R2.02_full.fits',field=(1,2)) #(Q,U)'''

    #mask = hp.reorder(hp.read_map('/Users/keir/Documents/spin_silc/masks/planck_pol_PSmask_nside2048.fits'),n2r=True) #0 where holes
    mask = hp.read_map('/Users/keir/Documents/spin_silc/masks/planck_pol_PSmask_nside512_05thresh.fits') #0 where holes
    
    new_map = cp.deepcopy(old_map)
    hole_pixs = np.where(mask == 0)[0]
    hole_neighbours = hp.get_all_neighbours(nside,hole_pixs) #8 x N_pix_holes

    new_map[:,hole_pixs] = np.mean(old_map,axis=-1)[:,None] #Start with mean map value
    for i in xrange(niter):
        print "Iteration no.", i+1, "/", niter
        new_map[:,hole_pixs] = np.mean(new_map[:,hole_neighbours],axis=1) #Average over neighbouring pixs

    resid_map = new_map - old_map
    T_map = np.zeros_like(new_map[0])

    hp.write_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_fwhm50_planck_pol_diffusePS_deconv_917_hybridTestD_6_1_recon_QUmaps_diffPaint.fits',(T_map,new_map[0],new_map[1]))

    '''hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_30_pr2_nside1024.fits',(T_map,new_map[0],new_map[1]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_44_pr2_nside1024.fits',(T_map,new_map[2],new_map[3]))'''

    '''hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_70_pr2_nside2048.fits',(T_map,new_map[0],new_map[1]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_100_pr2_nside2048.fits',(T_map,new_map[2],new_map[3]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_143_pr2_nside2048.fits',(T_map,new_map[4],new_map[5]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_217_pr2_nside2048.fits',(T_map,new_map[6],new_map[7]))
    hp.write_map('/Users/keir/Documents/spin_silc/maps/PR2/diffuseMaps/planck_pol_diffuse_353_pr2_nside2048.fits',(T_map,new_map[8],new_map[9]))'''