import healpy as hp
import spin_silc_utilities as su

if __name__ == "__main__":
    ellmax = 917
    mask_name = 'PSfsky'
    mask_fn = '/Users/keir/Documents/spin_silc/masks/planck_pol_PSmask_nside512_05thresh.fits' #RING ordering
    outfits_root = '/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noN_planck_pol_deconv_917_hybridTestD_6_1_recon'
    final_maps = hp.read_map('/Users/keir/Documents/spin_silc/recon_maps/spin_silc_noN_planck_pol_deconv_917_hybridTestD_6_1_recon_QUmaps.fits',field=(1,2))

    cl_output = su.masked_cl(ellmax,mask_name,mask_fn,outfits_root,final_maps)