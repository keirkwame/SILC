import healpy as hp
import spin_silc_utilities as su

if __name__ == "__main__":
    ellmax = 2301 #L
    mask_name = 'fullSky'
    mask_fn = '/Users/keir/Documents/spin_silc/masks/full_sky_nside2048.fits' #RING ordering
    outfits_root = '/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_lmax2300' #Up to end of recon
    final_maps = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map.fits',field=(1,2))
    
    cl_output = su.masked_cl(ellmax,mask_name,mask_fn,outfits_root,final_maps)