import numpy as np
import healpy as hp
import math as mh

'''n1map = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon.fits')
inpmap = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map.fits')

n1cls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls.fits')[:3400]'''
inpcls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_cls_ilcbeam_lmax3399.fits')[:3400]*1.e-12 #Convert to K

'''pixrecip = 1. / hp.pixwin(2048)[:3400]
ilcbeam = hp.gauss_beam(mh.radians(5./60.),lmax=3399)
n1alms = hp.map2alm(n1map,lmax=3399)
hp.almxfl(n1alms,pixrecip,inplace=True)
inpalms = hp.map2alm(inpmap,lmax=3399)
hp.almxfl(inpalms,ilcbeam*pixrecip,inplace=True)

residalms = n1alms - inpalms

biascls = hp.alm2cl(inpalms,alms2=residalms)
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_bias.fits',biascls)

debiascls = n1cls - (2.*biascls)
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_debias.fits',debiascls)'''

biascls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_bias.fits')
fracbias = biascls / inpcls
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_fracbias.fits',fracbias)

'''residcls = hp.alm2cl(residalms)
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_resid.fits',residcls)

deresidcls = n1cls - (2.*biascls) - residcls
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_deresid.fits',deresidcls)

#Check results
print np.sum(n1cls[2:] - inpcls[2:] - residcls[2:] - (2.*biascls[2:]))
print np.sum(n1cls[2:])'''
print np.mean(fracbias[2:1000])
print np.mean(fracbias[2:])