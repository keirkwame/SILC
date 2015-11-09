import numpy as np
import healpy as hp
import math as mh
import copy as cp
import matplotlib.pyplot as plt

#Latex font settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16.0)

#Plot N=1 and NILC maps
#n1map = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon.fits') * 1e6 #uK
#nilcmap_masked = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map.fits') * 1e6 #uK
'''psmask = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits')
n1map_masked = cp.deepcopy(n1map)
n1map_masked[psmask == 0] = 0. #None
nilcmap_masked = cp.deepcopy(nilcmap)
nilcmap_masked[psmask == 0] = 0. #None'''

ffp8map = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar1_ffp8_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon.fits') * 1e6 #uK
#inputmap = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map.fits') * 1e6 #uK
'''newmap = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon.fits')
psmask = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits')
newmap_masked = cp.deepcopy(newmap)
newmap_masked[psmask == 0] = 0. #None'''
'''ffp8map_masked = cp.deepcopy(ffp8map)
ffp8map_masked[psmask == 0] = 0. #None
inputmap_masked = cp.deepcopy(inputmap)
inputmap_masked[psmask == 0] = 0. #None'''

'''LIM = 300
hp.mollview(n1map_masked,unit=r'$\mu\mathrm{K}$',title=r'Directional ILC $(N=1)$',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/n1map.pdf')

hp.mollview(nilcmap_masked,unit=r'$\mu\mathrm{K}$',title=r'NILC',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/nilcmap.pdf')'''

#Plot difference map
'''n1alms = hp.map2alm(n1map,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(n1map))[:257]
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
hp.almxfl(n1alms,pixrecip*smoothbeam,inplace=True)
n1map_downgrade = hp.alm2map(n1alms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/hybrid_data/s2let_ilc_covar10_planck_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_lmax256_lmin2.fits',n1map_downgrade)'''

'''nilcalms = hp.map2alm(nilcmap,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(nilcmap))[:257]
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
hp.almxfl(nilcalms,pixrecip*smoothbeam,inplace=True)
nilcmap_downgrade = hp.alm2map(nilcalms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr2_lmax256.fits',nilcmap_downgrade)'''

ffp8alms = hp.map2alm(ffp8map,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(ffp8map))[:257]
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
hp.almxfl(ffp8alms,pixrecip*smoothbeam,inplace=True)
ffp8map_downgrade = hp.alm2map(ffp8alms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar1_ffp8_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_lmax256_lmin2.fits',ffp8map_downgrade)

'''inputalms = hp.map2alm(inputmap,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(inputmap))[:257]
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=256) #5 arcmin
hp.almxfl(inputalms,pixrecip*smoothbeam*ilcbeam,inplace=True)
inputmap_downgrade = hp.alm2map(inputalms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map_ilcbeam_lmax256_lmin2.fits',inputmap_downgrade)'''

#n1map_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_fp8_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_lmax256.fits')
#nilcmap_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr2_lmax256.fits')

#ffp8map_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_lmax256_lmin2.fits')
inputmap_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map_ilcbeam_lmax256_lmin2.fits')
ut78 = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr2_UT78_128_cut05.fits')
ffp8map_downgrade[ut78 == 0] = None
inputmap_downgrade[ut78 == 0] = None

'''LIM = 300
hp.mollview(n1map_downgrade,unit=r'$\mu\mathrm{K}$',title=r'Directional ILC $(N=1)$',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/n1map_dg.pdf')

hp.mollview(nilcmap_downgrade,unit=r'$\mu\mathrm{K}$',title=r'NILC',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/nilcmap_dg.pdf')'''

'''diffmap = n1map_downgrade - nilcmap_downgrade
LIM = 20
hp.mollview(diffmap,unit=r'$\mu\mathrm{K}$',title=r'Directional ILC $(N=1)$ - NILC',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/diffmap.pdf')'''

diffmap = ffp8map_downgrade - inputmap_downgrade
LIM = 15
hp.mollview(diffmap,unit=r'$\mu\mathrm{K}$',title=r'Directional ILC $(N=1)$ - Input',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/diffmap_covar1diffuse_ffp8masked15.pdf')

#Plot N=1 and NILC spectra
#ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=299) #5 arcmin

'''n1cls_masked = hp.anafast(n1map_masked,lmax=3599)
pixrecip = 1. / hp.pixwin(hp.get_nside(n1map))[:3600]
fsky = 1. - (float(np.sum(np.logical_not(psmask))) / len(psmask))
print fsky
n1cls_corrected = (n1cls_masked * pixrecip * pixrecip) / fsky
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_ffp8_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_cls_masked.fits',n1cls_corrected)'''

'''nilccls_masked = hp.anafast(nilcmap_masked,lmax=299)
pixrecip = 1. / hp.pixwin(hp.get_nside(nilcmap_masked))[:300]'''
'''fsky = 1. - (float(np.sum(np.logical_not(psmask))) / len(psmask))
print fsky'''
'''nilccls_corrected = nilccls_masked * pixrecip * pixrecip * ilcbeam * ilcbeam
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_cls_ilcbeam_lmax2999.fits',nilccls_corrected)'''

'''theorycls = hp.read_cl('/Users/keir/Software/camb/planck2015_2_scalCls.fits')[0][:3600] * 1e12 #(uK)^2
ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=3599) #5 arcmin
theorycls_corrected = theorycls * ilcbeam * ilcbeam'''

'''newcls_masked = hp.anafast(newmap_masked,lmax=299)
pixrecip = 1. / hp.pixwin(hp.get_nside(newmap))[:300]
fsky = 1. - (float(np.sum(np.logical_not(psmask))) / len(psmask))
print fsky
newcls_corrected = (newcls_masked * pixrecip * pixrecip) / fsky
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_cls_masked.fits',newcls_corrected)'''

#newcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_cls.fits')[:300] * 1e12
chanmaskcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar1_ffp8_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_cls.fits')[:300] * 1e12
inputcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_cls_ilcbeam_lmax299.fits') #Input
diffusecls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_cls.fits')[:300] * 1e12 #Diffuse & covar6
ell = np.arange(300)
invtwopi = 1. /(2.*mh.pi)

#Bin the data for visual clarity
binlen = 3 #Maybe 4?
chanmaskcls_binned = np.mean(np.reshape(chanmaskcls_corrected,(-1,binlen)),axis=-1)
inputcls_binned = np.mean(np.reshape(inputcls_corrected,(-1,binlen)),axis=-1)
#newcls_binned = np.mean(np.reshape(newcls_corrected,(-1,binlen)),axis=-1)
diffusecls_binned = np.mean(np.reshape(diffusecls,(-1,binlen)),axis=-1)
#theorycls_binned = np.mean(np.reshape(theorycls_corrected,(-1,binlen)),axis=-1)
ell_binned = np.mean(np.reshape(ell,(-1,binlen)),axis=-1)

f, ax = plt.subplots()
ax.plot(ell_binned,chanmaskcls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'Diffuse \& covar1') #Directional ILC $(N=1)$')
ax.plot(ell_binned,inputcls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'Input')
#ax.plot(ell_binned,newcls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'covar10')
ax.plot(ell_binned,diffusecls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'Diffuse \& covar15')
#ax.plot(ell_binned,theorycls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'Theory')
ax.set_xlim([0,300])
ax.set_xlabel(r'Multipole $\ell$')
ax.set_ylabel(r'$\ell (\ell + 1) C_{\ell} / 2 \pi$ $[{\mu\mathrm{K}}^2]$')
ax.legend(loc='lower right')
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/ffp8spec_chanmask7.pdf')


