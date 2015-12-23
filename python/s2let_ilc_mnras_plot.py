import numpy as np
import healpy as hp
import math as mh
import astropy.io.fits as af
import copy as cp
import matplotlib.pyplot as plt
import distinct_colours as dc

#Latex font settings
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=21.0) #21 for power spectra #17 for maps

#Increase line thicknesses
plt.rc('axes', linewidth=2.0)
plt.rc('xtick.major', width=2.0)
plt.rc('xtick.minor', width=2.0)
plt.rc('ytick.major', width=2.0)
plt.rc('ytick.minor', width=2.0)
plt.rc('lines', linewidth=1.5)

#Plot N=1 and NILC maps
n1map = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/scal_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_5_recon.fits') * 1e6 #uK
'''nilcmap = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/COM_CMB_IQU-smica-field-Int_2048_R2.00.fits') * 1e6 #uK'''
psmask = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits')
n1map_masked = cp.deepcopy(n1map)
n1map_masked[psmask == 0] = 0.
'''nilcmap_masked = cp.deepcopy(nilcmap)
nilcmap_masked[psmask == 0] = 0.'''

#ffp8map = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/scal_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_3_recon.fits') * 1e6 #uK
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
hp.mollview(n1map_masked,unit=r'$\mu\mathrm{K}$',title=r'SILC $(N=5)$',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/n5map4.pdf')'''

'''hp.mollview(nilcmap_masked,unit=r'$\mu\mathrm{K}$',title=r'FFP8 input',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/inputmap2_ffp8.pdf')'''

#Plot difference map
'''n1alms = hp.map2alm(n1map,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(n1map))[:257]
mapbeam = 1. / hp.gauss_beam(mh.radians(5./60.),lmax=256) #Remove existing map beam
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
hp.almxfl(n1alms,pixrecip*smoothbeam*mapbeam,inplace=True)
n1map_downgrade = hp.alm2map(n1alms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_dg.fits',n1map_downgrade)

nilcalms = hp.map2alm(nilcmap,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(nilcmap))[:257]
#mapbeam = 1. / af.getdata('/Users/keir/Documents/s2let_ilc_planck/COM_CMB_IQU-sevem-field-Int_2048_R2.01_full.fits',2).field('INT_BEAM')[:257]
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
#ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=256) #5 arcmin
hp.almxfl(nilcalms,smoothbeam*pixrecip,inplace=True)
nilcmap_downgrade = hp.alm2map(nilcalms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map_dg.fits',nilcmap_downgrade)'''

'''ffp8alms = hp.map2alm(ffp8map,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(ffp8map))[:257]
mapbeam = 1. / hp.gauss_beam(mh.radians(5./60.),lmax=256) #Remove existing map beam
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
hp.almxfl(ffp8alms,pixrecip*mapbeam*smoothbeam,inplace=True)
ffp8map_downgrade = hp.alm2map(ffp8alms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/scal_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_3_recon_lmin2_lmax256_fits',ffp8map_downgrade)'''

'''inputalms = hp.map2alm(inputmap,lmax=256)
pixrecip = 1. / hp.pixwin(hp.get_nside(inputmap))[:257]
smoothbeam = hp.gauss_beam(np.radians(80./60.),lmax=256) #80 arcmin smoothing
smoothbeam[:2] = 0. #Remove monopole and dipole
ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=256) #5 arcmin
hp.almxfl(inputalms,pixrecip*smoothbeam*ilcbeam,inplace=True)
inputmap_downgrade = hp.alm2map(inputalms,nside=128,pixwin=True)
hp.write_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map_ilcbeam_lmax256_lmin2.fits',inputmap_downgrade)'''

'''n1map_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_dg.fits')
nilcmap_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_map_dg.fits')'''

#ffp8map_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr2_lmin2_lmax256.fits')
#inputmap_downgrade = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/scal_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_1_recon_lmin2_lmax256_fits')
'''inputmap_downgrade = hp.remove_dipole(inputmap_downgrade) #Checking that no mono-/dipole'''
#ut78 = hp.ud_grade(hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits'),nside_out=128)
#ut78 = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr2_UT78.fits')
'''ffp8map_downgrade[ut78 < 1] = None
inputmap_downgrade[ut78 < 1] = None'''

'''LIM = 300
hp.mollview(n1map_downgrade,unit=r'$\mu\mathrm{K}$',title=r'Directional ILC $(N=1)$',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/n1map_dg.pdf')

hp.mollview(nilcmap_downgrade,unit=r'$\mu\mathrm{K}$',title=r'NILC',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/nilcmap_dg.pdf')'''

'''diffmap = n1map_downgrade - nilcmap_downgrade
print max(diffmap[ut78 == 1]), min(diffmap[ut78 == 1]), np.mean(diffmap[ut78 == 1]), np.std(diffmap[ut78 == 1])
print np.std(n1map_downgrade[ut78 == 1]), np.std(nilcmap_downgrade[ut78 == 1])
LIM = 15
hp.mollview(diffmap,unit=r'$\mu\mathrm{K}$',title=r'SILC $(N = 1)$ - input [FFP8]',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/diffmap6_n1_input.pdf')'''

'''diffmap = ffp8map_downgrade - inputmap_downgrade
print max(diffmap[ut78 == 1]), min(diffmap[ut78 == 1]), np.mean(diffmap[ut78 == 1]), np.std(diffmap[ut78 == 1])
print np.std(ffp8map_downgrade[ut78 == 1]), np.std(inputmap_downgrade[ut78 == 1])
LIM = 5
hp.mollview(diffmap,unit=r'$\mu\mathrm{K}$',title=r'$(N = 3)$ - $(N = 1)$',min=-1.*LIM,max=LIM)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/diffmap7_n3_n1_5.pdf')'''

#Plot N=1 and NILC spectra
#ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=299) #5 arcmin

n1cls_masked = hp.anafast(n1map_masked,lmax=3399)
pixrecip = 1. / hp.pixwin(hp.get_nside(n1map_masked))[:3400]
fsky = 1. - (float(np.sum(np.logical_not(psmask))) / len(psmask))
print fsky
chanmaskcls_corrected = (n1cls_masked * pixrecip * pixrecip) / fsky
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/scal_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_5_recon_masked.fits',chanmaskcls_corrected)
print "Calculated spectra"

'''nilccls_masked = hp.anafast(nilcmap_masked,lmax=3399)
#pixrecip = 1. / hp.pixwin(hp.get_nside(nilcmap_masked))[:3400]
fsky = 1. - (float(np.sum(np.logical_not(psmask))) / len(psmask))
print fsky
mapbeam = 1. / af.getdata('/Users/keir/Documents/s2let_ilc_planck/COM_CMB_IQU-smica-field-Int_2048_R2.00.fits',2).field('BEAM_WF')[:3400]
ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=3399) #5 arcmin
smicacls = (nilccls_masked * mapbeam * mapbeam * ilcbeam * ilcbeam) / fsky
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/smica_pr2_ilcbeam_lmax3399_cls_masked.fits',smicacls)'''

lmin = 2
lmax = 3392
theorycls = hp.read_cl('/Users/keir/Software/camb/planck2015_4_scalCls.fits')[0][lmin:lmax] * 1e12 #(uK)^2
ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=3391)[lmin:] #5 arcmin
theorycls_corrected = theorycls * ilcbeam * ilcbeam

'''newcls_masked = hp.anafast(newmap_masked,lmax=299)
pixrecip = 1. / hp.pixwin(hp.get_nside(newmap))[:300]
fsky = 1. - (float(np.sum(np.logical_not(psmask))) / len(psmask))
print fsky
newcls_corrected = (newcls_masked * pixrecip * pixrecip) / fsky
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_cls_masked.fits',newcls_corrected)'''

'''inputcls = hp.anafast(inputmap,lmax=3399)
pixrecip = 1. / hp.pixwin(hp.get_nside(inputmap))[:3400]
ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=3399) #5 arcmin
inputcls_corrected = inputcls * pixrecip * pixrecip * ilcbeam * ilcbeam
hp.write_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_cls_ilcbeam_lmax3399.fits',inputcls_corrected)'''

#newcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar10_ffp8_deconv_tapered_thresh_lmax3600_300_hybridC_0_1_recon_cls.fits')[:300] * 1e12

#FFP8
'''chanmaskcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_cls_masked.fits')[lmin:lmax] * 1e12
n2cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_2_recon_cls_masked.fits')[lmin:lmax] * 1e12
n3cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_3_recon_cls_masked.fits')[lmin:lmax] * 1e12
n4cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_masked.fits')[lmin:lmax] * 1e12
n5cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_5_recon_cls_masked.fits')[lmin:lmax] * 1e12
inputcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/ffp8_cmb_scl_000_full_cls_ilcbeam_lmax3399.fits')[lmin:lmax]'''

'''n1cls_debias = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_cls_debias.fits')[lmin:lmax] * 1e12
n2cls_debias = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_2_recon_cls_debias.fits')[lmin:lmax] * 1e12
n3cls_debias = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_3_recon_cls_debias.fits')[lmin:lmax] * 1e12
n4cls_debias = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_debias.fits')[lmin:lmax] * 1e12

n1cls_deresid = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_cls_deresid.fits')[lmin:lmax] * 1e12
n2cls_deresid = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_2_recon_cls_deresid.fits')[lmin:lmax] * 1e12
n3cls_deresid = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_3_recon_cls_deresid.fits')[lmin:lmax] * 1e12
n4cls_deresid = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/ffp8_pla_data/s2let_ilc_covar15_ffp8_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_deresid.fits')[lmin:lmax] * 1e12'''

#PR2
chanmaskcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/scal_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_6_1_recon_masked.fits')[lmin:lmax] #* 1e12
'''n2cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/diffuse_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_2_recon_cls_masked.fits')[lmin:lmax] * 1e12
n3cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/diffuse_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_3_recon_cls_masked.fits')[lmin:lmax] * 1e12
n4cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/diffuse_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_4_recon_cls_masked.fits')[lmin:lmax] * 1e12
n5cls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/diffuse_data/s2let_ilc_covar15_planck_diffuse_deconv_tapered_thresh_lmax3600_3600_hybridC_0_5_recon_cls_masked.fits')[lmin:lmax] * 1e12'''
inputcls_corrected = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/nilc_pr2_cls_masked.fits')[lmin:lmax] #* 1e12
smicacls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/smica_pr2_ilcbeam_lmax3399_cls_masked.fits')[lmin:lmax] #* 1e12'''
'''commcls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/commander_pr2_ilcbeam_lmax3400_cls.fits')[lmin:lmax] * 1e12
sevemcls = hp.read_cl('/Users/keir/Documents/s2let_ilc_planck/sevem_pr2_ilcbeam_lmax3100_cls.fits')[lmin:3092] * 1e12'''

#ilcbeam = hp.gauss_beam(np.radians(5./60.),lmax=1299) #5 arcmin
#diffusecls = hp.read_cl('/Users/keir/Software/camb/planck2015_2_scalCls.fits')[0][:1300] * 1e12 * ilcbeam * ilcbeam #Diffuse & covar6
#ell_long = np.arange(lmin,3392)
ell = np.arange(lmin,lmax)
#ell_short = np.arange(lmin,1292)
invtwopi = 1. /(2.*mh.pi)

chanmaskdls = chanmaskcls_corrected*ell*(ell+1)*invtwopi
'''n2dls = n2cls_corrected*ell*(ell+1)*invtwopi
n3dls = n3cls_corrected*ell*(ell+1)*invtwopi
n4dls = n4cls_corrected*ell*(ell+1)*invtwopi
n5dls = n5cls_corrected*ell*(ell+1)*invtwopi'''
inputdls = inputcls_corrected*ell*(ell+1)*invtwopi
smicadls = smicacls*ell*(ell+1)*invtwopi
'''commdls = commcls*ell*(ell+1)*invtwopi
sevemdls = sevemcls*ell_short*(ell_short+1)*invtwopi'''
theorydls = theorycls_corrected*ell*(ell+1)*invtwopi
'''n1debiasdls = n1cls_debias*ell*(ell+1)*invtwopi
n2debiasdls = n2cls_debias*ell*(ell+1)*invtwopi
n3debiasdls = n3cls_debias*ell*(ell+1)*invtwopi
n4debiasdls = n4cls_debias*ell*(ell+1)*invtwopi
n1deresiddls = n1cls_deresid*ell*(ell+1)*invtwopi
n2deresiddls = n2cls_deresid*ell*(ell+1)*invtwopi
n3deresiddls = n3cls_deresid*ell*(ell+1)*invtwopi
n4deresiddls = n4cls_deresid*ell*(ell+1)*invtwopi'''

#Bin the data for visual clarity
binlen = 10 #Maybe 4?
chanmaskdls_binned = np.mean(np.reshape(chanmaskdls,(-1,binlen)),axis=-1)
'''n2dls_binned = np.mean(np.reshape(n2dls,(-1,binlen)),axis=-1)
n3dls_binned = np.mean(np.reshape(n3dls,(-1,binlen)),axis=-1)
n4dls_binned = np.mean(np.reshape(n4dls,(-1,binlen)),axis=-1)
n5dls_binned = np.mean(np.reshape(n5dls,(-1,binlen)),axis=-1)'''
inputdls_binned = np.mean(np.reshape(inputdls,(-1,binlen)),axis=-1)
'''n1debiasdls_binned = np.mean(np.reshape(n1debiasdls,(-1,binlen)),axis=-1)
n2debiasdls_binned = np.mean(np.reshape(n2debiasdls,(-1,binlen)),axis=-1)
n3debiasdls_binned = np.mean(np.reshape(n3debiasdls,(-1,binlen)),axis=-1)
n4debiasdls_binned = np.mean(np.reshape(n4debiasdls,(-1,binlen)),axis=-1)
n1deresiddls_binned = np.mean(np.reshape(n1deresiddls,(-1,binlen)),axis=-1)
n2deresiddls_binned = np.mean(np.reshape(n2deresiddls,(-1,binlen)),axis=-1)
n3deresiddls_binned = np.mean(np.reshape(n3deresiddls,(-1,binlen)),axis=-1)
n4deresiddls_binned = np.mean(np.reshape(n4deresiddls,(-1,binlen)),axis=-1)'''
smicadls_binned = np.mean(np.reshape(smicadls,(-1,binlen)),axis=-1)
'''commdls_binned = np.mean(np.reshape(commdls,(-1,binlen)),axis=-1)
sevemdls_binned = np.mean(np.reshape(sevemdls,(-1,binlen)),axis=-1)'''
#newcls_binned = np.mean(np.reshape(newcls_corrected,(-1,binlen)),axis=-1)
#diffusecls_binned = np.mean(np.reshape(diffusecls,(-1,binlen)),axis=-1)
theorydls_binned = np.mean(np.reshape(theorydls,(-1,binlen)),axis=-1)
#ell_binned_long = np.mean(np.reshape(ell_long,(-1,binlen)),axis=-1)
ell_binned = np.mean(np.reshape(ell,(-1,binlen)),axis=-1)
#ell_binned_short = np.mean(np.reshape(ell_short,(-1,binlen)),axis=-1)

#Cosmic variance

cols = dc.get_distinct(4)

f, (ax0,ax1) = plt.subplots(2,sharex=True,figsize=(8,8))
ax0.plot(ell_binned,chanmaskdls_binned,color=cols[0],label = r'SILC $(N = 1)$')
'''ax0.plot(ell_binned,n2dls_binned,color=cols[1],label = r'$N=2$')
ax0.plot(ell_binned,n3dls_binned,color=cols[2],label = r'$N=3$')
ax0.plot(ell_binned,n4dls_binned,color=cols[3],label = r'$N=4$')
ax0.plot(ell_binned,n5dls_binned,color=cols[4],label = r'$N=5$')'''
#ax0.plot(ell_binned,inputdls_binned,ls='--',color='k',label = r'Input')
'''ax0.plot(ell_binned,n1debiasdls_binned,color=cols[0],label = r'$N=1$ (debiased)')
ax0.plot(ell_binned,n2debiasdls_binned,color=cols[1],label = r'$N=2$ (debiased)')
ax0.plot(ell_binned,n3debiasdls_binned,color=cols[2],label = r'$N=3$ (debiased)')
ax0.plot(ell_binned,n4debiasdls_binned,color=cols[3],label = r'$N=4$ (debiased)')'''
ax0.plot(ell_binned,inputdls_binned,color=cols[1],label = r'NILC')
ax0.plot(ell_binned,smicadls_binned,color=cols[3],label = r'SMICA')
'''ax0.plot(ell_binned,commdls_binned,color=cols[3],label = r'Commander')
ax0.plot(ell_binned_short,sevemdls_binned,color=cols[4],label = r'SEVEM')'''
#ax.plot(ell_binned,newcls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'covar10')
#ax.plot(ell_binned,diffusecls_binned*ell_binned*(ell_binned+1)*invtwopi,label = r'Theory')
ax0.plot(ell_binned,theorydls_binned,ls='--',color='k',label = r'Theory')
ax0.set_xlim([0,3400])
ax0.set_ylabel(r'$D_{\ell}$ $[{\mu\mathrm{K}}^2]$', labelpad = 3)
ax0.legend(prop={'size':21},frameon=False) #loc='lower right') #Usually 18 #Sometimes 20

#Cosmic variance

#Plot residuals
#FFP8
'''ax1.plot(ell_binned,chanmaskdls_binned - inputdls_binned,color=cols[0],label = r'$N=1$')
ax1.plot(ell_binned,n2dls_binned - inputdls_binned,color=cols[1],label = r'$N=2$')
ax1.plot(ell_binned,n3dls_binned - inputdls_binned,color=cols[2],label = r'$N=3$')
ax1.plot(ell_binned,n4dls_binned - inputdls_binned,color=cols[3],label = r'$N=4$')
ax1.plot(ell_binned,n5dls_binned - inputdls_binned,color=cols[4],label = r'$N=5$')'''
'''ax1.plot(ell_binned,n1debiasdls_binned - inputdls_binned,color=cols[0],label = r'$N=1$ (debiased)')
ax1.plot(ell_binned,n2debiasdls_binned - inputdls_binned,color=cols[1],label = r'$N=2$ (debiased)')
ax1.plot(ell_binned,n3debiasdls_binned - inputdls_binned,color=cols[2],label = r'$N=3$ (debiased)')
ax1.plot(ell_binned,n4debiasdls_binned - inputdls_binned,color=cols[3],label = r'$N=4$ (debiased)')
ax1.plot(ell_binned,n1deresiddls_binned - inputdls_binned,ls='-.',color=cols[0],label = r'$N=1$ (debiased)')
ax1.plot(ell_binned,n2deresiddls_binned - inputdls_binned,ls='-.',color=cols[1],label = r'$N=2$ (debiased)')
ax1.plot(ell_binned,n3deresiddls_binned - inputdls_binned,ls='-.',color=cols[2],label = r'$N=3$ (debiased)')
ax1.plot(ell_binned,n4deresiddls_binned - inputdls_binned,ls='-.',color=cols[3],label = r'$N=4$ (debiased)')'''

#PR2
ax1.plot(ell_binned,chanmaskdls_binned - theorydls_binned,color=cols[0],label = r'$N=1$')
'''ax1.plot(ell_binned,n2dls_binned - theorydls_binned,color=cols[1],label = r'$N=2$')
ax1.plot(ell_binned,n3dls_binned - theorydls_binned,color=cols[2],label = r'$N=3$')
ax1.plot(ell_binned,n4dls_binned - theorydls_binned,color=cols[3],label = r'$N=4$')
ax1.plot(ell_binned,n5dls_binned - theorydls_binned,color=cols[4],label = r'$N=5$')'''
ax1.plot(ell_binned,inputdls_binned - theorydls_binned,color=cols[1],label = r'NILC')
ax1.plot(ell_binned,smicadls_binned - theorydls_binned,color=cols[3],label = r'SMICA')
'''ax1.plot(ell_binned,commdls_binned - theorydls_binned,color=cols[3],label = r'Commander')
ax1.plot(ell_binned_short,sevemdls_binned - theorydls_binned[:len(sevemdls_binned)],color=cols[4],label = r'SEVEM')'''

#ax1.plot(ell_binned,inputdls_binned - theorydls_binned,color='g',label = r'NILC')
ax1.axhline(y=0.,ls='--',color='k')
ax1.set_xlim([0,3400])
#ax1.set_ylim([-200,200])
ax1.set_xlabel(r'Multipole $\ell$', labelpad = 10)
#xlab.set_position((0.5, 0.1))
ax1.set_ylabel(r'$\Delta D_{\ell}$ $[{\mu\mathrm{K}}^2]$', labelpad = -4)

yticks = ax1.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-1].label1.set_visible(False)

f.subplots_adjust(hspace=0,right=0.99)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

#ax.set_yscale('log')
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/smicaspec_n1_resids5.pdf')

#Panel 3
f2, ax2 = plt.subplots(1,figsize=(8,4))
botlim = 150

#FFP8
'''ax2.plot(ell_binned[:botlim],chanmaskdls_binned[:botlim] - inputdls_binned[:botlim],color=cols[0],label = r'$N=1$')
ax2.plot(ell_binned[:botlim],n2dls_binned[:botlim] - inputdls_binned[:botlim],color=cols[1],label = r'$N=2$')
ax2.plot(ell_binned[:botlim],n3dls_binned[:botlim] - inputdls_binned[:botlim],color=cols[2],label = r'$N=3$')
ax2.plot(ell_binned[:botlim],n4dls_binned[:botlim] - inputdls_binned[:botlim],color=cols[3],label = r'$N=4$')
ax2.plot(ell_binned[:botlim],n5dls_binned[:botlim] - inputdls_binned[:botlim],color=cols[4],label = r'$N=5$')'''
'''ax2.plot(ell_binned[:botlim],n1debiasdls_binned[:botlim] - inputdls_binned[:botlim],color=cols[0],label = r'$N=1$ (debiased)')
ax2.plot(ell_binned[:botlim],n2debiasdls_binned[:botlim] - inputdls_binned[:botlim],color=cols[1],label = r'$N=2$ (debiased)')
ax2.plot(ell_binned[:botlim],n3debiasdls_binned[:botlim] - inputdls_binned[:botlim],color=cols[2],label = r'$N=3$ (debiased)')
ax2.plot(ell_binned[:botlim],n4debiasdls_binned[:botlim] - inputdls_binned[:botlim],color=cols[3],label = r'$N=4$ (debiased)')
ax2.plot(ell_binned[:botlim],n1deresiddls_binned[:botlim] - inputdls_binned[:botlim],ls='-.',color=cols[0],label = r'$N=1$ (debiased)')
ax2.plot(ell_binned[:botlim],n2deresiddls_binned[:botlim] - inputdls_binned[:botlim],ls='-.',color=cols[1],label = r'$N=2$ (debiased)')
ax2.plot(ell_binned[:botlim],n3deresiddls_binned[:botlim] - inputdls_binned[:botlim],ls='-.',color=cols[2],label = r'$N=3$ (debiased)')
ax2.plot(ell_binned[:botlim],n4deresiddls_binned[:botlim] - inputdls_binned[:botlim],ls='-.',color=cols[3],label = r'$N=4$ (debiased)')'''

#PR2
ax2.plot(ell_binned[:botlim],chanmaskdls_binned[:botlim] - theorydls_binned[:botlim],color=cols[0],label = r'$N=1$')
'''ax2.plot(ell_binned[:botlim],n2dls_binned[:botlim] - theorydls_binned[:botlim],color=cols[1],label = r'$N=2$')
ax2.plot(ell_binned[:botlim],n3dls_binned[:botlim] - theorydls_binned[:botlim],color=cols[2],label = r'$N=3$')
ax2.plot(ell_binned[:botlim],n4dls_binned[:botlim] - theorydls_binned[:botlim],color=cols[3],label = r'$N=4$')
ax2.plot(ell_binned[:botlim],n5dls_binned[:botlim] - theorydls_binned[:botlim],color=cols[4],label = r'$N=5$')'''
ax2.plot(ell_binned[:botlim],inputdls_binned[:botlim] - theorydls_binned[:botlim],color=cols[1],label = r'NILC')
ax2.plot(ell_binned[:botlim],smicadls_binned[:botlim] - theorydls_binned[:botlim],color=cols[3],label = r'SMICA')
'''ax2.plot(ell_binned[:botlim],commdls_binned[:botlim] - theorydls_binned[:botlim],color=cols[3],label = r'Commander')
ax2.plot(ell_binned_short[:botlim],sevemdls_binned[:botlim] - theorydls_binned[:botlim],color=cols[4],label = r'SEVEM')'''

ax2.axhline(y=0.,ls='--',color='k')
ax2.set_xlim([0,1500])
#ax2.set_xlabel(r'Multipole $\ell$')
ax2.set_ylabel(r'$\Delta D_{\ell}$ $[{\mu\mathrm{K}}^2]$', labelpad = -4)
ax2.xaxis.tick_top()
yticks2 = ax2.yaxis.get_major_ticks()
yticks2[0].label1.set_visible(False)
yticks2[-1].label1.set_visible(False)
f2.subplots_adjust(right=0.99)
plt.savefig('/Users/keir/Documents/s2let_ilc_latex/s2let_ilc_papers/s2let_ilc_temp/smicaspec_n1_resids5_botpan.pdf')


