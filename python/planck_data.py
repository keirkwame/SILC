import numpy as np
import healpy as hp
import math as mh
import multiprocessing as mg
import matplotlib.pyplot as plt

def planck(i):
    planckfits = planckdir + planckroot + planckcode[i] + planckend
    planckmaps[i] = hp.read_map(planckfits)
    #pixrecip = 1. / hp.pixwin(hp.get_nside(planckmaps[i]))
    planckcls[i] = hp.anafast(planckmaps[i],lmax=3998) #* pixrecip * pixrecip
    planckclsfits = planckdir + planckroot + planckcode[i] + planckclsend
    hp.write_cl(planckclsfits,planckcls[i])
    #ax.plot(planckcls[i]*ell*(ell+1)*(10**12)*invtwopi,label=r"Planck 2015 data: %s GHz" % planckcode[i])

def ffp6(i):
    ffp6fits = ffp6dir + ffp6root + ffp6code[i] + ffp6end
    ffp6maps[i] = hp.read_map(ffp6fits)
    #pixrecip = 1. / hp.pixwin(hp.get_nside(planckmaps[i]))
    ffp6cls[i] = hp.anafast(ffp6maps[i],lmax=3998) #* pixrecip * pixrecip
    ffp6clsfits = ffp6dir + ffp6root + ffp6code[i] + ffp6clsend
    hp.write_cl(ffp6clsfits,ffp6cls[i])
    #ax.plot(ffp6cls[i]*ell*(ell+1)*(10**12)*invtwopi,label=r"FFP6 simulations: %s GHz" % ffp6code[i])

planckdir = '/home/keir/s2let_ilc_data/' #'/Users/keir/Documents/s2let_ilc_planck/deconv_data/'
planckroot = 'planck_deconv_tapered_'
planckcode = ['30','44','70','100','143','217','353','545','857']
planckend = '_pr2.fits'
planckclsend = '_pr2_cls.fits'

ffp6dir = '/home/keir/s2let_ilc_data/ffp6_data_withPS/' #'/Users/keir/Documents/s2let_ilc_planck/ffp6_data/'
ffp6root = 'ffp6_fiducial_withPS_tapered_'
ffp6code = planckcode
ffp6end = '.fits'
ffp6clsend = '_cls.fits'

ell = np.arange(3999)
invtwopi = 1. / (2.*mh.pi)

nprocess = 9

'''pool = mg.Pool(nprocess)
planck_output = pool.map(planck,np.arange(9))
pool.close()
pool.join()'''

pool2 = mg.Pool(nprocess)
ffp6_output = pool.map(ffp6,np.arange(9))
pool2.close()
pool2.join()

'''fig,ax = plt.subplots(1,1)
ax.legend()
plt.show()'''