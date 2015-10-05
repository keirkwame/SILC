import numpy as np
import healpy as hp
import math as mh
import multiprocessing as mg
import matplotlib.pyplot as plt

def planck(i):
    print "Doing another map analysis"
    planckfits = planckdir + planckroot + planckcode[i] + planckend
    planckmaps = hp.read_map(planckfits)
    #pixrecip = 1. / hp.pixwin(hp.get_nside(planckmaps[i]))
    planckcls = hp.anafast(planckmaps,lmax=3998) #* pixrecip * pixrecip
    planckclsfits = planckdir + planckroot + planckcode[i] + planckclsend
    hp.write_cl(planckclsfits,planckcls)
    #ax.plot(planckcls[i]*ell*(ell+1)*(10**12)*invtwopi,label=r"Planck 2015 data: %s GHz" % planckcode[i])

def ffp6(i):
    ffp6fits = ffp6dir + ffp6dircode[i] + planckprefix[i] + ffp6root + ffp6code[i] + ffp6end
    ffp6maps = hp.read_map(ffp6fits)
    #pixrecip = 1. / hp.pixwin(hp.get_nside(planckmaps[i]))
    ffp6cls = hp.anafast(ffp6maps,lmax=3998) #* pixrecip * pixrecip
    ffp6clsfits = ffp6dir + ffp6dircode[i] + planckprefix[i] + ffp6root + ffp6code[i] + ffp6clsend
    hp.write_cl(ffp6clsfits,ffp6cls)
    #ax.plot(ffp6cls[i]*ell*(ell+1)*(10**12)*invtwopi,label=r"FFP6 simulations: %s GHz" % ffp6code[i])

planckdir = '/Users/keir/Documents/s2let_ilc_planck/maps/PR1/' #'/home/keir/s2let_ilc_data/'
planckprefix = ['LFI','LFI','LFI','HFI','HFI','HFI','HFI','HFI','HFI']
planckroot = 'planck_filled_' #'planck_filled_minusgaussps2_' #'_SkyMap_'
planckcode = ['30','44','70','100','143','217','353','545','857'] #['030_1024_R2.01','044_1024_R2.01','070_2048_R2.01','100_2048_R2.00','143_2048_R2.00','217_2048_R2.00','353_2048_R2.00','545_2048_R2.00','857_2048_R2.00']
planckend = '_pr1.fits' #'_full.fits'
planckclsend = '_pr1_cls.fits' #'_full_cls.fits'

ffp6dir = '/Users/keir/Documents/ffp6_pla_fiducial_maps/PR1/frequency_map/' #'/Users/keir/Documents/s2let_ilc_planck/ffp6_data/' #'/home/keir/s2let_ilc_data/ffp6_data_withPS/'
ffp6dircode = ['030/','044/','070/','100/','143/','217/','353/','545/','857/']
ffp6root = '_SimMap_' #'ffp6_fiducial_withPS_'
ffp6code = ['030-ffp6_1024','044-ffp6_1024','070-ffp6_1024','100-ffp6_2048','143-ffp6_2048','217-ffp6_2048','353-ffp6_2048','545-ffp6_2048','857-ffp6_2048']
ffp6end = '_R1.10_nominal.fits' #'_pr1.fits'
ffp6clsend = '_R1.10_nominal_cls.fits' #'_pr1_cls.fits'

ell = np.arange(3999)
invtwopi = 1. / (2.*mh.pi)

nprocess = 4

'''pool = mg.Pool(nprocess)
planck_output = pool.map(planck,np.arange(9))
pool.close()
pool.join()'''

'''pool2 = mg.Pool(nprocess)
ffp6_output = pool2.map(ffp6,np.arange(9))
pool2.close()
pool2.join()'''

fig,ax = plt.subplots(3,3)
ax = ax.ravel()
stepno = 10
planckcls = [None]*9
ffp6cls = [None]*9
oldffp6cls = [None]*9
for i in xrange(len(planckcls)):
    planckclsfits = planckdir + planckroot + planckcode[i] + planckclsend
    planckcls[i] = hp.read_cl(planckclsfits)
    ffp6clsfits = ffp6dir + ffp6dircode[i] + planckprefix[i] + ffp6root + ffp6code[i] + ffp6clsend
    ffp6cls[i] = hp.read_cl(ffp6clsfits)
    oldffp6clsfits = '/Users/keir/Documents/s2let_ilc_planck/ffp6_data/' + 'ffp6_fiducial_withPS_' + planckcode[i] + '_cls.fits'
    oldffp6cls[i] = hp.read_cl(oldffp6clsfits)
    ax[i].set_xlabel(r"Multipole $\ell$")
    ax[i].set_ylabel(r"Power spectrum $\frac{\ell (\ell + 1) C_{\ell}}{2 \pi}$ (uK^2)")
    if i < 9:
        #ax[i].set_ylim([0,1e5])
        ax[i].set_yscale('log')
        ax[i].plot(ell[::stepno],planckcls[i][::stepno]*ell[::stepno]*(ell[::stepno]+1)*(10**12)*invtwopi,label=r"Planck 2013 data: %s GHz" % planckcode[i])
        ax[i].plot(ell[::stepno],ffp6cls[i][::stepno]*ell[::stepno]*(ell[::stepno]+1)*(10**12)*invtwopi,label=r"FFP6 simulations (PLA): %s GHz" % planckcode[i])
        ax[i].plot(ell[::stepno],oldffp6cls[i][::stepno]*ell[::stepno]*(ell[::stepno]+1)*(10**12)*invtwopi,label=r"FFP6 simulations (NERSC): %s GHz" % planckcode[i])
    else:
        ax[i].plot(planckcls[i]*ell*(ell+1)*(10**12)*invtwopi,label=r"Planck 2013 data [original beams]: %s GHz" % planckcode[i])
        ax[i].plot(ffp6cls[i]*ell*(ell+1)*(10**12)*invtwopi,label=r"FFP6 simulations [original beams]: %s GHz" % planckcode[i])
    ax[i].legend(loc='center right') #loc='upper left')

#ax.set_ylim([0,10000])
plt.show()