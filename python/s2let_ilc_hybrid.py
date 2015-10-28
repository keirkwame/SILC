import numpy as np
import healpy as hp
import math as mh
import multiprocessing as mg
import multiprocessing.pool
import pys2let as ps
import random
import string
import itertools
import os

#From http://stackoverflow.com/questions/15118344/system-error-while-running-subprocesses-using-multiprocessing
### A helper for letting the forked processes use data without pickling.
_data_name_cands = (
    '_data_' + ''.join(random.sample(string.ascii_lowercase, 10))
    for _ in itertools.count())
class ForkedData(object):
    '''
    Class used to pass data to child processes in multiprocessing without
    really pickling/unpickling it. Only works on POSIX.

    Intended use:
        - The master process makes the data somehow, and does e.g.
          data = ForkedData(the_value)
        - The master makes sure to keep a reference to the ForkedData object
          until the children are all done with it, since the global reference
          is deleted to avoid memory leaks when the ForkedData object dies.
        - Master process constructs a multiprocessing.Pool *after*
          the ForkedData construction, so that the forked processes
          inherit the new global.
        - Master calls e.g. pool.map with data as an argument.
        - Child gets the real value through data.value, and uses it read-only.
    '''
    # TODO: does data really need to be used read-only? don't think so...
    # TODO: more flexible garbage collection options
    def __init__(self, val):
        g = globals()
        self.name = next(n for n in _data_name_cands if n not in g)
        g[self.name] = val
        self.master_pid = os.getpid()

    def __getstate__(self):
        if os.name != 'posix':
            raise RuntimeError("ForkedData only works on OSes with fork()")
        return self.__dict__

    @property
    def value(self):
        return globals()[self.name]

    def __del__(self):
        if os.getpid() == self.master_pid:
            del globals()[self.name]

#From http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
class NoDaemonProcess(mg.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(mg.pool.Pool):
    Process = NoDaemonProcess

def smoothworker(i): #(j,n,map_index1,map_index2,smoothing_lmax,scale_fwhm) [map_index2<=map_index1]
    print "Smoothing independent covariance element", i[2], ",", i[3]

    #Map loading within sub-process
    if i[0] >= 0: #Wavelet scales
        wav_fits1 = wav_fits_root[i[2]] + '_' + wavparam_code + str(i[0]) + '_n' + str(i[1]+1) + '_double.npy'
        wav_fits2 = wav_fits_root[i[3]] + '_' + wavparam_code + str(i[0]) + '_n' + str(i[1]+1) + '_double.npy'
    else: #Scaling function
        wav_fits1 = scal_fits[i[2]][:-4] + '_double.npy'
        wav_fits2 = scal_fits[i[3]][:-4] + '_double.npy'
    map1 = np.real(np.load(wav_fits1,mmap_mode='r')) #Throw away zero imaginary part
    map2 = np.real(np.load(wav_fits2,mmap_mode='r'))
    R = np.multiply(map1,map2) + 0.j #Add back in zero imaginary part
    del map1,map2

    alms = ps.lm2lm_hp(ps.map2alm_mw(R,i[4],spin),i[4]) #No pixwin correct. with MW - calc alms to smooth - come out in MW order - so converted to HPX order
    del R
    gausssmooth = hp.gauss_beam(i[5],lmax=i[4]-1)
    hp.almxfl(alms,gausssmooth,inplace=True) #Multiply by gaussian beam

    print "Synthesising smoothed covariance map for element", i[2], ",", i[3]
    Rsmooth = np.real(ps.alm2map_mw(ps.lm_hp2lm(alms,i[4]),i[4],spin)) #Throw away zero imaginary part - input alm's in MW order
    del alms
    
    #SAVE smoothed covariance
    if i[0] >= 0: #Wavelet scales
        R_fits = wav_outfits_root + '_' + wavparam_code + str(i[0]) + '_n' + str(i[1]+1) + '_Rsmooth' + str(i[2]) + str(i[3]) +'.npy'
    else: #Scaling function
        R_fits = scal_outfits[:-4] + '_Rsmooth' + str(i[2]) + str(i[3]) +'.npy'
    np.save(R_fits,Rsmooth)
    del Rsmooth
    
    return 0

def zeropad(i): #(alms[HPX],scale_lmax,smoothing_lmax) - COULD RE-WRITE FOR MW ALM'S???
    print "Zero-padding the alm's"

    #Pre-allocate array for enlarged alm's
    new_alms = np.zeros(((i[2]*(i[2]+1))/2.),dtype=complex)
    for em in xrange(i[1]):
        startindex = em*i[1] - .5*em*(em-1)
        new_startindex = em*i[2] - .5*em*(em-1)
        new_alms[new_startindex:(new_startindex+i[1]-em)] = i[0][startindex:(startindex+i[1]-em)]

    return new_alms

def doubleworker(i): #i = (j,n,map_index,scale_lmax,smoothing_lmax)
    print "Doubling l_max of input map", i[2]+1, "/", nmaps
    
    #Map loading within sub-process
    if i[0] >= 0: #Wavelet scales
        wav_fits = wav_fits_root[i[2]] + '_' + wavparam_code + str(i[0]) + '_n' + str(i[1]+1) + '.npy'
    else: #Scaling function
        wav_fits = scal_fits[i[2]]
    map = np.load(wav_fits,mmap_mode='r') #Map still only stored on disk

    alms = ps.map2alm_mw(map,i[3],spin) #alm's to l(j) - come out in MW order
    del map
    alms = zeropad((ps.lm2lm_hp(alms,i[3]),i[3],i[4])) #New alm's larger
    map = ps.alm2map_mw(ps.lm_hp2lm(alms,i[4]),i[4],spin) #Input alm's in MW order
    del alms
    
    #SAVE doubled map
    if i[0] >= 0: #Wavelet scales
        double_fits = wav_fits_root[i[2]] + '_' + wavparam_code + str(i[0]) + '_n' + str(i[1]+1) + '_double.npy'
    else: #Scaling function
        double_fits = scal_fits[i[2]][:-4] + '_double.npy'
    np.save(double_fits,map)
    del map
    
    return 0

def s2let_ilc(mapsextra): #mapsextra = (j,n)
    print "Running S2LET ILC on wavelet scale", mapsextra[0], "/", jmax, "direction", mapsextra[1]+1, "/", ndir, "\n"

    if mapsextra[0] >= 0: #Wavelet scales
        scale_lmax = wav_bandlims[mapsextra[0]]
    else: #Scaling function
        scale_lmax = int(scal_bandlims)
    smoothing_lmax = 2*(scale_lmax-1)+1

    #Doubling lmax for input maps with zero-padding
    #Serial version
    '''mapsdouble = np.zeros((nrows,ps.mw_size(smoothing_lmax)),dtype=np.complex128) #Pre-allocate array
    for i in xrange(nrows):
        mapsdouble[i,:] = doubleworker((mapsextra[0][i],mapsextra[1],smoothing_lmax,mapsextra[2]))'''
    #Parallel version
    '''mapsextra2 = [(mapsextra[0],mapsextra[1],i,scale_lmax,smoothing_lmax) for i in xrange(nmaps)]
    print "Forming pool"
    pool2 = mg.Pool(nprocess2)
    print "Farming out workers to run doubling function"
    double_output = pool2.map(doubleworker,mapsextra2)
    print "Have returned from doubling workers\n"
    pool2.close()
    pool2.join()
    del pool2'''

    #Calculate scale_fwhm for smoothing kernel
    nsamp = 1200.
    npix = hp.nside2npix(1<<(int(0.5*scale_lmax)-1).bit_length()) #Equivalent no. HEALPIX pixels
    scale_fwhm = 4. * mh.sqrt(nsamp / npix)
    
    #Smooth covariance matrices
    #Serial version
    '''Rsmoothflat = np.zeros_like(Rflat) #Pre-allocate array
    for i in xrange(nindepelems):
        Rsmoothflat[i,:] = smoothworker((Rflat[i],smoothing_lmax,mapsextra[2],gausssmooth,mapsextra[1],mapsextra[3],i,mapsextra[4]))
    del Rflat'''
    #Parallel version
    nindepelems = int(nmaps*(nmaps+1)*.5) #No. indep. elements in symmetric covariance matrix
    Rextra = [None]*nindepelems
    k=0
    for i in xrange(nmaps):
        for j in xrange(i+1):
            Rextra[k] = (mapsextra[0],mapsextra[1],i,j,smoothing_lmax,scale_fwhm)
            k+=1
    print "Forming pool"
    pool3 = mg.Pool(nprocess3)
    print "Farming out workers to run smoothing function"
    R_output = pool3.map(smoothworker,Rextra)
    print "Have returned from smoothing workers\n"
    pool3.close()
    pool3.join()
    del pool3

    #Load R maps and form matrices
    print "Pre-allocating memory for complete covariance tensor\n"
    Rsmooth = np.zeros((ps.mw_size(smoothing_lmax),nmaps,nmaps),dtype=np.float64) #Pre-allocate array
    for i in xrange(nmaps):
        for j in xrange(i+1):
            if mapsextra[0] >= 0: #Wavelet scales
                R_fits = wav_outfits_root + '_' + wavparam_code + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '_Rsmooth' + str(i + 9 - nmaps) + str(j + 9 - nmaps) +'.npy'
            else: #Scaling function
                R_fits = scal_outfits[:-4] + '_Rsmooth' + str(i + 9 - nmaps) + str(j + 9 - nmaps) +'.npy'
            Rsmooth[:,i,j] = np.load(R_fits)
            if i != j:
                Rsmooth[:,j,i] = Rsmooth[:,i,j]

    #Compute inverse covariance matrices
    print "Calculating inverse covariance matrices\n"
    Rinv = np.linalg.inv(Rsmooth) #Parallel vers. slower!?- LARGEST MEMORY COST: 2*9*9*(8000^2)*complex128=0.2TB
    del Rsmooth

    #Compute weights vectors (at each pixel)
    wknumer = np.sum(Rinv,axis=-1)
    del Rinv
    wkdenom = np.sum(wknumer,axis=-1)
    wk = wknumer / wkdenom[:,None]
    del wknumer,wkdenom

    #Map loading within sub-process
    mapsdouble = np.zeros((len(wk),len(wk[0])),dtype=np.float64) #Pre-allocate array
    for i in xrange(nmaps):
        if mapsextra[0] >= 0: #Wavelet scales
            wav_fits = wav_fits_root[i + 9 - nmaps] + '_' + wavparam_code + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '_double.npy'
        else: #Scaling function
            wav_fits = scal_fits[i + 9 - nmaps][:-4] + '_double.npy'
        mapsdouble[:,i] = np.real(np.load(wav_fits,mmap_mode='r')) #Throw away zero imaginary part

    #Dot weights with maps (at each small pixel) - at double l(j)
    finalmap = np.sum(np.multiply(wk,mapsdouble),axis=-1) + 0.j #Add back in zero imaginary part
    del wk,mapsdouble
    
    #Downgrade resolution of MW maps
    print "Downgrading resolution of CMB wavelet map"
    finalmapalms = ps.lm2lm_hp(ps.map2alm_mw(finalmap,smoothing_lmax,spin),smoothing_lmax) #Come out in MW order - so converted to HPX order
    del finalmap
    if mapsextra[0] >= 0: #Wavelet scales
        alms_fname = wav_outfits_root + '_' + wavparam_code + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '_alms.fits'
    else: #Scaling function
        alms_fname = scal_outfits[:-4] + '_alms.fits'
    hp.write_alm(alms_fname,finalmapalms,lmax=scale_lmax-1,mmax=scale_lmax-1)
    del finalmapalms
    finalmapalmstruncate = hp.read_alm(alms_fname)
    finalmaphalf = ps.alm2map_mw(ps.lm_hp2lm(finalmapalmstruncate,scale_lmax),scale_lmax,spin)
    del finalmapalmstruncate
    
    #Saving output map
    if mapsextra[0] >= 0: #Wavelet scales
        wav_outfits = wav_outfits_root + '_' + wavparam_code + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.npy'
    else: #Scaling function
        wav_outfits = scal_outfits[:-4] + '.npy'
    np.save(wav_outfits,finalmaphalf)
    del finalmaphalf

    return 0

def test_ilc(mapsextra): #mapsextra = (j,n) [Testing on 100 GHz map]
    #Loading input map
    wav_fits = wav_fits_root[3] + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.npy'
    if mapsextra[0] == -1: #Scaling function
        wav_fits = scal_fits[3]
    anal_map = np.load(wav_fits,mmap_mode='r') #Map still only stored on disk
    
    #Saving output map
    wav_outfits = wav_outfits_root + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.npy'
    if mapsextra[0] == -1:
        wav_outfits = scal_outfits
    np.save(wav_outfits,anal_map)
    del anal_map

    return 0

if __name__ == "__main__":
    ##Input
    #Set directory structure
    comp = 1
    
    if comp == 0: #Keir's iMac
        nprocess = 1 #For distributing directions within scale
        nprocess2 = 3 #For doubling maps
        nprocess3 = 4 #For smoothing covariance elements
        fitsdir = '/Users/keir/Documents/s2let_ilc_planck/hybrid_data/'
    elif comp == 1: #Hypatia
        nprocess = 1 #For distributing directions within scale
        nprocess2 = 9 #For doubling maps
        nprocess3 = 12 #For smoothing covariance elements
        fitsdir = '/home/keir/s2let_ilc_data/hybrid_data/'
    
    nmaps = 9 #No. maps (Planck = 9)
    ellmax = 3600
    jmin = 0
    lamdas = np.array([60,2,1.3,1.2])
    wavparam_code = 'C'
    l_transitions = np.array([61,513,2017])
    ndir = 2 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals

    fitsroot = 'planck_deconv_tapered_thresh_lmax3600_'
    fitscode = ['30','44','70','100','143','217','353','545','857']
    scal_fits = [None]*nmaps
    wav_fits_root = [None]*nmaps
    for i in xrange(len(scal_fits)):
        scal_fits[i] = fitsdir + fitsroot + fitscode[i] + '_scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_fits_root[i] = fitsdir + fitsroot + fitscode[i] + '_wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    outdir = fitsdir
    outroot = 's2let_ilc_' + fitsroot
    scal_outfits = outdir + outroot + 'scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_outfits_root = outdir + outroot + 'wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    #outputdir = 's2let_ilc_output/' #Subdirectory for secondary output

    #Construct valid hybrid tiling
    scal_tiles, wav_tiles, scal_bandlims, wav_bandlims, jmax, l_bounds = ps.construct_hybrid_tiling(ellmax,jmin,lamdas,l_transitions)
    res = ps.verify_tiling(ellmax,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    if res == False:
        raise ValueError('\nAn invalid wavelet tiling has been chosen.\n')
    else:
        print '\nA valid wavelet tiling has been chosen.\n'

    #Run ILC on scaling function map
    scal_output = s2let_ilc((-1,-1)) #(j,n) = (-1,-1) for scaling function

    #Run ILC on wavelet maps in PARALLEL
    jmin_real = jmax - 1
    jmax_real = jmax - 1
    ndir_min = 0
    ndir_max = 0

    for j in xrange(jmin_real,jmax_real+1): #Loop over scales
        k = 0
        mapsextra = [None]*(ndir_max+1-ndir_min)
        for n in xrange(ndir_min,ndir_max+1): #Loop over directions
            mapsextra[k] = (j,n) #Map loading within sub-process
            k+=1
        print "\nForming non-daemonic pool"
        pool = MyPool(nprocess)
        print "Farming out workers to run S2LET ILC on wavelet scales\n"
        wav_output = pool.map(s2let_ilc,mapsextra)
        pool.close()
        pool.join()


