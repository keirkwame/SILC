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
    print "Smoothing another independent covariance element", i[2], ",", i[3]

    #Map loading within sub-process
    wav_fits1 = wav_fits_root[i[2]] + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '_double.npy'
    wav_fits2 = wav_fits_root[i[3]] + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '_double.npy'
    map1 = np.real(np.load(wav_fits1,mmap_mode='r')) #Throw away zero imaginary part
    map2 = np.real(np.load(wav_fits2,mmap_mode='r'))
    R = np.multiply(map1,map2) + 0.j #Add back in zero imaginary part
    del map1,map2

    alms = ps.map2alm_mw(R,i[4],spin) #No pixwin correct. with MW - calc alms to smooth
    del R
    gausssmooth = hp.gauss_beam(i[5],lmax=i[4]-1)
    hp.almxfl(alms,gausssmooth,inplace=True) #Multiply by gaussian beam
    
    '''if i[5] != -1: #If n != -1 (i.e. the maps are directional)
        print "Picking out directional component of covariance map"
        #Testing convolving gaussian-smoothed covariance maps with directional wavelet
        jmin_min = 1
        #Need to truncate alm's to scale_lmax (Is this necessary?)
        alms_fname = 'alms_dirwav_' + str(i[7]) + '_' + str(i[5]) + '_' + str(i[6]) + '.fits'
        hp.write_alm(alms_fname,alms,lmax=i[4]-1,mmax=i[4]-1)
        del alms
        alms_truncate = hp.read_alm(alms_fname)
        print "Analysing covariance map" #Could increase wavparam!?
        wav,scal = ps.analysis_lm2wav(alms_truncate,wavparam,i[4],jmin_min,ndir,spin,upsample)
        del alms_truncate
        #Delete wrong directions by zero-ing them
        print "Deleting wrong directions"
        jmax_min = ps.pys2let_j_max(wavparam,i[4],jmin_min)
        for j in xrange(jmin_min,jmax_min+1):
            for n in xrange(0,ndir):
                if n != i[5]:
                    offset,new_scale_lmax,nelem,nelem_wav = ps.wav_ind(j,n,wavparam,i[4],ndir,jmin_min,upsample)
                    wav[offset:offset+nelem] = 0.
        print "Synthesising directional covariance map"
        alms = ps.synthesis_wav2lm(wav,scal,wavparam,i[4],jmin_min,ndir,spin,upsample)
        del wav,scal
        #Expand alm's with zero-padding
        print "Zero-padding the alm's"
        nzeros = i[1] - i[4] #No. zeros to pad
        new_alms_temp = np.concatenate((alms[:i[4]],np.zeros(nzeros)))
        for em in xrange(1,i[4]):
            startindex = em*i[4] - .5*em*(em-1)
            new_alms_temp = np.concatenate((new_alms_temp,alms[startindex:(startindex+i[4]-em)],np.zeros(nzeros)))
        del alms
        print "Temporary length of alm's =", len(new_alms_temp)
        nfinalzeros = hp.Alm.getsize(i[1]-1) - len(new_alms_temp)
        alms = np.concatenate((new_alms_temp,np.zeros(nfinalzeros)))
        del new_alms_temp
        print "Final length of alm's =", len(alms)'''

    print "Synthesising smoothed covariance map for element", i[2], ",", i[3]
    Rsmooth = np.real(ps.alm2map_mw(alms,i[4],spin)) #Throw away zero imaginary part
    del alms
    
    #SAVE smoothed covariance
    R_fits = wav_outfits_root + '_Rsmooth' + str(i[2]) + str(i[3]) +'.npy'
    np.save(R_fits,Rsmooth)
    del Rsmooth
    
    return 0

def zeropad(i): #(alms,scale_lmax,smoothing_lmax)
    print "Zero-padding the alm's"

    '''nzeros = i[2] - i[1] #No. zeros to pad at each m
    new_alms_temp = np.concatenate((i[0][:i[1]],np.zeros(nzeros)))
    for em in xrange(1,i[1]):
        startindex = em*i[1] - .5*em*(em-1)
        new_alms_temp = np.concatenate((new_alms_temp,i[0][startindex:(startindex+i[1]-em)],np.zeros(nzeros)))
    nfinalzeros = hp.Alm.getsize(i[2]-1) - len(new_alms_temp)
    new_alms_temp = np.concatenate((new_alms_temp,np.zeros(nfinalzeros)))
    print "Final length of alm's =", len(new_alms_temp)'''

    #Pre-allocate array for enlarged alm's
    new_alms = np.zeros(((i[2]*(i[2]+1))/2.),dtype=complex)
    for em in xrange(i[1]):
        startindex = em*i[1] - .5*em*(em-1)
        new_startindex = em*i[2] - .5*em*(em-1)
        new_alms[new_startindex:(new_startindex+i[1]-em)] = i[0][startindex:(startindex+i[1]-em)]

    return new_alms

def doubleworker(i): #i = (j,n,map_index,scale_lmax,smoothing_lmax)
    print "Doubling l_max of input map", i[2]+1, "/", nmaps
    
    #TESTING map loading within sub-process
    wav_fits = wav_fits_root[i[2]] + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '.npy'
    map = np.load(wav_fits,mmap_mode='r') #Map still only stored on disk

    alms = ps.map2alm_mw(map,i[3],spin) #alm's to l(j)
    del map
    alms = zeropad((alms,i[3],i[4])) #New alm's larger
    map = ps.alm2map_mw(alms,i[4],spin)
    del alms
    
    #SAVE doubled map
    double_fits = wav_fits_root[i[2]] + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '_double.npy'
    np.save(double_fits,map)
    del map
    
    return 0

def s2let_ilc_dir_para(mapsextra): #mapsextra = (j,n)
    print "\nRunning Directional S2LET ILC on wavelet scale", mapsextra[0], "/", jmax, "direction", mapsextra[1]+1, "/", ndir, "\n"

    scale_lmax = wavparam**(mapsextra[0]+1) #lambda^(j+1)
    if scale_lmax > ellmax:
        scale_lmax = ellmax
    smoothing_lmax = 2.*(scale_lmax-1.)+1

    #Doubling lmax for input maps with zero-padding
    #Serial version
    '''mapsdouble = np.zeros((nrows,ps.mw_size(smoothing_lmax)),dtype=np.complex128) #Pre-allocate array
    for i in xrange(nrows):
        mapsdouble[i,:] = doubleworker((mapsextra[0][i],mapsextra[1],smoothing_lmax,mapsextra[2]))'''
    #Parallel version
    mapsextra2 = [(mapsextra[0],mapsextra[1],i,scale_lmax,smoothing_lmax) for i in xrange(nmaps)]
    
    print "Forming pool"
    pool2 = mg.Pool(nprocess2)
    print "\nFarming out workers to run doubling function"
    double_output = pool2.map(doubleworker,mapsextra2)
    print "Have returned from doubling workers\n"
    pool2.close()
    pool2.join()
    del pool2

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
    print "\nFarming out workers to run smoothing function"
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
            R_fits = wav_outfits_root + '_Rsmooth' + str(i) + str(j) +'.npy'
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
        wav_fits = wav_fits_root[i] + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '_double.npy'
        mapsdouble[:,i] = np.real(np.load(wav_fits,mmap_mode='r')) #Throw away zero imaginary part

    #Dot weights with maps (at each small pixel) - at double l(j)
    finalmap = np.sum(np.multiply(wk,mapsdouble),axis=-1) + 0.j #Add back in zero imaginary part
    del wk,mapsdouble
    
    #Downgrade resolution of MW maps
    print "Downgrading resolution of CMB wavelet map"
    finalmapalms = ps.map2alm_mw(finalmap,smoothing_lmax,spin)
    del finalmap
    alms_fname = wav_outfits_root + '_alms_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.fits'
    hp.write_alm(alms_fname,finalmapalms,lmax=scale_lmax-1,mmax=scale_lmax-1)
    del finalmapalms
    finalmapalmstruncate = hp.read_alm(alms_fname)
    finalmaphalf = ps.alm2map_mw(finalmapalmstruncate,scale_lmax,spin)
    del finalmapalmstruncate
    
    #Saving output map
    wav_outfits = wav_outfits_root + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.npy'
    if mapsextra[0] == -1:
        wav_outfits = scal_outfits
    np.save(wav_outfits,finalmaphalf)
    del finalmaphalf

    return 0

def test_ilc(mapsextra):
    #Saving output map
    wav_outfits = wav_outfits_root + '_j' + str(mapsextra[2]) + '_n' + str(mapsextra[3]+1) + '.npy'
    if mapsextra[2] == -1:
        wav_outfits = scal_outfits
    np.save(wav_outfits,mapsextra[0].value[2])

    return 0

if __name__ == "__main__":
    ##Input
    nmaps = 9 #No. maps (WMAP = 5) (Planck = 9)
    ellmax = 256 #S2LET parameters - actually band-limits to 1 less
    wavparam = 2
    ndir = 1 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals
    upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
    jmin = 6
    jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

    fitsdir = '/Users/keir/Documents/s2let_ilc_planck/deconv_data/' #'/home/keir/s2let_ilc_data/'
    fitsroot = 'planck_deconv_tapered_' #'ffp6_combined_mc_0000_deconv_' #'simu_dirty_beam_wmap_9yr_' #'wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    fitscode = ['30','44','70','100','143','217','353','545','857'] #['k','ka','q','v','w']
    scal_fits = [None]*nmaps
    wav_fits_root = [None]*nmaps
    for i in xrange(len(scal_fits)):
        scal_fits[i] = fitsdir + fitsroot + fitscode[i] + '_scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_fits_root[i] = fitsdir + fitsroot + fitscode[i] + '_wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir)

    outdir = fitsdir
    outroot = 's2let_ilc_dir_hypatia_memeff_' + fitsroot
    scal_outfits = outdir + outroot + 'scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_outfits_root = outdir + outroot + 'wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir)

    #Load scaling function maps for each channel
    '''scal_maps = [None]*nmaps
    for i in xrange(len(scal_maps)):
        print "Loading scaling function maps for channel", i+1, "/", len(scal_maps)
        scal_maps[i] = np.load(scal_fits[i])
    scal_maps = np.array(scal_maps)

    #Run ILC on scaling function map
    nprocess2 = 9
    nprocess3 = 45

    scaling_lmax = wavparam**jmin
    print "\nRunning Directional S2LET ILC on scaling function"
    scal_output = s2let_ilc_dir_para((ForkedData(scal_maps),scaling_lmax,-1,-1,spin,0)) #j,n=-1 signifies scaling func.
    print "Finished running Directional S2LET ILC on scaling function"
    del scal_maps,scal_output'''

    #Run ILC on wavelet maps in PARALLEL
    nprocess = 1
    nprocess2 = 3
    nprocess3 = 4

    jmin_real = jmin
    jmax_real = jmax
    ndir_min = 0
    ndir_max = ndir - 1

    for j in xrange(jmin_real,jmax_real+1): #Loop over scales
        mapsextra = [None]*ndir
        for n in xrange(ndir_min,ndir_max+1): #Loop over directions
            #mapsextra[i] = (ForkedData(wav_maps[:,offset:offset+nelem]),scale_lmax,j,n,spin,i)
            mapsextra = [(j,n)] #TESTING map loading within sub-process
        print "\nForming non-daemonic pool"
        pool = MyPool(nprocess)
        print "Farming out workers to run Directional S2LET ILC on wavelet scales"
        wav_output = pool.map(s2let_ilc_dir_para,mapsextra)
        pool.close()
        pool.join()


