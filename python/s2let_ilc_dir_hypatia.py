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

def smoothworker(i): #(Rflat[i],smoothing_lmax,spin,gausssmooth,scale_lmax,n,i,j)
    print "Smoothing another independent covariance element"
    alms = ps.map2alm_mw(i[0],i[1],i[2]) #No pixwin correct. with MW - calc alms to smooth
    hp.almxfl(alms,i[3],inplace=True) #Multiply by gaussian beam
    
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

    print "Synthesising smoothed covariance map"
    Rsmoothflat = ps.alm2map_mw(alms,i[1],i[2]) #Smooth covar in MW - calc final map to scale
    del alms
    return Rsmoothflat

def zeropad(i): #(alms,scale_lmax,smoothing_lmax,spin)
    print "Zero-padding the alm's"
    nzeros = i[2] - i[1] #No. zeros to pad at each m
    new_alms_temp = np.concatenate((i[0][:i[1]],np.zeros(nzeros)))
    for em in xrange(1,i[1]):
        startindex = em*i[1] - .5*em*(em-1)
        new_alms_temp = np.concatenate((new_alms_temp,i[0][startindex:(startindex+i[1]-em)],np.zeros(nzeros)))
    #print "Temporary length of alm's =", len(new_alms_temp)
    nfinalzeros = hp.Alm.getsize(i[2]-1) - len(new_alms_temp)
    new_alms_temp = np.concatenate((new_alms_temp,np.zeros(nfinalzeros)))
    print "Final length of alm's =", len(new_alms_temp)
    return new_alms_temp

def doubleworker(i): #i = (maps[i],scale_lmax,smoothing_lmax,spin)
    print "Doubling l_max of another input map"
    alms = ps.map2alm_mw(i[0],i[1],i[3]) #alm's to l(j)
    alms = zeropad((alms,i[1],i[2],i[3]))
    mapsdouble = ps.alm2map_mw(alms,i[2],i[3])
    del alms
    return mapsdouble

def s2let_ilc_dir_para(mapsextra): #mapsextra = (maps,scale_lmax,spin,n,j,i)
    print "\nRunning Directional S2LET ILC on wavelet scale", mapsextra[4], "/", jmax, "direction", mapsextra[3]+1, "/", ndir, "\n"
    nrows = len(mapsextra[0]) #No. rows in covar. matrix
    smoothing_lmax = 2.*mapsextra[1] #=4.*nside(j)
    
    #Doubling lmax for input maps with zero-padding
    #Serial version
    '''mapsdouble = np.zeros((nrows,ps.mw_size(smoothing_lmax)),dtype=np.complex128) #Pre-allocate array
    for i in xrange(nrows):
        mapsdouble[i,:] = doubleworker((mapsextra[0][i],mapsextra[1],smoothing_lmax,mapsextra[2]))'''
    #Parallel version
    print "Started packing input data for doubling"
    mapsextra2 = [(mapsextra[0][i],mapsextra[1],smoothing_lmax,mapsextra[2]) for i in xrange(nrows)]
    print "Finished packing input data for doubling"
    nprocess2 = 9
    print "Here"
    pool2 = mg.Pool(nprocess2)
    print "Farming out workers to run doubling function"
    mapsdouble = np.array(pool2.map(doubleworker,mapsextra2))
    print "Have returned from doubling workers"
    pool2.close()
    pool2.join()
    del mapsextra2
    
    #Calculating covariance matrix (at each pixel)
    print "Calculating covariance matrices"
    R = np.zeros((len(mapsdouble),len(mapsdouble),len(mapsdouble[0])),dtype=np.complex128) #Pre-allocate array
    for i in xrange(len(mapsdouble)):
        R[i,:,:] = np.multiply(mapsdouble,np.roll(mapsdouble,-i,axis=0))

    #Calculate scale_fwhm for smoothing kernel
    nsamp = 1200.
    npix = 12*((0.5*mapsextra[1])**2) #Equivalent number of HEALPix pixels
    scale_fwhm = 4. * mh.sqrt(nsamp / npix)
    
    #Smooth covariance matrices
    nindepelems = int(nrows*(nrows+1)*.5) #No. indep. elements in symmetric covariance matrix
    R = np.reshape(R,(nrows*nrows,len(R[0,0]))) #Flatten first two axes
    Rflatlen = len(R)
    gausssmooth = hp.gauss_beam(scale_fwhm,smoothing_lmax-1)
    #Serial version
    '''Rsmoothflat = np.zeros_like(Rflat) #Pre-allocate array
    for i in xrange(nindepelems):
        Rsmoothflat[i,:] = smoothworker((Rflat[i],smoothing_lmax,mapsextra[2],gausssmooth,mapsextra[1],mapsextra[3],i,mapsextra[4]))
    del Rflat'''
    #Parallel version
    print "Started packing input data for smoothing"
    Rextra = [(R[i],smoothing_lmax,mapsextra[2],gausssmooth,mapsextra[1],mapsextra[3],i) for i in xrange(nindepelems)]
    print "Finished packing input data for doubling"
    del R
    nprocess3 = 23
    pool3 = mg.Pool(nprocess3)
    Rsmooth = np.array(pool3.map(smoothworker,Rextra))
    pool3.close()
    pool3.join()
    del Rextra

    #Rearranging and padding out elements of Rsmooth
    Rsmooth[:nrows] = 0.5*Rsmooth[:nrows] #Mult diag elems by half - not double-count
    Rsmooth = np.vstack((Rsmooth,np.zeros((Rflatlen-len(Rsmooth),len(Rsmooth[0]))))) #Zero-pad
    Rsmooth = np.reshape(Rsmooth,(nrows,nrows,len(Rsmooth[0]))) #Reshape Rsmooth as mat.
    for i in xrange(1,len(Rsmooth[0])):
        Rsmooth[:,i,:] = np.roll(Rsmooth[:,i,:],i,axis=0) #Now in correct order-but with gaps
    Rsmooth = Rsmooth + np.transpose(Rsmooth,axes=(1,0,2)) #Gaps filled in

    #Compute inverse covariance matrices
    print "Calculating inverse covariance matrices"
    Rinv = np.linalg.inv(np.transpose(Rsmooth,axes=(2,0,1))) #Parallel vers. slower!?
    del Rsmooth

    #Compute weights vectors (at each pixel)
    wknumer = np.sum(Rinv,axis=-1)
    del Rinv
    wkdenom = np.sum(wknumer,axis=-1)
    wk = wknumer / wkdenom[:,None]
    del wknumer,wkdenom
    
    #Dot weights with maps (at each small pixel) - at double l(j)
    finalmap = np.sum(np.multiply(wk,mapsdouble.T),axis=-1)
    del wk,mapsdouble
    
    #Downgrade resolution of MW maps
    print "Downgrading resolution of CMB wavelet map"
    finalmapalms = ps.map2alm_mw(finalmap,smoothing_lmax,mapsextra[2])
    del finalmap
    alms_fname = 'alms_' + str(mapsextra[5]) + '.fits'
    hp.write_alm(alms_fname,finalmapalms,lmax=mapsextra[1]-1,mmax=mapsextra[1]-1)
    del finalmapalms
    finalmapalmstruncate = hp.read_alm(alms_fname)
    finalmaphalf = ps.alm2map_mw(finalmapalmstruncate,mapsextra[1],mapsextra[2])
    del finalmapalmstruncate
    
    #Saving output map
    wav_outfits = wav_outfits_root + '_j' + str(mapsextra[4]) + '_n' + str(mapsextra[3]+1) + '.npy'
    if mapsextra[4] == -1:
        wav_outfits = scal_outfits
    np.save(wav_outfits,finalmaphalf)
    del finalmaphalf

    return 0

if __name__ == "__main__":
    ##Input
    nmaps = 9 #No. maps (WMAP = 5) (Planck = 9)
    ellmax = 3400 #S2LET parameters - actually band-limits to 1 less
    wavparam = 2
    ndir = 2 #No. directions for each wavelet scale
    spin = 0 #0 for temp, 1 for spin signals
    upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
    jmin = 6
    jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

    fitsdir =  '/home/keir/s2let_ilc_data/' #'/Users/keir/Documents/s2let_ilc_planck/deconv_data/'
    fitsroot = 'planck_deconv_' #'simu_dirty_beam_wmap_9yr_' #'wmap_deconv_nosource_smoothw_extrapolated_9yr_'
    fitscode = ['30','44','70','100','143','217','353','545','857'] #['k','ka','q','v','w']
    scal_fits = [None]*nmaps
    wav_fits = [None]*nmaps
    for i in xrange(len(scal_fits)):
        scal_fits[i] = fitsdir + fitsroot + fitscode[i] + '_scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_fits[i] = fitsdir + fitsroot + fitscode[i] + '_wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'

    outdir = fitsdir
    outroot = 's2let_ilc_dir_hypatia_' + fitsroot
    scal_outfits = outdir + outroot + 'scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_outfits_root = outdir + outroot + 'wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir)

    #Load scaling function maps for each channel
    '''scal_maps = [None]*nmaps
    for i in xrange(len(scal_maps)):
        print "Loading scaling function maps for channel", i+1, "/", len(scal_maps)
        scal_maps[i] = np.load(scal_fits[i])
    scal_maps = np.array(scal_maps)

    #Run ILC on scaling function map
    scaling_lmax = wavparam**jmin
    print "Running Directional S2LET ILC on scaling function"
    scal_output = s2let_ilc_dir_para((scal_maps,scaling_lmax,spin,-1,-1,0)) #n,j=-1 signifies scaling func.
    print "Finished running Directional S2LET ILC on scaling function"
    del scal_maps'''

    #Load wavelet maps for each channel
    wav_maps = [None]*nmaps
    for i in xrange(len(wav_maps)):
        print "Loading wavelet maps for map", i+1, "/", len(wav_maps)
        wav_maps[i] = np.load(wav_fits[i],mmap_mode='r') #Maps still only stored on disk
    wav_maps = np.array(wav_maps) #1.1Gb!!!

    #Run ILC on wavelet maps in PARALLEL
    jmin_real = 6
    jmax_real = 6

    mapsextra = [None]*(jmax_real+1-(jmin_real))*2
    i = 0
    for j in xrange(jmin_real,jmax_real+1): #Loop over scales
        for n in xrange(0,2): #Loop over directions
            offset,scale_lmax,nelem,nelem_wav = ps.wav_ind(j,n,wavparam,ellmax,ndir,jmin,upsample)
            print "Forming input data structure for scale", j, "direction", n+1
            mapsextra[i] = ForkedData((wav_maps[:,offset:offset+nelem],scale_lmax,spin,n,j,i))
            i += 1
    del wav_maps

    nprocess = 2
    print "Forming non-pickleable data objects"
    #mapsextra_forked = ForkedData(mapsextra)
    print "Finished forming non-pickleable data objects"
    pool = MyPool(nprocess)
    print "Farming out workers to run Directional S2LET ILC on wavelet scales"
    wav_output = pool.map(s2let_ilc_dir_para,mapsextra)
    pool.close()
    pool.join()
    #wav_output = s2let_ilc_dir_para(mapsextra)
    del mapsextra
    del mapsextra_forked


