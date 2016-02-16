import numpy as np
import healpy as hp
import math as mh
import multiprocessing as mg
import pys2let as ps

def smoothworker(i): #(j,n,map_index1,map_index2,smoothing_lmax,scale_fwhm) [map_index2<=map_index1]
    print "Smoothing independent covariance element", i[2], ",", i[3]

    #Map loading within sub-process
    if i[0] >= 0: #Wavelet scales
        wav_fits1 = wav_fits_root[i[2]] + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '.npy'
        wav_fits2 = wav_fits_root[i[3]] + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '.npy'
    else: #Scaling function
        wav_fits1 = scal_fits[i[2]]
        wav_fits2 = scal_fits[i[3]]
    map1 = np.load(wav_fits1,mmap_mode='r') #Complex spin-wavelet coefficients
    map2 = np.conjugate(np.load(wav_fits2,mmap_mode='r'))
    R = np.multiply(map1,map2) #W_c W_d^*
    del map1,map2

    alms = ps.lm2lm_hp(ps.map2alm_mw(R,i[4],0),i[4]) #No pixwin correct. with MW - calc alms to smooth - come out in MW order - so converted to HPX order
    del R
    gausssmooth = hp.gauss_beam(i[5],lmax=i[4]-1)
    hp.almxfl(alms,gausssmooth,inplace=True) #Multiply by gaussian beam

    print "Synthesising smoothed covariance map for element", i[2], ",", i[3]
    Rsmooth = ps.alm2map_mw(ps.lm_hp2lm(alms,i[4]),i[4],0) #Input alm's in MW order
    del alms
    
    #Save smoothed covariance
    if i[0] >= 0: #Wavelet scales
        R_fits = wav_outfits_root + '_j' + str(i[0]) + '_n' + str(i[1]+1) + '_Rsmooth' + str(i[2]) + str(i[3]) +'.npy'
    else: #Scaling function
        R_fits = scal_outfits[:-4] + '_Rsmooth' + str(i[2]) + str(i[3]) +'.npy'
    np.save(R_fits,Rsmooth)
    del Rsmooth
    
    return 0

def s2let_ilc(mapsextra): #mapsextra = (j,n)
    print "Running Spin-SILC on wavelet scale", mapsextra[0], "/", jmax, "direction", mapsextra[1]+1, "/", ndir, "\n"

    if mapsextra[0] >= 0: #Wavelet scales
        scale_lmax = wav_bandlims[mapsextra[0]] #smoothing_lmax
    else: #Scaling function
        scale_lmax = int(scal_bandlims)
    real_lmax = 0.5*(scale_lmax-1.) + 1

    #Calculate scale_fwhm for smoothing kernel
    nsamp = 1200.
    npix = ps.mw_size(real_lmax)
    scale_fwhm = 10. * mh.sqrt(nsamp / npix) #Original was factor of 50

    #Smooth covariance matrices
    nindepelems = int(nmaps*(nmaps+1)*.5) #No. indep. elements in Hermitian covariance matrix
    Rextra = [None]*nindepelems
    k=0
    for i in xrange(nmaps):
        for j in xrange(i+1):
            Rextra[k] = (mapsextra[0],mapsextra[1],i,j,scale_lmax,scale_fwhm)
            k+=1
    print "Forming pool"
    pool2 = mg.Pool(nprocess2)
    print "Farming out workers to run smoothing function"
    R_output = pool2.map(smoothworker,Rextra)
    print "Have returned from smoothing workers\n"
    pool2.close()
    pool2.join()
    del pool2

    #Load R maps and form matrices
    print "Pre-allocating memory for complete covariance tensor\n"
    Rsmooth = np.zeros((ps.mw_size(scale_lmax),nmaps,nmaps),dtype=np.complex128) #Pre-allocate array
    for i in xrange(nmaps):
        for j in xrange(i+1):
            if mapsextra[0] >= 0: #Wavelet scales
                R_fits = wav_outfits_root + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '_Rsmooth' + str(i) + str(j) +'.npy'
            else: #Scaling function
                R_fits = scal_outfits[:-4] + '_Rsmooth' + str(i) + str(j) +'.npy'
            Rsmooth[:,i,j] = np.load(R_fits) #Bottom triangle
            if i != j:
                Rsmooth[:,j,i] = np.conj(Rsmooth[:,i,j]) #Matrices are Hermitian

    #Compute inverse covariance matrices & weights
    print "Calculating inverse covariance matrices\n"
    Rinv = np.linalg.inv(Rsmooth) #OPTIMISE INVERSE OF HERMITIAN MATRIX
    del Rsmooth
    wknumer = np.sum(Rinv,axis=1) #Sum over index c
    del Rinv
    wkdenom = np.sum(wknumer,axis=-1) #Sum over index d
    wk = wknumer / wkdenom[:,None]
    del wknumer,wkdenom

    #Map loading
    mapsdouble = np.zeros((len(wk),len(wk[0])),dtype=np.complex128) #Pre-allocate array
    for i in xrange(nmaps):
        if mapsextra[0] >= 0: #Wavelet scales
            wav_fits = wav_fits_root[i] + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.npy'
        else: #Scaling function
            wav_fits = scal_fits[i]
        mapsdouble[:,i] = np.load(wav_fits,mmap_mode='r') #Complex spin-wavelet coefficients

    #Dot weights with maps (at each small pixel) - at double l(j)
    finalmap = np.sum(np.multiply(wk,mapsdouble),axis=-1)
    del wk,mapsdouble

    #Saving output map
    if mapsextra[0] >= 0: #Wavelet scales
        wav_outfits = wav_outfits_root + '_j' + str(mapsextra[0]) + '_n' + str(mapsextra[1]+1) + '.npy'
    else: #Scaling function
        wav_outfits = scal_outfits
    np.save(wav_outfits,finalmap)
    del finalmap

    return 0

if __name__ == "__main__":
    ##Input
    #Set directory structure
    comp = 0
    
    if comp == 0: #Keir's iMac
        #nprocess = 1 #For distributing directions within scale
        nprocess2 = 4 #For smoothing covariance elements
        fitsdir = '/Users/keir/Documents/spin_silc/wavelet_maps/'
        outdir = '/Users/keir/Documents/spin_silc/ilc_maps/'
    elif comp == 1: #Hypatia
        #nprocess = 1 #For distributing directions within scale
        nprocess2 = 28 #For smoothing covariance elements
        fitsdir = ''
        outdir = ''
    
    nmaps = 7
    ellmax = 917
    jmin = 6
    lamdas = np.array([2,1.3]) #60,2])
    l_transitions = np.array([513]) #61])
    wavparam_code = 'TestD'
    ndir = 1 #No. directions for each wavelet scale
    spin = 2

    #Construct valid hybrid tiling
    scal_tiles, wav_tiles, scal_bandlims, wav_bandlims, jmax, l_bounds = ps.construct_hybrid_tiling(ellmax,jmin,lamdas,l_transitions)

    #Remove truncated wavelets
    ntrunc = 1 #1 or 2
    ellmax = wav_bandlims[-1-1*ntrunc] #Lower total bandlimit to match smallest wavelet
    print "New bandlimit =", ellmax
    scal_tiles = scal_tiles[:ellmax]
    wav_tiles = wav_tiles[:ellmax,:-1*ntrunc]
    wav_bandlims = wav_bandlims[:-1*ntrunc]
    jmax = jmax - ntrunc
    print jmax

    #Zeropad wavelets to allow double resolution
    smoothing_lmax = 2*(ellmax-1) + 1
    print "l_max = %i, smoothing l_max = %i" %(ellmax,smoothing_lmax)
    scal_tiles = np.concatenate((scal_tiles,np.zeros(smoothing_lmax-ellmax)))
    print len(scal_tiles)
    wav_tiles = np.concatenate((wav_tiles,np.zeros((smoothing_lmax-ellmax,len(wav_tiles[0])))))
    print len(wav_tiles), len(wav_tiles[0])
    scal_bandlims = 2*(scal_bandlims-1) + 1
    print scal_bandlims
    wav_bandlims = 2*(wav_bandlims-1) + 1
    print wav_bandlims

    res = ps.verify_tiling(smoothing_lmax,scal_tiles,wav_tiles.T.ravel(),scal_bandlims,wav_bandlims)
    if res == False:
        print '\nCheck that above admissibility condition is consistent with zero-padding of wavelets\n'
        pass
    else:
        print '\nA valid wavelet tiling has been chosen.\n'

    fitsroot = 'planck_pol_diffuse_deconv_'
    fitscode = ['30','44','70','100','143','217','353']
    scal_fits = [None]*nmaps
    wav_fits_root = [None]*nmaps
    for i in xrange(len(scal_fits)):
        scal_fits[i] = fitsdir + fitsroot + fitscode[i] + '_scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
        wav_fits_root[i] = fitsdir + fitsroot + fitscode[i] + '_wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    outroot = 'spin_silc_fwhm10_' + fitsroot
    scal_outfits = outdir + outroot + 'scal_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_outfits_root = outdir + outroot + 'wav_' + str(ellmax) + '_hybrid' + wavparam_code + '_' + str(jmin) + '_' + str(ndir)

    #Run ILC on scaling function map
    scal_output = s2let_ilc((-1,-1)) #(j,n) = (-1,-1) for scaling function

    #Run ILC on wavelet maps in PARALLEL
    jmin_real = 0 #!= jmin
    jmax_real = jmax
    ndir_min = 0
    ndir_max = ndir - 1

    for j in xrange(jmin_real,jmax_real+1): #Loop over scales
        for n in xrange(ndir_min,ndir_max+1): #Loop over directions
            #mapsextra = [(j,n)] #Map loading within sub-process
            #print "\nForming non-daemonic pool"
            #pool = MyPool(nprocess)
            #print "Farming out workers to run Spin-SILC on wavelet scales\n"
            #silc_output = pool.map(s2let_ilc,mapsextra)
            silc_output = s2let_ilc((j,n))
            #pool.close()
            #pool.join()
