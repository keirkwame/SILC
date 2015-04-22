import numpy as np
import healpy as hp
import math as mh
import multiprocessing as mg
import pys2let as ps

def smoothworker(i):
    print "Smoothing another independent covariance element"
    alms = ps.map2alm_mw(i[0],i[1],i[2]) #No pixwin correct. with MW sampling - calc alms to smooth
    #del i[0] #Everything gets moved down one index
    hp.almxfl(alms,i[3],inplace=True) #Multiply by gaussian beam
    
    '''if i[5] != -1:
        print "Picking out directional component of covariance map"
        #Testing convolving gaussian-smoothed covariance maps with directional wavelet
        jmin_min = 1
        #Need to truncate alm's to scale_lmax
        alms_fname = 'alms_' + str(i[6]) + '.fits'
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
                    #wav[offset:offset+nelem] = 0.
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
    
    #hp.almxfl(alms,i[3],inplace=True) #Multiply by gaussian beam
    print "Synthesising smoothed covariance map"
    Rsmoothflat = ps.alm2map_mw(alms,i[1],i[2]) #Smooth covariance in MW - calc final map to scale
    del alms
    return Rsmoothflat

def doubleworker(i):
    print "Doubling l_max of another input map"
    alms = ps.map2alm_mw(i[0],i[1],i[3]) #alm's to l(j)
    #del i[0] #Everything gets moved down one index
    
    print "Zero-padding the alm's"
    nzeros = i[2] - i[1] #No. zeros to pad
    new_alms_temp = np.concatenate((alms[:i[1]],np.zeros(nzeros)))
    for em in xrange(1,i[1]):
        startindex = em*i[1] - .5*em*(em-1)
        new_alms_temp = np.concatenate((new_alms_temp,alms[startindex:(startindex+i[1]-em)],np.zeros(nzeros)))
    del alms
    print "Temporary length of alm's =", len(new_alms_temp)
    nfinalzeros = hp.Alm.getsize(i[2]-1) - len(new_alms_temp)
    new_alms_temp = np.concatenate((new_alms_temp,np.zeros(nfinalzeros)))
    print "Final length of alm's =", len(new_alms_temp)

    mapsdouble = ps.alm2map_mw(new_alms_temp,i[2],i[3])
    del new_alms_temp
    return mapsdouble

def s2let_ilc_dir_para(mapsextra): #mapsextra = (maps,scale_lmax,spin,n,j,i)
    print "\nRunning Directional S2LET ILC on wavelet scale", mapsextra[4], "/", jmax, "direction", mapsextra[3]+1, "/", ndir, "\n"
    nrows = len(mapsextra[0]) #No. rows in covar. matrix
    smoothing_lmax = 2.*mapsextra[1] #=4.*nside(j)
    
    #Doubling lmax for input maps with zero-padding
    '''pool = mg.Pool(nprocess)
    mapsextra = [(maps[i],scale_lmax,smoothing_lmax,spin) for i in xrange(nrows)]
    del maps
    mapsdouble = np.array(pool.map(doubleworker,mapsextra))
    del mapsextra'''
    #Serial version
    mapsdouble = np.zeros((nrows,ps.mw_size(smoothing_lmax)),dtype=np.complex128) #Pre-allocate array
    for i in xrange(nrows):
        mapsdouble[i,:] = doubleworker((mapsextra[0][i],mapsextra[1],smoothing_lmax,mapsextra[2]))
    #mapsdouble = np.array(mapsdouble)
    
    #Calculating covariance matrix (at each pixel)
    #R = [None]*len(mapsdouble)
    R = np.zeros((len(mapsdouble),len(mapsdouble),len(mapsdouble[0])),dtype=np.complex128) #Pre-allocate array
    for i in xrange(len(mapsdouble)):
        R[i,:,:] = np.multiply(mapsdouble,np.roll(mapsdouble,-i,axis=0))
    #R = np.array(R)
    
    #Calculate scale_fwhm & smoothing_lmax
    nsamp = 1200.
    npix = hp.nside2npix(0.5*mapsextra[1]) #Equivalent number of HEALPix pixels
    scale_fwhm = 4. * mh.sqrt(nsamp / npix)
    
    #Smooth covariance matrices
    nindepelems = int(nrows*(nrows+1)*.5) #No. independent elements in symmetric covariance matrix
    Rflat = np.reshape(R,(nrows*nrows,len(R[0,0]))) #Flatten first two axes
    del R #NEW!!!
    Rflatlen = len(Rflat)
    gausssmooth = hp.gauss_beam(scale_fwhm,smoothing_lmax-1)
    
    #Testing zero-ing gaussian smoothing beam
    gauss_lmax = mapsextra[1]
    gausssmooth[gauss_lmax:] = 0.
    
    '''alms = [None]*nindepelems
    #alms_hp = [None]*nindepelems
    #alms_smooth = [None]*nindepelems
    Rsmoothflat = [None]*nindepelems #Only really need to pre-allocate this
    for i in xrange(nindepelems): #PARALLELISE
        print "Smoothing independent covariance element", i+1, "/", nindepelems
        alms[i] = ps.map2alm_mw(Rflat[i],scale_lmax,spin) #No pixwin correct. with MW sampling
        #alms_hp[i] = ps.lm2lm_hp(alms[i],smoothing_lmax) #Now in healpy ordering
        hp.almxfl(alms[i],gausssmooth,inplace=True) #Multiply by gaussian beam
        #alms_smooth[i] = ps.lm_hp2lm(alms_hp[i],smoothing_lmax) #Back in MW ordering
        Rsmoothflat[i] = ps.alm2map_mw(alms[i],scale_lmax,spin) #Smooth covariance in MW
    Rsmoothflat = np.array(Rsmoothflat)'''
    #Parallel version
    '''pool = mg.Pool(nprocess)
    Rflatextra = [(Rflat[i],smoothing_lmax,spin,gausssmooth,scale_lmax,en,i) for i in xrange(nindepelems)]
    del Rflat
    Rsmoothflat = np.array(pool.map(smoothworker,Rflatextra))
    del Rflatextra'''
    #Serial version
    Rsmoothflat = np.zeros_like(Rflat) #Pre-allocate array
    for i in xrange(nindepelems):
        Rsmoothflat[i,:] = smoothworker((Rflat[i],smoothing_lmax,mapsextra[2],gausssmooth,mapsextra[1],mapsextra[3],i))
    del Rflat
    #Rsmoothflat = np.array(Rsmoothflat)

    #Rearranging and padding out elements of Rsmooth
    Rsmoothflat[:nrows] = 0.5*Rsmoothflat[:nrows] #Multiply diag elements by half- not double-count
    Rsmoothflat = np.vstack((Rsmoothflat,np.zeros((Rflatlen-len(Rsmoothflat),len(Rsmoothflat[0]))))) #Zero-pad
    Rsmoothfat = np.reshape(Rsmoothflat,(nrows,nrows,len(Rsmoothflat[0]))) #Reshape Rsmooth as mat.
    del Rsmoothflat
    for i in xrange(1,len(Rsmoothfat[0])):
        Rsmoothfat[:,i,:] = np.roll(Rsmoothfat[:,i,:],i,axis=0) #Now in correct order-but with gaps
    Rsmoothfat = Rsmoothfat + np.transpose(Rsmoothfat,axes=(1,0,2)) #Gaps filled in

    #Compute inverse covariance matrices
    Rinv = np.linalg.inv(np.transpose(Rsmoothfat,axes=(2,0,1))) #Parallel vers. actually slower!?
    del Rsmoothfat

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
    #hp.write_alm('alms.fits',finalmapalms,lmax=mapsextra[1]-1,mmax=mapsextra[1]-1)
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

    return

##Input
nprocess = 4
nmaps = 5 #No. maps (WMAP = 5) (Planck = 9)
ellmax = 1024 #S2LET parameters - actually band-limits to 1 less
wavparam = 2
ndir = 1 #No. directions for each wavelet scale
spin = 0 #0 for temp, 1 for spin signals
upsample = 0 #0 for multiresolution, 1 for all scales at full resolution
jmin = 6
jmax = ps.pys2let_j_max(wavparam,ellmax,jmin)

fitsdir = 'deconv_data/'
fitsroot = 'wmap_deconv_smoothw_extrapolated_9yr_' #'planck_deconv_'
fitscode = ['k','ka','q','v','w'] #['30','44','70','100','143','217','353','545','857']
scal_fits = [None]*nmaps
wav_fits = [None]*nmaps
for i in xrange(len(scal_fits)):
    scal_fits[i] = fitsdir + fitsroot + fitscode[i] + '_scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
    wav_fits[i] = fitsdir + fitsroot + fitscode[i] + '_wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'

outdir = fitsdir
outroot = 's2let_ilc_dir_para_gauss_' + fitsroot
scal_outfits = outdir + outroot + 'scal_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir) + '.npy'
wav_outfits_root = outdir + outroot + 'wav_' + str(ellmax) + '_' + str(wavparam) + '_' + str(jmin) + '_' + str(ndir)

#Load scaling function maps for each channel
scal_maps = [None]*nmaps
for i in xrange(len(scal_maps)):
    scal_maps[i] = np.load(scal_fits[i])
scal_maps = np.array(scal_maps)

#Run ILC on scaling function map
scaling_lmax = wavparam**jmin
print "Running Directional S2LET ILC on scaling function"
scal_output = s2let_ilc_dir_para((scal_maps,scaling_lmax,spin,-1,-1,0)) #n,j=-1 signifies scaling func.
del scal_maps

#Load wavelet maps for each channel
wav_maps = [None]*nmaps
for i in xrange(len(wav_maps)):
    wav_maps[i] = np.load(wav_fits[i])
wav_maps = np.array(wav_maps) #1.1Gb!!!

#Run ILC on wavelet maps in PARALLEL
mapsextra = [None]*(jmax+1-(jmin))*(ndir)
i = 0
for j in xrange(jmin,jmax+1): #Loop over scales
    for n in xrange(0,ndir): #Loop over directions
        offset,scale_lmax,nelem,nelem_wav = ps.wav_ind(j,n,wavparam,ellmax,ndir,jmin,upsample)
        mapsextra[i] = (wav_maps[:,offset:offset+nelem],scale_lmax,spin,n,j,i)
        i += 1
del wav_maps

pool = mg.Pool(nprocess)
wav_output = pool.map_async(s2let_ilc_dir_para,mapsextra)
del mapsextra
