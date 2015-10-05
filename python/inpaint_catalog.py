import numpy as np
import healpy as hp
import copy as cp
import astropy.io.fits as af

if __name__ == "__main__":
    ncata = 9 #No. maps (WMAP = 5) (Planck = 9) (PCCS & PCCS2E = 15)
    nside = 2048
    npix = hp.nside2npix(nside)
    
    fwhms = np.radians(np.array([32.239,27.005,13.252,9.651,7.248,4.990,4.818,4.682,4.325,9.651,7.248,4.990,4.818,4.682,4.325]) / 60.) #http://wiki.cosmos.esa.int/planckpla2015/index.php/Effective_Beams - mean FWHM (arcmin converted to radians)
    radii = cp.deepcopy(fwhms) #Radii of point sources at each channel
    radii[:3] = 1.15*radii[:3] #LFI
    radii[3:] = 1.27*radii[3:] #HFI

    #Point source catalogue TXT file
    #psc = 'foreground_data/wmap_ptsrc_catalog_9yr_v5p1.txt'
    psc = [None]*ncata
    '''psc[0] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_030_R1.30.fits'
    psc[1] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_044_R1.30.fits'
    psc[2] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_070_R1.30.fits'
    psc[3] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_100_R1.20.fits'
    psc[4] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_143_R1.20.fits'
    psc[5] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_217_R1.20.fits'
    psc[6] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_353_R1.20.fits'
    psc[7] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_545_R1.20.fits'
    psc[8] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR1/PCCS_1.0/source_lists/COM_PCCS_857_R1.20.fits'
    '''
    
    psc[0] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_030_R2.04.fits'
    psc[1] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_044_R2.04.fits'
    psc[2] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_070_R2.04.fits'
    psc[3] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_100_R2.01.fits'
    psc[4] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_143_R2.01.fits'
    psc[5] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_217_R2.01.fits'
    psc[6] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_353_R2.01.fits'
    psc[7] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_545_R2.01.fits'
    psc[8] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_857_R2.01.fits'
    '''psc[9] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_100-excluded_R2.01.fits'
    psc[10] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_143-excluded_R2.01.fits'
    psc[11] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_217-excluded_R2.01.fits'
    psc[12] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_353-excluded_R2.01.fits'
    psc[13] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_545-excluded_R2.01.fits'
    psc[14] = '/Users/keir/Documents/s2let_ilc_planck/catalogs/PR2/PCCS_2.0/source_lists/COM_PCCS_857-excluded_R2.01.fits'
    '''

    #Counting total number of sources
    nsources = 0
    for i in xrange(ncata):
        nsources = nsources + len(af.open(psc[i])[1].data)
    T1 = [None]*nsources #Will be list of arrays
    source_mask = np.zeros(npix,dtype=np.int8) #Map of zeros - binary mask (1 if point source)
    array_mask = np.zeros(npix,dtype=np.uint16) #Map of indices (k+1) pointing to which T1 hole [0,65535]

    k = 0 #Real index
    for i in xrange(ncata):
        full_catalog = af.open(psc[i])[1].data
        coords = np.array([full_catalog['GLON'],full_catalog['GLAT']]) #GLON (degrees),GLAT (degrees)
        coords[1] = 90. - coords[1] #phi,theta
        coords = np.radians(coords) #Convert to radians  - ordered by phi

        for j in xrange(len(coords[0])): #Loop over point sources
            k_comb = k
            print "\nForming T vectors for point source", j+1, "/", len(coords[0]), "map", i, "k =", k
            T1temp = hp.query_disc(nside,hp.ang2vec(coords[1,j],coords[0,j]),radii[i])
            if np.sum(source_mask[T1temp]) > 0: #I.e. current hole overlaps another hole
                print "There is an overlap"
                old_holes = np.unique(array_mask[T1temp]) - 1 #Real indices of overlapping holes
                print "Indices of old_holes =", old_holes
                if len(old_holes) == 1:
                    print "New hole completely within old hole"
                    k_comb = old_holes[0]
                elif len(old_holes) == 2:
                    k_comb = old_holes[1]
                    T1[k_comb] = np.union1d(T1[k_comb],T1temp) #Union 1st old hole + new
                else: #If more than one old hole
                    k_comb = old_holes[1]
                    T1[k_comb] = np.union1d(T1[k_comb],T1temp) #Union 1st old hole + new
                    for m in xrange(len(old_holes)-2): #Loop over additional old holes
                        T1[k_comb] = np.union1d(T1[k_comb],T1[old_holes[m+2]])
            else: #I.e. completely new hole
                print "There is a new hole"
                T1[k] = T1temp
            print "Updating masks"
            source_mask[T1temp] = 1
            array_mask[T1temp] = k_comb+1 #k_comb cos can be different from k if holes join up - saved as k + 1
            del T1temp
            k+=1

    fsky = (len(source_mask[source_mask == 0]) + 0.) / npix
    print "\nNumber of sources =", nsources
    print "f_sky =", fsky
    write_mask = np.zeros_like(source_mask)
    write_mask[source_mask == 0] = 1
    #hp.write_map('/Users/keir/Documents/s2let_ilc_planck/s2let_ilc_pccs2e_sourcemask.fits',write_mask)


