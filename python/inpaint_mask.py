import numpy as np
import healpy as hp
import math as mh
import copy as cp

def find_holes(mask,nside):
    #TESTING
    '''mask[1e7:] = 1
    mask = hp.reorder(mask,r2n=True)'''
    
    mask_zero = np.where(mask==0)[0] #Indices where mask == 0 (npix_0)
    #print "Number of 'hole' pixels =", len(mask_zero)
    pix_neighbours = hp.get_all_neighbours(nside,mask_zero) #Neighbours to each hole pixel (8 x npix_0)[-1 if no]
    hole_index = np.array([-1]*len(mask_zero)) #Initialising array that holds hole index for each hole pixel
    hole_index[0] = 0 #mask_zero[0] #Creating first hole index - using pixel numbers as index
    k = 1

    for i in xrange(1,len(mask_zero)): #Looping over hole pixels
        print i, '/', len(mask_zero)
        sorted_index = np.searchsorted(mask_zero,pix_neighbours[:,i]) #Where neighbours fit in
        yindex = np.take(np.arange(len(mask_zero)),sorted_index,mode="clip") #Clips indices
        holes_to_be_joined = np.unique(hole_index[yindex][mask_zero[yindex] == pix_neighbours[:,i]]) #Inc. '-1'
        holes_to_be_joined = np.delete(holes_to_be_joined,np.where(holes_to_be_joined == -1)[0]) #Remove '-1'
        
        for j in xrange(1,len(holes_to_be_joined)): #Joining holes to lowest index
            hole_index[hole_index == holes_to_be_joined[j]] = holes_to_be_joined[0]
        if len(holes_to_be_joined) == 0: #If brand new hole
            hole_index[i] = k #mask_zero[i]
            k+=1
        else: #Else join to newly-expanded hole
            hole_index[i] = holes_to_be_joined[0]

    np.save('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask_holes_ring.npy',np.vstack((mask_zero,hole_index)))

    return mask_zero, hole_index

def find_rims(holes,nside):
    hole_indices = np.unique(holes[1])
    nholes = len(hole_indices)
    circpixs = [None]*nholes
    rimpixs = [None]*nholes
    rimindex = [None]*nholes
    
    for hole in xrange(nholes):
        print "Forming rim to hole", hole+1, "/", nholes
        #Extracting hole pixels and rim pixels
        hole_index = hole_indices[hole]
        circpixs[hole] = holes[0,np.where(holes[1] == hole_index)[0]]
        rimpixs[hole] = rims[0,np.where(rims[1] == hole_index)[0]]
        pix_neighbours = hp.get_all_neighbours(nside,circpixs[hole]) #Inc. "-1"
        rimpixs[hole] = np.setdiff1d(np.concatenate(pix_neighbours),holes[0]) #Remove any existing hole pixels
        lim = len(rimpixs[hole])
        while lim < 80: #Iterate until have at least 80 border pixels (N_side = 2048)
            pix_neighbours = hp.get_all_neighbours(nside,rimpixs[hole])
            new_rimpixs = np.setdiff1d(np.concatenate(pix_neighbours),np.concatenate((holes[0],rimpixs[hole]))) # -holes & rim
            rimpixs[hole] = np.concatenate((rimpixs[hole],new_rimpixs))
            lim = len(rimpixs[hole])
        rimindex[hole] = np.array([hole_index]*lim)

        np.save('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask_rims_ring.npy',np.vstack((np.concatenate(rimpixs),np.concatenate(rimindex))))

    return np.concatenate(rimpixs), np.concatenate(rimindex)

def gauss_inpaint(T1,T2,T1_tilde,T2_tilde,sigma12,sigma22):
    diff2 = T2 - T2_tilde
    #sigma12 = np.outer(T1,T2)
    #sigma22 = np.outer(T2,T2)
    #Add noise to diagonal elements
    for i in xrange(len(sigma22)):
        sigma22[i,i] = sigma22[i,i] + 1.e-20
    sigma22_inv = np.linalg.inv(sigma22)
    sigma_constrain = np.dot(sigma12,sigma22_inv)
    T1_constrain = np.dot(sigma_constrain,diff2)
    T1_hat = T1_tilde + T1_constrain

    return T1_hat #diff2,sigma12,sigma22

if __name__ == "__main__":
    nside = 2048
    #mask = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits') #0 where holes
    #y = find_holes(mask,nside)

    bad_map = hp.read_map('/home/keir/s2let_ilc_data/1dot2/s2let_ilc_dir_hypatia_memeff_planck_deconv_tapered_3999_1dot2_25_1_recon.fits') #/home/keir/s2let_ilc_data/masks/planck2015_2_cmb_map_99.fits') #/Users/keir/Documents/s2let_ilc_planck/planck2015_2_cmb_map_51.fits') #/Users/keir/Documents/s2let_ilc_planck/deconv_data/s2let_ilc_dir_hypatia_memeff_planck_deconv_tapered_3999_1dot2_25_1_recon.fits')
    good_map = hp.read_map('/home/keir/s2let_ilc_data/masks/planck2015_2_cmb_map_100.fits') #/Users/keir/Documents/s2let_ilc_planck/planck2015_2_cmb_map_52.fits')
    new_map = cp.deepcopy(bad_map)
    holemap = cp.deepcopy(bad_map)

    #Using query_disc to form holes and rims
    '''k = 0
    circpixs = [None]*6
    rimpixs = [None]*6
    for i in [15,21,30]:
        for j in xrange(2):
            theta = np.random.random(1)*mh.pi #random
            phi = np.random.random(1)*mh.pi*2. #random
            circpixs[k] = hp.query_disc(nside,hp.ang2vec(theta[0],phi[0]),np.radians(i/60.)) #21 arcmin
            rimpixs[k] = np.setxor1d(hp.query_disc(nside,hp.ang2vec(theta[0],phi[0]),np.radians(i/60.)+(0.00015*mh.pi)),circpixs[k],assume_unique=True) #21 arcmin + a bit [0.00015*mh.pi]
            print len(circpixs[k]), len(rimpixs[k])
            k+=1
    nholes = len(circpixs)'''
    
    #Using NILC mask holes
    holes = np.load('/home/keir/s2let_ilc_data/masks/nilc_pr1_builtmask_holes_ring.npy') #/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask_holes_ring.npy') #Pix no,hole index
    rims = np.load('/home/keir/s2let_ilc_data/masks/nilc_pr1_builtmask_rims_ring.npy') #Pix no., hole index
    hole_indices = np.unique(holes[1])
    nholes = len(hole_indices)
    circpixs = [None]*nholes
    rimpixs = [None]*nholes
    rimindex = [None]*nholes

    #Loading random realisations for covariance estimation
    rand_realise = [None]*98
    for i in xrange(98):
        print '\n', i+1
        fname = '/home/keir/s2let_ilc_data/masks/planck2015_2_cmb_map_' + str(i+1) + '.fits' #'/Users/keir/Documents/s2let_ilc_planck/planck2015_2_cmb_map_' + str(i+1) + '.fits'
        rand_realise[i] = hp.read_map(fname,memmap=True)
    
    #Filling in holes
    #for hole in xrange(len(circpixs)): #Looping over holes
    for hole in xrange(nholes):
        print "Filling in hole", hole+1, "/", nholes
        #Extracting hole pixels and rim pixels
        hole_index = hole_indices[hole]
        circpixs[hole] = holes[0,np.where(holes[1] == hole_index)[0]]
        rimpixs[hole] = rims[0,np.where(rims[1] == hole_index)[0]]
        rimpixs[hole] = np.delete(rimpixs[hole],np.where(rimpixs[hole] == -1)[0]) #Remove '-1'
        
        T1 = bad_map[circpixs[hole]]
        T2 = bad_map[rimpixs[hole]]
        T1_tilde = good_map[circpixs[hole]]
        T2_tilde = good_map[rimpixs[hole]]
        sigma12_rand = [None]*98
        sigma22_rand = [None]*98
        for i in xrange(98):
            #print i+1
            T1_rand = rand_realise[i][circpixs[hole]]
            T2_rand = rand_realise[i][rimpixs[hole]]
            sigma12_rand[i] = np.outer(T1_rand,T2_rand)
            sigma22_rand[i] = np.outer(T2_rand,T2_rand)
        sigma12_rand = np.mean(np.array(sigma12_rand),axis=0)
        sigma22_rand = np.mean(np.array(sigma22_rand),axis=0)
        
        newvals = gauss_inpaint(T1,T2,T1_tilde,T2_tilde,sigma12_rand,sigma22_rand)
        new_map[circpixs[hole]] = newvals

    hp.write_map('/home/keir/s2let_ilc_data/masks/s2let_ilc_dir_hypatia_memeff_planck_deconv_tapered_3999_1dot2_25_1_recon_inpaint.fits',new_map)

    resid_map = (new_map - bad_map) / bad_map

    holemask = np.zeros(hp.nside2npix(nside))
    holemask[np.concatenate(circpixs[:])] = 1
    holemask[np.concatenate(rimpixs[:])] = 2
    holemap[np.concatenate(circpixs[:])] = np.nan

    print "Calculating power spectra"
    '''bad_cls = hp.anafast(bad_map,lmax=3999)
    new_cls = hp.anafast(new_map,lmax=3999)
    ell = np.arange(len(bad_cls))
    invtwopi = 1. / (2.*mh.pi)'''
    
    '''circpixs2 = np.union1d(hp.query_disc(nside,hp.ang2vec(.25*mh.pi,0.*mh.pi),.005*mh.pi),hp.query_disc(nside,hp.ang2vec(.25*mh.pi,.005*mh.pi),.005*mh.pi))
    testmask = np.ones(hp.nside2npix(nside))
    testmask[circpixs] = 0
    testmask[circpixs2] = 0

    y = find_holes(testmask,nside)
    holemask = np.zeros(hp.nside2npix(nside)) - 1.
    holemask[y[0]] = y[1]'''
