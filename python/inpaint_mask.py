import numpy as np
import healpy as hp
import math as mh
import copy as cp
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
    '''hole_sizes = np.load(holes_dir + 'nilc_pr1_builtmask_holes_ring_sizes.npy') #Hole index, hole size
    hole_indices = hole_sizes[0,np.argsort(hole_sizes[1])[::-1]] #Descending order of size'''
    nholes = len(hole_indices)
    #nholes = 100 #10 largest holes
    circpixs = [None]*nholes
    rimpixs = [None]*nholes
    rimindex = [None]*nholes
    
    for hole in xrange(nholes):
        print "Forming rim to hole", hole+1, "/", nholes
        #Extracting hole pixels
        hole_index = hole_indices[hole]
        circpixs[hole] = holes[0,holes[1] == hole_index]
        pix_neighbours = hp.get_all_neighbours(nside,circpixs[hole]) #Inc. "-1"
        rimpixs[hole] = np.setdiff1d(np.concatenate(pix_neighbours),np.concatenate((holes[0],np.array([-1])))) #Remove any existing hole pixels & "-1"
        #lim = len(rimpixs[hole])
        for iter in xrange(3): #Ensure a minimum rim thickness of 3 pixels
            pix_neighbours = hp.get_all_neighbours(nside,rimpixs[hole])
            new_rimpixs = np.setdiff1d(np.concatenate(pix_neighbours),np.concatenate((holes[0],np.array([-1]),rimpixs[hole]))) # - holes & rim & "-1"
            rimpixs[hole] = np.concatenate((rimpixs[hole],new_rimpixs))
            lim = len(rimpixs[hole])
        while lim < 80: #Iterate until have at least 80 border pixels (N_side = 2048)
            pix_neighbours = hp.get_all_neighbours(nside,rimpixs[hole])
            new_rimpixs = np.setdiff1d(np.concatenate(pix_neighbours),np.concatenate((holes[0],np.array([-1]),rimpixs[hole]))) # - holes & rim & "-1"
            rimpixs[hole] = np.concatenate((rimpixs[hole],new_rimpixs))
            lim = len(rimpixs[hole])
        rimindex[hole] = np.array([hole_index]*lim)

    np.save('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask_rims_ring_enlarged3.npy',np.vstack((np.concatenate(rimpixs),np.concatenate(rimindex))))

    return np.concatenate(rimpixs), np.concatenate(rimindex)

def gauss_inpaint(T2,T1_tilde,T2_tilde,sigma12,sigma22):
    diff2 = T2 - T2_tilde
    #Add noise to diagonal elements
    for i in xrange(len(sigma22)):
        sigma22[i,i] = sigma22[i,i] + 1.e-20
    sigma22_inv = np.linalg.inv(sigma22)
    sigma_constrain = np.dot(sigma12,sigma22_inv)
    T1_constrain = np.dot(sigma_constrain,diff2)
    T1_hat = T1_tilde + T1_constrain

    return T1_hat

def inpaint_para(hole_index):
    print "Filling in hole for hole_index =", hole_index
    circpixs = holes[0,holes[1] == hole_index]
    rimpixs = rims[0,rims[1] == hole_index]
    #rimpixs = np.delete(rimpixs,np.where(rimpixs == -1)[0]) #Remove '-1'
    print "Size of hole =", len(circpixs), "Size of rim =", len(rimpixs)
    
    T2 = bad_map[rimpixs]
    T1_tilde = good_map[circpixs]
    T2_tilde = good_map[rimpixs]
    sigma12_rand = [None]*nrand
    sigma22_rand = [None]*nrand
    for i in xrange(nrand):
        T1_rand = rand_realise[i][circpixs]
        T2_rand = rand_realise[i][rimpixs]
        sigma12_rand[i] = np.outer(T1_rand,T2_rand)
        sigma22_rand[i] = np.outer(T2_rand,T2_rand)
    sigma12_rand = np.mean(np.array(sigma12_rand),axis=0)
    sigma22_rand = np.mean(np.array(sigma22_rand),axis=0)
    
    newvals = gauss_inpaint(T2,T1_tilde,T2_tilde,sigma12_rand,sigma22_rand)
    print "Finished filling in hole for hole_index =", hole_index

    return np.vstack((circpixs,newvals))

if __name__ == "__main__":
    #mask = hp.read_map('/Users/keir/Documents/s2let_ilc_planck/nilc_pr1_builtmask.fits') #0 where holes

    #Set directory structure
    comp = 1
    
    if comp == 0: #Keir's iMac
        nprocess = 4
        bad_dir = '/Users/keir/Documents/s2let_ilc_planck/hybrid_data/'
        #bad_dir = '/Users/keir/Documents/s2let_ilc_planck/'
        good_dir = '/Users/keir/Documents/planck2015_2_cmb_realisations/'
        holes_dir = '/Users/keir/Documents/s2let_ilc_planck/holes/'
    elif comp == 1: #Hypatia
        nprocess = 24
        #bad_dir = '/home/keir/s2let_ilc_data/hybrid_data/'
        bad_dir = '/home/keir/s2let_ilc_data/masks/'
        good_dir = '/home/keir/s2let_ilc_data/masks/'
        holes_dir = good_dir
    
    nrand = 998
    ncmb = 2 #nrand + ncmb = 1000

    #bad_map = hp.read_map(bad_dir + 's2let_ilc_planck_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon.fits')
    bad_map = hp.read_map(good_dir + 'planck2015_2_cmb_map_1.fits')
    good_map = hp.read_map(good_dir + 'planck2015_2_cmb_map_2.fits')
    #outfits = 's2let_ilc_planck_deconv_tapered_thresh_lmax3600_3600_hybridC_0_1_recon_inpaint.fits'
    outfits = 'planck2015_2_cmb_map_1_inpaint.fits'
    nside = hp.get_nside(bad_map)
    new_map = cp.deepcopy(bad_map)
    hole_map = cp.deepcopy(bad_map)

    #Using NILC mask holes and rims
    holes = np.load(holes_dir + 'nilc_pr1_builtmask_holes_ring.npy') #Pix no, hole index
    #holes = holes[:,np.where(holes[1] < 200)[0]] #Limit no. holes for testing
    rims = np.load(holes_dir + 'nilc_pr1_builtmask_rims_ring.npy') #Pix no., hole index
    #rims = rims[:,np.where(rims[1] < 200)[0]]

    '''hole_sizes = np.load(holes_dir + 'nilc_pr1_builtmask_holes_ring_sizes.npy') #Hole index, hole size
    hole_indices = hole_sizes[0,hole_sizes[1]<900]'''
    hole_indices = np.unique(holes[1]) #Sorted and unique
    print 'No. holes =', len(hole_indices)

    #Using query_disc to form holes and rims
    '''k = 0
    ncirc = 2
    circpixs_disc = [None]*ncirc
    circpixs_disc_indices = [None]*ncirc
    rimpixs_disc = [None]*ncirc
    rimpixs_disc_indices = [None]*ncirc
    for i in [21]: #[15,21,30]:
        for j in xrange(ncirc):
            #theta = np.random.random(1)*mh.pi #random
            #phi = np.random.random(1)*mh.pi*2. #random
            #theta = 0.5*mh.pi
            #phi = (k/float(ncirc))*mh.pi
            theta = mh.pi / 6.
            phi = k*theta
            circpixs_disc[k] = hp.query_disc(nside,hp.ang2vec(theta,phi),np.radians(i/60.)) #21 arcmin
            circpixs_disc_indices[k] = np.array([k]*len(circpixs_disc[k]))
            rimpixs_disc[k] = np.setxor1d(hp.query_disc(nside,hp.ang2vec(theta,phi),np.radians(i/60.)+(0.0006*mh.pi)),circpixs_disc[k],assume_unique=True) #21 arcmin + a bit [0.00015*mh.pi]
            rimpixs_disc_indices[k] = np.array([k]*len(rimpixs_disc[k]))
            print len(circpixs_disc[k]), len(rimpixs_disc[k])
            k+=1
    holes_pixs = np.concatenate(circpixs_disc)
    holes_pixs_indices = np.concatenate(circpixs_disc_indices)
    holes = np.vstack((holes_pixs,holes_pixs_indices))
    rims_pixs = np.concatenate(rimpixs_disc)
    rims_pixs_indices = np.concatenate(rimpixs_disc_indices)
    rims = np.vstack((rims_pixs,rims_pixs_indices))
    hole_indices = np.unique(holes[1]) #Sorted and unique'''

    #Loading random realisations for covariance estimation
    rand_realise = [None]*nrand
    for i in xrange(nrand):
        print i + ncmb + 1
        rand_realise[i] = np.load(good_dir + 'planck2015_2_cmb_map_' + str(i+ncmb+1) + '.npy',mmap_mode='r')
    
    #Filling in holes
    pool = mg.Pool(nprocess)
    newvals_list = pool.map(inpaint_para,hole_indices) #Ordered by hole index
    pool.close()
    pool.join()
    newvals_complete = np.concatenate(newvals_list,axis=-1)
    new_map[newvals_complete[0].astype(int)] = newvals_complete[1]
    hp.write_map(bad_dir + outfits,new_map)
    #hp.write_map(good_dir + 'planck2015_2_cmb_map_2_inpaint_test1thickcirc_rand999.fits',new_map)

    resid_map = new_map - bad_map
    hole_map[holes[0]] = np.nan
    hole_mask = np.zeros(hp.nside2npix(nside))
    hole_mask[holes[0]] = 1
    hole_mask[rims[0]] = 2


