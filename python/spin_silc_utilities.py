import numpy as np
import healpy as hp

def masked_cl(ellmax,mask_name,mask_fn,outfits_root,final_maps):
    #Masked spectra
    print "Calculating masked C_l"
    mask = hp.read_map(mask_fn) #0 where holes #RING ordering
    f_sky = np.sum(mask) / len(mask)
    print "f_sky =", f_sky
    masked_maps = [final_maps[0]*mask,final_maps[1]*mask] #(Q*mask,U*mask)
    T_map = [np.zeros_like(final_maps[0]),]
    masked_cls = hp.anafast(T_map+masked_maps,lmax=ellmax-1) #(TT,EE,BB,TE,EB,TB)
    pixrecip = np.concatenate((np.ones(2),np.reciprocal(hp.pixwin(hp.get_nside(masked_maps[0]),pol=True)[1][2:ellmax]))) #P pixwin #Not defined for l < 2
    clscode = ['EE','BB','EB']
    clsidx = [1,2,4]
    for i in xrange(len(clscode)):
        cl_outfits = outfits_root + '_' + clscode[i] + 'cls_' + mask_name + '.fits'
        hp.write_cl(cl_outfits,masked_cls[clsidx[i]] * pixrecip * pixrecip / f_sky)

    return 0