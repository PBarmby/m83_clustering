import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
import os
import matplotlib.pyplot as plt

ubvi_list = ['336','438','555','814']
ubvi_zps = [23.46, 24.98, 25.81, 24.67]
ubvi_pos_offsets = [(0.0001375,0.0001375), (0.00015,0.000142), (0.00014,0.00021), (0.000137, 0.00022)]
match_tol = 1.0 # arcsec

def go(filt_list = ubvi_list, zps=ubvi_zps, ers = '../raw_data/hlsp_wfc3ers_hst_wfc3_m83_cat_all_v2-corrected.txt'):

    ers_dat = Table.read(ers, format = 'ascii.commented_header')
    
    for i,filt in enumerate(filt_list):
        # generate Sextractor call
        img = 'f{}.fits'.format(filt) 
        catname = 'f{}.cat'.format(filt) 
        sysstr = 'sex {} -CATALOG_NAME {} -MAG_ZEROPOINT {}'.format(img,catname,zps[i])
        # run Sextractor
        os.system(sysstr)

        # extract relevant bits of ERS catalog
        m1 = 'mag05_{}'.format(filt)
        m2 = 'mag3_{}'.format(filt)
        ers_cols = ['ra', 'dec', m1, m1+'_unc', m2, m2+'_unc']
        gooddat = np.logical_or(ers_dat[m1]>0, ers_dat[m2]>0) 
        ers_subcat = ers_dat[ers_cols][gooddat]
        ers_coo = SkyCoord(ra=ers_subcat['ra']*u.degree, dec=ers_subcat['dec']*u.degree) 
       
        # read SE results
        new_cat = Table.read(catname, format = 'ascii')
        new_coo = SkyCoord(ra=(new_cat['col3']+ubvi_pos_offsets[i][0])*u.degree, dec=(new_cat['col4']+ubvi_pos_offsets[i][1])*u.degree) 

        # match results with ERS catalog
        idx, sep2d, sep3d = match_coordinates_sky(new_coo, ers_coo) 
        matched = sep2d < match_tol*u.arcsec

        # compute magnitude/uncert offsets
        dm1 = ers_subcat[m1][idx] - new_cat['col5']
        dm2 = ers_subcat[m2][idx] - new_cat['col6']
        du1 = ers_subcat[m1+'_unc'][idx] - new_cat['col7']
        du2 = ers_subcat[m2+'_unc'][idx] - new_cat['col8']
        
#        fig, ax = plt.subplots()

    #end of loop over images
    return
