import numpy as np
from astropy.table import table
import os


ubvi_list = ['336','438','555','814']
ubvi_zps = [23.46, 24.98, 25.81, 24.67]
ubvi_pos_offsets = [(0.0001375,0.0001375), (0.00015,0.000142), (0.00014,0.00021), (0.000137, 0.00022)]

ers_cols = []

def go(filt_list = ubvi_list, zps=ubvi_zps, ers = '../raw_data/hlsp_wfc3ers_hst_wfc3_m83_cat_all_v2-corrected.txt'):

    ers_dat = Table.read(ers, format = 'ascii.commented_header')
    
    for filt,i in enumerate(filt_list):
        # generate Sextractor call
        img = 'f{}.fits'.format(filt) 
        catname = 'f{}.cat'.format(filt) 
        sysstr = 'sex {} -CATALOG_NAME {} -MAG_ZEROPOINT {}'.format(img,catname,zps[i])
        # run Sextractor
        print sysstr
#        os.system(sysstr)

        # match results with ERS catalog
        ers_cols = ['ra','dec','mag05_{}'.format(filt), 'mag05_{}_unc'.format(filt), 'mag3_{}'.format(filt), 'mag3_{}_unc'.format(filt)]
        ers_subcat = ers_dat[ers_cols]
        
        
        # compute magnitude/uncert offsets

    #end of loop over images
    return
