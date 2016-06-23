import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from scipy.stats import binned_statistic
import os
import matplotlib.pyplot as plt

ubvi_list = ['336','438','555','814']
ubvi_zps = [23.46, 24.98, 25.81, 24.67]
ubvi_pos_offsets = [(0.0001375,0.0001375), (0.00015,0.000142), (0.00014,0.00021), (0.000137, 0.00022)]
match_tol = 1.0 # arcsec

def go(filt_list = ubvi_list, zps=ubvi_zps, ers = '../raw_data/hlsp_wfc3ers_hst_wfc3_m83_cat_all_v2-corrected.txt', rerun_se=False):

    ers_dat = Table.read(ers, format = 'ascii.commented_header')
    
    for i,filt in enumerate(filt_list):
        # generate Sextractor call
        img = 'f{}.fits'.format(filt) 
        catname = 'f{}.cat'.format(filt) 
        if rerun_se or not os.path.isfile(catname):
        # run Sextractor
            sysstr = 'sex {} -CATALOG_NAME {} -MAG_ZEROPOINT {}'.format(img,catname,zps[i])
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
        print('{}: matched {} of {} objects'.format(catname, matched.sum(), gooddat.sum()))

#        # compute magnitude/uncert offsets
        m1_old = ers_subcat[m1][idx][matched]
        m2_old = ers_subcat[m2][idx][matched]
        dm1 = (ers_subcat[m1][idx] - new_cat['col5'])[matched]
        dm2 = (ers_subcat[m2][idx] - new_cat['col6'])[matched]
        u1 = ers_subcat[m1+'_unc'][idx][matched]
        u2 = ers_subcat[m2+'_unc'][idx][matched]

        # compute some statistics
        bin_edges = np.arange(18,30,0.5)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2
        med_dm1, junk1, junk2 = binned_statistic(m1_old, dm1, 'median', bin_edges)
        med_dm2, junk1, junk2 = binned_statistic(m2_old, dm2, 'median', bin_edges)
                
        # make some plots
        fig, ax = plt.subplots(2,1)
        ax[0].plot(m1_old, new_cat['col5'][matched], marker='.', ms=5, label= '0.5pix')
        ax[0].plot(m2_old, new_cat['col6'][matched], marker='.', ms=5,label= '3pix')        
        ax[0].plot(bin_centers, med_dm1, marker='None', ls='solid')
        ax[0].plot(bin_centers, med_dm2, marker='None', ls='solid')
        ax[0].set_xlim(18,30)
        ax[0].set_ylim(18,30)
        ax[0].grid()
        ax[0].set_xlabel('Catalog magnitude')
        ax[0].set_ylabel('New mag')
        ax[0].legend(loc='upper left',numpoints=1)
        ax[0].set_title(filt)
        ax[1].plot(u1, new_cat['col7'][matched], marker='.', ms=5,label= '0.5pix')
        ax[1].plot(u2, new_cat['col8'][matched], marker='.', ms=5,label= '3pix')        
        ax[1].set_xlim(-0.1,10)
        ax[1].set_ylim(-0.05,0.5)
        ax[1].set_xlabel('Catalog magnitude unc')
        ax[1].set_ylabel('New mag unc')
        fig.tight_layout()
        figname = 'comp_{}.png'.format(filt)
        plt.savefig(figname)
        
    #end of loop over images
    return
