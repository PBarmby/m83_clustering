import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from scipy.stats import binned_statistic
import os
import matplotlib.pyplot as plt

ubvi_list = ['336','438','555','814']
ubvi_zps = [23.46, 24.98, 25.81, 24.67]
ubvi_pos_offsets = [(0.0001375,0.0001375), (0.000154,0.000154), (0.000148,0.000215), (0.000137, 0.00022)]
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

        # match SE results with ERS catalog
        # NB: matching results cross-checked vs TOPCAT: same
        idx, sep2d, sep3d = match_coordinates_sky(new_coo, ers_coo) 
        matched = sep2d < match_tol*u.arcsec
        print('{}: matched {} of {} objects'.format(catname, matched.sum(), gooddat.sum()))

#        # pull out data for comparison
        m1_old = ers_subcat[m1][idx][matched]
        m2_old = ers_subcat[m2][idx][matched]
        m1_new = new_cat['col5'][matched]
        m2_new = new_cat['col6'][matched]
        u1_old = ers_subcat[m1+'_unc'][idx][matched]
        u2_old = ers_subcat[m2+'_unc'][idx][matched]
        u1_new = new_cat['col7'][matched]
        u2_new = new_cat['col8'][matched]

        # compute some binned statistics
        bin_edges = np.arange(18,30,0.5)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_c1 = bin_edges[1:] - bin_width/2
        med_m1, junk1, junk2 = binned_statistic(m1_old, m1_new, 'median', bin_edges)
        med_m2, junk1, junk2 = binned_statistic(m2_old, m2_new, 'median', bin_edges)
        bin_edges = np.arange(0,3,0.25)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_c2 = bin_edges[1:] - bin_width/2
        med_u1, junk1, junk2 = binned_statistic(u1_old, u1_new, 'median', bin_edges)
        med_u2, junk1, junk2 = binned_statistic(u2_old, u2_new, 'median', bin_edges)
                
        # make some plots
        fig, ax = plt.subplots(2,1)
        ax[0].plot(m2_old, m2_new, marker='.', ms=1,label= '3pix')        
        ax[0].plot(m1_old, m1_new, marker='.', ms=1, label= '0.5pix')
        ax[0].plot(bin_c1, med_m1, marker='None', ls='solid',color='lightgreen')
        ax[0].plot(bin_c1, med_m2, marker='None', ls='solid',color='cyan')
        ax[0].set_xlim(18,30)
        ax[0].set_ylim(18,30)
        ax[0].grid(ls='dotted')
        ax[0].set_xlabel('ERS catalog magnitude')
        ax[0].set_ylabel('New mag')
        ax[0].legend(loc='upper left',numpoints=1,markerscale=5)
        ax[0].set_title('F{}W comparison'.format(filt))
        ax[1].plot(u2_old, u2_new, marker='.', ms=1,label= '3pix')
        ax[1].plot(u1_old, u1_new, marker='.', ms=1,label= '0.5pix')
        ax[1].plot(bin_c2, med_u1, marker='None', ls='solid',color='lightgreen')
        ax[1].plot(bin_c2, med_u2, marker='None', ls='solid',color='cyan')
        ax[1].set_xlim(-0.1,4)
        ax[1].set_ylim(-0.02,0.4)
        ax[1].grid(ls='dotted')
        ax[1].set_xlabel('ERS catalog magnitude unc')
        ax[1].set_ylabel('New mag unc')
        fig.tight_layout()
        figname = 'comp_{}.png'.format(filt)
        plt.savefig(figname)
        plt.close(fig)
        
    #end of loop over images
    return
