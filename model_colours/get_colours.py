# extract predicted colours in WFC3 bands from FSPS population models
import numpy as np
import fsps

broad_filts = [ 'wfc3_uvis_f225w', 'wfc3_uvis_f336w', 'wfc3_uvis_f438w', 'wfc3_uvis_f555w', 'wfc3_uvis_f814w']
narrow_filts = [ 'wfc3_uvis_f373n', 'wfc3_uvis_f487n','wfc3_uvis_f502n', 'wfc3_uvis_f657n', 'wfc3_uvis_f673n']
our_filters = broad_filts + narrow_filts
our_filters.sort()

# values of 'zmet' parameter for default isochrones (MIST): key is [Fe/H] value
met_vals = {'-2': 2, '-1': 6, '0':10, '+0.5': 12}
#met_vals = {'0':10}

def extract_cols(filters=our_filters, outfile='mist_ssp_feh'):
    # construct header for output file
    file_header = 'logt_yr '
    for f in filters:
        file_header = file_header + f[-5:] + ' '
    # generate SSP - this can take a while
    sp = fsps.StellarPopulation(compute_vega_mags=True, add_neb_emission=True, cloudy_dust=True, sfh=0,  zmet=10)
    ages = sp.log_age
    # modify SSP with different metallicity values, extract mags, write to file
    for met in met_vals.keys():
        output_file = outfile+met
        mags = sp.get_mags(zmet = met_vals[met], bands = filters)
        alldat = np.vstack((ages.T,mags.T)).T # transposing makes the stacking work correctly
        np.savetxt(output_file, alldat, fmt = '%+7.4f', header=file_header)
    # end of loop over metallicity values
    return

