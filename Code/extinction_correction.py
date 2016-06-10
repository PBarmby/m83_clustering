import numpy as np 
from astropy.table import Table
from matplotlib import pyplot as plt

def extinction_correction(path_):
    data_ = Table.read('data.txt', format='ascii.commented_header', guess=False)
    trial = Table.read('experiments.txt', format='ascii.commented_header', guess=False)
    mw_ex = 0.58
    for k in range(0, len(trial)): 
        wave1, wave2, wave3, wave4, grd, c1, c2 = data(data_,
                                                       trial['band1'][k],
                                                       trial['band2'][k],
                                                       trial['band3'][k],
                                                       trial['band4'][k])
        correction = wave3[grd] - wave4[grd] - mw_ex*(wave1[grd] - wave2[grd])
        plt.plot(c1*correction, c2*correction, '.')
    return


def data(data, band1, band2, band3, band4):
    '''Select data for analysis'''
    ratio = 0.25
    data = data#[10000:30000]
    # Colour 1
    wave1 = data[band1]
    wave1_unc = data[band1+'_unc']
    wave2 = data[band2]
    wave2_unc = data[band2+'_unc']
    # Colour 2
    wave3 = data[band3]
    wave3_unc = data[band3+'_unc']
    wave4 = data[band4]
    wave4_unc = data[band4+'_unc']
    # Change parameters to match data_file
    # Remove data pieces with no value

    wave1_trim = np.logical_and(wave1 != -99, wave1_unc != -99)
    wave2_trim = np.logical_and(wave2 != -99, wave2_unc != -99)
    wave3_trim = np.logical_and(wave3 != -99, wave3_unc != -99)
    wave4_trim = np.logical_and(wave4 != -99, wave4_unc != -99)

    colour1_ratio = np.logical_and(wave1_unc < ratio,
                                   wave2_unc < ratio)
    colour2_ratio = np.logical_and(wave3_unc < ratio,
                                   wave4_unc < ratio)

    gooddata1 = np.logical_and(np.logical_and(wave1_trim, wave2_trim),
                               np.logical_and(wave3_trim, wave4_trim))

    # Remove data above certain magnitude
    gooddata2 = np.logical_and(colour1_ratio, colour2_ratio)

    # Only data that match criteria for both colours
    greatdata = np.logical_and(gooddata1, gooddata2)

    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]
    
    return(wave1, wave2, wave3, wave4, greatdata, colour1, colour2)