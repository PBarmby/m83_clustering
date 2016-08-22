'''Plot model colours from FSPS'''
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d as p3
from astropy.table import Table
w1 = 0  # colour 1
w2 = 1  # colour 1
w3 = 1  # colour 2
w4 = 2  # colour 2
w5 = 1  # colour 3
w6 = 3  # colour 3

def make_plots(file_name, band1, band2, band3, band4):
    load_file = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\model_colours\\{}'.format(file_name)
    bands = np.array([band1, band2, band3, band4])
    survey_file = Table.read('data_v3.txt', format='ascii.commented_header',
                             guess=False)
    wave1s, wave2s, wave3s, wave4s, final_data = data(survey_file, band1,
                                                      band2, band3, band4)
    sur_filt = np.vstack([wave1s[final_data], wave2s[final_data],
                                wave3s[final_data], wave4s[final_data]])
    print len(sur_filt[0])

    model = Table.read(load_file, format='ascii.commented_header', guess=False)
    mod_filt = np.vstack([model[band1], model[band2], model[band3],
                               model[band4]])
    print mod_filt[0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(sur_filt[0] - sur_filt[1], sur_filt[2] - sur_filt[3],
               marker='.', s=4)
    ax.plot(mod_filt[0] - mod_filt[1], mod_filt[2] - mod_filt[3])
    ax.set_xlabel(band1 + ' - ' + band2)
    ax.set_ylabel(band3 + ' - ' + band4)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(sur_filt[w1] - sur_filt[w2],
                sur_filt[w3] - sur_filt[w4],
                sur_filt[w5] - sur_filt[w6], marker='.', s=4)
    ax2.plot(mod_filt[w1] - mod_filt[w2],
             mod_filt[w3] - mod_filt[w4],
             mod_filt[w5] - mod_filt[w6])
    ax2.set_xlabel(bands[w1] + ' - ' + bands[w2])
    ax2.set_ylabel(bands[w3] + ' - ' + bands[w4])
    ax2.set_zlabel(bands[w5] + ' - ' + bands[w6])

    return


def data(data_table, band1, band2, band3, band4):
    ratio = 0.2
    wave1 = data_table[band1]
    wave1_unc = data_table[band1+'_unc']
    wave2 = data_table[band2]
    wave2_unc = data_table[band2+'_unc']
    wave3 = data_table[band3]
    wave3_unc = data_table[band3+'_unc']
    wave4 = data_table[band4]
    wave4_unc = data_table[band4+'_unc']
    
    wave1_trim = np.logical_and(np.logical_and(wave1 != -99, wave1_unc !=-99),
                                wave1_unc < ratio)
    wave2_trim = np.logical_and(np.logical_and(wave2 != -99, wave2_unc !=-99),
                                wave2_unc < ratio)
    
    wave3_trim = np.logical_and(np.logical_and(wave3 != -99, wave3_unc !=-99),
                                wave3_unc < ratio)
    wave4_trim = np.logical_and(np.logical_and(wave4 != -99, wave4_unc !=-99),
                                wave4_unc < ratio)

    colour1_trim = np.logical_and(wave1_trim, wave2_trim)
    colour2_trim = np.logical_and(wave3_trim, wave4_trim)
    final_data = np.logical_and(colour1_trim, colour2_trim)

    return(wave1, wave2, wave3, wave4, final_data)