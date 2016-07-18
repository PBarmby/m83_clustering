'''Test labels'''
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
import pylab
import os
import os.path

colours = ['b', 'g', 'r', 'y', 'k', 'm']

def plot():
    labels= ['HII', 'PN', 'RSG', 'SNR', 'starcluster', 'wr']
    plots = Table.read('plots.txt', format='ascii.commented_header',
                       guess=False)
    label_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\cat_data\\clustering_labels\\'
    
    for p in range(0, len(plots)):
        trial = plots[p]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sum_ = 0
        col = '{}-{}_vs_{}-{}'.format(trial['band1'], trial['band2'],
                                      trial['band3'], trial['band4'])
        print col
        save_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\results\\labels\\{}'.format(col)
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        for i in range(0, len(labels)):
            label_file = label_path + labels[i] + '_labels.txt'
            label_ = Table.read(label_file, format='ascii.commented_header',
                                guess=False)
            w1, w2, w3, w4, c1, c2, gooddata = data(label_, trial)

            ax.scatter(c1, c2, marker='o', c=colours[i])
            ax.set_title(labels[i] + ' Colour-Colour Plot')
            ax.set_xlabel(trial['band1'] + ' - ' + trial['band2'])
            ax.set_ylabel(trial['band3'] + ' - ' + trial['band4'])
            print labels[i] + ' objects:' + str(len(c1))
            plot_file = '{}-plot.png'.format(labels[i])
            pylab.savefig(os.path.join(save_path, plot_file))
            sum_ = sum_ + len(c1)
        print 'Total_obj: ' + str(sum_) + '\n'

    return


def data(label, exp):
    ratio = 0.2

    wave1 = label[exp['band1']]
    wave1_unc = label[exp['band1']+'_unc']
    wave2 = label[exp['band2']]
    wave2_unc = label[exp['band2']+'_unc']
    # Colour 2
    wave3 = label[exp['band3']]
    wave3_unc = label[exp['band3']+'_unc']
    wave4 = label[exp['band4']]
    wave4_unc = label[exp['band4']+'_unc']

    # Change parameters to match data_file
    # Remove data pieces with no value
    wave1_trim = np.logical_and(np.logical_and(wave1 != -99, wave1_unc != -99),
                                wave1_unc < ratio)
    wave2_trim = np.logical_and(np.logical_and(wave2 != -99, wave2_unc != -99),
                                wave2_unc < ratio)

    wave3_trim = np.logical_and(np.logical_and(wave3 != -99, wave3_unc != -99),
                                wave3_unc < ratio)
    wave4_trim = np.logical_and(np.logical_and(wave4 != -99, wave4_unc != -99),
                                wave4_unc < ratio)

    colour1_trim = np.logical_and(wave1_trim, wave2_trim)
    colour2_trim = np.logical_and(wave3_trim, wave4_trim)

    final_data = np.logical_and(colour1_trim, colour2_trim)

    # Make colours
    colour1 = wave1[final_data] - wave2[final_data]
    colour2 = wave3[final_data] - wave4[final_data]

    return(wave1, wave2, wave3, wave4, colour1, colour2, final_data)
