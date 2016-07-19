'''Test labels'''
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
import pylab
import os
import os.path

colours = ['b', 'g', 'y', 'r', 'c', 'm']

def plot():
    labels= ['starcluster', 'HII', 'PN', 'RSG', 'SNR', 'WR']
    plots = Table.read('plots.txt', format='ascii.commented_header',
                       guess=False)
    label_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\cat_data\\clustering_labels\\'
    survey = Table.read('data_v3.txt', format='ascii.commented_header',
                      guess=False)
    for p in range(0, len(plots)):
        trial = plots[p]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        sum_ = 0
        col = '{}-{}_vs_{}-{}'.format(trial['band1'], trial['band2'],
                                      trial['band3'], trial['band4'])
        print col
        save_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\results\\labels'
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        for i in range(0, len(labels)):
            label_file = label_path + labels[i] + '_labels.txt'
            label_ = Table.read(label_file, format='ascii.commented_header',
                                guess=False)
            c1, c2, gooddata, c1_s, c2_s = data(label_, trial, survey)
            if i == 0:
                ax.scatter(c1_s, c2_s, c='k', s=2, zorder=i)

            ax.scatter(c1, c2, marker='o', c=colours[i], s=13, zorder=i+1,
                       edgecolor='None', label=labels[i])
            print labels[i] + ' objects:' + str(len(c1))
            sum_ = sum_ + len(c1)
        ax.legend(loc='upper left', fontsize=8)
        ax.set_title(col)
        ax.set_xlabel(trial['band1'] + ' - ' + trial['band2'])
        ax.set_ylabel(trial['band3'] + ' - ' + trial['band4'])
            
        plot_file = '{}-plot.png'.format(col)
        pylab.savefig(os.path.join(save_path, plot_file))
        plt.close()
        print 'Total_obj: ' + str(sum_) + '\n'

    return


def data(label, exp, survey_data):
    ratio = 0.2

    wave1 = label[exp['band1']]
    wave1_data = survey_data[exp['band1']]
    wave1_unc = label[exp['band1']+'_unc']
    wave1_data_unc = survey_data[exp['band1']+'_unc']
    wave2 = label[exp['band2']]
    wave2_data = survey_data[exp['band2']]
    wave2_unc = label[exp['band2']+'_unc']
    wave2_data_unc = survey_data[exp['band2']+'_unc']
    # Colour 2
    wave3 = label[exp['band3']]
    wave3_data = survey_data[exp['band3']]
    wave3_unc = label[exp['band3']+'_unc']
    wave3_data_unc = survey_data[exp['band3']+'_unc']
    wave4 = label[exp['band4']]
    wave4_data = survey_data[exp['band4']]
    wave4_unc = label[exp['band4']+'_unc']
    wave4_data_unc = survey_data[exp['band4']+'_unc']

    # Change parameters to match data_file
    # Remove data pieces with no value
    wave1_trim = np.logical_and(np.logical_and(wave1 != -99, wave1_unc != -99),
                                wave1_unc < ratio)
    wave1_data_trim = np.logical_and(np.logical_and(wave1_data != -99,
                                                    wave1_data_unc != -99),
                                     wave1_data_unc < ratio)
    wave2_trim = np.logical_and(np.logical_and(wave2 != -99, wave2_unc != -99),
                                wave2_unc < ratio)
    wave2_data_trim = np.logical_and(np.logical_and(wave2_data != -99,
                                                    wave2_data_unc != -99),
                                     wave2_data_unc < ratio)

    wave3_trim = np.logical_and(np.logical_and(wave3 != -99, wave3_unc != -99),
                                wave3_unc < ratio)
    wave3_data_trim = np.logical_and(np.logical_and(wave3_data != -99,
                                                    wave3_data_unc != -99),
                                     wave3_data_unc < ratio)
    wave4_trim = np.logical_and(np.logical_and(wave4 != -99, wave4_unc != -99),
                                wave4_unc < ratio)
    wave4_data_trim = np.logical_and(np.logical_and(wave4_data != -99,
                                                    wave4_data_unc != -99),
                                     wave4_data_unc < ratio)

    colour1_trim = np.logical_and(wave1_trim, wave2_trim)
    colour1_data_trim = np.logical_and(wave1_data_trim, wave2_data_trim)
    colour2_trim = np.logical_and(wave3_trim, wave4_trim)
    colour2_data_trim = np.logical_and(wave3_data_trim, wave4_data_trim)

    final_data = np.logical_and(colour1_trim, colour2_trim)
    survey_data = np.logical_and(colour1_data_trim, colour2_data_trim)

    # Make colours
    colour1 = wave1[final_data] - wave2[final_data]
    colour1_survey = wave1_data[survey_data] - wave2_data[survey_data]
    colour2 = wave3[final_data] - wave4[final_data]
    colour2_survey = wave3_data[survey_data] - wave4_data[survey_data]

    return(colour1, colour2, final_data, colour1_survey, colour2_survey)
