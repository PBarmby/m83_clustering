'''Plots'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pylab
import os
import os.path
import argparse
import pandas as pd 
from pandas.tools.plotting import scatter_matrix

cluster_colours = ['y', 'g', 'b', 'r', 'c', 'm', 'k', 'w', 'brown', 'darkgray', 'orange', 'pink','gold', 'lavender', 'salmon', 'g', 'b',
                   'r', 'c', 'm', 'k', 'm', 'w', 'y', 'g', 'b', 'r', 'c', 'm',
                   'k','m','w','y','g','b','r','c','m','k','m','w','y','g','b',
                   'r','c','m','k','m','w','y','g','b','r','c','m','k','m','w',
                   'y','g','b','r','c','m','k','m','w','y','g','b','r','c','m',
                   'k','m','w','y','g','b','r','c','m','k','m','w','y','g','b',
                   'r','c','m','k','m','w','y','g','b','r','c','m','k','m','w',
                   'y','g','b','r','c','m','k','m','w']


def plotting(data_table, path, plots, threshold, survey_objects):
    plot_path = make_directory(path)
    data_, trial = load_data(data_table, 'experiments.txt')

    for i in range(0, len(trial)):

        wave1, wave2, wave3, wave4, grd, c1, c2 = data(data_,
                                                       trial['band1'][i],
                                                       trial['band2'][i],
                                                       trial['band3'][i],
                                                       trial['band4'][i])
        if 'bvb' in plots:
            band_v_band(wave1, wave2, wave3, wave4, grd, trial[i], plot_path)
        if 'w_hist' in plots:
            wave_histogram(wave1, wave2, wave3, wave4, grd, trial[i],
                           plot_path)
        if 'c_hist' in plots:
            colour_histogram(c1, c2, trial[i], plot_path)
        if 'cvc' in plots:
            colour_v_colour(c1, c2, trial[i], plot_path)
        if 'wvc' in plots:
            wave_v_colour(wave1, wave2, wave3, wave4, c1, c2, grd, trial[i],
                          plot_path)

    return


def make_directory(p_path):
    '''Save results of each analysis
            save_path: set path of where you would like results saved'''
    # Create new plots_folder
    pl_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(
              p_path)
    if not os.path.exists(pl_path):
        os.makedirs(pl_path)

    return(pl_path)


def load_data(surveyfile_, experiments):
    '''User upload data file'''

    d_file_name = str(surveyfile_)
    t_file_name = str(experiments)
    data = Table.read(d_file_name, format='ascii.commented_header',
                      guess=False)
    tests = Table.read(t_file_name, format='ascii.commented_header',
                       guess=False)

    return (data, tests)


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


def band_v_band(wave1, wave2, wave3, wave4, greatdata, labels, path):

    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(211)
    ax.scatter(wave1[greatdata], wave2[greatdata], color='k', marker='.')
    ax.set_xlabel(labels['band1'])
    ax.set_ylabel(labels['band2'])
    ax.set_title(labels['band1']+' vs. '+labels['band2'], fontsize=12)

    ax = fig.add_subplot(212)
    ax.scatter(wave3[greatdata], wave4[greatdata], color='r', marker='.')
    ax.set_xlabel(labels['band3'])
    ax.set_ylabel(labels['band4'])
    ax.set_title(labels['band3']+' vs. '+labels['band4'], fontsize=12)

    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'plot_{}_vs_{}-{}_vs_{}.png'.format(labels['band1'],
                                               labels['band2'],
                                               labels['band3'],
                                               labels['band4'])
    path_ = ('{}\\band_band_plts').format(path)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return


def wave_histogram(wave1, wave2, wave3, wave4, greatdata, trial, path):
    fig2 = plt.figure(figsize=(8, 8))

    ax = fig2.add_subplot(2, 2, 1)
    ax.hist(wave1[greatdata], bins=50, range=[min(wave1[greatdata]),
            max(wave1[greatdata])])
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_title(trial['band1']+' Magnitude Frequency',
                 fontsize=12)

    ax = fig2.add_subplot(2, 2, 2)
    ax.hist(wave2[greatdata], bins=50, range=[min(wave2[greatdata]),
            max(wave2[greatdata])])
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_title(trial['band2']+' Magnitude Frequency', fontsize=12)

    ax = fig2.add_subplot(2, 2, 3)
    ax.hist(wave3[greatdata], bins=50, range=[min(wave3[greatdata]),
            max(wave3[greatdata])])
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_title(trial['band3']+' Magnitude Frequency', fontsize=12)

    ax = fig2.add_subplot(2, 2, 4)
    ax.hist(wave4[greatdata], bins=50, range=[min(wave4[greatdata]),
            max(wave4[greatdata])])
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_title(trial['band4']+' Magnitude Frequency', fontsize=12)
    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'w_hist_{}-{}-{}-{}.png'.format(trial['band1'],
                                                trial['band2'],
                                                trial['band3'],
                                                trial['band4'])
    path_ = ('{}\\wave_histogram').format(path)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return


def colour_histogram(colour1, colour2, trial, path):

    fig3 = plt.figure(figsize=(8, 8))

    ax = fig3.add_subplot(211)
    ax.hist(colour1, bins=50, range=[min(colour1), max(colour1)])
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_title(trial['band1'] + ' - ' + trial['band2'] + ' Frequency')

    ax = fig3.add_subplot(212)
    ax.hist(colour2, bins=50, range=[min(colour2), max(colour2)])
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_title(trial['band3'] + ' - ' + trial['band4'] + ' Frequency')
    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'c_hist_{}-{}_vs_{}-{}.png'.format(trial['band1'],
                                                   trial['band2'],
                                                   trial['band3'],
                                                   trial['band4'])

    path_ = ('{}\\colour_histogram').format(path)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return


def colour_v_colour(colour1, colour2, trial, path):

    matrix_file = 'matrix_{}-{}vs{}-{}.png'.format(trial['band1'],
                                                   trial['band2'],
                                                   trial['band3'],
                                                   trial['band4'])

    df = pd.DataFrame(np.vstack([colour1, colour2]).T,
                      columns=[trial['band1'] + '-' + trial['band2'],
                               trial['band3'] + '-' + trial['band4']])
    scatter_matrix(df, figsize=(8, 8), alpha=0.5,  diagonal='kde')
    path_ = ('{}\\color_color').format(path)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, matrix_file))
    plt.close()
    return


def wave_v_colour(wave1, wave2, wave3, wave4, colour1, colour2, greatdata,
                  trial, path):
    fig5 = plt.figure(figsize=(8, 8))
    ax = fig5.add_subplot(2, 2, 1)
    ax.scatter(colour1, wave1[greatdata], color='b', marker='.')
    ax.set_xlabel(trial['band1']+' - '+trial['band2'])
    ax.set_ylabel(trial['band1'])
    plt.gca().invert_yaxis()

    ax = fig5.add_subplot(2, 2, 2)
    ax.scatter(colour1, wave2[greatdata], color='b', marker='.')
    ax.set_xlabel(trial['band1']+' - '+trial['band2'])
    ax.set_ylabel(trial['band2'])
    plt.gca().invert_yaxis()

    ax = fig5.add_subplot(2, 2, 3)
    ax.scatter(colour2, wave3[greatdata], color='b', marker='.')
    ax.set_xlabel(trial['band3']+' - '+trial['band4'])
    ax.set_ylabel(trial['band3'])
    plt.gca().invert_yaxis()

    ax = fig5.add_subplot(2, 2, 4)
    ax.scatter(colour2, wave4[greatdata], color='b', marker='.')
    ax.set_xlabel(trial['band3']+' - '+trial['band4'])
    ax.set_ylabel(trial['band4'])
    plt.gca().invert_yaxis()

    '''Display interactive figure if # removed, if not, figures saved'''
    plt.show
    file_name = 'wave_colour_plot_{}-{}vs{}-{}.png'.format(trial['band1'],
                                                           trial['band2'],
                                                           trial['band3'],
                                                           trial['band4'])

    path_ = ('{}\\wave_color_plts').format(path)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return


def inputs():

    inputs = argparse.ArgumentParser(description='Plotting functions')
    inputs.add_argument("data_file", help="Specify data file")
    inputs.add_argument("-pp", "--plot_path", help="Save directory",
                        default="results")
    inputs.add_argument("-p", "--plots", help="Select plots to make",
                        choices=['bvb', 'w_hist', 'c_hist', 'cvc',
                                 'wvc'], nargs='*')
    # Data trimming
    inputs.add_argument("-r", "--ratio", help="Set noise to signal ratio",
                        default=0.1)
    inputs.add_argument("-dp", "--data_points", help="Set objects from survey",
                        default=10000)

    criteria = inputs.parse_args()
    plotting(criteria.data_file, criteria.plot_path, criteria.plots,
             criteria.ratio, criteria.data_points)
    return

if __name__ == "__main__":
    inputs()
