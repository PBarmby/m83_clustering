'''Makes histograms of magnitude distribution and distribution with uncertainty
    limit
        - data: data_v3.txt
 Run from Spyder terminal
 Set limit variable in each function'''
import numpy as np
import os, os.path
import pylab as pb
from matplotlib import pyplot as plt
from astropy.table import Table


def mag_hist_cut(path): 
    data_ = load_data('data_v3.txt')
    plot_path = make_directory(path)
    limit = 0.4  # CHANGE
    filters = data_.colnames
    for i in range(11, len(filters), 2):
        # Load wave and uncertainty
        band = filters[i]
        wave = data_[band]
        wave_unc = data_[band + '_unc']
        # Remove bad data
        bad_data = np.logical_and(wave != -99, wave_unc != -99)
        wave_trim = wave[bad_data]
        wave_unc_trim = wave_unc[bad_data]
        # Set limited data
        wave_limit = wave_trim[wave_unc_trim < limit]
        wave_unc_limit = wave_unc_trim[wave_unc_trim < limit]
        # Make histograms of wave and wave_limit
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        ax.hist(wave_trim, bins=100, range=[min(wave_trim), max(wave_trim)],
                color='b')
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        ax.hist(wave_limit, bins=100, range=[min(wave_limit), max(wave_limit)],
                color='r')
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        ax.set_title(band + ' Magnitude Frequency', fontsize=12)
        ax.set_xlabel(band + ' Magnitude')
        ax.set_ylabel('Frequency')
        '''Display interactive figure if # removed, if not, figures saved'''
        file_name = 'mag_hist_{}-unc_limit-{}.png'.format(band, limit)
        path_ = ('{}\\unclimit_histogram').format(plot_path)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pb.savefig(os.path.join(path_, file_name))
        plt.show
        plt.close()
    return
    
def unc_hist_cut(path):
    data_ = load_data('data_v3.txt')
    plot_path = make_directory(path)
    limit = 0.4  # CHANGE
    filters = data_.colnames
    for i in range(11, len(filters), 2):
        # Load wave and uncertainty
        band = filters[i]
        wave = data_[band]
        wave_unc = data_[band + '_unc']
        # Remove bad data
        bad_data = np.logical_and(wave != -99, wave_unc != -99)
        wave_trim = wave[bad_data]
        wave_unc_trim = wave_unc[bad_data]
        # Set limited data
        wave_limit = wave_trim[wave_unc_trim < limit]
        wave_unc_limit = wave_unc_trim[wave_unc_trim < limit]
        # Make histograms of wave and wave_limit
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        ax.hist(wave_unc_trim, bins=100, range=[min(wave_unc_trim), max(wave_unc_trim)],
                color='b')
        ax.xaxis.set_major_locator(plt.MultipleLocator(5))
        ax.hist(wave_unc_limit, bins=5, range=[min(wave_unc_limit), max(wave_unc_limit)],
                color='r')
        ax.xaxis.set_major_locator(plt.MultipleLocator(5))
        ax.set_title(band + ' Magnitude Frequency', fontsize=12)
        ax.set_xlabel(band + ' Magnitude')
        ax.set_ylabel('Frequency')
        # Display interactive figure if # removed, if not, figures saved
        file_name = 'unc_hist_{}-unc_limit-{}.png'.format(band, limit)
        path_ = ('{}\\unclimit_histogram\\unc_histograms').format(plot_path)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pb.savefig(os.path.join(path_, file_name))
        plt.show
        plt.close()
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

def load_data(surveyfile_):
    '''User upload data file'''

    d_file_name = str(surveyfile_)
    data = Table.read(d_file_name, format='ascii.commented_header',
                      guess=False)

    return (data)