import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pylab
import os
import os.path


def unc_plot(path):
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    names = data.colnames
    s_path = make_directory(path)
    for i in xrange(11, len(names), 2):
        band = data[names[i]]
        band_unc = data[names[i+1]]
        data_trim = np.logical_and(band != -99, band_unc != -99)
        wave_uncertainty(band, band_unc, data_trim, names[i], s_path)

    print "finished"
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


def wave_uncertainty(wave, wave_unc, greatdata, band, save):

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.scatter(wave[greatdata], wave_unc[greatdata], color='k', marker='.')
    ax.set_xlabel(band)
    ax.set_ylabel(band+'unc')
    ax.set_title(band+' vs. '+band+'_unc', fontsize=12)

    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'uncertainty_{}.png'.format(band)

    pylab.savefig(os.path.join(save, file_name))
    plt.close()
    return
