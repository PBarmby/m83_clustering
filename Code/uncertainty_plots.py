import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pylab
import os
import os.path

'''Data is only trimmed by removing -99s. UNC limit can be set by limiting the
    y-axis of the plots'''


def unc_plot(path):
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    names = data.colnames
    s_path = make_directory(path)
    for i in xrange(11, len(names), 4):
        band_05 = data[names[i]]
        band_3 = data[names[i+2]]
        band_05_unc = data[names[i+1]]
        band_3_unc = data[names[i+3]]
        data_trim05 = np.logical_and(band_05 != -99, band_05_unc != -99)
        data_trim3 = np.logical_and(band_3 != -99, band_3_unc != -99)
        data_trim = np.logical_and(data_trim05, data_trim3)
        wave_uncertainty(band_05, band_3, band_05_unc, band_3_unc,
                         data_trim, names[i], names[i+2], s_path)

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


def wave_uncertainty(wave_05, wave_3, wave_05_unc, wave_3_unc, greatdata,
                     band_05, band_3, save):

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(211)
    ax.scatter(wave_05[greatdata], wave_05_unc[greatdata],
               color='k', marker='.')
    plt.setp(ax.get_xticklabels(), visible=False)
    #ax.set_xlabel(band_05)
    ax.set_ylabel(band_05+'unc')
    ax.set_title(band_05+' vs. '+band_05+'_unc', fontsize=12)
    
    ax2 = fig.add_subplot(212, sharex=ax)
    ax2.scatter(wave_3[greatdata], wave_3_unc[greatdata],
               color='k', marker='.')
    ax2.set_ylabel(band_3+'unc')
    ax2.set_xlabel(band_05 + ' (top) - ' + band_3 + ' (bottom)')
    plt.setp(ax2.get_xticklabels(), fontsize=6)
    plt.gca().set_ylim(0, 0.2)
    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'uncertainty_02_lim_{}.png'.format(band_05)

    pylab.savefig(os.path.join(save, file_name))
    plt.close()
    return
