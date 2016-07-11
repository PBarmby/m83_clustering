import numpy as np
from astropy.table import Table

def data():
    '''Select data for analysis'''
    data = Table.read('data_v3.txt', format='ascii.commented_header',
                      guess=False)
    exp = Table.read('experiments.txt', format='ascii.commented_header',
                     guess=False)
    # ratio = 0.2
    data = data[:10000]
    exp_col = exp.colnames
    n_bands = 0

    for i in range(0, len(exp_col)):  # Find number of bands used
        if 'band' in exp_col[i]:
            n_bands = i + 1

    n_cols = n_bands/2
    bands = np.empty((len(data), n_bands*2))  # make array of zeros with size of data, bands
    colours = np.arange(0, n_cols)  # Make array of zeros with size of data, colours

    for e in range(0, len(exp)):
        bands = exp[0]
        for b in range(0, n_bands*2, 2):
            bands[:,b] = float(data[exp[exp_col[b]][e]])
            bands[:,b+1] = float(data[exp[exp_col[b]][e] + '_unc'])
        remove = []
        for t in range(0, len(bands)):
            if -99. in bands[t]:
                remove.append(t)
        print bands
    return


def load_data(surveyfile_, experiments):
    '''User upload data file'''

    d_file_name = str(surveyfile_)
    t_file_name = str(experiments)
    data = Table.read(d_file_name, format='ascii.commented_header',
                      guess=False)
    tests = Table.read(t_file_name, format='ascii.commented_header',
                       guess=False)

    return (data, tests)