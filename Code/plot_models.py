'''Plot model colours from FSPS'''
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table

def make_plots(file_name, band1, band2, band3, band4):
    load_file = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\model_colours\\{}'.format(file_name)
    #survey_file = Table.read('data_v3.txt', format='ascii.commented_header',
     #                        guess=False)
    model = Table.read(load_file, format='ascii.commented_header', guess=False)
    colour1 = model[band1] - model[band2]
    colour2 = model[band3] - model[band4]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(colour1, colour2)
    ax.set_xlabel(band1 + ' - ' + band2)
    ax.set_ylabel(band3 + ' - ' + band4)
    return