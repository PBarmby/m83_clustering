import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.table import Table
import os

def wave_unc_hist(path_):
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    #for i in range(11, len(data.colnames), 2):
    band1 = 'mag3_657'
    wave1 = data[band1]
    wave1_unc = data[band1+'_unc']
    wave1_trim = np.logical_and(wave1 != -99, wave1_unc != -99)
    x = wave1[wave1_trim]
    y = wave1_unc[wave1_trim]

    fig, axScatter = plt.subplots(figsize=(8, 8))
    axScatter.scatter(x, y)

    axScatter.set_xlabel(band1+' Magnitude')
    axScatter.set_ylabel(band1+' Uncertainty')
    divider = make_axes_locatable(axScatter)
    axHistx = divider.append_axes("top", 1.2, pad=0.3, sharex=axScatter)
    axHisty = divider.append_axes("right", 1.2, pad=0.3, sharey=axScatter)

    # make some labels invisible
    plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
             visible=True)
         
    # now determine nice limits by hand:
    binwidth = 1.0
    xymax = np.max(x)
    lim = (int(xymax/binwidth) + 1)*binwidth
    bins = np.arange(0, lim + binwidth, binwidth)
    axHistx.hist(wave1[wave1_trim], bins=500, range=[min(wave1[wave1_trim]),
                 max(wave1[wave1_trim])])
    axHistx.set_title(band1+' Magnitude vs. Uncertainty')
    axHisty.hist(wave1_unc[wave1_trim], bins=500, range=[min(wave1_unc[wave1_trim]),
                 max(wave1_unc[wave1_trim])], orientation='horizontal')

    axHisty.set_xticklabels([0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000], rotation=270)
    # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.

    #axHistx.axis["bottom"].major_ticklabels.set_visible(False)
    for tl in axHistx.get_xticklabels():
        tl.set_visible(True)
    #axHistx.set_yticks([ 20, 40])


    #axHisty.axis["left"].major_ticklabels.set_visible(False)
    for tl in axHisty.get_yticklabels():
        tl.set_visible(True)

    file_name = "{}_uncertainty_distribution".format(band1)
    path = "C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}".format(path_)
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(os.path.join(path, file_name))
    plt.draw()
    plt.show()
    plt.close()
        
    return()