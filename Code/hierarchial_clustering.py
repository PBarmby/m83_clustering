import numpy as np 
import matplotlib.pyplot as plt
import pylab
from astropy.table import Table

def data(): 
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    m83_data = data[:5000]
    band1 = 'mag05_555'
    band2 = 'mag05_814'
    band3 = 'mag05_225'
    band4 = 'mag05_336'

    ratio = 0.05
    #eps = 0.05
    #min_samp = 10
# Colour 1
    wave1 = m83_data[band1]
    wave1_unc = m83_data[band1+'_unc']
    wave2 = m83_data[band2]
    wave2_unc = m83_data[band2+'_unc']
# Colour 2
    wave3 = m83_data[band3]
    wave3_unc = m83_data[band3+'_unc']
    wave4 = m83_data[band4]
    wave4_unc = m83_data[band4+'_unc']
# Change parameters to match data_file
# Remove data pieces with no value

    wave1_trim = np.logical_and(wave1 != -99, wave1_unc != -99)
    wave2_trim = np.logical_and(wave2 != -99, wave2_unc != -99)
    wave3_trim = np.logical_and(wave3 != -99, wave3_unc != -99)
    wave4_trim = np.logical_and(wave4 != -99, wave4_unc != -99)

    colour1_ratio = np.logical_and(wave1_unc/wave1 < ratio,
                                   wave2_unc/wave2 < ratio)
    colour2_ratio = np.logical_and(wave3_unc/wave3 < ratio,
                                   wave4_unc/wave4 < ratio)

    gooddata1 = np.logical_and(np.logical_and(wave1_trim, wave2_trim),
                               np.logical_and(wave3_trim, wave4_trim))

# Remove data above certain magnitude
    gooddata2 = np.logical_and(colour1_ratio, colour2_ratio)
# Only data that match criteria for both colours
    greatdata = np.logical_and(gooddata1, gooddata2)
    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]
    print "Objects: {}".format(len(colour1))

    X = np.vstack([colour1, colour2]).T
    #heatmap(X)
    return(X)

def heatmap(dm):
    link = 'complete'
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist, squareform
    
    D1 = dm #squareform(pdist(dm, metric='euclidean'))
    D2 = dm.T #squareform(pdist(dm.T, metric='euclidean'))
    
    f = plt.figure(figsize=(8, 8))

    # add first dendrogram
    ax1 = f.add_axes([0.09, 0.1, 0.2, 0.6])
    Y = linkage(D1, method=link)
    Z1 = dendrogram(Y, orientation='right',show_leaf_counts=False, no_labels=True)
    ax1.set_xticks([])
    ax1.set_yticks([])

    # add second dendrogram
    ax2 = f.add_axes([0.3, 0.71, 0.6, 0.2])
    Y = linkage(D2, method=link)
    Z2 = dendrogram(Y, show_leaf_counts=False, no_labels=True)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # add matrix plot
    axmatrix = f.add_axes([0.3, 0.1, 0.6, 0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D1[idx1, :]
    D = D[:, idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap='hot')#pylab.cm.YlGnBu)#'hot')
    #axmatrix.set_xticks(range(len(Z1['leaves'])))
    axmatrix.set_xticks([])
    #axmatrix.set_xticklabels(Z1['leaves'])
    axmatrix.set_yticks([])
    
    axcolor = f.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im,cax=axcolor)
    f.show()

    return Y, Z2#, D1,Z2"""
