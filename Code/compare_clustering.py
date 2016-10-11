'''Compare the results of two classifications
1. Run from Spyder terminal
    Argument 1 (dim): '2d' or '3d' or '2a3' - puts dimensions of each clustering in file name
    Argument 2 (df_1): first id_ file from clustering
    Argument 3 (df_2): second id_ file from clustering
    Argument 4 (d1_path): path to df_1
    Argument 5 (d2_path): path to df_2
    Argument 6 (save_path): set path to save figures
    Argument 7 (plots): choose plots ('corr' or 'broad_pos') ONLY CORR WORKS
        Corr creates plots with bubbles comparing the objects in each clustering
        broad_pos was going to plot the objects in [U-B] vs. [V-I] space
            (probably not useful)

IMPORTANT: This file gets the BAND names from the name of the id_ file in parameters() function
                Do not change the format of the id_ file titles'''

import os
import numpy as np
import pylab
from matplotlib import pyplot as plt
from astropy.table import Table, join
# df_1: name of id_file 
# df_2: name of id_file
# d(n)_path: path to data files
# save_path: path to save plots (in results directory)
# plots: 'corr', 'broad_pos'

def compare_class(dim, df_1, df_2, d1_path, d2_path, save_path, plots):
    survey = Table.read('data_v3.txt', format='ascii.commented_header',
                        guess=False)

    # Load parameters from file titles
    bands_1, n_clust_1, alg_1, bands_2, n_clust_2, alg_2 = parameters(df_1,
                                                                      df_2)

    # Load two clusterings
    path1 = load_path(d1_path)
    path2 = load_path(d2_path)

    # Create saving path
    s_path = load_path(save_path)
    mk_dir(s_path)  # make the saving directory if not made
    label_data, clust_1data, clust_2data = data(df_1, df_2, survey, path1,
                                                path2, plots)
    # Make plots
    if "corr" in plots:
        corr_plot(label_data, s_path, bands_1, n_clust_1, alg_1, bands_2,
                  n_clust_2, alg_2, dim)
    if "broad_pos" in plots:
        broad_plot(clust_1data, clust_2data, s_path)

    return


def parameters(file1, file2):
    # Check clustering in file1 and assign starting position for bands in string
    if 'kmeans' in file1:
        clustering_1 = 'kmeans'
        start1 = 9
    if 'meanshift' in file1:
        clustering_1 = 'meanshift'
        start1 = 12
    # Check clustering in file1 and assign starting position for bands in string
    if 'kmeans' in file2:
        clustering_2 = 'kmeans'
        start2 = 9
    if 'meanshift' in file2:
        clustering_2 = 'meanshift'
        start2 = 12

    # file 1 parameters
    n_clust1 = int(file1[start1+1])
    wave1_1 = file1[start1+5:start1+14]
    wave2_1 = file1[start1+15:start1+24]
    wave3_1 = file1[start1+26:start1+35]
    wave4_1 = file1[start1+36:start1+45]
    file1_waves = np.array([wave1_1, wave2_1, wave3_1, wave4_1])

    # file 2 parameters
    n_clust2 = int(file2[start2+1])
    wave1_2 = file2[start2+5:start1+14]
    wave2_2 = file2[start2+15:start1+24]
    wave3_2 = file2[start2+26:start1+35]
    wave4_2 = file2[start2+36:start1+45]
    file2_waves = np.array([wave1_2, wave2_2, wave3_2, wave4_2])

    return(file1_waves, n_clust1, clustering_1, file2_waves, n_clust2,
           clustering_2)
    
    
def data(clustering1, clustering2, survey_data, path_1, path_2, plts):
    clust1 = Table.read(path_1 + '\\' + clustering1, format='ascii.commented_header',
                        guess=False)
    clust2 = Table.read(path_2 + '\\' + clustering2, format='ascii.commented_header',
                        guess=False)
    clust1_data = join(clust1, survey_data, 'object_id')
    clust2_data = join(clust2, survey_data, 'object_id')
    clust_labels = join(clust1, clust2, 'object_id')

    return(clust_labels, clust1_data, clust2_data)


def mk_dir(path_):
    if not os.path.exists(path_):
        os.makedirs(path_)
    return


def load_path(path_):
    path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}'.format(path_)
    return(path)


def corr_plot(labels, path_, waves1, nclust1, alg1, waves2, nclust2, alg2,
              dimensions):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Loop over n_clust in first clustering
    for i in range(0, nclust1):
        # Loop over n_clust in second clustering
        for k in range(0, nclust2):
            # Find objects in the same clusters of clustering 1 and 2
            objects = len(labels[np.logical_and(labels['cluster_number_1'] == i,
                                                labels['cluster_number_2'] == k)])
            ax.scatter(i+1, k+1, s=objects, color='0.75')
            ax.annotate(str(objects), xy=(i+1, k+1), xytext=(i+1.05, k+1.05))
    # Format plot
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    x_label = '{}: {}cl {}-{} vs. {}-{} {} dim'.format(alg1, nclust1,
                                                       waves1[0], waves1[1],
                                                       waves1[2], waves1[3],
                                                       dimensions)
    ax.set_xlabel(x_label, fontsize=8)

    ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
    y_label = '{}: {}cl {}-{} vs. {}-{} {} dim'.format(alg2, nclust2,
                                                       waves2[0], waves2[1],
                                                       waves2[2], waves2[3],
                                                       dimensions)
    ax.set_ylabel(y_label, fontsize=8)

    plt.show()
    file_name = '{}-{}cl_{}-{}_vs_{}-{}cl_{}-{}_{}-{}_{}dim_compare.png'.format(alg1,
                                                                nclust1,
                                                                waves1[0],
                                                                waves1[1],
                                                                alg2,
                                                                nclust2,
                                                                waves2[0],
                                                                waves2[1],
                                                                waves2[2],
                                                                waves2[3],
                                                                dimensions)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return


def broad_plot():
    return

