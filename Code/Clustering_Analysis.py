# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 13:24:41 2015
@author: Owner
"""

'''To Do'''
# Fix results_summary function

import os
import os.path
import argparse
import numpy as np
import pylab as pylab
from matplotlib import pyplot as plt
from astropy.table import Table, Column
import shutil

# Kmeans imports
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn import metrics 
from sklearn.metrics import pairwise_distances, silhouette_samples
from matplotlib.patches import Ellipse
from scipy.stats import norm

# meanshift imports
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import preprocessing

# mst imports
from scipy import sparse
from sklearn.mixture import GMM
from astroML.clustering import HierarchicalClustering, get_graph_segments

# used for plots
cluster_colours = ['y', 'g', 'b', 'r', 'c', 'm', 'k', 'w', 'brown', 'darkgray', 'orange', 'pink','gold', 'lavender', 'salmon', 'g', 'b',
                   'r', 'c', 'm', 'k', 'm', 'w', 'y', 'g', 'b', 'r', 'c', 'm',
                   'k','m','w','y','g','b','r','c','m','k','m','w','y','g','b',
                   'r','c','m','k','m','w','y','g','b','r','c','m','k','m','w',
                   'y','g','b','r','c','m','k','m','w','y','g','b','r','c','m',
                   'k','m','w','y','g','b','r','c','m','k','m','w','y','g','b',
                   'r','c','m','k','m','w','y','g','b','r','c','m','k','m','w',
                   'y','g','b','r','c','m','k','m','w']

# need this so that output files always have the same number of columns
max_num_clusters = 60

# Choose analysis and output


def clustering(save_plots, save_results, analysis, kmeans_input, plots,
               id_list, results_, data_file, input_file='experiments.txt'):
    '''DESCRIBE PROCESS HERE'''

    # Create saving directories
    plot_path, results_path = make_save_directory(save_plots, save_results)

    # User criteria
    analysis_criteria = analysis
    kmeans_input = kmeans_input
    generate_plots = plots
    id_output = id_list
    generate_results_summary = results_

    # Load data file
    data = load_data_file(data_file)

    # Import experiments listed in experiments.txt file
    run = np.genfromtxt(input_file, dtype='str')
    run_results = 'no'
    # Run analysis
    for i in range(0, len(run)):
        '''Get Data'''
        colour1, colour2, greatdata, x_data, y_data, id_data = \
            organize_data(run[i], data)

        results_title = ' '
        if "mag05" in run[i][0]:
            results_title = '05aperture_results.txt'
        elif "mag3" in run[i][0]:
            results_title = '3aperture_results.txt'
        '''Run Analysis'''
        
        numberofclusters = 0
        total_obj = 0
        num_obj = 0
        silhouette_score = 0
        mst_scale = 0
        bandwidth = 0
        if "meanshift" in analysis_criteria:
            numberofclusters, bandwidth = do_meanshift(plot_path, run[i, 0],
                                                       run[i, 1], run[i, 2],
                                                       run[i, 3], colour1,
                                                       colour2, generate_plots)
            run_results = 'yes'
        if "kmeans" in analysis_criteria:
            if "experiments.txt" in kmeans_input:
                numberofclusters = int(run[i, 4])
            silhouette_score, num_obj = do_kmeans(plot_path, run[i, 0],
                                                  run[i, 1], run[i, 2],
                                                  run[i, 3], colour1, colour2,
                                                  greatdata, numberofclusters,
                                                  generate_plots, id_output,
                                                  x_data, y_data, id_data)
            run_results = 'yes'
            total_obj = num_obj.sum()
        if "mst" in analysis_criteria:
            mst_scale = mst_clustering(generate_plots, x_data, y_data,
                                       run[i, 0], run[i, 1], run[i, 2],
                                       run[i, 3], plot_path)
            run_results = 'yes'
        '''WRITE RESULTS FILE'''
        if run_results == 'yes':
            write_results(results_path, results_title, analysis_criteria,
                          run[i, 0], run[i, 1], run[i, 2], run[i, 3],
                          numberofclusters, silhouette_score, total_obj,
                          num_obj, mst_scale, generate_results_summary,
                          bandwidth)
    if run_results == 'yes':
        shutil.copy2('C:\Users\Owner\Documents\GitHub\m83_clustering\Code\experiments.txt',
                     plot_path+'\\colours.txt')
    if 'yes' in generate_results_summary:
        results_summary(results_path, results_title)
    return()


def make_save_directory(p_path, r_path):
    '''Save results of each analysis
            save_path: set path of where you would like results saved'''
    # Create new analysis_folder
    res_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(r_path)
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    pl_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(
              p_path)
    if not os.path.exists(pl_path):
        os.makedirs(pl_path)

    return(pl_path, res_path)


def load_data_file(file_):
    '''User upload data file'''

    file_name = str(file_)
    data = Table.read(file_name, format='ascii.commented_header',
                      guess=False)

    return (data)


def organize_data(band_combinations, data_file):
    '''Select data for analysis'''
    # Colour 1
    data_file = data_file[:10000]
    wave1 = data_file[band_combinations[0]]
    wave2 = data_file[band_combinations[1]]
    # Colour 2
    wave3 = data_file[band_combinations[2]]
    wave4 = data_file[band_combinations[3]]

    # Change parameters to match data_file
    # Remove data pieces with no value
    gooddata1 = np.logical_and(np.logical_and(wave1 != -99, wave2 != -99),
                               np.logical_and(wave3 != -99, wave4 != -99))
    # Remove data above certain magnitude
    gooddata2 = np.logical_and(np.logical_and(wave1 < 26, wave2 < 26),
                               np.logical_and(wave3 < 26, wave4 < 26))
    # Only data that match criteria for both colours
    greatdata = np.logical_and(gooddata1, gooddata2)

    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]

    x = data_file['x'][greatdata]
    y = data_file['y'][greatdata]
    id_ = data_file['id'][greatdata].astype(np.int32)

    return (colour1, colour2, greatdata, x, y, id_)


def do_meanshift(s_path, band1, band2, band3, band4, colour1, colour2,
                 make_plot):
    '''Meanshift clustering to determine the number of clusters in the
        data, which is passed to KMEANS function'''
    # Truncate data
    X = np.vstack([colour1, colour2]).T
    '''Compute clustering with MeanShift'''
    # Scale data because meanshift generates circular clusters
    X_scaled = preprocessing.scale(X)
    # The following bandwidth can be automatically detected using
    # the routine estimate_bandwidth(X). Bandwidth can also be set manually.
    bandwidth = estimate_bandwidth(X)
    #bandwidth = 0.65
    # Meanshift clustering
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False)
    ms.fit(X_scaled)
    labels_unique = np.unique(ms.labels_)

    objects = ms.labels_[ms.labels_ >= 0]
    n_clusters = len(labels_unique[labels_unique >= 0])
    # Make plot
    if "meanshift" in make_plot:
        make_ms_plots(s_path, colour1, colour2, n_clusters, X, ms,
                      band1, band2, band3, band4, objects)
    return(n_clusters, bandwidth)


def do_kmeans(s_path, band1, band2, band3, band4, colour1, colour2, greatdata,
              number_clusters, make_plots, output_cluster_id, x, y, id_data):

    '''do K-means clustering on colours constructed from HST photometry band1,
    band 2, band3, band4 --- ie, names of HST  filters'''

    # Put data in format for clustering
    clusterdata = np.vstack([colour1, colour2]).T

    '''Compute KMeans Clustering'''
    # Data pre-processing
    scaler = preprocessing.StandardScaler()
    clf = KMeans(number_clusters, random_state=10)
    clf.fit(scaler.fit_transform(clusterdata))
    cluster_number = clf.predict(scaler.fit_transform(clusterdata))

    # Output object and cluster IDs to ID.txt file
    if "yes" in output_cluster_id:
        id_catologue(number_clusters, cluster_number, band1, band2, band3,
                     band4, id_data, s_path)

    # Compute the silhouette score
    labels = clf.labels_
    score = 0
    if number_clusters > 1:
        score = metrics.silhouette_score(scaler.fit_transform(clusterdata),
                                         labels)
    else:
        score = 'N/A'
        #sample_score = silhouette_samples(scaler.fit_transform(clusterdata), labels)
        #print sample_score
        #print "Average: {}".format(score)
    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        objects_per_cluster[i] = len(x_cluster)
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort

    # Generate plots if necessary
    if "kmeans_color" in make_plots:
        colour_kmeans_plot(s_path, band1, band2, band3, band4, clf, scaler,
                           colour1, colour2, number_clusters)
    if "kmeans_xy" in make_plots:
        xy_plot(s_path, colour1, colour2, number_clusters, cluster_number,
                band1, band2, band3, band4)

    return(score, objects_per_cluster)


def id_catologue(number_clusters, cluster_number, band1, band2, band3, band4,
                 id_data, save_):
    '''Create file with list of object ID and cluster number'''
    file_name = 'ID_{}cl_{}-{}vs{}-{}.txt'.format(str(number_clusters),
                                                  band1, band2, band3, band4)
    id_path = os.path.join(save_, file_name)
    id_tab = Table(data=[id_data, cluster_number],
                   names=['object_id', 'cluster_number'])
    id_tab.write(id_path, format='ascii.commented_header')

    return()


def mst_clustering(make_plots, x, y, band1, band2, band3, band4, p_path):

    X = np.vstack([x, y]).T

    # Boundaries for plots
    xmin, xmax = (min(x), max(x))
    ymin, ymax = (min(y), max(y))

    '''Compute MST clustering'''

    n_neighbors = 5
    edge_cutoff = 0.9
    cluster_cutoff = 20
    model = HierarchicalClustering(n_neighbors=n_neighbors,
                                   edge_cutoff=edge_cutoff,
                                   min_cluster_size=cluster_cutoff)
    model.fit(X)
    scale = np.percentile(model.full_tree_.data, 100 * edge_cutoff)
    n_components = model.n_components_
    labels = model.labels_

    # Get the x, y coordinates of the beginning and end of each line segment
    T_x, T_y = get_graph_segments(model.X_train_, model.full_tree_)
    T_trunc_x, T_trunc_y = get_graph_segments(model.X_train_,
                                              model.cluster_graph_)

    # Fit a GMM to each individual cluster
    Nx = 100
    Ny = 250
    Xgrid = np.vstack(map(np.ravel, np.meshgrid(np.linspace(xmin, xmax, Nx),
                                                np.linspace(ymin,
                                                            ymax, Ny)))).T
    density = np.zeros(Xgrid.shape[0])
    for i in range(n_components):
        ind = (labels == i)
        Npts = ind.sum()
        Nclusters = min(12, Npts / 5)
        gmm = GMM(Nclusters).fit(X[ind])
        dens = np.exp(gmm.score(Xgrid))
        density += dens / dens.max()
    density = density.reshape((Ny, Nx))

    if "mst" in make_plots:
        mst_plots(X, ymin, ymax, xmin, xmax, T_x, T_y, T_trunc_x, T_trunc_y,
                  density, band1, band2, band3, band4, p_path)

    return(scale)


def make_ms_plots(path, colour1, colour2, n_clusters, X, ms,
                  band1, band2, band3, band4, cluster_number):
    ''' Plot the results of mean shift clustering if needed'''
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    # plot density
    H, C1_bins, C2_bins = np.histogram2d(colour1, colour2, 51)
    ax.imshow(H.T, origin='lower', interpolation='nearest', aspect='auto',
              extent=[C1_bins[0], C1_bins[-1], C2_bins[0], C2_bins[-1]],
              cmap=plt.cm.binary)

    # plot clusters
    for i in range(n_clusters):
        Xi = X[ms.labels_ == i]
        H, b1, b2 = np.histogram2d(Xi[:, 0], Xi[:, 1], (C1_bins, C2_bins))
        bins = [0.1]
        ax.contour(0.5 * (C1_bins[1:] + C1_bins[:-1]), 0.5 * (C2_bins[1:]
                                                              + C2_bins[:-1]),
                   H.T, bins, colors='r')

    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(band1+' - '+band2)
    ax.set_ylabel(band3+' - '+band4)
    ax.set_title('mean-shift: '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                 fontsize=16)

    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'meanshift_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                         band1, band2,
                                                         band3, band4)
    pylab.savefig(os.path.join(path, file_name))
    plt.close()

    return()


def colour_kmeans_plot(path, band1, band2, band3, band4, clf, scaler, colour1,
                       colour2, number_clusters):

    '''Plot cluster data for KMEANS clustering'''
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    # Compute 2D histogram of the input
    H, C1_bins, C2_bins = np.histogram2d(colour1, colour2, 50)

    # Plot Density

    ax.imshow(H.T, origin='lower', interpolation='nearest', aspect='auto',
              extent=[C1_bins[0], C1_bins[-1],
                      C2_bins[0], C2_bins[-1]], cmap=plt.cm.binary)

    # Plot Cluster centers
    cluster_centers = scaler.inverse_transform(clf.cluster_centers_)
    for i in range(0, number_clusters):
        ax.scatter(cluster_centers[i, 0], cluster_centers[i, 1], s=40,
                   c=cluster_colours[i], edgecolors='k')

    # Plot cluster centers
    C1_centers = 0.5 * (C1_bins[1:] + C1_bins[:-1])
    C2_centers = 0.5 * (C2_bins[1:] + C2_bins[:-1])
    clusterdatagrid = np.meshgrid(C1_centers, C2_centers)
    clusterdatagrid = np.array(clusterdatagrid).reshape((2, 50 * 50)).T
    H = clf.predict(scaler.transform(clusterdatagrid)).reshape((50, 50))

    # Plot boundries
    for i in range(number_clusters):
        Hcp = H.copy()
        flag = (Hcp == i)
        Hcp[flag] = 1
        Hcp[~flag] = 0
        ax.contour(C1_centers, C2_centers, Hcp, [-0.5, 0.5], linewidths=1,
                   c='k')

    # Set axes for each plot
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))

    ax.set_xlabel(band1+' - '+band2)
    ax.set_ylabel(band3+' - '+band4)
    ax.set_title('k-means: '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                 fontsize=12)
    file_name = 'k_means_{}cl_{}-{}vs{}-{}.png'.format(str(number_clusters),
                                                       band1, band2, band3,
                                                       band4)
    pylab.savefig(os.path.join(path, file_name))
    plt.close()

    return ()


def xy_plot(path, x, y, number_clusters, cluster_number, band1, band2, band3,
            band4):
    '''Plot xy positions of objects in different clusters'''
    fig2 = plt.figure(figsize=(8, 8))
    ax2 = fig2.add_subplot(111)

    # Plot XY coordinates of each object and label with colour that
    # identifies cluster number
    ax2.set_xlim(min(x), max(x))
    ax2.set_ylim(min(y), max(y))
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        y_cluster = y[cluster_number == i]
        ax2.scatter(x_cluster, y_cluster, label=i, c=cluster_colours[i])

    ax2.set_xlabel(band1+' - '+band2)
    ax2.set_ylabel(band3+' - '+band4)
    ax2.set_title('kmeans: '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                  fontsize=16)
    ax2.legend(loc='upper left')

    file_name = 'xy_{}cl_{}-{}vs{}-{}.png'.format(str(number_clusters),
                                                  band1, band2, band3, band4)
    pylab.savefig(os.path.join(path, file_name))
    plt.close()

    return ()


def mst_plots(X, ymin, ymax, xmin, xmax, T_x, T_y, T_trunc_x, T_trunc_y,
              density, band1, band2, band3, band4, path):
    # Plot the results
    fig = plt.figure(figsize=(8, 12))
    fig.subplots_adjust(hspace=0, left=0.1, right=0.95, bottom=0.1, top=0.9)

    ax = fig.add_subplot(311, aspect='equal')
    ax.scatter(X[:, 1], X[:, 0], s=1, lw=0, c='k')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_ylabel('(Pix)')

    ax = fig.add_subplot(312, aspect='equal')
    ax.plot(T_y, T_x, c='k', lw=0.5)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_ylabel('(Pix)')

    ax = fig.add_subplot(313, aspect='equal')
    ax.plot(T_trunc_y, T_trunc_x, c='k', lw=0.5)
    ax.imshow(density.T, origin='lower', cmap=plt.cm.hot_r, extent=[ymin, ymax,
                                                                    xmin,
                                                                    xmax])

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('(Mpc)')
    ax.set_ylabel('(Mpc)')

    filename = 'mst_{}-{}vs{}-{}.png'.format(band1, band2, band3, band4)
    pylab.savefig(os.path.join(path, filename))
    # plt.show()

    return()


def write_results(save_path, name, methods, band1, band2, band3, band4,
                  numberofclusters, silhouette_score, total_obj, num_obj,
                  mst_scale, results_summary_gen, ms_bandwidth):

    file_path = '{}'.format(save_path)
    test_path = '{}\\{}'.format(file_path, name)
    header = '# meanshift kmeans mst band1 band2 band3 band4 bandwidth mst_scale number_of_clusters silhouette_score total_objects c_1 c_2 c_3 c_4 c_5 c_6 c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40 c_41 c_42 c_43 c_44 c_45 c_46 c_47 c_48 c_49 c_50 c_51 c_52 c_53 c_54 c_55 c_56 c_57 c_58 c_59 c_60'
    if not os.path.exists(test_path):
        header_path = os.path.join(file_path, name)
        results_file = open(header_path, "a")
        results_file.write(header + '\n')
        results_file.close()

    ms = "YES"
    km = "YES"
    mst = "YES"
    if "meanshift" not in methods:
        ms = "NO"
        ms_bandwidth = "N/A"
    if "kmeans" not in methods:
        km = "NO"
        silhouette_score = "N/A"
        total_obj = "N/A"
        num_obj = "N/A"
    if "mst" not in methods:
        mst = "NO"
        mst_scale = "N/A"
    name_path = os.path.join(file_path, name)
    results_file = open(name_path, "a")
    if ms == "YES":
        inputs = '{} {} {} {} {} {} {} {:.4f} {} {}'.format(ms, km, mst, band1,
                                                    band2, band3, band4,
                                                    ms_bandwidth, mst_scale,
                                                    numberofclusters)
    else: 
        inputs = '{} {} {} {} {} {} {} {} {} {}'.format(ms, km, mst, band1,
                                                    band2, band3, band4,
                                                    ms_bandwidth, mst_scale,
                                                    numberofclusters)
    if km == "YES":
        outputs = ' {:.4f} {:5d} {}'.format(silhouette_score, total_obj,
                                            np.array_str(num_obj,
                                                         max_line_width=500)[1:-1])
    else:
        outputs = "N/A N/A {}".format(np.array_str(
                                      np.array([0]*max_num_clusters),
                                      max_line_width=500)[1:-1])
    results_file.write(inputs + ' ' + outputs + '\n')
    results_file.close()

    return()


def results_summary(path, input_file):
    '''Compute and plot summaries for clustering analysis'''

    results_file = os.path.join(path, input_file)

    # read in the data -- this is not an ideal way to do it since it requires
    # knowledge of file structure
    results_table = Table.read(results_file, format='ascii.commented_header',
                               guess=False)

    remove = []
    for i in range(0, len(results_table)):
        if results_table['silhouette_score'][i] == 'N/A':
            remove.append(i)
    results_table.remove_rows(remove)

    num_clust = results_table['number_of_clusters']
    score = results_table['silhouette_score']
    total_obj = results_table['total_objects'].astype('float')

    # add a column with the size of the smallest cluster
    # have to do some tricky stuff since column corresponding to smallest
    # cluster varies dep on number of clusters

    results_table.add_column(Column(name='size_smallest',
                                    data=np.zeros(len(results_table)),
                                    dtype=np.int16))

    for i in range(0, len(results_table)):
        lastcol = 'c_{}'.format(num_clust[i])
        results_table['size_smallest'][i] = results_table[lastcol][i]

    # compute fraction of objects in largest cluster
    biggest_clust_fract = results_table['c_1']/total_obj
    # compute fraction of objects in smallest cluster
    smallest_clust_fract = results_table['size_smallest']/total_obj

    fig = plt.figure(figsize=(12,5))
    ax = fig.add_subplot(121)
    ax.scatter(num_clust, score)
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Score')
    ax.set_title('Number of Clusters vs Silhouette Score', fontsize=11)

    ax = fig.add_subplot(122)
    ax.scatter(score, biggest_clust_fract, c='r', marker='o',
               label='Largest Cluster')
    ax.scatter(score, smallest_clust_fract, c='b', marker='o',
               label='Smallest Cluster')
    ax.legend(loc='upper left')
    ax.set_xlabel('Score')
    ax.set_ylabel('Fractional size')
    ax.set_title('Silhouette Score vs. Fractional Size of Largest/Smallest Cluster', fontsize=11)

    filename = 'silhouette_score_plots.png'
    pylab.savefig(os.path.join(path, filename))

    # Bar graph
    for i in range(1, max(num_clust)+1):
        numberoftrials = len(num_clust[num_clust == i])  # Find  trials
        if numberoftrials > 0:
            fig, ax = plt.subplots()
            X = np.arange(0, numberoftrials)
            for n in range(0, i):  # Create stacked bar graph
                if n == 0:
                    yprev = (0*X).astype('float')
                else:
                    yprev = Y + yprev
                colname = 'c_{}'.format(n+1)
                Y = results_table[colname][num_clust==i]/results_table['total_objects'][num_clust==i].astype('float')  # Compute percentage of objects in each cluster and graph
                ax.bar(X, Y, color=cluster_colours[n], bottom=yprev)
                ax.set_title('Proportional Cluster Size')
                ax.set_xlabel('Number of Trials')
                ax.set_ylabel('Fractional Cluster Size')
                ax.xaxis.set_major_locator(plt.MultipleLocator(1))
            filename = '{}_clusters_vs_FractionalSize.png'.format(i)
            pylab.savefig(os.path.join(path, filename))

    return()


def user_input():

    inputs = argparse.ArgumentParser(description="Inputs for clustering analysis.")

    # Mandatory arguments: data_file
    inputs.add_argument("data_file", help="Choose the data file for analysis")

    # Optional Arguments: Save paths for results text files and figure files
    inputs.add_argument("-pp", "--plot_path", help="Enter path for saving plot files", default="results")
    inputs.add_argument("-rp", "--results_path", help="Enter path for saving results files", default="results")

    # Optional Arguments: Choose the clustering method and inputs for kmeans
    inputs.add_argument("-a", "--analysis", help = "Choose the methods you would like to use",
                        choices=['meanshift', 'kmeans', 'mst'], nargs='*',
                        default=[])
    inputs.add_argument("-kmi", "--kmeans_input", help="Choose the number of clusters input for kmeans", 
                        choices=['meanshift', 'experiments.txt'],
                        default = ['experiments.txt'])
    
    # Optional Arguments: Choose plots
    inputs.add_argument("-p","--plots", help = "Choose the plots you would like to make", 
                        choices =['meanshift', 'kmeans_color', 'kmeans_xy', 'mst'], nargs='*', default=[])
                        
    # Optional Arguments: Choose other functions to run (id, results_summary)
    inputs.add_argument("-id", "--id_list", help = "Produces object id list", choices = ['yes','no'], default='no')
    inputs.add_argument("-rs", "--results_summary", help="Produces results summary", choices =['yes', 'no'], default='no')

    # Optional Argument: Choose file name of input and output files    
    inputs.add_argument("-fn", "--file_names", help="Specify input and ouput file names. Default: experiments.txt, results.txt", default = ['experiments.txt, results.txt'], nargs='*')
    
    criteria = inputs.parse_args()     
    
    clustering(criteria.plot_path, criteria.results_path, criteria.analysis,
               criteria.kmeans_input, criteria.plots, criteria.id_list,
               criteria.results_summary, criteria.data_file)                                                         

    return()

if __name__ == "__main__": 
    user_input()
    