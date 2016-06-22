# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 13:24:41 2015
@author: Owner
"""

'''To Do'''
# Change results to individual clustering ??? 
# results_summary to take specific method or multiple 


import os
import os.path
import argparse
import numpy as np
import pylab as pylab
from matplotlib import pyplot as plt
from astropy.table import Table
import shutil
from itertools import cycle

# Kmeans imports
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn import metrics 
from sklearn.metrics import pairwise_distances, silhouette_samples
from scipy.stats import norm

# meanshift imports
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import preprocessing

# dbscan imports
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN

#affinity propagation imports
from sklearn.cluster import AffinityPropagation

'''-----------------------Global Variables----------------------------------'''
# used for plots
cluster_colours = ['y', 'g', 'b', 'r', 'c', 'm', 'k', 'g', 'b',
                   'r', 'c', 'm', 'k', 'm', 'y', 'g', 'b', 'r', 'c', 'm',
                   'k','m', 'y','g','b','r','c','m','k','m', 'y','g','b',
                   'r','c','m','k','m','y','g','b','r','c','m','k','m',
                   'y','g','b','r','c','m','k','m','y','g','b','r','c','m',
                   'k','m','y','g','b','r','c','m','k','m','y','g','b',
                   'r','c','m','k','m','y','g','b','r','c','m','k','m',
                   'y','g','b','r','c','m','k','m']
# need this so that output files always have the same number of columns
max_num_clusters = 60
'''-------------------------------------------------------------------------'''


def clustering(save_plots, save_results, analysis, kmeans_input, bw_in, plots,
               id_list, data_file, input_file='experiments.txt'):
    '''DESCRIBE PROCESS HERE'''
    # Create saving directories
    plot_path, results_path = make_save_directory(save_plots, save_results)

    # Load data and trials
    data, experiments = load_data_file(data_file, input_file)

    '''Run Analysis'''
    for i in range(0, len(experiments)):
        '''-----------------Metrics reset every trial-----------------------'''
        ms_n_clusters = 0
        db_n_clusters = 0
        af_n_clusters = 0
        km_n_clusters = 0
        total_obj = 0
        num_obj = 0
        ms_score = 0
        db_score = 0
        af_score = 0
        km_score = 0
        bandwidth = 0
        results_title = ' '
        '''-----------------------------------------------------------------'''

        # Get Data
        cluster_data_, greatdata, x_data, y_data, id_data = \
            organize_data(experiments[i], data)

        # Set title of results.txt
        if "mag05" in experiments['band1'][i]:
            results_title = '05aperture_results.txt'
        elif "mag3" in experiments['band1'][i]:
            results_title = '3aperture_results.txt'

        # Do clustering

        if "meanshift" in analysis:
            if 'experiments.txt' in bw_in:
                b_width_input = experiments['bandwidth'][i]
            else:
                b_width_input = 'estimate'

            ms_n_clusters, bandwidth, ms_score, ms_obj, ms_obj_p_cluster = \
                meanshift(plot_path, experiments[i], cluster_data_,
                             plots, b_width_input, id_list, id_data)

            meanshift_results(results_path, results_title, experiments[i],
                              ms_n_clusters, ms_score, bandwidth, ms_obj,
                              ms_obj_p_cluster)

        if "dbscan" in analysis:
            eps = float(experiments['eps'][i])
            min_samples = float(experiments['min_samples'][i])
            db_score, db_n_clusters, db_obj, db_obj_p_cluster = \
                dbscan(plot_path, experiments[i], cluster_data_,
                       plots, eps, min_samples)
            dbscan_results(results_path, results_title, experiments[i],
                           db_n_clusters, db_score, eps, min_samples, db_obj,
                           db_obj_p_cluster)

        if "affinity" in analysis:
            damping = float(experiments['damping'][i])
            preferences = float(experiments['preferences'][i])
            af_n_clusters, af_score, af_obj, af_obj_p_cluster = \
                affinity_propagation(plot_path, experiments[i], cluster_data_,
                                     plots, damping, preferences, id_list,
                                     id_data)
            affinity_propagation_results(results_path, results_title,
                                         experiments[i], af_n_clusters,
                                         af_score, damping, preferences,
                                         af_obj, af_obj_p_cluster)

        if "kmeans" in analysis:
            if "experiments.txt" in kmeans_input:
                km_n_clusters = int(experiments['n_clusters'][i])
            elif "meanshift" in kmeans_input:
                km_n_clusters = ms_n_clusters
            elif "dbscan" in kmeans_input:
                km_n_clusters = db_n_clusters
            elif "affinity" in kmeans_input:
                km_n_clusters = af_n_clusters
            if km_n_clusters >= 4:
                high, low = 3, 2
            elif km_n_clusters == 3:
                high, low = 4, 1
            else:
                high, low = 5, 0

            for a in range(km_n_clusters - low, km_n_clusters + high):
                km_score, num_obj = kmeans(plot_path, experiments[i],
                                              cluster_data_, greatdata,
                                              a, plots, id_list,
                                              x_data, y_data, id_data)
                total_obj = num_obj.sum()
                kmeans_results(results_path, results_title, experiments[i],
                               kmeans_input, a, km_score, total_obj, num_obj)

    # Copy experiments.txt to plots directory
    shutil.copy2('C:\Users\Owner\Documents\GitHub\m83_clustering\Code\experiments.txt',
                 plot_path+'\\inputs.txt')

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


def load_data_file(d_file, e_file):
    '''User upload data file'''

    d_file_name = str(d_file)
    e_file_name = str(e_file)
    data = Table.read(d_file_name, format='ascii.commented_header',
                      guess=False)
    trial = Table.read(e_file_name, format='ascii.commented_header',
                       guess=False)

    return (data, trial)


def organize_data(exp, data_file):
    '''Select data for analysis'''
    data = data_file
    ratio = 0.2

    wave1 = data[exp['band1']]
    wave1_unc = data[exp['band1']+'_unc']
    wave2 = data[exp['band2']]
    wave2_unc = data[exp['band2']+'_unc']
    # Colour 2
    wave3 = data[exp['band3']]
    wave3_unc = data[exp['band3']+'_unc']
    wave4 = data[exp['band4']]
    wave4_unc = data[exp['band4']+'_unc']
    # Colour 3
    wave5 = data[exp['band5']]
    wave5_unc = data[exp['band5']+'_unc']
    wave6 = data[exp['band6']]
    wave6_unc = data[exp['band6']+'_unc']
    # Colour 4
    wave7 = data[exp['band7']]
    wave7_unc = data[exp['band7']+'_unc']
    wave8 = data[exp['band8']]
    wave8_unc = data[exp['band8']+'_unc']
    # Colour 5
    wave9 = data[exp['band9']]
    wave9_unc = data[exp['band9']+'_unc']
    wave10 = data[exp['band10']]
    wave10_unc = data[exp['band10']+'_unc']
    # Colour 6
    wave11 = data[exp['band11']]
    wave11_unc = data[exp['band11']+'_unc']
    wave12 = data[exp['band12']]
    wave12_unc = data[exp['band12']+'_unc']
    #colour 7 
    wave13 = data[exp['band13']]
    wave13_unc = data[exp['band13']+'_unc']
    wave14 = data[exp['band14']]
    wave14_unc = data[exp['band14']+'_unc']
    # Change parameters to match data_file
    # Remove data pieces with no value
    
    wave1_trim = np.logical_and(np.logical_and(wave1 != -99, wave1_unc != -99), wave1_unc < ratio)
    wave2_trim = np.logical_and(np.logical_and(wave2 != -99, wave2_unc != -99), wave2_unc < ratio)

    wave3_trim = np.logical_and(np.logical_and(wave3 != -99, wave3_unc != -99), wave3_unc < ratio)
    wave4_trim = np.logical_and(np.logical_and(wave4 != -99, wave4_unc != -99), wave4_unc < ratio)

    wave5_trim = np.logical_and(np.logical_and(wave5 != -99, wave5_unc != -99), wave5_unc < ratio)
    wave6_trim = np.logical_and(np.logical_and(wave6 != -99, wave6_unc != -99), wave6_unc < ratio)

    wave7_trim = np.logical_and(np.logical_and(wave7 != -99, wave7_unc != -99), wave7_unc < ratio)
    wave8_trim = np.logical_and(np.logical_and(wave8 != -99, wave8_unc != -99), wave8_unc < ratio)

    wave9_trim = np.logical_and(np.logical_and(wave9 != -99, wave9_unc != -99), wave9_unc < ratio)
    wave10_trim = np.logical_and(np.logical_and(wave10 != -99, wave10_unc != -99), wave10_unc < ratio)

    wave11_trim = np.logical_and(np.logical_and(wave11 != -99, wave11_unc != -99), wave11_unc < ratio)
    wave12_trim = np.logical_and(np.logical_and(wave12 != -99, wave12_unc != -99), wave12_unc < ratio)

    wave13_trim = np.logical_and(np.logical_and(wave13 != -99, wave13_unc != -99), wave13_unc < ratio)
    wave14_trim = np.logical_and(np.logical_and(wave14 != -99, wave14_unc != -99), wave14_unc < ratio)

    colour1_trim = np.logical_and(wave1_trim, wave2_trim)
    colour2_trim = np.logical_and(wave3_trim, wave4_trim)
    colour3_trim = np.logical_and(wave5_trim, wave6_trim)
    colour4_trim = np.logical_and(wave7_trim, wave8_trim)
    colour5_trim = np.logical_and(wave9_trim, wave10_trim)
    colour6_trim = np.logical_and(wave11_trim, wave12_trim)
    colour7_trim = np.logical_and(wave13_trim, wave14_trim)

    gooddata1 = np.logical_and(colour1_trim, colour2_trim)
    gooddata2 = np.logical_and(colour3_trim, colour4_trim)
    gooddata3 = np.logical_and(colour5_trim, colour6_trim)

    # Only data that match criteria for both colours
    final_data = np.logical_and(np.logical_and(gooddata1, gooddata2),
                                np.logical_and(gooddata3, colour7_trim))

    # Make colours
    colour1 = wave1[final_data] - wave2[final_data]
    colour2 = wave3[final_data] - wave4[final_data]
    colour3 = wave5[final_data] - wave6[final_data]
    colour4 = wave7[final_data] - wave8[final_data]
    colour5 = wave9[final_data] - wave10[final_data]
    colour6 = wave11[final_data] - wave12[final_data]
    colour7 = wave13[final_data] - wave14[final_data]

    cluster_data = np.vstack([colour1, colour2, colour3, colour4, colour5,
                             colour6, colour7]).T
    x = data['x'][final_data]
    y = data['y'][final_data]
    id_ = data['id_'][final_data].astype(np.int32)

    return (cluster_data, final_data, x, y, id_)


def meanshift(s_path, bands, cluster_data, make_plot, bw_input,
              output_id, id_data):
    '''Meanshift clustering to determine the number of clusters in the
        data, which is passed to KMEANS function'''
    bandwidth = 0
    if 'experiments.txt' not in bw_input:
        # The following bandwidth can be automatically detected using
        bandwidth = estimate_bandwidth(cluster_data)
    else:
        bandwidth = bw_input

    X_scaled = preprocessing.scale(cluster_data)
    ms = MeanShift(bandwidth=bandwidth, cluster_all=True)
    ms.fit(X_scaled)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, n_clusters_):
        ith_cluster = X_scaled[labels == i]
        objects_per_cluster[i] = len(ith_cluster)
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    total_objects = len(X_scaled)

    # Compute silhouette_score for entire cluster or individuals
    if n_clusters_ > 1:
        average_score = metrics.silhouette_score(cluster_data, labels)
    else:
        average_score = -99
    # sample_score = metrics.silhouette_samples(X, labels)

    if "yes" in output_id:
        id_catologue(n_clusters_, labels, bands[0], bands[1], bands[2],
                     bands[3], id_data, s_path)
    # Make plot
    if "meanshift_density" in make_plot:
        meanshift_density(s_path, cluster_data, n_clusters_, cluster_data, ms,
                          bands[0], bands[1], bands[2], bands[4], labels,
                          cluster_centers)
    if "meanshift_colour" in make_plot:
        meanshift_colour(s_path, X_scaled, n_clusters_, labels,
                         cluster_centers, bands)

    return(n_clusters_, bandwidth, average_score, total_objects,
           objects_per_cluster)


def dbscan(s_path, band1, band2, band3, band4, colour1, colour2, make_plots,
           eps, min_samp):

    X = np.vstack([colour1, colour2]).T
    X = StandardScaler().fit_transform(X)

    # Compute DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samp).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    unique_labels = set(labels)

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, n_clusters_):
        ith_cluster = X[labels == i]
        objects_per_cluster[i] = len(ith_cluster)
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    total_objects = len(X)

    # Silhouette Score
    db_score = metrics.silhouette_score(X, labels)

    if 'dbscan' in make_plots:
        db_plot(unique_labels, labels, X, core_samples_mask, n_clusters_,
                band1, band2, band3, band4, s_path)
    return(db_score, n_clusters_, total_objects, objects_per_cluster)


def affinity_propagation(s_path, bands, cluster_data,
                         make_plots, damp, pref, output_id, id_data):

    # Compute Affinity Propagation
    af = AffinityPropagation(preference=pref, damping=damp).fit(cluster_data)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    n_clusters_ = len(cluster_centers_indices)

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, n_clusters_):
        ith_cluster = cluster_data[labels == i]
        objects_per_cluster[i] = len(ith_cluster)
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    total_objects = len(cluster_data)

    # Silhouette Score
    ap_score = metrics.silhouette_score(cluster_data, labels)
    if "yes" in output_id:
        id_catologue(n_clusters_, labels, bands, id_data, s_path)

    if 'affinity' in make_plots:
        affinity_plot(labels, cluster_centers_indices, cluster_data, n_clusters_, bands, s_path)

    return(n_clusters_, ap_score, total_objects, objects_per_cluster)


def kmeans(s_path, bands, cluster_data, greatdata,
              number_clusters, make_plots, output_cluster_id, x, y, id_data):

    '''do K-means clustering on colours constructed from HST photometry band1,
    band 2, band3, band4 --- ie, names of HST  filters'''

    # Put data in format for clustering
    clusterdata = cluster_data

    '''Compute KMeans Clustering'''
    # Data pre-processing
    scaler = preprocessing.StandardScaler()
    clf = KMeans(number_clusters, init='random')
    clf.fit(clusterdata) #change variable - data with labels
    #clf.fit(clusterdata)
    cluster_number = clf.predict(clusterdata)
    centers = clf.cluster_centers_
    #cluster_number = clf.predict(clusterdata)

    # Output object and cluster IDs to ID.txt file
    if "yes" in output_cluster_id:
        id_catologue(number_clusters, cluster_number, bands[0], bands[1], bands[2],
                     bands[3], id_data, s_path)

    # Compute the silhouette score
    labels = clf.labels_
    score = 0
    if number_clusters > 1:
        score = metrics.silhouette_score(clusterdata, labels)
    else:
        score = -99
        #sample_score = silhouette_samples(scaler.fit_transform(clusterdata), labels)
    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        objects_per_cluster[i] = len(x_cluster)
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort

    # Generate plots if necessary
    if "kmeans_density" in make_plots:
        kmeans_density(s_path, bands[0], bands[1], bands[2], bands[3], clf,
                       scaler, clusterdata, number_clusters, cluster_number)
    if "kmeans_colour" in make_plots:
        kmeans_colour(s_path, clusterdata, number_clusters, cluster_number,
                      bands, centers)

    return(score, objects_per_cluster)


def id_catologue(number_clusters, cluster_number, band1, band2, band3, band4,
                 id_data, save_):
    '''Create file with list of object ID and cluster number'''
    file_name = 'id_{}cl_{}-{}vs{}-{}.txt'.format(str(number_clusters),
                                                  band1, band2, band3, band4)
    id_path = os.path.join(save_, file_name)
    id_tab = Table(data=[id_data, cluster_number],
                   names=['object_id', 'cluster_number'])
    id_tab.write(id_path, format='ascii.commented_header')

    return()


def meanshift_density(path, cluster_data, n_clusters, X, ms,
                      band1, band2, band3, band4, labels_, centers):
    ''' Plot the results of mean shift clustering if needed'''
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    # plot density
    H, C1_bins, C2_bins = np.histogram2d(cluster_data[0], cluster_data[1], 51)
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
    ax.set_title('mean-shift ' + str(n_clusters) + ' : '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                 fontsize=14)

    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'ms_density_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                          band1, band2,
                                                          band3, band4)
    colors = ('{}-{}_{}-{}').format(band1, band2, band3, band4)
    path_ = ('{}\\meanshift-dens_{}').format(path, colors)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    return()


def meanshift_xy(path, X, n_clusters, labels_, centers, band1, band2, band3,
                 band4):
    
    return

def meanshift_colour(path, X, n_clusters, labels_, centers, bands):

    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    markers = cycle('ooooooo**')#*****<<<<<<>>>>>>>^^^^^^^')
    for i in range(1, 7):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for k, col, mark in zip(range(n_clusters), colors, markers):
            # cluster_score = sample_score[labels == k]
            my_members = labels_ == k
            cluster_center = centers[k]
            ax.scatter(X[my_members, 0], X[my_members, i], color=col,
                       marker=mark, label=k)
            ax.scatter(cluster_center[0], cluster_center[i], marker=mark,
                       color=col, edgecolor='k', s=100)

        # Format plot
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        ax.set_xlabel(bands[0] + ' - ' + bands[1])
        ax.set_ylabel(bands[i*2]+' - '+bands[i*2+1])
        ax.set_title('mean-shift ' + str(n_clusters) + ' : '+bands[0]+'-'+bands[1]+' vs. '+bands[i*2]+'-'+bands[i*2+1],
                     fontsize=14)
        ax.legend(loc='lower right')

        '''Display interactive figure if # removed, if not, figures saved'''
        # plt.show
        file_name = 'meanshift_color_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                                bands[0], bands[1],
                                                                bands[i*2], bands[i*2+1])
        colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[i*2], bands[i*2+1])
        path_ = ('{}\\{}').format(path, colours)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pylab.savefig(os.path.join(path_, file_name))
        plt.close()
    return()


def kmeans_density(path, band1, band2, band3, band4, clf, scaler, colour1,
                       colour2, number_clusters, cluster_number):

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
    
    file_name = 'kmeans_dens_{}cl_{}-{}vs{}-{}.png'.format(str(number_clusters),
                                                           band1, band2, band3,
                                                           band4)
    colors = ('{}-{}_{}-{}').format(band1, band2, band3, band4)
    path_ = ('{}\\kmeans-dens_{}').format(path, colors)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    return ()


def kmeans_colour(path, cluster_data, number_clusters, cluster_number, bands,
                  cluster_centers):
    '''Plot xy positions of objects in different clusters'''

    
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    markers = cycle('ooooooo**')#*****<<<<<<>>>>>>>^^^^^^^')
    for i in range(1, 7):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for k, col, mark in zip(range(number_clusters), colors, markers):
            class_members = cluster_number == k
            cluster_center = cluster_centers[k]
            ax.scatter(cluster_data[class_members, 0], cluster_data[class_members, i],
                       color=col, marker=mark, label=k)
            ax.scatter(cluster_center[0], cluster_center[i], marker=mark, color=col,
                       edgecolor='k', s=100)

        ax.set_xlabel(bands[0]+' - '+bands[1])
        ax.set_ylabel(bands[i*2]+' - '+bands[i*2+1])
        ax.set_title('kmeans: '+bands[0]+'-'+bands[1]+' vs. '+bands[i*2]+'-'+bands[i*2+1],
                     fontsize=16)
        ax.legend(loc='lower right')

        file_name = 'kmeans_xy_{}cl_{}-{}vs{}-{}.png'.format(str(number_clusters),
                                                             bands[0], bands[1], bands[i*2],
                                                             bands[i*2+1])
        colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[i*2], bands[i*2+1])
        path_ = ('{}\\{}').format(path, colours, number_clusters)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pylab.savefig(os.path.join(path_, file_name))
        plt.close()
    return ()


def db_plot(u_labels, labels_, X, core_samples_mask_, n_clusters, band1, band2,
            band3, band4, s_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(u_labels)))
    for k, col in zip(u_labels, colors):
        if k == -1:
        # Black used for noise.
            col = 'k'

        class_member_mask = (labels_ == k)
        xy = X[class_member_mask & core_samples_mask_]
        ax.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=14)

        xy = X[class_member_mask & ~core_samples_mask_]
        ax.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=6)

    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_title('dbscan ' + str(n_clusters) +': '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                  fontsize=16)
    ax.set_xlabel(band1 + ' - ' + band2)
    ax.set_ylabel(band1 + ' - ' + band2)    
    file_name = 'dbscan_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                  band1, band2, band3, band4)
    pylab.savefig(os.path.join(s_path, file_name))
    plt.close()
       
    return()


def affinity_plot(labels, cluster_center_indices, X, n_clusters, bands,
                  s_path):
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    markers = cycle('ooooooo**')#*****<<<<<<>>>>>>>^^^^^^^')
    for i in range(0, 6):
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)    
        for k, col, mark in zip(range(n_clusters), colors, markers):
            class_members = labels == k
            cluster_center = X[cluster_center_indices[k]]
            ax.scatter(X[class_members, 0], X[class_members, i], color=col,
                       marker=mark, label=k)
            ax.scatter(cluster_center[0], cluster_center[i], marker=mark,
                       color=col, edgecolor='k', s=100)
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        ax.set_title('affinity ' + str(n_clusters) +': '+bands[0]+'-'+bands[1]+' vs. '+bands[i*2]+'-'+bands[i*2+1],
                     fontsize=16)
        ax.set_xlabel(bands[0] + ' - ' + bands[1])
        ax.set_ylabel(bands[i*2] + ' - ' + bands[i*2+1])  
        file_name = 'affinity_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                          bands[0], bands[1], bands[i*2],
                                                          bands[i*2+1])
        colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[i*2],
                                         bands[i*2+1])
        path_ = ('{}\\{}\\affinity').format(s_path, colours)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pylab.savefig(os.path.join(path_, file_name))
    return()


def meanshift_results(save_path, name, bands, n_clusters,
                      s_score, b_width, total_obj, obj_per_cluster):

    # Create MS results file if it doesn't exist. If it does add to it.
    test_path = '{}\\{}'.format(save_path, name)
    header = '#clustering band1 band2 band3 band4 b_width eps min_samp damp '\
             'pref km_in score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 ' \
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40 c_41 c_42 c_43 c_44 c_45 '\
             'c_46 c_47 c_48 c_49 c_50 c_51 c_52 c_53 c_54 c_55 c_56 c_57 c_58 '\
             'c_59 c_60'
    if not os.path.exists(test_path):
        create_path = os.path.join(save_path, name)
        ms_results_file = open(create_path, "a")
        ms_results_file.write(header + '\n')
        ms_results_file.close()
    ms_results_path = os.path.join(save_path, name)
    ms_results_file = open(ms_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {:.4f} {} {} {} {} {}'.format('meanshift', bands[0],
                                                           bands[1], bands[2], bands[3],
                                                           b_width, 'N/A',
                                                           'N/A', 'N/A', 'N/A',
                                                           'N/A')
    outputs = '{:.4f} {} {} {}'.format(s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    ms_results_file.write(inputs + ' ' + outputs + '\n')
    ms_results_file.close()
    return()


def kmeans_results(save_path, name, bands, input_,
                   n_clusters, s_score, total_obj, obj_per_cluster):
    # Create KM results file if it doesn't exist. If it does add to it.
    test_path = '{}\\{}'.format(save_path, name)
    header = '#clustering band1 band2 band3 band4 b_width eps min_samp damp '\
             'pref km_in score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 ' \
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40 c_41 c_42 c_43 c_44 c_45 '\
             'c_46 c_47 c_48 c_49 c_50 c_51 c_52 c_53 c_54 c_55 c_56 c_57 c_58 '\
             'c_59 c_60'
    if not os.path.exists(test_path):
        create_path = os.path.join(save_path, name)
        km_results_file = open(create_path, "a")
        km_results_file.write(header + '\n')
        km_results_file.close()
    km_results_path = os.path.join(save_path, name)
    km_results_file = open(km_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {} {} {} {} {} {}'.format('kmeans', bands[0],
                                                       bands[1], bands[2],
                                                       bands[3], 'N/A',
                                                       'N/A', 'N/A', 'N/A',
                                                       'N/A', input_)
    outputs = '{:.4f} {} {} {}'.format(s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    km_results_file.write(inputs + ' ' + outputs + '\n')
    km_results_file.close()
    return()


def dbscan_results(save_path, name, band1, band2, band3, band4, n_clusters,
                   s_score, eps, min_samp, total_obj, obj_per_cluster):
    # Create KM results file if it doesn't exist. If it does add to it.
    test_path = '{}\\{}'.format(save_path, name)
    header = '#clustering band1 band2 band3 band4 b_width eps min_samp damp '\
             'pref km_in score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 ' \
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40 c_41 c_42 c_43 c_44 c_45 '\
             'c_46 c_47 c_48 c_49 c_50 c_51 c_52 c_53 c_54 c_55 c_56 c_57 c_58 '\
             'c_59 c_60'
    if not os.path.exists(test_path):
        create_path = os.path.join(save_path, name)
        db_results_file = open(create_path, "a")
        db_results_file.write(header + '\n')
        db_results_file.close()
    db_results_path = os.path.join(save_path, name)
    db_results_file = open(db_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {} {} {} {} {} {}'.format('dbscan', band1, band2,
                                                       band3, band4, 'N/A',
                                                       eps, min_samp, 'N/A',
                                                       'N/A', 'N/A')
    outputs = '{:.4f} {} {} {}'.format(s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    db_results_file.write(inputs + ' ' + outputs + '\n')
    db_results_file.close()
    return()


def affinity_propagation_results(save_path, name, bands,
                                 n_clusters, s_score, damp, pref, total_obj,
                                 obj_per_cluster):
    # Create KM results file if it doesn't exist. If it does add to it.
    test_path = '{}\\{}'.format(save_path, name)
    header = '#clustering band1 band2 band3 band4 b_width eps min_samp damp '\
             'pref km_in score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 ' \
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40 c_41 c_42 c_43 c_44 c_45 '\
             'c_46 c_47 c_48 c_49 c_50 c_51 c_52 c_53 c_54 c_55 c_56 c_57 c_58 '\
             'c_59 c_60'
    if not os.path.exists(test_path):
        create_path = os.path.join(save_path, name)
        af_results_file = open(create_path, "a")
        af_results_file.write(header + '\n')
        af_results_file.close()
    af_results_path = os.path.join(save_path, name)
    af_results_file = open(af_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {} {} {} {} {} {}'.format('affinity', bands[0], bands[1],
                                                       bands[1], bands[3], 'N/A',
                                                       'N/A', 'N/A', damp,
                                                       pref,  'N/A')
    outputs = '{:.4f} {} {} {}'.format(s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    af_results_file.write(inputs + ' ' + outputs + '\n')
    af_results_file.close()
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
                        choices=['meanshift', 'kmeans', 'dbscan', 'affinity'],
                        nargs='*', default=[])

    # Optional Arguments: Parameter secification 
    inputs.add_argument("-kmi", "--kmeans_input", help="Choose the number of clusters input for kmeans", 
                        choices=['meanshift', 'dbscan', 'affinity',
                                 'experiments.txt'],
                        default = ['experiments.txt'])
    inputs.add_argument("-bwi", "--bandwidth_input", help="Choose bandwidth for meanshift custering",
                        choices=['experiments.txt', 'estimate'],
                        default=['estimate'])
    
    # Optional Arguments: Choose plots
    inputs.add_argument("-p","--plots", help = "Choose the plots you would like to make", 
                        choices=['meanshift_density','meanshift_colour',
                                 'kmeans_density', 'kmeans_colour', 'dbscan',
                                 'affinity'], nargs='*', default=[])
                        
    # Optional Arguments: Choose other functions to run (id, results_summary)
    inputs.add_argument("-id", "--id_list", help = "Produces object id list",
                        choices=['yes', 'no'],
                        default=[])

    criteria = inputs.parse_args()     
    
    clustering(criteria.plot_path, criteria.results_path, criteria.analysis,
               criteria.kmeans_input, criteria.bandwidth_input, criteria.plots,
               criteria.id_list, criteria.data_file)                                                         

    return()

if __name__ == "__main__": 
    user_input()
    