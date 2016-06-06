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
from astropy.table import Table, Column
import shutil
from itertools import cycle

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

    # User criteria
    analysis_criteria = analysis
    kmeans_input = kmeans_input
    generate_plots = plots
    id_output = id_list
    b_width_input = bw_in
    
    # Load data and trials
    data, experiments = load_data_file(data_file, input_file)

    # Import experiments listed in experiments.txt file
    #run = np.genfromtxt(input_file, dtype='str')

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
        b1 = experiments['band1'][i]
        b2 = experiments['band2'][i]
        b3 = experiments['band3'][i]
        b4 = experiments['band4'][i]
        '''-----------------------------------------------------------------'''

        # Get Data
        colour1, colour2, greatdata, x_data, y_data, id_data = \
            organize_data(experiments[i], data)

        # Set title of results.txt
        if "mag05" in experiments['band1'][i]:
            results_title = '05aperture_results.txt'
        elif "mag3" in experiments['band1'][i]:
            results_title = '3aperture_results.txt'

        # Do clustering

        if "meanshift" in analysis_criteria:
            if 'experiments.txt' in b_width_input:
                b_width_input = experiments['bandwidth'][i]

            ms_n_clusters, bandwidth, ms_score, ms_obj, ms_obj_p_cluster = \
                do_meanshift(plot_path, b1, b2, b3, b4, colour1, colour2,
                             generate_plots, b_width_input, id_list, id_data)

            meanshift_results(results_path, results_title, b1, b2, b3, b4,
                              ms_n_clusters, ms_score, bandwidth, ms_obj,
                              ms_obj_p_cluster)

        if "dbscan" in analysis_criteria: 
            eps = float(experiments['eps'][i])
            min_samples = float(experiments['min_samples'][i])
            db_score, db_n_clusters, db_obj, db_obj_p_cluster = \
                dbscan(plot_path, b1, b2, b3, b4, colour1, colour2,
                       generate_plots, eps, min_samples)
            dbscan_results(results_path, results_title, b1, b2, b3, b4, 
                           db_n_clusters, db_score, eps, min_samples, db_obj,
                           db_obj_p_cluster)

        if "affinity" in analysis_criteria:
            damping = float(experiments['damping'][i])
            preferences = float(experiments['preferences'][i])
            af_n_clusters, af_score, af_obj, af_obj_p_cluster = \
                affinity_propagation(plot_path, b1, b2, b3, b4, colour1,
                                     colour2, generate_plots, damping,
                                     preferences)
            affinity_propagation_results(results_path, results_title, b1, b2,
                                         b3, b4, af_n_clusters, af_score,
                                         damping, preferences, af_obj,
                                         af_obj_p_cluster)

        if "kmeans" in analysis_criteria:
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
                km_score, num_obj = do_kmeans(plot_path, b1, b2, b3, b4, colour1,
                                              colour2, greatdata,
                                              a, generate_plots, id_output,
                                              x_data, y_data, id_data)
                total_obj = num_obj.sum()
                kmeans_results(results_path, results_title, b1, b2, b3, b4,
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


def organize_data(band_combinations, data_file):
    '''Select data for analysis'''
    data = data_file
    ratio = 0.05
    # Colour 1
    wave1 = data[band_combinations['band1']]
    wave1_unc = data[band_combinations['band1']+'_unc']
    wave2 = data[band_combinations['band2']]
    wave2_unc = data[band_combinations['band2']+'_unc']
    # Colour 2
    wave3 = data[band_combinations['band3']]
    wave3_unc = data[band_combinations['band3']+'_unc']
    wave4 = data[band_combinations['band4']]
    wave4_unc = data[band_combinations['band4']+'_unc']

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

    # Remove data above given noise/signal ratio
    gooddata2 = np.logical_and(colour1_ratio, colour2_ratio)

    # Only data that match criteria for both colours
    greatdata = np.logical_and(gooddata1, gooddata2)
    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]

    x = data['x'][greatdata]
    y = data['y'][greatdata]
    id_ = data['id'][greatdata].astype(np.int32)

    return (colour1, colour2, greatdata, x, y, id_)


def do_meanshift(s_path, band1, band2, band3, band4, colour1, colour2,
                 make_plot, bw_input, output_id, id_data):
    '''Meanshift clustering to determine the number of clusters in the
        data, which is passed to KMEANS function'''
    bandwidth = 0
    
    # Truncate data
    X = np.vstack([colour1, colour2]).T
    
    if 'experiments.txt' not in bw_input:
        # The following bandwidth can be automatically detected using
        bandwidth = estimate_bandwidth(X) #, quantile=0.1, n_samples=len(colour1))
    else:
        bandwidth = bw_input

    X_scaled = preprocessing.scale(X)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=True)
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
        average_score = metrics.silhouette_score(X, labels)
    else: 
        average_score = -99
    #sample_score = metrics.silhouette_samples(X, labels)

    if "yes" in output_id:
        id_catologue(n_clusters_, labels, band1, band2, band3,
                     band4, id_data, s_path)
    # Make plot
    if "meanshift_density" in make_plot:
        meanshift_density(s_path, colour1, colour2, n_clusters_, X, ms,
                          band1, band2, band3, band4, labels, cluster_centers)
    if "meanshift_colour" in make_plot: 
        meanshift_colour(s_path, X_scaled, n_clusters_, labels, cluster_centers,
                     band1, band2, band3, band4, )

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


def affinity_propagation(s_path, band1, band2, band3, band4, colour1, colour2,
                         make_plots, damp, pref, output_id, id_data):

    X = np.vstack([colour1, colour2]).T
    # Compute Affinity Propagation
    af = AffinityPropagation(preference=pref, damping=damp).fit(X)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    n_clusters_ = len(cluster_centers_indices)

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, n_clusters_):
        ith_cluster = X[labels == i]
        objects_per_cluster[i] = len(ith_cluster)
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    total_objects = len(X)

    # Silhouette Score
    ap_score = metrics.silhouette_score(X, labels)
    if "yes" in output_id:
        id_catologue(n_clusters_, labels, band1, band2, band3,
                     band4, id_data, s_path)

    if 'affinity' in make_plots:
        affinity_plot(labels, cluster_centers_indices, X, n_clusters_, band1,
                      band2, band3, band4, s_path)

    return(n_clusters_, ap_score, total_objects, objects_per_cluster)


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
    clf.fit(scaler.fit_transform(clusterdata)) #change variable - data with labels
    #clf.fit(clusterdata)
    cluster_number = clf.predict(scaler.fit_transform(clusterdata))
    #cluster_number = clf.predict(clusterdata)

    # Output object and cluster IDs to ID.txt file
    if "yes" in output_cluster_id:
        id_catologue(number_clusters, cluster_number, band1, band2, band3,
                     band4, id_data, s_path)

    # Compute the silhouette score
    labels = clf.labels_
    score = 0
    if number_clusters > 1:
        score = metrics.silhouette_score(scaler.fit_transform(clusterdata), labels)
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
        kmeans_density(s_path, band1, band2, band3, band4, clf, scaler,
                       colour1, colour2, number_clusters, cluster_number)
    if "kmeans_colour" in make_plots:
        kmeans_colour(s_path, colour1, colour2, number_clusters, cluster_number,
                band1, band2, band3, band4)

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


def meanshift_density(path, colour1, colour2, n_clusters, X, ms,
                      band1, band2, band3, band4, labels_, centers):
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

def meanshift_colour(path, X, n_clusters, labels_, centers, band1, band2, band3,
                 band4): 
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

    for k, col in zip(range(n_clusters), colors):
        #cluster_score = sample_score[labels == k]
        my_members = labels_ == k
        cluster_center = centers[k]
        ax.plot(X[my_members, 0], X[my_members, 1], col + '.', label = k,
                markersize=8)
        #print ("objects in cluster {}: {}").format(k+1, len(X[my_members]))
        #print("    Cluster {} Average Score: {:.4f}").format(k+1,
              #np.average(cluster_score))
        ax.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=10)
    
    # Format plot
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(band1 + ' - ' + band2)
    ax.set_ylabel(band3+' - '+band4)
    ax.set_title('mean-shift ' + str(n_clusters) + ' : '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                 fontsize=14)
    ax.legend(loc='upper left')

    '''Display interactive figure if # removed, if not, figures saved'''
    # plt.show
    file_name = 'meanshift_color_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                            band1, band2,
                                                            band3, band4)
    colors = ('{}-{}_{}-{}').format(band1, band2, band3, band4)
    path_ = ('{}\\meanshift-{}').format(path, colors)
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


def kmeans_colour(path, x, y, number_clusters, cluster_number, band1, band2, band3,
            band4):
    '''Plot xy positions of objects in different clusters'''
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    colors = ('{}-{}_{}-{}').format(band1, band2, band3, band4)
    # Plot XY coordinates of each object and label with colour that
    # identifies cluster number
    #ax2.set_xlim(min(x), max(x))
    #ax2.set_ylim(min(y), max(y))
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        y_cluster = y[cluster_number == i]
        ax2.plot(x_cluster, y_cluster, cluster_colours[i] + '.', label=i)

    ax2.set_xlabel(band1+' - '+band2)
    ax2.set_ylabel(band3+' - '+band4)
    ax2.set_title('kmeans: '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                  fontsize=16)
    ax2.legend(loc='upper left')

    file_name = 'kmeans_xy_{}cl_{}-{}vs{}-{}.png'.format(str(number_clusters),
                                                         band1, band2, band3,
                                                         band4)
    path_ = ('{}\\kmeans-{}').format(path, colors)
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


def affinity_plot(labels, cluster_center_indices, X, n_clusters, band1, band2,
                  band3, band4, s_path):
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)

    colors = cycle('bgrcmykwbgrcmykwbgrcmykwbgrcmykw')
    for k, col in zip(range(n_clusters), colors):
        class_members = labels == k
        cluster_center = X[cluster_center_indices[k]]
        ax.plot(X[class_members, 0], X[class_members, 1], col + '.')
        ax.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=14)
        #for x in X[class_members]:
        #   plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_title('affinity ' + str(n_clusters) +': '+band1+'-'+band2+' vs. '+band3+'-'+band4,
                  fontsize=16)
    ax.set_xlabel(band1 + ' - ' + band2)
    ax.set_ylabel(band1 + ' - ' + band2)  
    file_name = 'affinity_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                      band1, band2, band3,
                                                      band4)
    colors = ('{}-{}_{}-{}').format(band1, band2, band3, band4)
    path_ = ('{}\\affinity-{}').format(s_path, colors)
    if not os.path.exists(path_):
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    return()


def meanshift_results(save_path, name, band1, band2, band3, band4, n_clusters,
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
    inputs = '{} {} {} {} {} {:.4f} {} {} {} {} {}'.format('meanshift', band1,
                                                           band2, band3, band4,
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


def kmeans_results(save_path, name, band1, band2, band3, band4, input_,
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
    inputs = '{} {} {} {} {} {} {} {} {} {} {}'.format('kmeans', band1, band2,
                                                       band3, band4, 'N/A',
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


def affinity_propagation_results(save_path, name, band1, band2, band3, band4,
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
    inputs = '{} {} {} {} {} {} {} {} {} {} {}'.format('affinity', band1, band2,
                                                       band3, band4, 'N/A',
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
    