# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 13:24:41 2015
@author: Owner
"""
'''------------------------- Important Info -----------------------------------
MAKE SURE YOU CHANGE THE PATH IN MAKE_SAVE_DIRECTORY
Currently formatted for 2 colours.
Must be changed for different n_dimensions:
    - organize_data filters must be changed
    - headers in results files must be changed
        - New files for different dimesions
----------------------------------------------------------------------------'''

import os
import os.path
import argparse
import numpy as np
import pylab as pylab
from matplotlib import pyplot as plt
from astropy.table import Table
import shutil
from itertools import cycle
import time as time

# Kmeans imports
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn import metrics
from sklearn.metrics import silhouette_samples

# meanshift imports
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import preprocessing
from scipy.spatial import distance

# affinity propagation imports
from sklearn.cluster import AffinityPropagation

# Hierarchical imports
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph

'''-----------------------Global Variables----------------------------------'''
# used for plots
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',
          'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',
          'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*', '*',
           '.', '.', '.', '.', '.', '.', '.', '>', '>', '>', '>', '>', '>', '>',
           '+', '+', '+', '+', '+', '+', '+', '<', '<', '<', '<', '<', '<', '<',]

# need this so that output files always have the same number of columns
max_num_clusters = 40

# Set the base path and directory symbol for MAC or PC OS
# base_path = '/home/akiar/m83_clustering/'  # MAC
base_path = '/home/akiar/m83_clustering/'  # Sharcnet
# base_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\'  # PC
figure_save_symbol = '/'  # MAC Sharcnet
# figure_save_symbol = '\\'  # PC

# defined functions
from numpy import mean as avg
from numpy import square, sqrt
from numpy import std as stdev
from sklearn.metrics.pairwise import pairwise_distances as pdist
'''-------------------------------------------------------------------------'''


def clustering(save_plots, save_results, analysis, kmeans_input, bw_in, plots,
               id_list, data_file, write_res, ds9_cat, af_in, hc_input,
               input_file='experiments.txt'):
    '''DESCRIBE PROCESS HERE'''
    # Create saving directories
    plot_path, results_path = make_save_directory(save_plots, save_results)

    # Load data and trials
    data, experiments = load_data_file(data_file, input_file)

    '''Run Analysis'''
    for i in range(0, len(experiments)):
        '''-----------------Metrics reset every trial-----------------------'''
        ms_n_clusters = 0
        af_n_clusters = 0
        km_n_clusters = 0
        total_obj = 0
        num_obj = 0
        ms_score = 0
        af_score = 0
        km_score = 0
        bandwidth = 0
        results_title = ' '
        damping = 0
        preferences = 0
        n = 0
        st = time.time()
        '''-----------------------------------------------------------------'''

        # Get Data
        cluster_data_, greatdata, x_data, y_data, id_data, wave1, wave2, wave3, wave4 = \
            organize_data(experiments[i], data)

        # Set title of results.txt
        if "mag05" in experiments['band1'][i]:
            results_title = '05aperture_results_{}d.txt'.format(len(cluster_data_[0]))
        elif "mag3" in experiments['band1'][i]:
            results_title = '3aperture_results.txt'

        # Do clustering

        if "meanshift" in analysis:
            if 'experiments.txt' in bw_in:
                b_width_input = experiments['b_width'][i]
            else:
                b_width_input = estimate_bandwidth(cluster_data_)

            ms_n_clusters, bandwidth, ms_score, ms_obj, ms_obj_p_cluster, col = \
                meanshift(plot_path, experiments[i], cluster_data_,
                          plots, b_width_input, id_list, id_data, x_data,
                          y_data, ds9_cat, n, wave1, wave2, wave3, wave4)
            if "yes" in write_res:
                meanshift_results(results_path, results_title, experiments[i],
                                  ms_n_clusters, ms_score, bandwidth, ms_obj,
                                  ms_obj_p_cluster, col)

        if "hms" in analysis:
            if 'experiments.txt' in bw_in:
                b_width_input = experiments['b_width'][i]
            else:
                b_width_input = 'estimate'

            hms_n_clusters, bandwidth, hms_score, hms_obj, hms_obj_p_cluster =\
                hms(plot_path, experiments[i], cluster_data_, plots,
                    b_width_input, id_list, id_data, x_data, y_data, ds9_cat)

        if "affinity" in analysis:
            if 'experiments.txt' in af_in: 
                damping = float(experiments['damping'][i])
                preferences = float(experiments['preferences'][i])
            else: 
                damping = 0.95
                preferences = -len(cluster_data_)*0.15
            af_n_clusters, af_score, af_obj, af_obj_p_cluster, col = \
                affinity_propagation(plot_path, experiments[i], cluster_data_,
                                     plots, damping, preferences, id_list,
                                     id_data, x_data, y_data, ds9_cat, n,
                                     wave1, wave2, wave3, wave4)
            if "yes" in write_res:
                affinity_propagation_results(results_path, results_title,
                                             experiments[i], af_n_clusters,
                                             af_score, damping, preferences,
                                             af_obj, af_obj_p_cluster, col)
        
        if "hc" in analysis:
            if "experiments.txt" in hc_input:
                hc_n_clusters = int(experiments['n_clusters'][i])
            elif "meanshift" in kmeans_input:
                hc_n_clusters = ms_n_clusters
            elif "affinity" in kmeans_input:
                hc_n_clusters = af_n_clusters
            hc_score, num_obj, col = hierarchical(plot_path, experiments[i],
                                           cluster_data_, greatdata,
                                           hc_n_clusters, plots, id_list,
                                           x_data, y_data, id_data, ds9_cat,
                                           n)
            total_obj = num_obj.sum()
            if "yes" in write_res:
                hc_results(results_path, results_title, experiments[i],
                           hc_input, hc_n_clusters, hc_score, total_obj,
                           num_obj, col)

        if "kmeans" in analysis:
            if "experiments.txt" in kmeans_input:
                km_n_clusters = int(experiments['n_clusters'][i])
            elif "meanshift" in kmeans_input:
                km_n_clusters = ms_n_clusters
            elif "affinity" in kmeans_input:
                km_n_clusters = af_n_clusters
            if "experiments.txt" not in kmeans_input:
                if km_n_clusters >= 4:
                    high, low = 3, 2
                elif km_n_clusters == 3:
                    high, low = 4, 1
                else:
                    high, low = 5, 0
            else:
                high, low = 1, 0

            for a in range(km_n_clusters - low, km_n_clusters + high):
                km_score, num_obj, inertia, col = kmeans(plot_path, experiments[i],
                                           cluster_data_, greatdata,
                                           a, plots, id_list,
                                           x_data, y_data, id_data, ds9_cat,
                                           n, wave1, wave2, wave3, wave4)
                total_obj = num_obj.sum()
                if "yes" in write_res:
                    kmeans_results(results_path, results_title, experiments[i],
                                   kmeans_input, a, km_score, total_obj,
                                   num_obj, inertia, col)

        if "center_test" in analysis:
            n_clusters = experiments['n_clusters'][i]
            for n in range(1, 21):
                km_scor, num_obj, inertia, col = kmeans(plot_path, experiments[i],
                                          cluster_data_, greatdata,
                                          n_clusters, plots, id_list,
                                          x_data, y_data, id_data, ds9_cat,
                                          n, wave1, wave2, wave3, wave4)
                
        print "Finished Clustering"
        ed = time.time()
        print "Time: {}".format(ed - st)
    # Copy experiments.txt to plots directory
    shutil.copy2('experiments.txt',
                 plot_path + figure_save_symbol + 'inputs.txt')

    return()


def make_save_directory(p_path, r_path):
    '''Save results of each analysis
            save_path: set path of where you would like results saved'''
    # Create new analysis_folder
    res_path = '{}{}'.format(base_path, r_path)
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    pl_path = '{}{}'.format(base_path, p_path)
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

    # Change parameters to match data_file
    # Remove data pieces with no value
    wave1_trim = np.logical_and(np.logical_and(wave1 != -99, wave1_unc != -99),
                                wave1_unc < ratio)
    wave2_trim = np.logical_and(np.logical_and(wave2 != -99, wave2_unc != -99),
                                wave2_unc < ratio)

    wave3_trim = np.logical_and(np.logical_and(wave3 != -99, wave3_unc != -99),
                                wave3_unc < ratio)
    wave4_trim = np.logical_and(np.logical_and(wave4 != -99, wave4_unc != -99),
                                wave4_unc < ratio)

    colour1_trim = np.logical_and(wave1_trim, wave2_trim)
    colour2_trim = np.logical_and(wave3_trim, wave4_trim)

    final_data = np.logical_and(colour1_trim, colour2_trim)

    # Make colours
    colour1 = wave1[final_data] - wave2[final_data]
    colour2 = wave3[final_data] - wave4[final_data]
    wave1_fin = wave1[final_data]
    wave2_fin = wave2[final_data]
    wave3_fin = wave3[final_data]
    wave4_fin = wave4[final_data]

    cluster_data = np.vstack([colour1, colour2]).T 

    x = data['x'][final_data]
    y = data['y'][final_data]
    id_ = data['id_'][final_data].astype(np.int32)

    return (cluster_data, final_data, x, y, id_, wave1_fin, wave2_fin,
            wave3_fin, wave4_fin)


def meanshift(s_path, bands, cluster_data, make_plot, bw_input,
              output_id, id_data, x, y, ds9_cat, cent_test, w1, w2, w3, w4):
    '''---------------------------------------------------------------------'''
    '''Meanshift clustering to determine the number of clusters in the data,
    which can be passed to KMEANS function'''
    '''---------------------------------------------------------------------'''
    # TODO: SHOULD THE DATA BE SCALED????
    # X = cluster_data
    X_scaled = cluster_data  # preprocessing.scale(cluster_data)  

    # Compute clustering
    ms = MeanShift(bandwidth=bw_input, bin_seeding=False, cluster_all=True)
    ms.fit(X_scaled)

    # Cluster statistics
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    total_objects = len(X_scaled)

    # Compute silhouette_score for entire cluster or individuals
    if len(X_scaled) > n_clusters_ > 1:
        average_score = metrics.silhouette_score(X_scaled, labels)
        sample_score = metrics.silhouette_samples(X_scaled, labels)
    else:
        average_score = np.array(-99.0, dtype=float)
        sample_score = np.array(-99.0, dtype=float)
    
    # Set up save_path
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2], bands[3])
    
    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, n_clusters_):
        ith_cluster = X_scaled[labels == i]
        objects_per_cluster[i] = len(ith_cluster)
        if -99.0 in sample_score:
            write_cluster_stats(s_path, 'meanshift', n_clusters_, i,
                                objects_per_cluster[i], cluster_centers[i],
                                ith_cluster, average_score,
                                sample_score, cent_test, colours)
        else:
            write_cluster_stats(s_path, 'meanshift', n_clusters_, i,
                                objects_per_cluster[i], cluster_centers[i],
                                ith_cluster, average_score,
                                sample_score[labels == i], cent_test, colours)

    # Format for writing results
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    # Make plots
    if "meanshift_colour" in make_plot:
        meanshift_colour(s_path, X_scaled, n_clusters_, labels,
                         cluster_centers, bands, w1, w2, w3, w4)
    # Make catalgoues
    if "yes" in output_id:
        id_catologue('meanshift', n_clusters_, labels, bands, id_data, x, y,
                     s_path)
    if "yes" in ds9_cat:
        ds9_catalogue('meanshift', n_clusters_, labels, bands, x, y, s_path)

    return(n_clusters_, bw_input, average_score, total_objects,
           objects_per_cluster, colours)


def hms(s_path, bands, cluster_data, make_plot, bw_input, output_id, id_data,
        x, y, ds9_cat):
    # TODO: Figure out the best way to save and pass the results
            # Write function
    '''---------------------------------------------------------------------'''
    '''Perform hierarchical meanshift by creating incremental banwidth
    levels'''
    cluster_centers = []
    b_width = 0.0
    incriment = 3.0
    '''---------------------------------------------------------------------'''

    if "estimate" in bw_input:
        b_width = min(distance.pdist(cluster_data))

    while(len(cluster_centers) != 1):
        ms = MeanShift(bandwidth=b_width, bin_seeding=False, cluster_all=True)
        ms.fit(cluster_data)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters = len(labels_unique)
        if max_num_clusters > n_clusters > 1:
            hms_score = metrics.silhouette_score(cluster_data, labels)
            sample_score = metrics.silhouette_samples(cluster_data, labels)
        else:
            hms_score = -99.0
            sample_score = -99.0
        # Identify which cluster each object belongs to
        if n_clusters < max_num_clusters:
            objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
            for i in range(0, n_clusters):
                ith_cluster = cluster_data[labels == i]
                objects_per_cluster[i] = len(ith_cluster)
                hms_cluster_stats(s_path, 'hms', n_clusters, i,
                                  objects_per_cluster[i], cluster_centers[i],
                                  ith_cluster, hms_score,
                                  sample_score[labels == i], b_width)
                # Format for writing results
                objects_per_cluster.sort()  # sort from smallest to largest
                objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
            # Make Plots
            if "hms" in make_plot:
                hms_colour(s_path, cluster_data, n_clusters, labels,
                           cluster_centers, bands)
        # Make catalgoues
        if "yes" in output_id:
            id_catologue('hms', n_clusters, labels, bands, id_data, x, y,
                         s_path)
        if "yes" in ds9_cat:
            ds9_catalogue('hms', n_clusters, labels, bands, x, y, s_path)

        # Increase Bandwidth
        if b_width < 0.25:
            b_width = b_width*incriment
        else:
            b_width = b_width*1.1

    return(n_clusters, b_width, hms_score, len(cluster_data),
           objects_per_cluster)


def affinity_propagation(s_path, bands, cluster_data, make_plots, damp, pref,
                         output_id, id_data, x, y, ds9_cat, cent_test, w1, w2,
                         w3, w4):
    '''---------------------------------------------------------------------'''
    '''Perform Affinity Propagation clustering to determine the number of
    clusters in a given set of colours'''
    ap_score = 0
    '''---------------------------------------------------------------------'''

    # TODO: SHOULD THE DATA BE SCALED???
    # X_scaled = preprocessing.scale(cluster_data)

    # Compute similarities
    similarities = cluster_data
    # Compute Affinity Propagation
    af = AffinityPropagation(preference=pref, damping=damp).fit(similarities)

    # Calculate clustering statistics
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    n_clusters_ = len(cluster_centers_indices)
    total_objects = len(cluster_data)
    if n_clusters_ > 1:
        ap_score = metrics.silhouette_score(cluster_data, labels)
        sample_score = silhouette_samples(cluster_data, labels)
    else:
        ap_score = np.array(-99.0, dtype=float)
        sample_score = np.array(-99.0, dtype=float)

    # Set up save_path
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2], bands[3])

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, n_clusters_):
        ith_cluster = cluster_data[labels == i]
        objects_per_cluster[i] = len(ith_cluster)
        cluster_center = cluster_data[cluster_centers_indices[i]]
        if -99.0 in sample_score:
            write_cluster_stats(s_path, 'affinity', n_clusters_, i,
                                objects_per_cluster[i], cluster_center,
                                ith_cluster, ap_score,
                                sample_score, cent_test, colours)
        else:
            write_cluster_stats(s_path, 'affinity', n_clusters_, i,
                                objects_per_cluster[i], cluster_center,
                                ith_cluster, ap_score,
                                sample_score[labels == i], cent_test, colours)

    # Format for writing results
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort

    # Make plot
    if 'affinity' in make_plots:
        affinity_plot(labels, cluster_centers_indices, cluster_data,
                      n_clusters_, bands, s_path, w1, w2, w3, w4)

    # Write ID file
    if "yes" in output_id:
        id_catologue('affinity', n_clusters_, labels, bands, id_data, x, y,
                     s_path)
    if "yes" in ds9_cat:
        ds9_catalogue('affinity', n_clusters_, labels, bands, x, y, s_path)

    return(n_clusters_, ap_score, total_objects, objects_per_cluster, colours)


def kmeans(s_path, bands, cluster_data_, greatdata, number_clusters, make_plots,
           output_cluster_id, x, y, id_data, ds9_cat, cent_test, w1, w2, w3,
           w4):
    '''---------------------------------------------------------------------'''
    '''Perform K-means clustering on colours constructed from HST photometry
    using various combinations of filters'''
    score = 0
    '''---------------------------------------------------------------------'''
    # TODO: SHOULD THE DATA BE SCALED???
    # X_scaled = preprocessing.scale(cluster_data)
    # Data pre-processing
    cluster_data = cluster_data_ # preprocessing.scale(cluster_data_)  
    # Compute K-Means clustering
    km = KMeans(number_clusters, init='random')
    km.fit(cluster_data)

    # Compute clustering statistics
    sum_of_squares = km.inertia_
    centers = km.cluster_centers_
    labels = km.labels_
    if number_clusters > 1:
        score = metrics.silhouette_score(cluster_data, labels)
        sample_score = silhouette_samples(cluster_data, labels)
    else:
        score = np.array(-99.0, dtype=float)
        sample_score = np.array(-99.0, dtype=float)

    # Set up save_path
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2], bands[3])

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, number_clusters):
        x_cluster = cluster_data[labels == i]
        objects_per_cluster[i] = len(x_cluster)
        if -99.0 in sample_score:
            write_cluster_stats(s_path, 'kmeans', number_clusters, i,
                                objects_per_cluster[i], centers[i],
                                x_cluster, score,
                                sample_score, cent_test, colours)
        else:
            write_cluster_stats(s_path, 'kmeans', number_clusters, i,
                                objects_per_cluster[i], centers[i],
                                x_cluster, score,
                                sample_score[labels == i], cent_test, colours)

    # Format for writing results
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    # Generate plots
    if "kmeans_colour" in make_plots:
        kmeans_colour(s_path, cluster_data, number_clusters, labels,
                      bands, centers, cent_test, w1, w2, w3, w4)

    # Output object and cluster IDs to ID.txt file
    if "yes" in output_cluster_id:
        id_catologue('kmeans', number_clusters, labels, bands, id_data, x, y,
                     s_path)
    if "yes" in ds9_cat:
        ds9_catalogue('kmeans', number_clusters, labels, bands, x, y, s_path)
    
    return(score, objects_per_cluster, sum_of_squares, colours)


def hierarchical(s_path, bands, cluster_data_, greatdata, number_clusters, make_plots,
                 output_cluster_id, x, y, id_data, ds9_cat, cent_test):
    '''Perform Agglomerative Hierarchical Clustering'''
    # Define the structure A of the data. Here a 10 nearest neighbors
    
    connectivity = kneighbors_graph(cluster_data_, n_neighbors=20)

    # Compute clustering
    hc = AgglomerativeClustering(n_clusters=number_clusters, connectivity=connectivity,
                                 linkage='ward').fit(cluster_data_)
    
    # Compute statistics
    labels = hc.labels_
    if number_clusters > 1:
        score = metrics.silhouette_score(cluster_data_, labels)
        sample_score = silhouette_samples(cluster_data_, labels)
    else:
        score = np.array(-99.0, dtype=float)
        sample_score = np.array(-99.0, dtype=float)

    # Set up save_path
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2], bands[3])

    # Identify which cluster each object belongs to
    objects_per_cluster = np.zeros(max_num_clusters, dtype=np.int16)
    for i in range(0, number_clusters):
        x_cluster = cluster_data_[labels == i]
        objects_per_cluster[i] = len(x_cluster)
        if -99.0 in sample_score:
            write_cluster_stats(s_path, 'kmeans', number_clusters, i,
                                objects_per_cluster[i], 'N/A',
                                x_cluster, score,
                                sample_score, cent_test, colours)
        else:
            write_cluster_stats(s_path, 'kmeans', number_clusters, i,
                                objects_per_cluster[i], 'N/A',
                                x_cluster, score,
                                sample_score[labels == i], cent_test, colours)

    # Format for writing results
    objects_per_cluster.sort()  # sort from smallest to largest
    objects_per_cluster = objects_per_cluster[::-1]  # reverse sort
    # Generate plots
    if "h_colour" in make_plots:
        h_colour(s_path, cluster_data_, number_clusters, labels, bands)

    # Output object and cluster IDs to ID.txt file
    if "yes" in output_cluster_id:
        id_catologue('HC', number_clusters, labels, bands, id_data, x, y,
                     s_path)
    if "yes" in ds9_cat:
        ds9_catalogue('HC', number_clusters, labels, bands, x, y, s_path)
    
    return(score, objects_per_cluster, colours)


def id_catologue(clustering, number_clusters, cluster_number, waves, id_data,
                 x, y, save_):
    '''Create file with list of object ID and cluster number'''
    file_name = 'id_{}_{}cl_{}-{}vs{}-{}.txt'.format(clustering,
                                                     str(number_clusters),
                                                     waves[0], waves[1],
                                                     waves[2], waves[3])                                
    colours = ('{}-{}_{}-{}').format(waves[0], waves[1], waves[2], waves[3])
    path_ = ('{}{}{}{}{}').format(save_, figure_save_symbol, colours,
                                  figure_save_symbol, 'id_')
    if not os.path.exists(path_):
        os.makedirs(path_)

    id_path = os.path.join(path_, file_name)
    id_tab = Table(data=[id_data, x, y, cluster_number],
                   names=['object_id', 'x_pix', 'y_pix', 'cluster_number'])
    id_tab.write(id_path, format='ascii.commented_header')

    return()


def ds9_catalogue(clustering, n_clust, cluster_num, waves, x, y, save_):
    '''Create file with list of object x-y positions for each cluster'''
    colours = ('{}-{}_{}-{}').format(waves[0], waves[1], waves[2], waves[3])
    path_ = ('{}{}{}{}{}').format(save_, figure_save_symbol, colours,
                                  figure_save_symbol, 'ds9_' + clustering)
    if not os.path.exists(path_):
        os.makedirs(path_)

    ds_col = ['red', 'green', 'blue', 'cyan', 'magenta', 'black', 'white',
              'yellow']
    for i in range(0, n_clust):
        file_name = 'ds9_{}cl_cluster-{}_{}-{}vs{}-{}.reg'.format(str(n_clust),
                                                          str(i+1),
                                                          waves[0], waves[1],
                                                          waves[2], waves[3])
        x_coord = np.array(x[cluster_num == i])
        y_coord = np.array(y[cluster_num == i])
        ds9_path = os.path.join(path_, file_name)
        ds9_file = open(ds9_path, "w")
        for w in range(0, len(x_coord)): 
            coordinate_string = "{:.2f},{:.2f},".format(x_coord[w], y_coord[w])
            if i > 7:
                ds9_file.write("CIRCLE(" + coordinate_string + '15) # color = ' + ds_col[i-7] + '\n')
            else:
                ds9_file.write("CIRCLE(" + coordinate_string + '15) # color = ' + ds_col[i] + '\n')
        ds9_file.close()
        
    return()


def meanshift_colour(path, X, n_clusters, labels_, centers, bands, w1, w2, w3,
           w4):
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2], bands[3])
    path_ = ('{}{}{}').format(path, figure_save_symbol, colours)
    if not os.path.exists(path_):
        os.makedirs(path_)

    for i in range(1, 2):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for k in range(0, n_clusters):
            clust_col = plt.cm.jet(float(k) / np.max(labels_ + 1))
            my_members = labels_ == k
            cluster_center = centers[k]
            ax.scatter(X[my_members, 0], X[my_members, i], color=clust_col,
                       marker=markers[k], label=k, s=2, zorder=1)
            ax.scatter(cluster_center[0], cluster_center[i], marker=markers[k],
                       color=clust_col, edgecolor='k', s=100, zorder=2)

        # Format plot
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
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
        pylab.savefig(os.path.join(path_, file_name))
        plt.close()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    for c in range(0, n_clusters):
        my_members = labels_ == c
        clust_col = plt.cm.jet(float(c) / np.max(labels_ + 1))
        ax2.scatter(X[my_members, 1], w4[labels_ == c], color=clust_col, marker='.', s=4, label=c)
    ax2.legend(loc='lower right', fontsize=8)
    ax2.set_xlabel(bands[2]+' - '+bands[3])
    ax2.set_ylabel(bands[3])
    plt.gca().invert_yaxis()

    '''Display interactive figure if # removed, if not, figures saved'''
    file_name = 'meanshift_CMD_{}cl_{}-{}vs{}.png'.format(str(n_clusters),
                                                          bands[2], bands[3],
                                                          bands[3])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return()


def hms_colour(path, c_data, n_clusters, labels, centers, bands):

    for i in range(1, 2):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for k in range(0, n_clusters):
            my_members = labels == k
            cluster_center = centers[k]
            ax.scatter(c_data[my_members, 0], c_data[my_members, i],
                       color=colors[k], marker=markers[k], s=2, label=k,
                       zorder=1)
            ax.scatter(cluster_center[0], cluster_center[i], marker=markers[k],
                       color=colors[k], edgecolor='k', s=100, zorder=2)
        # Format plot
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        ax.set_xlabel(bands[0] + ' - ' + bands[1])
        ax.set_ylabel(bands[i*2]+' - '+bands[i*2+1])
        ax.set_title('mean-shift hierarchy' + str(n_clusters) + ' : '+bands[0]+'-'+bands[1]+' vs. '+bands[i*2]+'-'+bands[i*2+1],
                     fontsize=12)
        if n_clusters < 20:
            ax.legend(loc='lower right')
        '''Display interactive figure if # removed, if not, figures saved'''
        # plt.show
        file_name = 'HMeanShift_color_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                                bands[0], bands[1],
                                                                bands[i*2], bands[i*2+1])
        colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[i*2], bands[i*2+1])
        path_ = ('{}{}{}').format(path, figure_save_symbol, colours)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pylab.savefig(os.path.join(path_, file_name))
        plt.close()
    return


def kmeans_colour(path, cluster_data, number_clusters, cluster_number, bands,
                  cluster_centers, cent_test, w1, w2, w3, w4):
    '''Plot colour-colour diagrams for each colour space'''
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2], bands[3])
    if cent_test == 0:
        file_name = 'kmeans_col_{}cl_{}-{}vs{}-{}.png'.format(str(number_clusters),
                                                                  bands[0], bands[1], bands[2],
                                                                  bands[3])
        path_ = ('{}{}{}').format(path, figure_save_symbol, colours)
    else:
        file_name = 'kmeans_{}_col_{}cl_{}-{}vs{}-{}.png'.format(str(cent_test),
                                                                  str(number_clusters),
                                                                  bands[0], bands[1], bands[2],
                                                                  bands[3])
        path_ = ('{}{}{}{}{}').format(path, figure_save_symbol, colours,
                                      figure_save_symbol, 'cent_test_'+str(number_clusters))
    if not os.path.exists(path_):
        os.makedirs(path_)

    for i in range(1, 2):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for k in range(0, number_clusters):
            class_members = cluster_number == k
            cluster_center = cluster_centers[k]
            ax.scatter(cluster_data[class_members, 0],
                       cluster_data[class_members, i], color=colors[k],
                       marker=markers[k], label=k, s=2, zorder=1)
            ax.scatter(cluster_center[0], cluster_center[i], marker=markers[k],
                       color=colors[k], edgecolor='k', s=100, zorder=2)

        ax.set_xlabel(bands[0]+' - '+bands[1])
        ax.set_ylabel(bands[i*2]+' - '+bands[i*2+1])
        ax.set_title('kmeans: '+bands[0]+'-'+bands[1]+' vs. '+bands[i*2]+'-'+bands[i*2+1],
                     fontsize=16)
        ax.legend(loc='lower right')

        pylab.savefig(os.path.join(path_, file_name))
        plt.close()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    for c in range(0, number_clusters):
        my_members = cluster_number == c
        clust_col = plt.cm.jet(float(c) / np.max(cluster_number + 1))
        ax2.scatter(cluster_data[my_members, 1], w4[cluster_number == c],
                    color=clust_col, marker='.', s=4, label=c)
    ax2.legend(loc='lower right', fontsize=8)
    ax2.set_xlabel(bands[2]+' - '+bands[3])
    ax2.set_ylabel(bands[3])
    plt.gca().invert_yaxis()

    '''Display interactive figure if # removed, if not, figures saved'''
    file_name = 'kmeans_CMD_{}cl_{}-{}vs{}.png'.format(str(number_clusters),
                                                          bands[2], bands[3],
                                                          bands[3])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    return ()


def affinity_plot(labels, cluster_center_indices, X, n_clusters, bands,
                  s_path, w1, w2, w3, w4):
    for i in range(1, 2):
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        for k in range(0, n_clusters):
            class_members = labels == k
            cluster_center = X[cluster_center_indices[k]]
            ax.scatter(X[class_members, 0], X[class_members, i],
                       color=colors[k], marker=markers[k], label=k, s=2, zorder=1)
            ax.scatter(cluster_center[0], cluster_center[i], marker=markers[k],
                       color=colors[k], edgecolor='k', s=100, zorder=2)
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        ax.set_title('affinity ' + str(n_clusters) +': '+bands[0]+'-'+bands[1]+' vs. '+bands[i*2]+'-'+bands[i*2+1],
                     fontsize=16)
        ax.set_xlabel(bands[0] + ' - ' + bands[1])
        ax.set_ylabel(bands[i*2] + ' - ' + bands[i*2+1])
        ax.legend(loc='lower right')
        file_name = 'affinity_{}cl_{}-{}vs{}-{}.png'.format(str(n_clusters),
                                                            bands[0], bands[1],
                                                            bands[i*2],
                                                            bands[i*2+1])
        colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[i*2],
                                         bands[i*2+1])
        path_ = ('{}{}{}').format(s_path, figure_save_symbol, colours)
        if not os.path.exists(path_):
            os.makedirs(path_)
        pylab.savefig(os.path.join(path_, file_name))
    return()


def meanshift_results(save_path, name, bands, n_clusters,
                      s_score, b_width, total_obj, obj_per_cluster, colours):

    # Create MS results file if it doesn't exist. If it does add to it.
    colour_path = '{}{}{}'.format(save_path, figure_save_symbol, colours)
    if not os.path.exists(colour_path):
        os.makedirs(colour_path)
    test_path = '{}{}{}{}{}'.format(save_path, figure_save_symbol, colours,
                                    figure_save_symbol, name)
    header = '#clustering band1 band2 band3 band4 '\
             'b_width damp '\
             'pref km_in inertia score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 '\
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40'
    if not os.path.exists(test_path):
        create_path = os.path.join(colour_path, name)
        ms_results_file = open(create_path, "a")
        ms_results_file.write(header + '\n')
        ms_results_file.close()
    ms_results_path = os.path.join(colour_path, name)
    ms_results_file = open(ms_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {:.4f} {} {} {}'.format('meanshift', bands[0], 
                                                     bands[1], bands[2],
                                                     bands[3], # bands[4],
                                                     # bands[5],  # bands[6],
                                                       # bands[7], bands[8],
                                                     # bands[9], bands[10],
                                                     # bands[11]
                                                     b_width,
                                                     'N/A',
                                                     'N/A', 'N/A')
    outputs = '{} {} {} {} {}'.format('N/A', s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    ms_results_file.write(inputs + ' ' + outputs + '\n')
    ms_results_file.close()
    return()


def kmeans_results(save_path, name, bands, input_, n_clusters,
                   s_score, total_obj, obj_per_cluster, sos, colours):
    # Create KM results file if it doesn't exist. If it does add to it.
    colour_path = '{}{}{}'.format(save_path, figure_save_symbol, colours)
    if not os.path.exists(colour_path):
        os.makedirs(colour_path)
    test_path = '{}{}{}{}{}'.format(save_path, figure_save_symbol, colours,
                                    figure_save_symbol, name)
    header = '#clustering band1 band2 band3 band4 '\
             'b_width damp '\
             'pref km_in inertia score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 '\
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40'
    if not os.path.exists(test_path):
        create_path = os.path.join(colour_path, name)
        km_results_file = open(create_path, "a")
        km_results_file.write(header + '\n')
        km_results_file.close()
    km_results_path = os.path.join(colour_path, name)
    km_results_file = open(km_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {} {} {} {}'.format('kmeans', bands[0],
                                                 bands[1], bands[2],
                                                 bands[3], # bands[4],
                                                 # bands[5], #, bands[6],
                                                 # bands[7], bands[8],
                                                 # bands[9], bands[10],
                                                 # bands[11],
                                                 'N/A', 'N/A',
                                                 'N/A', input_)
    outputs = '{} {} {} {} {}'.format(str(sos), s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    km_results_file.write(inputs + ' ' + outputs + '\n')
    km_results_file.close()
    return()


def affinity_propagation_results(save_path, name, bands,
                                 n_clusters, s_score, damp, pref, total_obj,
                                 obj_per_cluster, colours):
    # Create KM results file if it doesn't exist. If it does add to it.
    colour_path = '{}{}{}'.format(save_path, figure_save_symbol, colours)
    if not os.path.exists(colour_path):
        os.makedirs(colour_path)
    test_path = '{}{}{}{}{}'.format(save_path, figure_save_symbol, colours,
                                    figure_save_symbol, name)
    header = '#clustering band1 band2 band3 band4 '\
             'b_width damp '\
             'pref km_in inertia score n_clust total_objects c_1 c_2 c_3 c_4 c_5 c_6 '\
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32 '\
             'c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40'
    if not os.path.exists(test_path):
        create_path = os.path.join(colour_path, name)
        af_results_file = open(create_path, "a")
        af_results_file.write(header + '\n')
        af_results_file.close()
    af_results_path = os.path.join(colour_path, name)
    af_results_file = open(af_results_path, "a")

    # Create strings for file
    inputs = '{} {} {} {} {} {} {} {} {}'.format('affinity', bands[0], bands[1],
                                                 bands[1], bands[3], # bands[4],
                                                 # bands[5],
                                                 # bands[6],
                                                 # bands[7], bands[8],
                                                 # bands[9], bands[10],
                                                 # bands[11],
                                                 'N/A', damp, pref,  'N/A')
    outputs = '{} {} {} {} {}'.format('N/A', s_score, n_clusters, total_obj,
                                       np.array_str(obj_per_cluster,
                                                    max_line_width=500)[1:-1])

    # Write results
    af_results_file.write(inputs + ' ' + outputs + '\n')
    af_results_file.close()
    return()


def write_cluster_stats(save_path, clustering, n_clust, cluster, n_obj,
                        center, c_data, t_score, c_score, cen_test, cols):
    header = '# clustering total_clust clust_num n_obj t_scr c_scr rms avg_dist '\
             'max_dist min_dist stdev cen_1 cen_2' \
             ' avg_col_1 avg_col_2'
    if cen_test == 0:
        name = 'cluster_statistics.txt'
    else:
        name = '{}cl_cent_test_statistics.txt'.format(n_clust)

    colour_path = '{}{}{}'.format(save_path, figure_save_symbol, cols)
    if not os.path.exists(colour_path):
        os.makedirs(colour_path)
    test_path = '{}{}{}{}{}'.format(save_path, figure_save_symbol, cols,
                                    figure_save_symbol, name)

    if not os.path.exists(test_path):
        create_path = os.path.join(colour_path, name)
        cluster_statistics = open(create_path, "a")
        cluster_statistics.write(header + '\n')
        cluster_statistics.close()
    cluster_statistics_path = os.path.join(colour_path, name)
    cluster_statistics = open(cluster_statistics_path, "a")
    
    # Calculate distance metrics
    if len(c_data) > 1:
        rms = sqrt(avg(square(c_data)))
        c_dist = pdist(c_data)
        cluster_dist = c_dist[c_dist != 0.0]
        max_dist = max(cluster_dist)
        min_dist = min(cluster_dist)
        mean_dist = avg(cluster_dist)
        distances = '{:.4f} {:.4f}  {:.4f} {:.4f} {:.4f}'.format(float(rms),
                                                         float(mean_dist),
                                                         float(max_dist),
                                                         float(min_dist),
                                                         float(stdev(c_data)))
    else:
        rms = 'N/A'
        max_dist = 'N/A'
        min_dist = 'N/A'
        mean_dist = 'N/A'
        distances = '{} {} {} {} {:.4f}'.format(rms, (mean_dist),
                                             (max_dist),
                                             (min_dist),
                                             float(stdev(c_data)))
    # Write stats
    inputs = '{} {} {} {} {} {}'.format(clustering, n_clust,
                                                cluster+1, n_obj, t_score,
                                                avg(c_score))

    centers = '{} {} '.format(center[0], center[1])
                                         # center[3], center[4], center[5])

    a = '{:.4f} {:.4f}'.format(avg(c_data[:,0]), avg(c_data[:,1]))
                               # avg(c_data[:,2]))  #, avg(c_data[:,3]),
                                      # avg(c_data[:,4]), avg(c_data[:,5]))
    cluster_statistics.write(inputs + ' ' + distances + ' ' + centers + ' ' + a + '\n')
    cluster_statistics.close()
    return


def hms_cluster_stats(save_path, clustering, n_clust, cluster, n_obj,
                      center, c_data, t_score, c_score, bandwidth):
    header = '# clustering bw total_clust clust_num n_obj t_scr c_scr avg_dist '\
             'max_dist min_dist stdev cen_1 cen_2 '\
             ' avg_col_1 avg_col_2 '
    name = 'hms_statistics.txt'
    test_path = '{}{}{}'.format(save_path, figure_save_symbol, name)
    if not os.path.exists(test_path):
        create_path = os.path.join(save_path, name)
        cluster_statistics = open(create_path, "a")
        cluster_statistics.write(header + '\n')
        cluster_statistics.close()
    cluster_statistics_path = os.path.join(save_path, name)
    cluster_statistics = open(cluster_statistics_path, "a")

    # Calculate distance metrics
    if len(c_data) > 1:
        rms = sqrt(avg(square(c_data)))
        c_dist = pdist(c_data)
        cluster_dist = c_dist[c_dist != 0.0]
        max_dist = max(cluster_dist)
        min_dist = min(cluster_dist)
        mean_dist = avg(cluster_dist)
        distances = '{:.4f} {:.4f}  {:.4f} {:.4f} {:.4f}'.format(float(rms),
                                                         float(mean_dist),
                                                         float(max_dist),
                                                         float(min_dist),
                                                         float(stdev(c_data)))
    else:
        max_dist = 'N/A'
        min_dist = 'N/A'
        mean_dist = 'N/A'
        rms = 'N/A'
        distances = '{} {} {} {} {:.4f}'.format(rms, mean_dist,
                                                max_dist,
                                                min_dist,
                                                float(stdev(c_data)))
    
    # Write stats
    inputs = '{} {:.4f} {} {} {} {:.4f} {:.4f}'.format(clustering, bandwidth,
                                                       n_clust, cluster+1,
                                                       n_obj, t_score,
                                                       avg(c_score))

    centers = '{:.4f} {:.4f} '.format(center[0], center[1])
    a = '{:.4f} {:.4f}'.format(avg(c_data[:,0]), avg(c_data[:,1]))
    cluster_statistics.write(inputs + ' ' + distances + ' ' + centers + ' ' + a + '\n')
    cluster_statistics.close()
    return


def user_input():

    inputs = argparse.ArgumentParser(description="Inputs for analysis.")

    # Mandatory arguments: data_file
    inputs.add_argument("data_file", help="Choose the data file for analysis")

    # Optional Arguments: Save paths for results text files and figure files
    inputs.add_argument("-pp", "--plot_path",
                        help="Enter path for saving plot files",
                        default="results")
    inputs.add_argument("-rp", "--results_path",
                        help="Enter path for saving results files",
                        default="results")

    # Optional Arguments: Choose the clustering method and inputs for kmeans
    inputs.add_argument("-a", "--analysis",
                        help="Choose the methods you would like to use",
                        choices=['meanshift', 'hms', 'kmeans', 'affinity',
                        'center_test', 'hc'], nargs='*', default=[])

    # Optional Arguments: Parameter secification
    inputs.add_argument("-kmi", "--kmeans_input",
                        help="Choose the number of clusters input for kmeans",
                        choices=['meanshift', 'affinity', 'experiments.txt'],
                        default=['experiments.txt'])
    inputs.add_argument("-hci", "--hierarchical_input",
                        help="Choose the number of clusters input for hc",
                        choices=['meanshift', 'affinity', 'experiments.txt'],
                        default=['experiments.txt'])
    inputs.add_argument("-bwi", "--bandwidth_input",
                        help="Choose bandwidth for meanshift custering",
                        choices=['experiments.txt', 'estimate'],
                        default=['estimate'])
    inputs.add_argument("-afp", "--affinity_input",
                        help="Choose parameters for affinity propagation",
                        choices=['experimments.txt', 'estimate'],
                        default=['estimate'])

    # Optional Arguments: Choose plots
    inputs.add_argument("-p", "--plots",
                        help="Choose the plots you would like to make",
                        choices=['meanshift_density', 'meanshift_colour',
                                 'kmeans_density', 'kmeans_colour',
                                 'affinity', 'hms', 'hc'], nargs='*',
                                 default=[])
    inputs.add_argument("-wr", "--write_results", help="Write main results",
                        choices=['yes', 'no'], default=['yes'])

    # Optional Arguments: Catalogue functions
    inputs.add_argument("-id", "--id_list", help="Produces object id list",
                        choices=['yes', 'no'], default=['no'])
    inputs.add_argument("-ds", "--ds_cat", help="Produces ds9 catalogue",
                        choices=['yes', 'no'], default=['no'])

    criteria = inputs.parse_args()

    clustering(criteria.plot_path, criteria.results_path, criteria.analysis,
               criteria.kmeans_input, criteria.bandwidth_input, criteria.plots,
               criteria.id_list, criteria.data_file, criteria.write_results,
               criteria.ds_cat, criteria.affinity_input,
               criteria.hierarchical_input)

    return()

if __name__ == "__main__":
    user_input()
