# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 13:24:41 2015

@author: Owner
"""
import numpy as np 
import pylab as pylab
from matplotlib import pyplot as plt 

from astropy.table import Table, Column


from sklearn.cluster import KMeans
from sklearn import preprocessing 
from sklearn import metrics 
from sklearn.metrics import pairwise_distances


from matplotlib.patches import Ellipse
from scipy.stats import norm

from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import preprocessing

# correspondance between photometric measurement and column name in data file
# (needed because we can't get ascii.commented_header reader to work)
band_names = {'05_225': 11, '3_225': 13, '05_336': 15,  '3_336': 17, '05_373': 19, '3_373': 21, '05_438':23,    '3_438':25,  '05_487':27, '3_487': 29, '05_502':31, '3_502':33, '05_555': 35, '3_555': 37, '05_657': 39 ,'3_657':41 ,'05_673':43, '3_673':45 , '05_814':47 , '3_814':49 }

# used for plots
cluster_colours = ['y','g','b','r','c','m','k','m','w','y','g','b','r','c','m','k','m','w','y','g','b','r','c','m','k','m','w'] 

# need this so that output files always have the same number of columns
max_num_clusters = 20

def do_everything(input_file = 'experiments.txt', output_file = 'results.txt', mp=True, oci=False, rs = True):
    '''Automate clustering process
       input: input_file:  a 5-column text file with 1 line per clustering run
                           each line lists the 4 filters to be used to construct colours, plus number of clusters
              mp: make output plots
              oci: output cluster IDs for each object
       output: output_file: a text file listing input+results from each clustering run'''
    
    run = np.genfromtxt(input_file, dtype='str')

    # TODO: check whether results file already exists; if not, open it and print a header line
    # if it does already exist, just open it
    results = open(output_file, 'a') 
    
    for i in range(0, len(run)):
        
        
        #New MEANSHIFT
        numberofclusters = do_meanshift (run[i,0], run[i,1], run[i,2], run[i,3], mp)
        print "Estimated number of clusters: ", numberofclusters 
        
        input_str =  '{} {}'.format(np.array_str(run[i][:-1])[1:-1],numberofclusters) # list of input parameters: bands and num of clusters
        
               
        
        score, num_obj =  do_kmeans(run[i,0], run[i,1], run[i,2], run[i,3], numberofclusters, make_plots=mp, output_cluster_id=oci)
        total_obj = num_obj.sum()
        output_str = ' {:.4f} {:5d} {}'.format(score, total_obj, np.array_str(num_obj)[1:-1])
        
        results.write(input_str + ' ' + output_str + '\n')
        
       

    results.close()
    
    if rs: 
        results_summary()
    
    return

def do_meanshift (band1, band2, band3, band4, make_plots):
    
    #----------------------------------------------------------------------
    if band1 == band2 or band3 == band4: 
        print "Not a good idea to use the same band in one colour, try again"
        return
    for band in [band1, band2, band3, band4]:
        if band not in band_names.keys():
            print "Can't find %s in band_name list" %band
            return
    #------------------------------------------------------------
    # Get the data
    data = np.loadtxt('hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1.txt')
    
    #Import 4 different wavelengths
    #Colour 1: 05_mag
    wave1 = data[:, band_names[band1]]
    wave2 = data[:, band_names[band2]]
    
    #Colour 2: 05_mag
    wave3 = data[:, band_names[band3]]
    wave4 = data[:, band_names[band4]]
    
    gooddata1 = np.logical_and(np.logical_and(wave1!=-99, wave2!=-99), np.logical_and(wave3!=-99, wave4!=-99)) # Remove data pieces with no value 
    gooddata2 = np.logical_and(np.logical_and(wave1<25, wave2<25), np.logical_and(wave3<25, wave4<25))
    greatdata = np.logical_and(gooddata1, gooddata2)
    
    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]
    
    x = data[:,0][greatdata]
    y = data[:,1][greatdata]
    r = data[:,2][greatdata]
    d = data[:,3][greatdata]

    # cut out some additional strange outliers
    # data = data[~((data['alphFe'] > 0.4) & (data['FeH'] > -0.3))]

    X = np.vstack([colour1, colour2]).T

    #----------------------------------------------------------------------
    # Compute clustering with MeanShift
    #
    # We'll work with the scaled data, because MeanShift finds circular clusters

    X_scaled = preprocessing.scale(X)

    # The following bandwidth can be automatically detected using
    # the routine estimate_bandwidth().  Because bandwidth estimation
    # is very expensive in memory and computation, we'll skip it here.

    bandwidth = estimate_bandwidth(X)
    #bandwidth = 0.7

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False)
    ms.fit(X_scaled)

    labels_unique = np.unique(ms.labels_)
    n_clusters = len(labels_unique[labels_unique >= 0])
    
    if make_plots: 
        make_ms_plots(colour1, colour2, n_clusters, X, ms, band1, band2, band3, band4)
    
    return(n_clusters)
    #print labels_unique
    #print bandwidth
    #print "number of estimated clusters : %d" % n_clusters

    #------------------------------------------------------------
    
    
def make_ms_plots(colour1, colour2, n_clusters, X, ms, band1, band2, band3, band4):    
    # Plot the results
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    # plot density
    H, C1_bins, C2_bins = np.histogram2d(colour1, colour2, 51)

    ax.imshow(H.T, origin='lower', interpolation='nearest', aspect='auto',
              extent=[C1_bins[0], C1_bins[-1],
                  C2_bins[0], C2_bins[-1]],
              cmap=plt.cm.binary)

    # plot clusters
    #colors = ['b', 'g', 'r', 'k']

    for i in range(n_clusters):
        Xi = X[ms.labels_ == i]
        H, b1, b2 = np.histogram2d(Xi[:, 0], Xi[:, 1], (C1_bins, C2_bins))

        bins = [0.1]

        ax.contour(0.5 * (C1_bins[1:] + C1_bins[:-1]),
               0.5 * (C2_bins[1:] + C2_bins[:-1]),
               H.T, bins, colors='r')

    ax.xaxis.set_major_locator(plt.MultipleLocator(0.3))
    #ax.set_xlim(-1.101, 0.101)
    #ax.set_ylim(C2_bins[0], C2_bins[-1])
    ax.set_xlabel(band1+' - '+band2)
    ax.set_ylabel(band3+' - '+band4)
    
    plt.show()
    
    return ()
    
def do_kmeans(band1, band2, band3, band4, number_clusters, make_plots, output_cluster_id):

    '''do K-means clustering on colours constructed from HST photometry band1, band 2,
    band3, band4 are keys from band_names --- ie, names of HST  filters'''
    
    #do some input checking
    
    if band1 == band2 or band3 == band4: 
        print "Not a good idea to use the same band in one colour, try again"
        return
    for band in [band1, band2, band3, band4]:
        if band not in band_names.keys():
            print "Can't find %s in band_name list" %band
            return
    
    #Import Data from WFC# ERS M83 Data Products: .txt file is catalog of all sources
    #Assign data sets for clustering
    
    data = np.loadtxt('hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1.txt')
    
    #Import 4 different wavelengths
    #Colour 1: 05_mag
    wave1 = data[:,band_names[band1]]
    wave2 = data[:,band_names[band2]]
    
    #Colour 2: 05_mag
    wave3 = data[:,band_names[band3]]
    wave4 = data[:,band_names[band4]]
    
    gooddata1 = np.logical_and(np.logical_and(wave1!=-99, wave2!=-99), np.logical_and(wave3!=-99, wave4!=-99)) # Remove data pieces with no value 
    gooddata2 = np.logical_and(np.logical_and(wave1<25, wave2<25), np.logical_and(wave3<25, wave4<25))
    greatdata = np.logical_and(gooddata1, gooddata2)
    
    colour1 = wave1[greatdata] - wave2[greatdata]
    colour2 = wave3[greatdata] - wave4[greatdata]
    
    x = data[:,0][greatdata]
    y = data[:,1][greatdata]
    r = data[:,2][greatdata]
    d = data[:,3][greatdata]

    id = data[:,4][greatdata].astype(np.int32)

    #Put data in the right format for clustering
    
    clusterdata = np.vstack([colour1, colour2]).T
    
    #Truncate data for speed
    #clusterdata = clusterdata[::5]
    #--------------------------------------------------------------------------
    
    #Compute KMeans Clustering
    
    scaler = preprocessing.StandardScaler() 
    clf = KMeans(number_clusters, random_state = 10)
    clf.fit(scaler.fit_transform(clusterdata))
    
    cluster_number = clf.predict(scaler.fit_transform(clusterdata))

#    print cluster_number

# output object and cluster IDs to a file
    if output_cluster_id:
        file_name = 'ID_'+str(number_clusters)+'cl_'+band1+'-'+band2+'vs'+band3+'-'+band4+'.txt'
        tmptab = Table([id,cluster_number])
        tmptab.write(file_name, format='ascii.no_header')
    
    #Compute the score
    # kmeans_model = KMeans(n_clusters = 3, random_state = 1).fit(clusterdata)
        
    labels = clf.labels_
    score = metrics.silhouette_score(scaler.fit_transform(clusterdata), labels, metric = 'euclidean')
#    print >> results, 'Silhouette score for bands %s %s %s %s: %f' % (band1, band2, band3, band4, score)
    
    #Print number of objects per cluster and send to txt file 'results'
#    print >> results, 'Cluster# #objects'
        
    objects_per_cluster = np.zeros(max_num_clusters,dtype=np.int16)
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        objects_per_cluster[i] = len(x_cluster)
#       print >> results, i, len(x_cluster)

    objects_per_cluster.sort() # sort from smallest to largest 
    objects_per_cluster = objects_per_cluster[::-1] # reverse sort
    
    if make_plots: 
        colour_kmeans_plot (band1, band2, band3, band4, clf, scaler, colour1, colour2, number_clusters)
        xy_plot (x, y, number_clusters, cluster_number, band1, band2, band3, band4)
            
    return(score, objects_per_cluster)
    # Do some stats 
    # (TBD) -- results of this should perhaps be included in return()

    #--------------------------------------------------------------------------
    
    #Visualize results -- should probably be a separate function
def colour_kmeans_plot(band1, band2, band3, band4, clf, scaler, colour1, colour2, number_clusters):
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
        
    # Compute 2D histogram of the input 
        
    H, C1_bins, C2_bins = np.histogram2d(colour1, colour2, 50)
        
    #Plot Density 
        
    ax = plt.axes()
    ax.imshow(H.T, origin = 'lower', interpolation ='nearest', aspect='auto', 
                  extent = [C1_bins[0], C1_bins[-1],
                            C2_bins[0], C2_bins[-1]],
                  cmap=plt.cm.binary)
                  
    #Plot Cluster centers
    cluster_centers = scaler.inverse_transform(clf.cluster_centers_)
        
    for i in range(0, number_clusters):
        ax.scatter(cluster_centers[i, 0], cluster_centers[i, 1],
                       s=40, c=cluster_colours[i], edgecolors='k')
        
        #Plot cluster boundaries 
        
    C1_centers = 0.5 * (C1_bins[1:] + C1_bins[:-1])
    C2_centers = 0.5 * (C2_bins[1:] + C2_bins[:-1])
        
    clusterdatagrid = np.meshgrid(C1_centers, C2_centers)
    clusterdatagrid = np.array(clusterdatagrid).reshape((2, 50 * 50)).T
        
    H = clf.predict(scaler.transform(clusterdatagrid)).reshape((50,50))
        

    for i in range(number_clusters):
        Hcp = H.copy()
        flag = (Hcp == i)
        Hcp[flag] = 1 
        Hcp[~flag] = 0
        
        ax.contour(C1_centers, C2_centers, Hcp, [-0.5, 0.5],

                       linewidths=1, c='k')
                       
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.3))
        
    ax.set_xlabel(band1+' - '+band2)
    ax.set_ylabel(band3+' - '+band4)
        
    file_name = 'Cluster_'+str(number_clusters)+'cl_'+band1+'-'+band2+'vs'+band3+'-'+band4+'.png'
    pylab.savefig(file_name)
        
    return ()

def xy_plot (x, y, number_clusters, cluster_number, band1, band2, band3, band4):
        # plot xy positions of objects in different clusters
   

    fig2 = plt.figure(figsize=(5,5))
    ax2 = fig2.add_subplot(111)
    
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        y_cluster = y[cluster_number == i]
        ax2.scatter(x_cluster, y_cluster, label = i, c=cluster_colours[i])

    #ax2.title('Clustering in colours '+band1+' - '+band2+' vs '+band3+' - '+band4)
    ax2.set_xlabel('X [pixels]')
    ax2.set_ylabel('Y [pixels]')
    ax2.legend()
    #plt.show()
    

    filename = 'XY_'+str(number_clusters)+'cl_'+band1+'-'+band2+'vs'+band3+'-'+band4+'.png'
    pylab.savefig(filename)
    
    return ()

    
def results_summary(input_file = 'results.txt'):
    '''compute and plot summaries for clustering analysis'''

    # read in the data -- this is not an ideal way to do it since it requires knowledge of file structure
    results_table = Table.read(input_file, format='ascii.no_header')
    num_clust = results_table['col5']
    score = results_table['col6']
    total_obj = results_table['col7'].astype('float')

    # add a column with the size of the smallest cluster
    # have to do some tricky stuff since column corresponding to smallest cluster varies dep on number of clusters
    last_clust_col = 7+ num_clust
    results_table.add_column(Column(name='size_smallest', data=np.zeros(len(results_table)),dtype=np.int16))
    for i in range(0,len(results_table)):
        lastcol = 'col{}'.format(last_clust_col[i])
        results_table['size_smallest'][i] = results_table[lastcol][i]

    biggest_clust_fract = results_table['col8']/total_obj # compute fraction of objects in largest cluster
    smallest_clust_fract = results_table['size_smallest']/total_obj # compute fraction of objects in smallest cluster

    fig, ax = plt.subplots()
    ax.scatter(num_clust, score)
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Score')
    fig.show()
    
    filename = 'n_clusters_vs_Score.png'
    #pylab.savefig(filename)
    
    #N = results_table['col5']
    
    
   
    results_3 = results_table[results_table['col5']==3]
    number_3trials = len(results_3)
    X = np.arange(0,number_3trials)
    A =  results_3['col8']/results_3['col7'].astype('float')
    B = results_3['col9']/results_3['col7'].astype('float')
    C = results_3['col10']/results_3['col7'].astype('float')
    
    fig, ax = plt.subplots()
    ax.bar (X, A, color = 'b')
    ax.bar (X, B, color = 'r', bottom = A)
    ax.bar (X, C, color = 'g', bottom = B+A)
    ax.set_xlabel('Trial Number')
    ax.set_ylabel('Fractional Size per Cluster')
    fig.show()
    
    results_4 = results_table[results_table['col5']==4]
    number_4trials = len(results_4)
    X = np.arange(0,number_4trials)
    A =  results_4['col8']/results_4['col7'].astype('float')
    B = results_4['col9']/results_4['col7'].astype('float')
    C = results_4['col10']/results_4['col7'].astype('float')
    D = results_4['col11']/results_4['col7'].astype('float')
    
    fig, ax = plt.subplots()
    ax.bar (X, A, color = 'b')
    ax.bar (X, B, color = 'r', bottom = A)
    ax.bar (X, C, color = 'g', bottom = B+A)
    ax.bar (X, D, color = 'k', bottom = A+B+C)
    ax.set_xlabel('Trial Number')
    ax.set_ylabel('Fractional Size per Cluster')
    fig.show()
    
    results_5 = results_table[results_table['col5']==5]
    number_5trials = len(results_5)
    X = np.arange(0,number_5trials)
    A =  results_5['col8']/results_5['col7'].astype('float')
    B = results_5['col9']/results_5['col7'].astype('float')
    C = results_5['col10']/results_5['col7'].astype('float')
    D = results_5['col11']/results_5['col7'].astype('float')
    E = results_5['col12']/results_5['col7'].astype('float')    
    
    fig, ax = plt.subplots()
    ax.bar (X, A, color = 'b')
    ax.bar (X, B, color = 'r', bottom = A)
    ax.bar (X, C, color = 'g', bottom = B+A)
    ax.bar (X, D, color = 'k', bottom = A+B+C)
    ax.bar (X, E, color = 'm', bottom = A+B+C+D)
    ax.set_xlabel('Trial Number')
    ax.set_ylabel('Fractional Size per Cluster')
    fig.show()
    
    results_6 = results_table[results_table['col5']==6]
    number_6trials = len(results_6)
    X = np.arange(0,number_6trials)
    A =  results_6['col8']/results_6['col7'].astype('float')
    B = results_6['col9']/results_6['col7'].astype('float')
    C = results_6['col10']/results_6['col7'].astype('float')
    D = results_6['col11']/results_6['col7'].astype('float')
    E = results_6['col12']/results_6['col7'].astype('float')    
    F = results_6['col13']/results_6['col7'].astype('float')
    
    fig, ax = plt.subplots()
    ax.bar (X, A, color = 'b')
    ax.bar (X, B, color = 'r', bottom = A)
    ax.bar (X, C, color = 'g', bottom = B+A)
    ax.bar (X, D, color = 'k', bottom = A+B+C)
    ax.bar (X, E, color = 'm', bottom = A+B+C+D)
    ax.bar (X, F, color = 'y', bottom = A+B+C+D+E)
    ax.set_xlabel('Trial Number')
    ax.set_ylabel('Fractional Size per Cluster')
    fig.show()
    
    results_7 = results_table[results_table['col5']==7]
    
    number_7trials = len(results_7)
    X = np.arange(0,number_7trials)
    A =  results_7['col8']/results_7['col7'].astype('float')
    B = results_7['col9']/results_7['col7'].astype('float')
    C = results_7['col10']/results_7['col7'].astype('float')
    D = results_7['col11']/results_7['col7'].astype('float')
    E = results_7['col12']/results_7['col7'].astype('float')    
    F = results_7['col13']/results_7['col7'].astype('float')
    G = results_7['col14']/results_7['col7'].astype('float')
    
    fig, ax = plt.subplots()
    ax.bar (X, A, color = 'b')
    ax.bar (X, B, color = 'r', bottom = A)
    ax.bar (X, C, color = 'g', bottom = B+A)
    ax.bar (X, D, color = 'k', bottom = A+B+C)
    ax.bar (X, E, color = 'm', bottom = A+B+C+D)
    ax.bar (X, F, color = 'y', bottom = A+B+C+D+E)
    ax.bar (X, G, color = 'c', bottom = A+B+C+D+E+F)
    ax.set_xlabel('Trial Number')
    ax.set_ylabel('Fractional Size per Cluster')
    fig.show()
   
    results_8 = results_table[results_table['col5']==8]
    number_8trials = len(results_8)
    
    X = np.arange(0,number_8trials)
    A =  results_8['col8']/results_8['col7'].astype('float')
    B = results_8['col9']/results_8['col7'].astype('float')
    C = results_8['col10']/results_8['col7'].astype('float')
    D = results_8['col11']/results_8['col7'].astype('float')
    E = results_8['col12']/results_8['col7'].astype('float')    
    F = results_8['col13']/results_8['col7'].astype('float')
    G = results_8['col14']/results_8['col7'].astype('float')
    H = results_8['col15']/results_8['col7'].astype('float')
    
    fig, ax = plt.subplots()
    ax.bar (X, A, color = 'b')
    ax.bar (X, B, color = 'r', bottom = A)
    ax.bar (X, C, color = 'g', bottom = B+A)
    ax.bar (X, D, color = 'k', bottom = A+B+C)
    ax.bar (X, E, color = 'm', bottom = A+B+C+D)
    ax.bar (X, F, color = 'y', bottom = A+B+C+D+E)
    ax.bar (X, G, color = 'c', bottom = A+B+C+D+E+F)
    ax.bar (X, H, color = 'w', bottom = A+B+C+D+E+F+G)
    ax.set_xlabel('Trial Number')
    ax.set_ylabel('Fractional Size per Cluster')
    fig.show()
    
    
    
        
    #fig, ax = plt.subplots()
    #ax.scatter(num_clust, biggest_clust_fract)
    #ax.scatter(num_clust, smallest_clust_fract)
    #ax.set_xlabel('Number of clusters')
    #ax.set_ylabel('Fractional size')
    #fig.show()

    filename = 'n_clusters_vs_FractionalSize.png'
    #pylab.savefig(filename)

    fig, ax = plt.subplots()
    ax.scatter(score, biggest_clust_fract)
    ax.scatter(score, smallest_clust_fract)
    ax.set_xlabel('Score')
    ax.set_ylabel('Fractional size')
    fig.show()
    
    filename = 'Score_vs_FractionalSize.png'
    #pylab.savefig(filename)

    return()

