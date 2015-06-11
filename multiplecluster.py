# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 13:24:41 2015

@author: Owner
"""
import numpy as np 
import pylab as pylab
from matplotlib import pyplot as plt 

from sklearn.cluster import KMeans
from sklearn import preprocessing 
from sklearn import metrics 
from sklearn.metrics import pairwise_distances


band_names = {'05_225': 11, '3_225': 13, '05_336': 15,  '3_336': 17, '05_373': 19, '3_373': 21, '05_438':23,    '3_438':25,  '05_487':27, '3_487': 29, '05_502':31, '3_502':33, '05_555': 35, '3_555': 37, '05_657': 39 ,'3_657':41 ,'05_673':43, '3_673':45 , '05_814':47 , '3_814':49 }

cluster_colours = ['y','g','b','r','c','m','k']

results = open('results.txt', 'a')

def do_everything():
    '''Automate clustering process'''
    
    run = np.genfromtxt('experiments.txt', dtype='str')
    
    for i in range(0, len(run)):
        
        title = run[i,0]+'-'+run[i,1]+'vs'+run[i,2]+'-'+run[i,3]
        print >> results, title
        
        print >> results, do_cluster(run[i,0], run[i,1], run[i,2], run[i,3], int(run[i,4]))

    return
  
    
def do_cluster(band1, band2, band3, band4, number_clusters):
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
    print cluster_number
    
    #results = open('results.txt', 'a')  #Send data to separate text file 'results'
    
    #Compute the score
    # kmeans_model = KMeans(n_clusters = 3, random_state = 1).fit(clusterdata)
    
        
    labels = clf.labels_
    score = metrics.silhouette_score(scaler.fit_transform(clusterdata), labels, metric = 'euclidean')
    print >> results, 'Silhouette score for bands %s %s %s %s: %f' % (band1, band2, band3, band4, score)
    
    #Print number of objects per cluster and send to txt file 'results'
    
    print >> results, 'Cluster# #objects'
        
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        print >> results, i, len(x_cluster)
    
    
    # Do some stats 
    
    
    
    
    
    
    
    #--------------------------------------------------------------------------
    
    #Visualize results 
    
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
    
    file_name = 'Cluster '+str(number_clusters)+'# '+band1+'-'+band2+'vs'+band3+'-'+band4+'.png'
    pylab.savefig(file_name)
    
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
    
    filename = 'XY '+str(number_clusters)+'# '+band1+'-'+band2+'vs'+band3+'-'+band4+'.png'
    pylab.savefig(filename)
    
    return   
    
    
