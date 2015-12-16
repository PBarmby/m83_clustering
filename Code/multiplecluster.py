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
# Can be changed based on file headers 
band_names = {'05_225': 11, '3_225': 13, '05_336': 15,  '3_336': 17, '05_373': 19, '3_373': 21, '05_438':23,    '3_438':25,  '05_487':27, '3_487': 29, '05_502':31, '3_502':33, '05_555': 35, '3_555': 37, '05_657': 39 ,'3_657':41 ,'05_673':43, '3_673':45 , '05_814':47 , '3_814':49 }

# used for plots
cluster_colours = ['y','g','b','r','c','m','k','m','w','y','g','b','r','c','m','k','m','w','y','g','b','r','c','m','k','m','w'] 

#Input file, columns correspond to different wavelengths 

inputdata  = 'hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1.txt'
data = Table.read(inputdata, format = 'ascii.commented_header', guess = False)
   
# need this so that output files always have the same number of columns
max_num_clusters = 20

# Choose analysis and output

def userinput(): 
    
    print "**Please enter a space between each input**" 
    analysis = raw_input("What analysis would you like to perform? (meanshift, kmeans, mst): ")
    
    if analysis in ['kmeans']: 
        KMEANS_number_cluster = raw_input("Would you like to use the number of clusters estimated by MEANSHIFT as input to KMEANS? (Yes/No): ")
    
    make_plots = raw_input("What plots would you like(cluster, color, xy, mst): ")
    classification_ID = raw_input("Would you like the objects to be catalogued (Yes/No): ")
    make_results_summary = raw_input("Would you like a results summary (Yes/No): ")

    #Save the inputs 
    
    inputs = [analysis, KMEANS_number_cluster, make_plots, classification_ID, make_results_summary]

    #Check Inputs

    for i in range (0,3):
    
        while inputs[i] not in ['meanshift', 'kmeans', 'mst', 'cluster', 'color', 'xy', 'mst', 'Yes', 'yes', 'No', 'no']: 
       
            print "Invalid input, please try again."
            print "**Please enter a space between each input**"
            analysis = raw_input("What analysis would you like to perform? (meanshift, kmeans, mst): ")
            make_plots = raw_input("What plots would you like(cluster, color, xy, mst): ")
            classification_ID = raw_input("Would you like the objects to be catalogued (Yes/No): ")
            make_results_summary = raw_input("Would you like a results summary (Yes/No): ")
    
    confirm = raw_input("Start clustering now? (Yes/Change my inputs): ")
    
    while confirm not in ['Yes', 'yes']:
    
        print "**Please enter a space between each input**" 
        analysis = raw_input("What analysis would you like to perform? (meanshift, kmeans, mst): ")
        make_plots = raw_input("What plots would you like(cluster, color, xy, tree): ")
        classification_ID = raw_input("Would you like the objects to be catalogued (Yes/No): ")
        make_results_summary = raw_input("Would you like a results summary (Yes/No): ")
        
        confirm = raw_input("Start clustering now? (Yes/Change my inputs): ")
    
    #do_everything(inputs)   # Pass user inputs to clustering functions 

    return

def do_everything(input_file = 'experiments.txt', output_file = 'results.txt', mp_MS=True, mp_KM_COLOUR=False, mp_KM_XY=True, oci=False, rs = False):
    '''Automate clustering process
       input: input_file:  a 4-column text file with 1 line per clustering run
                           each line lists the 4 filters to be used to construct colours
              mp: make output plots
              oci: output cluster IDs for each object
              rs: make graphs to summarize results
       output: output_file: a text file listing input+results from each clustering run'''
    
    #Run experiments listed in experiments.txt file 
    #Perform analysis based on user input
    
    run = np.genfromtxt(input_file, dtype='str')
       # perform_analysis = user_inputs  
    
    # check whether results file already exists; if not, open it and print a header line
    # if it does already exist, just open it
    
    results = open(output_file, 'a') 
    
    #Run clustering 
    
    for i in range(0, len(run)):
        
        '''Get Data'''
        #Colour 1 
        wave1 = data[run[i,0]]
        wave2 = data[run[i,1]]
    
        #Colour 2
        wave3 = data[run[i,2]]
        wave4 = data[run[i,3]]
    
        gooddata1 = np.logical_and(np.logical_and(wave1!=-99, wave2!=-99), np.logical_and(wave3!=-99, wave4!=-99)) # Remove data pieces with no value 
        gooddata2 = np.logical_and(np.logical_and(wave1<25, wave2<25), np.logical_and(wave3<25, wave4<25))  #Remove data above certain magnitude
        greatdata = np.logical_and(gooddata1, gooddata2)    #Only data that match criteria for both colours
    
        colour1 = wave1[greatdata] - wave2[greatdata]
        colour2 = wave3[greatdata] - wave4[greatdata]
    
        
    
    
    
    
        # MEANSHIFT to find the appropriate number of clusters
        numberofclusters = do_meanshift (run[i,0], run[i,1], run[i,2], run[i,3], colour1, colour2, mp_MS)
        
          
        
        
        
        
        
        
        input_str =  '{} {}'.format(np.array_str(run[i][:])[1:-1],numberofclusters) # list of input parameters: bands and num of clusters
        
               
        # run KMEANS clustering based on the number of clusters found using MEANSHIFT
        score, num_obj =  do_kmeans(run[i,0], run[i,1], run[i,2], run[i,3], colour1, colour2, greatdata, numberofclusters, mp_KM_COLOUR, mp_KM_XY, output_cluster_id=oci)
        total_obj = num_obj.sum()
        
        
        
        
        
        
        
        
        
        
        output_str = ' {:.4f} {:5d} {}'.format(score, total_obj, np.array_str(num_obj, max_line_width = 100)[1:-1])
        
        
        
        
        
        results.write(input_str + ' ' + output_str + '\n')
    
    results.close()
    
    
    
    
    
    
    
    
    if rs: 
        results_summary()
    
    return













def do_meanshift (band1, band2, band3, band4, colour1, colour2, make_plots):
    '''Does meanshift clustering to determine a number of clusters in the 
        data, which is passed to KMEANS function'''
        
    #Input Checking
    #if band1 == band2 or band3 == band4: 
        #print "Not a good idea to use the same band in one colour, try again"
        #return
    #for band in [band1, band2, band3, band4]:
        #if band not in band_names.keys():
            #print "Can't find %s in band_name list" %band
            #return
        
  
    #Truncate data
    X = np.vstack([colour1, colour2]).T

   
    # Compute clustering with MeanShift
   
    #Scale data because meanshift generates circular clusters 
    X_scaled = preprocessing.scale(X)

    # The following bandwidth can be automatically detected using
    # the routine estimate_bandwidth(). Bandwidth can also be set
    # as a value.

    bandwidth = estimate_bandwidth(X)
    
    # Meanshift clustering 
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False)
    ms.fit(X_scaled)

    labels_unique = np.unique(ms.labels_)
    n_clusters = len(labels_unique[labels_unique >= 0])
    
    #Make plot of clusters if needed
    if make_plots: 
        make_ms_plots(colour1, colour2, n_clusters, X, ms, band1, band2, band3, band4)
    
    return(n_clusters)
    

    
    
    

def do_kmeans(band1, band2, band3, band4, colour1, colour2, greatdata, number_clusters, make_colour, make_xy, output_cluster_id):

    '''do K-means clustering on colours constructed from HST photometry band1, band 2,
    band3, band4 are keys from band_names --- ie, names of HST  filters'''
   
    x = data['x'][greatdata]
    y = data['y'][greatdata]
    id = data['id'][greatdata].astype(np.int32)

    #Put data in the right format for clustering
    
    clusterdata = np.vstack([colour1, colour2]).T
    
     
    #Compute KMeans Clustering
    #Data pre-processing
    scaler = preprocessing.StandardScaler() 
    clf = KMeans(number_clusters, random_state = 10)
    
    clf.fit(scaler.fit_transform(clusterdata))
    
    cluster_number = clf.predict(scaler.fit_transform(clusterdata))

    # output object and cluster IDs to a file

    if output_cluster_id:
        file_name = 'ID_'+str(number_clusters)+'cl_'+band1+'-'+band2+'vs'+band3+'-'+band4+'.txt'
        tmptab = Table([id,cluster_number])
        tmptab.write(file_name, format='ascii.no_header')
    
    #Compute the silhouette score
    
        
    labels = clf.labels_
    score = metrics.silhouette_score(scaler.fit_transform(clusterdata), labels, metric = 'euclidean')

    #Identify which cluster each object belongs to 
    objects_per_cluster = np.zeros(max_num_clusters,dtype=np.int16)
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        objects_per_cluster[i] = len(x_cluster)

    objects_per_cluster.sort() # sort from smallest to largest 
    objects_per_cluster = objects_per_cluster[::-1] # reverse sort
    
    if make_colour: 
        colour_kmeans_plot (band1, band2, band3, band4, clf, scaler, colour1, colour2, number_clusters)
    
    if make_xy:        
        xy_plot (x, y, number_clusters, cluster_number, band1, band2, band3, band4)
            
    return(score, objects_per_cluster)
    
    
    
    
    
    
def make_ms_plots(colour1, colour2, n_clusters, X, ms, band1, band2, band3, band4):    
    ''' Plot the results of mean shift clustering if needed''' 
    
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    # plot density
    H, C1_bins, C2_bins = np.histogram2d(colour1, colour2, 51)

    ax.imshow(H.T, origin='lower', interpolation='nearest', aspect='auto',
              extent=[C1_bins[0], C1_bins[-1],
                  C2_bins[0], C2_bins[-1]],
              cmap=plt.cm.binary)

    # plot clusters
    

    for i in range(n_clusters):
        Xi = X[ms.labels_ == i]
        H, b1, b2 = np.histogram2d(Xi[:, 0], Xi[:, 1], (C1_bins, C2_bins))

        bins = [0.1]

        ax.contour(0.5 * (C1_bins[1:] + C1_bins[:-1]),
               0.5 * (C2_bins[1:] + C2_bins[:-1]),
               H.T, bins, colors='r')

    ax.xaxis.set_major_locator(plt.MultipleLocator(0.3))
   
    ax.set_xlabel(band1+' - '+band2)
    ax.set_ylabel(band3+' - '+band4)
    
    plt.show()
    
    return ()
    
    
    
    
    
    
    
    
    
    
    
    
def colour_kmeans_plot(band1, band2, band3, band4, clf, scaler, colour1, colour2, number_clusters):
    '''Plot cluster data for KMEANS clustering'''
    
    # Visualize results
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
        
    #Plot cluster centers
        
    C1_centers = 0.5 * (C1_bins[1:] + C1_bins[:-1])
    C2_centers = 0.5 * (C2_bins[1:] + C2_bins[:-1])
        
    clusterdatagrid = np.meshgrid(C1_centers, C2_centers)
    clusterdatagrid = np.array(clusterdatagrid).reshape((2, 50 * 50)).T
        
    H = clf.predict(scaler.transform(clusterdatagrid)).reshape((50,50))
        
    #Plot boundries 
    for i in range(number_clusters):
        Hcp = H.copy()
        flag = (Hcp == i)
        Hcp[flag] = 1 
        Hcp[~flag] = 0
        
        ax.contour(C1_centers, C2_centers, Hcp, [-0.5, 0.5],

                       linewidths=1, c='k')
                      
                      
    #Set axes for each plot
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.3))
        
    ax.set_xlabel(band1+' - '+band2)
    ax.set_ylabel(band3+' - '+band4)
        
    file_name = 'Cluster_'+str(number_clusters)+'cl_'+band1+'-'+band2+'vs'+band3+'-'+band4+'.png'
    pylab.savefig(file_name)
        
    return ()












def xy_plot (x, y, number_clusters, cluster_number, band1, band2, band3, band4):
    '''Plot xy positions of objects in different clusters'''
   
    fig2 = plt.figure(figsize=(5,5))
    ax2 = fig2.add_subplot(111)
    
    #Plot XY coordinates of each object and label with colour that identifies cluster number
    for i in range(0, number_clusters):
        x_cluster = x[cluster_number == i]
        y_cluster = y[cluster_number == i]
        ax2.scatter(x_cluster, y_cluster, label = i, c=cluster_colours[i])

    
    ax2.set_xlabel('X [pixels]')
    ax2.set_ylabel('Y [pixels]')
    ax2.legend()
    

    filename = 'XY_'+str(number_clusters)+'cl_'+band1+'-'+band2+'vs'+band3+'-'+band4+'.png'
    pylab.savefig(filename)
    
    return ()

    
    
    
    
    
    
    
    
    
    
def results_summary(input_file = 'results.txt'):
    '''compute and plot summaries for clustering analysis'''

    # read in the data -- this is not an ideal way to do it since it requires knowledge of file structure
    results_table = Table.read(input_file, format='ascii.commented_header', guess = False)
    num_clust = results_table['n_clusters']
    score = results_table['s_score']
    total_obj = results_table['total_objects'].astype('float')

    # add a column with the size of the smallest cluster
    # have to do some tricky stuff since column corresponding to smallest cluster varies dep on number of clusters
    last_clust_col = 7+ num_clust
    results_table.add_column(Column(name='size_smallest', data=np.zeros(len(results_table)),dtype=np.int16))
    for i in range(0,len(results_table)):
        lastcol = 'col{}'.format(last_clust_col[i])
        results_table['size_smallest'][i] = results_table[lastcol][i]

    biggest_clust_fract = results_table['c_1']/total_obj # compute fraction of objects in largest cluster
    smallest_clust_fract = results_table['size_smallest']/total_obj # compute fraction of objects in smallest cluster

    fig, ax = plt.subplots()
    ax.scatter(num_clust, score)
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Score')
    fig.show()
    
    filename = 'n_clusters_vs_Score.png'
    pylab.savefig(filename)
    
    #Bar graph
    
    for i in range (1,max(num_clust)+1):
        numberoftrials = len(num_clust[num_clust==i]) #Find number of trials in results.txt
         
        if numberoftrials > 0:
            fig, ax = plt.subplots()
                                  
            X = np.arange(0,numberoftrials)
            
            
            for n in range(0,i):   #Create stacked bar graph, add percentage of each cluster to the previous percentage
                if n == 0:
                    yprev = (0*X).astype('float')
                    
                else: 
                    yprev = Y+yprev
                    
                colname = 'col%d' % (n+8)
                Y = results_table[colname]/results_table['col7'].astype('float') #Compute percentage of objects in each cluster and graph
                ax.bar(X,Y,color = cluster_colours[n], bottom = yprev)  #Plot bar graph

    filename = 'n_clusters_vs_FractionalSize.png'
    pylab.savefig(filename)

    fig, ax = plt.subplots()
    ax.scatter(score, biggest_clust_fract)
    ax.scatter(score, smallest_clust_fract)
    ax.set_xlabel('Score')
    ax.set_ylabel('Fractional size')
    fig.show()
    
    filename = 'Score_vs_FractionalSize.png'
    pylab.savefig(filename)

    return()

