print(__doc__)
import matplotlib.pyplot as plt
from itertools import cycle
import numpy as np
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs
from sklearn import metrics
from sklearn import preprocessing
from astropy.table import Table
from scipy.spatial import distance

data = Table.read('data_v3.txt', format='ascii.commented_header',
                  guess=False)
m83_data = data#[:]
aperture = '05'
band1 = 'mag{}_555'.format(aperture)
band2 = 'mag{}_814'.format(aperture)

band3 = 'mag{}_373'.format(aperture)
band4 = 'mag{}_555'.format(aperture)

band5 = 'mag{}_487'.format(aperture)
band6 = 'mag{}_555'.format(aperture)

band7 = 'mag{}_502'.format(aperture)
band8 = 'mag{}_555'.format(aperture)

band9 = 'mag{}_555'.format(aperture)
band10 = 'mag{}_657'.format(aperture)

band11 = 'mag{}_555'.format(aperture)
band12 = 'mag{}_673'.format(aperture)

ratio = 0.1
n_clust = 14
incriment = 3
#colours = ['b', 'y', 'r', 'g', 'm', 'c', 'k', 'b', 'y', 'r', 'g', 'm', 'c', 'k']
#markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*', '*', '^', '>', '<', '+']

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
# Colour 3
wave5 = m83_data[band5]
wave5_unc = m83_data[band5+'_unc']
wave6 = m83_data[band6]
wave6_unc = m83_data[band6+'_unc']
# Colour 4
wave7 = m83_data[band7]
wave7_unc = m83_data[band7+'_unc']
wave8 = m83_data[band8]
wave8_unc = m83_data[band8+'_unc']
# Colour 5
wave9 = m83_data[band9]
wave9_unc = m83_data[band9+'_unc']
wave10 = m83_data[band10]
wave10_unc = m83_data[band10+'_unc']
# Colour 6
wave11 = m83_data[band11]
wave11_unc = m83_data[band11+'_unc']
wave12 = m83_data[band12]
wave12_unc = m83_data[band12+'_unc']

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

colour1_trim = np.logical_and(wave1_trim, wave2_trim)
colour2_trim = np.logical_and(wave3_trim, wave4_trim)
colour3_trim = np.logical_and(wave5_trim, wave6_trim)
colour4_trim = np.logical_and(wave7_trim, wave8_trim)
colour5_trim = np.logical_and(wave9_trim, wave10_trim)
colour6_trim = np.logical_and(wave11_trim, wave12_trim)

gooddata1 = np.logical_and(colour1_trim, colour2_trim)
gooddata2 = np.logical_and(colour3_trim, colour4_trim)
gooddata3 = np.logical_and(colour5_trim, colour6_trim)

# Only data that match criteria for both colours
final_data = np.logical_and(np.logical_and(gooddata1, gooddata2), gooddata3)

colour1 = wave1[final_data] - wave2[final_data]
colour2 = wave3[final_data] - wave4[final_data]
colour3 = wave5[final_data] - wave6[final_data]
colour4 = wave7[final_data] - wave8[final_data]
colour5 = wave9[final_data] - wave10[final_data]
colour6 = wave11[final_data] - wave12[final_data]

x = data['x'][final_data]
y = data['y'][final_data]

X = np.vstack([colour1, colour2, colour3, colour4, colour5, colour6]).T
print "objects: {}".format(len(X))
#X_scaled = preprocessing.scale(X)
# The following bandwidth can be automatically detected using
#bandwidth = estimate_bandwidth(X)
cluster_centers = []
b_width = min(distance.pdist(X))

while(len(cluster_centers) != 1):
    ms = MeanShift(bandwidth=b_width, bin_seeding=False, cluster_all=True)
    ms.fit(X)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_
    #print("N_clust_centers: {}").format(len(cluster_centers))
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    #print("Total Objects: {}").format(len(colour1))
    print("bandwidth: {:.4f}").format(b_width)
    print("number of estimated clusters : %d" % n_clusters_)
    if 50 > n_clusters_ > 1: 
        print("Silhouette Coefficient: %0.3f"
            % metrics.silhouette_score(X, labels))

    #print cluster_centers
    #sample_score = metrics.silhouette_samples(X, labels)
    ###############################################################################
    # Plot result

    #plt.clf()
    markers = cycle('ooooooo.......*******+++++++<<<<<<>>>>>>>^^^^^^^')
    colors = cycle('byrkmgcbyrkmgcbyrkmgcbyrkmgcbyrkmgc')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for k, col, mark in zip(range(n_clusters_), colors, markers):
        #cluster_score = sample_score[labels == k]
        my_members = labels == k
        cluster_center = cluster_centers[k]
            
        x_clust = x[labels == k]
        y_clust = y[labels ==k]
        ax.scatter(X[my_members, 0], X[my_members, 1], color=col, marker=mark,
                   s=4)
        #print ("objects in cluster {}: {} {}").format(k+1, len(X[my_members]), col)
        #print("    Cluster {} Average Score: {:.4f}").format(k+1,
        #     np.average(cluster_score))
        ax.plot(cluster_center[0], cluster_center[1], marker=mark, color=col,
                markersize=10)
    ax.set_xlabel("Colour 1")
    ax.set_ylabel("Colour " + str(i+1))
    ax.set_title('Estimated number of clusters: %d' % n_clusters_)
    plt.show()
    if b_width < 0.25:
        b_width = b_width*incriment
    else: 
        b_width = b_width*1.1
    # X = np.vstack(cluster_centers)
