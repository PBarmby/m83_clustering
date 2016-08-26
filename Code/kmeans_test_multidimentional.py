'''Set bands and aperture
Performs K-Means clustering in 6 colour space and creates interactive plots
Choose colours to plot in plotting area by changing index of X (array of all colours)'''

import numpy as np 
from astropy.table import Table
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn import metrics, preprocessing
from itertools import cycle
import os, os.path
ratio = 0.2
n_clust = 7
colours = ['b', 'y', 'r', 'g', 'm', 'c', 'k', 'b', 'y', 'r', 'g', 'm', 'c', 'k', 'b', 'y', 'r', 'g', 'm', 'c', 'k', 'b', 'y', 'r', 'g', 'm', 'c', 'k']
markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*', '*', '^', '^','^','^','^','^','^','^','^','^','>', '<', '+']

m83_data = Table.read('data_v3.txt', format='ascii.commented_header',
                  guess=False)

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

colour1 = wave1[final_data] - wave2[final_data] #
colour2 = wave3[final_data] - wave4[final_data] #
colour3 = wave5[final_data] - wave6[final_data] #
colour4 = wave7[final_data] - wave8[final_data] #
colour5 = wave9[final_data] - wave10[final_data]
colour6 = wave11[final_data] - wave12[final_data]

x = m83_data['x'][final_data]
y = m83_data['y'][final_data]
id_ = m83_data['id_'][final_data]
X = np.vstack([colour1, colour2, colour3, colour4, colour5, colour6]).T #, colour4, colour5, colour6, colour7]).T

print ("objects: {}").format(len(colour1))            

def kmeans(cluster_data):
    # Data pre-processing
    scaler = preprocessing.StandardScaler()
    km = KMeans(n_clust, init='random')#random_state=10)
    km.fit(cluster_data) #change variable - data with labels
    #clf.fit(clusterdata)
    #cluster_number = km.predict(scaler.fit_transform(clusterdata))
    #cluster_number = clf.predict(clusterdata)
    centers = km.cluster_centers_
    # Compute the silhouette score
    labels = km.labels_
    score = metrics.silhouette_score(cluster_data, labels)
    sample_score = metrics.silhouette_samples(cluster_data, labels)
    sum_of_squares = km.inertia_
    return(labels, score, sample_score, centers, sum_of_squares)

labels_, score, cluster_score, center, sos = kmeans(X)
print ("Sum-of-Squares: {}").format(sos)
#for k in range(0, 4):
fig = plt.figure()
ax = fig.add_subplot(111)


'''Plotting'''
# Choose colours by changing X index
colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
markers = cycle('ooooooo*******')
for i, col, mark in zip(range(n_clust), colors, markers):
    class_members = labels_ == i
    cluster_center = center[i]
    ax.scatter(X[class_members, 0], X[class_members, 1], color=col,
             marker=mark, s=2)
    ax.plot(cluster_center[0], cluster_center[1], marker=mark,
            markerfacecolor=col, markeredgecolor='k', markersize=14)
plt.show()
print "Score: {}".format(score)
