# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import norm

from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import preprocessing



#----------------------------------------------------------------------

#------------------------------------------------------------
# Get the data
data = np.loadtxt('hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1.txt')
    
    #Import 4 different wavelengths
    #Colour 1: 05_mag
wave1 = data[:,11]
wave2 = data[:,15]
    
    #Colour 2: 05_mag
wave3 = data[:,19]
wave4 = data[:,23]
    
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

#bandwidth = estimate_bandwidth(X)
bandwidth = 0.4

ms = MeanShift(bandwidth=bandwidth, bin_seeding=True, cluster_all=False)
ms.fit(X_scaled)

labels_unique = np.unique(ms.labels_)
n_clusters = len(labels_unique[labels_unique >= 0])
print labels_unique
print bandwidth
print "number of estimated clusters : %d" % n_clusters

#------------------------------------------------------------
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
colors = ['b', 'g', 'r', 'k']

for i in range(n_clusters):
    Xi = X[ms.labels_ == i]
    H, b1, b2 = np.histogram2d(Xi[:, 0], Xi[:, 1], (C1_bins, C2_bins))

    bins = [0.1]

    ax.contour(0.5 * (C1_bins[1:] + C1_bins[:-1]),
               0.5 * (C2_bins[1:] + C2_bins[:-1]),
               H.T, bins, colors='w')

ax.xaxis.set_major_locator(plt.MultipleLocator(0.3))
ax.set_xlim(-1.101, 0.101)
ax.set_ylim(C2_bins[0], C2_bins[-1])
ax.set_xlabel(r'$\rm [Fe/H]$')
ax.set_ylabel(r'$\rm [\alpha/Fe]$')

plt.show()