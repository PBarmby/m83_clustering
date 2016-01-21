# -*- coding: utf-8 -*-
"""
Created on Mon Dec 07 20:59:40 2015

@author: Owner
"""
# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
import numpy as np
from matplotlib import pyplot as plt

from scipy import sparse
from sklearn.mixture import GMM

from astroML.clustering import HierarchicalClustering, get_graph_segments
from astropy.table import Table, Column


 

#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=False)

#------------------------------------------------------------
# get data

inputdata = 'hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1.txt'

data = Table.read(inputdata, format = 'ascii.commented_header', guess = False)

wave1 = data['mag05_225']
wave2 = data['mag05_814']
    
#Colour 2
wave3 = data['mag05_336']
wave4 = data['mag05_555']
    
gooddata1 = np.logical_and(np.logical_and(wave1!=-99, wave2!=-99), np.logical_and(wave3!=-99, wave4!=-99)) # Remove data pieces with no value 
gooddata2 = np.logical_and(np.logical_and(wave1<25, wave2<25), np.logical_and(wave3<25, wave4<25))  #Remove data above certain magnitude
greatdata = np.logical_and(gooddata1, gooddata2)    #Only data that match criteria for both colours

x = data['x'][greatdata].astype('float')
y = data['y'][greatdata].astype('float')

X = data['x','y']
print X[0]
print X[1]

xmin, xmax = (0, 5000)
ymin, ymax = (0, 5000)

#------------------------------------------------------------
# Compute the MST clustering model
n_neighbors = 10
edge_cutoff = 0.9
cluster_cutoff = 10
model = HierarchicalClustering(n_neighbors=n_neighbors, edge_cutoff=edge_cutoff, min_cluster_size=cluster_cutoff)
model.fit(X)


scale_ = np.percentile(model.full_tree_.data, 100 * edge_cutoff)
print " scale: %2g Mpc" % scale_ 

print " scale: %2g Mpc" % np.percentile(model.full_tree_.data, 100 * edge_cutoff)

n_components = model.n_components_
labels = model.labels_

#------------------------------------------------------------
# Get the x, y coordinates of the beginning and end of each line segment
T_x, T_y = get_graph_segments(model.X_train_, model.full_tree_)

T_trunc_x, T_trunc_y = get_graph_segments(model.X_train_, model.cluster_graph_)

#------------------------------------------------------------
# Fit a GMM to each individual cluster
Nx = 100
Ny = 250
Xgrid = np.vstack(map(np.ravel, np.meshgrid(np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny)))).T
density = np.zeros(Xgrid.shape[0])

for i in range(n_components):
    ind = (labels == i)
    Npts = ind.sum()
    Nclusters = min(12, Npts / 5)

    gmm = GMM(Nclusters).fit(X[ind])
    dens = np.exp(gmm.score(Xgrid))
    density += dens / dens.max()

density = density.reshape((Ny, Nx))

#----------------------------------------------------------------------
# Plot the results
fig = plt.figure(figsize=(5, 6))
fig.subplots_adjust(hspace=0, left=0.1, right=0.95, bottom=0.1, top=0.9)

ax = fig.add_subplot(311, aspect='equal')
ax.scatter(X[:, 1], X[:, 0], s=1, lw=0, c='k')
ax.set_xlim(ymin, ymax)
ax.set_ylim(xmin, xmax)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.set_ylabel('(Mpc)')

ax = fig.add_subplot(312, aspect='equal')
ax.plot(T_y, T_x, c='k', lw=0.5)
ax.set_xlim(ymin, ymax)
ax.set_ylim(xmin, xmax)
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.set_ylabel('(Mpc)')

ax = fig.add_subplot(313, aspect='equal')
ax.plot(T_trunc_y, T_trunc_x, c='k', lw=0.5)
ax.imshow(density.T, origin='lower', cmap=plt.cm.hot_r,
          extent=[ymin, ymax, xmin, xmax])

ax.set_xlim(ymin, ymax)
ax.set_ylim(xmin, xmax)
ax.set_xlabel('(Mpc)')
ax.set_ylabel('(Mpc)')

plt.show()