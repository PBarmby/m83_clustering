'''Final plots for paper''' 
import matplotlib.pyplot as plt
import pylab
import os, os.path
from astropy.table import Table, Column


path_373 = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\broad_narrow\\U_OII\\clustering\\mag05_336-mag05_373_mag05_373-mag05_438_mag05_373-mag05_814\\05aperture_results_3d.txt'
results_373 = Table.read(path_373, format='ascii.commented_header', guess=False)
results_225 = Table.read('C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\broad_narrow\\UVW_U\\clustering\\mag05_225-mag05_336_mag05_555-mag05_814\\05aperture_results_2d.txt', format='ascii.commented_header', guess=False)

# Remove no scores in 373
remove = []
for i in range(0, len(results_373)):
    if results_373['score'][i] == -99.0:
        remove.append(i)
results_373.remove_rows(remove)

# Remove no scores in 373
remove = []
for j in range(0, len(results_225)):
    if results_225['score'][j] == -99.0:
        remove.append(j)
results_225.remove_rows(remove)

num_clust_373 = results_373['n_clust']
s_score373 = results_373['score']
num_clust225 = results_225['n_clust']
s_score225 = results_225['score']

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
# plot 373 
#ax.scatter(num_clust_373 [results_373['clustering'] == 'kmeans'],
 #          results_373['score'][results_373['clustering']=='kmeans'],
  #         c='r', s=40, label='F373n_kmeans')
ax.scatter(num_clust_373 [results_373['clustering'] == 'kmeans'],
           results_373['score'][results_373['clustering']=='kmeans'], marker='o',
           c='r', s=40, label='F373N K-Means')
ax.scatter(num_clust_373 [results_373['clustering'] == 'meanshift'],
           results_373['score'][results_373['clustering']=='meanshift'],
           c='r', marker='s', s=40, label='F373n Mean Shift')
    
# plot 225
ax.scatter(num_clust225[results_225['clustering'] == 'kmeans'],
           results_225['score'][results_225['clustering']=='kmeans'],
           facecolors='none', edgecolors='b', s=40, label='F225W K-Means')
ax.scatter(num_clust225[results_225['clustering'] == 'meanshift'],
           results_225['score'][results_225['clustering']=='meanshift'],
           facecolors='none', edgecolors='b', marker='s', s=40, label='F225W Mean Shift')

ax.legend(loc='best', fontsize=11, scatterpoints=1)
ax.set_xlabel('Number of clusters', fontweight='bold', fontsize=14)
ax.set_ylabel('Score', fontweight='bold', fontsize=14)
# ax.set_title('Number of Clusters vs Silhouette Score', fontsize=11)

filename = 'score_vs_nclust.png'
pylab.savefig(os.path.join('C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\Paper\\figs\\final_figures', filename))

# Plot bandwidth vs. score and n_clust on same axis
remove = []
for i in range(0, len(results_373)):
    if results_373['b_width'][i] == 'N/A':
        remove.append(i)
results_373.remove_rows(remove)
    
fig1 = plt.figure(figsize=(8,8))
ax1 = fig1.add_subplot(211)

ax2 = plt.subplot(211)
plt.plot(results_373['b_width'], results_373['score'], marker='o',
            color='r')
plt.setp(ax1.get_xticklabels(), fontsize=6)
ax2.set_ylabel('Score', fontweight='bold', fontsize=14)
plt.ylim(0.45, 0.54)
ax2.yaxis.set_major_locator(plt.MultipleLocator(0.01))

# share x only
ax3 = plt.subplot(212, sharex=ax2)
plt.plot(results_373['b_width'], results_373['n_clust'], marker='o',
            color='k')
# make these tick labels invisible
plt.setp(ax2.get_xticklabels(), visible=False)
ax3.set_ylabel('Number of Clusters', fontweight='bold', fontsize=14)
ax3.set_xlabel('Bandwidth', fontweight='bold', fontsize=14)
plt.ylim(4, 11)

filename = 'bandwidth_plot.png'
pylab.savefig(os.path.join('C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\Paper\\figs\\final_figures', filename))
