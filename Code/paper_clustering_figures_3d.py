'''Final paper clustering figures''' 
import numpy as np
import os
import pylab as pylab
from astropy.table import Table, join
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d as p3

# Dictionary for plot formatting
band_names = {'mag05_225':'F225W', 'mag05_336':'F336W', 'mag05_373':'F373N',
              'mag05_438':'F438W', 'mag05_487':'F487N', 'mag05_502':'F502N',
              'mag05_555':'F555W', 'mag05_657':'F657N', 'mag05_673':'F673N',
              'mag05_814':'F814W'}
save_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\Paper\\figs\\final_colour_figures\\'
survey_data = Table.read('data_v3.txt', format='ascii.commented_header',
                         guess=False)

# SET VALUES
band_1 = 'mag05_502'
band_2 = 'mag05_555'    # BASE FILTER
band_3 = 'mag05_336'
band_4 = 'mag05_438'
file_name = 'figure_10b_GS.png'
n_clust = 6
alg = 'kmeans'   # Change center colours to normal for meanshift, otherwise '0.65' for GS
col = 'k'
symbol = ['s', 'o', '^', '*', 'o', '>', '<']
plot_colours = ['0.4', '0.85', '0.5', '0.7', 'k', '0.175', '0.85']
#plot_colours =['m', 'y', 'b', 'g', 'k', 'c', 'k']
specific_path = 'broad_narrow\\OIII_V\\clustering'
colours = ('{}-{}_{}-{}_{}-{}').format(band_1, band_2, band_3, band_1,
                                       band_4, band_1)

# MODEL COLOURS
model_file = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\model_colours\\mist_ssp_feh+0.5.txt'
model = Table.read(model_file, format='ascii.commented_header', guess=False)
mod_filt = np.vstack([model[band_1], model[band_2], # 0 1
                      model[band_3], model[band_1], # 2 3
                      model[band_4], model[band_1]]) # 4 5

#LOAD ID FILE AND CLUSTER CENTERS
id_file_name = 'id_{}_{}cl_{}-{}vs{}-{}.txt'.format(alg, n_clust,
                                                    band_1, band_2,
                                                    band_3, band_4)
id_file_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}\\{}\\{}\\{}'.format(specific_path, colours, 'id_', id_file_name)
id_file = Table.read(id_file_path, format='ascii.commented_header',
                     guess=False)
cen_file_name = '3d_cluster_statistics.txt'
cen_file_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}\\{}\\{}'.format(specific_path, colours, cen_file_name)
centers = Table.read(cen_file_path, format='ascii.commented_header',
                     guess=False)

# GET CLUSTERING DATA FROM ID FILES AND SURVEY
cluster_data_3d = join(id_file, survey_data, 'object_id')

#CREATE COLOURS
wave1 = cluster_data_3d[band_1]
wave2 = cluster_data_3d[band_2]     #COMMON WAVE
wave3 = cluster_data_3d[band_3]
wave4 = cluster_data_3d[band_4]

colour1 = wave1 - wave2
mod_col_1 = mod_filt[0] - mod_filt[1]
colour2 = wave3 - wave1
mod_col_2 = mod_filt[2] - mod_filt[3]
colour3 = wave4 - wave1
mod_col_3 = mod_filt[4] - mod_filt[5]
base_colour = wave3 - wave4
base_mod_col = mod_filt[3] - mod_filt[5]    # CHANGE BASE MODEL COLOUR

# FIND CLUSTER CENTERS
clustering_data = centers[centers['clustering'] == alg]
num_clust_data = clustering_data[clustering_data['total_clust'] == n_clust]

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')    # CHANGE
'''
# Colour-Colour plot
for i in range(0, n_clust):
    clust_col = plot_colours[i]
    clust_cen_col = plot_colours[i]                     # Change'0.65' 
    my_members = cluster_data_3d['cluster_number'] == i
    cluster_center_1 = num_clust_data['cen_1'][i]
    cluster_center_2 = num_clust_data['cen_2'][i]
    cluster_center_3 = num_clust_data['cen_3'][i]
    base_cen = cluster_center_3 - cluster_center_2  # CHANGE
    ax.scatter(base_colour[my_members], colour1[my_members], color=clust_col,   #Change
               marker=symbol[i], label=i+1, s=2, zorder=1)
    if symbol[i]=='*':
        ax.scatter(base_cen, cluster_center_1, marker=symbol[i],    #Change
                   color=clust_cen_col, edgecolor=col, s=200, zorder=2)
    else:
        ax.scatter(base_cen, cluster_center_1, marker=symbol[i],    #Change
                   color=clust_cen_col, edgecolor=col, s=100, zorder=2)
# Plot model colours
ax.plot(base_mod_col, mod_col_1, color=col)    
ax.scatter(mod_filt[3][0] - mod_filt[5][0], mod_filt[0][0] - mod_filt[1][0],     #Change
           color=col, marker='o', s=70, zorder=2)

# Plot Chandar colours. Comment out if not UBVI combination
#Cluster Space
ax.plot([0.5, (2-0.33)/0.33], [0.5, 2], color=col)
ax.plot([0.5, -1], [0.5,-7], color=col)
ax.text(3, 0, "Cluster space")

# Star/Cluster space
ax.plot([-4, 0.24], [-0.8, -0.8], color=col)
ax.plot([0.24, -1], [-0.8, -7], color=col)
ax.text(-3, -2, "Star/Cluster space")

#Blue Star space
ax.plot([0.5, -4], [0.5, 0.5], color=col)
ax.text(-3, 0, "Blue Star space")

# Yellow star space
ax.text(-3, 1.25, "Yellow Star space")

# Plot limits - comment out if no Chandar colours used
plt.ylim(-3, 2)
plt.xlim(-4, 5)


ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.set_ylabel(band_names[band_1]+ ' - ' + band_names[band_2],   # Change
              fontweight='bold', fontsize=14)
ax.set_xlabel(band_names[band_3]+' - ' + band_names[band_4],    # Change
              fontweight='bold', fontsize=14)
'''


for x in range(0, n_clust):
    clust_col = plot_colours[x]
    cluster_center_1 = num_clust_data['cen_1'][x]
    cluster_center_2 = num_clust_data['cen_2'][x]
    cluster_center_3 = num_clust_data['cen_3'][x]
    members = cluster_data_3d['cluster_number'] == x
    ax.scatter(colour1[members], colour2[members], colour3[members],
               marker=symbol[x], color=clust_col, s=2, label=x+1)
    ax.scatter(cluster_center_1, cluster_center_2, cluster_center_3,
               color='0.65', #clust_col,             # CHANGE
               edgecolor='k', marker=symbol[x], s=100)
    # Plot model colours
ax.plot(mod_filt[0] - mod_filt[1], mod_filt[2] - mod_filt[3],
        mod_filt[4] - mod_filt[5], color='k')                          # Model colour
ax.scatter(mod_filt[0][0] - mod_filt[1][0],
           mod_filt[2][0] - mod_filt[3][0],
           mod_filt[4][0] - mod_filt[5][0], color='k',                 # Model colour
           marker='o', s=100, zorder=2)
# Format plot
ax.set_xlabel(band_names[band_1] + ' - ' + band_names[band_2],
              fontweight='bold', fontsize=14)
ax.set_ylabel(band_names[band_3] + ' - ' + band_names[band_1],
              fontweight='bold', fontsize=14)
ax.set_zlabel(band_names[band_4] + ' - ' + band_names[band_1],
              fontweight='bold', fontsize=14)

legend = ax.legend(loc='upper left', fontsize=11, scatterpoints=1)
for l in range (0, n_clust):
    legend.legendHandles[l]._sizes = [30]
pylab.savefig(os.path.join(save_path, file_name))
plt.show()
plt.close()


                                          