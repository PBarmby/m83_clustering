'''Final paper clustering figures''' 
import numpy as np
import os
import pylab as pylab
from astropy.table import Table, join
from matplotlib import pyplot as plt

# Dictionary for plot formatting
band_names = {'mag05_225':'F225W', 'mag05_336':'F336W', 'mag05_373':'F373N',
              'mag05_438':'F438W', 'mag05_487':'F487N', 'mag05_502':'F502N',
              'mag05_555':'F555W', 'mag05_657':'F657N', 'mag05_673':'F673N',
              'mag05_814':'F814W'}
save_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\Paper\\figs\\final_colour_figures\\'
survey_data = Table.read('data_v3.txt', format='ascii.commented_header',
                         guess=False)

# SET VALUES
band_1 = 'mag05_555'
band_2 = 'mag05_814'
band_3 = 'mag05_336'
band_4 = 'mag05_438'
n_clust = 5
alg = 'kmeans'   # Change center colours to normal for meanshift, otherwise '0.65' for GS
col = 'k'
symbol = ['s', 'o', '^', '*', 'o', '>', '<']
plot_colours = ['0.5', 'k', '0.85', '0.25', '0.65', '0.125', '0.85']
#plot_colours =['m', 'y', 'b', 'g', 'm', 'c', 'k']
file_name = 'figure_6_GS.png'
specific_path = 'broad_band\\U-B_V-I\\clustering'
colours = ('{}-{}_{}-{}').format(band_1, band_2, band_3, band_4)

# MODEL COLOURS
model_file = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\model_colours\\mist_ssp_feh+0.5.txt'
model = Table.read(model_file, format='ascii.commented_header', guess=False)
mod_filt = np.vstack([model[band_1], model[band_2],
                      model[band_3], model[band_4]])

#LOAD ID FILE AND CLUSTER CENTERS
id_file_name = 'id_{}_{}cl_{}-{}vs{}-{}.txt'.format(alg, n_clust,
                                                    band_1, band_2,
                                                    band_3, band_4)
id_file_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}\\{}\\{}\\{}'.format(specific_path, colours, 'id_', id_file_name)
id_file = Table.read(id_file_path, format='ascii.commented_header',
                     guess=False)
cen_file_name = 'cluster_statistics.txt'
cen_file_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}\\{}\\{}'.format(specific_path, colours, cen_file_name)
centers = Table.read(cen_file_path, format='ascii.commented_header',
                     guess=False)

# GET CLUSTERING DATA FROM ID FILES AND SURVEY
cluster_data_2d = join(id_file, survey_data, 'object_id')

#CREATE COLOURS
wave1 = cluster_data_2d[band_1]
wave2 = cluster_data_2d[band_2]
wave3 = cluster_data_2d[band_3]
wave4 = cluster_data_2d[band_4]
colour1 = wave1 - wave2
colour2 = wave3 - wave4

# FIND CLUSTER CENTERS
clustering_data = centers[centers['clustering'] == alg]
num_clust_data = clustering_data[clustering_data['total_clust'] == n_clust]

# Colour-Colour plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

for k in range(0, n_clust):
    clust_col = plot_colours[k]
    my_members = cluster_data_2d['cluster_number'] == k

    cluster_center_1 = num_clust_data['cen_1'][k]
    cluster_center_2 = num_clust_data['cen_2'][k]
    ax.scatter(colour1[my_members], colour2[my_members], color=clust_col,   #Change
               marker=symbol[k], label=k+1, s=2, zorder=1)
    if symbol[k]=='*':
        ax.scatter(cluster_center_1, cluster_center_2, marker=symbol[k],    #Change
                   color='0.65', edgecolor=col, s=200, zorder=2)
    else:
        ax.scatter(cluster_center_1, cluster_center_2, marker=symbol[k],    #Change
                   color='0.65', edgecolor=col, s=100, zorder=2)

# Plot model colours
ax.plot(mod_filt[0] - mod_filt[1], mod_filt[2] - mod_filt[3], color=col)    #Change
                                                     # Model Colour
ax.scatter(mod_filt[0][0] - mod_filt[1][0], mod_filt[2][0] - mod_filt[3][0],    #Change
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


# Format plot
ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
ax.set_xlabel(band_names[band_1]+ ' - ' + band_names[band_2],   # Change
              fontweight='bold', fontsize=14)
ax.set_ylabel(band_names[band_3]+' - ' + band_names[band_4],    # Change
              fontweight='bold', fontsize=14)

legend = ax.legend(loc='lower right', fontsize=11, scatterpoints=1)
for s in range (0, n_clust):
    legend.legendHandles[s]._sizes = [30]

pylab.savefig(os.path.join(save_path, file_name))
plt.show()
plt.close()



                                          