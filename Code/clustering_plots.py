'''Make clustering plots from id.txt files generated by clustering_(x)d.py
1. Set bands to create in plots.txt or plots_3d.txt
2. Set base1, base2, b_wave to reflect the "base colour" if using plots_3d.txt
3. Run from Spyder terminal
    Argument 1 (dimensions): Specify '2d' or '3d'
    Argument 2 (model_name): Select colour model to overlay
        - Change header in any other file than 'mist_ssp_feh+0.5.txt' so that
            filters match header from data_v3.txt
    Argument 3 (f_path): path to id_ file from specific clustering
    Argument 4 (s_path): path to save new figure'''

import numpy as np
import os
import pylab as pylab
from astropy.table import Table, join
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d as p3

'''-------------------------------------------------------------------------'''
# Set base colour
base1 = 3
base2 = 5
b_wave = 5

# Dictionary for plot formatting
band_names = {'mag05_225':'f225w', 'mag05_336':'f336w', 'mag05_373':'f373n',
              'mag05_438':'f438w', 'mag05_487':'f487n', 'mag05_502':'f502n',
              'mag05_555':'f555w', 'mag05_657':'f657n', 'mag05_673':'f673n',
              'mag05_814':'f814w'}
'''-------------------------------------------------------------------------'''

def plotting(dimensions, model_name, f_path, s_path):
    model_file = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\model_colours\\{}'.format(model_name)
    survey_data = Table.read('data_v3.txt', format='ascii.commented_header',
                             guess=False)
    save_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\{}\\'.format(s_path)
    model = Table.read(model_file, format='ascii.commented_header', guess=False)
    if '2d' in dimensions:  # Check if 2D or 3D 
        plots = Table.read('plots.txt', format='ascii.commented_header',
                           guess=False)  # Load bands to plot

        for i in range(0, len(plots)):  # Loop over all combos
            waves = plots[i]
            colours = ('{}-{}_{}-{}').format(waves[0], waves[1], waves[2],
                                             waves[3])
            mod_filt_2d = np.vstack([model[waves[0]], model[waves[1]],
                                     model[waves[2]], model[waves[3]]])
            n_clusters = plots['n_clust'][i]  # Find the number of clusters
            algorithm = plots['clustering'][i]  # Which algorithm

            # Load id file generated from clustering_2d.py
            id_file_2d, centers_2d = load_files(n_clusters, algorithm,
                                                waves, dimensions, colours,
                                                f_path)
            
            cluster_data_2d = join(id_file_2d, survey_data, 'object_id')
            colour1_2d, colour2_2d, wave_2, wave_4 = load_2d_cluster_data(cluster_data_2d,
                                                          waves)
            make_2d_plots(colour1_2d, colour2_2d, waves, n_clusters,
                          algorithm, cluster_data_2d, centers_2d,
                          save_path, wave_2, wave_4, mod_filt_2d)

    if '3d' in dimensions:
        plots_3d = Table.read('plots_3d.txt', format='ascii.commented_header',
                              guess=False)
        for j in range(0, len(plots_3d)):
            waves = plots_3d[j]
            colours = ('{}-{}_{}-{}_{}-{}').format(waves[0], waves[1],
                                                   waves[2], waves[3],
                                                   waves[4], waves[5])
            mod_filt_3d = np.vstack([model[waves[0]], model[waves[1]],
                                     model[waves[2]], model[waves[3]],
                                     model[waves[4]], model[waves[5]]])
            n_clusters = plots_3d['n_clust'][j]
            algorithm = plots_3d['clustering'][j]
            id_file_3d, centers_3d = load_files(n_clusters, algorithm,
                                                waves, dimensions,
                                                colours, f_path)
            cluster_data_3d = join(id_file_3d, survey_data, 'object_id')
            colour1, colour2, colour3, base, wave_base1, wave_base2 = load_3d_cluster_data(cluster_data_3d,
                                                             waves)
            make_3d_plots(colour1, colour2, colour3, waves, n_clusters,
                          algorithm, cluster_data_3d, centers_3d,
                          save_path, base, wave_base1, wave_base2, mod_filt_3d)
    return()


def load_files(n_clust, clustering, bands, dim, cols, cen_path):
    if '2d' in dim:
        id_file_name = 'id_{}_{}cl_{}-{}vs{}-{}.txt'.format(clustering,
                                                            n_clust, bands[0],
                                                            bands[1], bands[2],
                                                            bands[3])
        cen_file_name = 'cluster_statistics.txt'
    if '3d' in dim:
        id_file_name = 'id_{}_{}cl_{}-{}vs{}-{}.txt'.format(clustering,
                                                            n_clust, bands[0],
                                                            bands[1],
                                                            bands[base1],
                                                            bands[base2])
        cen_file_name = '3d_cluster_statistics.txt'
    id_file_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}\\{}\\{}\\{}'.format(cen_path, cols, 'id_', id_file_name)
    id_file = Table.read(id_file_path, format='ascii.commented_header',
                         guess=False)
    cen_file_path = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\results\\{}\\{}\\{}'.format(cen_path, cols, cen_file_name)
    cluster_center_data = Table.read(cen_file_path,
                                     format='ascii.commented_header',
                                     guess=False)
    return(id_file, cluster_center_data)


def load_2d_cluster_data(filter_data, bands):
    wave1 = filter_data[bands[0]]
    wave2 = filter_data[bands[1]]
    wave3 = filter_data[bands[2]]
    wave4 = filter_data[bands[3]]
    colour1 = wave1 - wave2
    colour2 = wave3 - wave4
    
    return(colour1, colour2, wave2, wave4)


def load_3d_cluster_data(filter_data, bands):
    wave1 = filter_data[bands[0]]
    wave2 = filter_data[bands[1]]
    wave3 = filter_data[bands[2]]
    wave4 = filter_data[bands[3]]
    wave5 = filter_data[bands[4]]
    wave6 = filter_data[bands[5]]
    colour1 = wave1 - wave2
    colour2 = wave3 - wave4
    colour3 = wave5 - wave6
    base_colour = filter_data[bands[base1]] - filter_data[bands[base2]]
    base_wave = filter_data[bands[b_wave]]
    return(colour1, colour2, colour3, base_colour, base_wave, wave2)


def make_2d_plots(c1, c2, bands, n_clust, alg, c_data, centers, path,
                  base_wave1, base_wave2, model_data):
    colours = ('{}-{}_{}-{}').format(bands[0], bands[1], bands[2],
                                           bands[3])
    # Find cluster centers
    clustering_data = centers[centers['clustering'] == alg]
    num_clust_data = clustering_data[clustering_data['total_clust'] == n_clust]

    # Colour-Colour plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for k in range(0, n_clust):
        clust_col = plt.cm.jet(float(k) / np.max(c_data['cluster_number'] + 1))
        my_members = c_data['cluster_number'] == k
        cluster_center_1 = num_clust_data['cen_1'][k]
        cluster_center_2 = num_clust_data['cen_2'][k]
        ax.scatter(c1[my_members], c2[my_members], color=clust_col,
                   marker='o', label=k+1, s=2, zorder=1)
        ax.scatter(cluster_center_1, cluster_center_2, marker='o',
                   color=clust_col, edgecolor='k', s=100, zorder=2)
    # Plot model colours
    ax.plot(model_data[0] - model_data[1], model_data[2] - model_data[3],
            color='r')
    # Format plot
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(band_names[bands[0]] + ' - ' + band_names[bands[1]])
    ax.set_ylabel(band_names[bands[2]]+' - '+band_names[bands[3]])
    ax.set_title(alg + ' ' + str(n_clust) + ' : '+band_names[bands[0]]+'-'+band_names[bands[1]]+' vs. '+band_names[bands[2]]+'-'+band_names[bands[3]],
                 fontsize=14)
    ax.legend(loc='lower right', fontsize=8)
    plt.show()
    file_name = '{}_color_{}cl_{}-{}vs{}-{}.png'.format(alg, str(n_clust),
                                                        band_names[bands[0]],
                                                        band_names[bands[1]],
                                                        band_names[bands[2]],
                                                        band_names[bands[3]])
    path_ = ('{}{}{}').format(path, '\\', colours)
    if not os.path.exists(path_):
        print "New_path"
        os.makedirs(path_)
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    # CMD1
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for a in range(0, n_clust):
        my_members = c_data['cluster_number'] == a
        clust_col = plt.cm.jet(float(a) / np.max(c_data['cluster_number'] + 1))
        ax1.scatter(c1[my_members], base_wave1[my_members],
                    color=clust_col, marker='.', s=4, label=a+1)
    ax1.legend(loc='lower right', fontsize=8)
    ax1.set_xlabel(band_names[bands[0]]+' - '+band_names[bands[1]])
    ax1.set_ylabel(band_names[bands[1]])
    plt.gca().invert_yaxis()
    plt.show()
    '''Display interactive figure if # removed, if not, figures saved'''
    file_name = '{}_CMD_{}cl_{}-{}vs{}.png'.format(alg, str(n_clust),
                                                   band_names[bands[0]],
                                                   band_names[bands[1]],
                                                   band_names[bands[1]])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    # CMD2
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    for c in range(0, n_clust):
        my_members = c_data['cluster_number'] == c
        clust_col = plt.cm.jet(float(c) / np.max(c_data['cluster_number'] + 1))
        ax2.scatter(c2[my_members], base_wave2[my_members],
                    color=clust_col, marker='.', s=4, label=c+1)
    ax2.legend(loc='lower right', fontsize=8)
    ax2.set_xlabel(band_names[bands[2]]+' - '+band_names[bands[3]])
    ax2.set_ylabel(band_names[bands[3]])
    plt.gca().invert_yaxis()
    plt.show()
    '''Display interactive figure if # removed, if not, figures saved'''
    file_name = '{}_CMD_{}cl_{}-{}vs{}.png'.format(alg, str(n_clust),
                                                   band_names[bands[2]],
                                                   band_names[bands[3]],
                                                   band_names[bands[3]])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    return


def make_3d_plots(c1, c2, c3, bands, n_clust, alg, c_data, centers,
                  path, base_colour, base_wave, wave2, model_data):
    colours = ('{}-{}_{}-{}_{}-{}').format(bands[0], bands[1], bands[2],
                                           bands[3], bands[4], bands[5])
    # Find cluster centers
    clustering_data = centers[centers['clustering'] == alg]
    num_clust_data = clustering_data[clustering_data['total_clust'] == n_clust]

    # Plot each colour against the original narrow-broad
    for i in range(1, 3):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for k in range(0, n_clust):
            clust_col = plt.cm.jet(float(k) / np.max(c_data['cluster_number'] + 1))
            my_members = c_data['cluster_number'] == k
            cluster_center_1 = num_clust_data['cen_1'][k]
            cluster_center_2 = num_clust_data['cen_2'][k]
            cluster_center_3 = num_clust_data['cen_3'][k]
            if i == 1:
                ax.scatter(c1[my_members], c2[my_members], color=clust_col,
                           marker='o', label=k+1, s=2, zorder=1)
                ax.scatter(cluster_center_1, cluster_center_2, marker='o',
                           color=clust_col, edgecolor='k', s=100, zorder=2)
            else:
                ax.scatter(c1[my_members], c3[my_members], color=clust_col,
                            marker='o', label=k+1, s=2, zorder=1)
                ax.scatter(cluster_center_1, cluster_center_3, marker='o',
                            color=clust_col, edgecolor='k', s=100, zorder=2)
        # Plot model colours
        ax.plot(model_data[0] - model_data[1],
                model_data[i*2] - model_data[i*2+1], color='r')
        # Format plot
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
        ax.set_xlabel(band_names[bands[0]] + ' - ' + band_names[bands[1]])
        ax.set_ylabel(band_names[bands[i*2]]+' - '+band_names[bands[i*2+1]])
        ax.set_title(alg + ' ' + str(n_clust) + ': '+band_names[bands[0]]+'-'+band_names[bands[1]]+' vs. '+band_names[bands[i*2]]+'-'+band_names[bands[i*2+1]],
                     fontsize=14)
        ax.legend(loc='lower left', fontsize=8)
        plt.show
        file_name = '{}_3d_projection_{}cl_{}-{}vs{}-{}.png'.format(alg, str(n_clust),
                                                               band_names[bands[0]],
                                                               band_names[bands[1]],
                                                               band_names[bands[i*2]],
                                                               band_names[bands[i*2+1]])
        path_ = ('{}{}{}').format(path, '\\', colours)
        if not os.path.exists(path_):
            print "made path"
            os.makedirs(path_)
        pylab.savefig(os.path.join(path_, file_name))
        plt.close()

    # 3d plot
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    for x in range(0, n_clust):
        clust_col = plt.cm.jet(float(x) / np.max(c_data['cluster_number'] + 1))
        cluster_center_1 = num_clust_data['cen_1'][x]
        cluster_center_2 = num_clust_data['cen_2'][x]
        cluster_center_3 = num_clust_data['cen_3'][x]
        members = c_data['cluster_number'] == x
        ax1.scatter(c1[members], c2[members], c3[members],
                    marker='o', color=clust_col, s=2, label=x+1)
        ax1.scatter(cluster_center_1, cluster_center_2, cluster_center_3,
                    c=clust_col, edgecolor='k', marker='o', s=100)
    # Plot model colours
    ax1.plot(model_data[0] - model_data[1], model_data[2] - model_data[3],
             model_data[4] - model_data[5], color='r')
    # Format plot
    ax1.set_xlabel(band_names[bands[0]] + ' - ' + band_names[bands[1]])
    ax1.set_ylabel(band_names[bands[2]] + ' - ' + band_names[bands[3]])
    ax1.set_zlabel(band_names[bands[4]] + ' - ' + band_names[bands[5]])
    ax1.legend(loc='upper left', fontsize=8)
    file_name = '{}_3d_{}cl_{}-{}vs{}-{}vs{}-{}.png'.format(alg, str(n_clust),
                                                            band_names[bands[0]],
                                                            band_names[bands[1]],
                                                            band_names[bands[2]],
                                                            band_names[bands[3]],
                                                            band_names[bands[4]],
                                                            band_names[bands[5]])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    # Base colour plot
    fig2 = plt.figure(figsize=(12,12))
    ax2 = fig2.add_subplot(111)
    for b in range(0, n_clust):
        clust_col = plt.cm.jet(float(b) / np.max(c_data['cluster_number'] + 1))
        objects = c_data['cluster_number'] == b
        cluster_center_1 = num_clust_data['cen_1'][b]
        cluster_center_2 = num_clust_data['cen_2'][b]
        cluster_center_3 = num_clust_data['cen_3'][b]
        base_cen = cluster_center_3 - cluster_center_2  # CHANGE
        ax2.scatter(base_colour[objects], c1[objects], marker='o',
                    color=clust_col, s=2, label=b+1, zorder=1)
        ax2.scatter(base_cen, cluster_center_1, marker='o',
                    color=clust_col, s=100, edgecolor='k', zorder=2)
    # Plot model colours
    ax2.plot(model_data[base1] - model_data[base2],
             model_data[0] - model_data[1], color='r')
    # Format Plot
    ax2.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax2.set_ylabel(band_names[bands[0]] + ' - ' + band_names[bands[1]])
    ax2.set_xlabel(band_names[bands[base1]]+' - ' + band_names[bands[base2]])
    plt.gca().invert_yaxis()
    ax2.set_title(alg + ' ' + str(n_clust) + ' : '+band_names[bands[0]]+'-'+band_names[bands[1]]+' vs. '+band_names[bands[base1]]+'-'+band_names[bands[base2]],
                  fontsize=14)
    ax2.legend(loc='lower right', fontsize=8)
    file_name = '{}_base_color_{}cl_{}-{}vs{}-{}.png'.format(alg, str(n_clust),
                                                             band_names[bands[0]],
                                                             band_names[bands[1]],
                                                             band_names[bands[base1]],
                                                             band_names[bands[base2]])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    # CMD1
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)
    for a in range(0, n_clust):
        my_members = c_data['cluster_number'] == a
        clust_col = plt.cm.jet(float(a) / np.max(c_data['cluster_number'] + 1))
        ax4.scatter(c1[my_members], wave2[my_members],
                    color=clust_col, marker='.', s=4, label=a+1)
    ax4.legend(loc='lower right', fontsize=8)
    ax4.set_xlabel(band_names[bands[0]]+' - '+band_names[bands[1]])
    ax4.set_ylabel(band_names[bands[1]])
    plt.gca().invert_yaxis()

    '''Display interactive figure if # removed, if not, figures saved'''
    file_name = '{}_CMD_{}cl_{}-{}vs{}.png'.format(alg, str(n_clust),
                                                          band_names[bands[0]],
                                                          band_names[bands[1]],
                                                          band_names[bands[1]])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()

    # CMD2
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    for c in range(0, n_clust):
        my_members = c_data['cluster_number'] == c
        clust_col = plt.cm.jet(float(c) / np.max(c_data['cluster_number'] + 1))
        ax3.scatter(base_colour[my_members], base_wave[my_members],
                    color=clust_col, marker='.', s=4, label=c+1)
    ax3.legend(loc='lower right', fontsize=8)
    ax3.set_xlabel(band_names[bands[base1]]+' - '+band_names[bands[base2]])
    ax3.set_ylabel(band_names[bands[base2]])
    plt.gca().invert_yaxis()

    '''Display interactive figure if # removed, if not, figures saved'''
    file_name = '{}_CMD_{}cl_{}-{}vs{}.png'.format(alg, str(n_clust),
                                                          band_names[bands[base1]],
                                                          band_names[bands[base2]],
                                                          band_names[bands[b_wave]])
    pylab.savefig(os.path.join(path_, file_name))
    plt.close()
    return
