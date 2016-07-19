'''Results Summary '''
''' Create various plots and statistics from the results of
clustering trials'''
import numpy as np
import matplotlib.pyplot as plt
import pylab
import os, os.path
from astropy.table import Table, Column

'''-------------------------------------------------------------------------'''
# Used for plots
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',
          'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',
          'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*', '*',
           '.', '.', '.', '.', '.', '.', '.', '>', '>', '>', '>', '>', '>', '>',
           '+', '+', '+', '+', '+', '+', '+', '<', '<', '<', '<', '<', '<', '<',]
max_num_clusters = 40
'''-------------------------------------------------------------------------'''


def results(file_name, general_path, save_path, plots):
    '''Main function of results.py'''
    # general_results_file: used for bar_graphs and paramaters vs. silhouette score
    general_results_file = file_name
    gen_path = make_path(general_path)
    plot_path = make_directory(save_path)

    general_results_data = load_data(gen_path, general_results_file)

    if general_results_file == '05aperture_results_2d.txt':
        clean_gen_res_data, n_clust = organize_data(general_results_data)
    else:
        clean_gen_res_data = general_results_data

    if "gen_score" in plots:
        silhouette_vs_nclust(clean_gen_res_data, plot_path)
    if "ms_score" in plots:
        bandwidth_vs_score(clean_gen_res_data, plot_path)
    if "af_score" in plots:
        affinity_par_vs_score(clean_gen_res_data, plot_path)
    if "bar" in plots:
        bar_graph(clean_gen_res_data, n_clust, plot_path)
    if "clust_cent" in plots:
        cluster_centers_plot(clean_gen_res_data, plot_path)
    if "inertia" in plots:
        inertia_fluctuations(clean_gen_res_data, plot_path)
    if "percentile" in plots:
        object_cluster_dist(clean_gen_res_data, n_clust, plot_path)

    return


def load_data(sum_path, file_):
    '''Load all results files from m83_clustering\results'''
    summary_results = Table.read(os.path.join(sum_path, file_),
                                 format='ascii.commented_header', guess=False)
    return(summary_results)


def organize_data(general_data):
    '''Select the type of clustering etc.'''

    num_clust = general_data['n_clust']

    # add a column with the size of the smallest cluster
    # have to do some tricky stuff since column corresponding to smallest
    # cluster varies dep on number of clusters
    general_data.add_column(Column(name='size_smallest',
                                   data=np.zeros(len(general_data)),
                                   dtype=np.int16))
    for i in range(0, len(general_data)):
        lastcol = 'c_{}'.format(num_clust[i])
        general_data['size_smallest'][i] = general_data[lastcol][i]

    return(general_data, num_clust)


def make_directory(p_path):
    '''Save each result
        - set path of where you would like results saved'''
    pl_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(
              p_path)
    if not os.path.exists(pl_path):
        os.makedirs(pl_path)

    return(pl_path)


def make_path(sum_path):
    '''Set up file paths'''
    g_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(sum_path)
    return(g_path)


def silhouette_vs_nclust(results_table, path):
    '''Create plots of features vs. the silhouette_score'''
    # Remove no scores
    remove = []
    for i in range(0, len(results_table)):
        if results_table['score'][i] == -99.0:
            remove.append(i)
    results_table.remove_rows(remove)

    num_clust = results_table['n_clust']
    s_score = results_table['score']
    total_obj = results_table['total_objects'].astype('float')
    # compute fraction of objects in largest cluster
    biggest_clust_fract = results_table['c_1']/total_obj
    # compute fraction of objects in smallest cluster
    smallest_clust_fract = results_table['size_smallest']/total_obj

    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(121)
    ax.scatter(num_clust[results_table['clustering'] == 'kmeans'],
               results_table['score'][results_table['clustering']=='kmeans'],
               c='r', label='kmeans')
    ax.scatter(num_clust[results_table['clustering'] == 'meanshift'],
               results_table['score'][results_table['clustering']=='meanshift'],
               c='b', label='meanshift')
    ax.scatter(num_clust[results_table['clustering'] == 'affinity'],
               results_table['score'][results_table['clustering']=='affinity'],
               c='y', label='affinity')
    ax.legend(loc='best', fontsize=11)
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Score')
    ax.set_title('Number of Clusters vs Silhouette Score', fontsize=11)

    ax = fig.add_subplot(122)
    ax.scatter(biggest_clust_fract, s_score, c='r', marker='o',
               label='Largest Cluster')
    ax.scatter(smallest_clust_fract, s_score, c='b', marker='o',
               label='Smallest Cluster')
    ax.legend(loc='best', fontsize=11)
    ax.set_xlabel('Fractional size')
    ax.set_ylabel('Score')
    ax.set_title('Silhouette Score vs. Fractional Size of Cluster',
                 fontsize=11)

    filename = 'silhouette_score_plots.png'
    pylab.savefig(os.path.join(path, filename))
    return


def bandwidth_vs_score(results_table, path):
    '''Create a plot of the bandwidth from mean-shift vs. the
    silhouette_score'''
    # Remove now bandwidth or score
    remove = []
    for i in range(0, len(results_table)):
        if results_table['b_width'][i] == 'N/A':
            remove.append(i)
        if results_table['score'][i] == 'N/A':
            remove.append(i)
    results_table.remove_rows(remove)

    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(121)
    ax.scatter(results_table['b_width'], results_table['score'], marker='o',
               s=6)
    ax.set_xlabel('Bandwidth', fontsize=11)
    ax.set_ylabel('Score', fontsize=11)
    ax.set_title('Bandwidth vs. Score', fontsize=12)

    ax = fig.add_subplot(122)
    ax.scatter(results_table['b_width'], results_table['n_clust'], marker='o',
               s=6)
    ax.set_xlabel('Bandwidth', fontsize=11)
    ax.set_ylabel('N_clusters', fontsize=11)
    ax.set_title('Bandwidth vs. N_clusters', fontsize=12)

    filename = 'meanshift_parameters.png'
    pylab.savefig(os.path.join(path, filename))
    return


def affinity_par_vs_score(results_table, path):
    '''Plot affinity parameters vs. the silhouette score'''
    remove = []
    for i in range(0, len(results_table)):
        if results_table['clustering'][i] != 'affinity':
            remove.append(i)
        if results_table['score'][i] == 'N/A':
            remove.append(i)
    results_table.remove_rows(remove)

    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(121)
    ax.scatter(results_table['damp'], results_table['score'], marker='o')
    ax.set_xlabel('Damping Factor', fontsize=11)
    ax.set_ylabel('Score', fontsize=11)
    ax.set_title('Damping vs. Score', fontsize=12)

    ax = fig.add_subplot(122)
    ax.scatter(results_table['pref'], results_table['score'], marker='o')
    ax.set_xlabel('Preferences', fontsize=11)
    ax.set_ylabel('Score', fontsize=11)
    ax.set_title('Preferences vs. Score', fontsize=12)

    filename = 'afpar_vs_score.png'
    pylab.savefig(os.path.join(path, filename))

    fig2 = plt.figure(figsize=(12, 5))
    ax2 = fig2.add_subplot(121)
    ax2.scatter(results_table['damp'], results_table['n_clust'], marker='o',
                s=6)
    ax2.set_xlabel('Damping Factor', fontsize=11)
    ax2.set_ylabel('N_clusters', fontsize=11)
    ax2.set_title('Damping vs. N_clusters', fontsize=12)

    ax2 = fig2.add_subplot(122)
    ax2.scatter(results_table['pref'], results_table['n_clust'], marker='o',
                s=6)
    ax2.set_xlabel('Preferences', fontsize=11)
    ax2.set_ylabel('N_clusters', fontsize=11)
    ax2.set_title('Preferences vs. N_clusters', fontsize=12)

    filename = 'afpar_vs_nclust.png'
    pylab.savefig(os.path.join(path, filename))

    return


def rms_fluctuations(cluster_data):
    '''Make a plot of the root mean square value for a given cluster'''
    return


def inertia_fluctuations(results_table, path):
    '''Make a plot of the inertia value vs n_clust for kmeans clustering'''
    num_clust = results_table['n_clust'][results_table['clustering'] == 'kmeans']
    inertia = results_table['inertia'][results_table['clustering'] == 'kmeans']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.plot(num_clust, inertia, 'k', lw=0.5, zorder=1)
    ax.scatter(num_clust, inertia, marker='o', s=20, c='k', zorder=2)
    ax.set_xlabel('N_clusters', fontsize=11)
    ax.set_ylabel('Inertia (Sum of Squares)', fontsize=11)
    ax.set_title('N_clusters vs. Inertia')

    file_name = 'inertia_plot.png'
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    return


def cluster_centers_plot(cluster_data, path):
    '''Plot the centers of each cluster after a center_test'''
    # Find the number of trials
    n_trials = np.arange(0, len(cluster_data[cluster_data['clust_num'] == 1]))
    # Make table of cluster centers and cluster number
    # cen_1: colour 1 center coordinate  
    cluster_centers = np.array(cluster_data['clust_num', 'cen_1', 'cen_2'])
    # Find the number of clusters imposed
    n_clusters = max(cluster_data['clust_num'])
    fig = plt.figure(figsize=(12, 8))
    for c in range(1, 3):  # Loop over the number of colours used
        ax = fig.add_subplot(1, 2, c)  # Create subplot for each colour
        for i in range(0, len(n_trials)):  # Loop over each trial in each colour
            for k in range(0, len(cluster_centers)):  # Loop over every center point
                if cluster_centers['clust_num'][k] == 1:  # Find the first cluster of a clustering
                    for n in range(k, k+n_clusters):  # Loop over each cluster
                        ax.scatter(i, cluster_centers['cen_' + str(c)][n],
                                   s=10, c=colors[n-k])  # Plot the center of each cluster
        ax.set_ylabel('Colour ' + str(c) + ' Cluster Centers', fontsize=10)
        ax.set_xlabel('Trial', fontsize=10)
    plt.suptitle(str(len(n_trials)) + ' Trials vs. Cluster Centers in each Colour')
    filename = '{}_cl_cluster_centers.png'.format(str(n_clusters))
    pylab.savefig(os.path.join(path, filename))
    plt.show()

    return


def bar_graph(results_table, num_clust, path):
    '''Create bar graph for each number of clusters vs. fractional size of
    each cluster'''
    # Bar graph
    for i in range(1, max(num_clust)+1):
        numberoftrials = len(num_clust[num_clust == i])  # Find  trials
        if numberoftrials > 0:
            fig, ax = plt.subplots()
            X = np.arange(0, numberoftrials)
            for n in range(0, i):  # Create stacked bar graph
                if n == 0:
                    yprev = (0*X).astype('float')
                else:
                    yprev = Y + yprev
                colname = 'c_{}'.format(n+1)
                # Compute percentage of objects in each cluster and graph
                Y = results_table[colname][num_clust==i]/results_table['total_objects'][num_clust==i].astype('float')
                ax.bar(X, Y, color=colors[n], bottom=yprev)
                ax.set_title('Proportional Size - ' + str(i) + ' Clusters')
                ax.set_xlabel('Number of Trials')
                ax.set_ylabel('Fractional Cluster Size')
                ax.xaxis.set_major_locator(plt.MultipleLocator(1))
            filename = '{}_clusters_vs_FractionalSize.png'.format(i)
            pylab.savefig(os.path.join(path, filename))
    return


def object_cluster_dist(results_table, num_clust, path):
    '''Determine which clusters hold ~99% of the objects and save percentage
    in cluster_percentage.txt'''
    filename = 'cluster_percentage.txt'
    header = '# clustering total_obj c_1 c_2 c_3 c_4 c_5 c_6 '\
             'c_7 c_8 c_9 c_10 c_11 c_12 c_13 c_14 c_15 c_16 c_17 c_18 c_19 '\
             'c_20 c_21 c_22 c_23 c_24 c_25 c_26 c_27 c_28 c_29 c_30 c_31 c_32'\
             ' c_33 c_34 c_35 c_36 c_37 c_38 c_39 c_40'

    test_path = '{}\\{}'.format(path, filename)
    if not os.path.exists(test_path):
        create_path = os.path.join(path, filename)
        cluster_percentage = open(create_path, "w")
        cluster_percentage.write(header + '\n')
        cluster_percentage.close()
    cluster_percentage_path = os.path.join(path, filename)
    

    for i in range(0, len(num_clust)):
        cluster_percentage = open(cluster_percentage_path, "a")
        cumulative_percentage_cluster = np.zeros(max_num_clusters,
                                                 dtype=float)
        total_obj = float(results_table['total_objects'][i])
        n_clust = results_table['n_clust'][i]
        for n in range(1, n_clust+1):
            if n == 1:
                P = float(results_table['c_'+str(n)][i])/float(total_obj)
            else:
                percentage = float(results_table['c_'+str(n)][i])/total_obj
                P = P + percentage
            cumulative_percentage_cluster[n-1] = float(P)
            inputs = ('{} {} {}').format(results_table['clustering'][i],  
                                         total_obj,
                                         np.array_str(np.around(cumulative_percentage_cluster,
                                                                decimals=4),
                                                      max_line_width=500)[1:-1])
        cluster_percentage.write(inputs + '\n')
        cluster_percentage.close()

    return
