'''Results Summary '''
''' Create various plots and statistics from the results of
clustering trials'''
import numpy as np
import matplotlib.pyplot as plt
import pylab
import os, os.path
import argparse
from astropy.table import Table, Column

'''-------Enter the names of the results files you would like to use--------'''
# general_results_file: used for bar_graphs and paramaters vs. silhouette score
general_results_file = '05aperture_results.txt'

# clustering_results_file: used for clustering specific plots
#   - cluster center movement etc.
clustering_results_file = ' '

# Used for plots
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',
          'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k',
          'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k']
markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o', '*', '*', '*', '*', '*', '*', '*',
           '.', '.', '.', '.', '.', '.', '.', '>', '>', '>', '>', '>', '>', '>',
           '+', '+', '+', '+', '+', '+', '+', '<', '<', '<', '<', '<', '<', '<',]
'''-------------------------------------------------------------------------'''


def results(general_path, save_path, plots):
    '''Main function of results.py'''
    gen_path = make_path(general_path)
    plot_path = make_directory(save_path)

    general_results_data = load_data(gen_path)
    clean_gen_res_data, n_clust = organize_data(general_results_data)
    if "s_score" in plots:
        silhouette_plots(clean_gen_res_data, plot_path)
    if "bar" in plots:
        bar_graph(clean_gen_res_data, n_clust, plot_path)
    return


def load_data(sum_path):
    '''Load all results files from m83_clustering\results'''
    summary_results = Table.read(os.path.join(sum_path, general_results_file),
                                 format='ascii.commented_header', guess=False)
    return(summary_results)


def organize_data(general_data):
    '''Select the type of clustering etc.'''

    remove = []
    for i in range(0, len(general_data)):
        if general_data['score'][i] == 'N/A':
            remove.append(i)
    general_data.remove_rows(remove)

    num_clust = general_data['n_clust']
    # s_score = general_data['score']
    # total_obj = general_data['total_objects'].astype('float')

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


def silhouette_plots(results_table, path):
    '''Create plots of features vs. the silhouette_score'''

    num_clust = results_table['n_clust']
    s_score = results_table['score']
    total_obj = results_table['total_objects'].astype('float')
    # compute fraction of objects in largest cluster
    biggest_clust_fract = results_table['c_1']/total_obj
    # compute fraction of objects in smallest cluster
    smallest_clust_fract = results_table['size_smallest']/total_obj

    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax.scatter(num_clust, s_score)
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Score')
    ax.set_title('Number of Clusters vs Silhouette Score', fontsize=11)

    ax = fig.add_subplot(122)
    ax.scatter(biggest_clust_fract, s_score, c='r', marker='o',
               label='Largest Cluster')
    ax.scatter(smallest_clust_fract, s_score, c='b', marker='o',
               label='Smallest Cluster')
    ax.legend(loc='best')
    ax.set_xlabel('Fractional size')
    ax.set_ylabel('Score')
    ax.set_title('Silhouette Score vs. Fractional Size of Cluster',
                 fontsize=11)

    filename = 'silhouette_score_plots.png'
    pylab.savefig(os.path.join(path, filename))
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
                ax.set_title('Proportional Cluster Size')
                ax.set_xlabel('Number of Trials')
                ax.set_ylabel('Fractional Cluster Size')
                ax.xaxis.set_major_locator(plt.MultipleLocator(1))
            filename = '{}_clusters_vs_FractionalSize.png'.format(i)
            pylab.savefig(os.path.join(path, filename))
    return




'''-------------------------------------------------------------------------'''
# Old code from Clustering_Analysis.py
# Results plots


def results_plots():
    

    return()

def results_summary(path, input_file):
    '''Compute and plot summaries for clustering analysis'''

    results_file = os.path.join(path, input_file)

    # read in the data -- this is not an ideal way to do it since it requires
    # knowledge of file structure
    results_table = Table.read(results_file, format='ascii.commented_header',
                               guess=False)

    remove = []
    for i in range(0, len(results_table)):
        if results_table['score'][i] == 'N/A':
            remove.append(i)
    results_table.remove_rows(remove)

    num_clust = results_table['n_clust']
    score = results_table['score']
    total_obj = results_table['total_objects'].astype('float')

    # add a column with the size of the smallest cluster
    # have to do some tricky stuff since column corresponding to smallest
    # cluster varies dep on number of clusters

    results_table.add_column(Column(name='size_smallest',
                                    data=np.zeros(len(results_table)),
                                    dtype=np.int16))

    for i in range(0, len(results_table)):
        lastcol = 'c_{}'.format(num_clust[i])
        results_table['size_smallest'][i] = results_table[lastcol][i]

    # compute fraction of objects in largest cluster
    biggest_clust_fract = results_table['c_1']/total_obj
    # compute fraction of objects in smallest cluster
    smallest_clust_fract = results_table['size_smallest']/total_obj

    fig = plt.figure(figsize=(12,5))
    ax = fig.add_subplot(121)
    ax.scatter(num_clust, score)
    ax.set_xlabel('Number of clusters')
    ax.set_ylabel('Score')
    ax.set_title('Number of Clusters vs Silhouette Score', fontsize=11)

    ax = fig.add_subplot(122)
    ax.scatter(biggest_clust_fract, score, c='r', marker='o',
               label='Largest Cluster')
    ax.scatter(smallest_clust_fract, score, c='b', marker='o',
               label='Smallest Cluster')
    ax.legend(loc='best')
    ax.set_xlabel('Fractional size')
    ax.set_ylabel('Score')
    ax.set_title('Silhouette Score vs. Fractional Size of Largest/Smallest Cluster', fontsize=11)

    filename = 'silhouette_score_plots.png'
    pylab.savefig(os.path.join(path, filename))

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
                ax.set_title('Proportional Cluster Size')
                ax.set_xlabel('Number of Trials')
                ax.set_ylabel('Fractional Cluster Size')
                ax.xaxis.set_major_locator(plt.MultipleLocator(1))
            filename = '{}_clusters_vs_FractionalSize.png'.format(i)
            pylab.savefig(os.path.join(path, filename))

    return()
    
def inputs():

    inputs = argparse.ArgumentParser(description='Plotting functions')
    inputs.add_argument("data_file", help="Specify data file")
    inputs.add_argument("-pp", "--plot_path", help="Save directory",
                        default="results")
    inputs.add_argument("-p", "--plots", help="Select plots to make",
                        choices=['sco_v_clust', 'bw_v_clust', '', 'cvc',
                                 'wvc'], nargs='*')
    # Data trimming
    inputs.add_argument("-r", "--ratio", help="Set noise to signal ratio",
                        default=0.1)
    inputs.add_argument("-dp", "--data_points", help="Set objects from survey",
                        default=10000)

    criteria = inputs.parse_args()
    plotting(criteria.data_file, criteria.plot_path, criteria.plots,
             criteria.ratio, criteria.data_points)
    return

if __name__ == "__main__":
    inputs()