import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import pylab
import os
import os.path
import argparse

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
        if results_table['silhouette_score'][i] == 'N/A':
            remove.append(i)
    results_table.remove_rows(remove)

    num_clust = results_table['number_of_clusters']
    score = results_table['silhouette_score']
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
    ax.scatter(score, biggest_clust_fract, c='r', marker='o',
               label='Largest Cluster')
    ax.scatter(score, smallest_clust_fract, c='b', marker='o',
               label='Smallest Cluster')
    ax.legend(loc='upper left')
    ax.set_xlabel('Score')
    ax.set_ylabel('Fractional size')
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
                Y = results_table[colname][num_clust==i]/results_table['total_objects'][num_clust==i].astype('float')  # Compute percentage of objects in each cluster and graph
                ax.bar(X, Y, color=cluster_colours[n], bottom=yprev)
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