'''Compare the results of two classifications'''
import os
import numpy as np
import pylab
from matplotlib import pyplot as plt
from astropy.table import Table, join


def compare_class(df_1, df_2, d1_path, d2_path, save_path, plots):
    survey = Table.read('data_v3.txt', format='ascii.commented_header',
                        guess=False)
    # Load two clusterings
    path1 = load_path(d1_path)
    path2 = load_path(d2_path)
    save_dir = mk_dir(save_path)  # make the saving directory
    label_data, clust_1data, clust_2data = data(df_1, df_2, survey, path1,
                                                path2)

    if "corr" in plots:
        corr_plot(label_data, save_dir)
    if "broad_pos" in plots:
        broad_plot(clust_1data, clust_2data, save_dir)

    return


def data(clustering1, clustering2, survey_data, path_1, path_2, plots):
    clust1 = Table.read(path_1 + clustering1, format='ascii.commented_header',
                        guess=False)
    clust2 = Table.read(path_2 + clustering2, format='ascii.commented_header',
                        guess=False)
    clust1_data = join(clust1, survey_data, 'id_')
    clust2_data = join(clust2, survey_data, 'id_')
    clust_labels = join(clust1, clust2, 'id_')

    return(clust_labels, clust1_data, clust2_data)


def mk_dir(path_):
    if not os.path.exists(path_):
        os.makedirs(path_)
    return


def load_path(path_):
    path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(path_)
    return(path)


def corr_plot(labels, path_):
    unique_labels = np.unique(labels['cluster_number_1'])
    n_clust = len(unique_labels)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(0, n_clust):
        for k in range(0, n_clust):
            objects = len(labels[np.logical_and(labels['cluster_number_1'] == i,
                                                labels['cluster_number_2'] == k)])
            ax.scatter(i+1, k+1, s=objects)
            ax.annotate(str(objects), xy=(i+1, k+1), xytext=(i+1.05, k+1.05))
            ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
            ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
    plt.show()
    pylab.savefig(path_ + 'similar_clusters.png')
    plt.close()
    return


def broad_plot():
    return

