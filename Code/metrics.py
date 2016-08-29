import numpy as np
from matplotlib import pyplot as plt
from numpy import var as var
from numpy import sqrt as sqrt


def sd_index(data_in_clust, survey_data, n_clust, cluster_centers):
    dis = []
    clust_var = np.arange(0, n_clust)
    data_var = abs(var(survey_data))

    for s in range(0, n_clust):
        clust_var[s] = abs(var(data_in_clust))

    scat = 1.0 / n_clust * sum(clust_var) / data_var

    for i in range(0, n_clust):
        for j in range(0, n_clust):
            if j != i:
                dis.append(distance(cluster_centers[i], cluster_centers[j]))
    dis_max = max(dis)
    dis_min = min(dis)
    dis_all = dis_max / (dis_min * sum(dis))
    
    sd_indx = scat*dis_m + dis_all

    return(sd_index)


def distance(x, y):
    dist = sqrt((x[1] - x[0])**2 + (y[1] - y[0] ** 2))
    return(dist)
