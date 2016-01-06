import numpy as np 
from astropy.table import Table, join
import matplotlib.pyplot as plt
from itertools import izip

def matched_tabs(tab1, tab2, tab_out, id1 = 'ID', id2='ID'):
    ''' match tab1 and tab2 on columns id1,id2
       return inner join with all columns'''

    # read data
    t1 = Table.read(tab1)
    t2 = Table.read(tab2)

    # rename ID columns if necessary
    if id1 != 'ID':
        t1.rename_column(id1, 'ID')
    if id2 != 'ID':
        t2.rename_column(id2, 'ID')

    # check for unique IDs
    if len(unique(t1, keys='ID')) != len(t1):
       print 'Table 1 has some non-unique IDs, quitting'
       return
    if len(unique(t2, keys='ID')) != len(t2):
       print 'Table 2 has some non-unique IDs, quitting'
       return

    # do the join and return the matched table
    match_tab = join(t1, t2, keys = 'ID')
    return(match_tab)

def plot_classes(match_tab, class1, class2, class1_lab = 'NED', class2_lab='K-means'):
    # group data by class1
    grouped_tab = match_tab.group_by(class1)

    class1_types = grouped_tab.groups.keys
    class1_ind = range(0,class1_types)

    # make plot
    fig,ax = plt.subplots()
# needs fix, see http://docs.astropy.org/en/stable/table/operations.html#grouped-operations
    for key, group in izip(class1_types, grouped_tab.groups):
        ax.plot(c1, group[c2])
    ax.set_xlabel(class1_lab)
    ax.set_ylabel(class2_lab)
    return
