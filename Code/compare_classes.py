# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:29:05 2016

@author: Owner
"""

import numpy as np 
from astropy.table import Table, join, unique
import matplotlib.pyplot as plt

def matched_tabs(tab1, tab2, id1 = 'ID', id2='ID'):
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
       print('Table 1 has some non-unique IDs, quitting')
       return
    if len(unique(t2, keys='ID')) != len(t2):
       print('Table 2 has some non-unique IDs, quitting')
       return

    # do the join and return the matched table
    match_tab = join(t1, t2, keys = 'ID')
    return(match_tab)

# although this works technically, it's not very easy to interpret visually
def plot_classes_v1(match_tab, class1, class2, class1_lab = 'NED', class2_lab='K-means'):
    # group data by class1
    grouped_tab = match_tab.group_by(class1)
    class1_types = grouped_tab.groups.keys[class1]
    ngroups_class1 = len(class1_types)
    ngroups_class2 = len(unique(match_tab,class2))

    # make plot
    fig,ax = plt.subplots()
    # loop over types
    for i, group in enumerate(grouped_tab.groups):
        ax.plot(0*group[class2]+i, group[class2], label=class1_types[i], ls='None',marker='s')
    ax.set_xlabel(class1_lab, fontweight='bold',fontsize=14)
    ax.set_ylabel(class2_lab, fontweight='bold',fontsize=14)
    ax.set_xlim(-1,ngroups_class1+1)
    ax.set_ylim(-1,ngroups_class2+1)
    ax.legend()
    fig.show()
    return

# intention: a set of bar charts, where the horizontal axis is class1
#  and the vertical axis counts the different classes in class2 
def plot_classes_v2(match_tab, class1, class2, class1_lab = 'NED', class2_lab='K-means'):
    # figure out how many different classes
    ngroups_class1 = len(unique(match_tab, class1))
    keys2 = unique(match_tab,class2)[class2]
    ngroups_class2 = len(keys2)
    n2in1 = np.zeros((ngroups_class1, ngroups_class2),dtype='int')              

    # group by class2, figure out how many of each class2 in each class1 group
    # surely there's a better way to do this, but I couldn't work it out
    grouped_tab = match_tab.group_by(class1)
    for i, grp in enumerate(grouped_tab.groups): 
        for j in range(0,ngroups_class2):
            n2in1[i,j] = len(grp[grp[class2]==keys2[j]]) 
            
    # make plot
    fig,ax = plt.subplots()
    c = ['b','g','r','c','y','m','0.5','k','blueviolet','chartreuse','crimson','darkslateblue','goldenrod','pink','brown','orange'] # set colour cycle
    xind = np.arange(ngroups_class1) # set number of classes

    # loop over class2
    for j in range(0, ngroups_class2):
        if j>0:
            bot += n2in1[:,j-1] # sum up the bottom of the bars
            ax.bar(xind, n2in1[:,j], bottom = bot, label=keys2[j],color=c[j])
        else:
            ax.bar(xind, n2in1[:,j], label=keys2[j],color=c[j])            
            bot = n2in1[:,j]*0
    ax.set_xlabel(class1_lab, fontweight='bold',fontsize=14)
    ax.set_ylabel('Number', fontweight='bold',fontsize=14)
    plt.xticks(xind + 0.4,grouped_tab.groups.keys[class1])
    ax.set_ylim(0,n2in1.sum(axis=1).max()*1.1)
    ax.legend()
    fig.show()

    return

