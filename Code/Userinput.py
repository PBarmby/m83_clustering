# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:19:08 2015

@author: Owner
"""
import argparse

    
inputs = argparse.ArgumentParser()
inputs.add_argument("data_file", help="Choose the data file for analysis")
inputs.add_argument("analysis", help = "Choose the methods you would like to use", choices=['meanshift', 'kmeans', 'mst'])
inputs.add_argument("-msi", "--kmeans_input", help="Choose the number of clusters input for kmeans", choices=['meanshift', 'experiments.txt'])
inputs.add_argument("-p","--plots", help = "Choose the plots you would like to make", choices =['meanshift', 'kmeans_color', 'kmeans_xy', 'mst'], default ="none")
inputs.add_argument("-id", "--id_list", help = "Produces object id list", choices = ['yes','no'], default="no")
inputs.add_argument("-rs", "--results_summary", help="Produces results summary", choices =['yes','no'], default="no")
inputs.add_argument("-s", "--save_results", help="Enter path for saving output files", default="no")    
    
criteria = inputs.parse_args()    

if 'meanshift' in criteria.analysis: 
    print 'meanshift' 
elif 'kmeans' in criteria.analysis:
    print 'kmeans' 
elif 'mst' in criteria.analysis:
    print 'mst' 
else: 
    print criteria.data_file
    





