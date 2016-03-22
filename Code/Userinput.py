# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:19:08 2015

@author: Owner
"""
import argparse

def user_input():
    
    inputs = argparse.ArgumentParser(description = "Inputs for clustering analysis.")
    inputs.add_argument("data_file", help="Choose the data file for analysis", default = [])
    inputs.add_argument("analysis", help = "Choose the methods you would like to use", 
                        choices=['meanshift', 'kmeans', 'mst'], nargs='*', default = [])
                        
    inputs.add_argument("-kmi", "--kmeans_input", help="Choose the number of clusters input for kmeans", 
                        choices=['meanshift', 'experiments.txt'], default = [])
    inputs.add_argument("-p","--plots", help = "Choose the plots you would like to make", 
                        choices =['meanshift', 'kmeans_color', 'kmeans_xy', 'mst'], nargs='*', default = [])
    inputs.add_argument("-id", "--id_list", help = "Produces object id list", choices = ['yes','no'], default='no')
    inputs.add_argument("-rs", "--results_summary", help="Produces results summary", choices =['yes','no'], default="no")
    inputs.add_argument("-fn", "--file_names", help="Specify input and ouput file names. Default: experiments.txt, results.txt", default = ['experiments.txt, results.txt'], nargs='*')
    inputs.add_argument("-s", "--save_results", help="Enter path for saving output files", default="no")
                                                                
    criteria = inputs.parse_args()    
                
    if 'meanshift' in criteria.analysis: 
        print 'Success!!'
    if 'kmeans' in criteria.analysis:
        print 'Success!!' 
    if 'mst' in criteria.analysis:
        print 'Success!!'
            
    return()

if __name__ == "__main__": 
    user_input()
    