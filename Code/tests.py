# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 14:49:49 2016

@author: Owner
"""
import os 
import os.path
import numpy as np
import fileinput

from matplotlib import pyplot as plt
from astropy.table import Table, Column

def hello():
    print "Hello World" 
    
    return 
    
def save_analysis(path):
    '''Save results of each analysis
            working_path: set path of Clustering_Analysis folder
            save_path: set path of where you would like results saved'''
            
    name_of_file = "test.txt"
    save_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(path) 
    file_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}\\{}'.format(path, name_of_file)
    header = '# clustering band1 band2 band3 band4 number_of_clusters'
    
    completeName = os.path.join(save_path, name_of_file) 
       
    if not os.path.exists(save_path):
        os.makedirs(save_path)
        if not os.path.exists(file_path):
            file1 = open(completeName, "a")
            file1.write(header + '\n')
            file1.close()
    
        
    #file_ = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}\\{}'.format(path, name_of_file)
    
    #if not os.path.exists(file_):
        
    
    # Create new analysis_folder
    # Save analysis results
         
    
    file1 = open(completeName, "a")
    file1.write("save results here"+'\n')
    file1.close()
    
    
       
    return()
    

    