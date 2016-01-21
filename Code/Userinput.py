# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:19:08 2015

@author: Owner
"""

def main(): 
    
    print "**Please enter a space between each input**" 
    analysis = raw_input("What analysis would you like to perform? (meanshift, kmeans, mst): ")
    make_plots = raw_input("What plots would you like(cluster, color, xy, mst): ")
    classification_ID = raw_input("Would you like the objects to be catalogued (Yes/No): ")
    make_results_summary = raw_input("Would you like a results summary (Yes/No): ")
    
    if analysis in ['kmeans']:
        meanshift_as_input_to_KMEANS = raw_input("Would you like to use the number of clusters estimated by MEANSHIFT as input to KMEANS? (Yes/No): ")
    
    #Save the inputs 
    
    inputs = [analysis, make_plots, classification_ID, make_results_summary]
    print inputs[0]
    #Check Inputs

    for i in range (0,3):
    
        while inputs[i] not in ['meanshift', 'kmeans', 'mst', 'cluster', 'color', 'xy', 'mst', 'Yes', 'yes', 'No', 'no']: 
       
            print "Invalid input, please try again."
            print "**Please enter a space between each input**"
            analysis = raw_input("What analysis would you like to perform? (meanshift, kmeans, mst): ")
            make_plots = raw_input("What plots would you like(cluster, color, xy, mst): ")
            classification_ID = raw_input("Would you like the objects to be catalogued (Yes/No): ")
            make_results_summary = raw_input("Would you like a results summary (Yes/No): ")
    
    confirm = raw_input("Start clustering now? (Yes/Change my inputs): ")
    
    while confirm not in ['Yes', 'yes']:
    
        print "**Please enter a space between each input**" 
        analysis = raw_input("What analysis would you like to perform? (meanshift, kmeans, mst): ")
        make_plots = raw_input("What plots would you like(cluster, color, xy, mst): ")
        classification_ID = raw_input("Would you like the objects to be catalogued (Yes/No): ")
        make_results_summary = raw_input("Would you like a results summary (Yes/No): ")
        
        confirm = raw_input("Start clustering now? (Yes/Change my inputs): ")
    

    
    return




