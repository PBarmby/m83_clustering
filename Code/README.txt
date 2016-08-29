# m83_clustering/Code directory README
# Contains all code written for figure generation, clustering, and computations

------------- .txt files -----------------

3d_experiments.txt 					input for: clustering_3d.py	File for selecting bands to cluster in three dimensions
all_labels.txt						File containing the magnitude data for all potential labels from cat_data/clustering_labels directory
data_test.txt 						First 10-20 rows of data_v1.txt used for testing code
data_v1.txt 						Original data file without flags removed or zero point/uncertainty corrections - DO NOT USE FOR CLUSTERING OR PLOT GENERATION
data_v2.txt 						generated with: new_catalogue.py	data_v1.txt with f657n 0.5/3 zero point corrected and flags removed - DO NOT USE FOR CLUSTERING OR PLOT GENERATION
data_v3.txt 						generated with: new_catalogue.py	data_v2.txt with all filter uncertanties corrected (BROADunc / 10, NARROWunc / 15) - USE FOR CLUSTERING AND PLOT GENERATION
experiments.txt 					input for: clustering_2d.py	File for selecting bands to cluster in two dimensions
plots.txt							input for: plots.py clustering_plots.py	File for selecting bands to create colour, CMD, histogram plots in two dimensions / recreate clustering plots in two dimensions with specific method and parameters
plots_3d.txt 						input for: clustering_plots.py 	File for selecting bands to recreate clustering plots in three dimensions
run.log							output file for all sharcnet submitted jobs
stats_experiments.txt					input for: wave_statistics.py	Used to select two bands to compute statistics on

------------- .py files ------------------

3d_plots.py						creates interactive 3d plot. Bands are set in file, run from Spyder
affinity_propagation.py				perform affinity propagation clustering and produces interactive plots. Bands are set in file, run from Spyder. Set number of colours by changing np.logical_and() lines.
Cluster.py							OLD DO NOT USE. original clustering code from Data Mining in Astronomy textbook. Performs K-Means clustering and produces plots from Data Mining in Astronomy textbook.
clustering_2d.py 					working clustering code for 2 dimensions. Run from command terminal. inputs: experiments.txt 	outputs: select from command line -h to list commands
clustering_3d.py					working clustering code for 3 dimensions, same operation as clustering_2d.py
									*Make sure to change the Global Variables and Function Variables outlined in IMPORTANT INFO at the top of the file
									*Make sure to change the base_path so the figures will save properly
Clustering_Analysis.py				OLD DO NOT USE. original clustering code from Data Mining in Astronomy textbook capable of running from terminal.
clustering_plots.py					creates plots from id_ files generated from clustering_2d.py and clustering_3d.py.	inputs: plots.txt or plots_3d.txt	Run from: Spyder terminal
compare_classes.py					Written by PB. not exactly sure what this does - I think it compares the NED and SINBAD data
compare_clustering.py 				creates bubble plots which compare the number of objects in each cluster in two clusterings. Run from Spyder terminal.
double_histogram.py					creates a histogram of all the magnitudes or uncertainties with a second histogram on top with an uncertainty limit.
extinction_correction.py				NOT WORKING. function to correct data file for Milky Way extinction.
hierarchical_clustering.py				NOT WORKING. Nathalie's code to create a heatmap after performing hierarchical clustering.
kmeans_test_multidimensional.py			performs k-means clustering on bands set in the function. Run from Spyder
labe_plot.py						creates colour-colour plots with different labels on top of data_v3.txt distribution. Inputs: plots.txt 	Run from Spyder
meanshift.py						performs meanshift on bands set in file
mult_clust_old.py					OLD DO NOT USE. first iteration of running multiple clusterings
multiplecluster.py					OLD DO NOT USE. iteration of running multiple clusterings
new_catalogue.py					creates a new catalogue from previous version of data_v(n).txt or create new label catalogue. Modify corrections to be made in mk_cat() function.
plot_models.py						plot colour models with data_v3.txt distribution. Set base colours in w1 - 6 variables.
plots.py							plot wave histogram, colour histogram, colour-colour, CMD, band v band. Input: plots.txt 	Run from: terminal -h for list of arguments
									* Make sure path is set properly for where you would like the plots to be saved (automatically adds on directory with colour name)
