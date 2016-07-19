README File for the RESULTS directories of the m83_clustering project.
All programs used to generate results are kept in m83_clustering/Code directory.

INPUT FILES
These files are used to define the colour combinations as input for the various .py files:
clustering.py: experiments.txt - change band combinations for each colour and parameters for clustering techniques.
	- Parameters are toggled from command line by -kmi -bwi -afp arguments (estimate/experiments.txt)
plots.py: plots.txt - change band combinations for each colour
wave_statistics.py: stats_experiments.txt - change bands for colour combination
For each .txt file, band1 & band2 are colour 1, band3 & band4 are colour2.

DIRECTORIES

results/band_positions:
- plots of the xy coordinates of various cuts in the data set

results/band_statistics:

- /id_files: list of object ids and magnitude for each wavelength and aperture
	- Code/wave_statistics.py {band_id(path)}. Run from python console.

- /mag_uncertainty_vs_unc-cut: histograms of magnitude distributions and the distribution with the uncertainty cut imposed ontop
	- Code/double_histogram.py {mag_hist_cut(path)} {unc_hist_cut(path)}. Run from python console.

- /uncertainty_dist_w_histogram: distribution of mag vs. uncertainty with histograms on shared axes
	- Code/wave_histogram.py {wave_unc_hist(path)}. Run from python console.

- /uncertainty-3_vs_05: unc vs. mag distribution of both apertures sharing an axis
	- Code/uncertainty_plots.py {unc_plot(path)}. Run from python console.

-band_unc_limit_statistics.txt: object counts for each band at various uncertainty limits
	- Code/wave_statistics.py {band_unc_limit(path)}. Run from python console.

-filter_detections.txt: detections per filter without uncertainty limit. Only removed -99 values.
	- Code/wave_statistics.py {n_detections(path)}. Run from python console.

-filter_statistics.txt: statistics calculated for each filter in both apertures (mean, median, standard deviation etc.)
	- Code/wave_statistics.py {band_stats(path)}. Run from python console.

-flag_objects.txt: id of objects flagged from original Chandar catalogue.
	- Since only 9 flagged objects, done by hand no code. Objects are removed from data_v2.txt on.

-mag3_555_col_statistics.txt: statistics calculated for all possible colour combinations with the 555 filter in 3 aperture. Repeated for 0.5 aperture: mag05_555_col_statistics.txt
	- Code/wave_statistics.py {col_stats(path)}. Run from python console.

-mag3_555_unc_limit_statistics.txt: number of detections at various uncertainty limits for colour combinations with the 555 filter. Repeated for 0.5 aperture: mag05_555_unc_limit_statistics.txt
	- Code/wave_statistics.py {colour_unc(path)}. Run from python console. 
 
results/broad_band:
- location for all of the results of the [broad-broad, broad-broad] clustering

results/broad_narrow: 
location for all the results of the [broad-narrow, broad-broad] clustering
/[broad_narrow]: each directory contains the results of the [broad-narrow] colour with all possible [broad-broad] combinations

- /[broad_narrow]/band_band_plts: plots of the band1 vs band2 of each colour
	- Code/plots.py {band_v_band}. Run from command line: python plots.py data_v3.txt -pp (path) -p bvb 

- /[broad_narrow]/clustering: directories for each of the possible [broad-broad] combinations for the given [broad-narrow] combination
- /[broad_narrow]/clustering/inputs.txt: last clustering run for this combination
	- Code/clustering_2d.py. Run from command line: python clustering_2d.py data_v3.txt -rp (path) -pp (path) -a (clusterings) -p (plots) -kmi (Kmeans input) -bwi (MS input) -afp (affinity input) -wr (output results files: Y/N)

- /[broad_narrow]/clustering/broad-narrow_broad-broad: plots of all the clusterings done for this combination
	- cluster_statistics.txt: file containing the statistics for each cluster of each clustering performed for this combination
		- Code/clustering_2d.py. Run from command line: python clustering_2d.py data_v3.txt -rp (path) -pp (path) -a (clusterings) -p (plots) -kmi (Kmeans input) -bwi (MS input) -afp (affinity input) -wr (output results files: Y/N)
	- 05aperture_results_2d.txt: file containing summary of each clustering 
		- Code/clustering_2d.py. Run from command line: python clustering_2d.py data_v3.txt -rp (path) -pp (path) -a (clusterings) -p (plots) -kmi (Kmeans input) -bwi (MS input) -afp (affinity input) -wr (output results files: YES)

- /[broad_narrow]/clustering/broad-narrow_broad-broad/results: plots of differnet relationships between clustering parameters and cluster metrics 
	- n_clusters_vs_Fractional_Size.png: bar graph of the fractional size of each cluster in a clustering for each trial with that n_clusters
		- Code/results.py {bar_graph(results_table, n_clust, path)}: Run from python console.
	- inertia_plot.png: plot of the sum of squares for all k-means trials for a combination
		- Code/results.py {inertia_fluctuations(results_table, path)}: Run from python console.
	- meanshift_parameters.png: plot of the bandwidth parameter against the silhouette score and n_clusters for all trials
		Code/results.py {bandwidth_vs_score(results_table, path)}: Run from python console.
	- silhouette_score.png: plot of the silhouette score against the number of clusters and fractional size of the largest and smallest cluster
		- Code/results.py {silhouette_vs_nclust(results_table, path)}: Run from python console.
	- cluster_percentage.txt: file containing the cumulative fractional size of each cluster in a clustering
		- Code/results.py {object_cluster_distribution(results_table, n_clust, path)}: Run from python console.

- /[broad_narrow]/colour_colour: plots of the colour-colour distributions for the combinations
	- Code/plots.py {colour_v_colour}. Run from command line: python plots.py data_v3.txt -pp (path) -p cvc

- /[broad_narrow]/colour_histogram: histogram of the colour distribution for the combinations
	- Code/plots.py {colour_histogram}. Run from command line: python plots.py data_v3.txt -pp (path) -p c_hist
 
- /[broad_narrow]/wave_colour_plts: plots of the colour-magnitude distributions for each band in combination
	- Code/plots.py {wave_v_colour}. Run from command line: python plots.py data_v3.txt -pp (path) -p wvc

- /[broad_narrow]/wave_histogram: histograms of magnitude distribution for each wave in combinations
	- Code/plots.py {wave_histogram}. Run from command line: python plots.py data_v3.txt -pp (path) -p w_hist

- /[broad_narrow]/inputs.txt: last trial run in plots.py

results/old_trials: 
- all old combinations of colours (not useful as they shared a broad band between them)

results/test_scale: 
- directory to test any new clusterings or code that produces results

					