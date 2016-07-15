README File for the RESULTS directories of the m83_clustering project.

results/.txt files and figures:
- various results plots that will be moved to each colour directory (fractional size, score relationships, parameter relationships)
- final_col_statistics.txt: colour statistics for all colours determined during the week of July 15th (subject to change as a "final" usually does)

results/band_positions:
- plots of the xy coordinates of various cuts in the data set

results/band_statistics:
- /id_files: list of object ids and magnitude for each wavelength and aperture
- /mag_uncertainty_vs_unc-cut: histograms of magnitude distributions and the distribution with the uncertainty cut imposed ontop
- /uncertainty_dist_w_histogram: distribution of mag vs. uncertainty with histograms on shared axes
- /uncertainty-3_vs_05: unc vs. mag distribution of both apertures sharing an axis
- .txt files: various files containing: object counts for uncertainties
					detections per filter
					statistics per filter
					flagged objects from data_V1.txt(m83_clustering\Code)
					colour statistics with _555 base 
results/broad_band:
- location for all of the results of the [broad-broad, broad-broad] clustering

results/broad_narrow: 
- location for all the results of the [broad-narrow, broad-broad] clustering
- /[broad_narrow]: each directory contains the results of the [broad-narrow] colour with all possible [broad-broad] combinations
- /[broad_narrow]/band_band_plts: plots of the bands of each colour
- /[broad_narrow]/clustering: directories for each of the possible [broad-broad] combinations for the given [broad-narrow] combination
- /[broad_narrow]/clustering/inputs.txt: last clustering run for this combination
- /[broad_narrow]/clustering/broad-narrow_broad-broad: plots of all the clusterings done for this combination
						       cluster_statistics.txt: file containing the statistics for each cluster of each clustering performed for this combination
- /[broad_narrow]/colour_colour: plots of the colour-colour distributions for the combinations
- /[broad_narrow]/colour_histogram: histogram of the colour distribution for the combinations
- /[broad_narrow]/wave_colour_plts: plots of the colour-magnitude distributions for each band in combination
- /[broad_narrow]/wave_histogram: histograms of magnitude distribution for each wave in combinations
- /[broad_narrow]/inputs.txt: last trial run in plots.py

results/old_trials: 
- all old combinations of colours (not useful as they shared a broad band between them)

results/test_scale: 
- directory to test any new clusterings or code that produces results

					