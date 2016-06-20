import numpy as np 
import os
from astropy.table import Table
from matplotlib import pyplot as plt
from itertools import cycle
'''Data is now trimmed with unc < 0.2 and removed all -99s
    - Exception: ID file which records all object detections'''

def band_unc_limit(path_):
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    band_names = data.colnames
    unc_limit_path = make_directory(path_)
    file_name = "band_unc_limit_statistics.txt"
    file_path = os.path.join(unc_limit_path, file_name)
    header = "# band limit n_objects"
    if not os.path.exists(file_path):
        unc_limit_file = open(file_path, "a")
        unc_limit_file.write(header + '\n')
        unc_limit_file.close()
    for i in range (11, len(data.colnames), 2):
        unc_limit_file = open(file_path, "a")
        band = band_names[i]
        band_mag = data[band]
        band_unc = data[band+'_unc']
        trim = np.logical_and(band_mag != -99, band_unc != -99)
        band_mag_trim = band_mag[trim]
        band_unc_trim = band_unc[trim]
        for k in range(10, 100, 10):
            limit = float(k)/100
            gooddata = band_mag_trim[band_unc_trim < limit]
            limit_string = "{} {:.4f} {}".format(band, limit, len(gooddata))
            unc_limit_file.write(limit_string + '\n')
        unc_limit_file.close()
    return()


def colour_unc(path_):
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    trials = Table.read('stats_experiments.txt',
                        format='ascii.commented_header', guess=False)
    col_limit_path = make_directory(path_)
    file_name = "{}-colour_unc_limit_statistics.txt".format(trials['band2'][0])
    file_path = os.path.join(col_limit_path, file_name)
    header = "# band1 band2 limit n_objects"

    if not os.path.exists(file_path):
        col_unc_file = open(file_path, "a")
        col_unc_file.write(header + '\n')
        col_unc_file.close()

    for i in range(0, len(trials)):
        col_unc_file = open(file_path, "a")
        # Load Band1
        band1 = trials['band1'][i]
        band1_mag = data[band1]
        band1_unc = data[band1 + '_unc']
        # Remove Band1 no detections
        band1_trim = np.logical_and(band1_mag != -99, band1_unc != -99)
        # Load Band2
        band2 = trials['band2'][i]
        band2_mag = data[band2]
        band2_unc = data[band2 + '_unc']
        # Remove Band2 no detections
        band2_trim = np.logical_and(band2_mag != -99, band2_unc != -99)
        # Remove no detections in one band
        colour_trim = np.logical_and(band1_trim, band2_trim)
        band1_mag_clean = band1_mag[colour_trim]
        band2_mag_clean = band2_mag[colour_trim]
        band1_unc_clean = band1_unc[colour_trim]
        band2_unc_clean = band2_unc[colour_trim]
        # Compute Color and uncertainty
        col_mag = band1_mag_clean - band2_mag_clean
        col_unc = np.sqrt(band1_unc_clean*band1_unc_clean+band2_unc_clean*band2_unc_clean)
        for k in range(1, 10):
            limit = float(k)/10
            gooddata = col_mag[col_unc < limit]
            limit_string = "{} {} {} {}".format(band1, band2, limit,
                                                len(gooddata))
            col_unc_file.write(limit_string + "\n")
        col_unc_file.close()
    return


def band_id(path_):
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    band_names = data.colnames
    id_path = make_directory(path_)
    #header = "# object_id magnitude"
    for i in range(11, len(band_names), 2):
        band = band_names[i]
        wave = data[band]
        object_id = data['id_']
        file_name = "{}_object_id.txt".format(band)
        file_path = os.path.join(id_path, file_name)
        #if not os.path.exists(file_path):
         #   id_file = open(file_path, "a")
          #  id_file.write(header + '\n')
           # id_file.close()
        #id_file = open(file_path, "a")
        observations = object_id[wave != -99]
        id_tab = Table(data=[observations, wave[wave != -99]],
                       names=['object_id', 'magnitude'])
        id_tab.write(file_path, format='ascii.commented_header')        
    return()


def col_stats(path_): 
    '''Not filtered based on uncertainty - only removed -99'''
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    trials = Table.read('stats_experiments.txt',
                        format='ascii.commented_header', guess=False)
    col_path = make_directory(path_)
    file_name = "05_col_statistics.txt"
    file_path = os.path.join(col_path, file_name)
    header = "# band1 band2 colour_mean colour_median colour_std colour_var colour_min colour_max num_obj unc_mean unc_median unc_std unc_var unc_min unc_max"
    if not os.path.exists(file_path):
        col_file = open(file_path, "a")
        col_file.write(header + '\n')
        col_file.close()
    for i in range(0, len(trials)):
        col_file = open(file_path, "a")
        band1 = trials['band1'][i]
        band1_mag = data[band1]
        band1_unc = data[band1+'_unc']
        band2 = trials['band2'][i]
        band2_mag = data[band2]
        band2_unc = data[band2+'_unc']
        remove_bad = np.logical_and(np.logical_and(band1_mag != -99, band1_unc != -99),
                              np.logical_and(band2_mag != -99, band2_unc != -99))
        limit = np.logical_and(band1_unc < 0.2, band2_unc < 0.2)
        trim = np.logical_and(remove_bad, limit)
        band1_mag_trim = band1_mag[trim]
        band1_unc_trim = band1_unc[trim]
        band2_mag_trim = band2_mag[trim]
        band2_unc_trim = band2_unc[trim]
        colour_mag = band1_mag_trim - band2_mag_trim
        colour_unc = np.sqrt(band1_unc_trim*band1_unc_trim+band2_unc_trim*band2_unc_trim)
        # Compute Magnitude Statistics
        mag_mean = np.mean(colour_mag)
        mag_median = np.median(colour_mag)
        mag_stdev = np.std(colour_mag)
        mag_var = np.var(colour_mag)
        mag_min = np.min(colour_mag)
        mag_max = np.max(colour_mag)
        mag_obj = len(colour_mag)
        # Create string to write
        mag_string = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {}".format(mag_mean, mag_median, mag_stdev, mag_var, mag_min, mag_max, mag_obj)#, mag_bad_obs)
        # Compute Uncertainty Statistics
        unc_mean = np.mean(colour_unc)
        unc_median = np.median(colour_unc)
        unc_stdev = np.std(colour_unc)
        unc_var = np.var(colour_unc)
        unc_min = np.min(colour_unc)
        unc_max = np.max(colour_unc)
        unc_string = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(unc_mean, unc_median, unc_stdev, unc_var, unc_min, unc_max)
        col_file.write(band1 + ' ' + band2 + ' ' + mag_string + ' ' + unc_string + '\n')
        col_file.close()
    return


def band_stats(path_): 
    '''Not filtered based on uncertainty - only removed -99 - now removed <0.2 unc'''
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    band_names = data.colnames
    stats_path = make_directory(path_)
    file_name = "filter_statistics-with_limit.txt"
    file_path = os.path.join(stats_path, file_name)
    header = "# band wave_mean wave_median wave_std wave_var wave_min wave_max num_obj unc_mean unc_median unc_std unc_var unc_min unc_max"
    if not os.path.exists(file_path):
        stats_file = open(file_path, "a")
        stats_file.write(header + '\n')
        stats_file.close()
    for i in range (11, len(band_names), 2):
        stats_file = open(file_path, "a")
        band = band_names[i]
        band_data = data[band_names[i]]
        band_unc_data = data[band_names[i+1]]
        # Remove bad data 
        remove_bad_data = np.logical_and(np.logical_and(band_data != -99,
                                                        band_unc_data != -99),
                                                        band_unc_data < 0.2)
        band_data_trim = band_data[remove_bad_data]
        band_unc_trim = band_unc_data[remove_bad_data]
        # Compute Magnitude Statistics
        mag_mean = np.mean(band_data_trim)
        mag_median = np.median(band_data_trim)
        mag_stdev = np.std(band_data_trim)
        mag_var = np.var(band_data_trim)
        mag_min = np.min(band_data_trim)
        mag_max = np.max(band_data_trim)
        mag_obj = len(band_data_trim)
        # Create string to write
        mag_string = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {}".format(mag_mean, mag_median, mag_stdev, mag_var, mag_min, mag_max, mag_obj)#, mag_bad_obs)
        # Compute Uncertainty Statistics
        unc_mean = np.mean(band_unc_trim)
        unc_median = np.median(band_unc_trim)
        unc_stdev = np.std(band_unc_trim)
        unc_var = np.var(band_unc_trim)
        unc_min = np.min(band_unc_trim)
        unc_max = np.max(band_unc_trim)
        #unc_bad_obs = remove_bad_data.count('False')
        # Create string to write
        unc_string = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(unc_mean, unc_median, unc_stdev, unc_var, unc_min, unc_max)#, unc_bad_obs)
        stats_file.write(band + ' ' + mag_string + ' ' + unc_string + '\n')
        stats_file.close()
    return


def make_directory(p_path):
    '''Save results of each analysis
            save_path: set path of where you would like results saved'''
    # Create new plots_folder
    pl_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(
              p_path)
    if not os.path.exists(pl_path):
        os.makedirs(pl_path)

    return(pl_path)
