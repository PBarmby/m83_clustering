import numpy as np 
import os
from astropy.table import Table
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
        
def stats(path_): 
    data = Table.read('data.txt', format='ascii.commented_header', guess=False)
    band_names = data.colnames
    stats_path = make_directory(path_)
    file_name = "filter_statistics.txt"
    file_path = os.path.join(stats_path, file_name)
    header = "# band wave_mean wave_median wave_std wave_var wave_min wave_max unc_mean unc_median unc_std unc_var unc_min unc_max"
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
        remove_bad_data = np.logical_and(band_data != -99, band_unc_data != -99)
        band_data_trim = band_data[remove_bad_data]
        band_unc_trim = band_unc_data[remove_bad_data]
        # Compute Magnitude Statistics
        mag_mean = np.mean(band_data_trim)
        mag_median = np.median(band_data_trim)
        mag_stdev = np.std(band_data_trim)
        mag_var = np.var(band_data_trim)
        mag_min = np.min(band_data_trim)
        mag_max = np.max(band_data_trim)
        #mag_bad_obs = remove_bad_data.count('False')
        # Create string to write
        mag_string = "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(mag_mean, mag_median, mag_stdev, mag_var, mag_min, mag_max)#, mag_bad_obs)
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
    