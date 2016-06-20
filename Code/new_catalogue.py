'''Write new catalogue without FLAGS and corrected 657 filter'''
import numpy as np 
from astropy.table import Table
import os 
import os.path


def mk_cat(path):
    correction = 1.6
    data = Table.read('hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1a-corrected.txt', format='ascii.commented_header', guess=False)

    for i in range (0, len(data)):
        if data['mag05_657'][i] != -99: 
            data['mag05_657'][i] = data['mag05_657'][i] + correction
        if data['mag3_657'][i] != -99: 
            data['mag3_657'][i] = data['mag3_657'][i] + correction
    Table.write(data, 'data_v2.txt', format='ascii.commented_header')
    #write_catalogue(data, path)
    return()

def write_catalogue(catalogue, path_): 
    cat_path = make_directory(path_)
    '''Create file with list of object ID and cluster number'''
    file_name = 'data_v2.txt'
    file_path = '{}\\{}'.format(cat_path, file_name)
    header = catalogue.colnames
    if not os.path.exists(file_path):
        create_path = os.path.join(cat_path, file_name)
        new_catalogue = open(create_path, "a")
        new_catalogue.write('# ' + str(header) + '\n')
        new_catalogue.close()
    new_catalogue_path = os.path.join(cat_path, file_name)
    new_catalogue = open(new_catalogue_path, "a")
    new_catalogue.write(str(catalogue))
    new_catalogue.close()
    #id_path = os.path.join(cat_path, file_name)
    #id_tab = Table(data=[catalogue], names=[catalogue.colnames])
                   #names=['x', 'y', 'ra', 'dec', 'id_', 'ci_white', 'mag0.5_white', 'mag3_white', 'flux_white', 'msky_white', 'flag', 'mag05_225', 'mag05_225_unc mag3_225 mag3_225_unc mag05_336 mag05_336_unc mag3_336 mag3_336_unc mag05_373 mag05_373_unc mag3_373 mag3_373_unc mag05_438 mag05_438_unc mag3_438 mag3_438_unc mag05_487 mag05_487_unc mag3_487 mag3_487_unc mag05_502 mag05_502_unc mag3_502 mag3_502_unc mag05_555 mag05_555_unc mag3_555 mag3_555_unc mag05_657 mag05_657_unc mag3_657 mag3_657_unc mag05_673 mag05_673_unc mag3_673 mag3_673_unc mag05_814 mag05_814_unc mag3_814 mag3_814_unc
#])
    #id_tab.write(id_path, format='ascii.commented_header')
    
    return()

def make_directory(p_path):
    '''Save results of each analysis
            save_path: set path of where you would like results saved'''
    # Create new plots_folder
    pl_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(
              p_path)
    if not os.path.exists(pl_path):
        os.makedirs(pl_path)

    return(pl_path)
