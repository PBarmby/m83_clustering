''' Make region files for different types of objects in M83_ers_pubcat.txt '''
'''Date: July 7, 2016'''
''' line 53: if 'WR*' in cat.. doesn't work'''
''' Types of objects:  #need to check these are correct definitions
DS9 Colors: Black, white, red, green, blue, cyan, magenta, yellow

Cl*: star cluster
XrayS: Xray source
SNR: Supernova remnent
HII: H2 region
RadioS: radio source
WR*: 
VisS:
RSG: Red supergiant 
MolCld: Molecular cloud
Star: star
PN: 
IR: 
Neb: Nebula
Galaxy: galaxy
GroupG:
ISM: interstellar medium
Possible_G: 
GGroup: 
G:
PartofG: 
LMXB: low mass xray binary
ULX:
SN: supernova'''

import os, os.path
import numpy as np
from astropy.table import Table


def make_region_files(file_name, file_path):
    ''' Make a catalogue file for each type of object '''
    general_catalogue = load_catalogue_file(file_name, file_path)
    # Set key for individual file
    general_catalogue.sort('Type')
    catalogue_object_type = general_catalogue['Type']
    types_of_objects = array_of_type_names(catalogue_object_type)

    for i in range(0, len(types_of_objects)):
        print types_of_objects[i]
        if 'Cl*' not in types_of_objects[i]:
            Table.write(general_catalogue[general_catalogue['Type'] == types_of_objects[i]],
                        str(types_of_objects[i]) + '_catalogue.txt',
                        format='ascii.commented_header')
        if types_of_objects[i] == 'Cl*':
            Table.write(general_catalogue[general_catalogue['Type'] == types_of_objects[i]],
                        'starcluster_catalogue.txt',
                        format='ascii.commented_header')
        if types_of_objects[i] == 'WR*':
            Table.write(general_catalogue[general_catalogue['Type'] == types_of_objects[i]],
                        'wr_catalogue.txt', format='ascii.commented_header')
            
    return


def make_ds9_files(file_name, file_path):
    ''' Make ds9 compatibale region files'''
    catalogue = load_catalogue_file(file_name, file_path)
    object_id = catalogue['id_'] 
    data = Table.read('C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\Code\\data_v3.txt',
                      format='ascii.commented_header', guess=False)
    coordinates_x = data['ra']
    coordinates_y = data['dec']
    object_x_coordinate = np.arange(0, len(object_id), dtype=float)
    object_y_coordinate = np.arange(0, len(object_id), dtype=float)

    reg_file_name = '{}_ra_dec.reg'.format('starclust')
    reg_file_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(file_path)

    for i in range(0, len(object_id)):
        object_x_coordinate[i] = coordinates_x[data['id_'] == object_id[i]]
        object_y_coordinate[i] = coordinates_y[data['id_'] == object_id[i]]

    region_file_ = os.path.join(reg_file_path, reg_file_name)
    region_file = open(region_file_, "w")
    for w in range(0, len(object_id)):
        coordinate_string = "{:.10f},{:.10f},".format(object_x_coordinate[w],
                                                    object_y_coordinate[w])
        region_file.write("CIRCLE(" + coordinate_string + '15) # color = blue' + '\n')
    region_file.close()
    return


def ra_dec(file_name, file_path):
    '''Make ds9 compatable xy file with ra and dec coordinates'''
    catalogue = load_catalogue_file(file_name, file_path)
    object_type = catalogue['Type'][0]
    object_id = catalogue['id_'] 
    data = Table.read('C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\Code\\data_v3.txt', format='ascii.commented_header',
                      guess=False)
    coordinates_x = data['ra']
    coordinates_y = data['dec']
    object_x_coordinate = np.arange(0, len(object_id), dtype=float)
    object_y_coordinate = np.arange(0, len(object_id), dtype=float)

    reg_file_name = '{}_ra_dec.txt'.format('starcluster')
    reg_file_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}\\ra_dec'.format(file_path)

    for i in range(0, len(object_id)):
        object_x_coordinate[i] = coordinates_x[data['id_'] == object_id[i]]
        object_y_coordinate[i] = coordinates_y[data['id_'] == object_id[i]]

    region_file_ = os.path.join(reg_file_path, reg_file_name)
    region_file = open(region_file_, "w")
    for w in range(0, len(object_id)):
        coordinate_string = "{:.10f} {:.10f}".format(object_x_coordinate[w],
                                                    object_y_coordinate[w])
        region_file.write(coordinate_string + '\n')
    region_file.close()
    return()


def dbase_ra_dec(file_name, file_path):
    '''Make ds9 compatable xy file with ra and dec coordinates'''
    catalogue = load_catalogue_file(file_name, file_path)
    object_type = catalogue['Type'][0]
    object_id = catalogue['id_'] 
    data = Table.read('C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\Code\\data_v3.txt', format='ascii.commented_header',
                      guess=False)
    coordinates_x = catalogue['RA']
    coordinates_y = catalogue['Dec']
    reg_file_name = '{}_ra_dec.txt'.format('starcluster')
    reg_file_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}\\ra_dec'.format(file_path)
    region_file_ = os.path.join(reg_file_path, reg_file_name)
    ra_dec_tab = Table(data=[coordinates_x, coordinates_y])
    ra_dec_tab.write(region_file_, format='ascii.no_header')
    return()


def label_cat(file_name, file_path):
    ''' Make label catalogue'''
    catalogue = load_catalogue_file(file_name, file_path)
    object_id = catalogue['id_']
    object_type = catalogue['Type'][0]
    data = Table.read('C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\Code\\data_v3.txt',
                      format='ascii.commented_header', guess=False)

    reg_file_name = '{}_label_catalogue.txt'.format('starcluster')
    reg_file_path = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(file_path)
    region_file_ = os.path.join(reg_file_path, reg_file_name)
    region_file = open(region_file_, "w")

    for w in range(0, len(object_id)):
        region_file.write(str(data[data['id_'] == object_id[w]]))

    region_file.close()
    return


def load_catalogue_file(name, path):
    ''' Load in general catalogue '''
    path_ = 'C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\{}'.format(path)
    cat_data = Table.read(os.path.join(path_, name),
                          format='ascii.commented_header', guess=False)
    return(cat_data)


def array_of_type_names(type_):
    t = []
    for i in range(1, len(type_)):
        if i == 1: 
            t.append(type_[i])
        if type_[i] != type_[i-1]:
            t.append(type_[i])
    return(t)
