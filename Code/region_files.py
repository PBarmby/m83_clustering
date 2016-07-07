''' Make region files for different types of objects in M83_ers_pubcat.txt '''
'''Date: July 7, 2016'''
''' line 53: if 'WR*' in cat.. doesn't work'''
''' Types of objects:  #need to check these are correct definitions

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
    ''' Make a region file for each type of object '''
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
