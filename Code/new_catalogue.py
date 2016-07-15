'''Write new catalogue without FLAGS and corrected 657 filter and corrected UNC
    - data_v2.txt has corrected 657 filter and removed flags 
    - data_v3.txt uses data_V2.txt to correct the mag_unc:
        - narrow filters: v2_unc/15
        - broad filters: v2_unc/10'''
import numpy as np
from astropy.table import Table
import os
import os.path


def label_cat():
    label_titles = ['PN', 'wr']
    for i in range(0, len(label_titles)):
        data = Table.read('data_v3.txt', format='ascii.commented_header',
                          guess=False)
        object_cat = Table.read('C:\\Users\\Owner\\Documents\\GitHub\\m83_clustering\\cat_data\\m83_ers_pubcat_by_object\\' + label_titles[i] + '_catalogue.txt',
                                format='ascii.commented_header', guess=False)
        remove = []
        object_id = object_cat['id_']
        for r in range(0, len(data)):
            if data['id_'][r] not in object_id:
                remove.append(r)
        data.remove_rows(remove)
        label_table = Table(data=data)
        label_table.write(label_titles[i] + '_labels.txt',
                          format='ascii.commented_header')
    return


def mk_cat():
    broad_correction = 10
    narrow_correction = 15
    data = Table.read('data_v2.txt', format='ascii.commented_header',
                      guess=False)
    broad = ['225', '336', '438', '555', '814']
    narrow = ['373', '487', '502', '657', '673']

    # Broad uncertanty correction
    for b in range(0, len(broad)):
        for i in range (0, len(data)):
            if data['mag05_' + broad[b]][i] != -99.0 and data['mag05_' + broad[b] + '_unc'][i] != -99.0:
                data['mag05_' + broad[b] + '_unc'][i] = data['mag05_' + broad[b] + '_unc'][i]/broad_correction
            if data['mag3_' + broad[b]][i] != -99.0 and data['mag3_' + broad[b] + '_unc'][i] != -99.0:
                data['mag3_' + broad[b] + '_unc'][i] = data['mag3_' + broad[b] + '_unc'][i]/broad_correction

    # Narrow uncertanty correction
    for n in range(0, len(narrow)):
        for j in range (0, len(data)):
            if data['mag05_' + narrow[n]][j] != -99.0 and data['mag05_' + narrow[n] + '_unc'][j] != -99.0:
                data['mag05_' + narrow[n] + '_unc'][j] = data['mag05_' + narrow[n] + '_unc'][j]/narrow_correction
            if data['mag3_' + narrow[n]][j] != -99.0 and data['mag3_' + narrow[n] + '_unc'][j] != -99.0:
                data['mag3_' + narrow[n] + '_unc'][j] = data['mag3_' + narrow[n] + '_unc'][j]/narrow_correction
    Table.write(data, 'data_v3.txt', format='ascii.commented_header')
    return()
