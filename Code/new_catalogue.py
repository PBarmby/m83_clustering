'''Write new catalogue without FLAGS and corrected 657 filter and corrected UNC
    - data_v2.txt has corrected 657 filter and removed flags 
    - data_v3.txt uses data_V2.txt to correct the mag_unc:
        - narrow filters: v2_unc/15
        - broad filters: v2_unc/10'''
import numpy as np
from astropy.table import Table
import os
import os.path


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
        for i in range (0, len(data)):
            if data['mag05_' + narrow[b]][i] != -99.0 and data['mag05_' + narrow[b] + '_unc'][i] != -99.0:
                data['mag05_' + narrow[b] + '_unc'][i] = data['mag05_' + narrow[b] + '_unc'][i]/narrow_correction
            if data['mag3_' + narrow[b]][i] != -99.0 and data['mag3_' + narrow[b] + '_unc'][i] != -99.0:
                data['mag3_' + narrow[b] + '_unc'][i] = data['mag3_' + narrow[b] + '_unc'][i]/narrow_correction
    Table.write(data, 'data_v3-test.txt', format='ascii.commented_header')
    return()
