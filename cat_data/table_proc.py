from astropy.table import Table
import numpy as np
import os

#usage:
#process_tab(ned_in, ned_out, type_col='Type',select_list = within_galaxy_types_ned)
#process_tab(simbad_in, simbad_out, type_col='OTYPE_S',select_list = within_galaxy_types_simbad)

within_galaxy_types_ned = ['*Cl','HII','PofG','Neb','SN','SNR', 'V*','WR*']
background_types_ned = ['G','GClstr','GGroup','QSO']
obs_types_ned = ['IrS','RadioS','UvES','UvS', 'VisS','XrayS']

within_galaxy_types_simbad = ['**','Assoc*','Candidate_SN*','Candidate_WR*','Cepheid','Cl*','HII','ISM','LMXB','MolCld','Nova','PN','PartofG','SN','SNR','SNR?','Star','ULX','V*','WR*','semi-regV*']
background_types_simbad = ['BLLac','ClG','EmG','Galaxy','GinCl','GroupG','Possible_G','StarburstG']
obs_types_simbad = ['EmObj','IR','Radio','Radio(sub-mm)','UV','X']



def process_tab(tab_in, tab_out, type_col, select_list = within_galaxy_types_ned):
    '''select specific object types from a table'''
    tab = Table.read(tab_in)
    tab['Type'] = np.char.strip(tab[type_col]) # remove whitespace -- helps with matching
    if type_col != 'Type': # change over to use Type for any table
        tab.delete_column(type_col)
    tg = tab.group_by('Type')
    mask = np.in1d(tg.groups.keys[type_col],select_list) # create mask for only wanted types
    wanted_types = tg.groups[mask]

    if tab_out != None:
        if os.path.exists(tab_out):
            os.unlink(tab_out)
        wanted_types.write(tab_out, format='fits')
    return(wanted_types)
    
