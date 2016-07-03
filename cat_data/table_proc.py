from astropy.table import Table, Column
import numpy as np
import os

#Reformat and combine NED and SIMBAD data tables
# reformat involves making object types and naming schemes match
# combine: match on RA/dec, then check for match on name/type
#  

#usage:
#table_proc.go(ned_in, simbad_in, combine_out)

within_galaxy_types_ned = ['*Cl','HII','PofG','Neb','SN','SNR', 'V*','WR*']
background_types_ned = ['G','GClstr','GGroup','QSO']
obs_types_ned = ['IrS','RadioS','UvES','UvS', 'VisS','XrayS']

within_galaxy_types_simbad = ['**','Assoc*','Candidate_SN*','Candidate_WR*','Cepheid','Cl*','HII','ISM','LMXB','MolCld','Nova','PN','PartofG','SN','SNR','SNR?','Star','ULX','V*','WR*','semi-regV*']
background_types_simbad = ['BLLac','ClG','EmG','Galaxy','GinCl','GroupG','Possible_G','StarburstG']
obs_types_simbad = ['EmObj','IR','Radio','Radio(sub-mm)','UV','X']

def go(ned_in, simbad_in, combine_out):
    
    # prelim processing
    # TODO: do we actually want to just select within-glx types?
    ned_proc = process_tab(ned_in, type_col='Type',select_list = within_galaxy_types_ned, rfmt_fn=rfmt_ned)
    simb_proc = process_tab(simbad_in, type_col='OTYPE_S',select_list = within_galaxy_types_simbad,rfmt_fn=rfmt_simbad)

    #(match output of these)
    matchtab = match(ned_proc, simb_proc)
    matchtab2 = process_match(matchtab, combine_out)
    
    return


def process_tab(tab_in, tab_out, type_col, select_list = within_galaxy_types_ned, rfmt_fn=None):
    '''select specific object types from a table'''
    tab = Table.read(tab_in)
    if type_col != 'Type':
        tab.rename_column(type_col, 'Type')
    tab['Type'] = np.char.strip(tab['Type']) # remove whitespace -- helps with matching
    tg = tab.group_by('Type')
    mask = np.in1d(tg.groups.keys['Type'],select_list) # create mask for only wanted types
    wanted_types = tg.groups[mask]

    # do some reformatting
    if rfmt_fn != None:
        wanted_types = rfmt_fn(wanted_types)

    if tab_out != None: # write to file
        if os.path.exists(tab_out):
            os.unlink(tab_out)
        wanted_types.write(tab_out, format='fits')
    return(wanted_types)


ned_replace_names = [("MESSIER 083:",""),("NGC 5236:",""), (" ", "")]
simbad_replace_names = [("M83-",""), (" ", "")]

ned_replace_types = [('*Cl','Cl*')]
simbad_replace_types = [("SNR?","SNR"), ("Cl*?", "*Cl")]

def reformat_cat(in_tab, old_name, new_name, old_type, new_type, replace_names, replace_types, remove_id):
    ''' reformat NED or SIMBAD catalog to make more intercompatible'''

    # change ID for name & type columns
    in_tab.rename_column(old_name, new_name)
    in_tab.rename_column(old_type, new_type)

    # reformat some object names
    for pair in replace_name:
        in_tab[new_name] = np.char.replace(pair[0], pair[1])

    # reformat some object types
    for pair in replace_types:
        in_tab[new_type] = np.char.replace(pair[0], pair[1])        

    # delete rows whose names are in remove_id
    idx = in_tab[new_name] in remove_id
    in_tab.remove_rows(idx)
    
    # all done
    return(in_tab)
        


def process_match(matched_tab_in, matched_tab_out=None):
    '''find secure matches btw NED and SIMBAD'''
    goodmatch = np.logical_or(name_match(matched_tab_in['Name_SIMBAD'],matched_tab_in['Name_NED']),
                              name_match(matched_tab_in['Type_1'],matched_tab_in['Type_2']))
    matched_tab_in.add_column(Column(goodmatch, name='Secure'))
    if matched_tab_out != None: # write to file
        if os.path.exists(matched_tab_out):
            os.unlink(matched_tab_out)
        matched_tab_in.write(matched_tab_out, format='fits')
    return(matched_tab_in)

def name_match(simbad_name, ned_name):
    matched = np.zeros(len(simbad_name),dtype='bool')
    matched[np.char.replace(simbad_name," ", "")==np.char.replace(ned_name," ", "")] = True
    ## TODO: account for cases where one name has leading zero in an ID (eg [WLC2001]03) and other doesn't

    return(matched)
