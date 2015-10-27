from astropy.table import Table, Column
import numpy as np
import os

#usage:
#table_proc.process_tab(ned_in, ned_out, type_col='Type',select_list = table_proc.within_galaxy_types_ned, rfmt_fn=table_proc.rfmt_ned)
#table_proc.process_tab(simbad_in, simbad_out, type_col='OTYPE_S',select_list = table_proc.within_galaxy_types_simbad,rfmt_fn=table_proc.rfmt_simbad)
#(match output of these with TOPCAT)
#matchtab = Table.read('ned_simbad_match.fits')
#matchtab2 = table_proc.process_match(matchtab,'ned_simbad_match2.fits')


within_galaxy_types_ned = ['*Cl','HII','PofG','Neb','SN','SNR', 'V*','WR*']
background_types_ned = ['G','GClstr','GGroup','QSO']
obs_types_ned = ['IrS','RadioS','UvES','UvS', 'VisS','XrayS']

within_galaxy_types_simbad = ['**','Assoc*','Candidate_SN*','Candidate_WR*','Cepheid','Cl*','HII','ISM','LMXB','MolCld','Nova','PN','PartofG','SN','SNR','SNR?','Star','ULX','V*','WR*','semi-regV*']
background_types_simbad = ['BLLac','ClG','EmG','Galaxy','GinCl','GroupG','Possible_G','StarburstG']
obs_types_simbad = ['EmObj','IR','Radio','Radio(sub-mm)','UV','X']



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

def rfmt_ned(in_tab):
    ''' strip prefixes from some NED names to better match SIMBAD
         obvs. very specific to this dataset'''
    in_tab['Object Name'] = np.char.replace(in_tab['Object Name'],"MESSIER 083:","")
    in_tab['Object Name'] = np.char.replace(in_tab['Object Name'],"NGC 5236:","")
    in_tab['Object Name'] = np.char.replace(in_tab['Object Name'], " ", "") # strip spaces
    in_tab.rename_column('Object Name','Name_NED')
    return(in_tab)

def rfmt_simbad(in_tab):
    ''' change some object types to better match NED
         obvs. very specific to this dataset'''
    in_tab['Type'] = np.char.strip(in_tab['Type'],"?") # only applies to "SNR?"
    in_tab['Type'][in_tab['Type']=='Cl*'] = '*Cl' 
    in_tab.rename_column('MAIN_ID','Name_SIMBAD')
    in_tab['Name_SIMBAD'] = np.char.replace(in_tab['Name_SIMBAD'],"M83-","")
    in_tab['Name_SIMBAD'] = np.char.replace(in_tab['Name_SIMBAD']," ", "") # strip spaces
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
