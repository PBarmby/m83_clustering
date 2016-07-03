from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
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

def go(ned_in, simbad_in, combine_out, match_tol = 1.0): # match_tol in arcsec
    
    # prelim processing
    ned_proc = reformat_cat(ned_in, old_name='Object Name', new_name='Name_N', old_type='Type', new_type='Type_N')
    sim_proc = reformat_cat(simbad_in, old_name='MAIN_ID', new_name='Name_S', old_type='OTYPE_S', new_type='Type_S')

    # construct coordinates needed for matching
    ned_coo = SkyCoord(ra=ned_proc['RA(deg)']*u.degree, dec=ned_proc['DEC(deg)']*u.degree)
    sim_coo = SkyCoord(ra=sim_proc['RA_d']*u.degree, dec=sim_proc['DEC_d']*u.degree) 

    idx_sim, sep2d_sim, sep3d = match_coordinates_sky(ned_coo, sim_coo) # idx is location in sim_coo for closest match to each ned_coo
    idx_ned, sep2d_ned, sep3d = match_coordinates_sky(sim_coo, ned_coo) # idx is location in ned_coo for closest match to each sim_coo


    matchtab = match(ned_proc, simb_proc)
    matchtab2 = process_match(matchtab, combine_out)
    
    return



ns_replace_names = [("MESSIER 083:",""),("NGC 5236:",""), ("M83-",""), (" ", "")]
ns_replace_types = [("SNR?","SNR"), ("Cl*?", "*Cl"), ('*Cl','Cl*')]
ns_remove_ids= ['NAMENGC5236Group', 'M83', 'MESSIER083', 'NGC5236GROUP']

def reformat_cat(in_tab, old_name, new_name, old_type, new_type, replace_names=ns_replace_names, replace_types=ns_replace_types, remove_id=ns_remove_ids):
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
    in_tab[new_type] = np.char.strip(tab[new_type]) # remove whitespace -- helps with matching
    
    # delete rows whose names are in remove_id
    # there's a non-loopy way to do this but I can't remember it
    remove_idx = []
    for i in range(0,len(in_tab)):
        if in_tab[i][new_name] in remove_id:
            remove_idx.append(i)
    in_tab.remove_rows(remove_idx)
                                
    # all done
    return(in_tab)

def process_match(matched_tab_in, matched_tab_out=None):
    '''find secure matches btw NED and SIMBAD'''
    goodmatch = np.logical_or(matched_tab_in['Name_SIMBAD']==matched_tab_in['Name_NED'],
                              matched_tab_in['Type_S']==matched_tab_in['Type_N'])
    matched_tab_in.add_column(Column(goodmatch, name='Secure'))
    if matched_tab_out != None: # write to file
        if os.path.exists(matched_tab_out):
            os.unlink(matched_tab_out)
        matched_tab_in.write(matched_tab_out, format='fits')
    return(matched_tab_in)

# not used    
def name_match(simbad_name, ned_name):
    matched = np.zeros(len(simbad_name),dtype='bool')
    matched[np.char.replace(simbad_name," ", "")==np.char.replace(ned_name," ", "")] = True
    ## TODO: account for cases where one name has leading zero in an ID (eg [WLC2001]03) and other doesn't

    return(matched)

# not using this
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
