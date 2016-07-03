from astropy.table import Table, Column, hstack, vstack
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

def go(ned_name, simbad_name, combine_out, match_tol = 1.0): # match_tol in arcsec

    ned_in = Table.read(ned_name)
    simbad_in = Table.read(simbad_name)
        
    # prelim processing
    ned_proc = reformat_cat(ned_in, old_name='Object Name', new_name='Name_N', old_type='Type', new_type='Type_N')
    sim_proc = reformat_cat(simbad_in, old_name='MAIN_ID', new_name='Name_S', old_type='OTYPE', new_type='Type_S')

    # construct coordinates needed for matching
    ned_coo = SkyCoord(ra=ned_proc['RA(deg)'], dec=ned_proc['DEC(deg)'])
    sim_coo = SkyCoord(ra=sim_proc['RA_d'], dec=sim_proc['DEC_d']) 

    # do the matching
    matched_ned, matched_sim, ned_only, sim_only = symmetric_match_sky_coords(ned_coo, sim_coo, match_tol*u.arcsec)

    # generate the matched table
    matchtab = hstack([ned_proc[matched_ned], sim_proc[matched_sim]], join_type = 'outer')

    # add on the unmatched objects
    matchtab2 = vstack([matchtab, ned_proc[ned_only], sim_proc[sim_only]],join_type = 'outer')

    # mark the really good matches
    finaltab = process_match(matchtab2)

    # save the result
    finaltab.write(combine_out, format='fits')
            
    return



ns_replace_names = [("MESSIER 083:",""),("NGC 5236:",""), ("M83-",""), ("NGC5236","")(" ", "")]
ns_replace_types = [('*Cl','Cl*'), ('PofG','Galaxy'),('X','XRayS'), ('Radio','RadioS')]
ns_remove_ids= ['NAMENGC5236Group', 'M83', 'MESSIER083', 'NGC5236GROUP']
ra_dec_cols = ['RA(deg)','DEC(deg)','RA_d','DEC_d']

def reformat_cat(in_tab, old_name, new_name, old_type, new_type, replace_names=ns_replace_names, replace_types=ns_replace_types, remove_id=ns_remove_ids):
    ''' reformat NED or SIMBAD catalog to make more intercompatible'''     
    # change units for RA/Dec
    for col in ra_dec_cols:
        if col in in_tab.colnames:
            in_tab[col].unit = u.degree
    
    # change ID for name & type columns
    in_tab.rename_column(old_name, new_name)
    in_tab.rename_column(old_type, new_type)

    # reformat some object names
    for pair in replace_names:
        in_tab[new_name] = np.char.replace(in_tab[new_name],pair[0], pair[1])

    # reformat some object types
    in_tab[new_type] = np.char.replace(in_tab[new_type],"?", "")
    for pair in replace_types:
        in_tab[new_type][in_tab[new_type]==pair[0]] = pair[1]        
    
    # delete rows whose names are in remove_id
    # there's a non-loopy way to do this but I can't remember it
    remove_idx = []
    for i in range(0,len(in_tab)):
        if in_tab[i][new_name] in remove_id:
            remove_idx.append(i)
    in_tab.remove_rows(remove_idx)
                                
    # all done
    return(in_tab)

def symmetric_match_sky_coords(coord1, coord2, tolerance):
    '''produce the symmetric match of coord1 to coord2
       output:
       index1_matched: index into coord1 for matched objects
       index2_matched: index into coord2 for matches of objects in index1_matched
       index1_unmatch: indices for unmatched objects in coord1
       index2_unmatch: indices for unmatched objects in coord2
    '''
    closest_2to1, sep2d_2to1, sep3d = match_coordinates_sky(coord1, coord2) # location in coord2 for closest match to each coord1. len = len(coord1)
    closest_1to2, sep2d_1to2, sep3d = match_coordinates_sky(coord2, coord1) # location in coord1 for closest match to each coord2. len = len(coord2)

    index1_matched = []
    index2_matched = []
    index1_unmatched = []
    index2_unmatched = []

    for i in range(0, len(coord1)): # doubtless there is a more Pythonic way to do this..
        # not sure this condition covers all of the possible corner cases. But good enough.
        if sep2d_1to2[i] < tolerance and closest_1to2[i] == closest_2to1[closest_1to2[i]]:     
            index1_matched.append(i)
            index2_matched.append(closest_2to1[i])
        else:
            index1_unmatched.append(i)

    for j in range(0, len(coord2)):
        if j not in index2_matched:
            index2_unmatched.append(j)
                        
    return(index1_matched, index2_matched, index1_unmatched, index2_unmatched)
    
    
def process_match(matched_tab_in):
    '''find secure matches btw NED and SIMBAD'''
    goodmatch = np.logical_or(matched_tab_in['Name_S']==matched_tab_in['Name_N'],
                              matched_tab_in['Type_S']==matched_tab_in['Type_N'])
    matched_tab_in.add_column(Column(goodmatch, name='Secure'))
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
