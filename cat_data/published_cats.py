from astropy.table import Table, Column, hstack, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
import numpy as np
import os

# published_cats:
#  tools for manipulating published catalogs of m83 objects, for comparison with Chandar et al catalog
#     - reformat and combine NED and SIMBAD data tables
#     - add in tables of data not included in NED or SIMBAD
#     - match to Chandar catalog to produce a list of best-matches with object types
#  

#usage:
#published_cats.ns_combine('ned-20160629.fits','simbad-20160629.fits','M83_NScomb.fits','M83_NSall.fits')
#published_cats.add_tables('M83_NSall.fits',['williams15_rsg.fits','kim12_wr.fits'],'M83_final.fits')
#published_cats.catalog_match('M83_final.fits', 'hlsp_wfc3ers_hst_wfc3_m83_cat_all_v2-corrected.txt','M83_ers_pubcat.txt')

# table columns to be renamed as part of processing: [I think (old, new)]
ned_rename = [('Name_N', 'Name'), ('RA(deg)', 'RA'), ('DEC(deg)', 'Dec'),('Type_N', 'Type')]
sim_rename = [('Name_S', 'Name'), ('RA_d', 'RA'), ('DEC_d', 'Dec'),('Type_S', 'Type')]

def ns_combine(ned_name, simbad_name, ns_combine, final_tab, match_tol = 1.0): # match_tol in arcsec

    ned_in = Table.read(ned_name)
    simbad_in = Table.read(simbad_name)
        
    # prelim processing
    ned_proc = reformat_cat(ned_in, old_name='Object Name', new_name='Name_N', old_type='Type', new_type='Type_N',
                            keepcols=['Object Name','RA(deg)','DEC(deg)','Type'])
    sim_proc = reformat_cat(simbad_in, old_name='MAIN_ID', new_name='Name_S', old_type='OTYPE', new_type='Type_S')

    # construct coordinates needed for matching
    ned_coo = SkyCoord(ra=ned_proc['RA(deg)'], dec=ned_proc['DEC(deg)'])
    sim_coo = SkyCoord(ra=sim_proc['RA_d'], dec=sim_proc['DEC_d']) 

    # do the matching
    matched_ned, matched_sim, ned_only, sim_only = symmetric_match_sky_coords(ned_coo, sim_coo, match_tol*u.arcsec)

    # generate the matched table
    # hstack is "horizontal stack", ie put the NED and SIMBAD columns for matched objects in the same row
    matchtab = hstack([ned_proc[matched_ned], sim_proc[matched_sim]], join_type = 'outer')
    
    # find the secure matches, save these as intermediate results
    matchtab2 = process_match(matchtab)
    matchtab2.write(ns_combine, format='fits')

    # process the secure match catalog
    keeplist = ['Name_N','RA(deg)','DEC(deg)','Type_N']
    matchtab3 = process_unmatch(Table(matchtab2[keeplist]), src='NS', rename_cols = ned_rename)

    #process the catalogs containing NED-only and SIMBAD-only objects
    nedcat = process_unmatch(ned_proc[ned_only], src = 'N', rename_cols= ned_rename)
    simcat = process_unmatch(sim_proc[sim_only], src = 'S', rename_cols = sim_rename)   
    
    # add the secure matches to the NED-only and SIMBAD-only catalogs
    # NB: I think this implies that the "insecure" matches just get thrown away - is that what we want?
    finaltab = vstack([matchtab3, nedcat, simcat],join_type = 'outer')
    
    # save the result
    finaltab.write(final_tab, format='fits')
            
    return

def add_tables(basetab_file, tab_file_list, outfile, jt = 'outer'):
   '''perform a vertical join of a list of tables to a base table
      input:
      basetab_file: filename for base table
      tab_file_list: list of filenames for additional tables
      outfile: name of output file
      jt: join type, should generaly be 'outer'
    '''    
    basetab = Table.read(basetab_file)
    tablist = [basetab]
    for filename in tab_file_list:
        tab = Table.read(filename)
        tablist.append(tab)
    stack = vstack(tablist, join_type= jt)
    stack.write(outfile)
    return

def catalog_match(pubcat_file, erscat_file, match_out_file, match_tol = 1.0):
    ''' matches combined NED/SIMBAD file to ERS source list
    '''
    pubcat = Table.read(pubcat_file, format = 'fits')
    erscat = Table.read(erscat_file, format='ascii.commented_header')

    # construct coordinates needed for matching
    pub_coo = SkyCoord(ra=pubcat['RA'], dec=pubcat['Dec'])
    ers_coo = SkyCoord(ra=erscat['ra']*u.degree, dec=erscat['dec']*u.degree) 

    # do the matching
#    closest_2to1, sep2d_2to1, sep3d = match_coordinates_sky(coord1, coord2) # location in coord2 for closest match to each coord1. len = len(coord1)
    closest, sep2d, sep3d = match_coordinates_sky(pub_coo, ers_coo) # location in coord2 for closest match to each coord1. len = len(coord1)
    matched  = sep2d < match_tol*u.arcsec
#    matched_ers, matched_pub, ers_only, pub_only = symmetric_match_sky_coords(ers_coo, pub_coo, match_tol*u.arcsec)

    # generate the matched table
    keeplist = ['id_','ra','dec']
    tmpcat = Table(erscat[keeplist])
    matchtab = hstack([tmpcat[closest][matched], pubcat[matched]], join_type = 'outer')

    # write the matched catalog to a file
    matchtab.write(match_out_file, format = 'ascii.commented_header')

    return

# notation to be replaced in object names: (old, new). Note this is galaxy-specific and really shoudn't be hardwired like this.
ns_replace_names = [(" ", ""), ("MESSIER083:",""),("NGC5236:",""), ("M83-",""), ("NGC5236","")]
# object types to be changed - makes NED more like SIMBAD
ns_replace_types = [('*Cl','Cl*'), ('PofG','Galaxy'),('X','XrayS'), ('Radio','RadioS')]
# objects to be removed from NED/SIMBAD list (e.g. the galaxy itself). Again galaxy-specific, shoudn't be hardwired like this.
ns_remove_ids= ['NAMENGC5236Group', 'M83', 'MESSIER083', 'NGC5236GROUP']
# names of RA & Dec columns - NED, SIMBAD respectively
ra_dec_cols = ['RA(deg)','DEC(deg)','RA_d','DEC_d']

def reformat_cat(in_tab, old_name, new_name, old_type, new_type, replace_names=ns_replace_names, replace_types=ns_replace_types, remove_id=ns_remove_ids, keepcols=None):
    ''' reformat NED or SIMBAD catalog to make more intercompatible'''     
    # only keep selected columns
    if keepcols!= None:
        in_tab = in_tab[keepcols]
        
    # change units for RA/Dec
    for col in ra_dec_cols:
        if col in in_tab.colnames:
            in_tab[col].unit = u.degree
    
    # change ID for name & type columns
    in_tab.rename_column(old_name, new_name)
    in_tab.rename_column(old_type, new_type)

    # reformat some object names: remove spaces, remove some parts of names that are repeated
    for pair in replace_names:
        in_tab[new_name] = np.char.replace(in_tab[new_name],pair[0], pair[1])

    # reformat some object types
    # remove spaces and question marks
    in_tab[new_type] = np.char.replace(in_tab[new_type],"?", "")
    in_tab[new_type] = np.char.replace(in_tab[new_type]," ", "")
    # rename some type indicators
    for pair in replace_types:
        in_tab[new_type][in_tab[new_type]==pair[0]] = pair[1]        
    
    # delete rows whose names are in list remove_id
    # there's a non-loopy way to do this but I can't remember it
    remove_idx = []
    for i in range(0,len(in_tab)):
        if in_tab[i][new_name] in remove_id:
            remove_idx.append(i)
    in_tab.remove_rows(remove_idx)
                                
    # all done
    return(in_tab)

def symmetric_match_sky_coords(coord1, coord2, tolerance):
    '''produce the symmetric match of coord1 to coord2, such that objects are only a match if
       distance is < tolerance AND
       each is the closest match in the other list  
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
        # (NB: should put in an assertion to check that tolerance is an angular unit)
        if sep2d_1to2[i] < tolerance and i == closest_2to1[closest_1to2[i]]:     
            index1_matched.append(i)
            index2_matched.append(closest_2to1[i])
        else:
            index1_unmatched.append(i)

    for j in range(0, len(coord2)):
        if j not in index2_matched:
            index2_unmatched.append(j)
                        
    return(index1_matched, index2_matched, index1_unmatched, index2_unmatched)
    
    
def process_match(matched_tab_in):
    '''take an already-matched table btw NED and SIMBAD, find secure matches
    where secure implies both name and object type match'''
    goodmatch = np.logical_or(matched_tab_in['Name_S']==matched_tab_in['Name_N'],
                              matched_tab_in['Type_S']==matched_tab_in['Type_N'])
    matched_tab_in.add_column(Column(goodmatch, name='Secure'))
    return(matched_tab_in)


def process_unmatch(tab_in, src, rename_cols):
    '''process unmatched catalogs: rename columns, add source column (ie, NED or SIMBAD)'''
    for pair in rename_cols:
        tab_in.rename_column(pair[0], pair[1])
    tab_in.add_column(Column(name='Source', length = len(tab_in), dtype='S2'))
    tab_in['Source'] = src
    return(tab_in)


# BELOW HERE IS OLD STUFF, not used    

within_galaxy_types_ned = ['*Cl','HII','PofG','Neb','SN','SNR', 'V*','WR*']
background_types_ned = ['G','GClstr','GGroup','QSO']
obs_types_ned = ['IrS','RadioS','UvES','UvS', 'VisS','XrayS']

within_galaxy_types_simbad = ['**','Assoc*','Candidate_SN*','Candidate_WR*','Cepheid','Cl*','HII','ISM','LMXB','MolCld','Nova','PN','PartofG','SN','SNR','SNR?','Star','ULX','V*','WR*','semi-regV*']
background_types_simbad = ['BLLac','ClG','EmG','Galaxy','GinCl','GroupG','Possible_G','StarburstG']
obs_types_simbad = ['EmObj','IR','Radio','Radio(sub-mm)','UV','X']

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
