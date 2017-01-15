from astropy.table import Table
from glob import glob
from pylatex import Document, Section, Subsection, Tabular, Math,  Figure, NewPage 
from pylatex.utils import NoEscape
#from StringIO import StringIO # Python 2.7
from io import StringIO # Python 3.x - need to correctly deal with this, probably by getting rid of py2.7..
import os
import numpy as np

# combine figures and tables for a given clustering run into a single document
# TODO: generalize to 2D
#  - resfile 05aperture_results_2d.txt 
#  - statfile='3d_cluster_statistics.txt'
# 
def doit_3d(outfile = 'summary',
                   resfile = '05aperture_results_3d.txt',
                   statfile='3d_cluster_statistics.txt'):

    # generate blank LaTeX doc
    geometry_options = {"tmargin": "1cm", "lmargin": "2cm"}
    doc = Document(geometry_options=geometry_options)

    # insert some introductory material
    dirname = os.getcwd()
    doc.append('Clustering output, found in')
    doc.append(dirname)

    # get the K-means section of results table and insert into the document
    with doc.create(Section('K-means results')):
        tmptex, ncl_list = get_stats(resfile, 'kmeans', oldcols_results, newcols_results, None, None)
        doc.append(tmptex)
        doc.append(NewPage())
    
        # add the individual K-means runs
        add_sections_to_doc(ncl_list, 'kmeans', doc, statfile)

    # now the intro material for meanshift
    # get the results table and insert into the document
    with doc.create(Section('Meanshift results')):
        tmptex, ncl_list = get_stats(resfile, 'meanshift', oldcols_res_ms, newcols_res_ms, None, None)
        doc.append(tmptex)
        doc.append(NewPage())

        # add the individual meanshift runs
        add_sections_to_doc(ncl_list, 'meanshift', doc, statfile)
            
    # turn the LaTex into a PDF
    doc.generate_tex(filepath=outfile)
#    doc.generate_pdf(outfile, clean_tex=False)
    
    # all done!
    return

# fix hard-wiring here?
oldcols_results = ['n_clust','inertia', 'score', 'total_objects', 'c_1', 'c_2', 'c_3', 'c_4', 'c_5', 'c_6', 'c_7', 'c_8']
newcols_results = ['Nclust','inertia', 'score', 'TotObj', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']

oldcols_res_ms = ['b_width', 'n_clust','score',  'total_objects', 'c_1', 'c_2', 'c_3', 'c_4', 'c_5', 'c_6', 'c_7', 'c_8','c_9','c_10']
newcols_res_ms = ['Bandw', 'Nclust','score', 'TotObj', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9','N10']

oldcols_2d_stats= ['clust_num','n_obj','t_scr','c_scr','rms','avg_dist','max_dist','min_dist','stdev','cen_1','cen_2','avg_col_1','avg_col_2']
newcols_2d_stats = ['Cluster','Nobj','tScore','cScore','rms','AvgDist','MaxDist','MinDist','Stdev','Cen1','Cen2','AvgCol1','AvgCol2']

oldcols_3d_stats= ['clust_num','n_obj','t_scr','c_scr','rms','avg_dist','max_dist','min_dist','stdev','cen_1','cen_2','cen_3','avg_col_1','avg_col_2','avg_col_3']
newcols_3d_stats = ['Cluster','Nobj','tScore','cScore','rms','AvgDist','MaxDist','MinDist','Stdev','Cen1','Cen2','Cen3','AvgCol1','AvgCol2','AvgCol3']

# extracts columns from text output files
def get_stats(statsfile, cluster_alg, oldcols, newcols, sel_cond2=None, sel_val2=None):

    # read the results table and select the right rows
    res_tab = Table.read(statsfile, format='ascii.commented_header')
    if sel_cond2 != None:
        sel_cond = np.logical_and(res_tab['clustering']== cluster_alg, res_tab[sel_cond2]==sel_val2)
    else:
        sel_cond = res_tab['clustering']== cluster_alg
    res_subtab = res_tab[sel_cond][oldcols]

    # figure out the number of clusters
    if sel_cond2 == None:
        ncl = res_subtab['n_clust'].data
        ncl_list = np.unique(ncl).tolist()
    else:
        ncl_list = []
        
    # rename some columns to get rid of underscores
    for i,col in enumerate(res_subtab.colnames):
        if col != newcols[i]:
            res_subtab.rename_column(col,newcols[i])

    # write the table to a string that pylatex can use
    tmptex = StringIO()
    res_subtab.write(tmptex,format='latex')
    tmpstr = NoEscape(tmptex.getvalue())
    tmptex.close()
    return(tmpstr, ncl_list)

def add_sections_to_doc(ncl_list, clustering_type, doc, statfile):
    for nclust in ncl_list:
        # generate a new subsection
        section_title = '%s results, Ncl=%d' % (clustering_type,nclust)

        # find the relevant plots
        glob_str = '%s*_%dcl*' % (clustering_type,nclust)
        all_plots = glob(glob_str)

        # extract the relevant parts of statfile and insert
        tmptex, junk = get_stats(statfile, clustering_type, oldcols_3d_stats, newcols_3d_stats, 'total_clust', nclust)

        # insert table and plots into the LaTex document
        with doc.create(Subsection(section_title)):
            doc.append(tmptex)
            for image_filename in all_plots:
                with doc.create(Figure()) as fig:
                    fig.add_image(image_filename)
                    fig.add_caption(image_filename)
            doc.append(NewPage())
    return

    
    # examples from documentation
#    with doc.create(Section('The simple stuff')):
#        doc.append('Some regular text and some')
#        doc.append(italic('italic text. '))
#        doc.append('\nAlso some crazy characters: $&#{}')
#        with doc.create(Subsection('Math that is incorrect')):
#            doc.append(Math(data=['2*3', '=', 9]))
#
#    with doc.create(Subsection('Table of something')):
#            with doc.create(Tabular('rc|cl')) as table:
#                table.add_hline()
#                table.add_row((1, 2, 3, 4))
#                table.add_hline(1, 2)
#                table.add_empty_row()
#                table.add_row((4, 5, 6, 7))

## examples from documentation
#    with doc.create(Subsection('Cute kitten pictures')):
#        with doc.create(Figure(position='h!')) as kitten_pic:
#            kitten_pic.add_image(image_filename, width='120px')
#            kitten_pic.add_caption('Look it\'s on its back')
