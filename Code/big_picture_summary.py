from astropy.Table import table
from glob import glob
from pylatex import Document, Section, Subsection, Tabular, Math,  Figure 
from pylatex.utils import NoEscape
from StringIO import StringIO
import os

# combine figures and tables for a given clustering run into a single document

def doit_kmeans_3d(outfile = 'summary',
                   resfile = '05aperture_results_3d.txt',
                   clust_stat='3d_cluster_statistics.txt'):

    # generate blank LaTeX doc
    geometry_options = {"tmargin": "1cm", "lmargin": "2cm"}
    doc = Document(geometry_options=geometry_options)

    # insert some introductory material
    dirname = os.getcwd()
    doc.append('K-means clustering output, found in')
    doc.append(dirname)

    # get the results table and insert into the document
    tmptex = get_results(resfile)
    doc.append(tmptex)

    doc.append(NewPage())
    
    for nclust in range(3,4):
        # generate a new section
        section_title = 'K-means colour-colour plots, K=%d' % nclust

        # find the relevant plots
        glob_str = 'kmeans*_%1dcl*' % nclust
        kmeans_plots = glob(glob_str)

        # extract the relevant parts of clust_stat  - NOT COMPLETE
        tmptex = get_stats(stattab)
         
        # insert table and plots into the LaTex document
        with doc.create(Section(section_title)):
            doc.append(tmptex)
            for image_filename in kmeans_plots:
                with doc.create(Figure()) as fig:
                    fig.add_image(image_filename)
                    fig.add_caption(image_filename)
            doc.append(NewPage())


    # turn the LaTex into a PDF
    doc.generate_tex(filepath=outfile)
#    doc.generate_pdf(outfile, clean_tex=False)
    
    # all done!
    return

collist= ['n_clust','inertia', 'score', 'total_objects', 'c_1', 'c_2', 'c_3', 'c_4', 'c_5', 'c_6', 'c_7', 'c_8']
newnames = ['Nclust','inertia', 'score', 'TotObj', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']
def get_results(resfile):
    # read the results table and convert it to LaTeX
    res_tab = Table.read(resfile, format='ascii.commented_header')
    res_tab_k = res_tab[res_tab['clustering']=='kmeans'][collist]

    # rename some columns to get rid of underscores
    for i,col in enumerate(res_tab_k.colnames):
        if col != newnames[i]:
            res_tab_k.rename_column(col,newnames[i])

    # write to a string that pylatex can use
    tmptex = StringIO()
    res_tab_k.write(tmptex,format='latex')
    tmpstr = NoEscape(tmptex.getvalue())
    tmptex.close()
    return(tmpstr)


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
