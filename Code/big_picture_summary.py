from astropy.Table import table
from glob import glob
from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Matrix
from pylatex.utils import italic
import os

# combine figures and tables for a given clustering run into a single document

texstart = "\begin{document}"

def doit_2d(outfile = 'summary'):
    # generate blank LaTeX doc
    geometry_options = {"tmargin": "1cm", "lmargin": "2cm"}
    doc = Document(geometry_options=geometry_options)

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

    # find the color-color plots
    kmeans_plots = glob('kmeans_color_*')
    
    # insert into a LaTex document

# examples from documentation
    with doc.create(Subsection('Cute kitten pictures')):
        with doc.create(Figure(position='h!')) as kitten_pic:
            kitten_pic.add_image(image_filename, width='120px')
            kitten_pic.add_caption('Look it\'s on its back')



    # read the clustering results table and insert that too
    
    # find the results plots

    # insert into a LaTex document

    # generate document name, delete old versions
    texdoc = outfile+'.tex'
    if os.file.is_file(texdoc):
        os.unlink(texdoc)

    # turn the LaTex into a PDF
    doc.generate_pdf(outfile, clean_tex=False)
    
    # all done!
    return

def twofigs(file1, file2):
    fig_string = "\begin{figure}\n \includegraphics[0.45\linewidth]{}\includegraphics[0.45\linewidth]{}"
    fig_string += "\end{figure}"
    fig_string += "\caption{Left: {}, Right: {}}" % use filename as caption?
    return(fig_string)

def onefig(file1):
    fig_string = "\begin{figure}\n \includegraphics{}"
    fig_string += "\end{figure}"
    fig_string += "\caption{}"
    return(fig_string)


