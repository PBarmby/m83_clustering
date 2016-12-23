from astropy.Table import table
from glob import glob
import os.file

# there is a pylatex package that might do a lot of the heavy lifting here..

texstart = "\begin{document}"

def doit_2d(outfile = 'summary'):
    # generate blank LaTeX doc
    texdoc = outfile+'.tex'
    if os.file.is_file(texdoc):
        os.unlink(texdoc)

    outfile = file.open(texdoc, "w")
    outfile.write(texstart)
    
    # find the color-color plots
    kmeans_plots = glob('kmeans_color_*')
    
    # insert into a LaTex document

    # read the clustering results table and insert that too
    
    # find the results plots

    # insert into a LaTex document

    # finish up the document
    outfile.write(texend)
    outfile.close()

    # turn the LaTex into a PDF
    os.sysstr('pdflatex %s' % outfile)
    
    # all done!
    return

def twofigs(file1, file2):
    fig_string = "\begin{figure}\n \includegraphics[0.45\linewidth]{}\includegraphics[0.45\linewidth]{}"
    fig_string += "\end{figure}"
    fig_string += "\caption{Left: {}, Right: {}}" % use filename as caption?
    return(fig_string)

def onefigs(file1):
    fig_string = "\begin{figure}\n \includegraphics[0.45\linewidth]{}\includegraphics[0.45\linewidth]{}"
    fig_string += "\end{figure}"
    fig_string += "\caption{Left: {}, Right: {}}"
    return(fig_string)


