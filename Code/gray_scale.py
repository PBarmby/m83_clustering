import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

fname = 'C:\\Users\\Alex\\Documents\\GitHub\\m83_clustering\\Paper\\figs\\final_colour_figures\\mag05_336-mag05_438_mag05_438-mag05_555_mag05_438-mag05_814\\meanshift_3d_projection_4cl_F336W-F438WvsF438W-F814W'

img = Image.open(fname+'.png').convert('LA')
img.save(fname+'_greyscale.png')