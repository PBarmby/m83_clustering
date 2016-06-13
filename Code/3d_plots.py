import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from astropy.table import Table

data = Table.read('data.txt', format='ascii.commented_header', guess=False)
aperture = '05'
band1 = 'mag{}_555'.format(aperture)
band2 = 'mag{}_814'.format(aperture)
band3 = 'mag{}_336'.format(aperture)
band4 = 'mag{}_438'.format(aperture)
'''Select data for analysis'''
ratio = 0.25
data = data#[10000:30000]
# Colour 1
wave1 = data[band1]
wave1_unc = data[band1+'_unc']
wave2 = data[band2]
wave2_unc = data[band2+'_unc']
# Colour 2
wave3 = data[band3]
wave3_unc = data[band3+'_unc']
wave4 = data[band4]
wave4_unc = data[band4+'_unc']
# Change parameters to match data_file
# Remove data pieces with no value

wave1_trim = np.logical_and(wave1 != -99, wave1_unc != -99)
wave2_trim = np.logical_and(wave2 != -99, wave2_unc != -99)
wave3_trim = np.logical_and(wave3 != -99, wave3_unc != -99)
wave4_trim = np.logical_and(wave4 != -99, wave4_unc != -99)
colour1_ratio = np.logical_and(wave1_unc < ratio,
                               wave2_unc < ratio)
colour2_ratio = np.logical_and(wave3_unc < ratio,
                               wave4_unc < ratio)

gooddata1 = np.logical_and(np.logical_and(wave1_trim, wave2_trim),
                           np.logical_and(wave3_trim, wave4_trim))

# Remove data above certain magnitude
gooddata2 = np.logical_and(colour1_ratio, colour2_ratio)

# Only data that match criteria for both colours
greatdata = np.logical_and(gooddata1, gooddata2)

colour1 = wave1[greatdata] - wave2[greatdata]
colour2 = wave3[greatdata] - wave4[greatdata]
x = data['x'][greatdata]
y = data['y'][greatdata]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(0, len(colour1)):
    if colour1[i] < 0:
        ax.scatter(x[i], colour1[i], y[i], c='b', marker='o')
    else:
        ax.scatter(x[i], colour1[i], y[i], c='r', marker='o')

ax.set_xlabel('X Pixel')
ax.set_ylabel(band1 + ' - ' + band2)
ax.set_zlabel('Y Pixel')

plt.show()
