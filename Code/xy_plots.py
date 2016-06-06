import numpy as np 
from astropy.table import Table
from matplotlib import pyplot as plt 

data = Table.read('data.txt', format='ascii.commented_header', guess=False)
ratio = 30
band1 = 'mag05_336'
band2 = 'mag05_502' 
band3 = 'mag05_555'
band4 = 'mag05_657' 

wave1 = data[band1]
wave1_unc = data[band1+'_unc']
wave2 = data[band2]
wave2_unc = data[band2+'_unc']
wave3 = data[band3]
wave3_unc = data[band3+'_unc']
wave4 = data[band4]
wave4_unc = data[band4+'_unc']

trim1 = np.logical_and(wave1 != -99, wave1_unc != -99)
trim1_1 = np.logical_and(wave1/wave1_unc > ratio, trim1)
wave1_trim = wave1[trim1_1]

trim2 = np.logical_and(wave2 != -99, wave2_unc != -99)
trim2_2 = np.logical_and(wave2/wave2_unc > ratio, trim2)
wave2_trim = wave2[trim2_2]

trim3 = np.logical_and(wave3 != -99, wave3_unc != -99)
trim3_3 = np.logical_and(wave3/wave3_unc > 70, trim3)
wave3_trim = wave3[trim3_3]

trim4 = np.logical_and(wave4 != -99, wave4_unc != -99)
trim4_4 = np.logical_and(wave4/wave4_unc > ratio, trim4)
wave4_trim = wave4[trim4_4]

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(2, 2, 1)
ax1.scatter(data['x'][trim1_1], data['y'][trim1_1])
ax1.set_title(band1+' XY')

ax1 = fig.add_subplot(2, 2, 2)
ax1.scatter(data['x'][trim2_2], data['y'][trim2_2])
ax1.set_title(band2+' XY')

ax1 = fig.add_subplot(2, 2, 3)
ax1.scatter(data['x'][trim3_3], data['y'][trim3_3])
ax1.set_title(band3+' XY')

ax1 = fig.add_subplot(2, 2, 4)
ax1.scatter(data['x'][trim4_4], data['y'][trim4_4])
ax1.set_title(band4+' XY')
