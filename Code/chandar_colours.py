'''Chandar colour segmentation''' 

import os
import numpy as np 
from astropy.table import Table
from matplotlib import pyplot as plt
import pylab

data = Table.read("data_v4.txt", format="ascii.commented_header", guess=False)

#Clean data
ratio = 0.2

#Colour 1: U - B
U = data["mag05_336"]
U_unc = data['mag05_336_unc']
B = data['mag05_438']
B_unc = data['mag05_438_unc']

#Colour 2: V - I 
V = data['mag05_555']
V_unc = data['mag05_555_unc']
I = data['mag05_814']
I_unc = data['mag05_814_unc']

U_trim = np.logical_and(np.logical_and(U != -99, U_unc != -99),
                        U_unc < ratio)
B_trim = np.logical_and(np.logical_and(B != -99, B_unc != -99),
                        B_unc < ratio)
V_trim = np.logical_and(np.logical_and(V != -99, V_unc != -99),
                        V_unc < ratio)
I_trim = np.logical_and(np.logical_and(I != -99, I_unc != -99),
                        I_unc < ratio)

colour1_trim = np.logical_and(U_trim, B_trim)
colour2_trim = np.logical_and(V_trim, I_trim)

final_data = np.logical_and(colour1_trim, colour2_trim)

# Colours
U_B = U[final_data] - B[final_data]
V_I = V[final_data] - I[final_data]

cluster_space = np.logical_and(5*V_I - U_B > 2, 0.33*V_I - U_B > -0.33)
star_w_cluster = np.logical_and(5*V_I-U_B < 2, U_B < -0.8)
blue_star = np.logical_and(np.logical_and(U_B > - 0.8, U_B < 0.5),
                           0.33*V_I - U_B < -0.33)
yellow_star = np.logical_and(U_B > 0.5, 0.33 * V_I - U_B < -0.33)

chandar_lines = plt.figure(figsize=(8,8))
ax = chandar_lines.add_subplot(111)
ax.scatter(V_I, U_B, marker='.', s=2)
col = "r"   # set line colour

#Cluster Space
ax.plot([0.5, (2-0.33)/0.33], [0.5, 2], color=col)
ax.plot([0.5, -1], [0.5,-7], color=col)
ax.text(3, 0, "Cluster space")

# Star/Cluster space
ax.plot([-4, 0.24], [-0.8, -0.8], color=col)
ax.plot([0.24, -1], [-0.8, -7], color=col)
ax.text(-3, -2, "Star/Cluster space")

#Blue Star space
ax.plot([0.5, -4], [0.5, 0.5], color=col)
ax.text(-3, 0, "Blue Star space")

# Yellow star space
ax.text(-3, 1.25, "Yellow Star space")

#format plot
ax.set_xlabel("V - I")
ax.set_ylabel("U - B")

plt.ylim(-3, 2)
plt.xlim(-4, 5)

# plt.gca().invert_yaxis()

