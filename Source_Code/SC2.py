from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
'''
    The code below creates the image for Figure 8 in the
    User Guide.
    Figure 8 is a comparison F1 temperatures in saturated S7 pixels seen
    through different colormaps.
'''
##############################################################################
#STEP1: First the relevant file is called on (the path to it is defined in the
#string fname below) and arrays are created to store the BT values from the
#S7 and F1 bands
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
F1 = Dataset(fname+'F1_BT_in.nc').variables['F1_BT_in'][:]
S7 = Dataset(fname+'S7_BT_in.nc').variables['S7_BT_in'][:]

##############################################################################
#STEP2: Arrays are defined to store the boolean masks of the F1/S7 image
#arrays. 'True' values in these arrays indicate the value is masked, and in
#the corresponding image array this value will be labelled 'masked' and will
#not be seen when plotting the image for that particular position
##############################################################################
F1mask, S7mask = np.ma.getmaskarray(F1), np.ma.getmaskarray(S7)

##############################################################################
#STEP3: An array is created which stores only the values from F1 which are 
#unmasked in F1 but masked in S7
##############################################################################
satpixels = np.where(np.logical_and(S7mask==True, F1mask==False), F1, 0)

##############################################################################
#STEP4: The array created from STEP3 is displayed with two seperate colormaps
#using the matplotlib.plt methods 'imshow' and various methods using 
#'subplots'
##############################################################################
fig, (ax0, ax1) = plt.subplots(ncols = 2, sharey=True)
im = ax0.imshow(satpixels,cmap='YlOrRd', vmin=100, vmax=400)
axes1 = plt.axes([0.09, 0.17, 0.37, 0.02])
fig.colorbar(im,cax=axes1, orientation='horizontal')
ax0.set_xlabel('a)')
ax0.set_xticks([])
ax0.set_yticks([])

im1 = ax1.imshow(satpixels,cmap='YlOrRd', vmin=300, vmax=350)
axes = plt.axes([0.55, 0.17, 0.37, 0.02])
fig.colorbar(im1,cax=axes, orientation='horizontal')
ax1.set_xlabel('b)')
ax1.set_xticks([])
ax1.set_yticks([])
fig.text(0.28, 0.06, 'Kelvin', ha='center', va='center')
fig.text(0.74, 0.06, 'Kelvin', ha='center', va='center')

fig.subplots_adjust(wspace=0.9)
plt.tight_layout(pad = 4)
plt.savefig('Figures/Fig8.png', dpi=1000)
plt.show()