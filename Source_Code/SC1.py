from netCDF4 import Dataset
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figure 7 in the
    User Guide.
    Figure 7 is a comparison of images between the S7 and F1 band.
'''
##############################################################################
#STEP1: First the relevant file is called on (the path to it is defined in the
#string fname below) and arrays are created to store the BT values from the
#S7 and F1 bands
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
S7 = Dataset(fname+'S7_BT_in.nc').variables['S7_BT_in'][:]
F1 = Dataset(fname+'F1_BT_in.nc').variables['F1_BT_in'][:]

##############################################################################
#STEP2: The minimum and maximum values of these image arrays are found and
#stored in lists
##############################################################################
mins = [S7.min(), F1.min()]
maxs = [S7.max(), F1.max()]
##############################################################################
#STEP3: These arrays are displayed using the matplotlib methods 'imshow' and
#various methods relating to subplotting. The plots are labelled and a 
#colorbar axis is adjusted to fit around them
##############################################################################
plt.subplot(121)
plt.imshow(S7, cmap='RdBu', vmin=min(mins), vmax=max(maxs))
plt.xlabel('a)')
plt.xticks([])
plt.yticks([])
plt.subplot(122)
plt.xlabel('b)')
plt.imshow(F1, cmap='RdBu', vmin=min(mins), vmax=max(maxs))
plt.xticks([])
plt.yticks([])

plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, wspace=0.3)
cax = plt.axes([0.85, 0.32, 0.02, 0.36])
cb = plt.colorbar(cax=cax)
cb.set_label('Kelvin')
plt.savefig('Figures/Fig7.png', dpi=1000)
plt.show()
