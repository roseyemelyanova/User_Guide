from netCDF4 import Dataset
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figure 14 in the
    User Guide.
    Figure 14 shows the distribution of orphan pixels by scan and pixel number.
'''
##############################################################################
#STEP1: Define the path to the relevant .SEN3/ file and store relevant 
#information from it in variables
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
coords = Dataset(fname+'indices_in.nc')
x = coords.variables['pixel_orphan_in'][:]
y = coords.variables['scan_orphan_in'][:]
orphval = Dataset(fname+'S8_BT_in.nc').variables['S8_BT_orphan_in'][:]

##############################################################################
#STEP4: Plot the information
##############################################################################
plt.xlim(x.min(),x.max())
plt.ylim(y.min(),y.max())
plt.xlabel('pixel number')
plt.ylabel('scan number')
plt.scatter(x,y,s=0.1,c=orphval, marker=',', cmap='RdBu')
plt.savefig('Figures/Fig14.png', dpi=1000)
plt.show()
