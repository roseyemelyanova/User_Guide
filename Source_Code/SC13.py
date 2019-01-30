from netCDF4 import Dataset
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figure 22 in the
    User Guide.
    Figure 22 is a plot of calibration offsets per pixel number in S6.
'''
##############################################################################
#STEP1: Define the path to the relevant .SEN3/ file and store useful 
#information in arrays
##############################################################################

fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
S6 = Dataset(fname+'S6_quality_an.nc')
caloffset = S6.variables['S6_cal_offset_an'][:]

##############################################################################
#STEP2: Plot the array
##############################################################################
plt.plot(caloffset[0,1])
plt.xlabel('pixel number')
plt.ylabel('calibration offset')
plt.tight_layout()
plt.savefig('Figures/Fig22.png', dpi=1000)