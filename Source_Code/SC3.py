from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
'''
    The code below creates the image for Figure 9 in the
    User Guide.
    Figure 9 is an RGB image.
'''
##############################################################################
#STEP1: A function is defined which will output an rgb image given a path
#to a .SEN3/ file
##############################################################################
def make_rgb(path):
    ##########################################################################
    #STEP1a: Arrays are created to store radiance information from 4 different
    #bands
    ##########################################################################
    s1 = Dataset(path+'S1_radiance_an.nc').variables['S1_radiance_an'][:]
    s2 = Dataset(path+'S2_radiance_an.nc').variables['S2_radiance_an'][:]
    s3 = Dataset(path+'S3_radiance_an.nc').variables['S3_radiance_an'][:]
    s5 = Dataset(path+'S5_radiance_an.nc').variables['S5_radiance_an'][:]
    ##########################################################################
    #STEP1b: The difference between s5-s1 arrays is stored in a new array.
    #new arrays are created to store the red, green and blue colors.
    ##########################################################################
    ndsi = s5 - s1
    r    = np.where(ndsi > 0, s5, s3)
    g    = np.where(ndsi > 0, s3, s2)
    b    = np.where(ndsi > 0, s2, s1)
    ##########################################################################
    #STEP1c: The red, green and blue arrays are stacked together into one
    #array, normalised and returned as a result of calling the function
    ##########################################################################
    rgb  = np.dstack((r / r.max(), g / g.max(), b / b.max()))
    return rgb

##############################################################################
#STEP2: define the path to the .SEN3/ file being used and call the function to 
#store its output in a new variable.
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
rgb1 = make_rgb(fname)

##############################################################################
#STEP3: create the plot using the matplotlib.plt method 'imshow'
##############################################################################
plt.yticks([])
plt.xticks([])
plt.imshow(rgb1)
plt.savefig('Figures/Fig9.png',dpi=1000)
plt.show()