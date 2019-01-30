from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
'''
    The code below creates the image for Figure 11 in the
    User Guide.
    Figure 11 is an image of F1 pixels which are saturated in S7 and unmasked
    by certain masks, illustrated from the results of SC4.
'''
##############################################################################
#STEP1: define a function which will remove selected masked pixels from the 
#image and plot the result
##############################################################################
def pltwom(fname,types,bitnos,band,mid,suffix):
    ##########################################################################
    #STEP1a: import relevant arrays and read in bitmask data, selecting only
    #unmasked values from the selected band
    ##########################################################################
    val = Dataset(fname+band+mid+suffix+'.nc').variables[band+mid+suffix][:]
    S7 = Dataset(fname+'S7_BT_in.nc').variables['S7_BT_in'][:]
    for i in range(len(types)):
        maskval = 1 << bitnos[i]
        flag = Dataset(fname+'flags'+suffix+'.nc').variables[types[i]+suffix][:]
        mask = (flag & maskval) >> bitnos[i]
        val = np.ma.where(mask==0,val,0)
    ##########################################################################
    #STEP1b: create boolean arrays containing information about masked pixels
    #in the selected band and S7. Then fill a new array with values from the
    #selected band which are unmasked within the band but masked in S7
    ##########################################################################
    valmask, S7mask = np.ma.getmaskarray(val), np.ma.getmaskarray(S7)
    satpix = np.where(np.logical_and(S7mask==True, valmask==False), val, 0)
    val = np.ma.where(satpix != 0, val,0)
    val = np.ma.masked_values(val,0)
    im = plt.imshow(val, cmap='RdBu', vmin=100,vmax=400)
    plt.xticks([])
    plt.yticks([])
    cbar = plt.colorbar(im)
    cbar.set_label('Kelvin')
    plt.savefig('Figures/Fig11.png', dpi=1000)
    plt.show()

##############################################################################
#STEP2: Define the path to the .SEN3/ file and call the function
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
pltwom(fname,['confidence','confidence','cloud'],[12,8,7],'F1','_BT','_in')