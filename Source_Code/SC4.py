from netCDF4 import Dataset
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figure 10 in the
    User Guide.
    Figure 10 is a comparison between different masks.
'''
##############################################################################
#STEP1: define a function which plots any mask in any flag file
##############################################################################
def pltmask(fname,axis, typeof,bitno,suffix):
    ##########################################################################
    #STEP1a: use boolean operators to read the bitmasked array, and use the
    #matplotlib.pyplot method 'imshow' to display the result
    ##########################################################################
    maskval = 1 << bitno
    flag = Dataset(fname+'flags'+suffix+'.nc').variables[typeof+suffix][:]
    mask = (flag & maskval) >> bitno
    axis.imshow(mask, cmap= 'YlOrRd',vmin=0.5,vmax=1,alpha=1)

##############################################################################
#STEP3: Define the path to the .SEN3/ file and create subplots for the images
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
#STEP2: create the plot
fig, (ax0,ax1,ax2) = plt.subplots(nrows=1, ncols=3, sharey= True)

im0= pltmask(fname,ax0,'cloud',7,'_in') # produces mask showing cloud pixels
im1= pltmask(fname,ax1,'confidence',12,'_in') # produces mask showing sun glint pixels
im2= pltmask(fname,ax2,'confidence',8,'_in') # produces mask showing cosmetically filled pixels

ax0.set_xlabel('a) Cloud')
ax0.set_xticks([])
ax0.set_yticks([])

ax1.set_xlabel('b) Sun Glint')
ax1.set_xticks([])
ax1.set_yticks([])

ax2.set_xlabel('c) Cosmetic')
ax2.set_xticks([])
ax2.set_yticks([])

plt.savefig('Figures/Fig10.png', dpi=1000)
plt.show()

