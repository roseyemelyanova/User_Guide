from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
'''
    The code below creates the images for Figure 16 in the
    User Guide.
    Figure 16 is a comparison between rounded and unrounded X/Y values.
'''
##############################################################################
#STEP1: Define the path to the relevant .SEN3/ file and create arrays to store
#relevant information
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
det = Dataset(fname+'indices_in.nc').variables['detector_in'][:]
deto = Dataset(fname+'indices_in.nc').variables['detector_orphan_in'][:]

cart = Dataset(fname+'cartesian_in.nc')
x, xo = cart.variables['x_in'][:], cart.variables['x_orphan_in'][:]
y, yo = cart.variables['y_in'][:], cart.variables['y_orphan_in'][:]
##############################################################################
#STEP2: Create the cosmetic pixel mask, cospix
##############################################################################
mask = 1 << 8
flag = Dataset(fname+'flags_in.nc').variables['confidence_in'][:]
cospix = (flag & mask) >> 8

##############################################################################
#STEP3: Only select valid indices, ie. where pixels are not cosmetically
#masked and are seen by one of the two available detectors. Orphan pixels
#do not have cosmetic masking
##############################################################################
cosidx = np.ma.where(np.logical_and(cospix==0,np.logical_or(det==0, det==1)))
idxo = np.ma.where(np.logical_or(deto==0,deto==1))
idx = np.ma.where(np.logical_or(det==0,det==1))

##############################################################################
#STEP4: divide by 1000 to get values in km, and add all x and y coordinates 
#together (i.e. orphan and non-orphan)
##############################################################################
x_unrounded = x[cosidx].astype('f4')/1000
x_orph_unrounded = xo[idxo].astype('f4')/1000
y_unrounded = y[cosidx].astype('f4')/1000
y_orph_unrounded = yo[idxo].astype('f4')/1000

allx_unrounded = np.append(x_unrounded,x_orph_unrounded)
ally_unrounded = np.append(y_unrounded,y_orph_unrounded)
##############################################################################
#STEP5: divide by 1000 to get values in km, round the values to the nearest km 
#and add all x and y rounded coordinates together (i.e. orphan and non-orphan)
##############################################################################

x_rounded = np.ma.round(x[cosidx].astype('f4')/1000).astype(int)
x_orph_rounded = np.ma.round(xo[idxo].astype('f4')/1000).astype(int)
y_rounded = np.ma.round(y[cosidx].astype('f4')/1000).astype(int)
y_orph_rounded = np.ma.round(yo[idxo].astype('f4')/1000).astype(int)

allx_rounded = np.append(x_rounded,x_orph_rounded)
ally_rounded = np.append(y_rounded,y_orph_rounded)

##############################################################################
#STEP6: Define the function to constrain the axes, as both resulting images
#will have the same axes
##############################################################################
def makeplt(axis):
    axis.set_xticks(np.arange(allx_rounded.min(), allx_rounded.max(), 1))
    axis.set_yticks(np.arange(ally_rounded.min(), ally_rounded.max(), 1))
    axis.xaxis.grid(True)
    axis.yaxis.grid(True)
    axis.set_xlim(allx_rounded.min(),allx_rounded.min()+10)
    axis.set_ylim(ally_rounded.min(),ally_rounded.min()+10)
    axis.set_aspect('auto')
    axis.set_xlabel('X')
    axis.set_ylabel('Y')
##############################################################################
#STEP7: Create the figures required
##############################################################################
fig = plt.figure()
ax = fig.gca()
makeplt(ax)
ax.scatter(allx_unrounded,ally_unrounded,marker='+', color='red')
plt.savefig('Figures/Fig16a.png', dpi=1000)
plt.show()


fig1 = plt.figure()
ax1 = fig1.gca()
makeplt(ax1)
ax1.scatter(allx_unrounded,ally_unrounded,marker='+', color='red')
ax1.scatter(allx_rounded,ally_rounded,marker='+', color='blue')
plt.legend(['rounded', 'unrounded'],loc='upper right', bbox_to_anchor=(1, 1.15), ncol=2)
plt.savefig('Figures/Fig16b.png', dpi=1000)
plt.show()

