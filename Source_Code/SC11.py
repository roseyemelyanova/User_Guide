from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import SC9
'''
    The code below creates the image for Figure 17 in the
    User Guide.
    Figure 17 is a comparison between a regridding attempt of the S8 band and
    the actual S8 image.
'''
##############################################################################
#STEP1: Define the variable to the relevant .SEN3/ file and import variables
#from SC9
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
det, deto = SC9.det, SC9.deto
S8, S8o = SC9.S8, SC9.S8o

##############################################################################
#STEP2: Create a cosmetic pixel mask
##############################################################################
mask = 1 << 8
flag = Dataset(fname+'flags_in.nc').variables['confidence_in'][:]
cospix = (flag & mask) >> 8

##############################################################################
#STEP3: fill arrays with X/Y information about all pixels - orphan and non 
#orphan
##############################################################################
cart = Dataset(fname+'cartesian_in.nc')
x, xo = cart.variables['x_in'][:], cart.variables['x_orphan_in'][:]
y, yo = cart.variables['y_in'][:], cart.variables['y_orphan_in'][:]

##############################################################################
#STEP4: find unmasked pixels (ie. where detectors are either 0 or 1 in this 
#case)
##############################################################################
idx = np.ma.where(np.logical_or(det==0, det==1))
idxo = np.ma.where(np.logical_or(deto ==0,deto ==1))

##############################################################################
#STEP5: Round the x and y values which are valid to the nearest km.
##############################################################################
x_rounded = np.ma.round(x[idx].astype('f4')/1000).astype(int)
x_orph_rounded = np.ma.round(xo[idxo].astype('f4')/1000).astype(int)
y_rounded = np.ma.round(y[idx].astype('f4')/1000).astype(int)
y_orph_rounded = np.ma.round(yo[idxo].astype('f4')/1000).astype(int)

##############################################################################
#STEP6: add all X/Y values together
##############################################################################
allx_rounded = np.append(x_rounded,x_orph_rounded)
ally_rounded = np.append(y_rounded,y_orph_rounded)

##############################################################################
#STEP7: Find unique values in a combination of X/Y, such that non-unique values 
#are orphan pixels
##############################################################################

temp_coords = allx_rounded + 10000*ally_rounded 
ucoords, indices,inv_idx,ucount= np.unique(temp_coords,return_index = True, return_inverse = True, return_counts=True)
S8all = np.append(S8[idx],S8o[idxo]).astype(int)

##############################################################################
#STEP8: Create a product grid and fill it with relevant pixel information
##############################################################################
new_grid = np.ma.zeros([ally_rounded.max()-ally_rounded.min()+1,allx_rounded.max()-allx_rounded.min()+1])
new_grid[ally_rounded[indices]-ally_rounded.min(),allx_rounded[indices]-allx_rounded.min()] = S8all[indices]
new_grid = np.ma.masked_values(new_grid,0)

fig, axes = plt.subplots(ncols = 2)
for i in range(2):
    axes[i].set_xticks([])
    axes[i].set_yticks([])

im = axes[0].imshow(new_grid[:,::-1],cmap = 'RdBu')
im1 = axes[1].imshow(S8, cmap='RdBu', aspect=1)

cax = plt.axes([0.15,0.16,0.73,0.02])
cb = plt.colorbar(im, cax, orientation='horizontal')
cb.set_label('Kelvin')
plt.subplots_adjust(wspace = 0.2)
plt.savefig('Figures/Fig17.png',dpi=1000)