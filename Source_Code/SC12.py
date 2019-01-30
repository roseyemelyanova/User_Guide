from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib
##############################################################################
#STEP0: Define path to relevant .SEN3/ file
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'


'''
    The code between lines 18 and 53 creates the images for Figure 18 in the
    User Guide.
    Figure 18 contains two images of the interpolated arrays in x and y from
    the tie point data over the product grid.
    Steps labelled A refer to those taken in this process.
'''

##############################################################################
#STEPA1: Import relevant filees and check/correct for the difference between
#start offsets of the product image and cartesian tie point grid. Use this
#to create interpolated X/Y coordinate grids, x_new and y_new
##############################################################################
coords = Dataset(fname+ 'cartesian_in.nc') 
vals = Dataset(fname+'S8_BT_in.nc')
image_tie_point = Dataset(fname+'cartesian_tx.nc')
if vals.start_offset == image_tie_point.start_offset:
    start_offset = 0.0
else:
    start_offset = vals.start_offset * float(vals.resolution.split()[2]) / float(image_tie_point.resolution.split()[2]) - image_tie_point.start_offset
interpx = (np.array(range(vals.dimensions['columns'].size)) - vals.track_offset) * float(vals.resolution.split()[1]) / float(image_tie_point.resolution.split()[1]) + image_tie_point.track_offset
interpy = np.array(range(vals.dimensions['rows'].size)) * float(vals.resolution.split()[2]) / float(image_tie_point.resolution.split()[2]) + start_offset
tie_x = np.array(range(image_tie_point.dimensions['columns'].size))
tie_y = np.array(range(image_tie_point.dimensions['rows'].size))

f = interpolate.RectBivariateSpline(tie_y, tie_x, image_tie_point.variables['x_tx'][:])
x_new = f(interpy, interpx).astype(int)
f1 = interpolate.RectBivariateSpline(tie_y,tie_x, image_tie_point.variables['y_tx'][:])
y_new = f1(interpy,interpx).astype(int)
##############################################################################
#STEPA2: Plot the interpolated arrays
##############################################################################

plt.imshow(x_new)
plt.xticks([])
plt.yticks([])
plt.savefig('Figures/Fig18a.png', dpi=1000)
plt.show()

plt.imshow(y_new)
plt.xticks([])
plt.yticks([])
plt.savefig('Figures/Fig18b.png', dpi=1000)
plt.show()

'''
    The code between lines 69 and 98 creates the image for Figure 19 in the
    User Guide.
    Figure 19 is an image showing the relationship between j (the column 
    number) and the interpolated x value, and i (the row number) and the
    interpolated y value. This is used as an argument for Figure 20, because
    it involves assuming the relationship is linear.
    Steps labelled B refer to those taken in this process.
'''
##############################################################################
#STEPB1: Create the subplot and set limits on the axes
##############################################################################
figure, ax = plt.subplots(ncols=2, figsize=(12,4))
labels = {0:['interpolated x value','j'], 1:['interpolated y value','i']}
ax[0].set_xlim(x_new.min(), x_new.max())
ax[0].set_ylim(0, 1499)
ax[0].set(xticks=([x_new.min(),0,x_new.max()]), yticks=([0,vals.track_offset,1499]))
ax[0].xaxis.grid(True)

##############################################################################
#STEPB2: Find the y intercept and gradient of the slope for y_new, then plot
#in the subplots created in STEPB1
##############################################################################
c2 = 1199.0/float(y_new.max()-y_new.min())
i_int = -c2*y_new.min()
ax[1].set_xlim(0,30000000)
ax[1].set_ylim(i_int,c2*30000000 + i_int)
ax[1].set(xticks=([0]), yticks=([i_int,0,c2*30000000 + i_int]))
ax[1].yaxis.grid(True)
for i in range(2):
    ax[i].set_xlabel(labels[i][0])
    ax[i].set_ylabel(labels[i][1], rotation=0)

ax[0].plot(x_new[0],np.arange(1500), 'b-')
ax[1].plot(y_new[:,0],np.arange(1200), 'b-')
ax[1].plot([0,y_new.min()],[i_int,0],'b--')
ax[1].plot([y_new.max(),30000000],[1199,c2*30000000 + i_int],'b--')

plt.subplots_adjust(wspace=0.8)
plt.tight_layout()
plt.savefig('Figures/Fig19.png', dpi=1000)
plt.show()

'''
    The code between lines 111 and 244 creates the image for Figure 20 in the
    User Guide.
    Figure 20 is the correctly regridded image, which is identical to the 
    original image. This can be checked by calling:
        S8 = Dataset(fname+'S8_BT_in.nc').variables['S8_BT_in']
        S8 - template
    at the end of this section or in the console after this script has been
    run.
    Steps labelled C refer to those taken in this process.
'''
##############################################################################
#STEPC1: create arrays containing information about the BT value (in this case, 
#because S8 is used) and X/Y coordinates for both orphan and non-orphan pixels, 
#as well as the detector value
##############################################################################
S8, S8_orphan = vals.variables['S8_BT_in'][:], vals.variables['S8_BT_orphan_in'][:]
x = coords.variables['x_in'][:]
y = coords.variables['y_in'][:]
x_orphan = coords.variables['x_orphan_in'][:]
y_orphan = coords.variables['y_orphan_in'][:]
det = Dataset(fname+'indices_in.nc').variables['detector_in'][:]

##############################################################################
#STEPC2: define i,j functions to calculate row and column number
##############################################################################
def j(xval):
    c1 = 1499.0/float(x_new.max()-x_new.min())
    return vals.track_offset - c1*xval

def i(yval):
    c2 = 1199.0/float(y_new.max()-y_new.min())
    y0 = c2*y_new.min()
    return c2*yval - y0

##############################################################################
#STEPC3: create cosmetic pixel mask, then select valid pixels only. Note that
#orphan pixel information is in 1D arrays as opposed to 2D, and restricts 
#the X/Y values to lie within the interpolated y_new. checks aren't in place 
#for x_new because that already has plenty of space to the sides due to the 
#presence of masked pixels, however an extension in Y will prove impossible to 
#fit into a 1200 by 1500 gridded image
##############################################################################
cospix = (Dataset(fname+'flags_in.nc').variables['confidence_in'][:] & 1 << 8) >> 8

x_in = np.where(np.logical_and(cospix==0,np.ma.getmaskarray(S8)==False), x, 0)
y_in = np.where(np.logical_and(cospix==0,np.ma.getmaskarray(S8)==False), y, 0)
S8_in =  np.where(np.logical_and(cospix==0,np.ma.getmaskarray(S8)==False), S8, 0)

x_cos = np.where(np.logical_and(cospix==1,np.ma.getmaskarray(S8)==False), x, 0)
y_cos= np.where(np.logical_and(cospix==1,np.ma.getmaskarray(S8)==False), y, 0)
S8_cos = np.where(np.logical_and(cospix==1,np.ma.getmaskarray(S8)==False), S8, 0)

x_in, x_cos = np.ma.masked_values(x_in,0), np.ma.masked_values(x_cos,0)
y_in, y_cos = np.ma.masked_values(y_in,0), np.ma.masked_values(y_cos,0)
S8_in, S8_cos = np.ma.masked_values(S8_in,0), np.ma.masked_values(S8_cos,0)

idx_orphan = np.where(y_orphan<y_new.max()) 
x_orphan_vals = x_orphan[idx_orphan]
y_orphan_vals = y_orphan[idx_orphan]
S8_orphan_vals = S8_orphan[idx_orphan]

##############################################################################
#STEPC4: find the i,j values of pixels, and round the i,j values for the 
#orphan pixels as these will need to be called as indices to the interpolated
#grids
##############################################################################
i_vals, j_vals = i(y_in), j(x_in)
i_cvals, j_cvals = i(y_cos), j(x_cos)
i_orphan_vals, j_orphan_vals = i(y_orphan_vals), j(x_orphan_vals)

idx_orphan_i = np.round(i_orphan_vals).astype(int)
idx_orphan_j = np.round(j_orphan_vals).astype(int)

##############################################################################
#STEPC5: create interpolated grids with 1 extra row, because interpolated grid
#positions represent the centers of pixels and some pixels may acceptably fall
#0.5km outside of this
##############################################################################
x_new_temp = np.vstack((x_new,x_new[1199,:]))
x_rolled = np.roll(x_new_temp,1)

y_new_temp = np.vstack((y_new,np.repeat((y_new[-1,0]+ float(y_new.max()-y_new.min())/1199.0).astype(int),1500)))
y_rolled = np.transpose(np.roll(np.transpose(y_new_temp),1))

##############################################################################
#STEPC6: find and apply corrections, note the difference in methods used 
#between orphan and non orphan pixels
##############################################################################
Dx, D_x_cos = x_new - x_in, x_new - x_cos
dx = x_new - x_rolled[0:1200,:]
dx[:,0] = 1000
Dj, Dj_cos = Dx.astype(float)/dx.astype(float), D_x_cos.astype(float)/dx.astype(float)

Dy, D_y_cos = y_new - y_in, y_new - y_cos
dy = y_new-y_rolled[0:1200,:]
dy[0,:] = 1004
Di, Di_cos = Dy.astype(float)/dy.astype(float), D_y_cos.astype(float)/dy.astype(float)

D_x_orphan = x_new_temp[idx_orphan_i,idx_orphan_j] - x_orphan_vals
dx = x_new_temp[idx_orphan_i, idx_orphan_j] - x_rolled[idx_orphan_i, idx_orphan_j]
Dj_orphan = D_x_orphan.astype(float)/dx.astype(float)

D_y_orphan = y_new_temp[idx_orphan_i,idx_orphan_j] - y_orphan_vals
dy = y_new_temp[idx_orphan_i,idx_orphan_j] - y_rolled[idx_orphan_i,idx_orphan_j]
Di_orphan = D_y_orphan.astype(float)/dy.astype(float)

##############################################################################
#STEPC7: apply the correction to calculated i,j values, then round them
##############################################################################
i_vals1, j_vals1 = i_vals + Di, j_vals + Dj
j_cvals1, i_cvals1 = j_cvals + Dj_cos, i_cvals + Di_cos
i_orphan_vals1, j_orphan_vals1 = i_orphan_vals + Di_orphan, j_orphan_vals + Dj_orphan

i_valsrounded = np.round(i_vals1).astype(int)
j_valsrounded = np.round(j_vals1).astype(int)

i_cvalsrounded = np.round(i_cvals1).astype(int)
j_cvalsrounded = np.round(j_cvals1).astype(int)

i_orphan_vals_rounded = np.round(i_orphan_vals1).astype(int)
j_orphan_vals_rounded = np.round(j_orphan_vals1).astype(int)

##############################################################################
#STEPC8: apply the correction to calculated i,j values, then round them and 
#create/fill the product grid, 'template' with their corresponding S8 values
##############################################################################

idx, idx_cos = np.where(np.logical_and(cospix==0,np.ma.getmaskarray(S8)==False)), np.where(np.logical_and(cospix==1,np.ma.getmaskarray(S8)==False))
all_i_vals = np.concatenate((i_orphan_vals_rounded,i_valsrounded[idx], i_cvalsrounded[idx_cos]))
all_j_vals = np.concatenate((j_orphan_vals_rounded,j_valsrounded[idx], j_cvalsrounded[idx_cos]))
all_S8 = np.ma.concatenate((S8_orphan_vals,S8_in[idx], S8_cos[idx_cos]))

template = np.zeros([1200,1500])
template[all_i_vals,all_j_vals] = all_S8
template = np.ma.masked_values(template,0)

im = plt.imshow(template, cmap='RdBu')
plt.xticks([])
plt.yticks([])
cb = plt.colorbar(im)

cb.set_label('Kelvin')
plt.tight_layout()
plt.savefig('Figures/Fig20.png',dpi=1000)

'''
    The code between lines 253 and 390 creates the image for Figure 21 in the
    User Guide.
    Figure 21 is a comparison between regridded images where the orphan pixels
    have been selected in different ways.
    Steps labelled D refer to those taken in this process.
'''
##############################################################################
#STEPD1: Create an orphan swapped image by simply changing how the arrays 
#called to make the tmeplate in lines 217-219 are organised
##############################################################################
all_i_vals = np.concatenate((i_valsrounded[idx], i_cvalsrounded[idx_cos],i_orphan_vals_rounded))
all_j_vals = np.concatenate((j_valsrounded[idx], j_cvalsrounded[idx_cos],j_orphan_vals_rounded))
all_S8 = np.ma.concatenate((S8_in[idx], S8_cos[idx_cos],S8_orphan_vals))

orph_swap = np.zeros([1200,1500])
orph_swap[all_i_vals,all_j_vals] = all_S8
orph_swap = np.ma.masked_values(orph_swap,0)
orph_swap[0,998] = template[0,998] #NOTE: by inspection this cell takes the fill value. 

##############################################################################
#STEPD2: Create the nearest neighbour image
#STEPD2a: Create a 1D array containing a combination of i and j such that 
#non-unique values indicate orphan pixels
##############################################################################
temp_values = all_i_vals + 10000000000*all_j_vals

##############################################################################
#STEPD2b: Find information about unique values of temp_array and use it to find 
#duplicate indices in temp_values, corresponding to orphans in original 1D 
#arrays all_i_vals, all_j_vals and all_S8
##############################################################################
temp_uniq, temp_uniq_idx ,temp_uniq_idv, temp_uniq_cnt = np.unique(temp_values, return_index = True, return_inverse=True, return_counts=True)
cnt_idx, = np.nonzero(temp_uniq_cnt > 1) 
idx_mask = np.in1d(temp_uniq_idv, cnt_idx)
idx_idx, = np.nonzero(idx_mask)
srt_idx = np.argsort(temp_uniq_idv[idx_mask])
dup_idx = np.split(idx_idx[srt_idx], np.cumsum(temp_uniq_cnt[temp_uniq_cnt > 1])[:-1]) 

##############################################################################
#STEPD2c: remove invalid values in dup_idx. This is a byproduct of the process
#in STEPD2b and can be improved upon
##############################################################################
i = 0
while i<len(dup_idx):
    if dup_idx[i].shape != (2L,):
        del dup_idx[i]
    i += 1

##############################################################################
#STEPD2d: Calculate the root mean squared (rms) distance to center of pixel in 
#X/Y and find the difference between these values for duplicate values
##############################################################################
#STEP18d: 
all_Dx = np.concatenate((D_x_orphan, Dx[idx],D_x_cos[idx_cos]))
all_Dy = np.concatenate((D_y_orphan, Dy[idx],D_y_cos[idx_cos]))
dist = np.sqrt(all_Dx**2 + all_Dy**2)

o_idx = np.concatenate(([dup_idx[h] for h in range(len(dup_idx))])) 
o_idx1, o_idx2 = o_idx[::2], o_idx[1::2]
compare = dist[o_idx1] - dist[o_idx2]

##############################################################################
#STEPD2e: create 1D arrays which minimise rms distances, fill a new array
#with these values as well as values of all unique i,j combinations
##############################################################################
c1 = np.ma.masked_values(np.where(compare>0,o_idx2,-1),-1)
c2 = np.ma.masked_values(np.where(compare<0,o_idx1,-1),-1)

new_i_vals = np.concatenate((all_i_vals[temp_uniq_idx],all_i_vals[c1], all_i_vals[c2]))
new_j_vals = np.concatenate((all_j_vals[temp_uniq_idx],all_j_vals[c1], all_j_vals[c2]))
new_S8_vals = np.ma.concatenate((all_S8[temp_uniq_idx],all_S8[c1], all_S8[c2]))

nearest = np.zeros([1200,1500])
nearest[new_i_vals, new_j_vals] = new_S8_vals
nearest = np.ma.masked_values(nearest,0)
nearest[0,998] = template[0,998]#NOTE: by inspection this cell takes the fill value. 

##############################################################################
#STEPD3: Create Weighted Pixel Image
#STEPD3a: Create 1D arrays from the corrections for non-orphans and cosmetic 
#pixels, containing values from idx and idx_cos
##############################################################################

D_i, D_j, D_i_cos, D_j_cos = Di[idx], Dj[idx], Di_cos[idx_cos], Dj_cos[idx_cos]

##############################################################################
#STEPD3b:calculate the rms distance to the center of the pixel in i/j and
#weight duplicate pixels based on the inverse of this value and use this
#to find a weighted S8 result for those pixels
##############################################################################
all_Di = np.concatenate((D_i,D_i_cos, Di_orphan))
all_Dj = np.concatenate((D_j,D_j_cos, Dj_orphan))
dist1 = np.sqrt(all_Di**2 + all_Dj**2)

m = np.where(dist1==0,1,1/dist1)

average_vals = ((m[o_idx1]*all_S8[o_idx1]) + (m[o_idx2]*all_S8[o_idx2]))/(m[o_idx1]+m[o_idx2])
uniq, indexs = np.unique(temp_values[o_idx], return_index = True) 

weighted_i_vals = np.concatenate((all_i_vals[temp_uniq_idx],all_i_vals[o_idx][indexs]))
weighted_j_vals = np.concatenate((all_j_vals[temp_uniq_idx],all_j_vals[o_idx][indexs]))
weighted_S8_vals = np.ma.concatenate((all_S8[temp_uniq_idx],average_vals))

weighted = np.zeros([1200,1500])
weighted[weighted_i_vals,weighted_j_vals] = weighted_S8_vals
weighted = np.ma.masked_values(weighted,0)

##############################################################################
#STEPD4: Plot the 4 different arrays in a given interval. The interval was
#detected by finding where the largest differences between the arrays were, 
#this process hasn't been included here
##############################################################################
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}

matplotlib.rc('font', **font)

fig, axes = plt.subplots(ncols=2, nrows=2, sharey= True, sharex = True, figsize = (5,5))
lbl = np.array([['a) SLSTR','b) Orphan'],['c) Neighbour','d) Weighted']])
for i,j in np.ndenumerate(axes):
    axes[i].set_xlim([1403,1405])
    axes[i].set_xticks([1403,1404,1405])
    axes[i].set(xlabel=lbl[i])
    axes[i].set_ylim([365,367])
    axes[i].set_yticks([365,366,367])


min_vals = [template[365:367,1403:1405].min(), orph_swap[365:367,1403:1405].min(), nearest[365:367,1403:1405].min(), weighted[365:367,1403:1405].min()]
max_vals = [template[365:367,1403:1405].max(), orph_swap[365:367,1403:1405].max(), nearest[365:367,1403:1405].max(), weighted[365:367,1403:1405].max()]
    
im0 = axes[0,0].imshow(template, vmin=min(min_vals), vmax= max(max_vals))
axes[0,1].imshow(orph_swap, vmin=min(min_vals), vmax= max(max_vals))
axes[1,0].imshow(nearest, vmin=min(min_vals), vmax= max(max_vals))
axes[1,1].imshow(weighted, vmin=min(min_vals), vmax= max(max_vals))
plt.subplots_adjust(left=0.15,right=0.77,bottom=0.2,top=0.9,wspace = 0.4)
cax = plt.axes([0.837,0.2,0.02,0.7])
cb = plt.colorbar(im0, cax)
cb.set_label('Kelvin')
fig.text(0.47, 0.08, 'x', ha='center', va='center')
fig.text(0.03, 0.53, 'y', ha='center', va='center')
plt.savefig('Figures/Fig21.png', dpi=1000)

plt.show()