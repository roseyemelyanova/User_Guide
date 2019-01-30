from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figure 29 in the
    User Guide.
    Figure 29 is a comparison between solar and satellite azimuth and zenith
    angles across the product grid dimensions/size.
'''
##############################################################################
#STEP1: Define the path to the relevant .SEN3/ file and assign necessary
#netcdf rootgroups to variables
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180822T200909_20180822T201209_20180822T221016_0179_035_028_2700_MAR_O_NR_002.SEN3/'
geomet = Dataset(fname+'geometry_tn.nc')
vals = Dataset(fname+'S8_BT_in.nc')

##############################################################################
#STEP2: Interpolate the geometric data over the range of the product grid.
#Note that line __ calls the integer '111' - this is because upon inspection,
#although the netCDF file used here has 130 columns, all columns from 111 
#onwards in all variables are masked. This has a big impact on the 
#interpolation method, and produces a null result
##############################################################################
if vals.start_offset == geomet.start_offset:
    start_offset = 0.0
else:
    start_offset = vals.start_offset * float(vals.resolution.split()[2]) / float(geomet.resolution.split()[2]) - geomet.start_offset


i_x = (np.array(range(vals.dimensions['columns'].size)) - vals.track_offset) * float(vals.resolution.split()[1]) / float(geomet.resolution.split()[1]) + geomet.track_offset
i_y = np.array(range(vals.dimensions['rows'].size)) * float(vals.resolution.split()[2]) / float(geomet.resolution.split()[2]) + start_offset
t_x = np.array(range(111))
t_y = np.array(range(geomet.dimensions['rows'].size))

##############################################################################
#STEP3: Define a function to interpolate over a specific variable in the 
#geometry rootgroup
##############################################################################
def grid(var):
    f = interpolate.RectBivariateSpline(t_y,t_x, geomet.variables[var][:, 0:111])
    arr = f(i_y,i_x)
    return arr

##############################################################################
#STEP4: Create the plot by using the function defined in STEP3 and using 
#subplotting
##############################################################################
fig, axes = plt.subplots(nrows = 2, ncols = 2)

solar_azimuth = grid('solar_azimuth_tn')
sat_azimuth = grid('sat_azimuth_tn')
solar_zenith = grid('solar_zenith_tn')
sat_zenith = grid('sat_zenith_tn')

idx = [[0,0],[0,1],[1,0],[1,1]]

axes[0,0].set_title('Solar Azimuth')
axes[0,0].set_xticks([])
axes[0,0].set_yticks([])
im0 = axes[0,0].imshow(solar_azimuth)

axes[0,1].set_title('Satellite Azimuth')
axes[0,1].set_xticks([])
axes[0,1].set_yticks([])
im1 = axes[0,1].imshow(sat_azimuth)

axes[1,0].set_title('Solar Zenith')
axes[1,0].set_xticks([])
axes[1,0].set_yticks([])
im2 = axes[1,0].imshow(solar_zenith)

axes[1,1].set_title('Satellite Zenith')
axes[1,1].set_xticks([])
axes[1,1].set_yticks([])
im3 = axes[1,1].imshow(sat_zenith)

fig.subplots_adjust(hspace = 0.8)

cax0= plt.axes([0.448,0.53,0.02,0.38])
cax2= plt.axes([0.448,0.05,0.02,0.375])
cax1= plt.axes([0.88,0.53,0.02,0.38])
cax3= plt.axes([0.88,0.05,0.02,0.375])

fig.colorbar(im0, cax = cax0)
fig.colorbar(im1, cax = cax1)
fig.colorbar(im2, cax = cax2)
fig.colorbar(im3, cax = cax3)
fig.text(0.98, 0.5, 'degrees', ha='center', va='center', rotation=270)
plt.tight_layout()
plt.savefig('Figures/Fig29.png', dpi=1000)