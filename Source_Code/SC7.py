from netCDF4 import Dataset
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
'''
    The code below creates the image for Figure 13 in the
    User Guide.
    Figure 13 is a combination of images on one portion of a 2d projection of
    the earth.
'''
##############################################################################
#STEP1: Define paths to relevant .SEN3/ files
##############################################################################
file_1 = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
file_2 = 'Data/S3A_SL_1_RBT____20180819T194643_20180819T194943_20180819T214913_0179_034_370_2700_MAR_O_NR_002.SEN3/'
file_3 = 'Data/S3A_SL_1_RBT____20180822T200909_20180822T201209_20180822T221016_0179_035_028_2700_MAR_O_NR_002.SEN3/'

##############################################################################
#STEP2: Create variables to store the rootgroup information for all relevant
#files listed in STEP1
##############################################################################
fnames = [file_1,file_2,file_3]
band, geo = 'S7_BT_in','geodetic_in'
F1_1, geo_1 = Dataset(fnames[0]+band+'.nc'),Dataset(fnames[0]+geo+'.nc')
F1_2, geo_2 = Dataset(fnames[1]+band+'.nc'),Dataset(fnames[1]+geo+'.nc')
F1_3, geo_3 = Dataset(fnames[2]+band+'.nc'),Dataset(fnames[2]+geo+'.nc')

##############################################################################
#STEP3: Use cartopy to create a projection of the portion of the earth's
#surface spanning the extent of the longitudes and latitudes of all relevant
#.SEN3/ files
##############################################################################
ax = plt.axes(projection=ccrs.PlateCarree(), extent = (-161.38,-128.06,7.62,23.74))
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.5,linestyle='--')
gl.xlabels_top, gl.ylabels_left = False, False
gl.xlines, gl.ylines = False, False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

##############################################################################
#STEP4: Overplot the image data onto this projection by looping over the number
#of images to be plotted. Ranges are supplied to ensure the images don't 
#overlap
##############################################################################
ranges = {0:[500,1200,500,1500], 1:[0,950,650,1500],2:[0,950,100,900]}
geos = [geo_1,geo_2,geo_3]
F1s = [F1_1,F1_2,F1_3]
for i in range(3):
    d = ranges[i]
    lats = geos[i].variables['latitude_in'][d[0]:d[1]:4,d[2]:d[3]:4]
    lons = geos[i].variables['longitude_in'][d[0]:d[1]:4,d[2]:d[3]:4]
    F1 = F1s[i].variables[band][d[0]:d[1]:4,d[2]:d[3]:4]
    if i==0:
        Fm = F1.min()
    im = plt.pcolormesh(lons,lats,F1,vmin=Fm , vmax = F1.max(), cmap = 'RdBu')

##############################################################################
#STEP5: Add a colorbar and save the image
##############################################################################
cb = plt.axes([0.125,0.1,0.775,0.02])
cbar = plt.colorbar(im, cax=cb, orientation='horizontal')
cbar.set_label('Kelvin')
plt.savefig('Figures/Fig13.png',dpi=1000)
