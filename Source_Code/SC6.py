from netCDF4 import Dataset
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figure 12 in the
    User Guide.
    Figure 12 is a comparison of images in different coordinate systems.
'''
##############################################################################
#STEP1: populate dictionaries which contain information about the file name 
#structure for each band (bands, suffix_dic) and a dictionary to indicate the 
#coordinate referenced in each coordinate system (cs)
##############################################################################
bands = {'F1':'_BT','F2':'_BT','S7':'_BT','S8':'_BT','S9':'_BT','S1':'_radiance', 'S2':'_radiance', 'S3':'_radiance','S4':'_radiance','S5':'_radiance','S6':'_radiance'}
cs={'Cartesian':['x','y'],'Geodetic':['longitude','latitude'],'Indices':['pixel','scan']}
suffix_dic = {'F1':['_in'],'F2':['_in'],'S7':['_in'],'S8':['_in'],'S9':['_in'],'S1':['_an'], 'S2':['_an'], 'S3':['_an'],'S4':['_an', '_bn','_cn'],'S5':['_an', '_bn','_cn'],'S6':['_an', '_bn','_cn']}

##############################################################################
#STEP2: create a function which takes a given axis and plots an image onto it, 
#given the coordinate system and the band
##############################################################################
def plot(fname,axis,band,system=None):
    ##########################################################################
    #STEP2a: Find the suffix by inspecting the dictionary, or asking for input.
    #Input is needed where the dictionary value is a list containing more than
    #one element. Import the image array to be plotted
    ##########################################################################
    
    if len(suffix_dic[band]) == 1:
        suffix = suffix_dic[band][0]
    else:
        suffix = raw_input('Please type the stripe and view in the form \'_<stripe><view>\': ')
    val = Dataset(fname+band+bands[band]+suffix+'.nc').variables[band+bands[band]+suffix][:]
    ##########################################################################
    #STEP2b: check if the coordinate system given is valid. If so, import 
    #relevant x/y coordinates and constrain the image to exist within those.
    #If the coordinate system is invalid or not provided, plot the regular
    #product grid
    ##########################################################################
    if system in cs:
        x = Dataset(fname+system+suffix+'.nc').variables[cs[system][0]+suffix][:]
        y = Dataset(fname+system+suffix+'.nc').variables[cs[system][1]+suffix][:]
        axis.set_xlim(x.min(),x.max())
        axis.set_ylim(y.min(),y.max())
        axis.set(xlabel=cs[system][0],ylabel=cs[system][1])
        axis.set_title(system)
        im = axis.contourf(x,y,val,cmap = 'RdBu')
    else:
        im = axis.contourf(val, cmap = 'RdBu')
        axis.set_title('Product')
    return im

##############################################################################
#STEP3: Define the path to the .SEN3/ file, create a subplot of images and fill 
#them using the function defined in STEP2
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'

fig, axes = plt.subplots(nrows=2, ncols=2)

im0 = plot(fname,axes[0,0],'F1','_in')
axes[0,0].set_xticks([])
axes[0,0].set_yticks([])
im1 = plot(fname,axes[0,1],'F1',system = 'Cartesian')
im2 = plot(fname,axes[1,0],'F1',system = 'Geodetic')
im3 = plot(fname,axes[1,1],'F1',system = 'Indices')

plt.subplots_adjust(wspace=0.5, hspace=0.65)
plt.savefig('Figures/Fig12.png', dpi=1000)