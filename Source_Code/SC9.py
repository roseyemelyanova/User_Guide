from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
'''
    The code below creates the images for Figures 15 and 2 in the
    User Guide.
    Figure 15 is a L0 image created from the S8 band information.
    Figure 2 is a L0 image using only one detector from the S8 band, and
    including a closer view to illustrate how scans and detector numbers are 
    related.
'''
##############################################################################
#STEP1: Define the path to the relevant .SEN3/ file and store relevant 
#information from it in variables
##############################################################################
fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'

S8 = Dataset(fname+'S8_BT_in.nc').variables['S8_BT_in'][:]
S8o = Dataset(fname+'S8_BT_in.nc').variables['S8_BT_orphan_in'][:]

scan = Dataset(fname+'indices_in.nc').variables['scan_in'][:]
pix = Dataset(fname+'indices_in.nc').variables['pixel_in'][:]
det = Dataset(fname+'indices_in.nc').variables['detector_in'][:]
scano = Dataset(fname+'indices_in.nc').variables['scan_orphan_in'][:]
pixo = Dataset(fname+'indices_in.nc').variables['pixel_orphan_in'][:]
deto = Dataset(fname+'indices_in.nc').variables['detector_orphan_in'][:]

##############################################################################
#STEP2: As this script is imported into SC11 to make use of variables above,
#set a condition which will only be met if this script is run
##############################################################################
if __name__ == "__main__":
    ##########################################################################
    #STEP3: Define a function to ungrid data from L1 information
    ##########################################################################
    
    l0_dms = [846,1199]
    def plotL0(det_no):
        ######################################################################
        #STEP3a: Create and populate a L0 template by splitting the detectors
        #and finding the seperate detector pixels and their values in S8.
        #l0_dms is defined in line 30 for simplicity. l0_S8 must exist outside
        #of the for loop to produce a full L0 image for both detectors, however
        #l0_dms definition in lines 47-49 doesn't allow this. Uncomment the
        #lines to see how it changes the image
        ######################################################################
        l0_S8 = np.zeros([l0_dms[0]*2,l0_dms[1]])
        for idet in range(det_no):
            ipix_valid = np.where(np.ma.getmaskarray(pix)==False)
            iscan_valid = np.where(np.logical_and(np.ma.getmaskarray(scan)==False,det==idet))
            iscan_valid_orphan = np.where(np.logical_and(np.ma.getmaskarray(deto)==False,deto==idet))
        
            s0 = scan[iscan_valid].min()
            #if idet == 0:
                #l0_dms = [scan[iscan_valid].max() - scan[iscan_valid].min() + 1,pix[ipix_valid].max() - pix[ipix_valid].min() + 1]
                #l0_S8 = np.zeros([l0_dms[0]*2,l0_dms[1]])
            l0_S8[((scan[iscan_valid]-s0)*2)+idet,pix[iscan_valid]] = S8[iscan_valid]
            l0_S8[((scano[iscan_valid_orphan]-s0)*2)+idet,pixo[iscan_valid_orphan]] = S8o[iscan_valid_orphan]
        return l0_S8
    ##########################################################################
    #STEP4: Call the function to produce the required plots
    ##########################################################################
    l0_S8 = plotL0(2)
    l0_S8 = np.ma.masked_values(l0_S8,0)
    im = plt.imshow(l0_S8, aspect=0.5, origin='lower', cmap='RdBu')
    plt.xlabel('pixel number')
    plt.ylabel('scan number')
    cax = plt.axes([0.9,0.12,0.02,0.76])
    plt.colorbar(im, cax)
    plt.savefig('Figures/Fig15.png',dpi=1000)
    plt.show()
    
    l0_det = plotL0(1)
    l0_det = np.ma.masked_values(l0_det,0)
    fig, axes = plt.subplots(ncols=2, nrows=1, gridspec_kw = {'width_ratios':[1, 1]})
    axes[0].imshow(l0_det, cmap='RdBu')
    axes[1].set_xlim(700,750)
    axes[1].set_ylim(625,700)
    axes[1].imshow(l0_det, cmap='RdBu', aspect=1)
    fig.text(0.02, 0.5, 'scan number', ha='center', va='center', rotation=90)
    fig.text(0.5, 0.02, 'pixel number', ha='center', va='center')
    plt.savefig('Figures/Fig2.png', dpi=1000)
    plt.show()
    