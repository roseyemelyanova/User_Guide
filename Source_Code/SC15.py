import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
'''
    The code below creates the image for Figures 26 and 28 in the
    User Guide.
    Figure 26 is a plot showing the detector response at different detector
    widths.
    Figure 28 is the area covered by each detector pixel along a scan, in a 
    1D model described in the related User Guide section.
'''
##############################################################################
#STEP1: Import the relevant arrays from relevant files
##############################################################################

fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'
time = Dataset(fname+ 'time_in.nc')

##############################################################################
#STEP2: The 'time_stamp_i' variable gives information about the time elapsed
#in microseconds from the beginning of the millenium until the time the scan
#begins. ts_corr contains only scan values from the beginning of the granule
##############################################################################
ts = time.variables['time_stamp_i'][:]
ts_corr = ts - ts[0]
##############################################################################
#NOTE: the 'i' grid contains 2 detectors, so every two entries in ts_corr are
#the same. Below attempts to find only the unique scan values
##############################################################################
ts_uniq = np.unique(ts_corr)

##############################################################################
#STEP3: Identify key variables. The time difference between the first and 
#second scan divided by the pixels in one scan yields the time elapsed per 
#pixel in that scan. t_start and t_int can be found in the quality files for 
#related bands -  Here, as 'i' is used so any TIR/TIRf bands will yield the 
#same values for these variables. H and R represent height of satellite from
#the ground and radius of the earth respectively (in km)
##############################################################################
pixel_width  = float(ts_uniq[1])/3670.0
t_start = 5
t_int   = 70
t_end = pixel_width- t_start - t_int

H = 814.5
R = 6371.0 

##############################################################################
#STEP4: Define function to return the length a detector sees on the ground for
#a given detector pixel along the swath
##############################################################################
def ang_length(input_pix):
    ##########################################################################
    #STEP4a: variables are defined to store angular and distance information:
    #pix   : pixel number. SSP has a pix number of 998, the track offset.
    #theta : angular width of one pixel from the satellite
    #Omega : angle subtended at earth's center with 'pix' as the arc length
    #c     : length of circle sector with radius R & arc length 2*pix
    #h     : height of circle sector with radius R & arc length 2*pix
    #omega : angle subtended at satellite between pixel and nadir direction
    #rho   : angle between lines connecting earth center, pixel and satellite
    #phi   : angle in circle sector with radius R & arc length pix between
    #        chord length and radius. This makes up an isosceles triangle
    #Phi   : angle between nadir direction and chord length of circle sector
    #        with radius R & arc length pix
    #d     : straight line distance from satellite to pix
    #lam   : angle between chord of same circle sector as Phi and horizontal
    #L     : actual length seen by satellite in angle theta - detector length
    #alpha : angle between chord of circle sector with radius R and length L
    #l     : length of chord of circle sector with radius R and length L
    ##########################################################################
    pix = np.abs(input_pix-998)
    theta = np.arctan(0.5/H)*2 
    Omega = pix/R
    c, h = 2*R*np.sin(Omega), R*(1-np.cos(Omega))
    omega = np.arctan(c/(2*(H+h)))
    rho = np.pi - omega - Omega
    phi = (np.pi - Omega)/2
    Phi = np.pi - omega - (rho-phi)
    if Omega != 0:
        d = np.sin(Phi) * (2*R*np.sin(Omega/2))/np.sin(omega)
        lam = np.pi/2 - phi
        alpha = rho - phi - lam
        l = np.sin(theta) * (d/np.sin(np.pi-theta-alpha))
        L = R*(2*np.arcsin(l/(2*R)))
    else:
        L = 1.
    return L

##############################################################################
#STEP5: Define function to calculate the area of the detector pixel
##############################################################################
def pix_area(m, img = False, last = False):
    n_resp = 500
    det_response = np.zeros(n_resp)
    grid_square = np.zeros(n_resp)
    if img == False:
        time_step = 2.25
    else:
        time_step = 1
    det_time     = np.array(range(n_resp))/time_step - (n_resp/2)
    idx = np.where((det_time >= -(pixel_width/2)) & (det_time<=(pixel_width/2)))
    idx1 = np.where((det_time >= -((pixel_width/2)*m)+t_start) & (det_time<=((pixel_width/2)*m)-t_end))
    det_response[idx1] = 1.0
    grid_square[idx] = 1.0
    overall = np.zeros(n_resp)
    for i in range(0,n_resp):
        offset = int(i-(n_resp/2))
        this_resp = np.roll(det_response,offset)
        combined= this_resp*grid_square
        overall[i] = overall[i] + np.sum(combined)
    if img == False:
        return np.count_nonzero(overall)
    else:
        img.plot(det_time,grid_square)
        img.plot(det_time,overall/overall.max())
        xleft, xright = img.get_xlim()
        ybottom, ytop = img.get_ylim() 
        img.set_aspect(abs((xright-xleft))*0.8)
        img.set_title('detector = {0} * grid square'.format(m))
        if last == True:
            img.legend(['grid pixel', 'detector pixel'], loc='upper right', bbox_to_anchor=(0.28, -0.2))
##############################################################################
#STEP6: Find the area and create desired plots
##############################################################################
figs, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

ax1 = pix_area(1,img= ax1)
ax2 = pix_area(3.014,img= ax2, last = True)

figs.text(0.5, 0.2, r'time $(\mu s)$', ha='center', va='center')
figs.text(0.06, 0.5, 'response', ha='center', va='center', rotation='vertical')

plt.savefig('Figures/Fig26.png', dpi=1000)
plt.show()

arr = np.array(range(1200))
for i in range(1200):
    length = ang_length(i)
    arr[i] = float(pix_area(length))

plt.ylabel('Area')
plt.xlabel('Pixel Number')
plt.plot(arr, 'r-')
plt.savefig('Figures/Fig28.png', dpi=1000)
plt.show()
