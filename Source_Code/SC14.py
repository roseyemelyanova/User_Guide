import metimg
'''
    The code below creates the images for Figures 23, 24 and 25 in the
    User Guide.
    Figure 23 is a plot of the atmospheric pressure as a function of height 
    form the satellite.
    Figure 24 shows predicted images of the specific humidity at these 
    atmospheric pressures.
    Figure 25 is the predicted wind speed at various times.
'''
##############################################################################
#THIS SCRIPT USES METIMG.PY - THE FOLLOWING ONLY CALLS ITS VARIABLES AND ITS
#METHODS
##############################################################################

fname = metimg.fname

metimg.tie_point_interp(fname,'p_atmos','yes','Figures/Fig23.png')
metimg.tie_point_interp(fname,'temperature_profile_tx','yes','Figures/Fig24.png')
metimg.tie_point_interp(fname,'u_wind_tx','yes','Figures/Fig25.png')