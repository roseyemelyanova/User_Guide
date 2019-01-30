#################################################################################################
FOR ALL USERS:
If you are reading this on your web browser, you need to checkout all folders and files from this
repository with svn if you want a hassle free way to recreate any data in the User Guide. If you 
are unfamiliar as to how this works, see:
https://www.tutorialspoint.com/svn/svn_checkout_process.htm
If you are using a windows machine it may be easier to use a package instead, such as TortoiseSVN:
https://tortoisesvn.net/

Once checked out onto your local machine, you are ready to begin coding!
#################################################################################################
HOW SC RELATES TO THE USER GUIDE:
Figures in the User Guide are labelled 'SC<num>' to indicate which source code document produced
them. You should refer to these documents to experiment with them or recreate their results with
your own datasets.
#################################################################################################
BREAKDOWN OF CONTENTS:
Everything is contained within the Source_Code folder.
---> 'Data' contains all .SEN3 files used throughout the User Guide
---> 'Figures' contains all Figures in the user guide, labelled 'Fig<num>' where num indicates 
      the figure number in the document.
NOTE: all figures are overwritten by the Source Code python scripts unless you change the name 
in the plt.savefig() method in each script.

There are 17 python scripts, labelled 'SC<num>' where num varies from 1 to 16. There is one 
script which is labelled 'metimg', which is imported into SC14 for testing purposes. SC11 also 
imports a script (SC9). 
#################################################################################################
GETTING STARTED:
Once you have checked out the repository, ensure to add it to your python path in whichever IDE 
you prefer to use. If you are new to python:
1. Install Anaconda (Python 2.7 is preferable) on your OS.
2. Run Spyder - you might have to find it through 'Anaconda Navigator' for the first time, 
especially if you are using a windows machine.
3. Install any relevant packages. netCDF4, cartopy and scipy are required for the scripts to run.
3. Add your checked out repository into your python path.
#################################################################################################
WHERE TO FIND OTHER DATA:
If you are a non expert user and/or do not know where to find SLSTR Level 1 products, EUMETSAT's 
CODA portal has an extensive database here: 
https://coda.eumetsat.int 

ESA also operates a data access scheme for marine data: 
https://sentinel.esa.int/web/sentinel/sentinel-data-access

NOTE: Coda works best in internet explorer to enable some functionalities (i.e. drawing a region
of interest across the world) on a windows machine. Elsewhere, chrome and firefox are preferable. 
Both of these websites require login access.
#################################################################################################
WHAT TO DO IF YOU FIND AN ISSUE OR PROBLEM, OR WANT TO SUGGEST SOMETHING:
Please report any issues/bugs to rose.yemelyanova@stfc.ac.uk

Happy coding!