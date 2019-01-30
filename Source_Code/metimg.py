from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import collections
import matplotlib

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}

matplotlib.rc('font', **font)

fname = 'Data/S3A_SL_1_RBT____20180817T185806_20180817T190106_20180817T205734_0179_034_341_2700_MAR_O_NR_002.SEN3/'

def plot_figures(figures, num, figname ,nrows = 1, ncols=1):

    fig, axeslist = plt.subplots(ncols=ncols, nrows=nrows, figsize=(num,num))
    if len(figures)>1:
        for ind,title in enumerate(figures):
            axeslist.ravel()[ind].imshow(figures[title], cmap='RdBu')
            axeslist.ravel()[ind].set_title(title)
            axeslist.ravel()[ind].set_axis_off()
    else:
        for ind,title in enumerate(figures):
            axeslist.imshow(figures[title], cmap='RdBu')
            axeslist.set_title(title)
            axeslist.set_axis_off()

    plt.savefig(figname,dpi=1000)

def tie_point_interp(fname, variable, mult, figname): # mult is either 'yes' or a number indicating which image to plot
    met = Dataset(fname+'met_tx.nc')
    vals = Dataset(fname+'S8_BT_in.nc')
    if vals.start_offset == met.start_offset:
        start_offset = 0.0
    else:
        start_offset = vals.start_offset * float(vals.resolution.split()[2]) / float(met.resolution.split()[2]) - met.start_offset
    lst = []
    dic = collections.OrderedDict()
    i_x = (np.array(range(vals.dimensions['columns'].size)) - vals.track_offset) * float(vals.resolution.split()[1]) / float(met.resolution.split()[1]) + met.track_offset
    i_y = np.array(range(vals.dimensions['rows'].size)) * float(vals.resolution.split()[2]) / float(met.resolution.split()[2]) + start_offset
    
    t_x = np.array(range(met.dimensions['columns'].size))
    t_y = np.array(range(met.dimensions['rows'].size))
    if len(met.variables[variable].shape) == 1:
        plt.plot(met.variables[variable][:])
        plt.xlabel('distance from satellite (arb)')
        plt.ylabel('pressure (Pa)')
        plt.tight_layout()
        plt.savefig(figname, dpi=1000)
    if len(met.variables[variable].shape) == 3:
        if met.variables[variable].shape[0] == 1:
            n = 0
            f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][n])
            arr = f(i_y,i_x)
            dic['{0}={1}'.format(met.variables[variable].dimensions[0],n)] = arr
            plot_figures(dic,5, figname)
        else:
            if mult != 'yes':
                f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][mult])
                arr = f(i_y,i_x)
                dic['{0}={1}'.format(met.variables[variable].dimensions[0],mult)] = arr
                plot_figures(dic,5, figname)
            if mult == 'yes':
                n=0
                while n<met.variables[variable].shape[0]:
                    f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][n])
                    arr = f(i_y,i_x)
                    dic['{0}={1}'.format(met.variables[variable].dimensions[0],n)] = arr
                    n += 1
                plot_figures(dic,8,figname, nrows=1, ncols=5)
    if len(met.variables[variable].shape) ==4:
        if met.variables[variable].shape[0] == 1:
            n1 = 0
        else:
            if mult != 'yes':
                n1 = mult
            else:
                n1=0
                while n1<met.variables[variable].shape[0]:
                    lst.append(n1)
                    n1 +=1
        if met.variables[variable].shape[1] == 1:
            n2 = 0
            if len(lst) == 0:
                f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][n1,n2])
                arr = f(i_y,i_x)
                plt.imshow(arr, cmap='RdBu')
            else:
                for i in lst:
                    f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][i,n2])
                    arr = f(i_y,i_x)
                    dic['{0}={1},{2}={3}'.format(met.variables[variable].dimensions[1],n2, met.variables[variable].dimensions[0],i)] = arr
                plot_figures(dic,13,figname,nrows=1,ncols=5)
        else:
            if mult != 'yes':
                n2 = mult
                f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][n1,n2])
                arr = f(i_y,i_x)
                plt.imshow(arr, cmap='RdBu')
            if mult == 'yes':
                n2 = 0
                while n2<met.variables[variable].shape[1]:
                    f = interpolate.RectBivariateSpline(t_y,t_x, met.variables[variable][n1,n2])
                    arr = f(i_y,i_x)
                    lst.append(arr)
                    dic['{0}={1},{2}={3}'.format(met.variables[variable].dimensions[1],n2, met.variables[variable].dimensions[0],n1)] = arr
                    n2 += 1
                fives = n2/5
                if n2%5 == 0:
                    plot_figures(dic,13,figname,fives,5)
                else:
                    plot_figures(dic,13,figname,fives+1,5)