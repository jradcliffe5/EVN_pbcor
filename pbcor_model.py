import astropy
from astropy.coordinates import SkyCoord
from astropy.modeling import models
from scipy import constants
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.wcs import WCS
import os
from astropy.coordinates import SkyCoord
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp2d
from matplotlib.colors import LogNorm
import pickle
## Set up observation information
plt.ioff()

def station_HPBW(station,frequency):
    HPBW = ((constants.c/frequency)/station)*(180/np.pi)
    return HPBW

def station_table():
  stations = {}
  # Effective aperture size in meters for each station, currently only L band
  stations["ARECIBO"] = 206 # Effective, physical diameter = 305
  stations["BADARY"] = 32 / 1.05
  stations["CAMBG32M"] = 32 / 1.05
  stations["EFLSBERG"] = 78 # 100
  stations["HART"] = 26/1.05
  stations["IRBENE"] = 30 / 1.05
  stations["JODRELL1"] = 67. #67 original value # 76
  stations["JODRELL2"] = 25/1.05
  stations["KUNMING"] = 40/1.05
  stations["MEDICINA"] = 25 #32/1.05
  stations["METSAHOV"] = 14/1.05
  stations["NOTO"] = 32 / 1.05
  stations["ONSALA60"] = 20/1.05
  stations["ONSALA85"] = 25/1.05
  stations["ONS_DBBC"] = stations["ONSALA85"]
  stations["ROBLEDO"] = 70/1.05
  stations["SHANGHAI"] = 22.5 # 25
  stations["SVETLOE"] = 32/1.05
  stations["TORUN"] = 32/1.05
  stations["URUMQI"] = 25/1.05
  stations["WSTRBORK"] = 25 / 1.05
  stations["YEBES40M"] = 40 / 1.05
  stations["ZELENCHK"] = 32 / 1.05
  stations["GBT_VLBA"] = 100. / 1.05
  return stations

def generate_psf(phase_centers,telescope,model,scale): ##generate voltage beams!
    if model == 'gaussian':
        for i in range(len(phase_centers)): #gaussian of form
            print telescope, station_HPBW(stations[telescope],frequency)*60
            xstd= station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2)))*u.degree*(1/np.cos((phase_centers[i].dec.radian)))
            ystd =station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2)))*u.degree
            if i == 0:
                EF = models.Gaussian2D(amplitude=scale[i], \
                x_mean=phase_centers[i].ra.degree,\
                y_mean=phase_centers[i].dec.degree,\
                x_stddev=xstd,\
                y_stddev=ystd, theta=0)
            else:
                EF = EF + models.Gaussian2D(amplitude=scale[i],\
                x_mean=phase_centers[i].ra.degree, \
                y_mean=phase_centers[i].dec.degree, \
                x_stddev=xstd, \
                y_stddev=ystd, theta=0)
    else:
        print 'Gaussians only please'
        sys.exit()
    return EF

#######################
######## Inputs #######
#######################
# Set up coordinate systems for observation
### Multiple pointings ###
Positions1_RA = ['12h37m20s','12h36m20s','12h36m20s','12h37m20s','12h36m50s']
Positions1_Dec = ['+62d16m28s','+62d16m28s','+62d09m28s','+62d09m28s','+62d12m58s']
scales = [1,1,1,1,1/16] ## relative observing time propto sqrt(time_per_pointing)
#--
### Central pointings ###
Positions2_RA = ['12h36m50s']
Positions2_Dec = ['+62d12m58s']

### central frequency of observations
frequency = 1.65849E9
### telescopes participating (in order of AIPS antenna table!)
telescopes = ['EFLSBERG','WSTRBORK','ROBLEDO','ONSALA60',"MEDICINA",'NOTO',"TORUN",\
              'SVETLOE',"BADARY","ZELENCHK","URUMQI","SHANGHAI","JODRELL1"]
outside_telescopes = ['EFLSBERG','JODRELL1']

#### Path to fitsfiles to extract corrdinates in EG078.npy is not here
path='../../wrong_model_MSSC/Tapered_weights/'

#######################
##-------------------##
#######################

## Generate skycoordinates of outside (c) and central (c1)
# outside pointings
c = SkyCoord(ra=Positions1_RA,\
            dec=Positions1_Dec, frame = 'icrs',unit=('hour','deg'))
# central pointings
c1 = SkyCoord(ra=Positions2_RA,dec=Positions2_Dec,frame ='icrs',unit=('hour','deg'))
print c.ra*u.degree
print c.dec*u.degree


## Set up beams as Gaussians for all positions
stations = station_table()

## Generate telescope single power patterns with sky projection included.
telescope_single_psf = {}
for i in range(len(telescopes)):
    if telescopes[i] in outside_telescopes:
        telescope_single_psf[telescopes[i]] = generate_psf(c,telescopes[i],'gaussian',scales)
    else:
        telescope_single_psf[telescopes[i]] = generate_psf(c1,telescopes[i],'gaussian',scales)


## Set the positions of observations
RA= []
DEC = []
filenames = []

if os.path.exists('./EG078B.npy'):
    RA = np.load('./EG078B.npy')[1]
    DEC = np.load('./EG078B.npy')[2]
    filenames = np.load('./EG078B.npy')[0]
else:
    for file in os.listdir(path):
        if file.endswith('.fits'):
            hdu = fits.open(path+file)
            print file
            #print 360+hdu[0].header['CRVAL1']
            RA.append(360+float(hdu[0].header['CRVAL1']))
            DEC.append(float(hdu[0].header['CRVAL2']))
            filenames.append(file[:8])
            np.save('EG078B.npy',[filenames,RA,DEC])
'''
## plot telescope_single_psfs power beams and non-power beams
import pandas as pd
#detections = pd.read_csv('VLBI_Catalogue_v10.csv',delimiter='\t')
#detections = SkyCoord(detections.VLBI_RA,detections.VLBI_Dec,frame='icrs',unit=('hour','deg'))
for i in range(len(telescopes)):
    x, y = np.mgrid[188.3:190.3:500j,61.7:62.7:500j]
    plot = telescope_single_psf[telescopes[i]]
    EG078B_positions = SkyCoord(RA,DEC,unit='deg',frame='icrs')
    plt.figure(i,figsize=(8, 8))
    plt.pcolormesh(x,y,np.sqrt(plot(x,y))/np.max(np.sqrt(plot(x,y))),cmap='magma',vmin=0.1,vmax=1)
    plt.contour(x,y,np.sqrt(plot(x,y))/np.max(np.sqrt(plot(x,y))),levels=[0.5],colors=('w'),linestyles=('-.'))
    #plt.colorbar()
    #plt.scatter(detections.ra,detections.dec)
    plt.gca().invert_xaxis()
    plt.grid(color='w')
    plt.savefig('PB_plots/'+telescopes[i]+'_single_voltage_beam.png',bbox_inches='tight',)
    plt.close('all')

for i in range(len(telescopes)):
    x, y = np.mgrid[188.3:190.3:500j,61.7:62.7:500j]
    plot = telescope_single_psf[telescopes[i]]
    EG078B_positions = SkyCoord(RA,DEC,unit='deg',frame='icrs')
    plt.figure(i,figsize=(8, 8))
    plt.pcolormesh(x,y,plot(x,y)/np.max(plot(x,y)),cmap='magma',vmin=0.1,vmax=1)
    plt.contour(x,y,plot(x,y)/np.max(plot(x,y)),levels=[0.5],colors=('w'),linestyles=('-.'))
    #plt.colorbar()
    #plt.scatter(detections.ra,detections.dec)
    plt.gca().invert_xaxis()
    plt.grid(color='w')
    plt.savefig('PB_plots/'+telescopes[i]+'_single_power_beam.png',bbox_inches='tight',)
    plt.close('all')
'''
## Generate voltage beam corrections for CLCOR by taking sqrt of power beam corrections
EG078B_CLCOR_corr = []

### Each filename need to grid the data and extract a value based upon the model
for i in range(len(filenames)):
    single_psf_corr = []
    print filenames[i]
    for j in range(len(telescopes)):
        x, y = np.mgrid[188.5:190:500j,61.9:62.6:500j]
        derived_corr_factor = [1/(telescope_single_psf[telescopes[j]](RA[i],DEC[i])/np.max(telescope_single_psf[telescopes[j]](x,y)))]
        single_psf_corr = single_psf_corr + derived_corr_factor
        print telescopes[j], derived_corr_factor
    EG078B_CLCOR_corr = EG078B_CLCOR_corr + [[filenames[i],RA[i],DEC[i],np.sqrt(single_psf_corr)]]

os.system('rm CLCOR_params.pckl')
f = open('CLCOR_params.pckl', 'wb')
pickle.dump(EG078B_CLCOR_corr, f)
f.close()

'''
##Generate baseline pair voltage beams
telescope_baseline_pairs = {}
key = telescope_single_psf.keys()
values = telescope_single_psf.values()

for i in range(len(values)):
    for j in range(len(values)):
        if key[i] != key[j]:
            if key[j]+key[i] not in telescope_baseline_pairs :
                telescope_baseline_pairs[key[i]+key[j]] = values[i] * values[j]
            else:
                print key[i]+' baseline '+key[j]+' already calculated'
        else:
            print key[i]+' is the same as '+key[j]
print telescope_baseline_pairs.keys()

#combine baseline beams to form primary beam
plot= []
for i in range(len(telescope_baseline_pairs.values())):
    if i == 0:
        plot = telescope_baseline_pairs.values()[i]
        print plot
    else:
        plot = plot+telescope_baseline_pairs.values()[i]


## Plots for combining the arrays
fig = plt.figure(1)
x, y = np.mgrid[188.5:190:500j,61.9:62.6:500j]
c1 = SkyCoord(RA,DEC,unit='deg',frame='icrs')
plt.figure(1)
plt.pcolormesh(x,y,plot(x,y)/np.max(plot(x,y)),cmap='magma',norm=LogNorm())
plt.colorbar()
plt.scatter(c1.ra,c1.dec)
plt.gca().invert_xaxis()
plt.show()

## Pwer beam plot
plt.figure(2)
plt.pcolormesh(x,y,(plot(x,y)/np.max(plot(x,y)))**2,cmap='viridis')
plt.colorbar()
plt.scatter(c1.ra,c1.dec)
plt.gca().invert_xaxis()
plt.show()
'''
