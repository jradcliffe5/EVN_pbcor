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
            if i == 0:
                EF = models.Gaussian2D(amplitude=scale[i], \
                x_mean=phase_centers[i].ra*u.degree,\
                y_mean=phase_centers[i].dec*u.degree,\
                x_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(np.log(2))))*u.degree*(1/np.cos(1.08)),\
                y_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(np.log(2))))*u.degree,theta=0)
            else:
                EF = EF + models.Gaussian2D(amplitude=scale[i], x_mean=phase_centers[i].ra*u.degree, y_mean=phase_centers[i].dec*u.degree, x_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(np.log(2))))*u.degree*(1/np.cos(1.08)),y_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(np.log(2))))*u.degree,theta=0)
    else:
        print 'Gaussians only please'
        sys.exit()
    return EF


# Set up coordinate systems for observation
c = SkyCoord(ra=['12h37m20s','12h36m20s','12h36m20s','12h37m20s','12h36m50s'],\
            dec=['+62d16m28s','+62d16m28s','+62d09m28s','+62d09m28s','+62d12m58s'], frame = 'icrs') # outside pointings
scales = [1,1,1,1,1/16]
c1 = SkyCoord(ra=['12h36m50s'],dec=['+62d12m58s'],frame = 'icrs') ## central pointings
print c.ra*u.degree
print c.dec*u.degree
## Set up beams as Gaussians for all positions
stations = station_table()

## Set up variables for the observations

frequency = 1.65849E9 ## central frequency of observations
telescopes = ['EFLSBERG','WSTRBORK','ROBLEDO','ONSALA60',"MEDICINA",'NOTO',"TORUN",\
              'SVETLOE',"BADARY","ZELENCHK","URUMQI","SHANGHAI","JODRELL1"] #telescopes participating

## Generate telescope single voltage patterns with sky projection included.
telescope_single_psf = {}
for i in range(len(telescopes)):
    if telescopes[i] == 'EFLSBERG' or telescopes[i] == 'JODRELL1':
        telescope_single_psf[telescopes[i]] = generate_psf(c,telescopes[i],'gaussian',scales)
    else:
        telescope_single_psf[telescopes[i]] = generate_psf(c1,telescopes[i],'gaussian',scales)


## Set the positions of observations
RA= []
DEC = []
filenames = []
path='/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/wrong_model_MSSC/Tapered_weights/'
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

## plot telescope_single_psf
'''
for i in range(len(telescopes)):
    x, y = np.mgrid[188.5:190:500j,61.9:62.6:500j]
    plot = telescope_single_psf[telescopes[i]]
    EG078B_positions = SkyCoord(RA,DEC,unit='deg',frame='icrs')
    plt.figure(i)
    plt.pcolormesh(x,y,plot(x,y)/np.max(plot(x,y)),cmap='magma')
    plt.colorbar()
    plt.scatter(EG078B_positions.ra,EG078B_positions.dec)
    plt.gca().invert_xaxis()
    #plt.show()
    plt.savefig(telescopes[i]+'_single_primary_beam.pdf',bbox_inches='tight')
    plt.close('all')
'''

## Generate voltage beam corrections for CLCOR
EG078B_CLCOR_corr = []

for i in range(len(filenames)):
    single_psf_corr = []
    for j in range(len(telescopes)):
        x, y = np.mgrid[188.5:190:500j,61.9:62.6:500j]
        single_psf_corr = single_psf_corr + [1/(telescope_single_psf[telescopes[j]](RA[i],DEC[i])/np.max(telescope_single_psf[telescopes[j]](x,y)))]
        print telescopes[j]
        print 1/(telescope_single_psf[telescopes[j]](RA[i],DEC[i])/np.max(telescope_single_psf[telescopes[j]](x,y)))
    EG078B_CLCOR_corr = EG078B_CLCOR_corr + [[filenames[i],RA[i],DEC[i],single_psf_corr]]

os.system('rm CLCOR_params.pckl')
f = open('CLCOR_params.pckl', 'wb')
pickle.dump(EG078B_CLCOR_corr, f)
f.close()

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

'''
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
