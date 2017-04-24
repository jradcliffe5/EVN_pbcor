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
## Set up observation information

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
  stations["SHANGHAI"] = 22.5 # 25
  stations["SVETLOE"] = 32/1.05
  stations["TORUN"] = 32/1.05
  stations["URUMQI"] = 25/1.05
  stations["WSTRBORK"] = 25 / 1.05
  stations["YEBES40M"] = 40 / 1.05
  stations["ZELENCHK"] = 32 / 1.05
  stations["GBT_VLBA"] = 100. / 1.05
  return stations

def generate_psf(phase_centers,telescope,model,scale):
    if model == 'gaussian':
        for i in range(len(phase_centers)):
            if i == 0:
                EF = models.Gaussian2D(amplitude=scale[i], x_mean=phase_centers[i].ra*u.degree, y_mean=phase_centers[i].dec*u.degree, x_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2))))*u.degree*(1/np.cos(1.08)),y_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2))))*u.degree,theta=0)
            else:
                EF = EF + models.Gaussian2D(amplitude=scale[i], x_mean=phase_centers[i].ra*u.degree, y_mean=phase_centers[i].dec*u.degree, x_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2))))*u.degree*(1/np.cos(1.08)),y_stddev=(station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2))))*u.degree,theta=0)
    else:
        print 'Gaussians only please'
        sys.exit()
    return EF

def convertAIPStoPythonImage(filename,outfilename):
    hdu_list = fits.open(filename)

    head = hdu_list['PRIMARY'].header
    head['CRVAL1'] = 360 - (head['CRVAL1']*-1)
    head['NAXIS'] = 2
    del head['NAXIS3']
    del head['NAXIS4']
    del head['CTYPE3'],  head['CRVAL3'], head['CDELT3'], head['CRPIX3'], head['CROTA3']
    del head['CTYPE4'], head['CRVAL4'], head['CDELT4'], head['CRPIX4'], head['CROTA4']
    #if filename.endswith('NA.fits') == False:
    #head.set('BMIN', float(hdu_list[1].data.field('BMIN')))
    #head.set('BMAJ', float(hdu_list[1].data.field('BMAJ')))
    #head.set('BPA', float(hdu_list[1].data.field('BPA')))
    image_data = hdu_list[0].data
    image_data = image_data[0,0,:,:]
    hdu = fits.PrimaryHDU(image_data, header=head)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outfilename,clobber=True)
    return outfilename


# Set up coordinate systems for observation
c = SkyCoord(ra=['12h37m20s','12h36m20s','12h36m20s','12h37m20s','12h36m50s'],dec=['+62d16m28s','+62d16m28s','+62d09m28s','+62d09m28s','+62d12m58s'], frame = 'icrs')
scales = [1,1,1,1,1./16.]
c1 = SkyCoord(ra=['12h37m20s'],dec=['+62d16m28s'])
print c.ra*u.degree
print c.dec*u.degree
## Set up beams as Gaussians for all positions
stations = station_table()

## Set up variables for the observations

frequency = 1.6E9
telescopes = ['EFLSBERG','JODRELL1','ONSALA60',"SHANGHAI","BADARY",'NOTO','SVETLOE']
telescope_single_psf = {}
for i in range(len(telescopes)):
    telescope_single_psf[telescopes[i]] = generate_psf(c,telescopes[i],'gaussian',scales)

print telescope_single_psf.values()

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

plot= []
for i in range(len(telescope_baseline_pairs.values())):
    if i == 0:
        plot = telescope_baseline_pairs.values()[i]
        print plot
    else:
        plot = plot+telescope_baseline_pairs.values()[i]

fig = plt.figure(1)
RA= []
DEC = []

## Set the positions of observations
if os.path.exists('./coordinatesDEC.npy') and os.path.exists('./coordinatesRA.npy'):
    RA = np.load('coordinatesRA.npy')
    DEC = np.load('coordinatesDEC.npy')
else:
    for file in os.listdir('./FITSImages/'):
        if file.endswith('.fits'):
            hdu = fits.open('FITSImages/'+file)
            print file
            #print 360+hdu[0].header['CRVAL1']
            RA.append(360+float(hdu[0].header['CRVAL1']))
            DEC.append(float(hdu[0].header['CRVAL2']))
            image_data = hdu[0].data
            image_data = image_data[0,0,:,:]
            np.save('coordinatesRA.npy',RA)
            np.save('coordinatesDEC.npy',DEC)

x, y = np.mgrid[188.5:190:500j,61.9:62.6:500j]
c1 = SkyCoord(RA,DEC,unit='deg',frame='icrs')
plt.pcolormesh(x,y,plot(x,y)/np.max(plot(x,y)),cmap='magma',norm=LogNorm())
plt.colorbar()
plt.scatter(c1.ra,c1.dec)
plt.gca().invert_xaxis()
plt.show()
