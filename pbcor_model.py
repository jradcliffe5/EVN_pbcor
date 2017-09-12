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
  stations["BADARY"] = 32 / 1.40
  stations["CAMBG32M"] = 32 / 1.40
  stations["EFLSBERG"] = 78 # 100
  stations["HART"] = 26/1.40
  stations["IRBENE"] = 30 / 1.40
  stations["JODRELL1"] = 67. #67 original value # 76
  stations["JODRELL2"] = 25/1.40
  stations["KUNMING"] = 40/1.40
  stations["MEDICINA"] = 25 #32/1.05
  stations["METSAHOV"] = 14/1.40
  stations["NOTO"] = 32 / 1.40
  stations["ONSALA60"] = 20/1.40
  stations["ONSALA85"] = 25/1.40
  stations["ONS_DBBC"] = stations["ONSALA85"]
  stations["ROBLEDO"] = 70/1.40
  stations["SHANGHAI"] = 22.5 # 25
  stations["SVETLOE"] = 32/1.40
  stations["TORUN"] = 32/1.40
  stations["URUMQI"] = 25/1.40
  stations["WSTRBORK"] = 25 / 1.40
  stations["YEBES40M"] = 40 / 1.40
  stations["ZELENCHK"] = 32 / 1.40
  stations["GBT_VLBA"] = 100. / 1.40
  return stations

def generate_psf(phase_centers,telescope,model): ##generate voltage beams!
    if model == 'gaussian':
        for i in range(len(phase_centers)): #gaussian of form
            print telescope, station_HPBW(stations[telescope],frequency)*60
            xstd= station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2)))*u.degree*(1/np.cos((phase_centers[i].dec.radian)))
            ystd =station_HPBW(stations[telescope],frequency)/(2*np.sqrt(2*np.log(2)))*u.degree
            if i == 0:
                EF = models.Gaussian2D(amplitude=1, \
                                       x_mean=phase_centers[i].ra.degree,\
                                       y_mean=phase_centers[i].dec.degree,\
                                       x_stddev=xstd,\
                                       y_stddev=ystd, theta=0)
            else:
                EF = EF + models.Gaussian2D(amplitude=1,\
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
multiple_pointing_names = ['P1','P2','P3','P4','HDFN']

scales = [1,1,1,1,1] ## relative observing time propto sqrt(time_per_pointing)
#--
### Central pointings ###
Positions2_RA = ['12h36m50s']
Positions2_Dec = ['+62d12m58s']

### central frequency of observations
frequency = 1.65849E9
### telescopes participating (in order of AIPS antenna table!)
## Variable telescopes
telescopes = ['EFLSBERG','WSTRBORK','ROBLEDO','ONSALA60','MEDICINA','NOTO','TORUN',\
              'SVETLOE','BADARY','ZELENCHK','URUMQI','SHANGHAI','JODRELL1']
AIPStelescope = [1,2,3,4,5,6,7,8,9,10,11,12,13]
outside_telescopes = ['EFLSBERG','JODRELL1']
outsideAIPStelescope = [1,13]

### Densities
plot_density = 500j  ## Normally don't have to change this but is sampling of the plot for gridding
pbcor_density = 1500j ## Sampling for the pbcor density, set to higher values if you have large FoV
## or rapidly changing PB. Normally ok to set to this value or lower

### PBCOR limits
## Set limits of the pbcor to be evaluated across. Will be automatic in future updates.
RA_lim = [188.3,190.3]
DEC_lim = [61.7,62.7]

#### Path to fitsfiles to extract coordinates in EG078.npy is not here
path='/net/10.0.6.249/volume1/data/radcliff/EG078B/MSSC_PBCOR/wrong_model_MSSC/Tapered_weights/'


#######################
##-------------------##
#######################
print 'GENERATING SKYCOORDINATES OF POINTINGS'
## Generate skycoordinates of outside (c) and central (c1)
# outside pointings
c = {}
for i in range(len(multiple_pointing_names)):
    c[multiple_pointing_names[i]] = SkyCoord(ra=[Positions1_RA[i]],\
                    dec=[Positions1_Dec[i]], frame ='icrs', unit=('hour','deg'))
### Need to include lists as SkyCoord now has len ability, pretty stupid bug

# central pointings
c1 = SkyCoord(ra=Positions2_RA,dec=Positions2_Dec,frame ='icrs',unit=('hour','deg'))

print 'SET UP STATIONS'
## Set up beams as Gaussians for all positions
stations = station_table()

print 'GENERATING TELSCOPE POWER BEAMS (SINGLE POINTING)'
## Generate telescope single power patterns with sky projection included.
telescope_single_psf = {}
telescope_multiple_psf = {}

### Make dictionary for single pointing telescopes
for i in range(len(telescopes)):
    if telescopes[i] not in outside_telescopes:
        telescope_single_psf[telescopes[i]] = generate_psf(c1,telescopes[i],'gaussian')

print 'GENERATING TELSCOPE POWER BEAMS (MULTIPLE POINTINGS)'

### Make dictionary for multiple pointings
for i in range(len(multiple_pointing_names)):
    telescope_multiple_psf[multiple_pointing_names[i]] = {} ### Need to initialise nested dictionary
    for j in range(len(outside_telescopes)):
        telescope_multiple_psf[multiple_pointing_names[i]][outside_telescopes[j]]= generate_psf(c[multiple_pointing_names[i]],outside_telescopes[j],'gaussian')

print 'GENERATING PHASE CENTRE POSITIONS'
## Set the positions of observations
RA= []
DEC = []
filenames = []

if os.path.exists('./pbcor.npy'):
    print 'pbcor.npy has been found, extracting RA, Dec and sourcenames'
    print 'if you have more phase centres to add, delete this and rerun'
    print 'remember to set the path in the inputs'
    RA = np.load('./pbcor.npy')[1]
    DEC = np.load('./pbcor.npy')[2]
    filenames = np.load('./pbcor.npy')[0]
else:
    print 'pbcor.npy not found, RA, DEC, and filenames extracted from %s' % path
    for file in os.listdir(path):
        if file.endswith('.fits'):
            hdu = fits.open(path+file)
            print file
            RA.append(360+float(hdu[0].header['CRVAL1']))
            DEC.append(float(hdu[0].header['CRVAL2']))
            filenames.append(file[:8])
            np.save('pbcor.npy',[filenames,RA,DEC])


## plot telescope_single_psfs power beams and non-power beams

### First plot the telescopes on the central pointing
print 'PLOTTING THE PRIMARY BEAMS OF TELESCOPES ON SINGLE POINTING'
for i in range(len(telescopes)):
    if telescopes[i] not in outside_telescopes:
        x, y = np.mgrid[RA_lim[0]:RA_lim[1]:plot_density,DEC_lim[0]:DEC_lim[1]:plot_density]
        plot = telescope_single_psf[telescopes[i]]
        EG078B_positions = SkyCoord(RA,DEC,unit='deg',frame='icrs')
        plt.figure(i,figsize=(8, 8))
        plt.pcolormesh(x,y,np.sqrt(plot(x,y))/np.max(np.sqrt(plot(x,y))),cmap='magma',vmin=0.1,vmax=1)
        plt.contour(x,y,np.sqrt(plot(x,y))/np.max(np.sqrt(plot(x,y))),levels=[0.5],colors=('w'),\
        linestyles=('-.'))
        #plt.colorbar()
        #plt.scatter(detections.ra,detections.dec)
        plt.gca().invert_xaxis()
        plt.grid(color='w')
        plt.savefig('PB_plots/'+telescopes[i]+'_single_voltage_beam.png',bbox_inches='tight',)
        plt.close('all')
for i in range(len(telescopes)):
    if telescopes[i] not in outside_telescopes:
        x, y = np.mgrid[RA_lim[0]:RA_lim[1]:plot_density,DEC_lim[0]:DEC_lim[1]:plot_density]
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

print 'PLOTTING THE PRIMARY BEAMS OF TELESCOPES ON MULTIPLE POINTING'
for j in range(len(multiple_pointing_names)):
    for i in range(len(outside_telescopes)):
        x, y = np.mgrid[RA_lim[0]:RA_lim[1]:plot_density,DEC_lim[0]:DEC_lim[1]:plot_density]
        plot = telescope_multiple_psf[multiple_pointing_names[j]][outside_telescopes[i]]
        EG078B_positions = SkyCoord(RA,DEC,unit='deg',frame='icrs')
        plt.figure(i,figsize=(8, 8))
        plt.pcolormesh(x,y,np.sqrt(plot(x,y))/np.max(np.sqrt(plot(x,y))),cmap='magma',vmin=0.1,vmax=1)
        plt.contour(x,y,np.sqrt(plot(x,y))/np.max(np.sqrt(plot(x,y))),levels=[0.5],colors=('w'),\
        linestyles=('-.'))
        #plt.colorbar()
        #plt.scatter(detections.ra,detections.dec)
        plt.gca().invert_xaxis()
        plt.grid(color='w')
        plt.savefig('PB_plots/%s_%s_single_voltage_beam.png' %  (multiple_pointing_names[j],outside_telescopes[i]),bbox_inches='tight',)
        plt.close('all')
for j in range(len(multiple_pointing_names)):
    for i in range(len(outside_telescopes)):
        x, y = np.mgrid[RA_lim[0]:RA_lim[1]:plot_density,DEC_lim[0]:DEC_lim[1]:plot_density]
        plot = telescope_multiple_psf[multiple_pointing_names[j]][outside_telescopes[i]]
        EG078B_positions = SkyCoord(RA,DEC,unit='deg',frame='icrs')
        plt.figure(i,figsize=(8, 8))
        plt.pcolormesh(x,y,plot(x,y)/np.max(plot(x,y)),cmap='magma',vmin=0.1,vmax=1)
        plt.contour(x,y,plot(x,y)/np.max(plot(x,y)),levels=[0.5],colors=('w'),linestyles=('-.'))
        #plt.colorbar()
        #plt.scatter(detections.ra,detections.dec)
        plt.gca().invert_xaxis()
        plt.grid(color='w')
        plt.savefig('PB_plots/%s_%s_single_power_beam.png' %  (multiple_pointing_names[j],outside_telescopes[i]),bbox_inches='tight',)
        plt.close('all')


print 'GRIDDING THE PRIMARY BEAMS TO EXTRACT VALUES FOR CORRECTIONS'
## Generate voltage beam corrections for CLCOR by taking sqrt of power beam corrections
corr_params = []

### Each filename need to grid the data and extract a value based upon the model
for i in range(len(filenames)):
    single_psf_corr = []
    print filenames[i]
    for j in range(len(telescopes)):
        if telescopes[j] not in outside_telescopes:
            print RA[i],DEC[i]
            if float(RA[i]) > 360.:
                RA_e = float(RA[i]) - 360.
            else:
                RA_e = float(RA[i])
            derived_corr_factor= {AIPStelescope[j]:np.sqrt(1/(telescope_single_psf[telescopes[j]](RA_e,DEC[i])))}
            single_psf_corr = single_psf_corr + [derived_corr_factor]
            print telescopes[j], derived_corr_factor
    corr_params = corr_params + [[filenames[i],RA[i],DEC[i],single_psf_corr]]

os.system('rm central_pointing_params.pckl')
f = open('central_pointing_params.pckl', 'wb')
pickle.dump(corr_params, f)
f.close()


### Each filename need to grid the data and extract a value based upon the model
def multiple_pointings_params(filenames,outside_telescopes,pointing_name,RA,DEC,telescope_multiple_psf,pbcor_density,RA_lim,DEC_lim):
    EG078B_CLCOR_corr =[]
    for i in range(len(filenames)): ### go through each phase center
        single_psf_corr = []
        print filenames[i]
        for j in range(len(outside_telescopes)): ## go through each telescope
            print RA[i],DEC[i]
            if float(RA[i]) > 360.:
                RA_e = float(RA[i]) - 360.
            else:
                RA_e = float(RA[i])
            derived_corr_factor = [{outsideAIPStelescope[j]:np.sqrt(1/(telescope_multiple_psf[pointing_name][outside_telescopes[j]](RA_e,DEC[i])))}]
            single_psf_corr = single_psf_corr + derived_corr_factor
            print outside_telescopes[j], derived_corr_factor
        EG078B_CLCOR_corr = EG078B_CLCOR_corr + [[filenames[i],RA[i],DEC[i],single_psf_corr]]
    return EG078B_CLCOR_corr

x = {}
for i in multiple_pointing_names:
    x[i] = multiple_pointings_params(filenames,outside_telescopes,i,RA,DEC,telescope_multiple_psf,\
    pbcor_density,RA_lim,DEC_lim)


os.system('rm outside_pointing_params.pckl')
f = open('outside_pointing_params.pckl', 'wb')
pickle.dump(x, f)
f.close()
print 'COMPLETE'
print 'SHOULD FIND PICKLE FILES FOR USE WITH apply_clcor_Parseltongue.py'
