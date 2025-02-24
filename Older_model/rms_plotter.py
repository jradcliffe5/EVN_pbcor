from astropy.io import fits
from astropy.wcs import WCS
import os
from astropy.coordinates import SkyCoord, Angle
import numpy as np
import matplotlib
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from regions import CircleSkyRegion
from matplotlib import *
from matplotlib import rc
from matplotlib import rcParams

rc('font', **{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)
rcParams['mathtext.default'] = 'regular'
figsize = plt.rcParams["figure.figsize"]
figsize[1]=9
figsize[0]=9
plt.rcParams["figure.figsize"]=figsize
matplotlib.rcParams.update({'font.size': 22})
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
    hdulist.writeto(outfilename,overwrite=True)
    return outfilename

def grid(x, y, z, resX=1000, resY=1000):
    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi,interp='linear')
    Z2 = ndimage.gaussian_filter(Z, sigma=1.0, order=1)
    #Z = interp2d(x, y, z, kind='cubic')
    #Z2 = Z(xi,yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z

## Create array of ra, dec, & intensity
RA= []
DEC = []
rms = []
os.system('rm *Py.fits')
os.system('rm *Py.fits')
f = []
for file in os.listdir('./'):
    if file.endswith('.fits'):
        hdu = fits.open(file)
        print file
	f.append(file)
        #print 360+hdu[0].header['CRVAL1']
        RA.append(360+float(hdu[0].header['CRVAL1']))
        DEC.append(float(hdu[0].header['CRVAL2']))
        image_data = hdu[0].data
        image_data = image_data[0,0,:,:]
        rms.append(np.sqrt(np.mean(np.square(image_data[150:1200,150:1200])))*1E6)
#print RA, DEC, rms
#print rms
#print np.nonzero(rms > 10E-6)
#print rms[np.argmax(rms>100E-6)]
#print np.mean(rms), np.std(rms)
c = SkyCoord(RA,DEC,unit='deg',frame='icrs')

### Create axes & scatter plot
filename = 'HDFC0155_PBCOR_NA_IM.fits' ## input fits file for correct wcs definition
outfilename = 'HDFC0155_NA_PBCOR_IMPy.fits'
convertAIPStoPythonImage(filename,outfilename)
hdu = fits.open(outfilename)
wcs = WCS(hdu[0].header)
fig = plt.figure()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
#ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
lim = 500000
shiftx = -50000
shifty = -10000
ax.set_xlim((lim*-1)+shiftx, lim+shiftx)
ax.set_ylim((lim*-1)+shifty, lim+shifty)
#ax1.set_xlim(-500000, 500000)HDFC0155_PBCOR_IM.fits
#ax1.set_ylim(-500000, 500000)
X, Y, Z = grid(c.ra.degree, c.dec.degree, rms)
Z_gauss = np.nan_to_num(np.array(gaussian_filter(Z,12)))
Z_gauss[Z_gauss==0] = 50
im = ax.pcolormesh(X,Y,Z_gauss,transform=ax.get_transform('icrs'),cmap='magma',vmin=np.amin(rms), vmax=40, alpha=1)
ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'), marker='+',color='g',norm=matplotlib.colors.LogNorm(),alpha=1)
CS = ax.contour(X,Y,gaussian_filter(Z,10),levels=np.linspace(np.amin(rms),40,6),colors='w',transform=ax.get_transform('icrs'), interpolation='none')
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0.00,axes_class=matplotlib.axes.Axes)
cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax,format='%d',ticks=np.linspace(np.amin(rms),40,6),extend='max')
cb.add_lines(CS)
cb.ax.xaxis.set_ticks_position('top')
region1 = CircleSkyRegion(SkyCoord('12h36m50.0s', '+62d12m58.00s', frame='icrs'), Angle(0.0625, 'deg'))
pixel_region1 = region1.to_pixel(wcs)
region2 = CircleSkyRegion(SkyCoord('12h37m20.0s', '+62d16m28.00s', frame='icrs'), Angle(0.0625, 'deg'))
pixel_region2 = region2.to_pixel(wcs)
region3 = CircleSkyRegion(SkyCoord('12h36m20.0s', '+62d16m28.00s', frame='icrs'), Angle(0.0625, 'deg'))
pixel_region3 = region3.to_pixel(wcs)
region4 = CircleSkyRegion(SkyCoord('12h36m20.0s', '+62d09m28.00s', frame='icrs'), Angle(0.0625, 'deg'))
pixel_region4 = region4.to_pixel(wcs)
region5 = CircleSkyRegion(SkyCoord('12h37m20.0s', '+62d09m28.00s', frame='icrs'), Angle(0.0625, 'deg'))
pixel_region5 = region5.to_pixel(wcs)
pixel_region1.plot(ax,facecolor='none', edgecolor='red',linestyle='--')
pixel_region2.plot(ax,facecolor='none', edgecolor='red',linestyle='--')
pixel_region3.plot(ax,facecolor='none', edgecolor='red',linestyle='--')
pixel_region4.plot(ax,facecolor='none', edgecolor='red',linestyle='--')
pixel_region5.plot(ax,facecolor='none', edgecolor='red',linestyle='--')
cax.set_xlabel("1$\sigma$ r.m.s. sensitivity ($\mu$Jy/bm)", labelpad=-80)
#print np.amin(rms), np.amax(rms)
fig.savefig('rms_pbcor.png',bbox_inches='tight',dpi=fig.dpi,format="png")

np.save('rms_pbcor.npy',[f,c.ra.degree,c.dec.degree,rms])
#plt.show()
'''
c = SkyCoord(RA,DEC,unit='deg',frame='icrs')
hdu = fits.open(outfilename)
wcs = WCS(hdu[0].header)
rms = 1/(np.array(rms)/5.5e-06)
print np.amin(rms), np.amax(rms)
fig = plt.figure()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
#ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
ax.set_xlim(-1100000, 1100000)
ax.set_ylim(-1100000, 1100000)
#ax1.set_xlim(-500000, 500000)HDFC0155_PBCOR_IM.fits
#ax1.set_ylim(-500000, 500000)
X, Y, Z = grid(c.ra.degree, c.dec.degree, rms)
im = ax.pcolormesh(X,Y,Z,transform=ax.get_transform('icrs'),cmap='magma',vmin=0.3, vmax=1)
ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'), marker='+',color='g',norm=matplotlib.colors.LogNorm())
#CS = ax.contour(X,Y,gaussian_filter(Z,10),levels=np.linspace(0,1,8),colors='w',transform=ax.get_transform('icrs'), interpolation='none')
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0.00,axes_class=matplotlib.axes.Axes)
cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax,format='%.1e',ticks=np.linspace(0.5,1,8),extend='min')
#cb.add_lines(CS)
cb.ax.xaxis.set_ticks_position('top')
cax.set_xlabel("1$\sigma$ r.m.s. sensitivity ($\mu$Jy/bm)", labelpad=-60)

#fig.savefig('rms_taper_cal_weights_PB.pdf',bbox_inches='tight',dpi=50,format="pdf")
plt.show()
'''
#above code tranforms AIPS into python and then sets out an axis based upon the wcs coords of the fits file
# need to get a fits file that is in the center of the field to replace HDFA....fits
# also need to set up an array of coordinates & rms of the fitsfiles to feed into an
# ax.scatter and set up the colour scale to reflect r.m.s.
# Also need to interpolate look at http://cgcooke.github.io/Scattered-Interpolation/
