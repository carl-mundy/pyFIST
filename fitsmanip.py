import numpy as np, os, pprint
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS 
from astropy.coordinates import SkyCoord
import pyregion

class fits_manip(object):

	def __init__(self, fits_path, pix_scale, hdu_idx=0, outfile='edited.fits'):
		""" Here we set up what we will need for the other functions.

			Args:
				fits_path (str): path to the FITS image file to open.
				pix_scale (float): pixel scale with an astropy unit attached. For example, if the 
					pixel scale is 0.1 ""/pix then supply 0.1*u.arcsec.
				hdu_idx (int): FITS hdu-index of the data one wishes to load.
				outfile (str): path to output image
		"""

		self.fits_path 		= fits_path
		self.fits_hdulist 	= fits.open(self.fits_path, memmap=True)
		self.fits_hdu 		= self.fits_hdulist[hdu_idx]
		self.fits_wcs 		= WCS(self.fits_hdu.header)
		self.pix_scale 		= pix_scale
		self.outfile 		= outfile
		self.mask 			= False

	def info(self):

		pprint.pprint(self.fits_hdu.header)
		pprint.pprint(self.fits_wcs.wcs.info())

	def makeMask(self):

		self.mask = np.ones(self.fits_hdu.data.shape, dtype=np.int8)

	def maskRegion(self, region_file):

		if self.mask is False:
			self.makeMask()

		region = pyregion.open(region_file)
		self.mask = np.logical_and(self.mask, 
			np.invert(region.get_mask(self.fits_hdu)))

		del region

	def maskStats(self):
		mask_total = np.product(self.mask.shape, dtype=float)
		print mask_total
		mask_good = np.sum(self.mask, dtype=float)
		print mask_good
		print mask_good / mask_total

	def cropImage(self, crop_coords=False, crop_coords_unit=u.degree, crop_radius=1*u.arcmin):
		""" Crop the FITS image, centered on some coordinates.

			Args:
				crop_coords: Coordinate string, or list to pass to astropy.coordinates.SkyCoord
					function.
				crop_coords_unit: Astropy unit specifiying the units of crop_coords.
				crop_radius: Astropy quantity specifying radius around object to crop.
		"""

		if crop_coords is False:
			# crop_coords = [self.fits_hdu.header['CRVAL1'], self.fits_hdu.header['CRVAL2']]
			crop_coords = np.array(self.fits_wcs.wcs_pix2world(self.fits_hdu.data.shape[0]/2., 
				self.fits_hdu.data.shape[1]/2., 1))
			print crop_coords
			crop_c = SkyCoord(crop_coords[0], crop_coords[1], unit=crop_coords_unit)
		else:
			crop_c = SkyCoord(crop_coords, unit=crop_coords_unit)

		crop_c_pix = np.array(self.fits_wcs.wcs_world2pix(crop_c.ra.degree, crop_c.dec.degree, 1))
		crop_radius_pixels = crop_radius.to(u.arcsec) / self.pix_scale.to(u.arcsec)

		x1 = np.clip(crop_c_pix[0]-crop_radius_pixels, 0, self.fits_hdu.data.shape[0]-1)
		x2 = np.clip(crop_c_pix[0]+crop_radius_pixels, 0, self.fits_hdu.data.shape[0]-1)
		y1 = np.clip(crop_c_pix[1]-crop_radius_pixels, 0, self.fits_hdu.data.shape[1]-1)
		y2 = np.clip(crop_c_pix[1]+crop_radius_pixels, 0, self.fits_hdu.data.shape[1]-1)

		self.fits_hdu.data = self.fits_hdu.data[y1:y2, x1:x2]

		if self.mask is not False:
			self.mask = self.mask[y1:y2, x1:x2]

		# Update WCS information
		self.fits_hdu.header['CRPIX1'] = (self.fits_hdu.data.shape[0]-0.5)/2.
		self.fits_hdu.header['CRPIX2'] = (self.fits_hdu.data.shape[1]-0.5)/2.
		self.fits_hdu.header['CRVAL1'] = crop_c.ra.degree
		self.fits_hdu.header['CRVAL2'] = crop_c.dec.degree
		self.fits_wcs = WCS(self.fits_hdu.header)

	def padImage(self, padsize):
		""" Pad a FITS image on all sides.

			Args:
				padsize: Astropy quantity detailing size of padding to apply.
		"""

		# Get required header information
		CRPIX1, CRPIX2 = self.fits_hdu.header['CRPIX1'], self.fits_hdu.header['CRPIX2']
		CRVAL1, CRVAL2 = self.fits_hdu.header['CRVAL1'], self.fits_hdu.header['CRVAL2']

		# Generate new image
		padpixels = (padsize.to(u.arcsec) / self.pix_scale.to(u.arcsec)).value
		newdims = (self.fits_hdu.data.shape[0]+(2*padpixels), 
			self.fits_hdu.data.shape[1]+(2*padpixels))
		paddedimg = np.zeros(newdims)
		paddedimg[padpixels:-padpixels, padpixels:-padpixels] = self.fits_hdu.data

		if self.mask is not False:
			newmask = np.zeros(newdims)
			newmask[padpixels:-padpixels, padpixels:-padpixels] = self.mask
			self.mask = newmask
			del newmask

		# Update WCS information
		self.fits_hdu.data = paddedimg
		self.fits_hdu.header['CRPIX1'] = CRPIX1 + padpixels
		self.fits_hdu.header['CRPIX2'] = CRPIX2 + padpixels
		self.fits_wcs = WCS(self.fits_hdu.header)

	def writeImage(self, outfile=False):

		if outfile is False:	
			outfile = self.outfile

		new_hdu = fits.PrimaryHDU(self.fits_hdu.data, header=self.fits_hdu.header)
		new_hdu.writeto(outfile, clobber=True)
		del new_hdu

	def writeMask(self, outfile=False):

		if outfile is False:
			outfile = self.outfile+'.mask.fits'

		new_hdu = fits.PrimaryHDU(self.mask, header=self.fits_hdu.header)
		new_hdu.writeto(outfile, clobber=True)


### TESTING ################
test = fits_manip('ngc1068.fits', 
			pix_scale=0.0996*u.arcsec)
test.makeMask()
test.maskRegion('test_region.reg')
# test.padImage(10.*u.arcsec)
test.cropImage()
# test.cropImage("02 42 40.771 -00 00 47.84", 
			# (u.hourangle, u.degree), 
			# crop_radius=0.5*u.arcmin)
test.writeImage()
test.maskStats()
# test.writeMask()
