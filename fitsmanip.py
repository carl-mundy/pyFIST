import numpy as np, os, pprint, time
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
try:
	import pyregion
except:
	'pyregion not available! - please check'

class fist(object):

	def __init__(self, fits_path, pix_scale, hdu_idx=0, outfile='edited.fits',
			verbose=False):
		""" Here we set up what we will need for the other functions.

			Args:
				fits_path (str): path to the FITS image file to open.
				pix_scale (float): pixel scale with an astropy unit attached. For example, if the 
					pixel scale is 0.1 ""/pix then supply 0.1*u.arcsec.
				hdu_idx (int): FITS hdu-index of the data one wishes to load.
				outfile (str): path to output image
				verbose (bool): False if no terminal output wanted, True if otherwise.
		"""

		self.fits_path 		= fits_path
		self.fits_hdulist 	= fits.open(self.fits_path, memmap=True)
		self.fits_hdu 		= self.fits_hdulist[hdu_idx]
		self.fits_wcs 		= WCS(self.fits_hdu.header)
		self.pix_scale 		= pix_scale
		self.outfile 		= outfile
		self.mask 			= False
		self.verbose		= verbose

	def makeMask(self):

		self.mask = np.ones(self.fits_hdu.data.shape, dtype=int)

	def maskRegion(self, region_file):
		""" Mask (i.e. set to zero) regions defined in a DS9 .reg region file in a 
			mask image of the same dimensions as the input image.

			Args:
				region_file (str): path to DS9 region file
		"""

		if self.mask is False:
			self.makeMask()

		region = pyregion.open(region_file)
		self.mask = np.array(np.logical_and(self.mask, 
			np.invert(region.get_mask(self.fits_hdu))), dtype=int)

		del region

	def maskXY(self, x, y, r, X=False, Y=False):
		""" Mask pixel coordinates.

			Args:
				x: image x-coordinate at center of masking area.
				y: image y-coordinate at center of masking area.
				r: radius within which to mask the image, in units of pixels.
				X: meshgrid X array. False if none passed and will be created.
				Y: meshgrid Y array. False if none passed and will be created.
		"""

		if (X is False) or (Y is False):
			X, Y = np.meshgrid(np.arange(self.fits_hdu.data.shape[0]-1), 
				np.arange(self.fits_hdu.data.shape[1]-1), sparse=True)

		self.mask[(((X-x)**2. + (Y-y)**2.) < r**2.)] = 0

		del X, Y, x, y, r

	def maskPosition(self, ra, dec, r, r_pixels=False, X=False, Y=False):
		""" Mask a WCS position in the mask image.

			Args:
				ra: astropy quantity representing the RA.
				dec: astropy quantity representing the Dec.
				r: astropy quantity or float representing the radius to mask out.
				r_pixels (bool): True if r in pixels, False if r is astropy quantity.
				X: meshgrid X array. False if none passed and will be created.
				Y: meshgrid Y array. False if none passed and will be created.
		"""

		if r_pixels is False:
			r = r.to(u.arcsec) / self.pix_scale.to(u.arcsec)

		img_xy = np.array(self.fits_wcs.wcs_world2pix(ra, dec, 0))
		self.maskXY(img_xy[0], img_xy[1], r, X, Y)

	def maskCatalog(self, cat_path, ra_col, dec_col, radii, col_units=u.degree, 
			r_pixels=False, cat_type='ascii'):
		""" Use a catalogue of positions to mask objects with radius in the list 
			of radii. 

			Args:
				cat_path (str): path to catalogue to load
				ra_col (str): column name for RA 
				dec_col (str): column name for Dec
				radii (list): list of radii for each object to mask out
				col_units: astropy units that the column data are in. Either 
					singular or a tuple of astropy units
				r_pixels (bool): True if radii are in pixels, False otherwise
				cat_type (str): astropy.table catalogue format identifier
		"""

		catalog = Table.read(cat_path, format=cat_type)
		ra, dec = catalog[ra_col], catalog[dec_col]

		if self.mask is False:
			self.makeMask()

		pos_c = SkyCoord(ra, dec, unit=col_units)

		# Make this here and pass it to other functions for speed
		X, Y = np.meshgrid(np.arange(self.fits_hdu.data.shape[0]-1), 
			np.arange(self.fits_hdu.data.shape[1]-1), sparse=True)

		t0 = time.clock()
		for obj in range(len(pos_c)):
			self.maskPosition(pos_c.ra[obj].degree, pos_c.dec[obj].degree, 
				radii[obj], r_pixels, X, Y)

		if self.verbose is True:
			print '#'*79
			print 'Finished masking {0} objects in catalog in {1:1.1f}s ({2:1.2f}s / object)'.format(
				obj+1, time.clock()-t0, (time.clock()-t0)/(obj+1))
			print '#'*79

		del X, Y, pos_c, ra, dec, catalog

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
				self.fits_hdu.data.shape[1]/2., 0))
			print crop_coords
			crop_c = SkyCoord(crop_coords[0], crop_coords[1], unit=crop_coords_unit)
		else:
			crop_c = SkyCoord(crop_coords, unit=crop_coords_unit)

		crop_c_pix = np.array(self.fits_wcs.wcs_world2pix(crop_c.ra.degree, crop_c.dec.degree, 0))
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

		new_hdu = fits.PrimaryHDU(np.array(self.mask, dtype=np.uint8), header=self.fits_hdu.header)
		new_hdu.writeto(outfile, clobber=True)

	def checkMaskXY(self, x, y):
		""" Check and return the mask value at image pixel (x,y):

			Args:
				x: x-coordinate in the image plane
				y: y-coordinate in the image plane

			Returns:
				maskvals: list of mask values at positions (x,y). 0 if outside image.
		"""

		if self.mask is False:
			self.makeMask()

		if not isinstance(x, (list, np.ndarray)):
			x = [x]
		if not isinstance(y, (list, np.ndarray)):
			y = [y]

		y = np.array(y, dtype=int)
		x = np.array(x, dtype=int)

		goodmask = (x < self.mask.shape[0])*(y < self.mask.shape[1])*(x >= 0)*(y >= 0)
		maskvals = np.zeros((x.shape), dtype=np.uint8)
		maskvals[goodmask] = self.mask[x[goodmask], y[goodmask]]

		return maskvals


	def checkMaskPosition(self, ra, dec, units=u.degree):

		""" Check mask values at a particular WCS position.

			Args:
				ra: ra astropy.unit quantity
				dec: dec astropy.unit quantity

			Returns: maskvals: list of mask values at positions (x,y). 0 if outside image.
		"""
		w_pos = SkyCoord(ra, dec, unit=units)
		px, py = self.fits_wcs.wcs_world2pix(w_pos.ra.degree, w_pos.dec.degree, 0)

		return self.checkMaskXY(px, py)

	def checkMaskCatalog(self, catpath, xxcol, yycol, units=u.degree, 
		inpixels=False, cat_type='ascii'):

		""" Check mask values at particular positions in a catalogue.

			Args:
				catpath: path to catalogue file to load in.
				xxcol: column name for ra (or x-coordinate).
				yycol: column name for dec (or y-coordinate).
				inpixels: True if xxcol, yycol values are in pixel coordinates. 
					False otherwise.
				cat_type: astropy.table format indicator. Default is ascii.

			Returns:
				maskvals: list of mask values at positions (x,y). 0 if outside image.
		"""

		catalog = Table.read(catpath, format=cat_type)
		if inpixels is False:
			return self.checkMaskPosition(catalog[xxcol], catalog[yycol], units=units)
		else:
			return self.checkMaskXY(catalog[xxcol], catalog[yycol])


### TESTING ################
test_w_pos = ["02 42 40.771 -00 00 47.84","02 42 40.771 -00 00 47.84"]
test_w = SkyCoord(test_w_pos, unit=(u.hourangle,u.degree))

test = fist('example/ngc1068.fits', 
		pix_scale=0.0996*u.arcsec,
		verbose=True)
# test.makeMask()
# radii = 5.*np.random.rand(63)
# test.maskCatalog('example/test_cat.tsv', 
# 		'RAJ2000', 
# 		'DEJ2000', 
# 		radii*u.arcsec)
# test.maskRegion('example/test_region.reg')
# test.padImage(10.*u.arcsec)
# test.cropImage()
# test.cropImage("02 42 40.771 -00 00 47.84", 
			# (u.hourangle, u.degree), 
			# crop_radius=0.5*u.arcmin)
# test.writeImage()
# test.maskStats()
# test.writeMask()
# print test.checkMaskPosition(test_w.ra.degree, test_w.dec.degree)
# print test.checkMaskXY(100,100)
print test.checkMaskPosition(40.66987917, -0.01328889)
# print test.checkMaskCatalog('example/test_cat.tsv', 
	# 'RAJ2000', 'DEJ2000',
	# cat_type='ascii')




