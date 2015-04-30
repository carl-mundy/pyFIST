from fitsmanip import fist
import astropy.units as u

# Open up the FITS image and define the pixel scale
test = fist('ngc1068.fits', pix_scale=0.0996*u.arcsec)
# Pad the image with zeros by 10"
test.padImage(10.*u.arcsec)
# Crop image centered on the location of NGC1068 with a box size of 0.5'
test.cropImage("02 42 40.771 -00 00 47.84", 
			(u.hourangle, u.degree), 
			crop_radius=0.5*u.arcmin)
# Write the manipulated FITS image to file
test.writeImage()
# Mask regions defined in a DS9 .reg file
test.maskRegion('test_region.reg')
# Write mask to file
test.writeMask()
