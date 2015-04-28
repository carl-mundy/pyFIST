# pyFIST
Python FITS Image Special Transformations

Package to help manipulate and transform astronomical FITS images while keeping WCS information intact, including
- Pad astronomical images
- Crop astronomical images
- Produce image masks using DS9 region files
- Produce image masks using catalogues of positions and radii

Package dependencies are
- `numpy` 1.9.2
- `astropy` 1.0.1
- `pyregion` 1.1.4 (optional)

NOTE: These functions have *not* been tested extensively yet.

## Planned Functionality
1. Cropping images
    - [x] at center of image
    - [x] at specified WCS location
    - [ ] at specfied pixel coordinates
2. Padding images
    - [x] on all sides
    - [ ] on each side separately
3. Masking images
    - [x] using DS9 region files
    - [x] using catalogue WCS coordinates and radii
    - [x] using pixel coordinates and radii
    - [ ] unmasking areas
4. Writing output
    - [x] writes cropped/padded original image
    - [x] writes cropped/padded mask image
    - [ ] writed masked/cropped/padded original image
  
## Usage
The functions are used from within the `fist` class and so must be imported first. When creating the `fist` object, the path to the FITS image must be passed as well as the image's pixel scale. An example would be:
```python
from fitsmanip import fist
import astropy.units as u

ngc1068 = fist('example/ngc1068.fits', 0.0996*u.arcsec)
```
From this point onwards, several functions can be applied including cropping, padding and masking the image. Whenever an image is cropped or padded, the same transformations are applied to the mask image as well.

### Cropping
The `.cropImage()` method can be used to crop the image (or, alternatively, create a cutout) at a selected WCS location, or by default at the center of the image.
As an example, if I wanted to crop the FITS image at the location of galaxy NGC1068 (as determined by [SIMBAD](http://simbad.u-strasbg.fr/simbad/)) with a box size of 2 arcminutes, one would use the following:
```python
ngc1068.cropImage("02 42 40.771 -00 00 47.84", 
			(u.hourangle, u.degree), 
			crop_radius=1*u.arcmin)
```
By default, calling `ngc1068.cropImage()` will crop the image at its center with a box side size of 2 arcminutes.

### Padding
The `.padImage()` method can be used to pad the FITS image with zeros to some specified size. Continuing our example, to pad the image on all sides by 50 arcseconds one would use:
```python
ngc1068.padImage(50*u.arcsec)
```

### Masking
One can generate a binary mask image detailing areas that are 'masked' out. Within this script, mask images are by default defined as 1 (unity) everywhere. Areas eventually masked out are given the value of 0 (zero).

One can mask WCS regions parsed from a [DS9](http://ds9.si.edu/site/Home.html) region file through the use of the [`pyregion`](http://pyregion.readthedocs.org/en/latest/) python module as follows:
```python
ncg1068.maskRegion('example/test_region.reg')
```
This takes all the regions within the .reg file and sets their areas within the mask image to 0 (zero).

pyFIST can also take a list of RA and Dec positions and radii (in either astropy.unit or pixel formats) and then mask each object in turn. This can either be accomplished on an object by object basis, or through openining a catalogue. The former can done done with
```python
ngc1068.maskPosition(ra, dec, r)
```
where `ra`, `dec` and `r` are all astropy.unit quantities (e.g. `ra = 40.2*u.degree`). To mask objects in a catalogue, one could use the following code:
```python
ngc1068.maskCatalog('example/test_cat.tsv', 
     'RAJ2000', 
     'DEJ2000', 
     radii*u.arcsec)
```
where `radii` is a list/array of radii values in arcseconds computed somewhere else in the script.

Finally, if the image pixel coordinates are already known, these can be masked using `ngc1068.maskXY(x, y, r)` where all variables are in pixels and not astropy.unit quantities.
