"""
aia_mkmovie
============

Provides 
    1. A way to download aia level1 data from the jsoc archive using drms
    2. A wrapper combining aia files into high quality movies
    3. A GUI which allows you to select a subregion of aia and set the color bar
    4. A module to create aia_images from a list
    5. A module to create a movie from a list of images

How to use the documentation
----------------------------
Documentation is available in two locations: documentation strings in the code
and on github at https://github.com/jprchlik/aia_mkmovie

The example docstring assume `aia_mkmovie` is imported as `am`:
  >>> import aia_mkmovie as am


View the function's documentation strings using the built-in ``help`` function:

  >>> help(am.aia_mkmovie)

Available subpackages
----------------------
aia_mkmovie
    Wrapper for downloading data, selecting plot regions, image creation, and movie creation
aia_select_cutout
    Only useful when called from within aia_mkmovie. It sets the color tables and location for plotting
aia_download_files
    Downloads files in a given time range with a given cadence at a given wavelength
aia_mkimage
    Makes images from a list of images
SMEARpy
    Retrieves AIA data from local SDO archives
grab_goes_xray_flux
    Downloads and processes goes x-ray fluxes and solar wind data from NOAA
make_movie
    Makes a movie with the given dimensions and cadence

"""
from aia_mkmovie import aia_mkmovie,aia_mkimage
