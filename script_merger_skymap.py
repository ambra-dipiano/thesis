import healpy as hp
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import norm
from astroquery.ned import Ned
from astropy.utils.data import download_file
from astropy.cosmology import default_cosmology
cosmo = default_cosmology.get()


filename= 'run0406_MergerID000126_skymap.fits.gz'
hpx = hp.read_map(filename)
sort=sorted(hpx, reverse=True)
cumsum = np.cumsum(sort)
index90,value90=min(enumerate(cumsum), key=lambda x: abs(x[1]-0.9))
index50,value50=min(enumerate(cumsum), key=lambda x: abs(x[1]-0.5))
npix = len(hpx)
nside = hp.npix2nside(npix)
sky_area = 4 * 180**2 / np.pi
area_pixel=sky_area / npix

# max prob. position
ipix_max = np.argmax(hpx)
theta, phi = hp.pix2ang(nside, ipix_max)
ra_max = np.rad2deg(phi)
dec_max = np.rad2deg(0.5 * np.pi - theta)
print('ra_max =',ra_max)
print('dec_max =',dec_max)


# compute area that encloses probability defined
area_probability90=area_pixel * index90
area_probability50=area_pixel * index50
print('area_probability90', area_probability90)
print('area_probability50', area_probability50)
