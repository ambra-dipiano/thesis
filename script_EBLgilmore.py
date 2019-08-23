# =========================== #
# ADD EBL GILMORE TO TEMPLATE #
# =========================== #

# IMPORTS ---!
from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import seaborn as sns
import sys
from sys import argv

from module_analysis import fits_ebl, add_ebl

workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/run0406_test/' 
template = workdir + 'run0406_ID000126.fits'
template_ebl = workdir + 'run0406_ID000126_ebl.fits'

pathebl = '/mnt/nvme0n1p1/piano_analysis/working-dir/'
fiducial = pathebl + 'gilmore_tau_fiducial.csv'

fits_ebl(template, template_ebl, fiducial, zfetch=False, z='0.10', plot=False)

with fits.open(template_ebl) as hdul :
  print(hdul.info())
  print(hdul[4].header)

print('done')
