# ======================== #
# GRB SIMULATION - ARCHIVE #
# ======================== #

# IMPORTS ---!
import gammalib
from class_analysis import *
from class_xml import *
import numpy as np
import csv
import os
import sys
from sys import argv

# --------------------------------- SETUP --------------------------------- !!!

# initialize global count ---!
chunk = int(sys.argv[1])  # global count
trials = int(sys.argv[2])  # number of trials
count = int(sys.argv[3])  # starting count

# ctools/cscripts parameters ---!
caldb = 'prod3b'
irf = 'South_z40_average_5h'

texp = [20000]  # exposure times (s)
texp.sort()
tint = len(texp)
tmin = 0  # slewing time (s)
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])

elow = 0.03  # simulation minimum energy (TeV)
ehigh = 150.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 0.5  # selection maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)

# pointing with off-axis equal to max prob GW ---!
offmax = (-1.475, -1.371)  # (deg)
trueRa = 33.057  # (deg)
trueDec = -51.841  # (deg)
pointRA = trueRa + offmax[0]  # (deg)
pointDEC = trueDec + offmax[1]  # (deg)

# others ---!
checks = True
if_fits = False
if_cut = False
if_ebl = True
skip_exist = False
extract_spec = False
fileroot = 'run0406_'
cfg_file = '/config_archive.xml'

# --------------------------------- INITIALIZE --------------------------------- !!!

cfg = xmlConfig(cfg_file)
p = cfgMng_xml(cfg)
# setup trials obj ---!
tObj = analysis(cfg_file)
tObj.pointing = [pointRA, pointDEC]
tObj.roi = roi
tObj.e = [elow, ehigh]
tObj.tmax = tmax
tObj.model = p.getWorkingDir() + 'run0406_ID000126.xml'
# add EBL to template ---!
if if_fits is True and chunk-1 == 0:
  tObj.template = p.getWorkingDir() + 'run0406_ID000126.fits' # nominal ---!
  new_template = p.getWorkingDir() + 'run0406_ID000126_ebl.fits' # absorbed ---!
  tObj.table = os.path.dirname(__file__) + 'gilmore_tau_fiducial.csv' # fiducial table ---!
  tObj.zfetch = True
  tObj.if_ebl = False
  tObj.fits_ebl(new_template)
  if_ebl = True
# assign template ---!
if if_ebl is True:
  template = p.getWorkingDir() + 'run0406_ID000126_ebl.fits'
  tObj.if_ebl = if_ebl
else :
  template = p.getWorkingDir() + 'run0406_ID000126.fits'
  tObj.if_ebl = False
tObj.template = template
print('!!! check ---- template=', tObj.template) if checks is True else None
# load template ---!
tObj.if_ebl = if_ebl
tObj.extract_spec = extract_spec
tbin_stop = tObj.load_template()
print('!!! check ---- tbin_stop=', tbin_stop) if checks is True else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

count += 1
tObj.seed = count
print('!!! check ---- seed=', tObj.seed) if checks is True else None
# attach ID to fileroot ---!
if if_ebl is True :
  f = fileroot + 'ebl%06d' % (count)
else :
  f = fileroot + 'sim%06d' % (count)
print('!!! check ---- file=', f) if checks is True else None

# --------------------------------- SIMULATION --------------------------------- !!!

print(tObj.tmax)
event_bins = []
# simulate ---!
for i in range(tbin_stop):
  if if_ebl is False:
    tObj.model = p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i
    tObj.event = p.getSimDir() + f + "_tbin%02d.fits" % i
    print('!!! check ---- simulation without EBL') if checks is True else None
  else:
    tObj.model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i
    tObj.event = p.getSimDir() + f + "_ebl_tbin%02d.fits" % i
    print('!!! check ---- simulation with EBL') if checks is True else None
  event_bins.append(tObj.event)
  if skip_exist is True:
    if not os.path.isfile(tObj.event):
      tObj.eventSim()
  else:
    print('!!! check ---tmax', tObj.tmax)
    exit()
    tObj.eventSim()
print('!!! check ---- simulation=', tObj.event) if checks is True else None
# observation list ---!
tObj.event = event_bins
tObj.event_list = p.getSimDir() + 'obs_%s.xml' % f
if skip_exist is True:
  if not os.path.isfile(tObj.event_list):
    tObj.obsList(obsname=f)
else:
  tObj.obsList(obsname=f)
print('!!! check ---- obs list=', tObj.event_list) if checks is True else None
