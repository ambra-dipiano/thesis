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
caldb_degraded = caldb.replace('prod', 'degr')
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
roi = 5  # region of interest for simulation and selection (deg)

# pointing with off-axis equal to max prob GW ---!
offmax = (-1.475, -1.371)  # (deg)
trueRa = 33.057  # (deg)
trueDec = -51.841  # (deg)
pointRA = trueRa + offmax[0]  # (deg)
pointDEC = trueDec + offmax[1]  # (deg)

# conditions control ---!
checks = True
if_fits = False
if_ebl = True
extract_spec = True
irf_degrade = False
skip_exist = False

# files ---!
fileroot = 'run0406_'
cfg_file = '/config_archive.xml'
ebl_table = '/mnt/nvme0n1p1/piano_analysis/working-dir/gilmore_tau_fiducial.csv'
nominal_template = 'run0406_ID000126.fits'
ebl_template = 'run0406_ID000126_ebl.fits'
model_pl = 'run0406_ID000126.xml'
tcsv = 'time_slices.csv'

# --------------------------------- INITIALIZE --------------------------------- !!!

cfg = xmlConfig(cfg_file)
p = cfgMng_xml(cfg)
# setup trials obj ---!
tObj = analysis(cfg_file)
tObj.pointing = [pointRA, pointDEC]
tObj.roi = roi
tObj.e = [elow, ehigh]
tObj.tmax = tmax
tObj.model = p.getWorkingDir() + model_pl
tObj.caldb = caldb
tObj.irf = irf
# degrade IRF if required ---!
if irf_degrade:
  tObj.degradeIRF()
# add EBL to template ---!
if if_fits:
  tObj.template = p.getWorkingDir() + nominal_template # nominal ---!
  new_template = p.getWorkingDir() + ebl_template # absorbed ---!
  tObj.table = ebl_table # fiducial ---!
  tObj.zfetch = True
  tObj.if_ebl = False
  tObj.fits_ebl(new_template)
  if_ebl = True
# assign template ---!
if if_ebl:
  template = p.getWorkingDir() + ebl_template
  tObj.if_ebl = if_ebl
else :
  template = p.getWorkingDir() + nominal_template
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
if if_ebl:
  f = fileroot + 'ebl%06d' % (count)
else :
  f = fileroot + 'sim%06d' % (count)
print('!!! check ---- file=', f) if checks is True else None

# --------------------------------- SIMULATION --------------------------------- !!!

event_bins = []
tObj.table = p.getDataDir() + tcsv
time = tObj.getTimeSlices()  # methods which returns time slice edges
# simulate ---!
for i in range(tbin_stop):
  tObj.t = [time[i], time[i + 1]]
  if if_ebl:
    tObj.model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i
    tObj.event = p.getSimDir() + f + "_ebl_tbin%02d.fits" % i
    print('!!! check ---- simulation with EBL') if checks is True else None
  else:
    tObj.model = p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i
    tObj.event = p.getSimDir() + f + "_tbin%02d.fits" % i
    print('!!! check ---- simulation without EBL') if checks is True else None
  event_bins.append(tObj.event)
  if skip_exist:
    if not os.path.isfile(tObj.event):
      tObj.eventSim()
  else:
    print('!!! check ---tmax', tObj.tmax)
    tObj.eventSim()
print('!!! check ---- simulation=', tObj.event) if checks is True else None
# observation list ---!
tObj.event = event_bins
tObj.event_list = p.getSimDir() + 'obs_%s.xml' % f
if skip_exist:
  if not os.path.isfile(tObj.event_list):
    tObj.obsList(obsname=f)
else:
  tObj.obsList(obsname=f)
print('!!! check ---- obs list=', tObj.event_list) if checks is True else None




