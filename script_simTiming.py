# ======================== #
# GRB SIMULATION - ARCHIVE #
# ======================== #

# IMPORTS ---!
import gammalib
import ctools
import cscripts
from module_analysis import *
import csv
import os
import sys
from sys import argv

# initialize global count ---!
chunk = int(sys.argv[1])  # global count
trials = int(sys.argv[2])  # number of trials
count = int(sys.argv[3])  # starting count

# work with absolute paths ---!
workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/run0406/'
runpath = workdir + 'run0406_ID000126/'
simpath = runpath + 'sim/'
selectpath = runpath + 'selected_sim/'
datapath = runpath + 'data/'
csvpath = runpath + 'csv/'
detpath = runpath + 'detection_all/'

# inputs ---!
template = workdir + 'run0406_ID000126_ebl.fits'
model = workdir + 'run0406_ID000126.xml'

# ctools/cscripts parameters ---!
caldb = 'prod3b'
irf = 'South_z40_average_100s'

sigma = 5  # detection acceptance (Gaussian)
texp = [2e4]  # exposure times (s)
texp.sort()
tint = len(texp)
tmin = 30  # slewing time (s)
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])

elow = 0.03  # simulation minimum energy (TeV)
ehigh = 1.0  # simulation maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)

# pointing with off-axis equal to max prob GW ---!
offmax = (-1.475, -1.371)  # (deg)
trueRa = 33.057  # (deg)
trueDec = -51.841  # (deg)
pointRA = trueRa + offmax[0]  # (deg)
pointDEC = trueDec + offmax[1]  # (deg)

# others ---!
checks = False
if_ebl = True
if if_ebl is False:
  fileroot = 'run0406_'
else:
  fileroot = 'run0406ebl_'

# =====================
# !!! LOAD TEMPLATE !!!
# =====================

t, tbin_stop = load_template(template, tmax, extract_spec=False, model=model, pathout=datapath)
print('!!! check ---- tbin_stop=', tbin_stop) if checks is True else None

for k in range(trials):
  count += 1
  # attach ID to fileroot ---!
  f = fileroot + 'sim%06d' % (count)
  print('!!! check ---- file=', f) if checks is True else None

  # ====================
  # !!! SIMULATE GRB !!!
  # ====================

  event_bins = []

  # simulate ---!
  for i in range(tbin_stop):
    if if_ebl is False:
      model = datapath + 'run0406_ID000126_tbin%d.xml' % i  # !!
      event = simpath + f + "_tbin%02d.fits" % i
    else:
      model = datapath + 'run0406_ID000126_ebl_tbin%d.xml' % i  # !!
      event = simpath + f + "_ebl_tbin%02d.fits" % i
    event_bins.append(event)
    if not os.path.isfile(event):
      simulate_event(model=model, event=event, t=[t[i], t[1 + i]], e=[elow, ehigh], caldb=caldb, irf=irf,
                     pointing=[pointRA, pointDEC], seed=count)

      # observation list ---!
  eventList = simpath + 'obs_%s.xml' % f
  if not os.path.isfile(eventList):
    observation_list(event=event_bins, eventList=eventList, obsname=f)


