# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMpORTS ---!
from pkg_blindsearch import *
import numpy as np
import csv
import os
import sys

# --------------------------------- SETUP --------------------------------- !!!

# initialize global count ---!
chunk = int(sys.argv[1])  # global count
trials = int(sys.argv[2])  # number of trials
count = int(sys.argv[3])  # starting count

# ctools/cscripts parameters ---!
caldb = 'prod3b'
# caldb_degraded = caldb.replace('prod', 'degr')
irf = 'South_z40_average_100s'

texp = [1, 5, 10, 100]  # exposure times (s)
texp.sort()
tint = len(texp)
tmin = 0  # slewing time (s)
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 0.5  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 0.5  # selection maximum energy (TeV)
roi = 5  # region of interest (deg)

# conditions control ---!
checks = True  # prints checks info ---!
irf_degrade = False  # use degraded irf ---!
skip_exist = False  # if an output already exists it skips the step ---!
debug = False  # prints logfiles on terminal ---!
if_log = True  # saves logfiles ---!

print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! sim energy range: [', elow, ', ', ehigh, '] (TeV)')
print('!!! *** !!! selection energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi: ', roi, ' (deg)')

# files ---!
model_pl = 'grb.xml'
model_bkg = 'CTAIrfBackground.xml'
tcsv = 'time_slices.csv'
cfg = xmlConfig()
p = ConfigureXml(cfg)

# pointing with off-axis equal to max prob GW ---!
true_coord = (33.057, -51.841)  # true position of source RA/DEC (deg)
offmax = (-1.475, -1.370)  # off-axis RA/DEC (deg)
pointing = (true_coord[0] + offmax[0], true_coord[1] + offmax[1])  # pointing direction RA/DEC (deg)

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup trials obj ---!
tObj = Analysis()
tObj.pointing = pointing
tObj.roi = roi
tObj.e = [elow, ehigh]
tObj.tmax = max(tmax)
tObj.caldb = caldb
tObj.irf = irf
tObj.debug = debug
tObj.if_log = if_log
# degrade IRF if required ---!
if irf_degrade:
  if count == 0:
    tObj.degradeIrf()
  # tObj.caldb = caldb_degraded
print('!!! check ---- caldb:', tObj.caldb) if checks is True else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
  count += 1
  tObj.seed = count
  print('\n\n!!! ************ STARTING TRIAL %d ************ !!!' %count) if checks is True else None
  print('!!! check ---- seed=', tObj.seed) if checks is True else None

  # --------------------------------- SIMULATION --------------------------------- !!!

  # attach ID to fileroot ---!
  f = 'bkg%06d' % (count)
  if irf_degrade:
    f += 'irf'
  # simulate ---!
  model = os.path.dirname(__file__) + model_bkg
  tObj.model = model
  event = p.getSimDir() + f + ".fits"
  if not skip_exist:
    if os.path.isfile(event):
      os.remove(event)
    tObj.output = event
    tObj.eventSim()
    print('!!! check ---- simulation=', event) if checks is True else None

  # --------------------------------- 2° LOOP :: texp --------------------------------- !!!

  # --------------------------------- SELECTION --------------------------------- !!!

  tObj.e = [emin, emax]
  for i in range(tint):
    print('\n\n!!! ************ STARTING TEXP %d ************ !!!\n\n' % texp[i]) if checks is True else None
    tObj.t = [tmin, tmax[i]]
    event_selected = event.replace(p.getSimDir(), p.getSelectDir()).replace('bkg', 'texp%ds_bkg' %texp[i])
    prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
    if not skip_exist:
      if os.path.isfile(event_selected):
        os.remove(event_selected)
      tObj.input = event
      tObj.output = event_selected
      tObj.eventSelect(prefix=prefix)
      print('!!! check ---- selection: ', event_selected) if checks is True else None

    # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

    model = os.path.dirname(__file__) + model_pl
    mObj = ManageXml(model)
    mObj.prmsFreeFix()
    mObj.closeXml()
    likeXml = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.fits', '_like.xml')
    if not skip_exist:
      if os.path.isfile(likeXml):
        os.remove(likeXml)
      tObj.input = event
      tObj.model = model
      tObj.output = likeXml
      tObj.maxLikelihood()
    likeObj = ManageXml(likeXml)
    print('!!! check ---- max likelihood: ', likeXml) if checks is True else None

    # --------------------------------- BEST FIT TSV --------------------------------- !!!

    ts_list, ts = ([] for j in range(2))
    ts_list.append(likeObj.loadTs())

    # only first elem ---!
    ts.append(ts_list[0])

    # --------------------------------- CLOSE LIKE XML --------------------------------- !!!

    likeObj.closeXml()

    # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

    header = '#trial,texp,TS\n'
    ID = 'ID%06d' % count
    csvName = p.getCsvDir() + 'bkg_%ds_chunk%02d.csv' % (texp[i], chunk)

    row = []
    print('\n\n!!! ---------- check trial:', count) if checks is True else None
    print('!!! ----- check texp:', texp[i]) if checks is True else None
    print('!!! *** check ts:', ts[0][0]) if checks is True else None

    row.append([ID, texp[i], ts[0][0]])
    print('!!! check row: seed %d --- texp' %i, texp[i], 's =====', row) if checks is True else None
    if os.path.isfile(csvName):
      with open(csvName, 'a') as f:
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    else:
      with open(csvName, 'w+') as f:
        f.write(header)
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    print('!!! check ---- data file: ', csvName) if checks is True else None

  # --------------------------------- CLEAR SPACE --------------------------------- !!!

  print('!!! check ---- ', count, ') trial done...') if checks is True else None
  if count > 4:
    os.system('rm ' + p.getSimDir() + '*bkg%06d*' % count)
    os.system('rm ' + p.getSelectDir() + '*bkg%06d*' % count)
    os.system('rm ' + p.getDetDir() + '*bkg%06d*' % count)

print('\n\n\n\n\n\n\ndone\n\n\n\n\n\n\n\n')



