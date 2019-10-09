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
caldb_degraded = caldb.replace('prod', 'degr')
irf = 'South_z40_average_100s'

sigma = 5  # detection acceptance (Gaussian)
texp = [1, 5, 10, 100]  # exposure times (s)
texp.sort()
tint = len(texp)
tmin = 30  # slewing time (s)
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 1.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 0.5  # selection maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)
wbin = 0.02  # skymap bin width (deg)
confidence = (0.68, 0.95, 0.9973)  # confidence interval for asymmetrical errors (%)
maxSrc = 1  # max candidates
corr_rad = 0.1  # Gaussian

# pointing with off-axis equal to max prob GW ---!
offmax = (-1.475, -1.371)  # (deg)
trueRa = 33.057  # (deg)
trueDec = -51.841  # (deg)
pointRA = trueRa + offmax[0]  # (deg)
pointDEC = trueDec + offmax[1]  # (deg)

# conditions control ---!
checks = True
if_cut = False
if_ebl = False
extract_spec = True
irf_degrade = False
src_sort = False
skip_exist = False
ebl_fits = False
debug = False
if_log = False

# files ---!
fileroot = 'run0406_'
cfg_file = '/config.xml'
ebl_table = os.environ.get('MORGANA') + '/gilmore_tau_fiducial.csv'
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
tObj.debug = debug
tObj.if_log = if_log
# degrade IRF if required ---!
if irf_degrade:
  tObj.degradeIRF()
# add EBL to template ---!
if ebl_fits:
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
else :
  template = p.getWorkingDir() + nominal_template
tObj.if_ebl = if_ebl
tObj.template = template
print('!!! check ---- template=', tObj.template) if checks is True else None
# load template ---!
tObj.extract_spec = extract_spec
tbin_stop = tObj.load_template()
print('!!! check ---- tbin_stop=', tbin_stop) if checks is True else None
print('!!! check ---- caldb:', tObj.caldb)
# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
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
    tObj.t = [time[i], time[i+1]]
    if if_ebl:
      tObj.model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i
      tObj.event = p.getSimDir() + f + "_ebl_tbin%02d.fits" % i
      print('!!! check ---- simulation %d with EBL' %(i+1)) if checks is True else None
    else:
      tObj.model = p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i
      tObj.event = p.getSimDir() + f + "_tbin%02d.fits" % i
      print('!!! check ---- simulation %d without EBL' %(i+1)) if checks is True else None
    event_bins.append(tObj.event)
    if not skip_exist:
      if os.path.isfile(tObj.event):
        os.remove(tObj.event)
      tObj.eventSim()
  print('!!! check ---- simulation=', tObj.event) if checks is True else None
  # observation list ---!
  tObj.event = event_bins
  tObj.event_list = p.getSimDir() + 'obs_%s.xml' % f
  if not skip_exist:
    if os.path.isfile(tObj.event_list):
      os.remove(tObj.event_list)
    tObj.obsList(obsname=f)
  print('!!! check ---- obs list=', tObj.event_list) if checks is True else None

  # --------------------------------- 2° LOOP :: texp --------------------------------- !!!

  raDet, decDet, Ndet = ([[] for i in range(4)] for j in range(3))
  ts, Nsrc, raFit, decFit = ([[] for i in range(4)] for j in range(4))
  Pref, Index, Pivot = ([[] for i in range(4)] for j in range(3))
  flux_ph, flux_en = ([[] for i in range(4)] for j in range(2))
  Cutoff = [[] for i in range(4)] if if_cut is True else None

  # --------------------------------- SELECTION --------------------------------- !!!

  tObj.e = [emin, emax]
  for i in range(tint):
    tObj.t = [tmin, tmax[i]]
    tObj.event_selected = tObj.event_list.replace(p.getSimDir(), p.getSelectDir()).replace('obs_', 'texp%ds_' % texp[i])
    prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
    if not skip_exist:
      if os.path.isfile(tObj.event_selected):
        os.system(tObj.event_selected)
      tObj.eventSelect(prefix=prefix)
    print('!!! check ---- selection: ', tObj.event_selected) if checks is True else None

  # --------------------------------- SKYMAP --------------------------------- !!!

    tObj.skymap = tObj.event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.xml', '_skymap.fits')
    if not skip_exist:
      if os.path.isfile(tObj.skymap):
        os.remove(tObj.skymap)
      tObj.eventSkymap(wbin=wbin)
    print('!!! check ---- skymaps: ', tObj.skymap) if checks is True else None

  # --------------------------------- DETECTION & MODELING --------------------------------- !!!

    tObj.corr_rad = corr_rad
    tObj.maxSrc = maxSrc
    tObj.detectionXml = tObj.skymap.replace('_skymap.fits', '_det%dsgm.xml' %sigma)
    if not skip_exist:
      if os.path.isfile(tObj.detectionXml):
        os.remove(tObj.detectionXml)
      tObj.runDetection()
    detObj = xmlMng(tObj.detectionXml, cfg_file)
    detObj.sigma = sigma
    detObj.if_cut = if_cut
    detObj.modXml()
    print('!!! check ---- detection.............', texp[i], 's done') if checks is True else None

  # --------------------------------- DETECTION RA & DEC --------------------------------- !!!

    pos = []
    pos.append(detObj.loadRaDec())
    print('!!! check ---- coords:', pos[0]) if checks is True else None
    raDet[i].append(pos[0][0][0]) if len(pos[0][0]) > 0 else raDet[i].append(np.nan)
    decDet[i].append(pos[0][1][0]) if len(pos[0][0]) > 0 else decDet[i].append(np.nan)
    Ndet[i].append(len(pos[0][0]))
    print('!!! check ---- number of detections in trial', k + 1, ' ====== ', Ndet[i][0]) if checks is True else None
    print('!!! check ---- 1° src DETECTED RA:', raDet[i][0], ' and DEC:', decDet[i][0]) if checks is True else None

  # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

    tObj.likeXml = tObj.detectionXml.replace('_det%dsgm.xml' % tObj.sigma, '_like%dsgm.xml' % tObj.sigma)
    if Ndet[i][0] > 0:
      if not skip_exist:
        if os.path.isfile(tObj.likeXml):
          os.remove(tObj.likeXml)
        tObj.maxLikelihood()
      likeObj = xmlMng(tObj.likeXml, cfg_file)
      if src_sort:
        likeObj.sortSrcTS()
    print('!!! check ---- max likelihoods: ', tObj.likeXml) if checks is True else None

  # --------------------------------- BEST FIT TSV --------------------------------- !!!

    tsList = []
    if Ndet[i][0] > 0:
      tsList.append(likeObj.loadTSV())
    else:
      tsList.append([np.nan])
    print('!!! check ---- TSV List for all sources: ', tsList[0]) if checks is True else None

    # only first elem ---!
    ts[i].append(tsList[0][0])
    print('!!! check ---- SRC001 TSV for candidate:', ts[i][0]) if checks is True else None

  # --------------------------------- Nsrc FOR TSV THRESHOLD --------------------------------- !!!

    # count src with TS >= 9
    n = 0
    for j in range(len(tsList[0])):
      if float(tsList[0][j]) >= 9:
        n += 1

    Nsrc[i].append(n)
    print('!!! check ---- SRC NUMBER with TS > 9 for trial ', k + 1, ' ===== ', Nsrc[i][0], ' for texp ==== ',
          texp[i], 's') if checks is True else None

  # --------------------------------- BEST FIT RA & DEC --------------------------------- !!!

    raList = []
    decList = []
    if Ndet[i][0] > 0:
      coord = likeObj.loadRaDec()
      raList.append(coord[0])
      decList.append(coord[1])
    else:
      raList.append([np.nan])
      decList.append([np.nan])

    raFit[i].append(raList[0][0])
    decFit[i].append(decList[0][0])

    print('!!! check --- RA FIT for all sources: ', raList[0], '\n!!! check --- RA FIT for candidate:',
          raFit[i][0]) if checks is True else None
    print('!!! check --- DEC FIT for all sources: ', decList[0], '\n!!! check --- DEC FIT for candidate:',
          decFit[i][0]) if checks is True else None

  # --------------------------------- BEST FIT SPECTRAL --------------------------------- !!!

    pref = []
    index = []
    pivot = []
    cutoff = [] if if_cut is True else None
    if Ndet[i][0] > 0:
      if if_cut:
        likeObj.if_cut = if_cut
      spectral = likeObj.loadSpectral()
      index.append(spectral[0])
      pref.append(spectral[1])
      pivot.append(spectral[2])
      if if_cut:
        cutoff.append(spectral[3])
    else:
      pref.append([np.nan])
      index.append([np.nan])
      pivot.append([np.nan])
      if if_cut:
        cutoff.append([np.nan])

    Index[i].append(index[0][0])
    Pref[i].append(pref[0][0])
    Pivot[i].append(pivot[0][0])
    if if_cut:
      Cutoff[i].append(cutoff[0][0])
    print('!!! check ----- index:', Index[i][0]) if checks is True else None
    print('!!! check ----- prefactor:', Pref[i][0]) if checks is True else None
    print('!!! check ----- pivot:', Pivot[i][0]) if checks is True else None
    if if_cut:
      print('!!! check ----- cutoff:', Cutoff[i]) if checks is True else None

  # --------------------------------- INTEGRATED FLUX --------------------------------- !!!

    if Ndet[i][0] > 0:
      flux_ph[i].append(tObj.photonFlux_pl(Index[i][0], Pref[i][0], Pivot[i][0]))  # E (MeV)
      flux_en[i].append(tObj.energyFlux_pl(Index[i][0], Pref[i][0], Pivot[i][0]))  # E (erg)
    else:
      flux_ph[i].append(np.nan)
      flux_en[i].append(np.nan)

    print('!!! check ----- my flux [ph/cm2/s]:', flux_ph[i][0]) if checks is True else None
    print('!!! check ----- my flux [erg/cm2/s]:', flux_en[i][0], '(*)') if checks is True else None

    # MISSING THE CUT-OFF OPTION ---!!!

  # --------------------------------- CLOSE XMLs --------------------------------- !!!

    detObj.closeXml()
    if Ndet[i][0] > 0:
      likeObj.closeXml()

  # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

  header = '#trial,t exp,sigma,Ndet,Nsrc,RA Src001,DEC Src001,RA Fit,DEC Fit,flux ph,flux erg,TSV\n'
  ID = 'ID%06d' % count
  for i in range(tint):
    csvName = p.getCsvDir() + fileroot + '%ds_chunk%02d.csv' % (texp[i], chunk)

    row = []
    print('!!! *** check Ndet:', Ndet[0])
    print('!!! *** check Nsrc:', Nsrc[0])
    print('!!! *** check raDet:', raDet[0])
    print('!!! *** check decDet:', decDet[0])
    print('!!! *** check raFit:', raFit[0])
    print('!!! *** check decFit:', decFit[0])
    print('!!! *** check flux_ph:', flux_ph[0])
    print('!!! *** check flux_en:', flux_en[0])
    print('!!! *** check ts:', ts[0])

    row.append([ID, texp[i], sigma, Ndet[i][0], Nsrc[i][0], raDet[i][0], decDet[i][0], raFit[i][0], decFit[i][0], flux_ph[i][0], flux_en[i][0], ts[i][0]])
    print('!!! check row --- iter', i, '=====', row) if checks is True else None
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
  if count > 1:
    os.system('rm ' + p.getSimDir() + '*run0406*%06d*' % count)
    os.system('rm ' + p.getSelectDir() + '*run0406*%06d*' % count)
    os.system('rm ' + p.getDetDir() + '*run0406*%06d*' % count)
print('!!! check end\n\ndone......chunk ', chunk, 'sim id from ', trials * (chunk - 1) + 1, ' to ',
      count) if checks is True else None
print('!!! check end\n\ndone...... removed all files in sim, selected_sim and detection_all except seeds from 1 to 4') if checks is True else None

print('\n\n\n\n\n\n\ndone\n\n\n\n\n\n\n\n')



