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
corr_rad = 0.1  # Gaussian
confidence = (0.68, 0.95, 0.9973)  # confidence interval for asymmetrical errors (%)
max_src = 10  # max candidates
ts_threshold = 9  # TS threshold for reliable detection
reduce_flux = None  # flux will be devided by factor reduce_flux, if nominal then set to None ---!

# conditions control ---!
checks = True  # prints checks info ---!
if_ebl = False  # uses the EBL absorbed template ---!
if_cut = False  # adds a cut-off parameter to the source model ---!
ebl_fits = False  # generate the EBL absorbed template ---!
extract_spec = True  # generates spectral tables and obs definition models ---!
irf_degrade = False  # use degraded irf ---!
src_sort = False  # sorts scandidates from highest TS to lowest ---!
skip_exist = False  # if an output already exists it skips the step ---!
debug = False  # prints logfiles on terminal ---!
if_log = True  # saves logfiles ---!

# files ---!
fileroot = 'run0406_'
cfg_file = '/config.xml'
ebl_table = os.environ.get('MORGANA') + '/gilmore_tau_fiducial.csv'
merge_map = 'run0406_MergerID000126_skymap.fits'
nominal_template = 'run0406_ID000126.fits'
ebl_template = 'run0406_ID000126_ebl.fits'
model_pl = 'run0406_ID000126.xml'
tcsv = 'time_slices.csv'
cfg = xmlConfig(cfg_file)
p = ConfigureXml(cfg)

# pointing with off-axis equal to max prob GW ---!
true_coord = (33.057, -51.841)  # true position of source RA/DEC (deg)
offmax = (-1.475, -1.370)  # off-axis RA/DEC (deg)
pointing = (true_coord[0] + offmax[0], true_coord[1] + offmax[1])  # pointing direction RA/DEC (deg)
# true_coord, pointing, offmax = getPointing(None, p.getWorkingDir()+nominal_template)
# pointing with off-axis equal to max prob GW ---!
print(true_coord, pointing, offmax) if checks is True else None

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup trials obj ---!
tObj = Analysis(cfg_file)
tObj.pointing = pointing
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
  if count == 0:
    tObj.degradeIrf() # exts=1, factor=2 by default => halve only Aeff ---!
  # tObj.caldb = caldb_degraded
# add EBL to template ---!
if ebl_fits:
  tObj.template = p.getWorkingDir() + nominal_template # nominal ---!
  new_template = p.getWorkingDir() + ebl_template # absorbed ---!
  tObj.table = ebl_table # fiducial ---!
  tObj.zfetch = True
  tObj.if_ebl = False
  tObj.fitsEbl(new_template)
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
tbin_stop = tObj.loadTemplate()
print('!!! check ---- tbin_stop=', tbin_stop) if checks is True else None
print('!!! check ---- caldb:', tObj.caldb)

# --------------------------------- REDUCE TEMPLATE FLUX  --------------------------------- !!!

if reduce_flux != None:
  tObj.factor = reduce_flux
  tObj.makeFainter()
  print('!!! check ---- reduce flux by factor %s' %str(reduce_flux)) if checks is True else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
  count += 1
  tObj.seed = count
  print('\n\n!!! ************ STARTING TRIAL %d ************ !!!\n\n' %count)
  print('!!! check ---- seed=', tObj.seed) if checks is True else None
  # attach ID to fileroot ---!
  if if_ebl:
    f = fileroot + 'ebl%06d' % (count)
  else :
    f = fileroot + 'sim%06d' % (count)
  if irf_degrade:
    f += 'irf'
  print('!!! check ---- file=', f) if checks is True else None

  # --------------------------------- SIMULATION --------------------------------- !!!

  event_bins = []
  tObj.table = p.getDataDir() + tcsv
  tgrid = tObj.getTimeSlices()  # methods which returns time slice edges
  # simulate ---!
  for i in range(tbin_stop):
    tObj.t = [tgrid[i], tgrid[i+1]]
    print(tgrid[i], tgrid[i+1]) if checks is True else None
    if if_ebl:
      tObj.model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i
      event = p.getSimDir() + f + "_ebl_tbin%02d.fits" % i
      print('!!! check ---- simulation %d with EBL' %(i+1)) if checks is True else None
    else:
      tObj.model = p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i
      event = p.getSimDir() + f + "_tbin%02d.fits" % i
      print('!!! check ---- simulation %d without EBL' %(i+1)) if checks is True else None
    if reduce_flux != None:
      tObj.model = tObj.model.replace('_tbin', '_flux%s_tbin' %str(reduce_flux))
      event = event .replace('_tbin', '_flux%s_tbin' %str(reduce_flux))
    event_bins.append(event)
    if not skip_exist:
      if os.path.isfile(event):
        os.remove(event)
      tObj.output = event
      tObj.eventSim()
  # observation list ---!
  event = event_bins
  event_list = p.getSimDir() + 'obs_%s.xml' % f
  if reduce_flux != None:
    event_list = event_list.replace('obs_', 'obs_flux%s_' %str(reduce_flux))
  if not skip_exist:
    if os.path.isfile(event_list):
      os.remove(event_list)
    tObj.input = event
    tObj.output = event_list
    tObj.obsList(obsname=f)
    print('!!! check ---- obs list=', event_list) if checks is True else None

  # --------------------------------- 2° LOOP :: texp --------------------------------- !!!

  ra_det, dec_det, Ndet = ([[] for i in range(4)] for j in range(3))
  ts, Nsrc, ra_fit, dec_fit = ([[] for i in range(4)] for j in range(4))
  pref, index, pivot = ([[] for i in range(4)] for j in range(3))
  flux_ph, flux_en = ([[] for i in range(4)] for j in range(2))
  cutoff = [[] for i in range(4)] if if_cut is True else None

  # --------------------------------- SELECTION --------------------------------- !!!

  tObj.e = [emin, emax]
  for i in range(tint):
    print('\n\n!!! ************ STARTING TEXP %d ************ !!!\n\n' % texp[i])
    tObj.t = [tmin, tmax[i]]
    event_selected = event_list.replace(p.getSimDir(), p.getSelectDir()).replace('obs_', 'texp%ds_' % texp[i])
    prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
    if not skip_exist:
      if os.path.isfile(event_selected):
        os.remove(event_selected)
      tObj.input = event_list
      tObj.output = event_selected
      tObj.eventSelect(prefix=prefix)
      print('!!! check ---- selection: ', event_selected) if checks is True else None

    # --------------------------------- SKYMAP --------------------------------- !!!

    skymap = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.xml', '_skymap.fits')
    if not skip_exist:
      if os.path.isfile(skymap):
        os.remove(skymap)
      tObj.input = event_selected
      tObj.output = skymap
      tObj.eventSkymap(wbin=wbin)
    print('!!! check ---- skymaps: ', skymap) if checks is True else None

    # --------------------------------- DETECTION & MODELING --------------------------------- !!!

    tObj.corr_rad = corr_rad
    tObj.max_src = max_src
    detectionXml = skymap.replace('_skymap.fits', '_det%dsgm.xml' %sigma)
    if not skip_exist:
      if os.path.isfile(detectionXml):
        os.remove(detectionXml)
      tObj.input = skymap
      tObj.output = detectionXml
      tObj.runDetection()
    detObj = ManageXml(detectionXml, cfg_file)
    detObj.sigma = sigma
    detObj.if_cut = if_cut
    detObj.modXml()
    print('!!! check ---- detection.............', texp[i], 's done') if checks is True else None
    print('\n\n!!! check ---- det mod: ', detectionXml) if checks is True else None

    # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

    #detObj.prmsFreeFix()
    likeXml = detectionXml.replace('_det%dsgm' % tObj.sigma, '_like%dsgm' % tObj.sigma)
    if not skip_exist:
      if os.path.isfile(likeXml):
        os.remove(likeXml)
      tObj.input = event_selected
      tObj.model = detectionXml
      tObj.output = likeXml
      tObj.maxLikelihood()
    likeObj = ManageXml(likeXml, cfg_file)
    if src_sort:
      highest_ts_src = likeObj.sortSrcTs()[0]
    else:
      highest_ts_src = None
      print('!!! check ---- highest TS: ', highest_ts_src) if checks is True else None
    print('\n\n!!! check ---- max likelihoods: ', likeXml) if checks is True else None

    # --------------------------------- DETECTION RA & DEC --------------------------------- !!!

    pos = []
    pos.append(detObj.loadRaDec(highest=highest_ts_src))
    print('!!! check ---- coords:', pos[0]) if checks is True else None
    ra_det[i].append(pos[0][0][0]) if len(pos[0][0]) > 0 else ra_det[i].append(np.nan)
    dec_det[i].append(pos[0][1][0]) if len(pos[0][0]) > 0 else dec_det[i].append(np.nan)
    Ndet[i].append(len(pos[0][0]))

    # --------------------------------- CLOSE DET XML --------------------------------- !!!

    detObj.closeXml()

    # --------------------------------- BEST FIT TSV --------------------------------- !!!

    ts_list = []
    ts_list.append(likeObj.loadTs()) if Ndet[i][0] > 0 else ts_list.append([np.nan])

    # only first elem ---!
    ts[i].append(ts_list[0][0])

    # --------------------------------- Nsrc FOR TSV THRESHOLD --------------------------------- !!!

    # count src with TS >= 9
    n = 0
    for j in range(len(ts_list[0])):
      if float(ts_list[0][j]) >= ts_threshold:
        n += 1

    Nsrc[i].append(n)

    # --------------------------------- BEST FIT RA & DEC --------------------------------- !!!

    ra_list = []
    dec_list = []
    coord = likeObj.loadRaDec() if Ndet[i][0] > 0 else None
    ra_list.append(coord[0]) if Ndet[i][0] > 0 else ra_list.append([np.nan])
    dec_list.append(coord[1]) if Ndet[i][0] > 0 else dec_list.append([np.nan])

    # only first elem ---!
    ra_fit[i].append(ra_list[0][0])
    dec_fit[i].append(dec_list[0][0])

    # --------------------------------- BEST FIT SPECTRAL --------------------------------- !!!

    pref_list = []
    index_list = []
    pivot_list = []
    cutoff_list = [] if if_cut is True else None
    if if_cut:
      likeObj.if_cut = if_cut
    spectral = likeObj.loadSpectral()
    index_list.append(spectral[0]) if Ndet[i][0] > 0 else index_list.append([np.nan])
    pref_list.append(spectral[1]) if Ndet[i][0] > 0 else pref_list.append([np.nan])
    pivot_list.append(spectral[2]) if Ndet[i][0] > 0 else pivot_list.append([np.nan])
    if if_cut:
      cutoff_list.append(spectral[3]) if Ndet[i][0] > 0 else cutoff_list.append([np.nan])

    # only first elem ---!
    index[i].append(index_list[0][0])
    pref[i].append(pref_list[0][0])
    pivot[i].append(pivot_list[0][0])
    if if_cut:
      cutoff[i].append(cutoff_list[0][0])

    # --------------------------------- INTEGRATED FLUX --------------------------------- !!!

    if Ndet[i][0] > 0:
      flux_ph[i].append(tObj.photonFluxPowerLaw(index[i][0], pref[i][0], pivot[i][0]))  # E (MeV)
      flux_en[i].append(tObj.energyFluxPowerLaw(index[i][0], pref[i][0], pivot[i][0]))  # E (erg)
    else:
      flux_ph[i].append(np.nan)
      flux_en[i].append(np.nan)

    # MISSING THE CUT-OFF OPTION ---!!!

    # --------------------------------- CLOSE LIKE XML --------------------------------- !!!

    if Ndet[i][0] > 0:
      likeObj.closeXml()

    # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

    header = '#trial,texp,sigma,Ndet,Nsrc,RA_det,DEC_det,RA_fit,DEC_fit,flux_ph,flux_erg,TS\n'
    ID = 'ID%06d' % count
    csvName = p.getCsvDir() + fileroot + '%ds_chunk%02d.csv' % (texp[i], chunk)

    row = []
    print('\n\n!!! ---------- check trial:', count)
    print('!!! ----- check texp:', texp[i])
    print('!!! *** check Ndet:', Ndet[i][0])
    print('!!! *** check Nsrc:', Nsrc[i][0])
    print('!!! *** check ra_det:', ra_det[i][0])
    print('!!! *** check dec_det:', dec_det[i][0])
    print('!!! *** check ra_fit:', ra_fit[i][0])
    print('!!! *** check dec_fit:', dec_fit[i][0])
    print('!!! *** check flux_ph:', flux_ph[i][0])
    # print('!!! *** check flux_en:', flux_en[i][0])
    print('!!! *** check ts:', ts[i][0])

    row.append([ID, texp[i], sigma, Ndet[i][0], Nsrc[i][0], ra_det[i][0], dec_det[i][0], ra_fit[i][0], dec_fit[i][0],
                flux_ph[i][0], flux_en[i][0], ts[i][0]])
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
  if count > 1:
    os.system('rm ' + p.getSimDir() + '*run*%06d*' % count)
    os.system('rm ' + p.getSelectDir() + '*run*%06d*' % count)
    os.system('rm ' + p.getDetDir() + '*run*%06d*' % count)

print('\n\n\n\n\n\n\ndone\n\n\n\n\n\n\n\n')



