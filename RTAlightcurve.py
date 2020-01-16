# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMPORTS ---!
from pkg_blindsearch import *
from module_plot import *
import numpy as np
import csv
import os
import sys

# --------------------------------- SETUP --------------------------------- !!!

# initialize global count ---!
# chunk = int(sys.argv[1])  # global count
# trials = int(sys.argv[2])  # number of trials
# count = int(sys.argv[3])  # starting count
# compact initialisation ---!
trials = 1  # trials
count = 0  # starting count
# cpus ---!
nthreads = 2
os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
os.environ['MKL_NUM_THREADS'] = str(nthreads)

# ctools/cscripts parameters ---!
caldb = 'prod3b-v2'
# caldb_degraded = caldb.replace('prod', 'degr')
irf = 'South_z40_0.5h'

sigma = 5  # detection acceptance (Gaussian)
texp = (10, 100)  # exposure times (s)
tdelay = 30  # slewing time (s)
tmax = []
for i in range(len(texp)):
  tmax.append(tdelay + texp[i])
ttotal = 1e4 #1e6  # maximum tobs (4h at least) simulation total time (s)
add_hours = 7200 #7200  # +2h observation time added after first none detection (s)
run_duration = 1200  # 20min observation run time for LST in RTA (s) ---!
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 10.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 10.0  # selection maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)
wbin = 0.02  # skymap bin width (deg)
corr_rad = 0.1  # Gaussian
max_src = 5  # max candidates
ts_threshold = 25  # TS threshold for reliable detection
reduce_flux = None  # flux will be devided by factor reduce_flux, if nominal then set to None

# conditions control ---!
checks1 = True  # prints info
checks2 = False  # prints more info
if_ebl = True  # uses the EBL absorbed template
if_cut = False  # adds a cut-off parameter to the source model
ebl_fits = False  # generate the EBL absorbed template
extract_spec = True  # generates spectral tables and obs definition models
irf_degrade = False  # use degraded irf
compute_degr = False  # compute irf degradation
src_sort = True  # sorts scandidates from highest TS to lowest
skip_exist = False  # skips the step if ID exists in csv (issue: if True than add+2h will start anew from last csv tbin)
debug = False  # prints logfiles on terminal
if_log = True  # saves logfiles

# path configuration ---!
cfg = xmlConfig()
p = ConfigureXml(cfg)
# files ---!
fileroot = 'run0406_'
ebl_table = p.getRootDir() + '/ebl_tables/gilmore_tau_fiducial.csv'
merge_map = 'run0406_MergerID000126_skymap.fits'
nominal_template = 'run0406_ID000126.fits'
ebl_template = 'run0406_ID000126_ebl.fits'
model_pl = 'run0406_ID000126.xml'
tcsv = 'time_slices.csv'

# pointing with off-axis equal to max prob GW ---!
true_coord = (33.057, -51.841)  # true position of source RA/DEC (deg)
offmax = (-1.475, -1.370)  # off-axis RA/DEC (deg)
pointing = (true_coord[0] + offmax[0], true_coord[1] + offmax[1])  # pointing direction RA/DEC (deg)
# true_coord, pointing, offmax = getPointing(None, p.getWorkingDir()+nominal_template)
# pointing with off-axis equal to max prob GW ---!
print('coords true:', true_coord, 'point', pointing, 'off', offmax) if checks2 else None

# recap and dof ---!
dof, m2, m1 = getDof()
print('!!! *** !!! dof = ', m2, ' - ', m1, ' = ', dof)
print('!!! *** !!! EBL ABSORPTION:', if_ebl)
print('!!! *** !!! MODEL CUTOFF:', if_cut)
print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! nominal caldb:', caldb)
print('!!! *** !!! irf:', irf)
print('!!! *** !!! TS SORT:', src_sort)
print('!!! *** !!! FLUX REDUCED factor:', reduce_flux)
print('!!! *** !!! sim energy range: [', elow, ', ', ehigh, '] (TeV)')
print('!!! *** !!! selection energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi:', roi, ' (deg)')
print('!!! *** !!! pointing:', pointing, ' (deg)')
print('!!! *** !!! blind detection confidence:', sigma, ' sigmas')
print('!!! *** !!! detection confidence ts threshold:', ts_threshold)
print('!!! *** !!! total observation time:', ttotal, ' s')
print('!!! *** !!! additional observation time:', add_hours, ' s')
print('!!! *** !!! delay time:', tdelay, ' s\n')

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup trials obj ---!
tObj = Analysis()
tObj.nthreads = nthreads
tObj.pointing = pointing
tObj.roi = roi
tObj.e = [elow, ehigh]
tObj.tmax = ttotal
tObj.model = p.getWorkingDir() + model_pl
tObj.debug = debug
tObj.if_log = if_log
# degrade IRF if required ---!
tObj.caldb = caldb
tObj.irf = irf
if irf_degrade:
  if compute_degr:
    tObj.degradeIrf()
  tObj.caldb = caldb.replace('prod', 'degr')
# add EBL to template ---!
if ebl_fits:
  tObj.template = p.getWorkingDir() + nominal_template  # nominal ---!
  new_template = p.getWorkingDir() + ebl_template  # absorbed ---!
  tObj.table = ebl_table  # fiducial ---!
  tObj.zfetch = True
  tObj.if_ebl = False
  tObj.fitsEbl(new_template)
# assign template ---!
template = p.getWorkingDir() + ebl_template
print('with EBL') if checks2 else None
tObj.if_ebl = if_ebl
tObj.template = template
print('!!! check ---- template=', tObj.template) if checks2 else None
# load template ---!
tObj.extract_spec = extract_spec
tbin_stop = tObj.loadTemplate()
print('!!! check ---- tbin_stop=', tbin_stop) if checks2 else None
print('!!! check ---- caldb:', tObj.caldb) if checks2 else None

# --------------------------------- REDUCE TEMPLATE FLUX  --------------------------------- !!!

if reduce_flux != None:
  tObj.factor = reduce_flux
  tObj.makeFainter()
  print('!!! check ---- reduce flux by factor %s' % str(reduce_flux)) if checks2 else None

# --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

for k in range(trials):
  count += 1
  tObj.seed = count
  clocking = tdelay-min(texp)  # simulate flowing time (subsequent temporal windows of 1s)
  GTIf = [run_duration for i in range(len(texp))]  # LST runs are of 20mins chunks ---!
  num = [1 for i in range(len(texp))]  # count on LST-like run chunks ---!
  print('\n\n!!! ************ STARTING TRIAL %d ************ !!!\n\n' % count) if checks1 else None
  print('!!! check ---- seed=', tObj.seed) if checks2 else None
  # attach ID to fileroot ---!
  f = fileroot + 'ebl%06d' % (count)
  if irf_degrade:
    f += 'irf'
  print('!!! check ---- file=', f) if checks2 else None

  # --------------------------------- SIMULATION --------------------------------- !!!

  event_bins = []
  tObj.table = p.getDataDir() + tcsv
  tgrid = tObj.getTimeSlices()  # methods which returns time slice edges
  # simulate ---!
  for i in range(tbin_stop):
    tObj.t = [tgrid[i], tgrid[i + 1]]
    if if_ebl:
      tObj.model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i
    else:
      tObj.model = p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i
    event = p.getSimDir() + f + "_tbin%02d.fits" % i
    if reduce_flux != None:
      tObj.model = tObj.model.replace('_tbin', '_flux%s_tbin' % str(reduce_flux))
      event = event.replace('_tbin', '_flux%s_tbin' % str(reduce_flux))
    event_bins.append(event)
    tObj.output = event
    if os.path.isfile(event):
      os.remove(event)
    tObj.eventSim()

  # --------------------------------- APPEND EVENTS IN PH-LIST --------------------------------- !!!

  # event_all = p.getSimDir() + f + '.fits'
  # if reduce_flux != None:
  #   event_all = event_all.replace('.fits', '_flux%s.fits' % str(reduce_flux))
  # tObj.input = event_bins
  # tObj.output = event_all
  # if run_duration == ttotal:
  #   tObj.appendEventsSinglePhList()
  #   phlist = event_all
  #   num_max = 1
  # else:
  #   num_max, phlist = tObj.appendEventsMultiPhList(max_length=run_duration, last=ttotal)
  #   print('phlist root name', phlist) if checks2 else None

  #breakpoint()

  # --------------------------------- 2° LOOP :: tbins --------------------------------- !!!

  tObj.e = [emin, emax]
  twindows = [int((ttotal-tdelay) / texp[i]) for i in range(len(texp))]  # number of temporal windows per exposure time in total time ---!
  tlast = [ttotal+tmax[i] for i in range(len(texp))]  # maximum observation time from last detection (not exceeding ttotal) ---!
  for i, t in enumerate(tlast):
    if t > ttotal:
      tlast[i] = ttotal
  is_detection = [True for i in range(len(texp))]  # controls which avoid forwarding of tlast for subsequent non-detections ---!
  # looping through all light-curve time intervals ---!
  for j in range(int(max(twindows))):
    clocking += min(texp)  # passing time second by second ---!
    print(clocking, 'clock', tlast, is_detection) if checks1 else None

    # --------------------------------- CLOCKING BREAK --------------------------------- !!!

    # check tlast, if globally reached then stop current trial ---!
    if clocking > max(tlast):
      print('end analysis trial', count, ' at clocking', tlast)
      break
    current_twindows = []

    # --------------------------------- 3° LOOP :: texp in tbin --------------------------------- !!!

    for i in range(len(texp)):
      if j == 0:
        current_twindows.append(texp[i])
      else:
        current_twindows.append(texp[i]) if clocking % texp[i] == 0 else None
    # looping for all the texp for which the tbin analysis needs to be computed ---!
    for i in range(len(current_twindows)):

      # --------------------------------- CLOCKING SKIP --------------------------------- !!!

      # check tlast, if locally reached then skip current bin ---!
      index = texp.index(current_twindows[i])
      if clocking > tlast[index]:
        print('skip analysis texp', texp[index]) if checks1 else None
        continue
      # --------------------------------- CHECK SKIP --------------------------------- !!!

      tbin = clocking / current_twindows[i]  # temporal bin number of this analysis
      IDbin = 'tbin%09d' % tbin
      csvName = p.getCsvDir() + fileroot + 'ID%06d_%ds.csv' % (count, texp[index])
      if irf_degrade:
        csvName = csvName.replace('.csv', '_degraded.csv')
      if os.path.isfile(csvName):
        skip = checkTrialId(csvName, IDbin)
      else:
        skip = False
      if skip_exist is True and skip is True:
        continue

      # --------------------------------- SELECTION TIME --------------------------------- !!!

      # if first tbin of tepx then don't add clocking time to selection edges ---!
      if clocking < tdelay:
        continue
      elif clocking == tdelay:
        tObj.t = [tdelay, tmax[index]]
      elif clocking > tdelay and texp[index] == min(texp):
        tObj.t = [clocking, texp[index]+clocking]
      elif clocking > tdelay and texp[index] != min(texp):
        tObj.t = [tdelay + clocking, tmax[index] + clocking]
      if tObj.t[1] > ttotal:
        tObj.t[1] = ttotal

      # --------------------------------- OBSERVATION LIST --------------------------------- !!!

      # using runs ---!
      # if run_duration != ttotal:
      #   if (tObj.t[0] < GTIf[index] or tObj.t[0] == GTIf[index]) and (tObj.t[1] == GTIf[index] or tObj.t[1] < GTIf[index]):
      #     event_bins = [phlist.replace('.fits', '_n%03d.fits' %num[index])]
      #   elif (tObj.t[0] < GTIf[index] or tObj.t[0] == GTIf[index]) and (tObj.t[1] == GTIf[index] or tObj.t[1] > GTIf[index]):
      #     event_bins = [phlist.replace('.fits', '_n%03d.fits' %num[index])]
      #     if num[index]+1 < num_max or num[index]+1 == num_max:
      #       event_bins.append(phlist.replace('.fits', '_n%03d.fits' %(num[index]+1)))
      #     GTIf[index] += run_duration
      #     num[index] += 1
      #   else:
      #     event_bins = [phlist.replace('.fits', '_n%03d.fits' %(num[index]+1))]
      #     GTIf[index] += run_duration
      #     num[index] += 1
      # else:
      #   event_bins = [phlist]

      # using tbins ---!
      phlist = p.getSelectDir() + f + '.fits'
      bins = tObj.getTimeBins(GTI=tObj.t, tgrid=tgrid)
      event_bins = []
      for bin in bins:
        event_bins.append(p.getSimDir() + f + '_tbin%02d.fits' %bin)

      print('event bins', event_bins) if checks2 else None
      # actual computation of obs list ---!
      event_list = phlist.replace('run0406_ebl', 'obs_t%dt%d_ebl' %(tObj.t[0], tObj.t[1])).replace('.fits', '.xml')
      if os.path.isfile(event_list):
        os.remove(event_list)
      tObj.input = event_bins
      tObj.output = event_list
      tObj.obsList(obsname='run0406_ID000126')
      print('event_bins', event_bins) if checks2 else None
      print('current phlist', phlist) if checks2 else None
      print('event_list', event_list) if checks2 else None

      # --------------------------------- SELECTION --------------------------------- !!!

      event_selected = event_list.replace(p.getSimDir(), p.getSelectDir()).replace('obs_', 'texp%ds_' %texp[i])
      prefix = p.getSelectDir() + 'texp%ds_' %texp[i]
      # select events ---!
      if os.path.isfile(event_selected):
        os.remove(event_selected)
      tObj.input = event_list
      tObj.output = event_selected
      tObj.eventSelect(prefix=prefix)
      print('selection', tObj.output) if checks2 else None

      # --------------------------------- SKYMAP --------------------------------- !!!

      skymap = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.xml', '_skymap.fits')
      if os.path.isfile(skymap):
        os.remove(skymap)
      tObj.input = event_selected
      tObj.output = skymap
      tObj.eventSkymap(wbin=wbin)
      print('skymap', tObj.output) if checks2 else None

      # --------------------------------- DETECTION & MODELING --------------------------------- !!!

      tObj.corr_rad = corr_rad
      tObj.max_src = max_src
      detectionXml = skymap.replace('_skymap.fits', '_det%dsgm.xml' % sigma)
      if os.path.isfile(detectionXml):
        os.remove(detectionXml)
      tObj.input = skymap
      tObj.output = detectionXml
      tObj.runDetection()
      #showSkymap(skymap=skymap, reg=detectionXml.replace('.xml', '.reg'), show=False)
      print('detection', tObj.output) if checks2 else None
      detObj = ManageXml(detectionXml)
      detObj.sigma = sigma
      detObj.if_cut = if_cut
      detObj.modXml()
      detObj.prmsFreeFix()

      # --------------------------------- CANDIDATES NUMBER --------------------------------- !!!

      pos = []
      pos.append(detObj.loadRaDec())
      Ndet = len(pos[0][0])

      # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

      likeXml = detectionXml.replace('_det%dsgm' % tObj.sigma, '_like%dsgm' % tObj.sigma)
      if os.path.isfile(likeXml):
        os.remove(likeXml)
      tObj.input = event_selected
      tObj.model = detectionXml
      tObj.output = likeXml
      tObj.maxLikelihood()
      print('likelihood', tObj.output) if checks2 else None
      likeObj = ManageXml(likeXml)
      if src_sort and Ndet > 0:
        highest_ts_src = likeObj.sortSrcTs()[0]
        print('!!! check ---- highest TS: ', highest_ts_src) if checks1 else None
      else:
        highest_ts_src = None

      # --------------------------------- DETECTION RA & DEC --------------------------------- !!!

      pos, ra_det, dec_det = ([] for n in range(3))
      pos.append(detObj.loadRaDec(highest=highest_ts_src))
      ra_det.append(pos[0][0][0]) if len(pos[0][0]) > 0 else ra_det.append(np.nan)
      dec_det.append(pos[0][1][0]) if len(pos[0][0]) > 0 else dec_det.append(np.nan)
      Ndet = len(pos[0][0])

      # --------------------------------- CLOSE DET XML --------------------------------- !!!

      detObj.closeXml()

      # --------------------------------- BEST FIT TSV --------------------------------- !!!

      ts_list, ts = ([] for n in range(2))
      ts_list.append(likeObj.loadTs()) if Ndet > 0 else ts_list.append([np.nan])

      # only first elem ---!
      ts.append(ts_list[0][0])
      print('ts:', ts[0]) if checks1 else None

      # --------------------------------- Nsrc FOR TSV THRESHOLD --------------------------------- !!!

      # count src with TS >= 9
      n = 0
      for j in range(len(ts_list[0])):
        if float(ts_list[0][j]) >= ts_threshold:
          n += 1

      Nsrc = n

      # --------------------------------- +2h FROM LAST DETECTION --------------------------------- !!!

      if (float(ts[0]) < ts_threshold or str(ts[0]) == 'nan') and is_detection[index]:
        is_detection[index] = False
        # add 2hrs of obs time ---!
        tlast[index] = tObj.t[1] + add_hours  # +2h ---!
        print('+2h tlast = ', tlast[index], ' with texp = ', texp[index], 'at clocking', clocking)
        # only 4hrs of simulation avialable, if tlast exceeds them then reset to ttotal ---!
        if tlast[index] > ttotal:
          tlast[index] = ttotal
          print('reset tlast = ', tlast[index], ' with texp = ', texp[index])

      # --------------------------------- BEST FIT RA & DEC --------------------------------- !!!

      ra_list, ra_fit, dec_list, dec_fit = ([] for n in range(4))
      coord = likeObj.loadRaDec() if Ndet > 0 else None
      ra_list.append(coord[0]) if Ndet > 0 else ra_list.append([np.nan])
      dec_list.append(coord[1]) if Ndet > 0 else dec_list.append([np.nan])

      # only first elem ---!
      ra_fit.append(ra_list[0][0])
      dec_fit.append(dec_list[0][0])

      # --------------------------------- BEST FIT SPECTRAL --------------------------------- !!!

      pref_list, pref, index_list, index, pivot_list, pivot = ([] for n in range(6))
      pref_err_list, pref_err = ([] for n in range(2))
      likeObj.if_cut = if_cut
      spectral = likeObj.loadSpectral()
      index_list.append(spectral[0]) if Ndet > 0 else index_list.append([np.nan])
      pref_list.append(spectral[1]) if Ndet > 0 else pref_list.append([np.nan])
      pivot_list.append(spectral[2]) if Ndet > 0 else pivot_list.append([np.nan])
      error = likeObj.loadPrefError()
      pref_err_list.append(error) if Ndet > 0 else pref_err_list.append([np.nan])

      # only first elem ---!
      index.append(index_list[0][0])
      pref.append(pref_list[0][0])
      pivot.append(pivot_list[0][0])
      pref_err.append(pref_err_list[0][0])

      # eventually cutoff ---!
      if if_cut:
        cutoff_list, cutoff = ([] for n in range(2))
        cutoff_list.append(spectral[3]) if Ndet > 0 else cutoff_list.append([np.nan])
        cutoff.append(cutoff_list[0][0])

      # --------------------------------- INTEGRATED FLUX --------------------------------- !!!

      flux_ph, flux_ph_err = ([] for n in range(2))
      norm_factor = 1
      if Ndet > 0:
        flux_ph.append(tObj.photonFluxPowerLaw(index[0], pref[0], pivot[0], norm_factor=norm_factor))  # E (MeV)
        flux_ph_err.append(tObj.photonFluxPowerLaw(index[0], pref_err[0], pivot[0], norm_factor=norm_factor))  # E (MeV)
      else:
        flux_ph.append(np.nan)
        flux_ph_err.append(np.nan)

      # MISSING THE CUT-OFF OPTION ---!!!

      # --------------------------------- CLOSE LIKE XML --------------------------------- !!!

      likeObj.closeXml()

      # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

      header = '#tbin,tinit,tend,Ndet,Nsrc,RA_det,DEC_det,RA_fit,DEC_fit,flux_ph,TS\n'
      ID = 'ID%06d' % count
      IDbin = 'tbin%09d' % tbin

      row = []
      if checks1:
        print('!!! ---------- check trial:', count)
        print('!!! ----- check texp:', texp[i], 's between: [', tObj.t[0], ', ', tObj.t[1], ' ] s')
        print('!!! *** check Ndet:', Ndet)
        print('!!! *** check Nsrc:', Nsrc)
        print('!!! *** check ra_det:', ra_det[0])
        print('!!! *** check dec_det:', dec_det[0])
        print('!!! *** check ra_fit:', ra_fit[0])
        print('!!! *** check dec_fit:', dec_fit[0])
        print('!!! *** check flux_ph:', flux_ph[0])
        print('!!! *** check ts:', ts[0])
        print('!!! *** ---------------------------')

      row.append([IDbin, tObj.t[0], tObj.t[1], Ndet, Nsrc, ra_det[0], dec_det[0], ra_fit[0], dec_fit[0],
                  flux_ph[0], ts[0]])
      if os.path.isfile(csvName):
        with open(csvName, 'a') as csv_file:
          w = csv.writer(csv_file)
          w.writerows(row)
          csv_file.close()
      else:
        with open(csvName, 'w+') as csv_file:
          csv_file.write(header)
          w = csv.writer(csv_file)
          w.writerows(row)
          csv_file.close()

    # --------------------------------- CLEAR SPACE --------------------------------- !!!

      os.system('rm ' + p.getSelectDir() + '*ebl%06d*' % count)
      os.system('rm ' + p.getDetDir() + '*ebl%06d*' % count)

  if int(count) != 1:
    os.system('rm ' + p.getSimDir() + '*ebl%06d*' % count)

print('\n\n!!! ================== END ================== !!!\n\n')


