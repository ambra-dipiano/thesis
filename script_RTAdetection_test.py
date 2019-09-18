# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMpORTS ---!
import gammalib
import ctools
import cscripts
# from astropy.io import fits
# from module_plot import showSkymap
from module_analysis import *
from class_xml import *
import numpy as np
import csv
import os
import sys
from sys import argv
# from lxml import etree as ET
# import os.path

# ==============
# !!! SET Up !!!
# ==============

# initialize global count ---!
chunk = int(sys.argv[1])  # global count
trials = int(sys.argv[2])  # number of trials
count = int(sys.argv[3])  # starting count

# ctools/cscripts parameters ---!
caldb = 'prod3b'
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
nbin = int(roi / wbin)  # skymap x,y axes (pixel)
confidence = (0.68, 0.95, 0.9973)  # confidence interval for asymmetrical errors (%)

# pointing with off-axis equal to max prob GW ---!
offmax = (-1.475, -1.371)  # (deg)
trueRa = 33.057  # (deg)
trueDec = -51.841  # (deg)
pointRA = trueRa + offmax[0]  # (deg)
pointDEC = trueDec + offmax[1]  # (deg)

# others ---!
checks = False
if_ebl = True
if_fits = True
fileroot = 'run0406_'

# inputs ---!
cfg = loadConfig()
p = fsMng(cfg)
model = p.getWorkingDir() + 'run0406_ID000126.xml'
if if_ebl is True :
  template = p.getWorkingDir() + 'run0406_ID000126_ebl.fits'
else :
  template = p.getWorkingDir() + 'run0406_ID000126.fits'

# template ---!
if if_fits is True :
  nominal = p.getWorkingDir() + 'run0406_ID000126.fits'
  absorbed = p.getWorkingDir() + 'run0406_ID000126_ebl.fits'
  fiducial = '/home/ambra/Desktop/cluster-morgana/gilmore_tau_fiducial.csv'

  fits_ebl(nominal, absorbed, fiducial, zfetch=False, z='0.10', plot=False)

# =====================
# !!! LOAD TEMpLATE !!!
# =====================

t, tbin_stop = load_template(template, tmax, extract_spec=True, model=model, pathout=p.getDataDir())
print('!!! check ---- tbin_stop=', tbin_stop) if checks is True else None

for k in range(trials):
  count += 1
  # attach ID to fileroot ---!
  if if_ebl is True :
    f = fileroot + 'ebl%06d' % (count)
  else :
    f = fileroot + 'sim%06d' % (count)
  print('!!! check ---- file=', f) if checks is True else None

  # ====================
  # !!! SIMULATE GRB !!!
  # ====================

  event_bins = []

  # simulate ---!
  for i in range(tbin_stop):
    if if_ebl is False:
      model = p.getDataDir() + 'run0406_ID000126_tbin%d.xml' % i  # !!
      event = p.getSimDir() + f + "_tbin%02d.fits" % i
    else:
      model = p.getDataDir() + 'run0406_ID000126_ebl_tbin%d.xml' % i  # !!
      event = p.getSimDir() + f + "_ebl_tbin%02d.fits" % i
    event_bins.append(event)
    if not os.path.isfile(event):
      simulate_event(model=model, event=event, t=[t[i], t[1 + i]], e=[elow, ehigh], caldb=caldb, irf=irf,
                     pointing=[pointRA, pointDEC], seed=count)

      # observation list ---!
  eventList = p.getSimDir() + 'obs_%s.xml' % f
  if not os.path.isfile(eventList):
    observation_list(event=event_bins, eventList=eventList, obsname=f)

    # ============================
  # !!! EVENT LIST SELECTION !!!
  # ============================

  selectedEvents = []

  for i in range(tint):
    selectedEvents.append(eventList.replace(p.getSimDir(), p.getSelectDir()).replace('obs_', 'texp%ds_' % texp[i]))
    prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
    if not os.path.isfile(selectedEvents[i]):
      select_event(eventList=eventList, event_selected=selectedEvents[i], prefix=prefix, t=[tmin, tmax[i]],
                   e=[emin, emax])

  print('!!! check --- selection: ', selectedEvents) if checks is True else None

  # ========================
  # !!! SELECTION SKYMAp !!!
  # ========================

  skymapName = []

  for i in range(tint):
    skymapName.append(selectedEvents[i].replace(p.getSelectDir(), p.getDetDir()).replace('.xml', '_skymap.fits'))
    if not os.path.isfile(skymapName[i]):
      skymap_event(event_selected=selectedEvents[i], sky=skymapName[i], e=[emin, emax], caldb=caldb, irf=irf, wbin=wbin,
                   nbin=nbin)

  print('!!! check --- skymaps: ', skymapName) if checks is True else None

  # ==============================
  # !!! DETECTION AND MODELING !!!
  # ==============================

  detXml = []
  detReg = []
  pos = []

  for i in range(tint):
    # needs the analysis class as the prms init goes in there, not in the xml ---!
    detectionXml, detectionReg = runDetection(skymap=skymapName[i])
    detObj = xmlMng(detectionXml)
    detObj.sigma = sigma
    detObj.maxSrc = 1
    newXml = detObj.modXml()
    detXml.append(detectionXml)
    detReg.append(detectionReg)
    pos.append(detObj.loadTSV())

    print('\n\n==========\n\n!!! check --- detection.............', texp[i], 's done\n\n!!! coords:', pos,
          '\n\n ==========\n\n') if checks is True else None

  # =========================
  # !!! DETECTED RA & DEC !!!
  # =========================

  raDet = []
  decDet = []
  Ndet = []

  for i in range(tint):
    raDet.append(pos[i][0][0]) if len(pos[i][0]) > 0 else raDet.append(np.nan)
    decDet.append(pos[i][1][0]) if len(pos[i][1]) > 0 else decDet.append(np.nan)
    Ndet.append(len(pos[i][0]))
    print('!!! check number of detections in trial', k + 1, ' ====== ', Ndet) if checks is True else None

  raSrc001 = raDet
  decSrc001 = decDet
  print('!!! check ------- 1° src DETECTED RA:', raSrc001, ' and DEC:', decSrc001) if checks is True else None

  # ==================
  # !!! LIKELIHOOD !!!
  # ==================

  resultsName = []

  for i in range(tint):
    resultsName.append(detXml[i].replace('_det%dsgm.xml' % sigma, '_det%dsgm_results.xml' % sigma))
    if Ndet[i] > 0:
      if not os.path.isfile(resultsName[i]):
        max_likelihood(event_selected=selectedEvents[i], detection_model=detXml[i], results=resultsName[i], caldb=caldb,
                       irf=irf)
  print('!!! check --- max likelihoods: ', resultsName) if checks is True else None

  # ======================
  # !!! LIKELIHOOD TSV !!!
  # ======================

  tsList = []

  for i in range(tint):
    if Ndet[i] > 0:
      tsList.append(xmlMng.loadTSV(resultsName[i]))
    else:
      tsList.append([np.nan])

  print('!!! check --- TSV List for each texp: ', tsList) if checks is True else None

  # only first elem ---!
  tsv = []
  for i in range(len(tsList)):
    tsv.append(tsList[i][0])

  print('!!! check --- SRC001 TSV for each texp:', tsv) if checks is True else None

  # =================
  # !!! N SOURCES !!!
  # =================

  # count src with TS >= 9
  Nsrc = []
  for i in range(tint):
    n = 0
    for j in range(len(tsList[i])):
      if float(tsList[i][j]) >= 9:
        n += 1

    Nsrc.append(n)

    print('!!! check ---- SRC NUMBER with TS > 9 for trial ', k + 1, ' ===== ', Nsrc[i], ' for texp ==== ',
          texp[i]) if checks is True else None

  # ===========================
  # !!! ASYMMETRICAL ERRORS !!!
  # ===========================

  errorsName = []

  for i in range(tint):
    errorsName.append(resultsName[i].replace('_results.xml', '_errors.xml'))
    if Ndet[i] > 0:
      if not os.path.isfile(errorsName[i]):
        errors_conf = confidence_lv(event_selected=selectedEvents[i], results=resultsName[i], asym_errors=errorsName[i],
                                    caldb=caldb, irf=irf)
  print('!!! check --- asym errors: ', errors_conf) if checks is True else None

  # ===========================
  # !!! LIKELIHOOD RA & DEC !!!
  # ===========================

  raList = []
  decList = []
  prefErr = []
  for i in range(tint):
    if Ndet[i] > 0:
      raList.append(xmlMng.loadRaDec(resultsName[i])[0])
      decList.append(xmlMng.loadRaDec(resultsName[i])[1])
    else:
      raList.append([np.nan])
      decList.append([np.nan])

  raFit = []
  decFit = []
  raFit_err = []
  decFit_err = []
  for i in range(len(raList)):
    raFit.append(raList[i][0])
    decFit.append(decList[i][0])

  print('!!! check --- RA FIT for each texp: ', raList, '\n\nand RA SRC001 for each texp:',
        raFit) if checks is True else None
  print('!!! check --- DEC FIT for each texp: ', decList, '\n\nand DEC SRC001 for each texp:',
        decFit) if checks is True else None

  # ==================================
  # !!! LIKELIHOOD SpECTRAL pARAMS !!!
  # ==================================

  pref = []
  prefErr = []
  index = []
  indexErr = []
  pivot = []
  for i in range(tint):
    if Ndet[i] > 0:
      pref.append(xmlMng.loadSpectral(resultsName[i])[0][0])
      index.append(xmlMng.loadSpectral(resultsName[i])[0][1])
      pivot.append(xmlMng.loadSpectral(resultsName[i])[0][2])
      prefErr.append(xmlMng.loadSpectral(resultsName[i])[1][0])
      indexErr.append(xmlMng.loadSpectral(resultsName[i])[1][1])

    else:
      pref.append([np.nan])
      index.append([np.nan])
      pivot.append([np.nan])
      prefErr.append([np.nan])
      indexErr.append([np.nan])

  pref = []
  Index = []
  pivot = []
  pref_err = []
  Index_err = []
  for i in range(len(pref)):
    pref.append(pref[i][0])
    Index.append(index[i][0])
    pivot.append(pivot[i][0])
    pref_err.append(prefErr[i][0])
    Index_err.append(indexErr[i][0])

  print('!!! check ----- prefactor:', pref) if checks is True else None
  print('!!! check ----- index:', Index) if checks is True else None
  print('!!! check ----- pivot:', pivot) if checks is True else None

  # ====================
  # !!! COMpUTE FLUX !!!
  # ====================

  #  (flux, f_lerr, f_uerr) = ((intensity, lerr, uerr) * (energy ** 2)) / (6.4215 * 1e5)

  flux = []
  flux_ph = []
  flux_erg = []
  flux_en = []
  for i in range(tint):
    if Ndet[i] > 0:
      gamma = Index[i] + 1
      cost = pref[i] / gamma
      cost_erg = cost * 6.42 * 1e5
      integ = (emax * 1e6 / pivot[i]) ** gamma - (emin * 1e6 / pivot[i]) ** gamma
      integ_erg = (emax * 1e6)

      flux.append(cost * integ)  # convert E (MeV)
      print('!!! check ----- next is gammalib.plaw_photon_flux') if checks is True else None
      flux_ph.append(pref[i] * gammalib.plaw_photon_flux(emin, emax, pivot[i], Index[i]))  # takes E (TeV)
      flux_erg.append(cost_erg * integ_erg)  # convert E (erg)
      print('!!! check ----- next is gammalib.plaw_energy_flux') if checks is True else None
      flux_en.append(pref[i] * gammalib.plaw_energy_flux(emin, emax, pivot[i], Index[i]))  # takes E (TeV)
    else:
      flux.append(np.nan)
      flux_ph.append(np.nan)
      flux_erg.append(np.nan)
      flux_en.append(np.nan)

  # ==================
  # !!! CHECK FLUX !!!
  # ==================

  #  uplimName = []

  #  for i in range(tint) :
  #    if Ndet[i] > 0 :
  #      integrated_flux(selectedEvents[i], resultsName[i], caldb=caldb, irf=irf, conf=0.95, erange=[0.03, 0.5])
  #  print('!!! check --- max likelihoods: ', resultsName)

  # ================== #
  # !!! WRITE CSV FILE #
  # ================== #

  csvName = []
  header = '#trial,t exp,sigma,Ndet,Nsrc,RA Src001,DEC Src001,RA Fit,DEC Fit,TSV,flux1 ph,flux2 ph,flux1 erg,flux2 erg\n'
  ID = 'ID%06d' % count

  for i in range(tint):
    csvName.append(p.getCsvDir() + fileroot + '%ds_chunk%02d.csv' % (texp[i], chunk))

    row = []
    row.append([ID, texp[i], sigma, Ndet[i], Nsrc[i], raSrc001[i], decSrc001[i], raFit[i], decFit[i], tsv[i], flux[i],
                flux_ph[i], flux_erg[i], flux_en[i]])

    print('!!! check row --- iter', i, '=====', row) if checks is True else None

    if os.path.isfile(csvName[i]) == True:
      with open(csvName[i], 'a') as f:
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    else:
      with open(csvName[i], 'w+') as f:
        f.write(header)
        w = csv.writer(f)
        w.writerows(row)
        f.close()

  print('!!! check --- data file: ', csvName) if checks is True else None

  # ==================== #
  # !!! DELETE SIM FILES #
  # ==================== #

  print('!!! check --- ', count, ') trial done...') if checks is True else None

  if count > 4:
    os.system('rm ' + p.getSimDir() + '*sim%06d*' % count)
    os.system('rm ' + p.getSelectDir() + '*sim%06d*' % count)
    os.system('rm ' + p.getDetDir() + '*sim%06d*' % count)

print('!!! check end\n\ndone......chunk ', chunk, 'sim id from ', trials * (chunk - 1) + 1, ' to ',
      count) if checks is True else None
print(
  '!!! check end\n\ndone...... removed all files in sim, selected_sim and detection_all except seeds from 1 to 4') if checks is True else None

print('\n\n\n\n\n\n\ndone\n\n\n\n\n\n\n\n') if checks is True else None



