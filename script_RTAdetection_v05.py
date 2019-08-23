# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMPORTS ---!
import ctools
import gammalib
import cscripts
from astropy.io import fits
#from module_plot import showSkymap
from module_xml import *
import numpy as np
import csv
import os
import sys
from sys import argv


task = int(sys.argv[1])

caldb = 'prod3b'
irf = 'South_z40_average_100s'

# work with absolute paths ---!
workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/run0406/'
runpath = workdir + 'run0406_ID000126/'
simpath = runpath + 'sim/'
selectpath = runpath + 'selected_sim/'
datapath = runpath + 'data/'
csvpath = runpath + 'csv/'
detpath = runpath + 'detection_all/'
fileroot = 'run0406_'


# trials
trials = 10
# trials count
count = 10*(task-1)

# =====================
# !!! SET UP TRIALS !!!
# =====================

sigma = 5
texp = [1, 5, 10, 100]
texp.sort()
tint = len(texp)
tmin = 30
tmax = []
for i in range(tint):
  tmax.append(tmin + texp[i])

emin = 0.03  # 300 GeV
emax = 0.5  # 500 GeV
roi = 5


# ====================
# !!! SIMULATE GRB !!!
# ====================

# pointing with off-axis equal to max prob GW ---!
offmax = [31.582, -53.211]

# open template ---!
hdul = fits.open('%srun0406_ID000126.fits' % workdir)
# positions [deg] ---!
primary = hdul[0].header
trueRA = primary['RA']
trueDEC = primary['DEC']
# energybins [GeV] ---!
energy = np.array(hdul[1].data)
# timebins [s] ---!
time = np.array(hdul[2].data)
# spectra ---!
spectra = np.array(hdul[3].data)

Nt = len(time)
Ne = len(energy)
pointRA = offmax[0]
pointDEC = offmax[1]

# time grid ---!
t = [0.0 for x in range(Nt + 1)]
for i in range(Nt - 1):
  t[i + 1] = time[i][0] + (time[i + 1][0] - time[i][0]) / 2
# tmax in last bin ---!
t[Nt] = time[Nt - 1][0] + (time[Nt - 1][0] - t[Nt - 1])

tbin_stop = 0
for bin in range(len(t)) :
  if t[bin] <= 100 :
    tbin_stop += 1
  else :
    continue

print(tbin_stop)

# energy grid ---!
en = [1.0 for x in range(Ne + 1)]
for i in range(Ne - 1):
  en[i + 1] = energy[i][0] + (energy[i + 1][0] - energy[i][0]) / 2
# Emax in last bin ---!
en[Ne] = energy[Ne - 1][0] + (energy[Ne - 1][0] - en[Ne - 1])

# simulate N independent photon-lists of the same event ---!
for k in range(trials):
  count += k + 1
  # attach ID to fileroot
  f = fileroot + 'sim%06d' % (count)

  # photon lists for each bin ---!
  for i in range(tbin_stop):

    sim = ctools.ctobssim()
    sim["inmodel"] = datapath + 'run0406_ID000126_tbin%d.xml' % i
    sim["outevents"] = simpath + f + "_tbin%02d.fits" % i
    sim["caldb"] = "prod3b"
    sim["irf"] = "South_z40_average_100s"
    sim["ra"] = pointRA
    sim["dec"] = pointDEC
    sim["rad"] = 5.0
    sim["tmin"] = t[i]
    sim["tmax"] = t[i + 1]
    sim["emin"] = 0.03
    sim["emax"] = 1.0
    sim["seed"] = count
    sim["logfile"] = simpath + f + "_tbin%02d.log" % i
    sim.execute() if k > 0 else None

  # combine in observatiion list
  xml = gammalib.GXml()
  obslist = xml.append('observation_list title="observation library"')

  for i in range(tbin_stop):
    obs = obslist.append('observation name="run0406_sim%05d" id="%02d" instrument="CTA"' % (count,i))
    obs.append('parameter name="EventList" file="%s%s_tbin%02d.fits"' % (simpath, f, i))

  eventList = '%sobs_%s.xml' % (detpath, f)
  xml.save(eventList)
  print('!!! check --- eventList: ', eventList)

  # assign eventList to observation list ---!

  # =======================
  # !!! EVENT LIST INFO !!!
  # =======================

  # inspect observation definition files ---!
  # eventList = phList.replace(xmlpath, simpath)
  #pointingReg = eventList.replace('.xml', '.reg')
  #
  # # obs ---!
  #info = cscripts.csobsinfo()
  #info['inobs'] = eventList
  #info['outds9file'] = pointingReg
  #info['logfile'] = pointingReg.replace('.reg', '_pointing.log')
  #info['debug'] = bool('no')
  #info.execute()

  # ============================
  # !!! EVENT LIST SELECTION !!!
  # ============================

  #eventList = phList
  selectedEvents = []

  for i in range(tint):
    selectedEvents.append(eventList.replace('obs_', 'texp%ds_' % texp[i]))

    selection = ctools.ctselect()
    selection['inobs'] = eventList
    selection['outobs'] = selectedEvents[i]
    selection['usepnt'] = bool('yes')
    selection['prefix'] = selectpath + 'texp%ds_' % texp[i]
    selection['rad'] = roi
    selection['tmin'] = tmin
    selection['tmax'] = tmax[i]
    selection['emin'] = emin
    selection['emax'] = emax
    selection['logfile'] = selectedEvents[i].replace('.xml', '.log')
    selection['debug'] = bool('no')
    selection.execute() if k > 0 else None

  print('!!! check --- selection: ', selectedEvents)


  # ========================
  # !!! SELECTION SKYMAP !!!
  # ========================

  skymapName = []

  for i in range(tint):
    skymapName.append(selectedEvents[i].replace('.xml', '_skymap.fits'))

    skymap = ctools.ctskymap()
    skymap['inobs'] = selectedEvents[i]
    skymap['outmap'] = skymapName[i]
    skymap['irf'] = irf
    skymap['caldb'] = caldb
    skymap['emin'] = emin
    skymap['emax'] = emax
    skymap['usepnt'] = bool('yes')
    skymap['nxpix'] = 200
    skymap['nypix'] = 200
    skymap['binsz'] = 0.02
    skymap['coordsys'] = 'CEL'
    skymap['proj'] = 'CAR'
    skymap['bkgsubtract'] = 'IRF'
    skymap['logfile'] = skymapName[i].replace('.fits', '.log')
    skymap['debug'] = bool('no')
    skymap.execute() if k > 0 else None 

    #showSkymap(skymapName[i], show=False)

  print('!!! check --- skymaps: ', skymapName)

  # ==============================
  # !!! DETECTION AND MODELING !!!
  # ==============================

  detXml = []
  detReg = []
  pos = [[]]

  for i in range(tint):
    det, reg, coord = srcDetection_spcModeling(skymapName[i], sigma=sigma, maxSrc=10)
    detXml.append(det)
    detReg.append(reg)
    pos.append(coord)
    print('\n\n==========\n\n!!! check --- detection.............', texp[i],'s done\n\ncoords:', pos, '\n\n ==========\n\n')

  # print(detReg)
  # print(np.array(detXml).shape)
  # print(np.array(detReg).shape)
  # print(np.array(pos).shape)
  # print(pos)
  print('!!! check --- detection xml: ', detXml)

  # =========================
  # !!! DETECTED RA & DEC !!!
  # =========================

  raDet = []
  decDet = []
  raSrc001 = []
  decSrc001 = []
  Nsrc = []

  for i in range(tin):
    raDet.append(pos[i][0][0]) if len(pos[i][0]) > 0 else raDet.append(['NAN'])
    decDet.append(pos[i][1][0]) if len(pos[i][1]) > 0 else decDet.append(['NAN'])
    Nsrc.append(len(pos[i][0]))

  # flatten the lists ---!
  # raSrc001 = [float(el) for grp in raDet for el in grp]
  # decSrc001 = [float(el) for grp in decDet for el in grp]

  # already flattened lists ---!
  raSrc001 = raDet
  decSrc001 = decDet


  # =========================
  # !!! DETECTION SKYMAPS !!!
  # =========================

  # skymaps ---!
  # n = 0
  # for i in range(trials) :
  # 	if Nsrc[i] > 0 :
  # 	n += 1
  # 	showSkymap(skymapName[i], reg=detReg[i], suffix='detect%dsgm' % sigma, show=False)

  # print(n, 'detection skymaps')

  # ==================
  # !!! LIKELIHOOD !!!
  # ==================

  resultsName = []

  for i in range(tin):
    resultsName.append(detXml[i].replace('Mod', '_results'))

    if Nsrc[i] > 0:
      like = ctools.ctlike()
      like['inobs'] = selectedEvents[i]
      like['inmodel'] = detXml[i]
      like['outmodel'] = resultsName[i]
      like['caldb'] = caldb
      like['irf'] = irf
      like['fix_spat_for_ts'] = bool('yes')
      like['logfile'] = resultsName[i].replace('.xml', '.log')
      like['debug'] = bool('no')  # default
      like.execute()
      i += 1

  print('!!! check --- max likelihoos: ', resultsName)

  # ======================
  # !!! LIKELIHOOD TSV !!!
  # ======================

  tsList = []

  for i in range(tin):
    if Nsrc[i] > 0:
      tsList.append(getTSV(resultsName[i]))
    else:
      tsList.append(['NAN'])

  # flatten list ---!
  # tsv = [float(el) for grp in tsList for el in grp]

  # only first elem ---!
  tsv = []
  for i in range(len(tsList)):
    tsv.append(tsList[i][0])


  # ===========================
  # !!! LIKELIHOOD RA & DEC !!!
  # ===========================

  raList = []
  decList = []
  for i in range(tin):
    if Nsrc[i] > 0:
      raList.append(getRaDec(resultsName[i])[0])
      decList.append(getRaDec(resultsName[i])[1])
    else:
      raList.append(['NAN'])
      decList.append(['NAN'])

  # flatten lists ---!
  # raFit = [float(el) for grp in raList for el in grp]
  # decFit = [float(el) for grp in decList for el in grp]

  # already flattened lists ---!
  raFit = []
  decFit = []
  for i in range(len(raList)):
    raFit.append(raList[i][0])
    decFit.append(decList[i][0])


  # =================================
  # !!! LIKELIHOOD RA & DEC ERROR !!!
  # =================================

  raErr = []
  decErr = []
  for i in range(tin):
    if Nsrc[i] > 0:
      raErr.append(getRaDec_errors(resultsName[i])[0])
      decErr.append(getRaDec_errors(resultsName[i])[1])
    else:
      raErr.append(['NAN'])
      decErr.append(['NAN'])

  # flatten lists ---!
  # raFit_err = [float(el) for grp in raErr for el in grp]
  # decFit_err = [float(el) for grp in decErr for el in grp]

  # already flattened lists ---!
  raFit_err = []
  decFit_err = []
  mean_err = []

  for i in range(len(raErr)):
    raFit_err.append(raErr[i][0])
    decFit_err.append(decErr[i][0])
    mean_err.append(0.5 * (float(raErr[i][0]) + float(decErr[i][0])))


  # ================== #
  # !!! WRITE CSV FILE #
  # ================== #

  row = []
  csvName = []
  header = '#trial,t int,sigma,Nsrc,RA Src001,DEC Src001,RA Fit,DEC Fit,RA err,DEC err,MEAN err,TSV\n'

  for i in range(tin):
    csvName.append(csvpath + fileroot + '' + '_%03dsec.csv' % texp[i])
    row.append(['ID' + str(count), texp[i], sigma, Nsrc[i],
                raSrc001[i], decSrc001[i], raFit[i], decFit[i],
                raFit_err[i], decFit_err[i], mean_err[i], tsv[i]])


    open(csvName[i], 'w+').write(header).close()

    with open(csvName[i], 'w+') as f:
      w = csv.writer(f)
      w.writerows(row[i])

    f.close()

  print('!!! check --- data file: ', csvName)

  # ==================== #
  # !!! DELETE SIM FILES #
  # ==================== #

  #os.system('rm ' + simpath + '*')

  #print(count, 'trial done...')

print('done')
