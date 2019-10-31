# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMPORTS ---!
import gammalib
import ctools
import cscripts
from astropy.io import fits
#from module_plot import showSkymap
from module_xml_old import *
import numpy as np
import csv
import os
import sys
from sys import argv
from lxml import etree as ET
import os.path



# =========================== #
# !!! PROPER CODE STARTS HERE #
# =========================== #

# initialize global count ---!
chunk = int(sys.argv[1]) # global count
trials = int(sys.argv[2]) # number of trials
count = int(sys.argv[3]) # starting count  

caldb = 'prod3b'
irf = 'South_z40_average_100s'

# work with absolute paths ---!
workdir = '/home/ambra/Desktop/cluster-morgana/run0406_test/run0406/'
runpath = workdir + 'run0406_ID000126/'
simpath = runpath + 'sim/'
selectpath = runpath + 'selected_sim/'
datapath = runpath + 'data/'
csvpath = runpath + 'csv/'
detpath = runpath + 'detection_all/'
fileroot = 'run0406_'

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
offmax = [-1.475, -1.371]
trueRa = 33.057
trueDec = -51.841
pointRA = trueRa + offmax[0]
pointDEC = trueDec + offmax[1]

# open template ---!
hdul = fits.open('%srun0406_ID000126.fits' % workdir)
# energybins [GeV] ---!
energy = np.array(hdul[1].data)
# timebins [s] ---!
time = np.array(hdul[2].data)
# spectra ---!
spectra = np.array(hdul[3].data)

hdul.close()

Nt = len(time)
Ne = len(energy)


# time grid ---!
t = [0.0 for x in range(Nt + 1)]
for i in range(Nt - 1):
  t[i + 1] = time[i][0] + (time[i + 1][0] - time[i][0]) / 2
# tmax in last bin ---!
t[Nt] = time[Nt - 1][0] + (time[Nt - 1][0] - t[Nt - 1])

# stop the second after higher tmax ---!
tbin_stop = 1
for bin in range(len(t)) :
  if t[bin] <= max(tmax) :
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

print('!!! check start simulations ------ ')

# simulate N independent photon-lists of the same event ---!
for k in range(trials):
  count += 1
  # attach ID to fileroot
  f = fileroot + 'sim%06d' % (count)

  # photon lists for each bin ---!
  for i in range(tbin_stop):

    sim = ctools.ctobssim()
    sim["inmodel"] = datapath + 'run0406_ID000126_tbin%02d.xml' % i
    sim["outevents"] = simpath + f + "_tbin%02d.fits" % i
    sim["caldb"] = "prod3b"
    sim["irf"] = "South_z40_average_100s"
    sim["ra"] = pointRA
    sim["dec"] = pointDEC
    sim["rad"] = 5.0
    sim["tmin"] = t[i]
    sim["tmax"] = t[i + 1]
    sim["emin"] = 0.03
    sim["emax"] = 0.5
    sim["seed"] = count
    sim["logfile"] = simpath + f + "_tbin%02d.log" % i
    sim["debug"] = False
    sim.execute() 

  # combine in observatiion list
  xml = gammalib.GXml()
  obslist = xml.append('observation_list title="observation library"')

  for i in range(tbin_stop):
    obs = obslist.append('observation name="run0406_sim%06d" id="%02d" instrument="CTA"' % (count,i))
    obs.append('parameter name="EventList" file="%s%s_tbin%02d.fits"' % (simpath, f, i))

  eventList = '%sobs_%s.xml' % (detpath, f)
  xml.save(eventList)
  print('!!! check --- eventList: ', eventList)

  # assign eventList to observation list ---!

  # ============================
  # !!! EVENT LIST SELECTION !!!
  # ============================

  selectedEvents = []

  for i in range(tint) :

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
    selection['debug'] = False
    selection.execute()

  print('!!! check --- selection: ', selectedEvents)


  # ========================
  # !!! SELECTION SKYMAP !!!
  # ========================

  skymapName = []

  for i in range(tint) :
    skymapName.append(selectedEvents[i].replace(selectpath, detpath).replace('.xml', '_skymap.fits'))

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
    skymap['debug'] = False
    skymap.execute() 

  print('!!! check --- skymaps: ', skymapName)

  # ==============================
  # !!! DETECTION AND MODELING !!!
  # ==============================

  detXml = []
  detReg = []
  pos = []

  for i in range(tint) :
    det, reg, coord = srcDetection_spcModeling(skymapName[i], sigma=sigma, maxSrc=10)
    detXml.append(det)
    detReg.append(reg)
    pos.append(coord)

    print('\n\n==========\n\n!!! check --- detection.............', texp[i],'s done\n\n!!! coords:', pos, '\n\n ==========\n\n')


  # =========================
  # !!! DETECTED RA & DEC !!!
  # =========================

  raDet = []
  decDet = []
  Ndet = []

  for i in range(tint) :
    fixParams(detXml[i], spt_prms=['RA', 'DEC'])
    raDet.append(pos[i][0][0]) if len(pos[i][0]) > 0 else raDet.append(np.nan)
    decDet.append(pos[i][1][0]) if len(pos[i][1]) > 0 else decDet.append(np.nan)
    Ndet.append(len(pos[i][0]))
    print('!!! check number of detections in trial', k+1, ' ====== ', Ndet)

  print('!!! check RA-DET list:', raDet)

  # flatten the lists ---!
  #raSrc001 = [float(el) for grp in raDet for el in grp]
  #decSrc001 = [float(el) for grp in decDet for el in grp]

  # already flattened lists ---!
  raSrc001 = raDet
  decSrc001 = decDet
  print('!!! check ------- 1° src DETECTED RA:', raSrc001, ' and DEC:', decSrc001)


  # ==================
  # !!! LIKELIHOOD !!!
  # ==================

  resultsName = []

  for i in range(tint) :

    resultsName.append(detXml[i].replace('_det%dsgm_Mod.xml' % sigma, '_det%dsgm_results.xml' % sigma))

    if Ndet[i] > 0:
      like = ctools.ctlike()
      like['inobs'] = selectedEvents[i]
      like['inmodel'] = detXml[i]
      like['outmodel'] = resultsName[i]
      like['caldb'] = caldb
      like['irf'] = irf
      like['fix_spat_for_ts'] = False
      like['logfile'] = resultsName[i].replace('.xml', '.log')
      like['debug'] = False  # default
      like.execute()


  print('!!! check --- max likelihoods: ', resultsName)

  # ======================
  # !!! LIKELIHOOD TSV !!!
  # ======================

  tsList = []

  for i in range(tint):
    if Ndet[i] > 0:
      tsList.append(getTSV(resultsName[i]))
    else:
      tsList.append([np.nan])

  # flatten list ---!
  #tsv = [float(el) for grp in tsList for el in grp]

  print('!!! check --- TSV List for each texp: ', tsList)

  # only first elem ---!
  tsv = []
  for i in range(len(tsList)):
    tsv.append(tsList[i][0])

  print('!!! check --- SRC001 TSV for each texp:', tsv)

  # count src with TS >= 9
  Nsrc = []
  for i in range(tint):
    n = 0
    for j in range(len(tsList[i])):
      if float(tsList[i][j]) >= 9 :
        n += 1

    Nsrc.append(n) 
  
    print('!!! check ---- SRC NUMBER with TS > 9 for trial ', k+1, ' ===== ', Nsrc[i], ' for texp ==== ', texp[i])

  # ===========================
  # !!! LIKELIHOOD RA & DEC !!!
  # ===========================


  raList = []
  decList = []
  for i in range(tint):
    if Ndet[i] > 0:
      raList.append(getRaDec(resultsName[i])[0])
      decList.append(getRaDec(resultsName[i])[1])
    else:
      raList.append([np.nan])
      decList.append([np.nan])

  # flatten lists ---!
  # raFit = [float(el) for grp in raList for el in grp]
  # decFit = [float(el) for grp in decList for el in grp]

  # already flattened lists ---!
  raFit = []
  decFit = []
  for i in range(len(raList)):
    raFit.append(raList[i][0])
    decFit.append(decList[i][0])


  print('!!! check --- RA FIT for each texp: ', raList, '\n\nand SRC001 for each texp:', raFit)
  print('!!! check --- RA FIT for each texp: ', decList, '\n\nand SRC001 for each texp:', decFit)

  # =================================
  # !!! LIKELIHOOD RA & DEC ERROR !!!
  # =================================


  raErr = []
  decErr = []
  for i in range(tint):
    if Ndet[i] > 0:
      raErr.append(getRaDec_errors(resultsName[i])[0])
      decErr.append(getRaDec_errors(resultsName[i])[1])
    else:
      raErr.append([np.nan])
      decErr.append([np.nan])

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

  csvName = []
  header = '#trial,t exp,sigma,Ndet,Nsrc,RA Src001,DEC Src001,RA Fit,DEC Fit,RA err,DEC err,MEAN err,TSV\n'
  ID = 'ID%06d' % count
  
  for i in range(tint):
    csvName.append(csvpath + fileroot + 'v07_%ds_chunk%02d.csv' % (texp[i], chunk))

    row = []
    row.append([ID, texp[i], sigma, Ndet[i], Nsrc[i], raSrc001[i], decSrc001[i], raFit[i], decFit[i], raFit_err[i], decFit_err[i], mean_err[i], tsv[i]])

    print('!!! check row --- iter', i, '=====', row)

    if os.path.isfile(csvName[i]) == True :
      with open(csvName[i], 'a') as f:
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    else :
      with open(csvName[i], 'w+') as f:
        f.write(header)
        w = csv.writer(f)
        w.writerows(row)
        f.close()
    
  print('!!! check --- data file: ', csvName)

  # ==================== #
  # !!! DELETE SIM FILES #
  # ==================== #


  print('!!! check --- ', count, ') trial done...')
  if count > 1:
    os.system('rm ' + simpath + '*sim%06d*' % count)
    os.system('rm ' + selectpath + '*sim%06d*' % count)
    os.system('rm ' + detpath + '*sim%06d*' % count)

print('!!! check end\n\ndone......chunk ', chunk, 'sim id from ', trials*(chunk-1)+1 , ' to ', count)
print('!!! check end\n\ndone...... removed all files in sim, selected_sim and detection_all')
