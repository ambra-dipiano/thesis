#!/bin/python3.6

# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMPORTS ---!
import gammalib
import ctools
import cscripts
from astropy.io import fits
#from module_plot import showSkymap
from module_xml import *
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

# test with 10 sec ---!
texp = int(sys.argv[1]) # s
chunk = int(sys.argv[2]) # global count

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
csvName = csvpath + fileroot + 'v06_%ds.csv' % texp



# =====================
# !!! SET UP TRIALS !!!
# =====================

# trials
trials = 10
# trials count
count = trials*(chunk-1)

sigma = 5
tmin = 30
tmax = tmin + texp

emin = 0.03  # 300 GeV
emax = 0.5  # 500 GeV
roi = 5

tsList = []
raList = []
decList = []
raErr_list = []
decErr_list = []

# ====================
# !!! SIMULATE GRB !!!
# ====================

# pointing with off-axis equal to max prob GW ---!
offmax = [31.582, -53.211]



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

print('!!! check start simulations ------ ')

# simulate N independent photon-lists of the same event ---!
for k in range(trials):
  count += 1
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
    sim.execute() 

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

  # ============================
  # !!! EVENT LIST SELECTION !!!
  # ============================

  selectedEvents = eventList.replace('obs_', 'texp%ds_' % texp)

  selection = ctools.ctselect()
  selection['inobs'] = eventList
  selection['outobs'] = selectedEvents
  selection['usepnt'] = bool('yes')
  selection['prefix'] = selectpath + 'texp%ds_' % texp
  selection['rad'] = roi
  selection['tmin'] = tmin
  selection['tmax'] = tmax
  selection['emin'] = emin
  selection['emax'] = emax
  selection['logfile'] = selectedEvents.replace('.xml', '.log')
  selection['debug'] = bool('no')
  selection.execute()

  print('!!! check --- selection: ', selectedEvents)


  # ========================
  # !!! SELECTION SKYMAP !!!
  # ========================

  skymapName = selectedEvents.replace(selectpath, detpath).replace('.xml', '_skymap.fits')

  skymap = ctools.ctskymap()
  skymap['inobs'] = selectedEvents
  skymap['outmap'] = skymapName
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
  skymap['logfile'] = skymapName.replace('.fits', '.log')
  skymap['debug'] = bool('no')
  skymap.execute() 

  print('!!! check --- skymaps: ', skymapName)

  # ==============================
  # !!! DETECTION AND MODELING !!!
  # ==============================



  detXml, detReg, pos = srcDetection_spcModeling(skymapName, sigma=sigma, maxSrc=10, exclrad=0.2)
  #det, reg, coord = srcDetection_spcModeling(skymapName, sigma=sigma, maxSrc=10)
  #detXml.append(det)
  #detReg.append(reg)
  #pos.append(coord)

  print('!!! check --- detection xml: ', detXml)

  # =========================
  # !!! DETECTED RA & DEC !!!
  # =========================

  if len(pos[0]) > 0 :
    raSrc001 = pos[0][0]
    decSrc001 = pos[1][0]
  else :
    raSrc001 = 'NAN'
    decSrc001 = 'NAN'

  Ndet = len(pos[0])


  # ==================
  # !!! LIKELIHOOD !!!
  # ==================

  resultsName = detXml.replace('Mod', '_results')

  if Ndet > 0:
    like = ctools.ctlike()
    like['inobs'] = selectedEvents
    like['inmodel'] = detXml
    like['outmodel'] = resultsName
    like['caldb'] = caldb
    like['irf'] = irf
    like['fix_spat_for_ts'] = bool('yes')
    like['logfile'] = resultsName.replace('.xml', '.log')
    like['debug'] = bool('no')  # default
    like.execute()
    i += 1

  print('!!! check --- max likelihoos: ', resultsName)

  # ======================
  # !!! LIKELIHOOD TSV !!!
  # ======================




  if Ndet > 0:
    tsList.append(getTSV(resultsName))
  else:
    tsList.append('NAN')

  # only first elem ---!
  tsv = tsList[0][0]

  print('!!! check --- TSV List: ', tsList)

  # count src with TS >= 9
  Nsrc = 0
  for i in range(len(tsList)):
    if float(tsList[i][0]) >= 9 :
       Nsrc += 1 
  


  # ===========================
  # !!! LIKELIHOOD RA & DEC !!!
  # ===========================


  if Ndet > 0:
    raList.append(getRaDec(resultsName)[0])
    decList.append(getRaDec(resultsName)[1])
  else:
    raList.append(['NAN'])
    decList.append(['NAN'])

  raFit = raList[0][0]
  decFit = decList[0][0]


  # =================================
  # !!! LIKELIHOOD RA & DEC ERROR !!!
  # =================================



  if Ndet > 0:
    raErr_list.append(getRaDec_errors(resultsName)[0])
    decErr_list.append(getRaDec_errors(resultsName)[1])
  else:
    raErr_list.append(['NAN'])
    decErr_list.append(['NAN'])

  raFit_err = raErr_list[0][0]
  decFit_err = decErr_list[0][0]
  mean_err = 0.5 * (float(raFit_err) + float(decFit_err))


  # ================== #
  # !!! WRITE CSV FILE #
  # ================== #

  row = []
  header = '#trial,t int,sigma,Nsrc,RA Src001,DEC Src001,RA Fit,DEC Fit,RA err,DEC err,MEAN err,TSV\n'

  row.append(['ID%06d' % count, texp,sigma,Nsrc,raSrc001,decSrc001,raFit,decFit,raFit_err,decFit_err,mean_err,tsv])

  print("before file management")

  if os.path.isfile(csvName) == True :
    with open(csvName, 'a') as f:
      w = csv.writer(f)
      w.writerows(row)
      f.close()
  else :
    with open(csvName, 'w+') as f:
      f.write(header)
      w = csv.writer(f)
      w.writerows(row)
      f.close()

  print('!!! check --- data file: ', csvName)

  # ==================== #
  # !!! DELETE SIM FILES #
  # ==================== #


  print('!!! check --- ', count, ') trial done...')


#os.system('rm ' + simpath + '*')
#os.system('rm ' + selectpath + '*')
#os.system('rm ' + detpath + '*')
print('done......chunk ', chunk, 'sim id from ', 10*(chunk-1)+1 , ' to ', count)
