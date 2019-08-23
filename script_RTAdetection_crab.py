#!/bin/python3.6

# ===================== #
# TESTING WILKS THEOREM #
# ===================== #

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

# initialize global count ---!
chunk = int(sys.argv[1]) # global count
trials = int(sys.argv[2]) # number of trials
count = int(sys.argv[3]) # starting count  

caldb = 'prod3b'
irf = 'South_z40_average_100s'

# work with absolute paths ---!
workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/crab/'
runpath = workdir 
simpath = runpath + 'sim/'
selectpath = runpath + 'selected_sim/'
datapath = runpath + 'data/'
csvpath = runpath + 'csv/'
detpath = runpath + 'detection_all/'

fileroot = 'crab_'
crab = datapath + 'crab.xml'


# =====================
# !!! SET UP TRIALS !!!
# =====================


sigma = 5
#texp = [1, 5, 10, 100]
texp = [2000]
#texp.sort()
#tint = len(texp)
#tmin = 30
#tmax = []
#for i in range(tint):
#  tmax.append(tmin + texp[i])

emin = 0.03  # 300 GeV
emax = 100.0  # 500 GeV
roi = 5
wbin = 0.02
nbin = int(roi/wbin)

tstart = 0
tstop = 2000


# ====================
# !!! SIMULATE GRB !!!
# ====================

# pointing with off-axis equal to max prob GW ---!
offmax = [0, 0]
true = [83.633, 22.014]
pointRA = true[0] + offmax[0]
pointDEC = true[1] + offmax[1]


print('!!! check start simulations ------ ')

# simulate N independent photon-lists of the same event ---!
for k in range(trials):
  count += 1
  # attach ID to fileroot
  f = fileroot + 'sim%06d.fits' % (count)
  event = simpath + f 

  sim = ctools.ctobssim()
  sim["inmodel"] = crab
  sim["outevents"] = event
  sim["caldb"] = "prod3b"
  sim["irf"] = "South_z40_average_100s"
  sim["ra"] = pointRA
  sim["dec"] = pointDEC
  sim["rad"] = 5.0
  sim["tmin"] = tstart
  sim["tmax"] = tstop
  sim["emin"] = emin
  sim["emax"] = emax
  sim["seed"] = count
  sim["logfile"] = event.replace('.fits', '.log')
  sim.execute() 

  # combine in observatiion list
  xml = gammalib.GXml()
  obslist = xml.append('observation_list title="observation library"')

  obs = obslist.append('observation name="Crab" id="sim%06d" instrument="CTA"' % count)
  obs.append('parameter name="EventList" file="%s%s"' % (simpath, f))

  eventList = '%sobs_%s' % (detpath, f.replace('.fits', '.xml'))
  xml.save(eventList)
  print('!!! check --- eventList: ', eventList)

  # assign eventList to observation list ---!


  # ========================
  # !!! SELECTION SKYMAP !!!
  # ========================

  skymapName = event.replace(simpath, detpath).replace('.fits', '_skymap.fits')

  skymap = ctools.ctskymap()
  skymap['inobs'] = event
  skymap['outmap'] = skymapName
  skymap['irf'] = irf
  skymap['caldb'] = caldb
  skymap['emin'] = emin
  skymap['emax'] = emax
  skymap['usepnt'] = bool('yes')
  skymap['nxpix'] = nbin
  skymap['nypix'] = nbin
  skymap['binsz'] = wbin
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


  detXml, detReg, pos = srcDetection_spcModeling(skymapName, sigma=sigma, maxSrc=10)

  print('\n\n==========\n\n!!! check --- detection.............', texp[0],'s done\n\n!!! coords:', pos, '\n\n ==========\n\n')


  # ===========================
  # !!! FREE/FIX PARAMETERS !!!
  # ===========================


  freeParams(detXml, spc_prms=['Prefactor'])
  print('!!! check ------- freed params')
  fixParams(detXml, spc_prms=['PivotEnergy', 'Index'], spt_prms=['RA', 'DEC'])
  print('!!! check ------- fixed params')


  # ==================
  # !!! LIKELIHOOD !!!
  # ==================

  resultsName = event.replace(simpath, detpath).replace('.fits', '_results.xml')

  like = ctools.ctlike()
  like['inobs'] = event
  like['inmodel'] = crab
  like['outmodel'] = resultsName
  like['caldb'] = caldb
  like['irf'] = irf
  like['fix_spat_for_ts'] = bool('yes') # only if you didn't fix the bkg
  like['logfile'] = resultsName.replace('.xml', '.log')
  like['debug'] = bool('no')  # default
  like.execute()


  print('!!! check --- max likelihoods: ', resultsName)

  # ======================
  # !!! LIKELIHOOD TSV !!!
  # ======================

  tsList = getTSV(resultsName)

  print('!!! check --- TSV List for each texp: ', tsList)

  # only first elem ---!
  if len(tsList) > 0 :
    tsv = tsList[0] 
  else:
    tsv = np.nan

  print('!!! check --- SRC001 TSV for each texp:', tsv)

  # ===========================
  # !!! LIKELIHOOD RA & DEC !!!
  # ===========================


  raList = getRaDec(resultsName)[0]
  decList = getRaDec(resultsName)[1]

  # already flattened lists ---!
  if len(raList) > 0 :
    raFit = raList[0]
    decFit = decList[0]
  else :
    raFit = np.nan
    decList = np.nan


  print('!!! check --- RA FIT for each texp: ', raList, '\n\nand SRC001 for each texp:', raFit)
  print('!!! check --- RA FIT for each texp: ', decList, '\n\nand SRC001 for each texp:', decFit)

  # ================== #
  # !!! WRITE CSV FILE #
  # ================== #

  csvName = []
  header = '#trial,t exp,RA fit,DEC fit,TSV\n'
  ID = 'ID%06d' % count
  
  csvName.append(csvpath + fileroot + '%ds_chunk%02d.csv' % (texp[0], chunk))
  #csvName.append(csvpath + fileroot + 'v07_%ds.csv' % texp[i])


  row = []
  row.append([ID, texp[0], raFit, decFit, tsv])

  print('!!! check row --- iter', k, '=====', row)

  if os.path.isfile(csvName[0]) == True :
    with open(csvName[0], 'a') as f:
      w = csv.writer(f)
      w.writerows(row)
      f.close()
  else :
    with open(csvName[0], 'w+') as f:
      f.write(header)
      w = csv.writer(f)
      w.writerows(row)
      f.close()
    
  print('!!! check --- data file: ', csvName)

  # ==================== #
  # !!! DELETE SIM FILES #
  # ==================== #


  print('!!! check --- ', count, ') trial done...')

  if count > 4 :
    os.system('rm ' + simpath + '*sim%06d*' % count)
    os.system('rm ' + selectpath + '*sim%06d*' % count)
    os.system('rm ' + detpath + '*sim%06d*' % count)

print('!!! check end\n\ndone......chunk ', chunk, 'sim id from ', trials*(chunk-1)+1 , ' to ', count)
print('!!! check end\n\ndone...... removed all files in sim, selected_sim and detection_all except seeds 1-4')

