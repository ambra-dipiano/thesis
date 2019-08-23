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
workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/run0406_bkg/'
runpath = workdir + 'run0406_ID000126/'
simpath = runpath + 'sim/'
selectpath = runpath + 'selected_sim/'
datapath = runpath + 'data/'
csvpath = runpath + 'csv/'
detpath = runpath + 'detection_all/'

fileroot = 'run0406_'
bkg = workdir + 'model_CTAIrfBackground.xml'
grb = workdir + 'model_run0406_ID000126.xml'

# =====================
# !!! SET UP TRIALS !!!
# =====================


sigma = 5
texp = [1]
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
pointRA = offmax[0]
pointDEC = offmax[1]


#print('!!! check start simulations ------ ')

# simulate N independent photon-lists of the same event ---!
for k in range(trials):
  count += 1
  # attach ID to fileroot
  f = fileroot + 'bkg%09d.fits' % (count)
  event = simpath + f 

  sim = ctools.ctobssim()
  sim["inmodel"] = bkg
  sim["outevents"] = event
  sim["caldb"] = "prod3b"
  sim["irf"] = "South_z40_average_100s"
  sim["ra"] = pointRA
  sim["dec"] = pointDEC
  sim["rad"] = 5.0
  sim["tmin"] = 0
  sim["tmax"] = 35
  sim["emin"] = 0.03
  sim["emax"] = 1.0
  sim["seed"] = count
  sim["logfile"] = event.replace('.fits', '.log')
  sim.execute() 

  # combine in observatiion list
  xml = gammalib.GXml()
  obslist = xml.append('observation_list title="observation library"')

  obs = obslist.append('observation name="BKG" id="bkg%09d" instrument="CTA"' % count)
  obs.append('parameter name="EventList" file="%s%s"' % (simpath, f))

  eventList = '%sobs_%s' % (detpath, f.replace('.fits', '.xml'))
  xml.save(eventList)
#  print('!!! check --- eventList: ', eventList)

  # assign eventList to observation list ---!

  # ============================
  # !!! EVENT LIST SELECTION !!!
  # ============================

  selectedEvents =  eventList.replace('obs_', 'texp%ds_' % texp[0]).replace(simpath, selectpath)

  selection = ctools.ctselect()
  selection['inobs'] = eventList
  selection['outobs'] = selectedEvents
  selection['usepnt'] = bool('yes')
  selection['prefix'] = selectpath + 'texp%ds_' % texp[0]
  selection['rad'] = roi
  selection['tmin'] = tmin
  selection['tmax'] = tmax[0]
  selection['emin'] = emin
  selection['emax'] = emax
  selection['logfile'] = selectedEvents[i].replace('.xml', '.log')
  selection['debug'] = bool('no')
  selection.execute()

#  print('!!! check --- selection: ', selectedEvents)


  # ===========================
  # !!! FREE/FIX PARAMETERS !!!
  # ===========================

  if count == 0 :
    freeParams(grb, spc_prms=['Prefactor', 'Index'], spt_prm=['RA', 'DEC'], free_bkg=True, bkg_prms=['Prefactor', 'Index'])
#    print('!!! check ------- freed params')
    fixParams(grb, spc_prms=['PivotEnergy'], fix_bkg=True, bkg_prms=['PivotEnergy'])
#    print('!!! check ------- fixed params')


  # ==================
  # !!! LIKELIHOOD !!!
  # ==================

  resultsName = selectedEvents.replace('.xml', '_results.xml').replace(selectpath, simpath)

  like = ctools.ctlike()
  like['inobs'] = selectedEvents
  like['inmodel'] = grb
  like['outmodel'] = resultsName
  like['caldb'] = caldb
  like['irf'] = irf
  like['fix_spat_for_ts'] = bool('yes')
  like['logfile'] = resultsName.replace('.xml', '.log')
  like['debug'] = bool('no')  # default
  like.execute()


#  print('!!! check --- max likelihoods: ', resultsName)

  # ======================
  # !!! LIKELIHOOD TSV !!!
  # ======================

  tsList = getTSV(resultsName)

#  print('!!! check --- TSV List for each texp: ', tsList)

  # only first elem ---!
  if len(tsList) > 0 :
    tsv = tsList[0] 
  else:
    tsv = 0.0

#  print('!!! check --- SRC001 TSV for each texp:', tsv)



  # ================== #
  # !!! WRITE CSV FILE #
  # ================== #

  csvName = []
  header = '#trial,t exp,TSV\n'
  ID = 'BKG%09d' % count
  
  csvName.append(csvpath + fileroot + 'v07_%ds_chunk%02d.csv' % (texp[0], chunk))
  #csvName.append(csvpath + fileroot + 'v07_%ds.csv' % texp[i])


  row = []
  row.append([ID, texp[0], tsv])

#  print('!!! check row --- iter', i, '=====', row)

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
    
#  print('!!! check --- data file: ', csvName)

  # ==================== #
  # !!! DELETE SIM FILES #
  # ==================== #


#  print('!!! check --- ', count, ') trial done...')

  if count > 10 :
    os.system('rm ' + simpath + '*bkg%09d*' % count)
    os.system('rm ' + selectpath + '*bkg%09d*' % count)
    os.system('rm ' + detpath + '*bkg%09d*' % count)

#print('!!! check end\n\ndone......chunk ', chunk, 'sim id from ', trials*(chunk-1)+1 , ' to ', count)
#print('!!! check end\n\ndone...... removed all files in sim, selected_sim and detection_all')

