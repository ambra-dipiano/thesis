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

count=1

# =========================== #
# !!! PROPER CODE STARTS HERE #
# =========================== #

caldb = 'prod3b'
irf = 'South_z40_average_30m'

# work with absolute paths ---!
workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/run0406_lc/'
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
texp = 2000
tmin = 0
tmax = tmin+texp
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
  if t[bin] <= tmax :
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
f = fileroot + '%ds' % (texp)

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
    obs = obslist.append('observation name="run0406_sim%06d" id="%02d" instrument="CTA"' % (count,i))
    obs.append('parameter name="EventList" file="%s%s_tbin%02d.fits"' % (simpath, f, i))

  eventList = '%sobs_%s.xml' % (detpath, f)
  xml.save(eventList)
print('!!! check --- eventList: ', eventList)


# ==============
# !!! SKYMAP !!!
# ==============

skymapName = eventList.replace('obs_', '').replace('.xml', '_skymap.fits')

skymap = ctools.ctskymap()
skymap['inobs'] = eventList
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

detXml, detReg, pos = srcDetection_spcModeling(skymapName, sigma=sigma, maxSrc=10)

print('\n\n==========\n\n!!! check --- detection.............', texp,'s done\n\n!!! coords:', pos, '\n\n ==========\n\n')


# ==================
# !!! LIKELIHOOD !!!
# ==================


resultsName = eventList.replace(simpath, detpath).replace('.xml', '_results.xml' )

like = ctools.ctlike()
like['inobs'] = eventList
like['inmodel'] = detXml
like['outmodel'] = resultsName
like['caldb'] = caldb
like['irf'] = irf
like['fix_spat_for_ts'] = bool('yes')
like['logfile'] = resultsName.replace('.xml', '.log')
like['debug'] = bool('no')  # default
like.execute()


print('!!! check --- max likelihoods: ', resultsName)

# ===========================
# !!! LIKELIHOOD RA & DEC !!!
# ===========================


raFit = getRaDec(resultsName)[0]
decFit = getRaDec(resultsName)[1]


print('!!! check --- RA FIT: ', raFit)
print('!!! check --- DEC FIT: ', decFit)

# ==================
# !!! LIGHTCURVE !!!
# ==================

wbin = [1, 5, 10, 100]
nbin = []
for i in range(len(wbin)) :
  nbin.append(int(texp/wbin[i]))

lc = []

for i in range(len(wbin)) :
  lc.append(resultsName.replace('results', 'lightcurve_%ds' %wbin[i]).replace('.xml', '.fits'))
  print('!!! check start LC n.', i+1, 'name:', lc[i])

  lightcurve = cscripts.cslightcrv()
  lightcurve['inobs'] = eventList
  lightcurve['inmodel'] = resultsName
  lightcurve['srcname'] = 'Src001'
  lightcurve['caldb'] = caldb
  lightcurve['irf'] = irf
  lightcurve['outfile'] = lc[i]
  lightcurve['tbinalg'] = 'LIN' # <FILE|LIN|GTI>
  lightcurve['edisp'] = bool('yes')
  lightcurve['tmin'] = tmin
  lightcurve['tmax'] = tmax
  lightcurve['tbins'] = nbin[i]
  lightcurve['method'] = '3D'
  lightcurve['emin'] = emin
  lightcurve['emax'] = emax
  lightcurve['enumbins'] = 0 # 0 for unbinned 3D only
  lightcurve['xref'] = float(raFit[0])
  lightcurve['yref'] = float(decFit[0])
  lightcurve['logfile'] = lc[i].replace('.fits', '.log')
  lightcurve['debug'] = bool('no') 
  lightcurve.execute()

print('!!! check ------- lc n.', i+1, ') completed:', lc[i])
