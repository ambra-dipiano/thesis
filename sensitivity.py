# MIT License
# Copyright (c) 2019, 2020 Ambra Di Piano
# ---------------------------------------
# =======================
# !!! IRF SENSITIVITY !!!
# =======================

import pandas as pd

from module_plot import *
from pkg_blindsearch import *

# files and path ---!
model = '/home/ambra/Desktop/cluster-morgana/grb.xml'
cfg = xmlConfig('/config_irf.xml')
p = ConfigureXml(cfg=cfg)
path = p.getWorkingDir()
outpath = p.getRunDir()
pngpath = p.getPngDir()
# irf and caldb ---!
caldb_nom = 'prod3b-v2'
caldb_deg = caldb_nom.replace('prod', 'degr')
irf = 'South_z40_0.5h'
# setup ---!
nominal = True
degraded = True
compute = True
plot = True
sens_type = 'Differential'
print(sens_type, 'sensitivity')

caldb = []
if nominal:
  caldb.append(caldb_nom)
  print('Nominal cladb')
if degraded:
  caldb.append(caldb_deg)
  print('Degraded caldb')

e = [0.03, 150.0]
texp = [1, 5, 10, 100]
pointing = (33.057, -51.841) # pointing direction RA/DEC (deg) - centered on the source
nbins = 20 # energy bins for sensitivity computation
src_name = 'GRB'

# INITIALIZE ---!
if compute:
  for i in range(len(caldb)):
    for j in range(len(texp)):
      print('texp=', texp[j])
      event = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_phlist.fits'
      results = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_maxlike.xml'
      output = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_sens.csv'
      nObj = Analysis('/config_irf.xml')
      nObj.e = e
      nObj.t = [0, texp[j]]
      nObj.caldb = caldb[i]
      nObj.irf = irf
      nObj.model = model
      # NOMINAL SIM ---!
      nObj.output = event
      nObj.pointing = pointing  # (deg)
      nObj.eventSim()
      # NOMINAL MAX LIKELIHOOD ---!
      nObj.input = event
      nObj.output = results
      nObj.maxLikelihood()
      # NOMINAL SENS ---!
      nObj.sens_type = sens_type
      nObj.model = results
      nObj.output = output
      nObj.src_name = src_name
      nObj.eventSens(bins=nbins)

