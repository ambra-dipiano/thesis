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
sens_type = 'Integral'

caldb = []
if nominal:
  caldb.append(caldb_nom)
  print('compute nominal')
if degraded:
  caldb.append(caldb_deg)
  print('compute degraded')

e = [0.03, 150.0]
texp = [1, 5, 10, 100]
pointing = (33.057, -51.841) # pointing direction RA/DEC (deg)

# INITIALIZE ---!
if compute:
  for i in range(len(caldb)):
    print('caldb:', caldb[i])
    for j in range(len(texp)):
      print('texp=', texp[j])
      event = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_crab.fits'
      results = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_crab_like.xml'
      output = outpath + 'texp%ds_' %texp[j] + caldb[i] + '_crab_sens.csv'
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
      print('sim')
      # NOMINAL MAX LIKELIHOOD ---!
      nObj.input = event
      nObj.output = results
      nObj.maxLikelihood()
      print('max like')
      # NOMINAL SENS ---!
      nObj.sens_type = sens_type
      nObj.model = results
      nObj.output = output
      nObj.src_name = 'GRB'
      if sens_type.capitalize() == 'Integral':
        print('integral')
        nObj.eventSens(bins=0)
      else:
        print('differential')
        nObj.eventSens(bins=20)
      print('sens')

# ------------------------------------- PLOT --------------------------------------- !!!

if plot:
  csv = [[], []]
  savefig1, savefig2, savefig3 = [], [], []
  list_sens_nom, list_flux_nom, list_sens_deg, list_flux_deg = [], [], [], []
  for i in range(len(caldb)):
    for j in range(len(texp)):
      csv[i].append(outpath + 'texp%ds_' %texp[j] + caldb[i] + '_crab_sens.csv')
      pngroot = caldb[i] + '_texp%ds' %texp[j]
      if sens_type.capitalize() != 'Integral':
        savefig1.append(pngpath + pngroot + '_sensDiff.png')
        savefig2.append(pngpath + pngroot + '_sensDiff_phflux.png')
        savefig3.append(pngpath + pngroot + '_sensDiff_ratio.png')
      else:
        savefig1.append(pngpath + pngroot + '_sensInt.png')
        savefig2.append(pngpath + pngroot + '_sensInt_phflux.png')
        savefig3.append(pngpath + pngroot + '_sensInt_ratio.png')

  for j in range(len(texp)):
    title = caldb_nom + ': ' + irf.replace('_', '\_') + ' with texp=%ds' %texp[j]
    # nominal
    df_nom = pd.read_csv(csv[0][j])
    cols = list(df_nom.columns)
    energy_nom = np.array(df_nom[cols[0]])
    sens_nom = np.array(df_nom[cols[6]])
    flux_nom = np.array(df_nom[cols[4]])
    # degraded ---!
    df_deg = pd.read_csv(csv[1][j])
    energy_deg = np.array(df_deg[cols[0]])
    sens_deg = np.array(df_deg[cols[6]])
    flux_deg = np.array(df_deg[cols[4]])

    # check ---!
    print('texp%ds: nominal, degraded' %texp[j])
    # print('min:', sens_nom.min(), sens_deg.min())
    # print('max:', sens_nom.max(), sens_deg.max())
    # print('bin5:', sens_nom[5], sens_deg[5])
    print('ratios:', sens_nom/sens_deg)

    showSensitivity([10**energy_nom, 10**energy_deg], [sens_nom, sens_deg], savefig=savefig1[j], marker=['+', 'x'],
                    xlabel='energy (TeV)', ylabel='E$^2$F sensitivity (erg/cm$^2$/s)',
                    label=['full sens irf', 'degraded irf'], title=title, fontsize=12, show=False)

    showSensitivity([10**energy_nom, 10**energy_deg], [flux_nom, flux_deg], savefig=savefig2[j], marker=['+', 'x'],
                    xlabel='energy (TeV)', ylabel='ph flux (ph/cm$^2$/s)',
                    label=['full sens irf', 'degraded irf'], title=title, fontsize=12, show=False)

    # plot ratio ---!
    plt.figure()
    ax = plt.subplot(111, xscale='log')
    plt.title(caldb_nom + ': '+ irf.replace('_', '\_') + ' (texp=%ds)' %texp[j])
    plt.xlabel('energy (TeV)')
    plt.ylabel('ratio nominal/degraded sensitivity')
    plt.plot(10**energy_nom, sens_nom/sens_deg)
    plt.axhline(0.5, c='r', ls='-.')
    # plt.ylim(0., 1.)
    # plt.show()
    plt.savefig(savefig3[j])

    # plot 1 enery bin integral sens: x=texp, y1=sens, y2=phflux ---!
    x = [1, 5, 10, 100]
    y1 = [[6.358567480537711e-11, 2.1582881778139934e-11, 1.3048873052454111e-11, 2.7600318180325007e-12]]
    y2 = [[6.902673455761295e-09, 2.3429740362872147e-09, 1.4165472006464716e-09, 2.9962091973864583e-10]]
    y1.append([3.45734663098577e-11, 1.3048873052446157e-11,7.910828789733936e-12, 1.6811190279480429e-12])
    y2.append([3.753193606911769e-09, 1.4165472006456084e-09, 8.587762584435284e-10, 1.824973270427709e-10])
    l = ['degraded', 'nominal']
    # sens ---!
    plt.figure()
    ax = plt.subplot(111, xscale='log', yscale='log')
    plt.title(caldb_nom + ': ' + irf.replace('_', '\_'))
    plt.xlabel('exposure time (s)')
    plt.ylabel('E$^2$F sensitivity (erg/cm$^2$/s)')
    for i in range(len(y1)):
      plt.plot(x, y1[i], label=l[i])
      plt.scatter(x, y1[i])
    plt.legend()
    plt.savefig(pngpath + caldb_nom + '_sensInt_texp.png')
    # ph flux ---!
    plt.figure()
    ax = plt.subplot(111, xscale='log', yscale='log')
    plt.title(caldb_nom + ': ' + irf.replace('_', '\_'))
    plt.xlabel('exposure time (s)')
    plt.ylabel('ph flux (ph/cm$^2$/s)')
    for i in range(len(y2)):
      plt.plot(x, y2[i], label=l[i])
      plt.scatter(x, y2[i])
    plt.legend()
    plt.savefig(pngpath + caldb_nom + '_fluxInt_texp.png')