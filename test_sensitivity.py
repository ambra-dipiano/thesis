# ===============================
# !!! DEGRADE IRF SENSITIVITY !!!
# ===============================

import pandas as pd

from module_plot import *
from pkg_blindsearch import *

# files and path ---!
path = '/home/ambra/Desktop/cluster-morgana/irf_degraded/'
outpath = path + 'South_z20_average_100s/'
model = '$CTOOLS/share/models/crab.xml'
bkg = outpath + 'irf_bkg.xml'
# irf and caldb ---!
prod_n = 3
caldb_nom = 'prod%db' %prod_n
caldb_deg = 'degr%db' %prod_n
irf = 'South_z20_average_100s'
# nominal outputs ---!
event_nom = outpath+'crab.fits'
output_nom = outpath+'prod%db_sens.csv' %prod_n
results_nom = outpath+'prod%db_results.xml' %prod_n
# degraded outputs ---!
event_deg = outpath+'crab_degraded.fits'
output_deg = outpath+'degr%db_sens.csv' %prod_n
results_deg = outpath+'degr%db_results.xml' %prod_n
# setup ---!
nominal = False
degraded = False
plot = True

caldb, event, results, output = [], [], [], []
if nominal:
  caldb.append(caldb_nom)
  event.append(event_nom)
  results.append(results_nom)
  output.append(output_nom)
  print('compute nominal')
if degraded:
  caldb.append(caldb_deg)
  event.append(event_deg)
  results.append(results_deg)
  output.append(output_deg)
  print('compute degraded')
  irfObj = Analysis()
  irfObj.irf = irf
  irfObj.caldb = caldb
  irfObj.degradeIrf(bkg=False)
e = [0.03, 150.0]
t = [0, 150]
pointing = [83.63, 22.01]

# INITIALIZE ---!
for i in range(len(caldb)):
  nObj = Analysis('/config_irf.xml')
  nObj.e = e
  nObj.t = t
  nObj.caldb = caldb[i]
  nObj.irf = irf
  nObj.model = model
  # NOMINAL SIM ---!
  nObj.output = event[i]
  nObj.pointing = pointing  # (deg)
  nObj.eventSim()
  # NOMINAL MAX LIKELIHOOD ---!
  nObj.input = event[i]
  nObj.output = results[i]
  nObj.maxLikelihood()
  # NOMINAL SENS ---!
  nObj.model = results[i]
  nObj.output = output[i]
  nObj.src_name = 'Crab'
  nObj.eventSens(bins=20)

# ------------------------------------- PLOT --------------------------------------- !!!

if plot:
  pngroot = outpath + caldb_nom + '_' + irf.strip('_')
  savefig1 = pngroot + '_sensdiff.png'
  savefig2 = pngroot + 'sensdiff_phflux.png'
  savefig3 = pngroot + 'sensdiff_ratio.png'
  title = caldb_nom + ': ' + irf.replace('_', '\_')

  # nominal
  df_nom = pd.read_csv(output_nom)
  cols = list(df_nom.columns)
  energy_nom = np.array(df_nom[cols[0]])
  sens_nom = np.array(df_nom[cols[6]])
  flux_nom = np.array(df_nom[cols[4]])
  # degraded ---!
  df_deg = pd.read_csv(output_deg)
  energy_deg = np.array(df_deg[cols[0]])
  sens_deg = np.array(df_deg[cols[6]])
  flux_deg = np.array(df_deg[cols[4]])

  # check ---!
  print('nominal, degraded')
  print('min:', sens_nom.min(), sens_deg.min())
  print('max:', sens_nom.max(), sens_deg.max())
  print('bin5:', sens_nom[5], sens_deg[5])
  print('ratios:', sens_nom/sens_deg)

  showSensitivity([10**energy_nom, 10**energy_deg], [sens_nom, sens_deg], savefig1, marker=['+', 'x'],
                  xlabel='energy (TeV)', ylabel='E$^2$F sensitivity (erg/cm$^2$/s)',
                  label=['full sens irf', 'degraded irf'], title=title, fontsize=12)

  showSensitivity([10**energy_nom, 10**energy_deg], [flux_nom, flux_deg], savefig2, marker=['+', 'x'],
                  xlabel='energy (TeV)', ylabel='ph flux (ph/cm$^2$/s)',
                  label=['full sens irf', 'degraded irf'], title=title, fontsize=12)

  # plot ratio ---!
  plt.figure()
  ax = plt.subplot(111)
  plt.title('ratio '+ irf.replace('_', '\_'))
  plt.xlabel('energy (TeV)')
  plt.ylabel('ratio nominal/degraded')
  plt.plot(10**energy_nom, sens_nom/sens_deg)
  plt.axhline(0.5, c='r', ls='-.')
  #plt.ylim(0., 1.)
  plt.show()
  plt.savefig(savefig3)