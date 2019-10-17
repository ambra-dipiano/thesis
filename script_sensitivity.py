# ===============================
# !!! DEGRADE IRF SENSITIVITY !!!
# ===============================

from pkg_blindsearch import *
from module_plot import *
import pandas as pd

path = '/home/ambra/Desktop/cluster-morgana/irf_degraded/'
outpath = path + 'South_z20_average_100s/'
model = '$CTOOLS/share/models/crab.xml'
bkg = outpath + 'irf_bkg.xml'

prod_n = 3
caldb_nom = 'prod%db' %prod_n
caldb_deg = 'degr%db' %prod_n
irf = 'South_z20_average_100s'

event_nom = outpath+'crab.fits'
event_deg = outpath+'crab_degraded.fits'
output_nom = outpath+'prod%db_sens.csv' %prod_n
output_deg = outpath+'degr%db_sens.csv' %prod_n
results_nom = outpath+'prod%db_results.xml' %prod_n
results_deg = outpath+'degr%db_results.xml' %prod_n

# INITIALIZE ---!
nObj = analysis('/config_irf.xml')
nObj.e = [0.03, 150.0]
nObj.t = [0, 100]
nObj.caldb = caldb_nom
nObj.irf = irf
nObj.model = bkg
# NOMINAL SIM ---!
nObj.output = event_nom
nObj.pointing = [83.63, 22.01]  # (deg)
nObj.eventSim()
# NOMINAL MAX LIKELIHOOD ---!
nObj.input = event_nom
nObj.output = results_nom
nObj.maxLikelihood()
# NOMINAL SENS ---!
nObj.model = results_nom
nObj.output = output_nom
nObj.srcName = 'CTABackgroundModel'
nObj.eventSens()
print('crab nominal done')

# INITIALIZE ---!
dObj = analysis('/config_irf.xml')
dObj.e = [0.03, 150.0]
dObj.t = [0, 100]
dObj.caldb = caldb_nom
dObj.irf = irf
dObj.model = bkg
# DEGRADE IRF ---!
dObj.degradeIRF()
dObj.irf = 'North_z20_average_100s'
# DEGRADE SIM ---!
dObj.output = event_deg
dObj.pointing = [83.63, 22.01]  # (deg)
dObj.eventSim()
# NOMINAL MAX LIKELIHOOD ---!
dObj.input = event_deg
dObj.output = results_deg
dObj.maxLikelihood()
# DEGRADE SENS ---!
dObj.model = results_deg
dObj.output = output_deg
dObj.srcName = 'CTABackgroundModel'
dObj.eventSens()
print('Crab degraded done')

# PLOT ---!

savefig1 = outpath + 'sensitivity_differential.png'
savefig2 = outpath + 'sens_phflux_differential.png'

df_nom = pd.read_csv(output_nom)
df_deg = pd.read_csv(output_deg)
cols = list(df_nom.columns)
energy_nom = np.array(df_nom[cols[0]])
energy_deg = np.array(df_deg[cols[0]])
sens_nom = np.array(df_nom[cols[6]])
sens_deg = np.array(df_deg[cols[6]])

flux_nom = np.array(df_nom[cols[4]])
flux_deg = np.array(df_deg[cols[4]])

print(sens_nom.max(), sens_deg.max())
print(sens_nom[5], sens_deg[5])

# showSensitivity([10**energy_nom, 10**energy_deg], [sens_nom, sens_deg], savefig1, marker=['+', 'x'],
#                 xlabel='energy (TeV)', ylabel='E$^2$F sensitivity (erg/cm$^2$/s)',
#                 label=['full sens irf', 'degraded irf'], fontsize=12)
#
# showSensitivity([10**energy_nom, 10**energy_deg], [flux_nom, flux_deg], savefig2, marker=['+', 'x'],
#                 xlabel='energy (TeV)', ylabel='ph flux (ph/cm$^2$/s)',
#                 label=['full sens irf', 'degraded irf'], fontsize=12)