# ----------------------------- #
# Copyright 2020 Ambra Di Piano #
# ----------------------------- # -------------------------------------------------- #
# Redistribution and use in source and binary forms, with or without modification,   #
# are permitted provided that the following conditions are met:                      #
# 1. Redistributions of source code must retain the above copyright notice,          #
# this list of conditions and the following disclaimer.                              #
# 2. Redistributions in binary form must reproduce the above copyright notice,       #
# this list of conditions and the following disclaimer in the documentation and/or   #
# other materials provided with the distribution.                                    #
# 3. Neither the name of the copyright holder nor the names of its contributors      #
# may be used to endorse or promote products derived from this software without      #
# specific prior written permission.                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. #
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,   #
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,     #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,      #
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE    #
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  #
# OF THE POSSIBILITY OF SUCH DAMAGE.                                                 #
# ---------------------------------------------------------------------------------- #

# =======================
# !!! IRF SENSITIVITY !!!
# =======================

import pandas as pd
import time

from module_plot import *
from pkg_blindsearch import *

# files and path ---!
cfg = xmlConfig('./config.xml')
p = ConfigureXml(cfg=cfg)
model = p.getModelsDir() + 'crab.xml'
# irf and caldb ---!
caldb_nom = 'prod3b-v2'
caldb_deg = caldb_nom.replace('prod', 'degr')
irf = 'North_z60_0.5h'
# setup ---!
nominal = True
degraded = False
compute = True
plot = True
sens_type = 'integral'
nbins = 1  # energy bins for sensitivity computation

print('Compute', sens_type, 'sensitivity')

caldb = []
if nominal:
    caldb.append(caldb_nom)
    print('Use nominal cladb')
if degraded:
    caldb.append(caldb_deg)
    print('Use degraded caldb')

e = [0.6, 10.0]
tmin = 0
tmax = 30
texp = []
for i in range(tmax-tmin):
    texp.append(i+1)

#pointing = (33.057, -51.841) # pointing direction RA/DEC (deg) - centered on the source
pointing = (83.63, 22.51)
src_name = 'Crab'

# INITIALIZE ---!
if compute:
    for i in range(len(caldb)):
        for j in range(len(texp)):
            print('texp = ', texp[j], ' s')
            event = p.getObsDir() + 'texp%ds_' %texp[j] + caldb[i] + '_phlist.fits'
            results = p.getDetDir() + 'texp%ds_' %texp[j] + caldb[i] + '_maxlike.xml'
            output = p.getCsvDir() + 'texp%ds_' %texp[j] + caldb[i] + '_sens.csv'
            # INITIALIZE ---!
            nObj = Analysis()
            nObj.e = e
            nObj.t = [0, texp[j]]
            nObj.caldb = caldb[i]
            nObj.irf = irf
            nObj.model = model
            # SIM ---!
            nObj.output = event
            nObj.pointing = pointing  # (deg)
            clock = time.time()
            nObj.eventSim()
            print('sim time', time.time() - clock)
            # MAX LIKELIHOOD ---!
            nObj.input = event
            nObj.output = results
            clock = time.time()
            nObj.maxLikelihood()
            print('ctlike time', time.time() - clock)
            # SENS ---!
            nObj.sens_type = sens_type
            nObj.model = results
            nObj.output = output
            nObj.src_name = src_name
            clock = time.time()
            nObj.eventSens(bins=nbins)
            print('sens time', time.time() - clock)

#os.system('rm '+ p.getObsDir() + '*')
#os.system('rm '+ p.getDetDir() + '*maxlike*')

# ------------------------------------- PLOT --------------------------------------- !!!

if plot:
    csv = [[], []]
    savefig1, savefig2, savefig3 = [], [], []
    list_sens_nom, list_flux_nom, list_sens_deg, list_flux_deg = [], [], [], []
    for i in range(len(caldb)):
        for j in range(len(texp)):
            csv[i].append(p.getCsvDir() + 'texp%ds_' %texp[j] + caldb[i] + '_sens.csv')
            pngroot = caldb[i] + '_texp%ds' %texp[j]
            if sens_type.capitalize() != 'Integral':
                savefig1.append(p.getPngDir() + pngroot + '_sensDiff.png')
                savefig2.append(p.getPngDir() + pngroot + '_sensDiff_phflux.png')
                savefig3.append(p.getPngDir() + pngroot + '_sensDiff_vsTime.png')
            else:
                savefig1.append(p.getPngDir() + pngroot + '_sensInt.png')
                savefig2.append(p.getPngDir() + pngroot + '_sensInt_phflux.png')
                savefig3.append(p.getPngDir() + pngroot + '_sensInt_vsTime.png')

    fluxes_nom, fluxes_deg = [], []
    for j in range(len(texp)):
        title = caldb_nom + ': ' + irf.replace('_', '\_') + ' with texp=%ds' %texp[j]
        # nominal
        if nominal:
            df_nom = pd.read_csv(csv[0][j])
            cols = list(df_nom.columns)
            energy_nom = np.array(df_nom[cols[0]])
            sens_nom = np.array(df_nom[cols[6]])
            flux_nom = np.array(df_nom[cols[4]])
            fluxes_nom.append(flux_nom[0])
        # degraded ---!
        if degraded:
            df_deg = pd.read_csv(csv[1][j])
            cols = list(df_deg.columns)
            energy_deg = np.array(df_deg[cols[0]])
            sens_deg = np.array(df_deg[cols[6]])
            flux_deg = np.array(df_deg[cols[4]])
            fluxes_deg.append(flux_deg[0])

        if degraded and nominal:
            showSensitivity([10**energy_nom, 10**energy_deg], [sens_nom, sens_deg], savefig=savefig1[j], marker=['+', 'x'], xlabel='energy (TeV)', ylabel='E$^2$F sensitivity (erg/cm$^2$/s)', label=['full sens irf', 'degraded irf'], title=title, fontsize=12, show=False, tex=False)

            showSensitivity([10**energy_nom, 10**energy_deg], [flux_nom, flux_deg], savefig=savefig2[j], marker=['+', 'x'], xlabel='energy (TeV)', ylabel='ph flux (ph/cm$^2$/s)', label=['full sens irf', 'degraded irf'], title=title, fontsize=12, show=False, tex=False)

    showSensitivity([np.array(texp)], [np.array(fluxes_nom)], savefig=p.getPngDir() + 't5t20_sensInt_vsTime.png', marker=['+', 'x'], xlabel='time (s)', ylabel='ph flux (ph/cm$^2$/s)', label=['full sens irf', 'degraded irf'], title='', fontsize=12, ratio=False, xscale='log', show=True, tex=False)

