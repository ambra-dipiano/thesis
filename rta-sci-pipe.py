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

# ======================================================= #
# TESTING cssrcdetect FOR CTA-RTA PERFORMANCE EVALUATIONS #
# ======================================================= #

# IMPORTS ---!
from pkg_blindsearch import *
# from module_plot import *
import numpy as np
import csv
import os
import time

# --------------------------------- SETUP --------------------------------- !!!

# compact initialisation ---!
trials = 1  # trials
count = 0  # starting count
# cpus ---!
nthreads = 2
os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
os.environ['MKL_NUM_THREADS'] = str(nthreads)

# ctools/cscripts parameters ---!
caldb = 'prod3b-v2'
# caldb_degraded = caldb.replace('prod', 'degr')
irf = 'South_z40_0.5h'

sigma = 5  # detection acceptance (Gaussian)
texp = (10, 100)  # exposure times (s)
tdelay = 50  # slewing time (s)
tmax = []
for i in range(len(texp)):
    tmax.append(tdelay + texp[i])
ttotal = 600  # 1e4  # maximum tobs (4h at least) simulation total time (s)
add_hours = 10  # 7200  # +2h observation time added after first none detection (s)
run_duration = 600  # 1200  # 20min observation run time for LST in RTA (s) ---!
elow = 0.03  # simulation minimum energy (TeV)
ehigh = 150.0  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 150.0  # selection maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)
wbin = 0.02  # skymap bin width (deg)
corr_rad = 0.1  # Gaussian
max_src = 5  # max candidates
ts_threshold = 25  # TS threshold for reliable detection
highest_ts = ts_threshold  # minimum threshold to repoint (repoint will be for monotonically higher ts only)
reduce_flux = None  # flux will be devided by factor reduce_flux, if nominal then set to None

# conditions control ---!
checks1 = True  # prints info
checks2 = False  # prints more info
if_ebl = True  # uses the EBL absorbed template
if_cut = False  # adds a cut-off parameter to the source model
ebl_fits = False  # generate the EBL absorbed template
extract_spec = True  # generates spectral tables and obs definition models
irf_degrade = False  # use degraded irf
compute_degr = False  # compute irf degradation
src_sort = True  # sorts scandidates from highest TS to lowest
use_runs = True  # if True uses phlists of run_duration otherwise uese the template format
repoint = False  # repoint to source coords after positive detection
skip_exist = False  # skips the step if ID exists in csv
debug = False  # prints logfiles on terminal
if_log = True  # saves logfiles

# path configuration ---!
cfg = xmlConfig()
p = ConfigureXml(cfg)
# files ---!
targets = ['run0367_ID000227']
merger_maps = [p.getMergersDir() + target.replace('_ID', '_MergerID') + '_skymap.fits' for target in targets]

# recap and dof ---!
dof, m2, m1 = getDof()
print('!!! *** !!! dof = ', m2, ' - ', m1, ' = ', dof)
print('!!! *** !!! MODEL CUTOFF:', if_cut)
print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! nominal caldb:', caldb)
print('!!! *** !!! irf:', irf)
print('!!! *** !!! TS SORT:', src_sort)
print('!!! *** !!! selection energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi:', roi, ' (deg)')
print('!!! *** !!! blind detection confidence:', sigma, ' sigmas')
print('!!! *** !!! detection confidence ts threshold:', ts_threshold)
print('!!! *** !!! total observation time:', ttotal, ' s')
print('!!! *** !!! additional observation time:', add_hours, ' s')
print('!!! *** !!! delay time:', tdelay, ' s\n')
if use_runs:
    print('handle data in runs of', run_duration, 's\n')
else:
    print('handle data with template structure\n')

del dof, m1, m2

# --------------------------------- INITIALIZE --------------------------------- !!!

count += 1
for idx, target in enumerate(targets):
    p.setRunDir(target)
    print('------------------------------------------\n\nTarget:', target)
    if irf_degrade:
        phlist = p.getObsDir() + 'deg%06.fits' % count
    else:
        phlist = p.getObsDir() + 'ebl%06.fits' % count
    # pointing with off-axis equal to max prob GW ---!
    pointing = getPointingAlert(merger_map=merger_maps[idx])
    # pointing = [round(pointing[0], 3), round(pointing[1], 3)]
    print('Pointing coordinates:', pointing)

    # setup trials obj ---!
    observation = Analysis()
    observation.nthreads = nthreads
    observation.pointing = pointing
    observation.roi = roi
    # degrade IRF if required ---!
    observation.caldb = caldb
    observation.irf = irf
    if irf_degrade:
        if compute_degr:
            observation.degradeIrf()
        observation.caldb = caldb.replace('prod', 'degr')

    # --------------------------------- 1° LOOP :: trials  --------------------------------- !!!

    count += 1
    observation.seed = count
    clocking = tdelay - min(texp)  # simulate flowing time (subsequent temporal windows of 1s)
    GTIf = [run_duration for i in range(len(texp))]  # LST runs are of 20mins chunks ---!
    num = [1 for i in range(len(texp))]  # count on LST-like run chunks ---!
    print('\n\n!!! ************ STARTING TRIAL %d ************ !!!\n\n' % count) if checks1 else None
    print('!!! check ---- seed=', observation.seed) if checks2 else None
    # attach ID to fileroot ---!
    f = target + 'ebl%06d' % count
    #  if irf_degrade:
    #    f += 'irf'
    print('!!! check ---- obs:', target) if checks2 else None

    # --------------------------------- 2° LOOP :: tbins --------------------------------- !!!

    observation.e = [emin, emax]
    twindows = [int((ttotal - tdelay) / texp[i]) for i in
                range(len(texp))]  # number of temporal windows per exposure time in total time ---!
    tlast = [ttotal + tmax[i] for i in
             range(len(texp))]  # maximum observation time from last detection (not exceeding ttotal) ---!
    for i, t in enumerate(tlast):
        if t > ttotal:
            tlast[i] = ttotal
    is_detection = [True for i in
                    range(len(texp))]  # controls which avoid forwarding of tlast for subsequent non-detections ---!
    # looping through all light-curve time intervals ---!
    for j in range(int(max(twindows))):
        clocking += min(texp)  # passing time second by second ---!
        print(clocking, 'clock', tlast, is_detection) if checks1 else None

        # --------------------------------- CLOCKING BREAK --------------------------------- !!!

        # check tlast, if globally reached then stop current trial ---!
        if clocking >= max(tlast):
            print('end analysis trial', count, ' at clocking', tlast)
            break
        current_twindows = []

        # --------------------------------- 3° LOOP :: texp in tbin --------------------------------- !!!

        for i in range(len(texp)):
            if j == 0:
                current_twindows.append(texp[i])
            else:
                current_twindows.append(texp[i]) if clocking % texp[i] == 0 else None
        # looping for all the texp for which the tbin analysis needs to be computed ---!
        for i in range(len(current_twindows)):

            # --------------------------------- CLOCKING SKIP --------------------------------- !!!

            # check tlast, if locally reached then skip current bin ---!
            index = texp.index(current_twindows[i])
            if clocking > tlast[index]:
                print('skip analysis texp', texp[index]) if checks1 else None
                continue
            # --------------------------------- CHECK SKIP --------------------------------- !!!

            tbin = clocking / current_twindows[i]  # temporal bin number of this analysis
            IDbin = 'tbin%09d' % tbin
            csv_name = p.getCsvDir() + 'tesi_tdel%d_deg%s_%ds.csv' % (tdelay, str(irf_degrade), texp[index])
            if os.path.isfile(csv_name) and skip_exist:
                skip = checkTrialId(csv_name, IDbin)
            else:
                skip = False
            if skip_exist is True and skip is True:
                continue

            # --------------------------------- SELECTION TIME --------------------------------- !!!

            # if first tbin of tepx then don't add clocking time to selection edges ---!
            if clocking < tdelay:
                continue
            elif clocking == tdelay:
                observation.t = [tdelay, tmax[index]]
            elif clocking > tdelay and texp[index] == min(texp):
                observation.t = [clocking, texp[index] + clocking]
            elif clocking > tdelay and texp[index] != min(texp):
                observation.t = [tdelay + clocking, tmax[index] + clocking]
            if observation.t[1] > ttotal:
                observation.t[1] = ttotal

            # --------------------------------- OBSERVATION LIST --------------------------------- !!!

            event_list = p.getObsDir() + 'obs_' + target + '.xml'
            if os.path.isfile(event_list):
                os.remove(event_list)
            observation.input = phlist
            observation.output = event_list
            observation.obsList(obsname='run0406_ID000126')
            print('phlist:', phlist) if checks2 else None
            print('observation list:', event_list) if checks2 else None

            # --------------------------------- SELECTION --------------------------------- !!!

            event_selected = event_list.replace(p.getObsDir(), p.getSelectDir()).replace('obs_', 'texp%ds_' % texp[i])
            prefix = p.getSelectDir() + 'texp%ds_' % texp[i]
            # select events ---!
            if os.path.isfile(event_selected):
                os.remove(event_selected)
            observation.input = event_list
            observation.output = event_selected
            observation.eventSelect(prefix=prefix)
            print('selection', observation.output) if checks2 else None

            # --------------------------------- SKYMAP --------------------------------- !!!

            skymap = event_selected.replace(p.getSelectDir(), p.getDetDir()).replace('.xml', '_skymap.fits')
            if os.path.isfile(skymap):
                os.remove(skymap)
            observation.input = event_selected
            observation.output = skymap
            observation.eventSkymap(wbin=wbin)
            print('skymap', observation.output) if checks2 else None

            # --------------------------------- DETECTION & MODELING --------------------------------- !!!

            observation.corr_rad = corr_rad
            observation.max_src = max_src
            detectionXml = skymap.replace('_skymap.fits', '_det%dsgm.xml' % sigma)
            if os.path.isfile(detectionXml):
                os.remove(detectionXml)
            observation.input = skymap
            observation.output = detectionXml
            observation.runDetection()
            # showSkymap(file=skymap, reg=detectionXml.replace('.xml', '.reg'), show=False)
            print('detection', observation.output) if checks2 else None
            deobservation = ManageXml(detectionXml)
            deobservation.sigma = sigma
            deobservation.if_cut = if_cut
            deobservation.modXml()
            deobservation.prmsFreeFix()

            # --------------------------------- CANDIDATES NUMBER --------------------------------- !!!

            pos = [deobservation.loadRaDec()]
            Ndet = len(pos[0][0])

            # --------------------------------- MAX LIKELIHOOD --------------------------------- !!!

            start_time = time.time() if checks1 else None
            likeXml = detectionXml.replace('_det%dsgm' % observation.sigma, '_like%dsgm' % observation.sigma)
            if os.path.isfile(likeXml):
                os.remove(likeXml)
            observation.input = event_selected
            observation.model = detectionXml
            observation.output = likeXml
            observation.maxLikelihood()
            if checks1:
                end_time = time.time() - start_time
                print(end_time, 's with texp=%d s' % texp[i])
                print('likelihood', observation.output)
            likeObj = ManageXml(likeXml)
            if src_sort and Ndet > 0:
                highest_ts_src = likeObj.sortSrcTs()[0]
                print('!!! check ---- highest TS: ', highest_ts_src) if checks1 else None
            else:
                highest_ts_src = None

            # --------------------------------- DETECTION RA & DEC --------------------------------- !!!

            pos, ra_det, dec_det = ([] for n in range(3))
            pos.append(deobservation.loadRaDec(highest=highest_ts_src))
            ra_det.append(pos[0][0][0]) if len(pos[0][0]) > 0 else ra_det.append(np.nan)
            dec_det.append(pos[0][1][0]) if len(pos[0][0]) > 0 else dec_det.append(np.nan)
            Ndet = len(pos[0][0])

            # --------------------------------- CLOSE DET XML --------------------------------- !!!

            deobservation.closeXml()
            del deobservation

            # --------------------------------- BEST FIT TSV --------------------------------- !!!

            ts_list, ts = ([] for n in range(2))
            ts_list.append(likeObj.loadTs()) if Ndet > 0 else ts_list.append([np.nan])

            # only first elem ---!
            ts.append(ts_list[0][0])
            print('ts:', ts[0]) if checks1 else None

            # --------------------------------- Nsrc FOR TSV THRESHOLD --------------------------------- !!!

            # count src with TS >= 9
            n = 0
            for k in range(len(ts_list[0])):
                if float(ts_list[0][k]) >= ts_threshold:
                    n += 1

            Nsrc = n

            # --------------------------------- REPOINTING ---------------------------------- !!!

            # if positive detection has been achieved, use cource coordinates not original pointing
            if repoint:
                if float(ts[0]) > highest_ts:
                    pointing = (ra_det[0], dec_det[0])
                    highest_ts = float(ts[0])
                    print('repointing to', pointing, 'with TS:', ts[0])

            # --------------------------------- +2h FROM LAST DETECTION --------------------------------- !!!

            # if no detection set False and defined end of observation ---!
            if (float(ts[0]) < ts_threshold or str(ts[0]) == 'nan') and is_detection[index]:
                is_detection[index] = False
                # add 2hrs of obs time ---!
                tlast[index] = observation.t[1] + add_hours  # +2h ---!
                print('+2h; tlast = ', tlast[index], ' with texp = ', texp[index], 'at clocking', clocking)
                # only 4hrs of simulation avialable, if tlast exceeds them then reset to ttotal ---!
                if tlast[index] > ttotal:
                    tlast[index] = ttotal
                    print('reset tlast = ', tlast[index], ' with texp = ', texp[index])
            # if later detection then reset True ---!
            elif float(ts[0]) >= ts_threshold and not is_detection[index]:
                is_detection[index] = True

            # --------------------------------- BEST FIT RA & DEC --------------------------------- !!!

            ra_list, ra_fit, dec_list, dec_fit = ([] for n in range(4))
            coord = likeObj.loadRaDec() if Ndet > 0 else None
            ra_list.append(coord[0]) if Ndet > 0 else ra_list.append([np.nan])
            dec_list.append(coord[1]) if Ndet > 0 else dec_list.append([np.nan])

            # only first elem ---!
            ra_fit.append(ra_list[0][0])
            dec_fit.append(dec_list[0][0])

            # --------------------------------- BEST FIT SPECTRAL --------------------------------- !!!

            pref_list, pref, index_list, index, pivot_list, pivot = ([] for n in range(6))
            pref_err_list, pref_err = ([] for n in range(2))
            likeObj.if_cut = if_cut
            spectral = likeObj.loadSpectral()
            index_list.append(spectral[0]) if Ndet > 0 else index_list.append([np.nan])
            pref_list.append(spectral[1]) if Ndet > 0 else pref_list.append([np.nan])
            pivot_list.append(spectral[2]) if Ndet > 0 else pivot_list.append([np.nan])
            error = likeObj.loadPrefError()
            pref_err_list.append(error) if Ndet > 0 else pref_err_list.append([np.nan])

            # only first elem ---!
            index.append(index_list[0][0])
            pref.append(pref_list[0][0])
            pivot.append(pivot_list[0][0])
            pref_err.append(pref_err_list[0][0])

            # eventually cutoff ---!
            if if_cut:
                cutoff_list, cutoff = ([] for n in range(2))
                cutoff_list.append(spectral[3]) if Ndet > 0 else cutoff_list.append([np.nan])
                cutoff.append(cutoff_list[0][0])

            # --------------------------------- INTEGRATED FLUX --------------------------------- !!!

            flux_ph, flux_ph_err = ([] for n in range(2))
            if Ndet > 0:
                flux_ph.append(observation.photonFluxPowerLaw(index[0], pref[0], pivot[0]))  # E (MeV)
                flux_ph_err.append(observation.photonFluxPowerLaw(index[0], pref_err[0], pivot[0]))  # E (MeV)
            else:
                flux_ph.append(np.nan)
                flux_ph_err.append(np.nan)

            # MISSING THE CUT-OFF OPTION ---!!!

            # --------------------------------- CLOSE LIKE XML --------------------------------- !!!

            likeObj.closeXml()
            del likeObj

            # --------------------------------- RESULTS TABLE (csv) --------------------------------- !!!

            header = '#tbin,tinit,tend,Ndet,Nsrc,RA_det,DEC_det,RA_fit,DEC_fit,flux_ph,flux_ph_err,TS\n'
            ID = 'ID%06d' % count
            IDbin = 'tbin%09d' % tbin

            row = []
            if checks1:
                print('!!! ---------- check trial:', count)
                print('!!! ----- check texp:', texp[i], 's between: [', observation.t[0], ', ', observation.t[1],
                      ' ] s')
                print('!!! *** check Ndet:', Ndet)
                print('!!! *** check Nsrc:', Nsrc)
                print('!!! *** check ra_det:', ra_det[0])
                print('!!! *** check dec_det:', dec_det[0])
                print('!!! *** check ra_fit:', ra_fit[0])
                print('!!! *** check dec_fit:', dec_fit[0])
                print('!!! *** check flux_ph:', flux_ph[0])
                print('!!! *** check flux_ph_err:', flux_ph_err[0])
                print('!!! *** check ts:', ts[0])
                print('!!! *** ---------------------------')

            row.append(
                [IDbin, observation.t[0], observation.t[1], Ndet, Nsrc, ra_det[0], dec_det[0], ra_fit[0], dec_fit[0],
                 flux_ph[0], flux_ph_err[0], ts[0]])
            if os.path.isfile(csv_name):
                with open(csv_name, 'a') as csv_file:
                    w = csv.writer(csv_file)
                    w.writerows(row)
                    csv_file.close()
            else:
                with open(csv_name, 'w+') as csv_file:
                    csv_file.write(header)
                    w = csv.writer(csv_file)
                    w.writerows(row)
                    csv_file.close()

            # --------------------------------- CLEAR SPACE --------------------------------- !!!

            os.system('rm ' + p.getSelectDir() + '*ebl%06d*' % count)
            os.system('rm ' + p.getDetDir() + '*ebl%06d*' % count)

    del observation
print('\n\n!!! ================== END ================== !!!\n\n')
