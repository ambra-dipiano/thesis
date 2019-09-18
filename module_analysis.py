# ============================ #
# MODULE OF ANALYSIS FUNCTIONS #
# ============================ #

import gammalib
import ctools
import cscripts
from astropy.io import fits
import numpy as np
import os.path
import pandas as pd
from scipy.interpolate import interp1d, interp2d


# EXTRACT SPECTRUM ---!
def extract_spectrum(model, Nt, Ne, tbin_stop, energy, spectra, ebl=None, tau=None, if_ebl=False, pathout=None) :
  print('work in progress')

  for i in range(tbin_stop):
    if if_ebl is False:
      filename = pathout + 'spec_tbin' + str(i) + '.out'
    else:
      filename = pathout + 'spec_ebl_tbin' + str(i) + '.out'

    if os.path.isfile(filename):
      os.system('rm ' + filename)
    if not os.path.isfile(filename):
      os.system('touch ' + filename)
      out_file = open(filename, 'a')
      # out_file.write("E[MeV],Flux[fotoni/MeV/cm^2/s]"+"\n")
      out_file.close()

  # ebl ---!
  if if_ebl == True:
    for i in range(Nt):
      outfile = pathout + 'spec_ebl_tbin' + str(i) + '.out'
      out_file = open(outfile, 'a')
      for j in range(Ne):
        # write spectral data in E [MeV] and I [ph/cm2/s/MeV]
        if ebl is not None and tau is None:
          out_file.write(str(energy[j][0] * 1000) + ' ' + str(ebl[i][j] / 1000) + "\n")
        if tau is not None and ebl is None :
          out_file.write(str(energy[j][0] * 1000.0) + ' ' + str((spectra[i][j] / 1000.0) * np.exp(-tau[j])) + "\n")
      out_file.close()

      os.system('cp ' + model + ' ' + pathout + 'run0406_ID000126_ebl_tbin' + str(i) + '.xml')
      s = open(pathout + 'run0406_ID000126_ebl_tbin' + str(i) + '.xml').read()
      s = s.replace('data/spec', 'spec_ebl_tbin' + str(i))
      f = open(pathout + 'run0406_ID000126_ebl_tbin' + str(i) + '.xml', 'w')
      f.write(s)
      f.close()

  # no ebl ---!
  else:
    for i in range(Nt):
      outfile = pathout + 'spec_tbin' + str(i) + '.out'
      out_file = open(outfile, 'a')
      for j in range(Ne):
        # write spectral data in E [MeV] and I [ph/cm2/s/MeV]
        out_file.write(str(energy[j][0] * 1000.0) + ' ' + str(spectra[i][j] / 1000.0) + "\n")
      out_file.close()

      os.system('cp ' + model + ' ' + pathout + 'run0406_ID000126_tbin' + str(i))
      s = open(pathout + 'run0406_ID000126_tbin' + str(i) + '.xml').read()
      s = s.replace('spec', 'spec_tbin' + str(i))
      f = open(pathout + 'run0406_ID000126_tbin' + str(i) + '.xml', 'w')
      f.write(s)
      f.close()

  return

# LOAD TEMPLATE ---!
def load_template(template, tmax, extract_spec=False, model=None, pathout=None) :
  # open template ---!
  hdul = fits.open(template)
  # energybins [GeV] ---!
  energy = np.array(hdul[1].data)
  # timebins [s] ---!
  time = np.array(hdul[2].data)
  # spectra ---!
  spectra = np.array(hdul[3].data)
  # ebl ---!
  if len(hdul) == 5 :
    ebl = np.array(hdul[4].data) 
    if_ebl = True
  else :
    if_ebl = False 

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
  if tmax != None :
    tbin_stop = 1
    for bin in range(len(t)) :
      if t[bin] <= max(tmax) :
        tbin_stop += 1
      else :
        continue
  else :
    tbin_stop = None

  # energy grid ---!
  en = [1.0 for x in range(Ne + 1)]
  for i in range(Ne - 1):
    en[i + 1] = energy[i][0] + (energy[i + 1][0] - energy[i][0]) / 2
  # Emax in last bin ---!
  en[Ne] = energy[Ne - 1][0] + (energy[Ne - 1][0] - en[Ne - 1])

  if extract_spec is True and if_ebl is True :
    extract_spectrum(model, Nt, Ne, tbin_stop, energy=energy, spectra=spectra,
                     ebl=ebl, if_ebl=if_ebl, pathout=pathout)
  if extract_spec is True and if_ebl is False :
    # here you should call the add_ebl and/or fits_ebl ---!!!!!!!!!!
    extract_spectrum(model, Nt, Ne, tbin_stop, energy=energy, spectra=spectra,
                     if_ebl=if_ebl, pathout=pathout)
  return t, tbin_stop

# ADD EBL TO TEMPLATE ---!
def add_ebl(table, z, time, energy, spectra, plot=False) :

  df = pd.read_csv(table)
  cols = list(df.columns)
  df.dropna()
  # interpolate ---!
  tau = np.array(df[z])
  E = np.array(df[cols[0]])/1e3  # MeV --> GeV
  interp = interp1d(E, tau)
  tau_gilmore = np.array(interp(energy))
  ebl_gilmore = np.empty_like(spectra)
  # compute ---!
  for i in range(len(time)):
    for j in range(len(energy)): # QUI ERRORE
      ebl_gilmore[i][j] = spectra[i][j] * np.exp(-tau_gilmore[j])

  if plot is True :
    return ebl_gilmore, E, energy, tau, tau_gilmore
  else :
    return ebl_gilmore

# CREATE EBL FITS MODEL [WIP] ---!
def fits_ebl(template, template_ebl, table, zfetch=True, z=None, plot=False) :
  with fits.open(template) as hdul:
    # energybins [GeV] ---!
    energy = np.array(hdul[1].data)
    # timebins [s] ---!
    time = np.array(hdul[2].data)
    # spectra ---!
    spectra = np.array(hdul[3].data)
    # redshift [must approx to chose the column] ---!
    if zfetch is True :
      z = hdul[0].header['REDSHIFT']
    # retrive the ebl ---!
    if plot is True :
      ebl, x, x2, y, y2 = add_ebl(table, z, time, energy, spectra, plot=plot)
    else :
      ebl = add_ebl(table, z, time, energy, spectra, plot=plot)
    # update fits ---!
    hdu = fits.BinTableHDU(name='EBL Gilmore', data=ebl)
    header = hdu.header
    header.set('UNITS', 'ph/cm2/s/GeV', ' ')
    hdu = fits.BinTableHDU(name='EBL Gilmore', data=ebl, header=header)
    hdul.append(hdu)
    # save to new ---!
    os.system('rm ' + template_ebl)
    hdul.writeto(template_ebl, overwrite=False)

  if plot is True :
    return x, y, x2, y2
  else :
    return

# SIMULATE EVENT LIST ---!
def simulate_event(model, event, t=[0, 2000], e=[0.03, 150.0], caldb='prod2', irf='South_0.5h', roi=5, pointing=[83.63, 22.01], seed=1) :
  sim = ctools.ctobssim()
  sim["inmodel"] = model  # xml
  sim["outevents"] = event  # fits
  sim["caldb"] = caldb
  sim["irf"] = irf
  sim["ra"] = pointing[0]
  sim["dec"] = pointing[1]
  sim["rad"] = roi
  sim["tmin"] = t[0]
  sim["tmax"] = t[1]
  sim["emin"] = e[0]
  sim["emax"] = e[1]
  sim["seed"] = seed
  sim["logfile"] = event.replace('.fits', '.log')
  if os.path.isfile(event) is False:
    sim.execute()
  else :
    pass

  return

# CREATE OBSERVATION LIST ---!
def observation_list(event, eventList, obsname) :
  # combine in observatiion list
  xml = gammalib.GXml()
  obslist = xml.append('observation_list title="observation library"')

  for i in range(len(event)):
    obs = obslist.append('observation name="%s" id="%02d" instrument="CTA"' % (obsname,i))
    obs.append('parameter name="EventList" file="%s"' % event[i])

  xml.save(eventList)

  return

# EVENT SELECTION ---!
def select_event(eventList, event_selected, prefix, t=[0, 2000], e=[0.1, 100.0], roi=5) :
  selection = ctools.ctselect()
  selection['inobs'] = eventList
  selection['outobs'] = event_selected
  selection['usepnt'] = bool('yes')
  selection['prefix'] = prefix
  selection['rad'] = roi
  selection['tmin'] = t[0]
  selection['tmax'] = t[1]
  selection['emin'] = e[0]
  selection['emax'] = e[1]
  selection['logfile'] = event_selected.replace('.xml', '.log')
  selection['debug'] = bool('no')
  if os.path.isfile(event_selected) is False:
    selection.execute()
  else :
    pass

  return

# SKYMAP ---!
def skymap_event(event_selected, sky, e=[0.1, 100.0], caldb='prod2', irf='South_0.5h', wbin=0.02, nbin=200) :
  skymap = ctools.ctskymap()
  skymap['inobs'] = event_selected
  skymap['outmap'] = sky
  skymap['irf'] = irf
  skymap['caldb'] = caldb
  skymap['emin'] = e[0]
  skymap['emax'] = e[1]
  skymap['usepnt'] = bool('yes')
  skymap['nxpix'] = nbin
  skymap['nypix'] = nbin
  skymap['binsz'] = wbin
  skymap['coordsys'] = 'CEL'
  skymap['proj'] = 'CAR'
  skymap['bkgsubtract'] = 'IRF'
  skymap['logfile'] = sky.replace('.fits', '.log')
  skymap['debug'] = bool('no')
  if os.path.isfile(sky) is False:
    skymap.execute()
  else :
    pass

  return

# DETECTION ---!
def runDetection(skymap, sigma=5, maxSrc=10, bkgType='irf', exclrad=0.5, corr_rad=0.1) :
  detectionXml = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.xml' % sigma)
  detectionReg = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.reg' % sigma)

  detection = cscripts.cssrcdetect()
  detection['inmap'] = skymap
  detection['outmodel'] = detectionXml
  detection['outds9file'] = detectionReg
  detection['srcmodel'] = 'POINT'
  detection['bkgmodel'] = bkgType.upper()
  detection['threshold'] = int(sigma)
  detection['maxsrcs'] = maxSrc
  detection['exclrad'] = exclrad
  detection['corr_rad'] = corr_rad
  detection['corr_kern'] = 'GAUSSIAN'
  detection['logfile'] = detectionXml.replace('.xml', '.log')
  detection['debug'] = bool('no')
  detection.execute()

  return detectionXml, detectionReg
#   def __runDetection(self, skymap) :
#     self.__detectionXml = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.xml' % self.sigma)
#     self.__detectionReg = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.reg' % self.sigma)
#
#     detection = cscripts.cssrcdetect()
#     detection['inmap'] = skymap
#     detection['outmodel'] = self.__detectionXml
#     detection['outds9file'] = self.__detectionReg
#     detection['srcmodel'] = 'POINT'
#     detection['bkgmodel'] = self.bkgType.upper()
#     detection['threshold'] = int(self.sigma)
#     detection['maxsrcs'] = self.maxSrc
#     detection['exclrad'] = self.exclrad
#     detection['corr_rad'] = self.corr_rad
#     detection['corr_kern'] = 'GAUSSIAN'
#     detection['logfile'] = self.__detectionXml.replace('.xml', '.log')
#     detection['debug'] = bool('no')
#     detection.execute()
#
#     return

# CTLIKE ---!
def max_likelihood(event_selected, detection_model, results, caldb='prod2', irf='South_0.5h') :
  like = ctools.ctlike()
  like['inobs'] = event_selected
  like['inmodel'] = detection_model
  like['outmodel'] = results
  like['caldb'] = caldb
  like['irf'] = irf
  like['fix_spat_for_ts'] = bool('no')
  like['logfile'] = results.replace('.xml', '.log')
  like['debug'] = bool('no')  # default
  if os.path.isfile(results) is False:
    like.execute()
  else :
    pass

  return

# ASYM ERRORS ---!
def confidence_lv(event_selected, results, asym_errors, srcname='Src001', caldb='prod2', irf='South_0.5h', confidence_level=[0.6827, 0.9545, 0.9973]) :
  errors_conf = []
  for i in range(len(confidence_level)) :
    errors_conf.append(asym_errors.replace('_errors', '_%derrors' %(confidence_level[i]*100)))
    if not os.path.isfile(errors_conf[i]) :
      err = ctools.cterror()
      err['inobs'] = event_selected
      err['inmodel'] = results
      err['srcname'] = srcname
      err['outmodel'] = errors_conf[i]
      err['caldb'] = caldb
      err['irf'] = irf
      err['confidence'] = confidence_level[i]
      err['logfile'] = asym_errors.replace('.xml', '.log')
      err['debug'] = bool('no')  # default
      if os.path.isfile(asym_errors) is False:
        err.execute()
      else:
        pass

  return errors_conf

# UPPERLIM FLUX ---!
def integrated_flux(event_selected, results, srcname='Src001', caldb='prod2', irf='South_0.5h', conf=0.95, eref=1, erange=[1, 100], sgmrange=[0,10]) :
  uplim = ctools.ctulimit()
  uplim['inobs'] = event_selected
  uplim['inmodel'] = results
  uplim['srcname'] = srcname
  uplim['caldb'] = caldb
  uplim['irf'] = irf
  uplim['confidence'] = conf  # default
  uplim['sigma_min'] = sgmrange[0]  # default
  uplim['sigma_max'] = sgmrange[1]  # default
  uplim['eref'] = eref  # default reference energy for differential limit (in TeV)
  uplim['emin'] = erange[0]  # default minimum energy for integral flux limit (in TeV)
  uplim['emax'] = erange[1]  # default maximum energy for integral flux limit (in TeV)
  uplim['logfile'] = results.replace('results.xml', 'flux.log')
  uplim['debug'] =  bool('no')  # default
  uplim.execute()

  return

# DEGRADE IRF ---!
def degrade_IRF(irf, degraded_irf, factor=3) :
  extension = ['EFFECTIVE AREA', 'BACKGROUND']
  field = [4, 6]
  #  field = ['EFFAREA', 'BKG']
  inv = 1 / factor
  with fits.open(irf) as hdul:
    col = []
    for i in range(len(extension)):
      np.array(col.append(hdul[extension[i]].data.field(field[i])[:].astype(float)))

  a = np.where(np.array([i*inv for i in col[0]]) is np.nan, 0., np.array([i*inv for i in col[0]]))
  b = []
  for i in range(len(col[1][0])) :
    b.append(np.where(col[1][0][i] is np.nan, 0., col[1][0][i]) * inv)

  b = np.array(b)
  tmp = [a, b]

  with fits.open(degraded_irf, mode='update') as hdul:
    for i in range(len(extension)):
      hdul[extension[i]].data.field(field[i])[:] = tmp[i]
    # save changes ---!
    hdul.flush

  return

# CTA SENSITIVITY ---!
def sensitivity(model, event, output, caldb='prod2', irf='South_0.5h', t=100, e=[0.03, 150.0], roi=5, srcName='Crab',
                sigma=5, bins=20, npix=200, binsz=0.05, type='Differential') :
  sens = cscripts.cssens()
  sens['inobs'] = event
  sens['inmodel'] = model
  sens['srcname'] = srcName
  sens['caldb'] = caldb
  sens['irf'] = irf
  sens['outfile'] = output
  sens['duration'] = t
  sens['rad'] = roi
  sens['emin'] = e[0]
  sens['emax'] = e[1]
  sens['bins'] = bins
  sens['sigma'] = sigma
  sens['type'] = type
  sens['npix'] = npix
  sens['binsz'] = binsz
  sens['logfile'] = output.replace('', '.log')
  if os.path.isfile(output) is False:
    sens.execute()
  else :
    pass

  return