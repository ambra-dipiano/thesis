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
from scipy.interpolate import interp1d
import untangle
import csv
import re
import subprocess

def xmlConfig(cfg_file) :
  #Load configuration file
  cfgFile = os.path.dirname(__file__)+str(cfg_file)
  with open(cfgFile) as fd:
    cfg = untangle.parse(fd.read())
  return cfg.config

# --------------------------------- CLASS xml CONFIGURATION --------------------------------- !!!

class cfgMng_xml() :
  def __init__(self, cfg) :
    self.__initPath(cfg)

  def __initPath(self, cfg) :
    self.__cfg = cfg
    self.__workdir = self.__cfg.dir.workdir['path']
    self.__runpath = self.__cfg.dir.runpath['path']
    self.__datapath = self.__cfg.dir.datapath['path']
    self.__simpath = self.__cfg.dir.simpath['path']
    self.__selectpath = self.__cfg.dir.selectpath['path']
    self.__detpath = self.__cfg.dir.detpath['path']
    self.__csvpath = self.__cfg.dir.csvpath['path']

  def getWorkingDir(self) :
    return self.__workdir
  def setWorkingDir(self, workingDir) :
    self.__workdir = workingDir

  def getRunDir(self) :
    return self.__runpath.replace('${workdir}', self.getWorkingDir())
  def setRunDir(self, runDir) :
    self.__runpath = runDir

  def getDataDir(self) :
    return self.__datapath.replace('${runpath}', self.getRunDir())
  def setDataDir(self, dataDir) :
    self.__runpath = dataDir

  def getSimDir(self) :
    return self.__simpath.replace('${runpath}', self.getRunDir())
  def setSimDir(self, simDir) :
    self.__runpath = simDir

  def getSelectDir(self) :
    return self.__selectpath.replace('${runpath}', self.getRunDir())
  def setSelectDir(self, selectDir) :
    self.__runpath = selectDir

  def getDetDir(self) :
    return self.__detpath.replace('${runpath}', self.getRunDir())
  def setDetDir(self, detDir) :
    self.__runpath = detDir

  def getCsvDir(self) :
    return self.__csvpath.replace('${runpath}', self.getRunDir())
  def setCsvDir(self, csvDir) :
    self.__runpath = csvDir

# --------------------------------- CLASS ANALYSIS --------------------------------- !!!

class analysis() :
  def __init__(self, cfg_file):
    global p, CTOOLS
    CTOOLS = os.environ.get('CTOOLS')
    # conf ---!
    self.__cfg = xmlConfig(cfg_file)
    p = cfgMng_xml(self.__cfg)
    self.__pathout = p.getDataDir()
    self.seed = 1
    # files ---!
    self.model, self.template, self.table = (str() for i in range(3))
    self.event, self.event_list, self.event_selected, self.skymap = (str() for i in range(4))
    self.detectionXml, self.detectionReg, self.likeXml = (str() for i in range(3))
    self.sensCsv = str()
    self.caldb = 'prod2'
    self.irf = 'South_0.5h'
    self.irf_nominal =  CTOOLS + '/share/caldb/data/cta/%s/bcf/%s/irf_file.fits' % (self.caldb, self.irf)
    self.irf_degraded = self.irf_nominal.replace('prod', 'degr')
    # condition control ---!
    self.if_ebl = True
    self.extract_spec = False
    self.plot = False
    self.zfetch = False
    self.debug = False
    # data ---!
    self.z = 0.1
    self.t = [0, 2000]
    self.tmax = [1800]
    self.e = [0.03, 150.0]
    self.roi = 5
    self.pointing = [83.63, 22.01]
    self.sigma = 5
    self.maxSrc = 10
    # algorithm miscellaneous ---!
    self.z_ind = 1
    self.coord_sys = 'CEL'
    self.sky_subtraction = 'IRF'
    self.bkgType = 'irf'
    self.exclrad = 0.5
    self.corr_rad = 0.1
    self.src_type = 'POINT'
    self.corr_kern = 'GAUSSIAN'
    self.srcName = 'Src001'
    self.sgmrange = [0, 10]
    self.confidence = 0.95
    self.eref = 1
    # irf degradation ---!
    self.factor = 2
    self.sensType = 'Differential'

  def __openFITS(self):
    global hdul
    hdul = fits.open(self.template)
    return hdul
  def __closeFITS(self, hdul):
    hdul.close()

  @property # to be experimented ---!
  def __getFitsData(self):
    global energy, time, spectra
    hdul = self.__openFITS()
    energy = np.array(hdul[1].data)
    time = np.array(hdul[2].data)
    spectra = np.array(hdul[3].data)
    if len(hdul) == 5 :
      global ebl
      ebl = np.array(hdul[4].data)
      self.if_ebl = True
      self.__closeFITS(hdul)
      return energy, time, spectra, ebl
    else :
      self.if_ebl = False
      self.__closeFITS(hdul)
      return energy, time, spectra

  def __openCSV(self):
    df = pd.read_csv(self.table)
    df.dropna()
    return df

  def __getCsvData(self):
    df = self.__openCSV()
    cols = list(df.columns)
    tau_gilmore = np.array(df[cols[self.z_ind]])
    E = np.array(df[cols[0]]) / 1e3  # MeV --> GeV ---!
    return tau_gilmore, E

  def getTimeSlices(self):
    global time
    df = self.__openCSV()
    cols = list(df.columns)
    time = np.append(0, np.array(df[cols[1]]))
    for i in range(len(time)):
      if time[i] > max(self.tmax):
        time[i] = max(self.tmax)
        sliceObj = slice(0, i+1)
        break

    return time[sliceObj]

  def __add_ebl(self):
    energy, time, spectra = self.__getFitsData
    tau_gilmore, E = self.__getCsvData()
    # interpolate ---!
    interp = interp1d(E, tau_gilmore)
    global tau, ebl
    tau = np.array(interp(energy))
    ebl = np.empty_like(spectra)
    # compute ---!
    for i in range(len(time)):
      for j in range(len(energy)):
        ebl[i][j] = spectra[i][j] * np.exp(-tau[j])

    if self.plot is True:
      return E, tau_gilmore, energy, tau
    else:
      return

  def __zfetch(self):
    hdul = self.__openFITS()
    # fetch z and chose the table column with min distance from it ---!
    z = hdul[0].header['REDSHIFT']
    with open(self.table, 'r') as f:
      reader = csv.reader(f)
      hdr = next(reader)
    zlist = []
    for el in hdr:
      zlist.append(re.sub('[^0-9,.]', '', el))
    zlist.remove('')
    zlist = [float(i) for i in zlist]
    self.z = min(zlist, key=lambda x:abs(x-z))
    self.z_ind = zlist.index(self.z) +1
    return

  def fits_ebl(self, template_ebl):
    global x, y, x2, y2
    hdul = self.__openFITS() # open template ---!
    if self.zfetch is True:
      self.__zfetch()
    if self.plot is True:
      x, y, x2, y2 = self.__add_ebl()
    else:
      self.__add_ebl()
    # update fits ---!
    hdu = fits.BinTableHDU(name='EBL Gilmore', data=ebl)
    header = hdu.header
    header.set('UNITS', 'ph/cm2/s/GeV', ' ')
    hdu = fits.BinTableHDU(name='EBL Gilmore', data=ebl, header=header)
    hdul.append(hdu)
    # save to new ---!
    if os.path.isfile(template_ebl):
      os.system('rm ' + template_ebl)
    hdul.writeto(template_ebl, overwrite=False) # write new template ---!
    self.__closeFITS(hdul) # close template ---!
    if self.plot is True:
      return x, y, x2, y2
    else:
      return

  def __extractSpc(self):
    # time slices table ---!
    table = self.__pathout + 'time_slices.csv'
    if os.path.isfile(table):
      os.system('rm ' + table)
    if not os.path.isfile(table):
      os.system('touch ' + table)
      with open(table, 'a') as tab:
        tab.write('#bin,tmax_bin')

    # spectra and models ---!
    for i in range(Nt):
      if self.if_ebl is False:
        filename = self.__pathout + 'spec_tbin%02.out' % i
      else:
        filename = self.__pathout + 'spec_ebl_tbin%02d.out' % i

      if os.path.isfile(filename):
        os.system('rm ' + filename)
      if not os.path.isfile(filename):
        os.system('touch ' + filename)

      # time slices table ---!
      with open(table, 'a') as tab:
        tab.write('\n' + str(i) + ', ' + str(time[i][0]))

      # ebl ---!
      if self.if_ebl is True:
        with open(filename, 'a') as f:
          for j in range(Ne):
            # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
            if ebl is not None:
              f.write(str(energy[j][0] * 1000) + ' ' + str(ebl[i][j] / 1000) + "\n")
        # write bin models ---!
        os.system('cp ' + str(self.model) + ' ' + str(self.__pathout) + 'run0406_ID000126_ebl_tbin%02d.xml' % i)
        s = open(self.__pathout + 'run0406_ID000126_ebl_tbin%02d.xml' % i).read()
        s = s.replace('data/spec', 'spec_ebl_tbin%02d' % i)
        with open(self.__pathout + 'run0406_ID000126_ebl_tbin%02d.xml' % i, 'w') as f:
          f.write(s)
      # no ebl ---!
      else:
        with open(filename, 'a') as f:
          for j in range(Ne):
            # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
            f.write(str(energy[j][0] * 1000.0) + ' ' + str(spectra[i][j] / 1000.0) + "\n")
        # write bin models ---!
        os.system('cp ' + str(self.model) + ' ' + str(self.__pathout) + 'run0406_ID000126_tbin%02d.xml' % i)
        s = open(self.__pathout + 'run0406_ID000126_tbin%02d.xml' % i).read()
        s = s.replace('spec', 'spec_tbin%02d' % i)
        with open(self.__pathout + 'run0406_ID000126_tbin%02d.xml' % i, 'w') as f:
          f.write(s)

    return

  def load_template(self) :
    global energy, time, spectra, ebl
    if self.if_ebl is True :
      energy, time, spectra, ebl = self.__getFitsData
    if self.if_ebl is False:
      energy, time, spectra = self.__getFitsData

    global Nt, Ne, tbin_stop
    Nt = len(time)
    Ne = len(energy)

    # time grid ---!
    t = [0.0 for x in range(Nt + 1)]
    for i in range(Nt - 1):
      t[i + 1] = time[i][0] + (time[i + 1][0] - time[i][0]) / 2
    # tmax in last bin ---!
    t[Nt] = time[Nt - 1][0] + (time[Nt - 1][0] - t[Nt - 1])

    # stop the second after higher tmax ---!
    if self.tmax != None :
      tbin_stop = 1
      for bin in range(len(t)) :
        if t[bin] <= max(self.tmax) :
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

    if self.extract_spec is True:
      self.__extractSpc()

    return tbin_stop

  def eventSim(self) :
    sim = ctools.ctobssim()
    sim["inmodel"] = self.model
    sim["outevents"] = self.event
    sim["caldb"] = self.caldb
    sim["irf"] = self.irf
    sim["ra"] = self.pointing[0]
    sim["dec"] = self.pointing[1]
    sim["rad"] = self.roi
    sim["tmin"] = self.t[0]
    sim["tmax"] = self.t[1]
    sim["emin"] = self.e[0]
    sim["emax"] = self.e[1]
    sim["seed"] = self.seed
    sim["logfile"] = self.event.replace('.fits', '.log')
    sim["debug"] = self.debug
    sim.execute()

    return

  def obsList(self, obsname):
    # if '.fits' in str(self.event):
    #   self.event_list = 'obs_' + self.event.replace('.fits', '.xml')
    # else:
    #   self.event_list = 'obs_' + self.event
    xml = gammalib.GXml()
    obslist = xml.append('observation_list title="observation library"')

    for i in range(len(self.event)):
      obs = obslist.append('observation name="%s" id="%02d" instrument="CTA"' % (obsname, i))
      obs.append('parameter name="EventList" file="%s"' % self.event[i])
    xml.save(self.event_list)

    return

  def eventSelect(self, prefix):
    selection = ctools.ctselect()
    selection['inobs'] = self.event_list
    selection['outobs'] = self.event_selected
    selection['usepnt'] = True
    selection['prefix'] = prefix
    selection['rad'] = self.roi
    selection['tmin'] = self.t[0]
    selection['tmax'] = self.t[1]
    selection['emin'] = self.e[0]
    selection['emax'] = self.e[1]
    selection['logfile'] = self.event_selected.replace('.xml', '.log')
    selection['debug'] = self.debug
    selection.execute()

    return

  def eventSkymap(self, wbin=0.02):
    nbin = int(self.roi / wbin)
    skymap = ctools.ctskymap()
    skymap['inobs'] = self.event_selected
    skymap['outmap'] = self.skymap
    skymap['irf'] = self.irf
    skymap['caldb'] = self.caldb
    skymap['emin'] = self.e[0]
    skymap['emax'] = self.e[1]
    skymap['usepnt'] = True
    skymap['nxpix'] = nbin
    skymap['nypix'] = nbin
    skymap['binsz'] = wbin
    skymap['coordsys'] = self.coord_sys.upper()
    skymap['proj'] = 'CAR'
    skymap['bkgsubtract'] = self.sky_subtraction.upper()
    skymap['logfile'] = self.skymap.replace('.fits', '.log')
    skymap['debug'] = self.debug
    skymap.execute()

    return

  def runDetection(self) :
    self.detectionXml = '%s' % self.skymap.replace('_skymap.fits', '_det%ssgm.xml' % self.sigma)
    self.detectionReg = '%s' % self.skymap.replace('_skymap.fits', '_det%ssgm.reg' % self.sigma)

    detection = cscripts.cssrcdetect()
    detection['inmap'] = self.skymap
    detection['outmodel'] = self.detectionXml
    detection['outds9file'] = self.detectionReg
    detection['srcmodel'] = self.src_type.upper()
    detection['bkgmodel'] = self.bkgType.upper()
    detection['threshold'] = int(self.sigma)
    detection['maxsrcs'] = self.maxSrc
    detection['exclrad'] = self.exclrad
    detection['corr_rad'] = self.corr_rad
    detection['corr_kern'] = self.corr_kern.upper()
    detection['logfile'] = self.detectionXml.replace('.xml', '.log')
    detection['debug'] = self.debug
    detection.execute()

    return

  def maxLikelihood(self):
    like = ctools.ctlike()
    like['inobs'] = self.event_selected
    like['inmodel'] = self.detectionXml
    like['outmodel'] = self.likeXml
    like['caldb'] = self.caldb
    like['irf'] = self.irf
    like['fix_spat_for_ts'] = True
    like['logfile'] = self.likeXml.replace('.xml', '.log')
    like['debug'] = self.debug
    like.execute()

    return

  def confLevels(self, asym_errors):
    self.confidence_level=[0.6827, 0.9545, 0.9973]
    self.errors_conf = []
    for i in range(len(self.confidence_level)):
      self.errors_conf.append(asym_errors.replace('_errors', '_%2derr' % (self.confidence_level[i] * 100)))
      if not os.path.isfile(self.errors_conf[i]):
        err = ctools.cterror()
        err['inobs'] = self.event_selected
        err['inmodel'] = self.likeXml
        err['srcname'] = self.srcName
        err['outmodel'] = self.errors_conf[i]
        err['caldb'] = self.caldb
        err['irf'] = self.irf
        err['confidence'] = self.confidence_level[i]
        err['logfile'] = asym_errors.replace('.xml', '.log')
        err['debug'] = self.debug
        err.execute()

    return self.errors_conf

  def integrFlux(self):
    uplim = ctools.ctulimit()
    uplim['inobs'] = self.event_selected
    uplim['inmodel'] = self.likeXml
    uplim['srcname'] = self.srcName
    uplim['caldb'] = self.caldb
    uplim['irf'] = self.irf
    uplim['confidence'] = self.confidence
    uplim['sigma_min'] = self.sgmrange[0]
    uplim['sigma_max'] = self.sgmrange[1]
    uplim['eref'] = self.eref  # default reference energy for differential limit (in TeV)
    uplim['emin'] = self.e[0]  # default minimum energy for integral flux limit (in TeV)
    uplim['emax'] = self.e[1]  # default maximum energy for integral flux limit (in TeV)
    uplim['logfile'] = self.likeXml.replace('results.xml', 'flux.log')
    uplim['debug'] = self.debug
    uplim.execute()

    return

  def photonFlux_pl(self, gamma, k0, e0):
    e1 = self.e[0]*1e6
    e2 = self.e[1]*1e6
    delta = gamma + 1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

  def energyFlux_pl(self, gamma, k0, e0):
    k0 *= 1.60218e-6
    e0 *= 1.60218e-6
    e1 = self.e[0]*1e6 * 1.60218e-6
    e2 = self.e[1]*1e6 * 1.60218e-6
    delta = gamma+1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

  def degradeIRF(self):
    caldb_degr = self.caldb.replace('prod', 'degr')
    nominal_cal = CTOOLS + '/share/caldb/data/cta/' + self.caldb
    degraded_cal = CTOOLS + '/share/caldb/data/cta/' + caldb_degr
    # copy caldb if not ---!
    if not os.path.isdir(degraded_cal):
      os.system('cp -r %s %s' % (nominal_cal, degraded_cal))
    # permissions ---!
    if os.geteuid() == 0:
      print('!!! root')
      subprocess.run(['sudo', 'chmod', '-R', '777', degraded_cal], check=True)
    else:
      print('!!! not root')
      subprocess.run(['sudo', 'chmod', '-R', '777', degraded_cal], check=True)

    # degrade ---!
    extension = ['EFFECTIVE AREA', 'BACKGROUND']
    field = [4, 6]
    #  field = ['EFFAREA', 'BKG']
    inv = 1 / self.factor
    with fits.open(self.irf_nominal) as hdul:
      col = []
      for i in range(len(extension)):
        col.append(hdul[extension[i]].data.field(field[i])[:].astype(float))
      # np.array(col)

    a = np.where(np.array([i * inv for i in col[0]]) is np.nan, 0., np.array([i * inv for i in col[0]]))
    b = []
    for i in range(len(col[1][0])):
      b.append(np.where(col[1][0][i] is np.nan, 0., col[1][0][i]) * inv)

    b = np.array(b)
    tmp = [a, b]

    with fits.open(self.irf_degraded, mode='update') as hdul:
      for i in range(len(extension)):
        hdul[extension[i]].data.field(field[i])[:] = tmp[i]
        print('!!! DEGRADING IRF BY FACTOR %d !!!' %self.factor)
      # save changes ---!
      hdul.flush()
    # update and change permissions back ---!
    self.caldb.replace('prod', 'degr')
    # permissions ---!
    if os.geteuid() == 0:
      subprocess.run(['sudo', 'chmod', '-R', '755', degraded_cal], check=True)
    else:
      subprocess.run(['sudo', 'chmod', '-R', '755', degraded_cal], check=True)
    return

  def eventSens(self, bins=20, wbin=0.05):
    sens = cscripts.cssens()
    nbin = int(self.roi / wbin)
    sens['inobs'] = self.event
    sens['inmodel'] = self.likeXml
    sens['srcname'] = self.srcName
    sens['caldb'] = self.caldb
    sens['irf'] = self.irf
    sens['outfile'] = self.sensCsv
    sens['duration'] = self.t[1] - self.t[0]
    sens['rad'] = self.roi
    sens['emin'] = self.e[0]
    sens['emax'] = self.e[1]
    sens['bins'] = bins
    sens['sigma'] = self.sigma
    sens['type'] = self.sensType.capitalize()
    sens['npix'] = nbin
    sens['binsz'] = wbin
    sens['logfile'] = self.sensCsv.replace('', '.log')
    sens.execute()

    return