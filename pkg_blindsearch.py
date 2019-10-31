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
from lxml import etree as ET
import random

def xmlConfig(cfg_file) :
  # load configuration file ---!
  cfgFile = os.path.dirname(__file__)+str(cfg_file)
  # create configuration dictionary ---!
  with open(cfgFile) as fd:
    cfg = untangle.parse(fd.read())
  return cfg.config

def getTrueCoords(fits_file):
  with fits.open(fits_file) as hdul:
    ra = hdul[0].header['RA']
    dec = hdul[0].header['DEC']
  return (ra, dec)

def getPointing(merge_map, fits_file, roi=5):
  if merge_map==None:
    offaxis = (roi*random.random(), roi*random.random())
  else:
    with fits.open(merge_map) as hdul:
      # search max prob coords WIP ---!
      offaxis = (0,0)
  true_coord = getTrueCoords(fits_file)
  pointing = (true_coord[0] + offaxis[0], true_coord[1] + offaxis[1])
  return true_coord, pointing, offaxis

# --------------------------------- CLASS xml CONFIGURATION --------------------------------- !!!

class ConfigureXml() :
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

  # working directory ---!
  def getWorkingDir(self) :
    return self.__workdir
  def setWorkingDir(self, workingDir) :
    self.__workdir = workingDir

  # directory containing runs ---!
  def getRunDir(self) :
    return self.__runpath.replace('${workdir}', self.getWorkingDir())
  def setRunDir(self, runDir) :
    self.__runpath = runDir

  # directory storing runs data ---!
  def getDataDir(self) :
    return self.__datapath.replace('${runpath}', self.getRunDir())
  def setDataDir(self, dataDir) :
    self.__runpath = dataDir

  # target directory for simulations ---!
  def getSimDir(self) :
    return self.__simpath.replace('${runpath}', self.getRunDir())
  def setSimDir(self, simDir) :
    self.__runpath = simDir

  # target directory for selections ---!
  def getSelectDir(self) :
    return self.__selectpath.replace('${runpath}', self.getRunDir())
  def setSelectDir(self, selectDir) :
    self.__runpath = selectDir

  # target directory for pipeline products ---!
  def getDetDir(self) :
    return self.__detpath.replace('${runpath}', self.getRunDir())
  def setDetDir(self, detDir) :
    self.__runpath = detDir

  # target directory for output tables ---!
  def getCsvDir(self) :
    return self.__csvpath.replace('${runpath}', self.getRunDir())
  def setCsvDir(self, csvDir) :
    self.__runpath = csvDir

# --------------------------------- CLASS ANALYSIS --------------------------------- !!!

class analysis() :
  '''
  This class contains a mixture of ctools wrapper and pipeline methods. The former are used to easily access and set
  ctools while the latter handles all the analysis necessities: from handling the EBL absorption to degrading the IRFs,
  from extracting spectra to reading the template time bins needed for the simulations. Each public field (self.field)
  can be set accordingly to the user needs from a python script, while private fields (self.__field) is for internal
  usage. Equally, public methods (methodName()) can be invoked within a python script once the class is instanced while
  private methods (__methodName()) should only be used within the class itself.
  '''
  def __init__(self, cfgFile):
    # location of ctools ---!
    self.__CTOOLS = os.environ.get('CTOOLS')
    # path initial configuration ---!
    self.__cfg = xmlConfig(cfgFile)
    self.__p = ConfigureXml(self.__cfg)
    # files fields ---!
    self.model, self.template, self.table, self.sensCsv = (str() for i in range(4))
    self.output, self.input = (str() for i in range(2))
    self.caldb = 'prod2'  # production name in calibration database ---!
    self.irf = 'South_0.5h'  # irf ID name ---!
    # condition control ---!
    self.if_ebl = True  # set/unset EBL absorption feature ---!
    self.extract_spec = False  # set/unset spectra extraction feature ---!
    self.plot = False  # option for retrieving plotting values ---!
    self.zfetch = False  # set/unset automatic fetching of redshift ---!
    self.debug = False  # set/unset debug mode for ctools ---!
    self.if_log = True  # set/unset logfiles for ctools ---!
    # data fields ---!
    self.t = [0, 2000]  # time range (s/MJD) ---!
    self.tmax = [1800]  # maximum exposure time needed (s) ---!
    self.e = [0.03, 150.0]  # energy range (TeV) ---!
    self.roi = 5  # region of indeterest (deg) ---!
    self.pointing = [83.63, 22.01]  # RA/DEC or GLON/GLAT (deg) ---!
    self.sigma = 5  # Gaussian significance (sigmas) ---!
    self.max_src = 10  # Max number of candidates to list during blind-detection ---!
    # ctools miscellaneous ---!
    self.seed = 1  # MC seed ---!
    self.coord_sys = 'CEL'  # coordinate system <CEL|GAL> ---!
    self.sky_subtraction = 'IRF'  # skymap subtraction type <NONE|IRF|RING> ---!
    self.bkgType = 'irf'  # background model <Irf|Aeff|Racc> ---!
    self.src_type = 'POINT'  # source model type ---!
    self.srcName = 'Src001'  # name of source of interest ---!
    self.exclrad = 0.5  # radius around candidate to exclude from further search ---!
    self.corr_kern = 'GAUSSIAN'  # smoothing type ---!
    self.corr_rad = 0.1  # radius for skymap smoothing ---!
    self.sgmrange = [0, 10]  # range of gaussian sigmas ---!
    self.confidence = 0.95  # confidence level (%) ---!
    self.eref = 1  # energy reference for flux computation ---!
    self.sensType = 'Differential'  # sensitivity type <Integral|Differential> ---!
    # ebl specifics ---!
    self.z = 0.1  # redshift value ---!
    self.z_ind = 1  # redshift value index ---!
    # irf degradation & flux reduction ---!
    self.factor = 2
    # fits extension array ---!
    self.__time, self.__energy, self.__spectra, self.__ebl = (float() for i in range(4))

  # open and close the FITS files ---!
  def __openFITS(self):
    hdul = fits.open(self.template)
    return hdul
  def __closeFITS(self, hdul):
    hdul.close()
    return

  # retrive FITS data ---!
  def __getFitsData(self):
    hdul = self.__openFITS()
    self.__energy = np.array(hdul[1].data)
    self.__time = np.array(hdul[2].data)
    self.__spectra = np.array(hdul[3].data)
    if self.if_ebl:
      try:
        self.__ebl = np.array(hdul[4].data)
      except:
        raise IndexError('Template extensions out of range. Unable to load EBL absorbed spectra.')

    self.__closeFITS(hdul)
    return

  # load csv tabl in pandas DataFrame and drop NaN values---!
  def __openCSV(self):
    df = pd.read_csv(self.table)
    df.dropna()
    return df

  # retrive csv data ---!
  def __getCsvData(self):
    df = self.__openCSV()
    cols = list(df.columns)
    tau_gilmore = np.array(df[cols[self.z_ind]])
    E = np.array(df[cols[0]]) / 1e3  # MeV --> GeV ---!
    return tau_gilmore, E

  # retrive csv temporal bin grid of the template in use and return the necessary slice ---!
  def getTimeSlices(self):
    df = self.__openCSV()
    cols = list(df.columns)
    self.__time = np.append(0, np.array(df[cols[1]]))
    for i in range(len(self.__time)):
      if self.__time[i] > max(self.tmax):
        self.__time[i] = max(self.tmax)
        sliceObj = slice(0, i+1)
        break
    return self.__time[sliceObj]

  # compute the EBL absorption ---!
  def __addEbl(self):
    self.__getFitsData()
    tau_gilmore, E = self.__getCsvData()
    # interpolate linearly ---!
    interp = interp1d(E, tau_gilmore)
    tau = np.array(interp(self.__energy))
    self.__ebl = np.empty_like(self.__spectra)
    # compute absorption ---!
    for i in range(len(self.__time)):
      for j in range(len(self.__energy)):
        self.__ebl[i][j] = self.__spectra[i][j] * np.exp(-tau[j])

    if self.plot:
      return E, tau_gilmore, self.__energy, tau
    else:
      return

  # retrive the redshift and evaluate which table column is the nearest, then access its index ---!
  def __zfetch(self):
    hdul = self.__openFITS()
    # fetch z and chose the table column with min distance from it ---!
    z = hdul[0].header['REDSHIFT']
    with open(self.table, 'r') as f:
      reader = csv.reader(f)
      hdr = next(reader)
    zlist = []
    # load only the redshift columns ---!
    for el in hdr:
      zlist.append(re.sub('[^0-9,.]', '', el))
    zlist.remove('')
    zlist = [float(i) for i in zlist]
    # find nearest ---!
    self.z = min(zlist, key=lambda x:abs(x-z))
    self.z_ind = zlist.index(self.z) +1
    return

  # add EBL extension to a FITS template ---!
  def fitsEbl(self, template_ebl):
    hdul = self.__openFITS()
    if self.zfetch:
      self.__zfetch()
    if self.plot:
      x, y, x2, y2 = self.__addEbl()
    else:
      self.__addEbl()
    # update fits ---!
    hdu = fits.BinTableHDU(name='EBL Gilmore', data=self.__ebl)
    header = hdu.header
    header.set('UNITS', 'ph/cm2/s/GeV', ' ')
    hdu = fits.BinTableHDU(name='EBL Gilmore', data=self.__ebl, header=header)
    hdul.append(hdu)
    # save to new ---!
    if os.path.isfile(template_ebl):
      os.remove(template_ebl)
    hdul.writeto(template_ebl, overwrite=True)
    self.__closeFITS(hdul)
    if self.plot:
      return x, y, x2, y2
    else:
      return

  # extract template spectra, create xml model files and time slices csv file ---!
  def __extractSpec(self):
    # time slices table ---!
    table = self.__p.getDataDir() + 'time_slices.csv'
    if os.path.isfile(table):
      os.remove(table)
    with open(table, 'w+') as tab:
      tab.write('#bin,tmax_bin')

    # spectra and models ---!
    for i in range(self.__Nt):
      if self.if_ebl:
        filename = self.__p.getDataDir() + 'spec_ebl_tbin%02d.out' % i
      else:
        filename = self.__p.getDataDir() + 'spec_tbin%02d.out' % i
      if os.path.isfile(filename):
        os.remove(filename)

      # time slices table ---!
      with open(table, 'a') as tab:
        tab.write('\n' + str(i) + ', ' + str(self.__time[i][0]))

      # ebl ---!
      if self.if_ebl:
        with open(filename, 'a+') as f:
          for j in range(self.__Ne):
            # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
            if self.__ebl is not None:
              f.write(str(self.__energy[j][0] * 1000) + ' ' + str(self.__ebl[i][j] / 1000) + "\n")
        # write bin models ---!
        os.system('cp ' + str(self.model) + ' ' + str(self.__p.getDataDir()) + 'run0406_ID000126_ebl_tbin%02d.xml' % i)
        s = open(self.__p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i).read()
        s = s.replace('data/spec', 'spec_ebl_tbin%02d' % i)
        with open(self.__p.getDataDir() + 'run0406_ID000126_ebl_tbin%02d.xml' % i, 'w') as f:
          f.write(s)
      # no ebl ---!
      else:
        with open(filename, 'a+') as f:
          for j in range(self.__Ne):
            # write spectral data in E [MeV] and I [ph/cm2/s/MeV] ---!
            f.write(str(self.__energy[j][0] * 1000.0) + ' ' + str(self.__spectra[i][j] / 1000.0) + "\n")
        # write bin models ---!
        os.system('cp ' + str(self.model) + ' ' + str(self.__p.getDataDir()) + 'run0406_ID000126_tbin%02d.xml' % i)
        s = open(self.__p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i).read()
        s = s.replace('data/spec', 'spec_tbin%02d' % i)
        with open(self.__p.getDataDir() + 'run0406_ID000126_tbin%02d.xml' % i, 'w') as f:
          f.write(s)
    return

  # read template and return tbin_stop containing necessary exposure time coverage ---!
  def loadTemplate(self) :
    self.__getFitsData()

    self.__Nt = len(self.__time)
    self.__Ne = len(self.__energy)

    # time grid ---!
    t = [0.0 for x in range(self.__Nt + 1)]
    for i in range(self.__Nt - 1):
      t[i + 1] = self.__time[i][0] + (self.__time[i + 1][0] - self.__time[i][0]) / 2
    # tmax in last bin ---!
    t[self.__Nt] = self.__time[self.__Nt - 1][0] + (self.__time[self.__Nt - 1][0] - t[self.__Nt - 1])

    # stop the second after higher tmax ---!
    if self.tmax != None :
      tbin_stop = 1
      for bin in range(len(t)) :
        if t[bin] <= max(self.tmax) :
          tbin_stop += 1
        else :
          continue
    else :
      raise ValueError('Maximum exposure time (tmax) is larger than the template temporal evolution.')

    # energy grid ---!
    en = [1.0 for x in range(self.__Ne + 1)]
    for i in range(self.__Ne - 1):
      en[i + 1] = self.__energy[i][0] + (self.__energy[i + 1][0] - self.__energy[i][0]) / 2
    # Emax in last bin ---!
    en[self.__Ne] = self.__energy[self.__Ne - 1][0] + (self.__energy[self.__Ne - 1][0] - en[self.__Ne - 1])

    # extract spectrum if required ---!
    if self.extract_spec:
      self.__extractSpec()
    return tbin_stop

  # ctobssim wrapper ---!
  def eventSim(self) :
    sim = ctools.ctobssim()
    sim["inmodel"] = self.model
    sim["outevents"] = self.output
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
    sim["logfile"] = self.output.replace('.fits', '.log')
    sim["debug"] = self.debug
    if self.if_log:
      sim.logFileOpen()
    sim.execute()
    return

  # create observation list with gammalib ---!
  def obsList(self, obsname):
    xml = gammalib.GXml()
    obslist = xml.append('observation_list title="observation library"')

    for i in range(len(self.input)):
      obs = obslist.append('observation name="%s" id="%02d" instrument="CTA"' % (obsname, i))
      obs.append('parameter name="EventList" file="%s"' % self.input[i])
    xml.save(self.output)
    return

  # ctselect wrapper ---!
  def eventSelect(self, prefix):
    selection = ctools.ctselect()
    selection['inobs'] = self.input
    selection['outobs'] = self.output
    selection['usepnt'] = True
    selection['prefix'] = prefix
    selection['rad'] = self.roi
    selection['tmin'] = self.t[0]
    selection['tmax'] = self.t[1]
    selection['emin'] = self.e[0]
    selection['emax'] = self.e[1]
    selection['logfile'] = self.output.replace('.xml', '.log')
    selection['debug'] = self.debug
    if self.if_log:
      selection.logFileOpen()
    selection.execute()
    return

  # ctskymap wrapper ---!
  def eventSkymap(self, wbin=0.02):
    nbin = int(self.roi / wbin)
    skymap = ctools.ctskymap()
    skymap['inobs'] = self.input
    skymap['outmap'] = self.output
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
    skymap['logfile'] = self.output.replace('.fits', '.log')
    skymap['debug'] = self.debug
    if self.if_log:
      skymap.logFileOpen()
    skymap.execute()
    return

  # cssrcdetect wrapper ---!
  def runDetection(self) :
    self.detectionXml = '%s' % self.input.replace('_skymap.fits', '_det%ssgm.xml' % self.sigma)
    self.detectionReg = '%s' % self.input.replace('_skymap.fits', '_det%ssgm.reg' % self.sigma)

    detection = cscripts.cssrcdetect()
    detection['inmap'] = self.input
    detection['outmodel'] = self.output
    detection['outds9file'] = self.detectionReg
    detection['srcmodel'] = self.src_type.upper()
    detection['bkgmodel'] = self.bkgType.upper()
    detection['threshold'] = int(self.sigma)
    detection['maxsrcs'] = self.max_src
    detection['exclrad'] = self.exclrad
    detection['corr_rad'] = self.corr_rad
    detection['corr_kern'] = self.corr_kern.upper()
    detection['logfile'] = self.output.replace('.xml', '.log')
    detection['debug'] = self.debug
    if self.if_log:
      detection.logFileOpen()
    detection.execute()
    return

  # ctlike wrapper ---!
  def maxLikelihood(self):
    like = ctools.ctlike()
    like['inobs'] = self.input
    like['inmodel'] = self.model
    like['outmodel'] = self.output
    like['caldb'] = self.caldb
    like['irf'] = self.irf
    like['refit'] = True
    like['max_iter'] = 500
    like['fix_spat_for_ts'] = False
    like['logfile'] = self.output.replace('.xml', '.log')
    like['debug'] = self.debug
    if self.if_log:
      like.logFileOpen()
    like.execute()
    return

  # cterror wrapper ---!
  def confLevels(self, asym_errors):
    self.confidence_level=[0.6827, 0.9545, 0.9973]
    self.output = []
    for i in range(len(self.confidence_level)):
      self.output.append(asym_errors.replace('_errors', '_%2derr' % (self.confidence_level[i] * 100)))
      if not os.path.isfile(self.output[i]):
        err = ctools.cterror()
        err['inobs'] = self.input
        err['inmodel'] = self.model
        err['srcname'] = self.srcName
        err['outmodel'] = self.output[i]
        err['caldb'] = self.caldb
        err['irf'] = self.irf
        err['confidence'] = self.confidence_level[i]
        err['logfile'] = self.output[i].replace('.xml', '.log')
        err['debug'] = self.debug
        if self.if_log:
          err.logFileOpen()
        err.execute()
    return self.output

  # ctulimit wrapper ---!
  def integrFlux(self):
    uplim = ctools.ctulimit()
    uplim['inobs'] = self.input
    uplim['inmodel'] = self.model
    uplim['srcname'] = self.srcName
    uplim['caldb'] = self.caldb
    uplim['irf'] = self.irf
    uplim['confidence'] = self.confidence
    uplim['sigma_min'] = self.sgmrange[0]
    uplim['sigma_max'] = self.sgmrange[1]
    uplim['eref'] = self.eref  # default reference energy for differential limit (in TeV)
    uplim['emin'] = self.e[0]  # default minimum energy for integral flux limit (in TeV)
    uplim['emax'] = self.e[1]  # default maximum energy for integral flux limit (in TeV)
    uplim['logfile'] = self.model.replace('results.xml', 'flux.log')
    uplim['debug'] = self.debug
    if self.if_log:
      uplim.logFileOpen()
    uplim.execute()
    return

  # compute integral photon flux for PL model ---!
  def photonFluxPowerLaw(self, gamma, k0, e0):
    e1 = self.e[0]*1e6
    e2 = self.e[1]*1e6
    delta = gamma + 1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

  # compute integral energy flux for PL model ---!
  def energyFluxPowerLaw(self, gamma, k0, e0):
    k0 *= 1.60218e-6
    e0 *= 1.60218e-6
    e1 = self.e[0]*1e6 * 1.60218e-6
    e2 = self.e[1]*1e6 * 1.60218e-6
    delta = gamma+1
    factor = k0 / (e0**gamma * delta)
    flux = factor * (e2**delta - e1**delta)
    return flux

  # degrade IRFs via Effective Area and/or Background ---!
  def degradeIrf(self, exts=1):
    self.irf_nominal =  self.__CTOOLS + '/share/caldb/data/cta/%s/bcf/%s/irf_file.fits' % (self.caldb, self.irf)
    self.irf_degraded = self.irf_nominal.replace('prod', 'degr')
    # only if not existing ---!
    if not os.path.isfile(self.irf_degraded):
      caldb_degr = self.caldb.replace('prod', 'degr')
      nominal_cal = self.__CTOOLS + '/share/caldb/data/cta/' + self.caldb
      degraded_cal = self.__CTOOLS + '/share/caldb/data/cta/' + caldb_degr
      # permissions check through allowed id list ---!
      # open all ---!
      if os.geteuid() in [0, 1126]:
        subprocess.run(['chmod', '-R', '777', self.__CTOOLS + '/share/caldb/data/cta/'], check=True)
      else:
        subprocess.run(['sudo', 'chmod', '-R', '777', self.__CTOOLS + '/share/caldb/data/cta/'], check=True)
      # create degr caldb path if not existing ---!
      if not os.path.isdir(degraded_cal):
        os.mkdir(degraded_cal)
      if not os.path.isfile(degraded_cal+'/caldb.indx'):
        os.system('cp %s/caldb.indx %s/caldb.indx' %(nominal_cal, degraded_cal))
      if not os.path.isdir(degraded_cal+'/bcf'):
        os.mkdir(degraded_cal+'/bcf')
      if not os.path.isdir(degraded_cal+'/bcf/'+self.irf):
        os.mkdir(degraded_cal+'/bcf/'+self.irf)
      if os.path.isfile(self.irf_degraded):
        os.system('rm %s' %self.irf_degraded)
      if not os.path.isfile(self.irf_degraded):
        os.system('cp %s %s' %(self.irf_nominal, self.irf_degraded))
      # permissions check through allowed id list ---!
      # close all and open only degraded caldb ---!
      if os.geteuid() in [0, 1126]:
        subprocess.run(['chmod', '-R', '755', self.__CTOOLS + '/share/caldb/data/cta/'], check=True)
        subprocess.run(['chmod', '-R', '777', degraded_cal], check=True)
      else:
        subprocess.run(['sudo', 'chmod', '-R', '755', self.__CTOOLS + '/share/caldb/data/cta/'], check=True)
        subprocess.run(['sudo', 'chmod', '-R', '777', degraded_cal], check=True)

      # degrade ---!
      extension = ['EFFECTIVE AREA', 'BACKGROUND']
      field = [4, 6]
      #  field = ['EFFAREA', 'BKG']
      inv = 1 / self.factor
      with fits.open(self.irf_nominal) as hdul:
        col = []
        for i in range(exts):
          col.append(hdul[extension[i]].data.field(field[i])[:].astype(float))
      # effective area multiplied by inv of factor ---!
      a = np.where(np.array([i * inv for i in col[0]]) is np.nan, 0., np.array([i * inv for i in col[0]]))
      # background counts operate only where values is not NaN nor zero ---!
      b = []
      for i in range(len(col[1][0])):
        b.append(np.where(col[1][0][i] is np.nan, 0., col[1][0][i]) * inv) # left unchanged ---!

      b = np.array(b)
      tmp = [a, b]

      # save changes to degraded IRF ---!
      with fits.open(self.irf_degraded, mode='update') as hdul:
        for i in range(len(extension)):
          hdul[extension[i]].data.field(field[i])[:] = tmp[i]
        # save changes ---!
        hdul.flush()
      # update caldb ---!
      self.caldb.replace('prod', 'degr')
      # permissions check through allowed id list ---!
      # close degraded caldb ---!
      if os.geteuid() in [0, 1126]:
        subprocess.run(['chmod', '-R', '755', degraded_cal], check=True)
      else:
        subprocess.run(['sudo', 'chmod', '-R', '755', degraded_cal], check=True)
    return

  # cssens wrapper ---!
  def eventSens(self, bins=20, wbin=0.05):
    sens = cscripts.cssens()
    nbin = int(self.roi / wbin)
    sens['inobs'] = self.input
    sens['inmodel'] = self.model
    sens['srcname'] = self.srcName
    sens['caldb'] = self.caldb
    sens['irf'] = self.irf
    sens['outfile'] = self.output
    sens['duration'] = self.t[1] - self.t[0]
    sens['rad'] = self.roi
    sens['emin'] = self.e[0]
    sens['emax'] = self.e[1]
    sens['bins'] = bins
    sens['sigma'] = self.sigma
    sens['type'] = self.sensType.capitalize()
    sens['npix'] = nbin
    sens['binsz'] = wbin
    sens['logfile'] = self.output.replace('.csv', '.log')
    if self.if_log:
      sens.logFileOpen()
    sens.execute()
    return

  # reduce flux of template by given factor ---!
  def __reduceFluxSpec(self):
    spec_files = []
    # r: root, d: directories, f: files ---!
    for r, d, f in os.walk(self.__p.getDataDir()):
      for file in f:
        if self.if_ebl:
          if '.out' in file and 'ebl' in file and 'flux' not in file:
            spec_files.append(os.path.join(r, file))
        else:
          if '.out' in file and 'ebl' not in file and 'flux' not in file:
            spec_files.append(os.path.join(r, file))

    spec_files.sort()
    # new files with relative suffix ---!
    for i in range(len(spec_files)):
      if self.if_ebl:
        new_file = spec_files[i].replace('spec_ebl_tbin', 'spec_ebl_flux%d_tbin' %self.factor)
      else:
        new_file = spec_files[i].replace('spec_tbin', 'spec_flux%d_tbin' %self.factor)
      if os.path.isfile(new_file):
        os.remove(new_file)
      # modify by a given factor ---!
      with open(spec_files[i], 'r') as input, open(new_file, 'w+') as output:
        df = pd.read_csv(input, sep=' ', header=None)
        df.iloc[:,1] = df.iloc[:,1].apply(lambda x: float(x)/self.factor)
        df.to_csv(path_or_buf=output, sep=' ', index=False, header=None)
    return

  # replace path/to/spectrum/file.out in the xml model file ---!
  def __replaceSpecFile(self):
    xml_files = []
    # r: root, d: directories, f: files ---!
    for r, d, f in os.walk(self.__p.getDataDir()):
      for file in f:
        if self.if_ebl:
          if '.xml' in file and 'ebl' in file and 'flux' not in file:
            xml_files.append(os.path.join(r, file))
        else:
          if '.xml' in file and 'ebl' not in file and 'flux' not in file:
            xml_files.append(os.path.join(r, file))

    xml_files.sort()
    # replace ---!
    for i in range(len(xml_files)):
      if self.if_ebl:
        new_file = xml_files[i].replace('ID000126_ebl_tbin', 'ID000126_ebl_flux%d_tbin' %self.factor)
      else:
        new_file = xml_files[i].replace('ID000126_tbin', 'ID000126_flux%d_tbin' %self.factor)
      if os.path.isfile(new_file):
        os.remove(new_file)
      with open(xml_files[i], 'r') as input, open(new_file, 'w+') as output:
        content = input.read()
        if self.if_ebl:
          content = content.replace('spec_ebl_tbin', 'spec_ebl_flux%d_tbin' %self.factor)
        else:
          content = content.replace('spec_tbin', 'spec_flux%d_tbin' %self.factor)
        output.write(content)
    return

  # execute the flux reduction and consequent substitution of files ---!
  def makeFainter(self):
    self.__reduceFluxSpec()
    self.__replaceSpecFile()
    return

# --------------------------------- CLASS xml HANDLING --------------------------------- !!!

class ManageXml():
  '''
  This class contains all the methods which read, generate and modify xml files, as needed for the analysis.
  They are not comprehensive of all the parameters one could want to access though more could be added according
  to necessity. In the future this class will also handle the pipeline configuration files. 
  '''
  def __init__(self, xml, cfgFile):
    self.__xml = xml
    self.__cfg = xmlConfig(cfgFile)
    p = ConfigureXml(self.__cfg)
    self.file = open(self.__xml)
    self.srcLib = ET.parse(self.file)
    self.root = self.srcLib.getroot()
    self.tsvList = []
    self.pos = []
    self.err = []
    self.spectral = []
    self.sigma = 5
    self.default_model = True
    self.instr = 'CTA'
    self.bkgType = 'Irf'
    self.srcAtt = []
    self.bkgAtt = []
    self.tscalc = True
    self.if_cut = False

  def __getSrcObj(self):
    src = self.root.findall('source')
    return src

  def __skipNode(self, cfg):
    '''
    :retun true for skip node
    '''
    src = self.__getSrcObj()
    if src.attrib[cfg.get('idAttribute')] in cfg.get('skip'):
      return True

    for filter in cfg.get('filters'):
      if src.attrib[filter.get('attribute')] == filter.get('value'):
        return True

    if len(cfg.get('selectors')) == 0:
      return False

    for select in cfg.get('selectors'):
      if src.attrib[select.get('attribute')] == select.get('value'):
        return False

    return True

  def loadTs(self, highest=None):
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          tsv = src.attrib['ts']
          self.tsvList.append(tsv)
      else:
        if src.attrib['name'] == highest:
          tsv = src.attrib['ts']
          self.tsvList.append(tsv)

    return self.tsvList

  def loadRaDec(self, highest=None):
    posRaList, posDecList = ([] for i in range(2))
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          posRaList.append(ra)
          posDecList.append(dec)
      else:
        if src.attrib['name'] == highest:
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          posRaList.append(ra)
          posDecList.append(dec)
    self.pos = [posRaList, posDecList]
    return self.pos

  def loadConfInt(self, highest=None):
    raList, decList = ([] for i in range(2))
    for src in self.root.findall('source'):
      if highest == None:
        if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          raList.append(ra)
          decList.append(dec)
      else:
        if src.attrib['name'] == highest:
          ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
          dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
          raList.append(ra)
          decList.append(dec)
    self.err = [raList, decList]
    return self.err

  def loadSpectral(self, highest=None):
    indexList, prefList, pivotList = ([] for i in range(3))
    if self.if_cut is True :
      cutoffList = []

    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        if highest == None:
          index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
              src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
          pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
              src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
          pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
              src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
          indexList.append(index)
          prefList.append(pref)
          pivotList.append(pivot)
          if self.if_cut is True :
            cutoff = float(src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['scale'])
            cutoffList.append(cutoff)
        else:
          if src.attrib['name'] == highest:
            index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
            pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
            pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
                src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
            indexList.append(index)
            prefList.append(pref)
            pivotList.append(pivot)
            if self.if_cut is True:
              cutoff = float(src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['value']) * float(
                  src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['scale'])
              cutoffList.append(cutoff)

    if self.if_cut is False:
      self.spectral = [indexList, prefList, pivotList]
    else:
      self.spectral = [indexList, prefList, pivotList, cutoffList]
    return self.spectral

  def __saveXml(self):
    self.srcLib.write(self.__xml, encoding="UTF-8", xml_declaration=True,
                      standalone=False, pretty_print=True)
    return

  def __setModel(self):
    if self.default_model is True:
      Att_Prefactor = {'name': 'Prefactor', 'scale': '1e-16', 'value': '5.7', 'min': '1e-07', 'max': '1000.0', 'free': '1'}
      Att_Index = {'name': 'Index', 'scale': '-1', 'value': '2.4', 'min': '0', 'max': '5.0', 'free': '1'}
      Att_PivotEn = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '0.3', 'min': '1e-07', 'max': '1000.0', 'free': '0'}
      Bkg_Prefactor = {'name': 'Prefactor', 'scale': '1', 'value': '1', 'min': '1e-03', 'max': '1e+3', 'free': '1'}
      Bkg_Index = {'name': 'Index', 'scale': '1', 'value': '0.0', 'min': '-5', 'max': '+5.0', 'free': '1'}
      Bkg_PivotEn = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '0.3', 'min': '0.01', 'max': '1000.0', 'free': '0'}

      self.srcAtt = [Att_Prefactor, Att_Index, Att_PivotEn]
      self.bkgAtt = [Bkg_Prefactor, Bkg_Index, Bkg_PivotEn]
      if self.if_cut is True:
        Att_CutOff = {'name': 'CutoffEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '1'}
        self.srcAtt.append(Att_CutOff)

      return self.srcAtt, self.bkgAtt
    else:
      pass

  def modXml(self, overwrite=True):
    self.__setModel()
    # source ---!
    i = 0
    for src in self.root.findall('source'):
      i += 1
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        src.set('tscalc', '1') if self.tscalc is True else None
        # remove spectral component ---!
        rm = src.find('spectrum')
        src.remove(rm)
        # new spectrum ---!
        if self.if_cut:
          spc = ET.SubElement(src, 'spectrum', attrib={'type': 'ExponentialCutoffPowerLaw'})
        else:
          spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
        spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
        spc.tail = '\n\t\t'.replace('\t', ' ' * 2)
        src.insert(0, spc)
        # new spectral params ---!
        for j in range(len(self.srcAtt)):
          prm = ET.SubElement(spc, 'parameter', attrib=self.srcAtt[j])
          if prm.attrib['name'] == 'Prefactor' and i > 1:
            prm.set('value', str(float(prm.attrib['value']) / 2 ** (i - 1)))
          prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < len(self.srcAtt) else '\n\t\t'.replace('\t', ' ' * 2)
          spc.insert(j, prm)
      # background ---!
      else:
        # set bkg attributes ---!
        src.set('instrument', '%s' % self.instr.upper()) if self.instr.capitalize() != 'None' else None
        if self.bkgType.capitalize() == 'Aeff' or self.bkgType.capitalize() == 'Irf':
          src.set('type', 'CTA%sBackground' % self.bkgType.capitalize())
        if self.bkgType.capitalize() == 'Racc':
          src.set('type', 'RadialAcceptance')
        # remove spectral component ---!
        rm = src.find('spectrum')
        src.remove(rm)
        # new bkg spectrum ---!
        spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
        spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
        spc.tail = '\n\t'.replace('\t', ' ' * 2)
        # new bkg params ---!
        for j in range(len(self.bkgAtt)):
          prm = ET.SubElement(spc, 'parameter', attrib=self.bkgAtt[j])
          prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < len(self.bkgAtt) else '\n\t\t'.replace('\t', ' ' * 2)

    # instead of override original xml, save to a new one with suffix "_mod" ---!
    if not overwrite:
      self.__xml = self.__xml.replace('.xml', '_mod.xml')
    self.__saveXml()
    return

  def prmsFreeFix(self):
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        for prm in src.findall('*/parameter'):
          if prm.attrib['name'] not in self.__cfg.xml.bkg.free:
            prm.set('free', '0')
        for free in self.__cfg.xml.src.free:
          src.find('*/parameter[@name="%s"]' % free['prm']).set('free', '1') if free['prm'] != None else None
      else:
        for prm in src.findall('*/parameter'):
          if prm.attrib['name'] not in self.__cfg.xml.bkg.free:
            prm.set('free', '0')
        for free in self.__cfg.xml.bkg.free:
          src.find('*/parameter[@name="%s"]' % free['prm']).set('free', '1') if free['prm'] != None else None

    self.__saveXml()
    return

  def sortSrcTs(self):
    src = self.root.findall("*[@ts]")
    self.root[:-1] = sorted(src, key=lambda el: (el.tag, el.attrib['ts']), reverse=True)
    from_highest = []
    for src in self.root.findall("*[@ts]"):
      from_highest.append(src.attrib['name'])
    self.__saveXml()
    return from_highest

  def closeXml(self):
    self.file.close()
    return