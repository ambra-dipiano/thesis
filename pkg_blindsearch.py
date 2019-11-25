# ============================ #
# MODULE OF ANALYSIS FUNCTIONS #
# ============================ #

import gammalib
import ctools
import cscripts
import pyregion
from astropy.io import fits
import numpy as np
import os.path
import pandas as pd
from matplotlib.colors import SymLogNorm
from scipy.interpolate import interp1d
import untangle
import csv
import re
import subprocess
from lxml import etree as ET
import matplotlib
matplotlib.rcsetup.validate_backend('agg')
import matplotlib.pyplot as plt
import seaborn as sns

def xmlConfig(cfgfile='/config.xml') :
  '''
  This function loads a configuration xml file from the same directory where the code itself is stored.
  :param cfgfile: str:
            relative path of the configuration xml file wrt the directory of this file, default = "/config.xml"
  :return: cfg.config: dict
            dictionary containing the paths information
  '''
  # load configuration file ---!
  cfgfile = os.path.dirname(__file__)+str(cfgfile)
  # create configuration dictionary ---!
  with open(cfgfile) as fd:
    cfg = untangle.parse(fd.read())
  return cfg.config

def getDof(cfgfile='/config.xml'):
  cfg = xmlConfig(cfgfile)
  if type(cfg.xml.src.free) is list:
    src = len(cfg.xml.src.free)
  else:
    src = len([cfg.xml.src.free])
  bkg = len(cfg.xml.bkg.free)
  dof = src
  return dof, bkg+src, bkg

def getTrueCoords(fits_file):
  '''
  This function extract source position information from a fits table file.
  :param fits_file: str
            absolute path of the fits table file
  :return: (ra, dec): tuple
            source coordinates RA/DEC
  '''
  with fits.open(fits_file) as hdul:
    ra = hdul[0].header['RA']
    dec = hdul[0].header['DEC']
  return (ra, dec)

def getPointing(merge_map, fits_file, roi=5):
  '''
  This function extract source position information from a fits table file, max probability coordinates from a
  probability fits map and compute the corresponding off axis. If None prob map is given the offaxis is randomly
  computed and the pointing coordinates accordingly derived.
  :param merge_map: str
            absolute path of the probability map fits file
  :param fits_file: str
            absolute path of the fits table file
  :param roi: scalar
            region of interest
  :return: true_coord: tuple
            source coordinates RA/DEC
  :return: pointing: tuple
            pointing coordinates RA/DEC
  :return: offaxis: tuple
            offaxis angle between pointing and position RA/DEC
  '''
  true_coord = getTrueCoords(fits_file)
  if merge_map==None:
    offaxis = (np.random.uniform(0,5,1), np.random.uniform(0,5,1))
    pointing = (true_coord[0] + offaxis[0], true_coord[1] + offaxis[1])
  else:
    with fits.open(merge_map) as hdul:
      # search max prob coords WIP ---!
      offaxis = (0,0)
      pointing = (true_coord[0] + offaxis[0], true_coord[1] + offaxis[1])
  return true_coord, pointing, offaxis

def checkTrialId(file, id):
  '''
  This function checks if a trial ID is already existing within a data file.
  :param file: data file (str)
  :param id: trial ID (str)
  :return: True id the ID exists, False if it doesn't
  '''
  with open(file=file) as f:
    df = pd.read_csv(f)
    cols = list(df.columns)
    ids = df[cols[0]]
  if id in list(ids):
    skip = True
  else:
    skip= False
  return skip

# --------------------------------- CLASS xml CONFIGURATION --------------------------------- !!!

class ConfigureXml() :
  '''
  This class handles the configuration of paths for the analysis.
  '''
  def __init__(self, cfg) :
    self.__initPath(cfg)

  def __initPath(self, cfg) :
    self.__cfg = cfg
    self.__root = self.__cfg.dir.root['path']
    self.__workdir = self.__cfg.dir.workdir['path']
    self.__runpath = self.__cfg.dir.runpath['path']
    self.__datapath = self.__cfg.dir.datapath['path']
    self.__simpath = self.__cfg.dir.simpath['path']
    self.__selectpath = self.__cfg.dir.selectpath['path']
    self.__detpath = self.__cfg.dir.detpath['path']
    self.__csvpath = self.__cfg.dir.csvpath['path']
    self.__pngpath = self.__cfg.dir.pngpath['path']

  def __checkDir(self, dir):
    isdir = os.path.isdir(dir)
    return isdir
  def __makeDir(self, dir):
    if not self.__checkDir(dir=dir):
      os.mkdir(dir)

  def getRootDir(self):
    return self.__root

  # working directory ---!
  def getWorkingDir(self):
    self.__makeDir(self.__workdir.replace('${root}', self.__root))
    return self.__workdir.replace('${root}', self.__root)
  def setWorkingDir(self, workingDir):
    self.__workdir = workingDir
    self.__makeDir(workingDir)

  # directory containing runs ---!
  def getRunDir(self):
    self.__makeDir(self.__runpath.replace('${workdir}', self.getWorkingDir()))
    return self.__runpath.replace('${workdir}', self.getWorkingDir())
  def setRunDir(self, runDir):
    self.__runpath = runDir
    self.__makeDir(runDir)

  # directory storing runs data ---!
  def getDataDir(self):
    self.__makeDir(self.__datapath.replace('${runpath}', self.getRunDir()))
    return self.__datapath.replace('${runpath}', self.getRunDir())
  def setDataDir(self, dataDir):
    self.__runpath = dataDir
    self.__makeDir(dataDir)

  # target directory for simulations ---!
  def getSimDir(self):
    self.__makeDir(self.__simpath.replace('${runpath}', self.getRunDir()))
    return self.__simpath.replace('${runpath}', self.getRunDir())
  def setSimDir(self, simDir):
    self.__runpath = simDir
    self.__makeDir(simDir)

  # target directory for selections ---!
  def getSelectDir(self):
    self.__makeDir(self.__selectpath.replace('${runpath}', self.getRunDir()))
    return self.__selectpath.replace('${runpath}', self.getRunDir())
  def setSelectDir(self, selectDir):
    self.__runpath = selectDir
    self.__makeDir(selectDir)

  # target directory for pipeline products ---!
  def getDetDir(self):
    self.__makeDir(self.__detpath.replace('${runpath}', self.getRunDir()))
    return self.__detpath.replace('${runpath}', self.getRunDir())
  def setDetDir(self, detDir):
    self.__runpath = detDir
    self.__makeDir(detDir)

  # target directory for output tables ---!
  def getCsvDir(self):
    self.__makeDir(self.__csvpath.replace('${runpath}', self.getRunDir()))
    return self.__csvpath.replace('${runpath}', self.getRunDir())
  def setCsvDir(self, csvDir):
    self.__runpath = csvDir
    self.__makeDir(csvDir)

  def getPngDir(self):
    self.__makeDir(self.__pngpath.replace('${workdir}', self.getWorkingDir()))
    return self.__pngpath.replace('${workdir}', self.getWorkingDir())
  def setPngDir(self, pngDir):
    self.__pngpath = pngDir
    self.__makeDir(pngDir)

# --------------------------------- CLASS ANALYSIS --------------------------------- !!!

class Analysis() :
  '''
  This class contains a mixture of ctools wrapper and pipeline methods. The former are used to easily access and set
  ctools while the latter handles all the analysis necessities: from handling the EBL absorption to degrading the IRFs,
  from extracting spectra to reading the template time bins needed for the simulations. Each public field (self.field)
  can be set accordingly to the user needs from a python script, while private fields (self.__field) is for internal
  usage. Equally, public methods (methodName()) can be invoked within a python script once the class is instanced while
  private methods (__methodName()) should only be used within the class itself.
  '''
  def __init__(self, cfgfile='/config.xml'):
    # location of ctools ---!
    self.__CTOOLS = os.environ.get('CTOOLS')
    # path initial configuration ---!
    self.__cfg = xmlConfig(cfgfile)
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
    self.t = [0, 1800]  # time range (s/MJD) ---!
    self.tmax = 1800  # maximum exposure time needed (s) ---!
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
    self.src_name = 'Src001'  # name of source of interest ---!
    self.exclrad = 0.5  # radius around candidate to exclude from further search ---!
    self.corr_kern = 'GAUSSIAN'  # smoothing type ---!
    self.corr_rad = 0.1  # radius for skymap smoothing ---!
    self.sgmrange = [0, 10]  # range of gaussian sigmas ---!
    self.confidence = 0.95  # confidence level (%) ---!
    self.eref = 1  # energy reference for flux computation ---!
    self.sensType = 'Differential'  # sensitivity type <Integral|Differential> ---!
    self.nthreads = 1
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
      if self.__time[i] > self.tmax or self.__time[i] == self.tmax:
        self.__time[i] = self.tmax
        bin = i+1
        break
    sliceObj = slice(0, bin)
    return self.__time[sliceObj]

  # compute the EBL absorption ---!
  def __addEbl(self):
    self.__getFitsData()
    tau_gilmore, E = self.__getCsvData()
    # interpolate linearly handling NaNs/inf/zeroes ---!
    with np.errstate(invalid='raise'):
      interp = interp1d(E, tau_gilmore, bounds_error=False)
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
      tbin_stop = 0
      for bin in range(len(t)) :
        if t[bin] <= self.tmax :
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
    sim["nthreads"] = self.nthreads
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
    # selection["nthreads"] = self.nthreads
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
    # skymap["nthreads"] = self.nthreads
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
    # detection["nthreads"] = self.nthreads
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
    like["nthreads"] = self.nthreads
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
        err['srcname'] = self.src_name
        err['outmodel'] = self.output[i]
        err['caldb'] = self.caldb
        err['irf'] = self.irf
        err['confidence'] = self.confidence_level[i]
        err["nthreads"] = self.nthreads
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
    uplim['srcname'] = self.src_name
    uplim['caldb'] = self.caldb
    uplim['irf'] = self.irf
    uplim['confidence'] = self.confidence
    uplim['sigma_min'] = self.sgmrange[0]
    uplim['sigma_max'] = self.sgmrange[1]
    uplim['eref'] = self.eref  # default reference energy for differential limit (in TeV)
    uplim['emin'] = self.e[0]  # default minimum energy for integral flux limit (in TeV)
    uplim['emax'] = self.e[1]  # default maximum energy for integral flux limit (in TeV)
    uplim["nthreads"] = self.nthreads
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

  # initialize paths for caldb degradation: directories and files ---!
  def __initCaldbIrf(self):
    nominal_irf =  self.__CTOOLS + '/share/caldb/data/cta/%s/bcf/%s/irf_file.fits' % (self.caldb, self.irf)
    degraded_irf = nominal_irf.replace('prod', 'degr')
    caldb_degr = self.caldb.replace('prod', 'degr')
    folder = self.__CTOOLS + '/share/caldb/data/cta/'
    nominal_cal =  folder + self.caldb
    degraded_cal = folder + caldb_degr
    return folder, nominal_cal, nominal_irf, degraded_cal, degraded_irf

  # updates the degraded caldb index by replacing all "prod" references with "degr" ---!
  def __updateCaldbIndex(self, index):
    # read content ---!
    with open(index, 'r', encoding="ISO-8859-1") as f:
      filedata = f.read()
    # Replace the target keyword ---!
    filedata = filedata.replace('prod', 'degr').replace('PROD', 'DEGR')
    # Write the file out again ---!
    with open(index, 'w', encoding="ISO-8859-1") as f:
      f.write(filedata)
    return

  # create copy of caldb and corresponding caldb.inx file ---!
  def __mockNominalCaldb(self, nominal_cal, nominal_irf, degraded_cal, degraded_irf):
    if not os.path.isdir(degraded_cal):
      os.mkdir(degraded_cal)
    if not os.path.isfile(degraded_cal+'/caldb.indx'):
      os.system('cp %s/caldb.indx %s/caldb.indx' %(nominal_cal, degraded_cal))
      # update caldb.indx file ---!
      self.__updateCaldbIndex(degraded_cal+'/caldb.indx')
    if not os.path.isdir(degraded_cal+'/bcf'):
      os.mkdir(degraded_cal+'/bcf')
    if not os.path.isdir(degraded_cal+'/bcf/'+self.irf):
      os.mkdir(degraded_cal+'/bcf/'+self.irf)
    if os.path.isfile(degraded_irf):
      os.system('rm %s' %degraded_irf)
    if not os.path.isfile(degraded_irf):
      os.system('cp %s %s' %(nominal_irf, degraded_irf))
    return

  # change permission to 777 and ask for password if user id not in idlist param ---!
  def __openPermission(self, path, idlist=(0,1126)):
    if os.geteuid() in idlist:
      subprocess.run(['chmod', '-R', '777', path], check=True)
    else:
      subprocess.run(['sudo', 'chmod', '-R', '777', path], check=True)
    return

  # change permission to 755 and ask for password if user id not in idlist param ---!
  def __closePermission(self, path, idlist=(0,1126)):
    if os.geteuid() in idlist:
      subprocess.run(['chmod', '-R', '755', path], check=True)
    else:
      subprocess.run(['sudo', 'chmod', '-R', '755', path], check=True)
    return

  # degrade Aff by self.factor (either scalar or array of energy-bins dimension) ---!
  def __degrAeff(self, nominal_irf, degraded_irf, r=False):
    # initialise ---!
    inv = 1 / self.factor
    extension = 'EFFECTIVE AREA'
    field = 4
    with fits.open(nominal_irf) as hdul:
      elo = np.array(hdul[extension].data.field(0)[:].astype(float)[0])
      ehi = np.array(hdul[extension].data.field(1)[:].astype(float)[0])
      e = elo + 0.5*(ehi - elo)
      tlo = np.array(hdul[extension].data.field(2)[:].astype(float)[0])
      thi = np.array(hdul[extension].data.field(3)[:].astype(float)[0])
      theta = tlo + 0.5*(thi - tlo)
      aeff = np.array(hdul[extension].data.field(field)[:].astype(float)[0])
    # effective area multiplied by inv of factor ---!
    a = np.where(np.array([i * inv for i in aeff]) is np.nan, 0., np.array([i * inv for i in aeff]))
    # degrade and save new ---!
    with fits.open(degraded_irf, mode='update') as hdul:
      hdul[extension].data.field(field)[:] = a
      # save changes ---!
      hdul.flush()
    # return only if bkg counts must be degraded ---!
    if not r:
      return
    else:
      return aeff, a, theta, e

  # degrade bkg counts by normalise for aeff nominal and multiply times aeff degraded ---!
  def __degrBkg(self, nominal_irf, degraded_irf):
    # degrade Aeff and get its returns ---!
    aeff_nom, aeff_deg, theta, e_aeff = self.__degrAeff(nominal_irf=nominal_irf, degraded_irf=degraded_irf, r=True)
    # initialise ---!
    extension = 'BACKGROUND'
    field = 6
    with fits.open(nominal_irf) as hdul:
      xlo = np.array(hdul[extension].data.field(0)[:].astype(float)[0])
      xhi = np.array(hdul[extension].data.field(1)[:].astype(float)[0])
      x = xlo + 0.5*(xhi - xlo)
      ylo = np.array(hdul[extension].data.field(2)[:].astype(float)[0])
      yhi = np.array(hdul[extension].data.field(3)[:].astype(float)[0])
      y = ylo + 0.5*(yhi - ylo)
      elo = np.array(hdul[extension].data.field(4)[:].astype(float)[0])
      ehi = np.array(hdul[extension].data.field(5)[:].astype(float)[0])
      e_bkg = elo + 0.5*(ehi - elo)
      bkg = np.array(hdul[extension].data.field(field)[:].astype(float)[0])
    # spatial pixel/deg conversion factor ---!
    conv_factor = (xhi.max() - xlo.min()) / theta.max()
    # interpolated Aeff via energy grid ---!
    nominal_interp, degraded_interp = ([[]*i for i in range(len(theta))] for i in range(2))
    for i in range(len(theta)):
      fnom = interp1d(e_aeff[:], aeff_nom[i,:])
      nominal_interp[i].append(fnom(e_bkg[:]))
      fdeg = interp1d(e_aeff[:], aeff_deg[i,:])
      degraded_interp[i].append(fdeg(e_bkg[:]))

    # flatten list of theta interpolations (theta array of energy frames) ---!
    nominal_interp = np.array([item for sublist in nominal_interp for item in sublist])
    degraded_interp = np.array([item for sublist in degraded_interp for item in sublist])
    # empty copy of bkg tensor ---!
    b = np.empty_like(bkg)
    for idf, frame in enumerate(bkg[:,0,0]):
      for idx, xpix in enumerate(bkg[idf,:,0]):
        for idy, ypix in enumerate(bkg[idf,idx,:]):
          # find radius in degrees ---!
          r = np.sqrt((0 - xpix)**2 + (0 - ypix)**2)
          rdegree = r * conv_factor
          # find corresponding theta index ---!
          angle = min(theta, key=lambda x:abs(x-rdegree))
          idtheta = np.where(np.isin(theta[:], angle))
          # degrade the background count for frame/x/y point ---!
          if nominal_interp[idtheta,idf] == 0.:
            b[idf, idx, idy] = 0.
          else:
            b[idf,idx,idy] = bkg[idf,idx,idy] / nominal_interp[idtheta,idf] * degraded_interp[idtheta,idf]

    # save to new ---!
    with fits.open(degraded_irf, mode='update') as hdul:
      hdul[extension].data.field(field)[:] = b
      # save changes ---!
      hdul.flush()
    return

  # degrade IRFs via Effective Area and/or Background ---!
  def degradeIrf(self, bkg=True):
    # initialize ---!
    folder, nominal_cal, nominal_irf, degraded_cal, degraded_irf = self.__initCaldbIrf()
    # open all folder permission ---!
    self.__openPermission(path=folder)
    # create degr caldb path if not existing ---!
    self.__mockNominalCaldb(nominal_cal=nominal_cal, nominal_irf=nominal_irf,
                            degraded_cal=degraded_cal, degraded_irf=degraded_irf)
    # close all folder permission and open only degraded caldb permission ---!
    self.__closePermission(path=folder)
    self.__openPermission(path=degraded_cal)
    # degradation aeff ---!
    if not bkg:
      self.__degrAeff(nominal_irf=nominal_irf, degraded_irf=degraded_irf)
    # degradation bkg counts ---!
    else:
      self.__degrBkg(nominal_irf=nominal_irf, degraded_irf=degraded_irf)

    # close degraded caldb permission ---!
    self.__closePermission(degraded_cal)
    # update caldb ---!
    self.caldb = self.caldb.replace('prod', 'degr')
    return

  # cssens wrapper ---!
  def eventSens(self, bins=20, wbin=0.05, enumbins=0):
    sens = cscripts.cssens()
    nbin = int(self.roi / wbin)
    sens['inobs'] = self.input
    sens['inmodel'] = self.model
    sens['srcname'] = self.src_name
    sens['caldb'] = self.caldb
    sens['irf'] = self.irf
    sens['outfile'] = self.output
    sens['duration'] = self.t[1] - self.t[0]
    sens['rad'] = self.roi
    sens['emin'] = self.e[0]
    sens['emax'] = self.e[1]
    sens['bins'] = bins
    if enumbins != 0:
      sens['enumbins'] = enumbins
      sens['npix'] = nbin
      sens['binsz'] = wbin
    sens['sigma'] = self.sigma
    sens['type'] = self.sensType.capitalize()
    sens["nthreads"] = self.nthreads
    sens['logfile'] = self.output.replace('.csv', '.log')
    sens['debug'] = self.debug
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

  # returns a random total delay time (slew time + gw latency) within given ranges ---!
  def totalDelay(self, slew=(0,50), gw_latency=(0,36000)):
    tslew = np.random.uniform(slew[0], slew[1], 1)
    tgw = np.random.uniform(gw_latency[0], gw_latency[1])
    delay = tslew + tgw
    return delay

# --------------------------------- CLASS xml HANDLING --------------------------------- !!!

class ManageXml():
  '''
  This class contains all the methods which read, generate and modify xml files, as needed for the analysis.
  They are not comprehensive of all the parameters one could want to access though more could be added according
  to necessity. In the future this class will also handle the pipeline configuration files. Each public field (self.field)
  can be set accordingly to the user needs from a python script, while private fields (self.__field) is for internal
  usage. Equally, public methods (methodName()) can be invoked within a python script once the class is instanced while
  private methods (__methodName()) should only be used within the class itself.
  '''
  def __init__(self, xml, cfgfile='/config.xml'):
    self.__xml = xml
    self.__cfg = xmlConfig(cfgfile)
    self.__p = ConfigureXml(self.__cfg)
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

  def setTsTrue(self):
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
        src.set('tscalc', '1')
    self.__saveXml()
    return

  def modXml(self, overwrite=True):
    self.__setModel()
    #self.setTsTrue() if self.tscalc is True else None
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

# --------------------------------- CLASS PIPELINE GRAPHICS --------------------------------- !!!
class Graphics():
  '''
  This class handles all the graphic methods connected to the analysis pipeline. Each public field (self.field)
  can be set accordingly to the user needs from a python script, while private fields (self.__field) is for internal
  usage. Equally, public methods (methodName()) can be invoked within a python script once the class is instanced while
  private methods (__methodName()) should only be used within the class itself.
  '''

  def __init__(self, cfgfile='/config.xml'):
    # init paths ---!
    self.__cfg = xmlConfig(cfgfile)
    self.__p = ConfigureXml(self.__cfg)
    # plotting fields ---!
    self.title = ['title'] # title for
    self.xlabel = ['x']
    self.ylabel = ['y']
    self.grphlabel = ['data']
    self.cbarlabel = ['counts']
    self.colors = ['b']
    self.markers = ['+']
    self.fontsize = 12
    self.cmap = 'jet'
    # axis and images fields ---!
    self.subplots = (1, 1)  # (rows, cols)
    # conditional control ---!
    self.show = False
    self.usetex = True
    self.sharex = 'col'
    self.sharey = 'row'
    # files ---!
    self.imgname = 'image.png'

  # handles SD9 regions WIP ---!
  def __handleReg(self, file):
    with pyregion.open(file) as r:
      r[0].attr[1]['color'] = self.colors
    # NEED TO SAVE CHANGES ---!
    return r

  # plots skymaps ---!
  def showSkymap(self, skyfile, rfile=None, rcolor='k'):
    # add path ---!
    self.imgname = self.__p.getPngDir() + self.imgname
    for id, file in enumerate(skyfile):
      skyfile[id] = self.__p.getDetDir() + file
    for id, file in enumerate(rfile):
      rfile[id] = self.__p.getDetDir() + file
    # set graph ---!
    self.colors = rcolor
    print('before')
    fig, axs = plt.subplots(self.subplots[0], self.subplots[1], sharex=self.sharex, sharey=self.sharey)
    print('after')
    plt.rc('text', usetex=self.usetex)
    # grid of plots ---!
    num = 0
    for cols in range(self.subplots[0]):
      for rows in range(self.subplots[1]):
        num += 1
        print(num)
        # data ---!
        with fits.open(skyfile[num-1]) as hdulist:
          data = hdulist[0].data
          # plot skymap ---!
          axs[cols, rows].imshow(data[num-1], cmap=self.cmap, norm=SymLogNorm(1), interpolation='gaussian',
                     label=self.grphlabel[num-1])
          # plot regions ---!
          if rfile != None:
            r = pyregion.open(rfile[num-1]).as_imagecoord(hdulist[0].header)
            for i in range(len(r)):
              r[i].attr[1]['color'] = self.colors
              patch_list, text_list = r.get_mpl_patches_texts()
              for p in patch_list:
                axs[cols, rows].add_patch(p)
              for t in text_list:
                axs[cols, rows].add_artist(t)

    # decorate plot ---!
    fig.xlabel(self.xlabel)
    fig.ylabel(self.ylabel)
    fig.suptitle(self.title) if self.title != None else None
    fig.colorbar().set_label(self.cbarlabel)
    # save fig ---!
    plt.savefig(self.imgname)
    # show fig ---!
    plt.show() if self.show else None
    plt.close()
    return
