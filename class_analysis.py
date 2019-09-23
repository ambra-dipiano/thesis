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
import json
import xmltodict
import untangle

def xmlConfig() :
  #Load configuration file
  cfgFile = os.path.dirname(__file__)+'/config.xml'
  with open(cfgFile) as fd:
    cfg = untangle.parse(fd.read())
  return cfg.config

class cfgMng() :
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

  def getSelectDir(self) :
    return self.__cfg.dir.selectpath.replace('${runpath}', self.getRunDir())

  def getDataDir(self) :
    return self.__cfg.dir.datapath.replace('${runpath}', self.getRunDir())

  def getCsvDir(self) :
    return self.__cfg.dir.csvpath.replace('${runpath}', self.getRunDir())

  def getDetDir(self) :
    return self.__cfg.dir.detpath.replace('${runpath}', self.getRunDir())

class analysis() :
  def __init__(self):
    global p
    self.__cfg = xmlConfig()
    p = cfgMng(self.__cfg)
    self.model = None
    self.template = None
    self.table = None
    self.tmax = 1800
    self.if_ebl = False
    self.extract_spec = False
    self.__pathout = p.getDataDir()
    self.plot = False
    self.z = 0.1
    self.zfetch = False

  def __openFITS(self):
    global hdul
    hdul = fits.open(self.template)
    return
  def __closeFITS(self):
    hdul.close()

  def __getFitsData(self):
    global energy, time, spectra, ebl
    hdul = self.__openFITS()
    energy = np.array(hdul[1].data)
    time = np.array(hdul[2].data)
    spectra = np.array(hdul[3].data)
    if len(hdul) == 5 :
      ebl = np.array(hdul[4].data)
      self.if_ebl = True
    else :
      self.if_ebl = False
    return

  def __openCSV(self):
    df = pd.read_csv(self.table)
    df.dropna()
    return df

  def __getCsvData(self):
    df = self.__openCSV()
    cols = list(df.columns)
    tau = np.array(df[self.z])
    E = np.array(df[cols[0]]) / 1e3  # MeV --> GeV ---!
    return  tau, E

  def __add_ebl(self):
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
    z = hdul[0].header['REDSHIFT']
    self.z = min([0.5, 0.9, 0.1, 0.09, 1.2], key=lambda x:abs(x-0.097)) # fetch min distance values from given one ---!
    return

  def fits_ebl(self, template_ebl):

    self.__getFitsData()
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
    os.system('rm ' + template_ebl)
    hdul.writeto(template_ebl, overwrite=False)

    if self.plot is True:
      return x, y, x2, y2
    else:
      return

  def __extractSpc(self):
    for i in range(tbin_stop):
      if self.if_ebl is False:
        filename = self.__pathout + 'spec_tbin' + str(i) + '.out'
      else:
        filename = self.__pathout + 'spec_ebl_tbin' + str(i) + '.out'

      if os.path.isfile(filename):
        os.system('rm ' + filename)
      if not os.path.isfile(filename):
        os.system('touch ' + filename)
        out_file = open(filename, 'a')
        # E[MeV], Flux[fotoni/MeV/cm^2/s]
        out_file.close()

    # ebl ---!
    if self.if_ebl is True:
      for i in range(Nt):
        outfile = self.__pathout + 'spec_ebl_tbin' + str(i) + '.out'
        out_file = open(outfile, 'a')
        for j in range(Ne):
          # write spectral data in E [MeV] and I [ph/cm2/s/MeV]
          if ebl is not None and tau is None:
            out_file.write(str(energy[j][0] * 1000) + ' ' + str(ebl[i][j] / 1000) + "\n")
          if tau is not None and ebl is None:
            out_file.write(str(energy[j][0] * 1000.0) + ' ' + str((spectra[i][j] / 1000.0) * np.exp(-tau[j])) + "\n")
        out_file.close()

        os.system('cp ' + self.model + ' ' + self.__pathout + 'run0406_ID000126_ebl_tbin' + str(i) + '.xml')
        s = open(self.__pathout + 'run0406_ID000126_ebl_tbin' + str(i) + '.xml').read()
        s = s.replace('data/spec', 'spec_ebl_tbin' + str(i))
        f = open(self.__pathout + 'run0406_ID000126_ebl_tbin' + str(i) + '.xml', 'w')
        f.write(s)
        f.close()

    # no ebl ---!
    else:
      for i in range(Nt):
        outfile = self.__pathout + 'spec_tbin' + str(i) + '.out'
        out_file = open(outfile, 'a')
        for j in range(Ne):
          # write spectral data in E [MeV] and I [ph/cm2/s/MeV]
          out_file.write(str(energy[j][0] * 1000.0) + ' ' + str(spectra[i][j] / 1000.0) + "\n")
        out_file.close()

        os.system('cp ' + self.model + ' ' + self.__pathout + 'run0406_ID000126_tbin' + str(i))
        s = open(self.__pathout + 'run0406_ID000126_tbin' + str(i) + '.xml').read()
        s = s.replace('spec', 'spec_tbin' + str(i))
        f = open(self.__pathout + 'run0406_ID000126_tbin' + str(i) + '.xml', 'w')
        f.write(s)
        f.close()

    return

  def load_template(self) :
    self.__getFitsData()
    self.__closeFITS()

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

    if self.extract_spec is True and self.if_ebl is True :
      self.__extractSpc()
    if self.extract_spec is True and self.if_ebl is False :
      # here you should call the add_ebl and/or fits_ebl ---!!!!!!!!!!
      self.__extractSpc()

    return
