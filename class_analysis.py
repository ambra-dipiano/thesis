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

def loadConfig() :
  #Load configuration file
  cfgFile = os.path.dirname(__file__)+'/config.json'
  with open(cfgFile) as json_cfg_file:
      cfg = json.load(json_cfg_file)
  return cfg

class fsMng() :
  def __init__(self, cfg) :
    self.__initPath(cfg.dir)

  def __initPath(self, cfg) :
    self.__cfg = cfg.dir
    self.__workdir = cfg.dir.workdir
    self.__runpath = cfg.dir.runpath
    self.__datapath = cfg.dir.datapath
    self.__simpath = cfg.dir.simpath
    self.__selectpath = cfg.dir.selectpath
    self.__detpath = cfg.dir.detpath
    self.__csvpath = cfg.dir.csvpath
    if self.__workdir.endswith(".") == False :
      self.__workdir+'/'
    if self.__runpath.endswith(".") == False :
      self.__runpath+'/'

  def getWorkingDir(self) :
    return self.__workdir
  def setWorkingDir(self, workingDir) :
    self.__workdir = workingDir
    if self.__workdir.endswith(".") == False :
      self.__workdir+'/'

  def getRunDir(self) :
    return self.__runpath.replace('${workdir}', self.getWorkinDir())
  def setRunDir(self, runDir) :
    self.__runpath = runDir
    if self.__runpath.endswith(".") == False :
      self.__runpath+'/'

  def getSelectDir(self) :
    return self.__cfg.dir.selectpath.replace('${runpath}', self.getRunDir())

  def getDataDir(self) :
    return self.__cfg.dir.datapath.replace('${runpath}', self.getRunDir())

  def getCsvDir(self) :
    return self.__cfg.dir.csvpath.replace('${runpath}', self.getRunDir())

  def getDetDir(self) :
    return self.__cfg.dir.detpath.replace('${runpath}', self.getRunDir())