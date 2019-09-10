# ====================== #
# CLASS FOR XML HANDLING #
# ====================== #


# IMPORT ---!
#import xml.etree.ElementTree as ET
#import lxml.etree as ET
from lxml import etree as ET
import numpy as np
import cscripts
import json
import os

def loadConfig() :
  #Load configuration file
  cfgFile = os.path.dirname(__file__)+'/config.json'
  with open(cfgFile) as json_cfg_file:
      cfg = json.load(json_cfg_file)
  return cfg
  
class fsMng:
  def __init__(self, cfg) :
    self.__initPath(cfg.dir)
    
  def __initPath(self,cfg) :
    self.__cfg = cfg.dir
    self.__workdir = cfg.dir.workdir
    self.__runpath = cfg.dir.runpath
    self.selectpath = cfg.dir.selectpath
    self.datapath = cfg.dir.datapath
    self.csvpath = cfg.dir.csvpath
    self.detpath = cfg.dir.detpath
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
  
#  def resolvePath(self.path) :

class xmlMng :
  def __init__(self, xml) :
    fs = fsMng()
    self.file = open(fs.resolvePath(xml))
    self.srcLib = ET.parse(self.file)
    self.root = self.srcLib.getroot()
    self.tsvList = []
    self.pos = [[], []]
    self.err = [[], []]
    self.spectral = [[], [], []]
    self.sigma = 5
    self.bkgType = 'Irf'
    self.maxSrc = 10
    self.corr_rad = 0.1
    self.exclrad = 0.5
    self.srcAtt = [[], [], []]
    self.bkgAtt = [[], [], []]

  def __getSrcName(self) :
    src = self.root.findall('source')
    return src

  def __skipNode(self,node,cfg, src) :
    '''
    :node: node for skip chacks
    :retun true for skip node
    '''
    if src.attrib[cfg.idAttribute] in cfg.skip:
      return True

    for filter in cfg.filters :
      if src.attrib[filter.attribute] == filter.value :
        return True

    if len(cfg.selectors) == 0 :
      return False
      
    for select in cfg.selectors :
      if src.attrib[select.attribute] == select.value:
        return False

    return True # shouldn't not return anything here? Or a logical sum of the previous?

  # RETRIVES TSV v01 ---!
  def loadTSV(self) :
    for src in self.root.findall('source') :
      if self.__skipNode(src, self.__cfg.xml.tsv) :
        continue

      # source candidates ---!
      self.tsvList.append(src.attrib['ts'])
    return self.tsvList

  def loadRaDec(self) :
    posRaList, posDecList = []
    for src in self.root.findall('source'):
      if self.__skipNode(src, self.__cfg.xml.RaDec):
        continue
        
      # source candidates ---!
      ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
      dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
      posRaList.append(ra.attrib['value'])
      posDecList.append(dec.attrib['value'])

    self.pos = [posRaList, posDecList]
    return self.pos

  def loadConfInt(self) :
    raList, decList = []
    for src in self.root.findall('source'):
      if self.__skipNode(src, self.__cfg.xml.ConfInt) :
        continue

      ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
      dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
      raList.append(ra)
      decList.append(dec)

    self.err = [raList, decList]
    return self.err

  def loadSpectral(self) :
    indexList, prefList, pivotList = []
    for src in self.root.findall('source'):
      if self.__skipNode(src, self.__cfg.xml.ConfInt):
        continue

      index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
          src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
      pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
          src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
      pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
          src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
      indexList.append(index)
      prefList.append(pref)
      pivotList.append(pivot)

    self.spectral = [indexList, prefList, pivotList]
    return self.spectral

  def __runDetection(self, skymap) :
    detectionXml = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.xml' % self.sigma)
    detectionReg = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.reg' % self.sigma)

    detection = cscripts.cssrcdetect()
    detection['inmap'] = skymap
    detection['outmodel'] = detectionXml
    detection['outds9file'] = detectionReg
    detection['srcmodel'] = 'POINT'
    detection['bkgmodel'] = self.bkgType.upper()
    detection['threshold'] = int(self.sigma)
    detection['maxsrcs'] = self.maxSrc
    detection['exclrad'] = self.exclrad
    detection['corr_rad'] = self.corr_rad
    detection['corr_kern'] = 'GAUSSIAN'
    detection['logfile'] = detectionXml.replace('.xml', '.log')
    detection['debug'] = bool('no')
    detection.execute()

    return

  def modXml(self, skymap) :
    self.__runDetection(skymap)

    Att_Prefactor = {'name':'Prefactor', 'scale':'1e-16', 'value':'5.7', 'min':'1e-07', 'max':'1000.0', 'free':'1'}
    Att_Index = {'name':'Index', 'scale':'-1', 'value':'2.4', 'min':'0', 'max':'5.0', 'free':'1'}
    Att_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'1e-07', 'max':'1000.0', 'free':'0'}
    self.srcAtt = [Att_Prefactor, Att_Index, Att_PivotEn]
    Bkg_Prefactor = {'name':'Prefactor', 'scale':'1', 'value':'1', 'min':'1e-03', 'max':'1e+3', 'free':'1'}
    Bkg_Index = {'name':'Index', 'scale':'1', 'value':'0.0', 'min':'-5', 'max':'+5.0', 'free':'1'}
    Bkg_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'0.01', 'max':'1000.0', 'free':'0'}
    self.bkgAtt = [Bkg_Prefactor, Bkg_Index, Bkg_PivotEn]


    return