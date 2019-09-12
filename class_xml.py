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
  with open(cfgFile) as json_cfg_file :
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
    self.instr = 'CTA'
    self.maxSrc = 10
    self.corr_rad = 0.1
    self.exclrad = 0.5
    self.default_model = True
    self.srcAtt = [[], [], []]
    self.bkgAtt = [[], [], []]
    self.tscalc = True

  def __getSrcObj(self) :
    src = self.root.findall('source')
    return src

  def __skipNode(self,node,cfg,src) :
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
      if src.attrib[select.attribute] == select.value :
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
    posRaList, posDecList = ([] for i in range(2))
    for src in self.root.findall('source') :
      if self.__skipNode(src, self.__cfg.xml.RaDec) :
        continue
        
      # source candidates ---!
      ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
      dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
      posRaList.append(ra.attrib['value'])
      posDecList.append(dec.attrib['value'])

    self.pos = [posRaList, posDecList]
    return self.pos

  def loadConfInt(self) :
    raList, decList = ([] for i in range(2))
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
    indexList, prefList, pivotList = ([] for i in range(3))
    for src in self.root.findall('source'):
      if self.__skipNode(src, self.__cfg.xml.ConfInt) :
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
    self.__detectionXml = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.xml' % self.sigma)
    self.__detectionReg = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.reg' % self.sigma)

    detection = cscripts.cssrcdetect()
    detection['inmap'] = skymap
    detection['outmodel'] = self.__detectionXml
    detection['outds9file'] = self.__detectionReg
    detection['srcmodel'] = 'POINT'
    detection['bkgmodel'] = self.bkgType.upper()
    detection['threshold'] = int(self.sigma)
    detection['maxsrcs'] = self.maxSrc
    detection['exclrad'] = self.exclrad
    detection['corr_rad'] = self.corr_rad
    detection['corr_kern'] = 'GAUSSIAN'
    detection['logfile'] = self.__detectionXml.replace('.xml', '.log')
    detection['debug'] = bool('no')
    detection.execute()

    return

  def __saveXml(self) :
    self.srcLib.write(self.__detectionXml, encoding="UTF-8", xml_declaration=True,
                 standalone=False, pretty_print=True)
    return

  def __setModel(self) :
    if self.default_model is True :
      Att_Prefactor = {'name':'Prefactor', 'scale':'1e-16', 'value':'5.7', 'min':'1e-07', 'max':'1000.0', 'free':'1'}
      Att_Index = {'name':'Index', 'scale':'-1', 'value':'2.4', 'min':'0', 'max':'5.0', 'free':'1'}
      Att_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'1e-07', 'max':'1000.0', 'free':'0'}
      Bkg_Prefactor = {'name':'Prefactor', 'scale':'1', 'value':'1', 'min':'1e-03', 'max':'1e+3', 'free':'1'}
      Bkg_Index = {'name':'Index', 'scale':'1', 'value':'0.0', 'min':'-5', 'max':'+5.0', 'free':'1'}
      Bkg_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'0.01', 'max':'1000.0', 'free':'0'}

      self.srcAtt = [Att_Prefactor, Att_Index, Att_PivotEn]
      self.bkgAtt = [Bkg_Prefactor, Bkg_Index, Bkg_PivotEn]
      return self.srcAtt, self.bkgAtt
    else :
      pass

  def modXml(self, skymap) :
    self.__runDetection(skymap)
    self.__setModel()

    i = 0
    for src in self.root.findall('source') :
      if self.__skipNode(src, self.__cfg.xml.src) :
        continue

      i += 1
      src.set('tscalc', '1') if self.tscalc is True else None
      # remove spectral component ---!
      rm = src.find('spectrum')
      src.remove(rm)
      # new spectrum ---!
      spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
      spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
      spc.tail = '\n\t\t'.replace('\t', ' ' * 2)
      src.insert(0, spc)
      # new spectral params ---!
      for j in range(3):
        prm = ET.SubElement(spc, 'parameter', attrib=self.srcAtt[j])
        if prm.attrib['name'] == 'Prefactor' and i > 1 :
          prm.set('value', str(float(prm.attrib['value']) / 2 ** (i - 1)))
          prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < 2 else '\n\t\t'.replace('\t', ' ' * 2)
          spc.insert(j, prm)

    for src in self.root.findall('source[@name]') :
      if src.attrib[self.__cfg.idAttribute] in self.__cfg.skip :
        # set bkg attributes ---!
        src.set('instrument', '%s' % self.instr.upper()) if self.instr.capitalize() != 'None' else None
        if self.bkgType.capitalize() == 'Aeff' or self.bkgType.capitalize() == 'Irf' :
          src.set('type', 'CTA%sBackground' % self.bkgType.capitalize())
        if self.bkgType.capitalize() == 'Racc' :
          src.set('type', 'RadialAcceptance')
        # remove spectral component ---!
        rm = src.find('spectrum')
        src.remove(rm)
        # new bkg spectrum ---!
        spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
        spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
        spc.tail = '\n\t'.replace('\t', ' ' * 2)
        # new bkg params ---!
        for j in range(3):
          prm = ET.SubElement(spc, 'parameter', attrib=self.bkgAtt[j])
          prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < 2 else '\n\t\t'.replace('\t', ' ' * 2)

    self.__saveXml()
    return self.__detectionXml, self.__detectionReg

  def FreeFixPrms(self) :
    for src in self.root.findall('source') :
      if self.__skipNode(src, self.__cfg.xml.src) :
        continue

      for free in self.__cfg.xml.src.free :
        src.find('spatialModel/parameter[@name="%s"]' % free.value).set('free', '1')
      for fix in self.__cfg.xml.src.fix :
        src.find('spatialModel/parameter[@name="%s"]' % fix.value).set('free', '0')

    for src in self.root.findall('source') :
      if self.__skipNode(src, self.__cfg.xml.src) :
        for bkg in self.__cfg.xml.src.bkg :
          src.find('spatialModel/parameter[@name="%s"]' % bkg.value).set('free', '1')
        for prm in src not in self.__cfg.xml.src.bkg :
          prm.set('free', '0')

    self.__saveXml()
    return

  def closeXml(self) :
    self.file.close()
    return



