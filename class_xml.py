# ====================== #
# CLASS FOR XML HANDLING #
# ====================== #


# IMPORT ---!
from lxml import etree as ET
from operator import attrgetter
import json
import os
from class_analysis import (xmlConfig, cfgMng_xml)


def jsonConfig():
  # Load configuration file
  cfgFile = os.path.dirname(__file__) + '/config.json'
  with open(cfgFile) as json_cfg_file:
    cfg = json.load(json_cfg_file)
  return cfg

# --------------------------------- CLASS json CONFIGURATION --------------------------------- !!!

class cfgMng_json():
  cfg = jsonConfig()

  def __init__(self, cfg):
    self.__initPath(cfg.get('dir'))

  def __initPath(self, cfg):
    self.__cfg = cfg
    self.__workdir = self.__cfg.get('workdir')
    self.__runpath = self.__cfg.get('runpath')
    self.datapath = self.__cfg.get('datapath')
    self.simpath = self.__cfg.get('simpath')
    self.selectpath = self.__cfg.get('selectpath')
    self.detpath = self.__cfg.get('detpath')
    self.csvpath = self.__cfg.get('csvpath')
    if self.__workdir.endswith(".") == False:
      self.__workdir = self.__workdir + '/'
    if self.__runpath.endswith(".") == False:
      self.__runpath = self.__runpath + '/'

  def getWorkingDir(self):
    return self.__workdir

  def setWorkingDir(self, workingDir):
    self.__workdir = workingDir
    if self.__workdir.endswith(".") == False:
      self.__workdir = self.__workdir + '/'

  def getRunDir(self):
    return self.__runpath.replace('${workdir}', self.getWorkingDir())

  def setRunDir(self, runDir):
    self.__runpath = runDir
    if self.__runpath.endswith(".") == False:
      self.__runpath = self.__runpath + '/'

  def getDataDir(self):
    return self.datapath.replace('${runpath}', self.getRunDir())

  def getSimDir(self):
    return self.simpath.replace('${runpath}', self.getRunDir())

  def getSelectDir(self):
    return self.selectpath.replace('${runpath}', self.getRunDir())

  def getDetDir(self):
    return self.detpath.replace('${runpath}', self.getRunDir())

  def getCsvDir(self):
    return self.csvpath.replace('${runpath}', self.getRunDir())

  def resolvePath(self):
    return

# --------------------------------- CLASS xml HANDLING --------------------------------- !!!

class xmlMng():

  def __init__(self, xml):
    self.__xml = xml
    self.__cfg = jsonConfig()
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

  def __skipNode(self, src, cfg):
    '''
    :retun true for skip node
    '''
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

  def loadTSV(self):
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
      # if self.__skipNode(src, self.__cfg.get('xml').get('tsv')):
      #   continue
        tsv = src.attrib['ts']
        self.tsvList.append(tsv)
    return self.tsvList

  def loadRaDec(self):
    posRaList, posDecList = ([] for i in range(2))
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
      # if self.__skipNode(src, self.__cfg.get('xml').get('RaDec')):
      #   continue
        ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
        dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
        posRaList.append(ra)
        posDecList.append(dec)
    self.pos = [posRaList, posDecList]
    return self.pos

  def loadConfInt(self):
    raList, decList = ([] for i in range(2))
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
      # if self.__skipNode(src, self.__cfg.get('xml').get('ConfInt')):
      #   continue
        ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
        dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
        raList.append(ra)
        decList.append(dec)
    self.err = [raList, decList]
    return self.err

  def loadSpectral(self):
    global cutoffList
    indexList, prefList, pivotList = ([] for i in range(3))
    if self.if_cut is True :
      cutoffList = []

    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
      # if self.__skipNode(src, self.__cfg.get('xml').get('src')):
      #   continue
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

    if self.if_cut is False:
      self.spectral = [indexList, prefList, pivotList]
    else:
      self.spectral = [indexList, prefList, pivotList, cutoffList]
    return self.spectral

  def __saveXml(self):
    self.srcLib.write(self.__xml, encoding="UTF-8", xml_declaration=True,
                      standalone=False, pretty_print=True)
    return self.__xml

  def __setModel(self):
    if self.default_model is True:
      Att_Prefactor = {'name': 'Prefactor', 'scale': '1e-16', 'value': '5.7', 'min': '1e-07', 'max': '1000.0', 'free': '1'}
      Att_Index = {'name': 'Index', 'scale': '-1', 'value': '2.4', 'min': '0', 'max': '5.0', 'free': '1'}
      Att_PivotEn = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1', 'min': '1e-07', 'max': '1000.0', 'free': '0'}
      Bkg_Prefactor = {'name': 'Prefactor', 'scale': '1', 'value': '1', 'min': '1e-03', 'max': '1e+3', 'free': '1'}
      Bkg_Index = {'name': 'Index', 'scale': '1', 'value': '0.0', 'min': '-5', 'max': '+5.0', 'free': '1'}
      Bkg_PivotEn = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1', 'min': '0.01', 'max': '1000.0', 'free': '0'}

      self.srcAtt = [Att_Prefactor, Att_Index, Att_PivotEn]
      self.bkgAtt = [Bkg_Prefactor, Bkg_Index, Bkg_PivotEn]
      if self.if_cut is True:
        Att_CutOff = {'name': 'CutoffEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '1'}
        self.srcAtt.append(Att_CutOff)

      return self.srcAtt, self.bkgAtt
    else:
      pass

  def modXml(self):
    self.__setModel()
    # source ---!
    i = 0
    for src in self.root.findall('source'):
      i += 1
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
      # if self.__skipNode(src=src, cfg=self.__cfg.get('xml').get('src')):
      #   continue
        src.set('tscalc', '1') if self.tscalc is True else None
        # remove spectral component ---!
        rm = src.find('spectrum')
        src.remove(rm)
        # new spectrum ---!
        if self.if_cut is True:
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
      # for src in self.root.findall('source[@name]'):
      #   if self.__skipNode(src=src, cfg=self.__cfg.get('xml').get('src')):
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

    self.__xml = self.__saveXml()
    return

  def FreeFixPrms(self):
    for src in self.root.findall('source'):
      if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
      # if self.__skipNode(src, self.__cfg.xml.src):
      #   continue
        for free in self.__cfg.xml.src.free:
          src.find('spatialModel/parameter[@name="%s"]' % free.value).set('free', '1')
        for fix in self.__cfg.xml.src.fix:
          src.find('spatialModel/parameter[@name="%s"]' % fix.value).set('free', '0')
      else:
      # if self.__skipNode(src, self.__cfg.xml.src):
        for bkg in self.__cfg.xml.src.bkg:
          src.find('spatialModel/parameter[@name="%s"]' % bkg.value).set('free', '1')
        for prm in src not in self.__cfg.xml.src.bkg:
          prm.set('free', '0')

    self.__saveXml()
    return

  def sortSrcTS(self):
    src = self.root.findall("*[@ts]")
    self.root[:-1] = sorted(src, key=lambda el: (el.tag, el.attrib['ts']), reverse=True)
    self.__saveXml()
    return

  def closeXml(self):
    self.file.close()
    return

