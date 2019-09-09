# ================================ #
# MODULE OF XML HANDLING FUNCTIONS #
# ================================ #


# IMPORT ---!
#import xml.etree.ElementTree as ET
#import lxml.etree as ET
from lxml import etree as ET
import numpy as np
import cscripts
import json
import os

def loadConfig();
  #Load configuration file
  cfgFile = os.path.dirname(a_module.__file__)+'/config.json';
  with open(cfgFile) as json_cfg_file:
      cfg = json.load(json_cfg_file);
  return cfg;
  
class fsMng:
  def __init__(self):
    self.__initPath(cfg.dir);
    
  def __initPath(self,cfg)
    self.__cfg = cfg.dir;
    self.__workdir = cfg.dir.workdir;
    self.__runpath = cfg.dir.runpath;
    self.selectpath = cfg.dir.selectpath;
    self.datapath = cfg.dir.datapath;
    self.csvpath = cfg.dir.csvpath;
    self.detpath = cfg.dir.detpath;
    if(!self.__workdir.endswith(".")) self.__workdir+'/';
    if(!self.__runpath.endswith(".")) self.__runpath+'/';
  
  def getWorkingDir(self)
    return self.__workdir;
  def setWorkingDir(self, workingDir)
    self.__workdir = workingDir;
    if(!self.__workdir.endswith(".")):
      self.__workdir+'/';

  def getRunDir(self)
    return self.__runpath.replace('${workdir}',self.getWorkinDir());
  def setRunDir(self, runDir)
    self.__runpath = runDir;
    if(!self.__runpath.endswith(".")):
      self.__runpath+'/';

  def getSelectDir(self)
    return self.__cfg.dir.selectpath.replace('${runpath}',self.getRunDir());

  def getDataDir(self)
    return self.__cfg.dir.datapath.replace('${runpath}',self.getRunDir());

  def getCsvDir(self)
    return self.__cfg.dir.csvpath..replace('${runpath}',self.getRunDir());

  def getDetDir(self)
    return self.__cfg.dir.detpath.replace('${runpath}',self.getRunDir());
  
  def resolvePath(self.path)
  
'''
:param likeXml: likelihood result model definition xml file (str)
'''
class xmlMng:
  def __init__(self, likeXml):
    fs = fsMng();
    self.file = open(fs.resolvePath(likeXml));
    self.srcLib = ET.parse(file);
    self.root = srcLib.getroot();
    self.tsvList = [];
    self.pos = [[], []];
    self.err = [[], []];

  def __skipNode(self,node,cfg)
    '''
    :node: node for skip chacks
    :retun true for skip node
    '''
    if src.attrib[cfg.idAttribute] in cfg.skip:
      return true;

    for filter in cfg.filters
      if src.attrib[filter.attribute] == filter.value:
        return true;

    if count(cfg.selectors)==0:
      return false;
      
    for select in cfg.selectors
      if src.attrib[select.attribute] == select.value:
        return false;

    return true;

  # RETRIVES TSV v01 ---!
  def loadTSV(self) :
    for src in self.root.findall('source'):
      if self.__skipNode(src,cfg.xml.tsv):
        continue;

      # source candidates ---!
      self.tsvList.append(src.attrib['ts'])

  def loadRaDec(self) :
    posRaList = [];
    posDecList = [];
    errRaList = [];
    errDecList = [];
    for src in root.findall('source'):
      if self.__skipNode(src,cfg.xml.raDec):
        continue;
        
      # source candidates ---!
      posRaList.append(ra.attrib['value'])
      posDecList.append(dec.attrib['value'])
      errRaList.append(ra.attrib['error'])
      errDecList.append(dec.attrib['error'])

    self.pos = [posRaList,posDecList]
    self.err = [errRaList,errDecList]
    
    
