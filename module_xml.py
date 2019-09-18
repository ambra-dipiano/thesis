# ================================ #
# MODULE OF XML HANDLING FUNCTIONS #
# ================================ #

# IMPORT ---!
#import xml.etree.ElementTree as ET
#import lxml.etree as ET
from lxml import etree as ET
import numpy as np
import cscripts


# RETRIVES TSV v01 ---!
def getTSV(likeXml) :
  '''

  :param likeXml: likelihood result model definition xml file (str)
  :return: tsv: list of TSV for each src (list of float)
  '''

  file = open(likeXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  i = 0
  tsvList = []
  for src in root.findall('source'):
    i += 1
    # source candidates ---!
    if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel' :
      tsv = src.attrib['ts'] if src.attrib['tscalc'] == '1' else None
      tsvList.append(tsv)

  file.close()
  return tsvList



# RETRIVES RA & DEC v01 ---!
def getRaDec(likeXml) :
  '''
  to use if max likelihood fit has been ran with RA & DEC as free.
  :param likeXml: likelihood result model definition xml file (str)
  :return: pos: source positions after likelihood (list of float)
  '''

  file = open(likeXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  i = 0
  raList = []
  decList = []
  for src in root.findall('source'):
    i += 1

    # source candidates ---!
    if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel' :
      ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
      dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
      raList.append(ra)
      decList.append(dec)

  pos = [raList, decList]
  file.close()

  return pos



# RETRIVES RA & DEC ERRORS v01 ---!
def getRaDec_errors(likeXml) :
  '''
  to use if max likelihood fit has been ran with RA & DEC as free.
  :param likeXml: likelihood result model definition xml file (str)
  :return: pos: source positions after likelihood (list of float)
  '''

  file = open(likeXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  i = 0
  raList = []
  decList = []
  for src in root.findall('source'):
    i += 1

    # source candidates ---!
    if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel' :
      ra = src.find('spatialModel/parameter[@name="RA"]').attrib['error']
      dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['error']
      raList.append(ra)
      decList.append(dec)

  err = [raList, decList]
  file.close()

  return err



# RETRIVES CONFIDENCE LEVEL ERRORS v01 ---!
def getConfInt_gauss(errors) :
  '''
  '''

  file = open(errors)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  i = 0
  raList = []
  decList = []
#  prefList = []
#  indexList = []
  for src in root.findall('source'):
    i += 1

    # source candidates ---!
    if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel' :
      ra = src.find('spatialModel/parameter[@name="RA"]').attrib['value']
      dec = src.find('spatialModel/parameter[@name="DEC"]').attrib['value']
#      prefactor = src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']
#      index = src.find('spectrum/parameter[@name="Index"]').attrib['value']
      raList.append(ra)
      decList.append(dec)


  err = [raList, decList]
  file.close()

  return err


# RETRIVES SPECTRAL VALUES v01 ---!
def getSpectral(likeXml) :
  '''
  to use if max likelihood fit has been ran with RA & DEC as free.
  :param likeXml: likelihood result model definition xml file (str)
  :return: pos: source positions after likelihood (list of float)
  '''

  file = open(likeXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  i = 0
  indexList = []
  prefactorList = []
  pivotList = []
  indexList_error = []
  prefactorList_error = []
  for src in root.findall('source'):
    i += 1

    # source candidates ---!
    if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel' :
      prefactor = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value'])*float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
      index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value'])*float(src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
      pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value'])*float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
      prefactor_err = src.find('spectrum/parameter[@name="Prefactor"]').attrib['error'] 
      index_err = src.find('spectrum/parameter[@name="Index"]').attrib['error']
      prefactorList.append(prefactor)
      indexList.append(index)
      pivotList.append(pivot)
      prefactorList_error.append(prefactor_err)
      indexList_error.append(index_err)

  val = [prefactorList, indexList, pivotList]
  err = [prefactorList_error, indexList_error]
  file.close()

  return val, err


# RUN DETECTION AND MODEL SPECTRAL COMPONENTS v09 ---!
def srcDetection_spcModeling(skymap, sigma=5, instr='CTA', bkgType='Irf', src_attrib='none', bkg_attrib='none', tsv=True, maxSrc=20, exclrad = 0.5, if_ebl=True) :
  ''''
  Runs cssrcdetect tool to detected src candidates in a counts map. The listing file is modeled in its spectral components.
  :param:
  skymap = name of counts map file (str)
  sigma = Gaussian threshold for detection acceptance, takes 5 as default (int)
  instr = <CTA|HESS|MAGIC|VERITAS> takes CTA as default (str)
  bkgType = <None|Irf|Aeff> takes Irf as default (str)
  src_attrib = list of dicts containing source spectral model params, takes 'none' as default
                and builds a default model (list of dicts or dict)
  bkg_attrib = list of dicts containing bkg spectral model params, takes 'none' as default
                and builds a default model (list of dicts or dict)
  tsv = boolean param for Test Statistic Value setting, takes True as default (bool)
  maxsrc = max number of src candidates to be detected, takes 20 as default (int)
  exclrad = radius from detection to exclude in further search, takes 1deg as default (float)
  :return:
  srcLib = detection modified (ET obj)
  detectionXml_mod = detection xml modified file (str)
  detectionReg = DS9 region for detected candidates (str)
  posList = Ra & DEC of detected src candidates (2d list of floats)
  '''

  # filenames for detection model definition XML and region DS9 files ---!
  detectionXml = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.xml' % sigma)
  detectionReg = '%s' % skymap.replace('_skymap.fits', '_det%ssgm.reg' % sigma)

  # RUN DETECTION ---!
  detection = cscripts.cssrcdetect()
  detection['inmap'] = skymap
  detection['outmodel'] = detectionXml
  detection['outds9file'] = detectionReg
  detection['srcmodel'] = 'POINT'
  detection['bkgmodel'] = bkgType.upper()
  detection['threshold'] = int(sigma)
  detection['maxsrcs'] = maxSrc
  detection['exclrad'] = exclrad
  detection['corr_rad'] = 0.1
  detection['corr_kern'] = 'GAUSSIAN'
  detection['logfile'] = detectionXml.replace('.xml', '.log')
  detection['debug'] = bool('no')
  detection.execute()

  # MODELING THE SPECTRAL COMPONENT OF DETECTED CANDIDATE SOURCES ---!

  file = open(detectionXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  # src ---!
  # default dicts of params attributes ---!
  if src_attrib == 'none' :
    Att_Prefactor = {'name':'Prefactor', 'scale':'1e-16', 'value':'5.7', 'min':'1e-07', 'max':'1000.0', 'free':'1'}
    Att_Index = {'name':'Index', 'scale':'-1', 'value':'2.4', 'min':'0', 'max':'5.0', 'free':'1'}
    Att_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'1e-07', 'max':'1000.0', 'free':'0'}
    srcAtt = [Att_Prefactor, Att_Index, Att_PivotEn]
    if if_ebl is True :
      Att_CutOff = {'name': 'CutoffEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '1'}
      srcAtt.append(Att_CutOff)

  # alternative input dicts of params attributes ---!
  else :
    srcAtt = src_attrib

  # bkg ---!
  # default dicts of params attributes ---!
  if bkg_attrib == 'none' :
    Bkg_Prefactor = {'name':'Prefactor', 'scale':'1', 'value':'1', 'min':'1e-03', 'max':'1e+3', 'free':'1'}
    Bkg_Index = {'name':'Index', 'scale':'1', 'value':'0.0', 'min':'-5', 'max':'+5.0', 'free':'1'}
    Bkg_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'0.01', 'max':'1000.0', 'free':'0'}
    bkgAtt = [Bkg_Prefactor, Bkg_Index, Bkg_PivotEn]

  # alternative input dicts of params attributes ---!
  else :
    bkgAtt = bkg_attrib

  i = 0
  raList = []
  decList = []
  for src in root.findall('source'):
    i += 1

    # source candidates ---!
    if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel' :
      # set TSV ---!
      src.set('tscalc', '1') if tsv == True else None
      # remove spectral component ---!
      rm = src.find('spectrum')
      src.remove(rm)
      # new spectrum ---!
      spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
      spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
      spc.tail = '\n\t\t'.replace('\t', ' ' * 2)
      src.insert(0, spc)

      # new spectral params ---!
      for j in range(len(srcAtt)):
        prm = ET.SubElement(spc, 'parameter', attrib=srcAtt[j])
        prm.set('value', str(float(prm.attrib['value']) / 2 ** (i - 1))) if prm.attrib['name'] == 'Prefactor' \
                                                                            and i > 1 else None
        prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < (len(srcAtt)-1) else '\n\t\t'.replace('\t', ' ' * 2)
        spc.insert(j, prm)

      # store detected src positions (RA & DEC) ---!
      raList.append(src.find('spatialModel/parameter[@name="RA"]').attrib['value'])
      decList.append(src.find('spatialModel/parameter[@name="DEC"]').attrib['value'])



    # background ---!
    else:
      # set bkg attributes ---!
      src.set('instrument', '%s' % instr.upper()) if instr.capitalize() != 'None' else None
      src.set('type', 'CTA%sBackground' % bkgType.capitalize()) if bkgType.capitalize() == 'Aeff' \
                                                                   or bkgType.capitalize() == 'Irf' else None
      src.set('type', 'RadialAcceptance') if bkgType.capitalize() == 'Racc' else None
      # remove spectral component ---!
      rm = src.find('spectrum')
      src.remove(rm)
      # new bkg spectrum ---!
      spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
      spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
      spc.tail = '\n\t'.replace('\t', ' ' * 2)
      # new bkg params ---!
      for j in range(len(bkgAtt)):
        prm = ET.SubElement(spc, 'parameter', attrib=bkgAtt[j])
        prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < 2 else '\n\t\t'.replace('\t', ' ' * 2)

  # store ra & dec in position list ---!
  posList = [raList, decList]

  # filename for modified detection model definition XML file ---!
  #detectionXml_mod = detectionXml.replace('.xml', 'Mod.xml')
  # write new file ---!
  srcLib.write(detectionXml, encoding="UTF-8", xml_declaration=True,
               standalone=False, pretty_print=True)

  file.close()

  return detectionXml, detectionReg, posList



# MODELING XML DETECTION FILE HES TYPE v01 ---!
def getAttribs_spcModel(instr='CTA', bkgType='Irf', srcType='Powerlaw'):
  ''''
  :param:
  instr = <NONE|CTA|HESS|MAGIC|VERITAS> takes CTA as default (str)
  bkgType = <None|Irf|Aeff|Cube> takes Irf as default (str)
  srcType = <PowerLaw>
  :return:
  bkgAtt = collection of params attributes for bkg modeling (list of dicts)
  srcAtt = collection of params attributes for src modeling (list of dicts)
  '''

  # SOURCE ---!
  if instr.upper() == 'HESS' and srcType.lower() == 'powerlaw' :

    # def use case for HESS, PointSource and Powerlaw
    Att_Prefactor = {'name': 'Prefactor', 'scale': '1e-16', 'value': '5.7', 'min': '1e-07', 'max': '1000.0', 'free': '1'}
    Att_Index = {'name': 'Index', 'scale': '-1', 'value': '2.4', 'min': '0', 'max': '+5.0', 'free': '1'}
    Att_PivotEn = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '0.3', 'min': '1e-07', 'max': '1000.0', 'free': '0'}

  else :

    # def default for CTA, PointSource and PowerLaw ---!
    Att_Prefactor = {'name':'Prefactor', 'scale':'1e-16', 'value':'5.7', 'min':'1e-07', 'max':'1000.0', 'free':'1'}
    Att_Index = {'name':'Index', 'scale':'-1', 'value':'2.4', 'min':'0', 'max':'5.0', 'free':'1'}
    Att_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'0.3', 'min':'1e-07', 'max':'1000.0', 'free':'0'}

  # BACKGROUND ---!
  if instr.upper() == 'HESS' and bkgType.capitalize() == 'Aeff' :

    # def use case for HESS, CTAAeffBackground and Powerlaw
    Bkg_Prefactor = {'name': 'Prefactor', 'scale': '1e-13', 'value': '1', 'min': '0', 'free': '1'}
    Bkg_Index = {'name': 'Index', 'scale': '-2.5', 'value': '1', 'min': '-4', 'max': '4', 'free': '1'}
    Bkg_PivotEn = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1', 'free': '0'}

  else :

    # def default for CTA, CTAIrfBackground and PowerLaw ---!
    Bkg_Prefactor = {'name':'Prefactor', 'scale':'1', 'value':'1', 'min':'1e-03', 'max':'1e+3', 'free':'1'}
    Bkg_Index = {'name':'Index', 'scale':'1', 'value':'0.0', 'min':'-5', 'max':'+5.0', 'free':'1'}
    Bkg_PivotEn = {'name':'PivotEnergy', 'scale':'1e6', 'value':'1', 'min':'0.01', 'max':'1000.0', 'free':'0'}


  srcAtt = [Att_Prefactor, Att_Index, Att_PivotEn]
  bkgAtt = [Bkg_Prefactor, Bkg_Index, Bkg_PivotEn]

  return bkgAtt, srcAtt



# ON/OFF MODEL v01 ---!
def writeXml_onoffModel(modelXml):
  '''
  Given a model definition XML file it modifies its background source instrument attributes in 'CTAOnOff'
  and saves on a new file.
  :param
  modelXml = model xml file to modify
  :return:
  srcLib = on/off model (ET obj)
  modelXml_mod = on/oof model file (str)
  '''

  # parse and getroot ---!
  file = open(modelXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  # search for bkg and set new instrument ---!
  bkg = root.find('source[@type="CTAIrfBackground"]')
  bkg.set('instrument', 'CTAOnOff')

  # write on new file ---!
  srcLib.write(modelXml.replace('model', 'onoff_model'), encoding="UTF-8", xml_declaration=True,
              standalone=False, pretty_print=True)
  modelXml_mod = modelXml.replace('model', 'onoff_model')
  #srcLib.write(modelXml_mod, encoding="UTF-8", xml_declaration=True)

  file.close()

  return srcLib, modelXml_mod



# FREE PARAMS v01 ---!
def freeParams(modelXml, spt_prms='none', spc_prms='none', free_bkg=False, bkg_prms='none') :
  '''

  :param modelXml: which file to parse (str)
  :param spc_prms, spt_prms: which parameters to free (list)
  :return:
  '''
  file = open(modelXml)
  srcLib = ET.parse(file)
  root = srcLib.getroot()

  for src in root.findall('source'):

    # free spatial prms ---!
    if (src.attrib['name'] != 'Background' or src.attrib['name'] != 'CTABackgroundModel') and spt_prms != 'none' :
      for i in range(len(spt_prms)) :
        src.find('spatialModel/parameter[@name="%s"]' % spt_prms[i]).set('free', '1')

    # free spectral prms ---!
    if (src.attrib['name'] != 'Background' or src.attrib['name'] != 'CTABackgroundModel') and spc_prms != 'none' :
      for i in range(len(spc_prms)) :
        src.find('spectrum/parameter[@name="%s"]' % spc_prms[i]).set('free', '1')

    # free bkg prms ---!
    if free_bkg == True :
      if (src.attrib['name'] == 'Background' or src.attrib['name'] == 'CTABackgroundModel') and bkg_prms != 'none' :
        for i in range(len(bkg_prms)) :
          src.find('spectrum/parameter[@name="%s"]' % bkg_prms[i]).set('free', '1')

  # save ---!
  srcLib.write(modelXml, encoding="UTF-8", xml_declaration=True,
               standalone=False, pretty_print=True)
  file.close()

  return



# FIX PARAMS v01 ---!
def fixParams(modelXml, spt_prms='none', spc_prms='none', fix_bkg=False, bkg_prms='none') :
  '''

  :param modelXml: which file to parse (str)
  :param spc_prms, spt_prms: which parameters to free (list)
  :return:
  '''

  file = open(modelXml)
  srcLib = ET.parse(modelXml)
  root = srcLib.getroot()

  for src in root.findall('source'):

    # free spatial prms ---!
    if (src.attrib['name'] != 'CTABackgroundModel' or src.attrib['name'] != 'Background') and spt_prms != 'none' :
      for i in range(len(spt_prms)) :
        src.find('spatialModel/parameter[@name="%s"]' % spt_prms[i]).set('free', '0')

    # free spectral prms ---!
    if (src.attrib['name'] != 'CTABackgroundModel' or src.attrib['name'] != 'Background') and spc_prms != 'none' :
      for i in range(len(spc_prms)) :
        src.find('spectrum/parameter[@name="%s"]' % spc_prms[i]).set('free', '0')

    # free bkg prms ---!
    if fix_bkg == True :
      if (src.attrib['name'] == 'CTABackgroundModel' or src.attrib['name'] == 'Background') and bkg_prms != 'none' :
        for i in range(len(bkg_prms)) :
          src.find('spectrum/parameter[@name="%s"]' % bkg_prms[i]).set('free', '0')  
  
  # save ---!
  srcLib.write(modelXml, encoding="UTF-8", xml_declaration=True,
               standalone=False, pretty_print=True)
  file.close()

  return



# BKG MODELING ---!
def bkgModeling(modelXml, bkgXml, spectral=True, spatial=False) :
  '''

  :param modelXml: model definition XML file to modify (str)
  :param bkgXml:  model definition XML file of backgroun (str)
  :param spectral, spatial: parameters to modify if True (bool)
  :return: bkgObj, bkgMod
  '''
  file_src = open(modelXml)
  file_bkg = open(bkgXml)
  srcLib = ET.parse(file_src)
  bkgLib = ET.parse(file_bkg)
  rootSrc = srcLib.getroot()
  rootBkg = bkgLib.getroot()

  # spectral ---!
  if spectral == True :
    spc_new = rootBkg.find('source/spectrum')
    spc_old = rootSrc.find('source[@name="Background"]')
    rm = spc_old.find('spectrum')
    spc_old.remove(rm)
    spc_old.append(spc_new)

  # spatial ---!
  if spatial == True :
    spt_new = rootBkg.find('source/spatialModel[@type="Multiplicative"]')
    spt_old = rootSrc.find('source[@name="Background"]')
    rm = spt_old.find('spatialModel')
    spt_old.remove(rm)
    spt_old.append(spt_new)

  # filename for modified detection model definition XML file ---!
  filename = modelXml.replace('Mod.xml', 'Bkg.xml')
  # write new file ---!
  srcLib.write(filename, encoding="UTF-8", xml_declaration=True,
               standalone=False, pretty_print=True)

  file_src.close()
  file_bkg.close()

  return filename




# CREATE MODEL v01 ---!
def writeXml_obsModel_bkgIrf(simNum, Nsrc=1):
  ''''
  Creates a new model definition XML file, the first source is set to have Crab parameters values.
  :param:
  simNum = prefix of filename to which append '_model.xml'
  Nsrc = number of sources in 5 deg squared centered around Crab
  :return:
  modelXml = model (ET obj)
  modelXml_name = model file (str)
  '''

  modelXml_name = '%s_model.xml' % simNum
  file = open(modelXml_name, 'w+').close()

  # define src attributes dictionaries and values lists ---!
  srcName = ['Src0%02d' % (i+1) for i in range(int(Nsrc))]
  prmName = ('Prefactor', 'Index', 'PivotEnergy', 'RA', 'DEC')
  prmScale = (1e-16, -1, 1e-16, 1, 1)
  PrefactorValue = [11.4/2**(i+1) for i in range(int(Nsrc))]
  RaValue = [np.random.uniform(81.13, 86.13) for i in range(int(Nsrc)-1)]
  RaValue.insert(0, 83.6331)
  DecValue = [np.random.uniform(24.51, 19.51) for i in range(int(Nsrc)-1)]
  DecValue.insert(0, 22.0145)
  prmValue = (1, 2.48, 0.3, 1, 1)
  prmMin = (1e-07, 0.0, 0.01, -360, -90)
  prmMax = (100.0, 5.0, 100.0, 360, 90)
  prmFree = (1, 1, 0, 0, 0)

  # build xml src elements ---!
  root = ET.Element('source_library', attrib={'title': 'source library'})
  root.text = '\n\t'.replace('\t', ' ' * 2)
  root.tail = '\n\t'.replace('\t', ' ' * 2)

  # insert nodes ---!
  for i in range(int(Nsrc)) :
    src = ET.SubElement(root, 'source', attrib={'title':tuple(srcName)[i], 'type': 'PointSource'} )
    src.text = '\n\t\t'.replace('\t', ' ' * 2)
    src.tail = '\n\t'.replace('\t', ' ' * 2)

    # spectral component ---!
    spectrum = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
    spectrum.text = '\n\t\t\t'.replace('\t', ' ' * 2)
    spectrum.tail = '\n\t\t'.replace('\t', ' ' * 2)

    for j in range(3) :
      param = ET.SubElement(spectrum, 'parameter', attrib={'name':str(prmName[j]), 'scale': str(prmScale[j]), 'value':str(prmValue[j]),
                                                           'min': str(prmMin[j]), 'max' : str(prmMax[j]), 'free':str(prmFree[j])})
      param.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < 2 else '\n\t\t'.replace('\t', ' ' * 2)

    # spatial component ---!
    spatial = ET.SubElement(src, 'spatialModel', attrib={'type': 'PointSource'})
    spatial.text = '\n\t\t\t'.replace('\t', ' ' * 2)
    spatial.tail = '\n\t'.replace('\t', ' ' * 2)

    for j in range(2) :
      param = ET.SubElement(spatial, 'parameter', attrib={'name':str(prmName[j+3]), 'scale': str(prmScale[j+3]), 'value':str(prmValue[j+3]),
                                                           'min': str(prmMin[j+3]), 'max' : str(prmMax[j+3]), 'free':str(prmFree[j+3])})
      param.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < 1 else '\n\t\t'.replace('\t', ' ' * 2)

    # change Prefactor, Ra, Dec values for each source ---!
    for prm in spectrum.findall('parameter[@name="Prefactor"]') :
      prm.set('value', str(PrefactorValue[i]))

    for prm in spatial.findall('parameter') :
      prm.set('value', str(RaValue[i]) if prm.attrib['name'] == 'RA' else str(DecValue[i]))

  # redefine attributes dictionaries for bkg ---!
  prmScale = (1.0, 1.0, 1e6)
  prmValue = (1.0, 0.0, 1.0)
  prmMin = (1e-3, -5.0, 0.01)
  prmMax = (1e3, 5.0, 1000.0)
  prmFree = (1, 1, 0)

  # insert bkg node ---!
  bkg = ET.SubElement(root, 'source', attrib={'name':'CTABackgroundModel', 'type':'CTAIrfBackground', 'instrument':'CTA'})
  bkg.text = '\n\t\t'.replace('\t', ' ' * 2)
  bkg.tail = '\n'

  # bkg spectral component ---!
  spectrum = ET.SubElement(bkg, 'spectrum', attrib={'type': 'PowerLaw'})
  spectrum.text = '\n\t\t\t'.replace('\t', ' ' * 2)
  spectrum.tail = '\n\t'.replace('\t', ' ' * 2)

  # bkg params ---!
  for j in range(3):
    param = ET.SubElement(spectrum, 'parameter', attrib={'name': str(prmName[j]), 'scale': str(prmScale[j]), 'value': str(prmValue[j]),
                                  'min' : str(prmMin[j]), 'max': str(prmMax[j]), 'free': str(prmFree[j])})
    param.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < 2 else '\n\t\t'.replace('\t', ' ' * 2)

  # create file obj ---!
  modelXml = ET.ElementTree(root)

  # save on new file ---!
  modelXml.write(modelXml_name, encoding="UTF-8", xml_declaration=True, standalone=False, pretty_print=True)
  #modelXml.close()

  return modelXml_name
