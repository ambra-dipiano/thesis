# =============== #
# CRAB LIGHTCURVE #
# =============== #

# IMPORTS ---!
import gammalib
import ctools
import cscripts
#from module_plot import *
from module_xml import *

caldb = 'prod3b'
irf = 'South_z20_average_30m'

# work with absolute paths ---!
workdir = '/mnt/nvme0n1p1/piano_analysis/working-dir/crab_lc/'
runpath = workdir 
simpath = runpath + 'sim/'
selectpath = runpath + 'selected_sim/'
datapath = runpath + 'data/'
csvpath = runpath + 'csv/'
detpath = runpath + 'detection_all/'
fileroot = 'crab'

# ==================
# !!! CRAB MODEL !!!
# ==================

#model = writeXml_obsModel_bkgIrf('%scrab' % datapath, Nsrc=1)
model = datapath + 'crab.xml'

# =====================
# !!! SET UP TRIALS !!!
# =====================

# trials
#trials = 10
# trials count
count = 0
sigma=5

texp = 2000
tmin = 0
tmax = tmin+texp
#for i in range(tint):
#  tmax.append(tmin + texp[i])

emin = 0.1  # 0.1 TeV
emax = 100.0  # 100 TeV
roi = 5


# ====================
# !!! SIMULATE CRAB !!!
# ====================

# pointing with off-axis equal to max prob GW ---!
offmax = [0., 0.]
#RA = 83.633212   
#DEC = 22.014460
RA = 83.633
DEC = 22.014
pointRA = RA + offmax[0]
pointDEC = DEC + offmax[1]

# time setting ---!
start = 0
stop = 2000

#print('!!! check start simulations ------ ')

# simulate N independent photon-lists of the same event ---!
count += 1
event = simpath + fileroot + ".fits"

print('!!! check start simulation ----')
sim = ctools.ctobssim()
sim["inmodel"] = model 
sim["outevents"] = event 
sim["caldb"] = caldb
sim["irf"] = irf
sim["ra"] = pointRA
sim["dec"] = pointDEC
sim["rad"] = 5.0
sim["tmin"] = start
sim["tmax"] = stop
sim["emin"] = emin
sim["emax"] = emax
sim["seed"] = count
sim["logfile"] = simpath + fileroot + ".log" 
sim.execute() 

print('!!! check end simulation ---')

# combine in observation list
xml = gammalib.GXml()
obslist = xml.append('observation_list title="observation library"')
obs = obslist.append('observation name="Crab" id="%06d" instrument="CTA"' % (count))
obs.append('parameter name="EventList" file="%s%s.fits"' % (simpath, fileroot))

eventList = '%sobs_%s.xml' % (detpath, fileroot)
xml.save(eventList)
print('!!! check --- eventList: ', eventList)

# ==============
# !!! SKYMAP !!!
# ==============

skymapName = event.replace(simpath, detpath).replace('.fits', '_skymap.fits')

skymap = ctools.ctskymap()
skymap['inobs'] = event
skymap['outmap'] = skymapName
skymap['irf'] = irf
skymap['caldb'] = caldb
skymap['emin'] = emin
skymap['emax'] = emax
skymap['usepnt'] = bool('yes')
skymap['nxpix'] = 200
skymap['nypix'] = 200
skymap['binsz'] = 0.02
skymap['coordsys'] = 'CEL'
skymap['proj'] = 'CAR'
skymap['bkgsubtract'] = 'IRF'
skymap['logfile'] = skymapName.replace('.fits', '.log')
skymap['debug'] = bool('no')
skymap.execute() 

print('!!! check --- skymaps: ', skymapName)

# ==============================
# !!! DETECTION AND MODELING !!!
# ==============================

detXml, detReg, pos = srcDetection_spcModeling(skymapName, sigma=sigma, maxSrc=10)

print('\n\n==========\n\n!!! check --- detection.............', texp,'s done\n\n!!! coords:', pos, '\n\n ==========\n\n')


# ==================
# !!! LIKELIHOOD !!!
# ==================


resultsName = event.replace(simpath, detpath).replace('.fits', '_results.xml' )

like = ctools.ctlike()
like['inobs'] = event
like['inmodel'] = detXml
like['outmodel'] = resultsName
like['caldb'] = caldb
like['irf'] = irf
like['fix_spat_for_ts'] = bool('yes')
like['logfile'] = resultsName.replace('.xml', '.log')
like['debug'] = bool('no')  # default
like.execute()


print('!!! check --- max likelihoods: ', resultsName)

# ===========================
# !!! LIKELIHOOD RA & DEC !!!
# ===========================


raFit = getRaDec(resultsName)[0]
decFit = getRaDec(resultsName)[1]


print('!!! check --- RA FIT: ', raFit)
print('!!! check --- DEC FIT: ', decFit)

# ==================
# !!! LIGHTCURVE !!!
# ==================

wbin = [1, 5, 10, 100]
nbin = []
for i in range(len(wbin)) :
  nbin.append(int(tmax/wbin[i]))

lc = []

for i in range(len(wbin)) :
  lc.append(resultsName.replace('results', 'lightcurve%ds_%ds' %(tmax, wbin[i])).replace('.xml', '.fits'))
  print('!!! check start LC n.', i+1, 'name:', lc[i])

  lightcurve = cscripts.cslightcrv()
  lightcurve['inobs'] = event
  lightcurve['inmodel'] = resultsName
  lightcurve['srcname'] = 'Src001'
  lightcurve['caldb'] = caldb
  lightcurve['irf'] = irf
  lightcurve['outfile'] = lc[i]
  lightcurve['tbinalg'] = 'LIN' # <FILE|LIN|GTI>
  lightcurve['edisp'] = bool('yes')
  lightcurve['tmin'] = tmin
  lightcurve['tmax'] = tmax
  lightcurve['tbins'] = nbin[i]
  lightcurve['method'] = '3D'
  lightcurve['emin'] = emin
  lightcurve['emax'] = emax
  lightcurve['enumbins'] = 0 # 0 for unbinned 3D only
  lightcurve['xref'] = float(raFit[0])
  lightcurve['yref'] = float(decFit[0])
  lightcurve['logfile'] = lc[i].replace('.fits', '.log')
  lightcurve['debug'] = bool('no') 
  lightcurve.execute()

print('!!! check ------- lc n.', i+1, ') completed:', lc[i])
