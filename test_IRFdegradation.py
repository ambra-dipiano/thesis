from pkg_blindsearch import *

# files ---!
caldb = 'prod3b'
irf = 'South_z40_average_100s'

# path configuration ---!
cfg = xmlConfig()
p = ConfigureXml(cfg)

# initialise ---!
irfObj = Analysis()
irfObj.irf = irf
irfObj.caldb = caldb
irfObj.degradeIrf(bkg=True)