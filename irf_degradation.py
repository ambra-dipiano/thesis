from pkg_blindsearch import Analysis
import os

# files ---!
caldb = 'prod3b'
irf = os.listdir(os.environ.get('CTOOLS') + '/share/caldb/data/cta/' + caldb + '/bcf/')
#print(irf)

# initialise ---!
for fits in irf:
  irfObj = Analysis()
  irfObj.irf = fits
  irfObj.caldb = caldb
  irfObj.degradeIrf(bkg=True)