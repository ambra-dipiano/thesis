# MIT License
# Copyright (c) 2019, 2020 Ambra Di Piano
# ---------------------------------------
# ================== #
#  CALDB DEGRADATION #
# ================== #

from pkg_blindsearch import Analysis
import os

# files ---!
caldb = 'prod3b-v2'
irf = os.listdir(os.environ.get('CTOOLS') + '/share/caldb/data/cta/' + caldb + '/bcf/')
#print(irf)

# initialise ---!
for fits in irf:
  irfObj = Analysis()
  irfObj.irf = fits
  irfObj.caldb = caldb
  irfObj.degradeIrf()