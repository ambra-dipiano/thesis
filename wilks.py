# ----------------------------- #
# Copyright 2020 Ambra Di Piano #
# ----------------------------- # -------------------------------------------------- #
# Redistribution and use in source and binary forms, with or without modification,   #
# are permitted provided that the following conditions are met:                      #
# 1. Redistributions of source code must retain the above copyright notice,          #
# this list of conditions and the following disclaimer.                              #
# 2. Redistributions in binary form must reproduce the above copyright notice,       #
# this list of conditions and the following disclaimer in the documentation and/or   #
# other materials provided with the distribution.                                    #
# 3. Neither the name of the copyright holder nor the names of its contributors      #
# may be used to endorse or promote products derived from this software without      #
# specific prior written permission.                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. #
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,   #
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,     #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,      #
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE    #
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  #
# OF THE POSSIBILITY OF SUCH DAMAGE.                                                 #
# ---------------------------------------------------------------------------------- #
# ===================================
# !!! SET UP PATH AND MERGE FILES !!!
# ===================================

import os
import sys
sys.path.append('/home/ambra/Desktop/cluster-morgana/')
from module_statistics import *

dof = 1
folder = 'tesi_bkg_1e6_nominal_%ddof_fullE' % dof
path = '/home/ambra/Desktop/cluster-morgana/archive_tests/tesi_02_bkg/' + folder + '/run0406_bkg/run0406_ID000126/csv/'
png_path = '/home/ambra/Desktop/cluster-morgana/archive_tests/tesi_02_bkg/' + folder + '/png/'
if not os.path.isdir(png_path):
  os.mkdir(png_path)

Nchunk = 20

texp = [1, 5, 10, 100]
sigma = [5]
chunk = [i + 1 for i in range(Nchunk)]

# csvName[texp][chunk]
csvName = [[] * i for i in range(len(texp))]
for i in range(len(chunk)):
  for j in range(len(texp)):
    csvName[j].append('bkg_%ds_chunk%02d.csv' % (texp[j], chunk[i]))

# merge files ---!
csvMerged = []
for j in range(len(texp)):
  csvMerged.append('bkg_%ds.csv' % texp[j])

  fout = open(path + csvMerged[j], 'w+')
  # first file ---!
  for line in open(path + csvName[j][0]):
    fout.write(line)
  # remaining files ---!
  for i in range(len(chunk) - 1):
    f = open(path + csvName[j][i + 1])
    next(f)  # skip the header ---!
    for line in f:
      fout.write(line)
    f.close()
  fout.close()

print('data files merge completed')

show = False
fontsize = 18
for n in range(len(texp)):
  filename = csvMerged[n]
  print('!======== texp = ', texp[n], '(s) ========!')
  # load DataFrame and column names ---!
  df = pd.read_csv(path + filename)
  cols = list(df.columns)
  trials = len(df[cols[0]])
  print('verify trials = ', trials)
  # drop duplicates ---!
  df.sort_values(cols[0], inplace=True)
  # dropping ALL duplicte values
  df.drop_duplicates(subset=cols[0], keep='last', inplace=True)
  trials = len(df[cols[0]])
  print('verify trials = ', trials)
  # drop NaN ---!
  # df = df.dropna()

  # set arrays ---!
  trial = np.array(df[cols[0]])
  tsv = np.array(df[cols[-1]])

  tsv.sort()

  wbin = 1
  nbin = int(tsv.max() / wbin)
  if nbin == 0:
    nbin = 1
  print('ts bin:', nbin)

  # -------------------------------- STATS ---!

  ts = []
  for i in range(trials):
    if tsv[i] < 0.0 or tsv[i] == np.nan:
      ts.append(0.0)
    else:
      ts.append(tsv[i])

  # chi2, chi2r = chi2_reduced(ts, trials, df=dof, nbin=nbin, width=wbin, var=False)
  # print('var=False; chi2=', chi2, '; chi2r=', chi2r)

  # -------------------------------- PLOT ---!

  fig, ax = ts_wilks(ts, trials, df=dof, nbin=nbin, width=wbin, figsize=(8, 12), fontsize=fontsize,
                     title='prod3b-v2: South\_z40\_0.5h (texp=%ds)' % texp[n], show=True,
                     filename=png_path + filename.replace('.csv', '_wilks.png'), ylim=(1e-7, 2e0))

  fig, ax = p_values(ts, trials, df=dof, nbin=nbin, width=wbin, figsize=(8, 12), fontsize=fontsize,
                     title='prod3b-v2: South\_z40\_0.5h (texp=%ds)' % texp[n], show=True,
                     filename=png_path + filename.replace('.csv', '_pvalues.png'), ylim=(1e-7, 2e0))

  fig, ax = ts_wilks_cumulative(ts, trials, df=dof, nbin=nbin, width=wbin, figsize=(12, 8),
                                fontsize=fontsize,
                                title='prod3b-v2: South\_z40\_0.5h (texp=%ds)' % texp[n],
                                filename=png_path + filename.replace('.csv', '_cumulative.png'))



