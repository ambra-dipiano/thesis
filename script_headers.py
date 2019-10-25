# ====================== #
# !!! CHANGE HEADERS !!! #
# ====================== #

import os
import fileinput
import csv

files = []
# r=root, d=directories, f = files
for r, d, f in os.walk(os.environ.get('MORGANA')+'/run0406_1e4detOnly_Ebl_noCut_Degr_noSort/run0406_det/run0406_ID000126/csv/'):
    for file in f:
        if '.csv' in file:
            files.append(os.path.join(r, file))

def line_prepender(filename, line):
  with open(filename, 'r+') as csvfile:
    next(csvfile)
    content = csvfile.read()
    csvfile.seek(0, 0)
    csvfile.write(line.rstrip('\r\n') + '\n' + content)

for f in files:
  header_list = ['#trials', 'texp', 'sigma', 'Ndet', 'RAdet', 'DECdet']
  hdr = ",".join(header_list)

  line_prepender(f, hdr)


print('done')