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

# for f in files:
#   with open(f, "w") as outfile:
#     for line in fileinput.input(f, inplace=False):
#       if fileinput.isfirstline():
#         outfile.write(['#trial,t exp,sigma,Ndet,RA Src001,DEC Src001\n'])
#       else:
#         outfile.write(line)

# for f in files:
#   with open(f) as file:
#       r = csv.reader(file)
#       w = csv.writer(file)
#       next(r, None)  # skip the first row from the reader, the old header
#       # write new header
#       w.writerow(['#trial,t exp,sigma,Ndet,RA Src001,DEC Src001\n'])
#       # copy the rest
#       for row in r:
#           w.writerow(row)

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