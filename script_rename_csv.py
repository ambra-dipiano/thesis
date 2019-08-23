#!/bin/python3.6

import os
import sys
from sys import argv

N=20
start = 1
texp = [1,5,10,100]
start_new = start+N

for i in range(N) :
  for j in range(len(texp)) :
    os.system('mv $WORK/run0406_bkg/run0406_ID000126/csv/run0406_v07_%ds_chunk%02d.csv $WORK/run0406_bkg/run0406_ID000126/csv/run0406_bkg_%ds_chunk%02d.csv' %(texp[j], start+i, texp[j], start_new+i))

