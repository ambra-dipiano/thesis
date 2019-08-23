#!/bin/python3.6

import os
import sys
from sys import argv

N=20
start = 4729454
stop = start+N

action=sys.argv[1]

if action=='rm' :
  for i in range(N) :
    os.system('rm slurm-%d.out' %(start+i))
elif action=='mv' :
  for i in range(N) :
    os.system('mv slurm-%d.out $WORK/run0406/.' %(start+i))
else :
  print('input either "rm" or "mv"')
