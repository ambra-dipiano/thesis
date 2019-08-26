#!/bin/python3.6

import sys
import os
from sys import argv

N = int(sys.argv[1])
trials = int(sys.argv[2])
start_count = int(sys.argv[3])
start_chunk = 0

task = ['run0406_task%03d.sh' % (start_chunk+i+1) for i in range(N)] 
jobs = 'job_run.sh'

for i in range(N):
  # write task ---!
  with open(task[i], 'w+') as f:
    f.write('#!/bin/bash')
    f.write('\n\n')
    f.write('#SBATCH -N1')
    f.write('\n#SBATCH --ntasks-per-node=1')
    f.write('\n#SBATCH --account=pianoambra@morgana')
    f.write('\n\n')
    f.write('python /mnt/nvme0n1p1/piano_analysis/working-dir/script_RTAdetection_v11.py %d %d %d' % (start_chunk+i+1, trials, start_count+trials*i))
    f.close()


# write jobs ---!
with open(jobs, 'w+') as f:
  f.write('#!/bin/bash')
  f.write('\n\n')
  for i in range(N):
    f.write('sbatch run0406_task%03d.sh\n' % (start_chunk+i+1))
  f.close()

os.system('chmod 777 job_run.sh')

