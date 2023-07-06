import numpy
import subprocess
import re
import os
import glob
import time

root_dir = os.getcwd()

folders = glob.glob('0*')
folders.sort()

for fldr in folders:

  os.chdir(fldr)

  dsets = glob.glob(fldr + '*.nc')
  print(dsets)

  for dset in dsets:
    f = open('../code/run_compress.sh')
    text = f.read()
    f.close()

    SF1 = re.compile('REPL1')

    new = SF1.sub('mpiexec -n 1 ./compress.py -f ' + dset,text)

    f = open('tmp', 'w')
    f.write(new)
    f.close()

    command = 'qsub tmp'
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()
    time.sleep(1)

  os.chdir(root_dir)
