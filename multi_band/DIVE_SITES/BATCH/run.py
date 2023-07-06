import subprocess
import glob
import os
import re

window_size = 6366200 #3183100
dive_site_window_size = 12732
dx = 2

root_dir = os.getcwd()

dirs = glob.glob('0*')
print(dirs)

for d in dirs:

  print('Working in directory ', d)
  command = 'cp ./code/* ./' + d
  process = subprocess.Popen(command, shell=True)
  process.wait()

  os.chdir(d)

  f = open(d + '_dive_sites.txt')
  nodes = f.readline().strip('\n')
  nodes = max(1, int(nodes))
  f.close()

  print(d)
  f = open('submit_tmp')
  text = f.read()
  f.close()
  SF = re.compile('REPL1')
  SF2 = re.compile('REPL2')
  new = SF.sub('#PBS -l select=' + str(nodes) + ':ncpus=36',text)
  new2 = SF2.sub('./main -nl ' + str(nodes) + ' --in_name=' + d + ' --window_size=' + str(window_size) + ' --dive_site_window_size=' + str(dive_site_window_size) + ' --dx=' + str(dx),new)

  g = open('submit_cheyenne', 'w')
  g.write(new2)
  g.close()

  command = "rm submit_tmp"
  process = subprocess.Popen(command,shell=True)
  process.wait()

  command = "qsub submit_cheyenne"
  subprocess.Popen(command,shell=True)

  os.chdir(root_dir)

