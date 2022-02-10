import os
import sys
import subprocess

path = '/eos/user/f/fernance/LLP_Analysis/UL/2018/TTTo2L2Nu/'
for _f,f in enumerate(os.listdir(path)):
    if not '.root' in f: continue
    filename = path + f
    print('launching' + filename)
    command = 'python /afs/cern.ch/work/f/fernance/private/Analysis-Utils/CMSSW_9_4_4/src/Analysis/CONDOR-Launcher/toCONDOR.py '
    command += 'microcentury ' 
    command += 'python $PWD/plotVertexPattern.py -t {0} -f {1}'.format(str(_f), filename)
    os.system(command)
    #subprocess.Popen(["/bin/bash", "-i", "-c", 'toCONDOR espresso python $PWD/plotVertexPattern.py -t {0} -f {1}'.format(str(_f), filename)])
    #sp.communicate()
