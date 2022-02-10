import os
import sys
import subprocess

path = '/eos/user/f/fernance/LLP_Analysis/UL/2017/EG_Run2017B/'

fullpaths = []

for era in 'BCDEF':
    for _f in os.listdir('/eos/user/f/fernance/LLP_Analysis/UL/2017/EG_Run2017'+era+'/'):
        fullpaths.append('/eos/user/f/fernance/LLP_Analysis/UL/2017/EG_Run2017'+era+'/' + _f)

print(fullpaths)   


for _f,f in enumerate(fullpaths):
    if not '.root' in f: continue
    #filename = path + f
    print('launching' + f)
    command = 'python /afs/cern.ch/work/f/fernance/private/Analysis-Utils/CMSSW_9_4_4/src/Analysis/CONDOR-Launcher/toCONDOR.py '
    command += 'microcentury ' 
    command += 'python $PWD/checkMass.py -t {0} -f {1}'.format(str(_f), f)
    os.system(command)
    #subprocess.Popen(["/bin/bash", "-i", "-c", 'toCONDOR espresso python $PWD/plotVertexPattern.py -t {0} -f {1}'.format(str(_f), filename)])
    #sp.communicate()
