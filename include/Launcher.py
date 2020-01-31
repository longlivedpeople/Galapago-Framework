import os
import time

class Launcher:

    def __init__(self, script, ID = '', output = 'withID', name = False, gridui = True):

        self.basename = name if name else 'Galapago'
        self.script = script
        self.ID = ID
        self.gridui = gridui

        # Identify workpath:
        self.workpath = ''
        for level in script.split('/')[:-1]: self.workpath += level + '/'

        # Identify cmssw release
        self.cmssw = ''
        for level in script.split('/'): 
            self.cmssw += level + '/'
            if 'CMSSW' in level:
                self.cmssw += 'src' 
                break

        # Define file names for job submission
        self.auxscript = self.workpath + '_auxScript'+str(ID)+'.py'
        self.auxsubmit = self.workpath + '_auxSubmit'+str(ID)+'.sh'
        self.condorfile = self.workpath + '_auxCondorFile' + str(ID)+'.sh'
        self.condorsub = self.workpath + '_auxCondorSub' + str(ID)+'.sh'

        # Define qeue (CONDOR)
        self.qeue = 'espresso'
        self.logs = self.workpath + 'logs/'

        self.output = output


        # Execute:
        self.makeSubmitScript()

        if self.gridui:
            self.makeGriduiSubmitFile()
        else: # default condor
            self.makeCondorSubmitFiles()


    def makeGriduiSubmitFile(self):

        text = """
###################
# FLAG definition #
###################

# Request the Bourne Shell
#$ -S /bin/bash

# Change to the current working directory before starting the job
#$ -cwd

# Change the job name to "hello_world"

#$ -N {0}

# Resource request. We request 1MB of memory, and 60 seconds of wall
# clock time, that more than is enough for the test.
# -l mem_free=2G
# -l h_rt=4:0:00

# We are using the "l.tests" project for the examples
#$ -P l.gaes


#################
# Actual script #
#################

pushd /gpfs/users/fernanc/CMSSW_9_4_4/src/
eval `scramv1 runtime -sh`
pushd
python {1} --out {2}
"""
        _auxSubmit = open(self.auxsubmit, 'w')
        _auxSubmit.write(text.format(self.basename, self.auxscript, self.output))
        _auxSubmit.close()

    def makeCondorSubmitFiles(self):

        templateCONDOR = """#!/bin/bash
pushd {0}
eval `scramv1 runtime -sh`
pushd
python {1} --out {2}
"""

        _f = open(self.condorfile, 'w')
        _f.write(templateCONDOR.format(self.cmssw, self.auxscript, self.output))
        _f.close()


        templateCONDORsub = """
universe                = vanilla
executable              = $(filename)
output                  = {0}$(ClusterId).$(ProcId).out
error                   = {0}$(ClusterId).$(ProcId).err
log                     = {0}$(ClusterId).log
Notify_user             = fernance@cern.ch
+JobFlavour = "{1}" 
queue filename matching {2}
"""

        _fs = open(self.condorsub, 'w')
        _fs.write(templateCONDORsub.format(self.logs, self.qeue, self.condorfile))
        _fs.close()


    def makeSubmitScript(self):

        original = open(self.script, 'r')
        original_lines = original.readlines()
        copy = open(self.auxscript, 'w')

        for line in original_lines:
            if '############# Plotting' in line: break
            copy.write(line)

        original.close()
        copy.close()


    def addOrder(self, order):

        with open(self.auxscript, 'a') as _file:
            _file.write(4*' ' + order)
            _file.close()


    def launch(self):

        if self.gridui:
            os.system('chmod +x ' + self.auxsubmit)
            os.system('qsub -o '+ self.workpath + 'logs/gridui'+str(self.ID)+'.log -e ' + self.workpath + 'logs/gridui'+str(self.ID)+'.err ' + self.auxsubmit)
        else:
            os.system('chmod +x ' + self.condorfile)
            os.system('chmod +x ' + self.condorsub)
            os.system('condor_submit ' + self.condorsub)


    def clear(self):

        os.system('rm ' + self.auxsubmit)
        os.system('rm ' + self.auxscript)


