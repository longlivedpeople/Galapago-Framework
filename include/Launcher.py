import os
import time

class Launcher:

    def __init__(self, script, env, ID = '', output = 'random', name = False):

        self.basename = name if name else 'Galapago'
        self.script = script
        self.ID = ID
        self.env = env


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

        # Define queue (only in condorCONDOR)
        self.queue = 'microcentury'
        self.logs = self.workpath + 'logs/'

        # Create logs/ folder
        if not os.path.exists(self.logs):
            os.makedirs(self.logs)

        self.output = output


        # Execute:
        self.makeSubmitScript()

        if self.env == 'gridui':
            self.makeGriduiSubmitFile()
        if self.env == 'condor': 
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

    def makeCondorSubmitFiles(self, submitfile = False):

        if not submitfile:
            submitfile = self.auxscript

        templateCONDOR = """#!/bin/bash
pushd {0}
eval `scramv1 runtime -sh`
pushd
""".format(self.cmssw)

        templateCONDOR += 'python {0}'.format(submitfile)
        if self.output != '':
            templateCONDOR += ' --out {0}'.format(self.output)

        _f = open(self.condorfile, 'w')
        _f.write(templateCONDOR)
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
        _fs.write(templateCONDORsub.format(self.logs, self.queue, self.condorfile))
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


    def addOrder(self, order, submitfile = False):

        if not submitfile:
            submitfile = self.auxscript

        with open(submitfile, 'a') as _file:
            if order[0] != ' ': 
                order = 4*' ' + order # correct indentation (provisional)
            _file.write(order)
            _file.close()


    def launch(self, submitfile = False):

        if not submitfile: 
            submitfile = self.auxsubmit

        if self.env == 'gridui':
            os.system('chmod +x ' + submitfile)
            os.system('qsub -o '+ self.workpath + 'logs/gridui'+str(self.ID)+'.log -e ' + self.workpath + 'logs/gridui'+str(self.ID)+'.err ' + submitfile)
        else:
            os.system('chmod +x ' + self.condorfile)
            os.system('chmod +x ' + self.condorsub)
            os.system('condor_submit ' + self.condorsub)


    def clear(self):

        os.system('rm ' + self.auxsubmit)
        os.system('rm ' + self.auxscript)


