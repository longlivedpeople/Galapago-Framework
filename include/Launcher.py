import os
import time
from datetime import datetime

class Launcher:

    def __init__(self, script, env, ID = '', output = 'random', name = False):

        self.basename = name if name else 'Galapago'
        self.script = script
        self.ID = ID
        self.env = ''
        self.queue = ''

        if env == 'gridui':
            self.env = 'gridui' # deprecated (old gridui queue system)
        elif env == 'slurm':
            self.env = 'slurm'            
        elif env in ['espresso', 'microcentury', 'longlunch', 'workday', 'tomorrow', 'testmatch', 'nextweek']:
            self.env = 'condor'
            self.queue = env
        elif env == 'condor': # default (to be changed)
            self.env = 'condor'
            self.queue = 'microcentury' 


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
        self.bashfile = self.workpath + '_auxBashFile' + str(ID)+'.sh'
        self.condorsub = self.workpath + '_auxCondorSub' + str(ID)+'.sh'

        # Define queue (only in condorCONDOR)
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
            self.makeBashFile()
            self.makeCondorSubmitFile()
        if self.env == 'slurm':
            self.makeBashFile()


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

    def makeBashFile(self):

        templateBASH = """#!/bin/bash
pushd {0}
eval `scramv1 runtime -sh`
pushd
""".format(self.cmssw)

        templateBASH += 'python {0}'.format(self.auxscript)
        if self.output != '':
            templateBASH += ' --out {0}'.format(self.output)

        _f = open(self.bashfile, 'w')
        _f.write(templateBASH)
        _f.close()


    def makeCondorSubmitFile(self):

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
        _fs.write(templateCONDORsub.format(self.logs, self.queue, self.bashfile))
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
        elif self.env == 'condor':
            os.system('chmod +x ' + self.bashfile)
            os.system('chmod +x ' + self.condorsub)
            os.system('condor_submit ' + self.condorsub)
        elif self.env == 'slurm':
            os.system('chmod +x ' + self.bashfile)
            now = datetime.now()
            time_code = now.strftime("_%d%m_%H%M%S")
            os.system('sbatch -o ' + self.workpath +'logs/'+time_code+'.log -e '+self.workpath + 'logs/' + time_code + '.err --qos=gridui_sort --partition=cloudcms ' + self.bashfile)

    def clear(self):

        os.system('rm ' + self.auxsubmit)
        os.system('rm ' + self.auxscript)


