import os
import time

class Launcher:

    def __init__(self, script, ID = '', output = 'withID', name = False):

        self.basename = name if name else 'Galapago'
        self.script = script
        self.ID = ID

        self.workpath = ''
        for level in script.split('/')[:-1]: self.workpath += level + '/'
        self.auxscript = self.workpath + '_auxScript'+str(ID)+'.py'
        self.auxsubmit = self.workpath + '_auxSubmit'+str(ID)+'.sh'

        self.output = output

        self.makeSubmitScript()
        self.makeGriduiSubmitFile(out = self.output)


    def makeGriduiSubmitFile(self, out):

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
        _auxSubmit.write(text.format(self.basename, self.auxscript, out))
        _auxSubmit.close()

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

        os.system('chmod +x ' + self.auxsubmit)
        os.system('qsub -o '+ self.workpath + 'logs/gridui'+str(self.ID)+'.log -e ' + self.workpath + 'logs/gridui'+str(self.ID)+'.err ' + self.auxsubmit)

    def clear(self):

        os.system('rm ' + self.auxsubmit)
        os.system('rm ' + self.auxscript)


