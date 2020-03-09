import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect
from string import Template

import include.Sample as Sample
import include.helper as helper
import include.Canvas as Canvas

LUMI = 36.773
OUTPATH = os.path.abspath('./') + '/launchWithGridui/outputHistograms'
WORKDIR = os.path.abspath('./') + '/'

def createSubmissionFile(DATFILE, SAMPLE, OPTION):

    parameters = {}
    parameters['LUMI'] = LUMI
    parameters['OUTPATH'] = OUTPATH
    parameters['DATFILE'] = DATFILE
    parameters['SAMPLE'] = SAMPLE
    parameters['OPTION'] = OPTION
    
    temp_file = open('launchWithGridui/TEMPLATE_submitJob.txt', 'r')
    input_text = temp_file.read()
    output_text = Template(input_text).safe_substitute(parameters)    
    temp_file.close()

    sub_file = open('launchWithGridui/submit_' + SAMPLE + '.sh', 'w')
    sub_file.write(output_text)
    sub_file.close()



def getSubmissionCommand(SAMPLE):

    COMMAND = "chmod +x " + WORKDIR + 'launchWithGridui/submit_' + SAMPLE + '.sh'
    COMMAND = COMMAND + '\n'
    COMMAND = COMMAND + 'qsub -o ' + WORKDIR + 'launchWithGridui/logs/' + SAMPLE + '.log '
    COMMAND = COMMAND + '-e ' + WORKDIR + 'launchWithGridui/logs/' + SAMPLE + '.err '
    COMMAND = COMMAND + WORKDIR + 'launchWithGridui/submit_' + SAMPLE + '.sh'
    COMMAND = COMMAND + '\n'

    return COMMAND


if __name__=='__main__':


    ### Load the different trees:
    treeMC = Sample.Tree("dat/MC.dat", "MC", 0, "aux_MC.root")
    treeSI = Sample.Tree("dat/SI.dat", "SI", 0, "aux_SI.root")
    treeDATA = Sample.Tree("dat/DATA.dat", "DATA", 0, "aux_DATA.root")

    ### Create the run.sh file
    shFile = open("launchWithGridui/runCode.sh", "w")


    ### Loop over MC samples:
    for b in treeMC.blocks:
        for s in b.samples:
            createSubmissionFile("dat/MC.dat", s.name, "MC")
            command = getSubmissionCommand(s.name)
            shFile.write(command)

    ### Loop over SI samples:
    for b in treeSI.blocks:
        for s in b.samples:
            createSubmissionFile("dat/SI.dat", s.name, "SI")
            command = getSubmissionCommand(s.name)
            shFile.write(command)


    ### Loop over DATA samples:
    for b in treeDATA.blocks:
        for s in b.samples:
            createSubmissionFile("dat/DATA.dat", s.name, "DATA")
            command = getSubmissionCommand(s.name)
            shFile.write(command)

    shFile.close()

