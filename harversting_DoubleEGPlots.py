import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager
from include.Utils import *


################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C+')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()


    ############# Muon data definition
    DoubleEGB = ['DoubleEG_Run2016B']
    DoubleEGC = ['DoubleEG_Run2016C']
    DoubleEGD = ['DoubleEG_Run2016D']
    DoubleEGE = ['DoubleEG_Run2016E']
    DoubleEGF = ['DoubleEG_Run2016F']
    DoubleEGG = ['DoubleEG_Run2016G']
    DoubleEGH = ['DoubleEG_Run2016H']

    DoubleEG_list = []
    DoubleEG_list.append(DoubleEGB)
    DoubleEG_list.append(DoubleEGC)
    DoubleEG_list.append(DoubleEGD)
    DoubleEG_list.append(DoubleEGE)
    DoubleEG_list.append(DoubleEGF)
    DoubleEG_list.append(DoubleEGG)
    #DoubleEG_list.append(DoubleEGH)


    ############# Background definition
    Backgrounds = []
    Backgrounds.append('DYJetsToLL_M-50') 
    Backgrounds.append('DYJetsToLL_M-10to50') 
    Backgrounds.append('WW') 
    Backgrounds.append('WZ') 
    Backgrounds.append('ZZ') 
    Backgrounds.append('TT') 

    ############# Signal definition
    Signals = []
    #Signals.append('DisplacedSUSY_350_148_173')
    Signals.append('HXX_400_50_400mm')
    #Signals.append('HXX_400_50_40mm')
    #Signals.append('HXX_400_50_4mm')

    ############# Parameter definition
    lumiB = 5.79
    lumiC = 2.57
    lumiD = 4.25
    lumiE = 4.01
    lumiF = 3.10
    lumiG = 7.54
    lumiH = 8.61
    lumi_total = lumiB + lumiC + lumiD + lumiE + lumiF + lumiG + lumiH # luminosity
    lumi_list = [lumiB, lumiC, lumiD, lumiE, lumiF, lumiG, lumiH]
    

    filename = 'dat/Samples_cern_filling.dat'

    ############# Tree creation
    for i in range(0, len(DoubleEG_list)):
        dataset = DoubleEG_list[i]
        lumi = lumi_list[i]

        print(">>>>>> PLOT " + dataset[0] + " plots")

        treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
        treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Signals, 'SI'), name = 'SI', isdata = 0 )
        treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, dataset, 'DATA'), name = 'DATA', isdata = 1 )


        makeDataMCPlot(lumi, 'hEE_dPhi', True, treeMC, treeDATA, treeSI, WORKPATH + opts.input, outtag = dataset[0])
        makeDataMCPlot(lumi, 'hEE_mass', True, treeMC, treeDATA, treeSI, WORKPATH + opts.input, outtag = dataset[0])
        makeDataMCPlot(lumi, 'hEE_trackIxy', True, treeMC, treeDATA, treeSI, WORKPATH + opts.input, outtag = dataset[0])
        makeDataMCPlot(lumi, 'hEE_trackDxy', True, treeMC, treeDATA, treeSI, WORKPATH + opts.input, outtag = dataset[0])
        makeDataMCPlot(lumi, 'hEE_Lxy', True, treeMC, treeDATA, treeSI, WORKPATH + opts.input, outtag = dataset[0])
        makeDataMCPlot(lumi, 'hEE_Ixy', True, treeMC, treeDATA, treeSI, WORKPATH + opts.input, outtag = dataset[0])




