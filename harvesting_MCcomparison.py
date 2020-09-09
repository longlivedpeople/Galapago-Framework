import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time
import shutil

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager
#from include.Utils import *




def makeMCComparison(lumi, hnameA, hnameB, ylog, treeMC, treeDATA, inputdir, labelA = 'A', labelB = 'B', xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', DATAlabel = '', extralabel = ''):


    ### Get histograms
    luminosity = lumi
    SShnameA = hnameA.split('__')[0] + '_SS'
    SShnameB = hnameB.split('__')[0] + '_SS'

    hSS_A = treeDATA.getLoopTH1F(inputdir, SShnameA)
    hBKG_A = treeMC.getLoopTH1F(inputdir, hnameA)
    hSS_B = treeDATA.getLoopTH1F(inputdir, SShnameB)
    hBKG_B = treeMC.getLoopTH1F(inputdir, hnameB)
    
    hBKG_A.Add(hSS_A)
    hBKG_B.Add(hSS_B)


    ### Histogram tuning 
    hBKG_A.SetMarkerStyle(24)
    hBKG_A.SetMarkerSize(0.8)
    hBKG_A.SetMarkerColor(r.kRed)
    hBKG_B.SetMarkerStyle(24)
    hBKG_B.SetMarkerSize(0.8)
    hBKG_B.SetMarkerColor(r.kBlue)
    
    ### Normalized
    hBKG_A.Scale(1.0/hBKG_A.Integral())
    hBKG_B.Scale(1.0/hBKG_B.Integral())

    
    ### Get maximum
    maxValA = hBKG_A.GetMaximum()
    maxValB = hBKG_B.GetMaximum()
    maxVal = max(maxValA, maxValB)

    ### Set Maximum
    if not ylog:
        hBKG_A.SetMaximum(1.3*maxVal)
        hBKG_B.SetMaximum(1.3*maxVal)
    else:
        hBKG_A.SetMaximum(100.0*maxVal)
        hBKG_B.SetMaximum(100.0*maxVal)
        #hBKG_A.SetMinimum(0.1)
        #hBKG_B.SetMaximum(0.1)
 
    ### Canvas object
    plot = Canvas.Canvas('MCComparison_'+hnameA, 'png', 0.5, 0.79, 0.9, 0.87, 1)
    plot.addHisto(hBKG_A, 'P', labelA, 'p', r.kRed, 1, 0)
    plot.addHisto(hBKG_B, 'P, SAME', labelB, 'p', r.kBlue, 1, 1)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if extralabel:
        plot.addLatex(0.17, 0.76, extralabel, font = 42)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/MCComparisons_' + outtag + '/'
    plot.saveRatio(1, 0, ylog, luminosity, hBKG_A, hBKG_B, label="Ratio", outputDir = outdir)
    

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


    ############# EG data definition
    DoubleEGB = 'DoubleEG_Run2016B'
    DoubleEGC = 'DoubleEG_Run2016C'
    DoubleEGD = 'DoubleEG_Run2016D'
    DoubleEGE = 'DoubleEG_Run2016E'
    DoubleEGF = 'DoubleEG_Run2016F'
    DoubleEGG = 'DoubleEG_Run2016G'
    DoubleEGH = 'DoubleEG_Run2016H'

    DoubleEG_list = []
    DoubleEG_list.append(DoubleEGB)
    DoubleEG_list.append(DoubleEGC)
    DoubleEG_list.append(DoubleEGD)
    DoubleEG_list.append(DoubleEGE)
    DoubleEG_list.append(DoubleEGF)
    DoubleEG_list.append(DoubleEGG)
    DoubleEG_list.append(DoubleEGH)

    ############# Muon data definition
    DoubleMuonB = 'DoubleMuon_Run2016B'
    DoubleMuonC = 'DoubleMuon_Run2016C'
    DoubleMuonD = 'DoubleMuon_Run2016D'
    DoubleMuonE = 'DoubleMuon_Run2016E'
    DoubleMuonF = 'DoubleMuon_Run2016F'
    DoubleMuonG = 'DoubleMuon_Run2016G'
    DoubleMuonH = 'DoubleMuon_Run2016H'

    DoubleMuon_list = []
    DoubleMuon_list.append(DoubleMuonB)
    DoubleMuon_list.append(DoubleMuonC)
    DoubleMuon_list.append(DoubleMuonD)
    DoubleMuon_list.append(DoubleMuonE)
    DoubleMuon_list.append(DoubleMuonF)
    DoubleMuon_list.append(DoubleMuonG)
    DoubleMuon_list.append(DoubleMuonH)


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

    ############# Luminosity definition
    lumiB = 5.79
    lumiC = 2.57
    lumiD = 4.25
    lumiE = 4.01
    lumiF = 3.10
    lumiG = 7.54
    lumiH = 8.61

    lumiEG = {}
    lumiEG['DoubleEG_Run2016B'] = lumiB
    lumiEG['DoubleEG_Run2016C'] = lumiC
    lumiEG['DoubleEG_Run2016D'] = lumiD
    lumiEG['DoubleEG_Run2016E'] = lumiE
    lumiEG['DoubleEG_Run2016F'] = lumiF
    lumiEG['DoubleEG_Run2016G'] = lumiG
    lumiEG['DoubleEG_Run2016H'] = lumiH

    lumi_EG = 0.0
    for dataset in DoubleEG_list: lumi_EG += lumiEG[dataset]


    lumiMuon = {}
    lumiMuon['DoubleMuon_Run2016B'] = lumiB
    lumiMuon['DoubleMuon_Run2016C'] = lumiC
    lumiMuon['DoubleMuon_Run2016D'] = lumiD
    lumiMuon['DoubleMuon_Run2016E'] = lumiE
    lumiMuon['DoubleMuon_Run2016F'] = lumiF
    lumiMuon['DoubleMuon_Run2016G'] = lumiG
    lumiMuon['DoubleMuon_Run2016H'] = lumiH

    lumi_Muon = 0.0
    for dataset in DoubleMuon_list: lumi_Muon += lumiMuon[dataset]




    filename = 'dat/Samples_cern_filling.dat'

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )

    makeMCComparison(lumi = lumi_Muon, hnameA = 'hMM_trackIxy_highPU', hnameB = 'hMM_trackIxy_lowPU', ylog = True, treeMC = treeMC, treeDATA = treeDATA, inputdir = opts.input, labelA = 'High PU regime: nPV > 30', labelB = 'Low PU regime: nPV < 30', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', extralabel = '')
    
