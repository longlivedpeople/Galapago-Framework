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




def makeSignalBackgroundPlot(hname, ylog, treeMC, treeDATA, treeSI, inputdir, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', DATAlabel = '', extralabel = ''):


    SShname = hname.split('__')[0] + '_SS'

    hSS = treeDATA.getLoopTH1F(inputdir, SShname)
    hBKG = treeMC.getLoopTH1F(inputdir, hname)
    hBKG.Add(hSS)

    
    ### Signal histograms
    s_histos = []
    if treeSI:
    
        hSIS = treeSI.getLoopStack(inputdir, hname)
    
        for _i, _h in enumerate(hSIS.GetHists()):
            _h.Scale(1.0/_h.Integral())
            s_histos.append(copy.deepcopy(_h))
   
    hBKG.Scale(1.0/hBKG.Integral())


    ### Get maximum
    maxBKG = hBKG.GetMaximum()
    maxSI = 0 if not treeSI else max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxBKG, maxSI])
    
        ### Set Maximum
    if not ylog:
        hBKG.SetMaximum(1.15*maxVal)
        hBKG.SetMinimum(0.0)
        if treeSI:
            for _h in s_histos: 
                if not yshift:
                    _h.SetMaximum(1.15*maxVal)
                else:
                    _h.SetMaximum(yshift*maxVal)
                _h.SetMinimum(0.0)
    else:
        if not yshift:
            hBKG.SetMaximum(10.0*maxVal)
        else:
            hBKG.SetMaximum(yshift*maxVal)
        hBKG.SetMinimum(0.1)
        if treeSI:
            for _h in s_histos: 
                if not yshift:
                    _h.SetMaximum(10.0*maxVal)
                else:
                    _h.SetMaximum(yshift*maxVal)
                _h.SetMinimum(0.1)
    
    
    ### -> Canvas object
    plot = Canvas.Canvas(hname, 'png', 0.55, 0.69, 0.87, 0.87, 1)
    
    ### Add background:
    hBKG.SetLineColor(r.kAzure-7)
    hBKG.SetFillColorAlpha(r.kAzure, 0.5)
    hBKG.GetXaxis().SetLabelSize(0.04)
    hBKG.GetXaxis().SetTitleSize(0.045)
    hBKG.GetYaxis().SetLabelSize(0.04)
    hBKG.GetYaxis().SetTitleSize(0.045)
    #hBKG.SetFillStyle(3016)
    plot.addHisto(hBKG, 'HIST', 'Background', 'f', '', 1, 0) # Background
    
    ### Add signals:
    if treeSI:
        for i,_h in enumerate(s_histos):
            _h.SetLineWidth(2) # provisional
            plot.addHisto(_h, 'HIST, SAME', _h.GetTitle(), 'l', _h.GetFillColor(), 1, i+1) # Signal
    

    plot.addLatex(0.17, 0.85, 'Normalized to 1', font = 62)
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if extralabel:
        plot.addLatex(0.17, 0.76, extralabel, font = 42)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SIvsBKG_' + outtag + '/'
    plot.save(1, 0, ylog, '', '', outputDir = outdir)
    

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
    #DoubleMuon_list.append(DoubleMuonC)
    DoubleMuon_list.append(DoubleMuonD)
    #DoubleMuon_list.append(DoubleMuonE)
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
    Signals.append('HXX_1000_350_350mm')
    Signals.append('HXX_1000_150_100mm')
    Signals.append('HXX_400_150_400mm')
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

    ################################
    ######## DoubleEG Plots ########
    ################################
       
    ## Full luminosity
    treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename,Signals, 'MC'), name = 'SI', isdata = 0 )
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
    treeMuon = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )
    treeEG = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )



    makeSignalBackgroundPlot(hname = 'hMM_dPhi_full', ylog = False, treeMC = treeMC, treeDATA = treeMuon, treeSI = treeSI, inputdir = opts.input, xlabel = '', outtag = 'SIBKG_comp', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '')
    makeSignalBackgroundPlot(hname = 'hEE_dPhi_full', ylog = False, treeMC = treeMC, treeDATA = treeEG, treeSI = treeSI, inputdir = opts.input, xlabel = '', outtag = 'SIBKG_comp', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '')
