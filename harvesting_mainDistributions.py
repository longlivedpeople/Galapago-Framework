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




def makeDataMCPlot(lumi, hname, hnameDATA, ylog, treeMC, treeDATA, inputdir, treeSI = False, xlabel = '', outtag = '', yshift = 100.0, LLlabel = '', DATAlabel = '', extralabel = '', xlog = False):


    
    luminosity = lumi

    hSS = treeDATA.getLoopTH1F(inputdir, hnameDATA)
    hOS = treeMC.getLoopStack(inputdir, hname)

    # hSS tunning:
    hSS.SetLineColor(r.kBlack)
    hSS.SetFillColor(r.kMagenta+1)
    hSS.SetTitle('QCD, W+jets')
    #print(hSS.GetEntries())

    ### Combine SS + OS contributions:
    hBKG = r.THStack('hBKG_%s'%(hname), ';'+hOS.GetXaxis().GetTitle()+';'+hOS.GetYaxis().GetTitle()) # Background stacked
    hBKGtotal = copy.deepcopy(hSS) # Background total (for ratio)
    hBKG.Add(copy.deepcopy(hSS))

    for _h in hOS.GetHists():
        _h.Scale(lumi/35.87)
        hBKG.Add(copy.deepcopy(_h))
        hBKGtotal.Add(copy.deepcopy(_h))


    nBCK = len(hOS.GetHists()) + 1
    hBKGtotal.SetMarkerStyle(20) # Auxiliar to save the ratio correctly 
    
    
    ### Signal histograms
    s_histos = []
    if treeSI:
    
        hSIS = treeSI.getLoopStack(inputdir, hname)
    
        for _i, _h in enumerate(hSIS.GetHists()):
            s_histos.append(copy.deepcopy(_h))
   
        
    ### Data histogram
    hDATA = treeDATA.getLoopTH1F(inputdir, hname)
    hDATA.SetMarkerStyle(20)
    hDATA.SetMarkerSize(0.8)
        
    
    ### Get maximum
    maxValMC = hBKGtotal.GetMaximum()
    maxValSI = 0 if not treeSI else max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxValDATA = 0 if not treeDATA else hDATA.GetMaximum()
    maxVal = max([maxValMC, maxValSI, maxValDATA])
    
        ### Set Maximum
    if not ylog:
        hBKG.SetMaximum(1.3*maxVal)
        hBKG.SetMinimum(0.0)
        hBKGtotal.SetMaximum(1.3*maxVal)
        hBKGtotal.SetMinimum(0.0)
        if treeSI:
            for _h in s_histos: 
                if not yshift:
                    _h.SetMaximum(1.3*maxVal)
                else:
                    _h.SetMaximum(yshift*maxVal)
                _h.SetMinimum(0.0)
        if treeDATA:
            if not yshift:
                hDATA.SetMaximum(1.3*maxVal)
            else:
                hDATA.SetMaximum(yshift*maxVal)
            hDATA.SetMinimum(0.0)
    else:
        if not yshift:
            hBKG.SetMaximum(10.0*maxVal)
        else:
            hBKG.SetMaximum(yshift*maxVal)
        hBKG.SetMinimum(0.1)
        if not yshift:
            hBKGtotal.SetMaximum(10.0*maxVal)
        else:
            hBKGtotal.SetMaximum(yshift*maxVal)
        hBKGtotal.SetMinimum(0.1)
        if treeSI:
            for _h in s_histos: 
                if not yshift:
                    _h.SetMaximum(10.0*maxVal)
                else:
                    _h.SetMaximum(yshift*maxVal)
                _h.SetMinimum(0.1)
        if treeDATA:
            if not yshift:
                hDATA.SetMaximum(10.0*maxVal)
            else:
                hDATA.SetMaximum(yshift*maxVal)
            hDATA.SetMinimum(0.1)
    
    
    ### -> Canvas object
    if treeDATA:
        plot = Canvas.Canvas(hname, 'png', 0.6, 0.65, 0.85, 0.9, 1)
    else:
        plot = Canvas.Canvas(hname, 'png', 0.6, 0.79, 0.9, 0.9, 1)
    
    ### Add background:
    plot.addStack(hBKG, 'HIST', 1, 0) # Background
    
    ### Add signals:
    if treeSI:
        for i,_h in enumerate(s_histos):
            _h.SetLineWidth(2) # provisional
            plot.addHisto(_h, 'HIST, SAME', _h.GetTitle(), 'l', _h.GetFillColor(), 1, i+nBCK) # Signal
    
    ### Add DATA:
    if treeDATA:
        plot.addHisto(hDATA, 'P, SAME', 'Data', 'p', r.kBlack, 1, nBCK + len(s_histos))
    
    ### DATAlabel banner:
    if DATAlabel:
        plot.addLatex(0.17, 0.85, DATAlabel, font = 62)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if extralabel:
        plot.addLatex(0.17, 0.76, extralabel, font = 42)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/Plots_' + outtag + '/'
    if treeDATA:
        plot.saveRatio(1, 0, ylog, luminosity, hDATA, hBKGtotal, label="Data/BKG", outputDir = outdir, xlog = xlog)
    else:
        plot.save(1, 0, ylog, luminosity, '', outputDir = outdir)
    

################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-e', '--EGinput', action='store', type=str, dest='EGinput', default='', help='Target directory')
    parser.add_option('-m', '--Muoninput', action='store', type=str, dest='Muoninput', default='', help='Target directory')
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
    Signals.append('HXX_1000_150_100mm')
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
    """ 
    for i in range(0, len(DoubleEG_list)):
        dataset = DoubleEG_list[i]
        lumi = lumiEG[dataset]

        print(">>>>>> PLOT " + dataset + " plots")

        treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
        #treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Signals, 'SI'), name = 'SI', isdata = 0 )
        treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, [dataset], 'DATA'), name = 'DATA', isdata = 1 )


        makeDataMCPlot(lumiEG[dataset], 'hEE_dPhi', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_mass', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2' )
        makeDataMCPlot(lumiEG[dataset], 'hEE_trackIxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_trackDxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_Lxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_Ixy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_trackIxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_trackDxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_Lxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_Ixy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_leadingPt', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_subleadingPt', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hEE_normalizedChi2', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hE_eta', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hE_dxy', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hE_dxyError', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiEG[dataset], 'hE_charge', True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = dataset, LLlabel = 'EE', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
    """

    ##### Joint luminosity
    
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )

    makeDataMCPlot(lumi_EG, 'hEEOS0_dPhi', 'hEESS0_dPhi',True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = 'DoubleEG_Full2016', LLlabel = 'EE', DATAlabel = 'DoubleEG_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', yshift = 10000.0)
    makeDataMCPlot(lumi_EG, 'hEEOS0_mass', 'hEESS0_mass',True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = 'DoubleEG_Full2016', LLlabel = 'EE', DATAlabel = 'DoubleEG_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_EG, 'hEEOS0_trackIxy', 'hEESS0_trackIxy',True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = 'DoubleEG_Full2016', LLlabel = 'EE', DATAlabel = 'DoubleEG_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_EG, 'hEEOS0_Ixy', 'hEESS0_Ixy',True, treeMC, treeDATA, WORKPATH + opts.EGinput, outtag = 'DoubleEG_Full2016', LLlabel = 'EE', DATAlabel = 'DoubleEG_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    
    

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    """
    for i in range(0, len(DoubleMuon_list)):
        dataset = DoubleMuon_list[i]
        lumi = lumiMuon[dataset]

        print(">>>>>> PLOT " + dataset + " plots")

        treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
        treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Signals, 'SI'), name = 'SI', isdata = 0 )
        treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, [dataset], 'DATA'), name = 'DATA', isdata = 1 )


        makeDataMCPlot(lumiMuon[dataset], 'hMM_dPhi', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_mass', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_cosAlpha', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_trackIxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_trackDxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_Lxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_Ixy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_leadingPt', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_subleadingPt', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hMM_normalizedChi2', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hM_eta', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hM_dxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hM_dxyError', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hM_charge', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hM_numberOfValidHits', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
        makeDataMCPlot(lumiMuon[dataset], 'hM_normChi2', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = dataset, LLlabel = 'MM', DATAlabel = dataset, extralabel = '|Control region #Delta#Phi| > #pi/2')
    """
    
    ## Full luminosity
    """
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )

    makeDataMCPlot(lumi_Muon, 'hMM_nPU_weighted', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_nPU_unweighted', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_dPhi', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', yshift = 10000.0)
    makeDataMCPlot(lumi_Muon, 'hMM_mass', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_cosAlpha', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', yshift = 10000.0)
    makeDataMCPlot(lumi_Muon, 'hMM_trackIxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_trackDxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_Lxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_Ixy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_trackIxy_log', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', xlog = True)
    makeDataMCPlot(lumi_Muon, 'hMM_trackDxy_log', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', xlog = True)
    makeDataMCPlot(lumi_Muon, 'hMM_Lxy_log', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', xlog = True)
    makeDataMCPlot(lumi_Muon, 'hMM_Ixy_log', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', xlog = True)
    makeDataMCPlot(lumi_Muon, 'hMM_leadingPt', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_subleadingPt', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hMM_normalizedChi2', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hM_eta', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', yshift = 10000.0)
    makeDataMCPlot(lumi_Muon, 'hM_dxy', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hM_dxyError', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    makeDataMCPlot(lumi_Muon, 'hM_charge', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', yshift = 10000.0)
    makeDataMCPlot(lumi_Muon, 'hM_numberOfValidHits', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2', yshift = 1000.0)
    makeDataMCPlot(lumi_Muon, 'hM_normChi2', True, treeMC, treeDATA, WORKPATH + opts.Muoninput, outtag = 'DoubleMuon_Full2016', LLlabel = 'MM', DATAlabel = 'DoubleMuon_Full2016', extralabel = '|Control region #Delta#Phi| > #pi/2')
    """


