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




def makeBlindedPlot(lumi, hname_SR, hname_CR, ylog, treeDATA, inputdir, treeSI = False, hname_SI = '', rebin = False, limit = 0.0, lines = [], xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = ''):


    ### Get histograms
    luminosity = lumi

    hSR_ = treeDATA.getLoopTH1F(inputdir, hname_SR)
    hCR_ = treeDATA.getLoopTH1F(inputdir, hname_CR)
    
    ### rebinins:
    if type(rebin) != bool:
        if len(rebin) > 1:
            hSR = hSR_.Rebin(len(rebin)-1, hSR_.GetName() + '_rebined', rebin)
            hCR = hCR_.Rebin(len(rebin)-1, hCR_.GetName() + '_rebined', rebin)
        else:
            hSR = hSR_.Rebin(rebin)
            hCR = hCR_.Rebin(rebin)
    else:
        hSR = hSR_.Clone()
        hCR = hCR_.Clone()

    ### Initialize blinded histogram
    hSRblinded = hSR.Clone()

    Bover = 0.0
    for n in range(hSR.GetNbinsX() + 1):
        if limit and hSR.GetBinLowEdge(n) >= limit:
            hSRblinded.SetBinContent(n, 0.0)
            hSRblinded.SetBinError(n, 0.0)
            Bover += hCR.GetBinContent(n)
        else:
            hSRblinded.SetBinContent(n, hSR.GetBinContent(n))
            hSRblinded.SetBinError(n, hSR.GetBinError(n))
    print('B:', Bover)

    hCR.SetFillColorAlpha(r.kCyan-6, 0.8) 
    hCR.SetLineColor(r.kCyan-6) 
    hCR.GetXaxis().SetTitleSize(0.045)
    hCR.GetYaxis().SetTitleSize(0.045)
    hSRblinded.SetMarkerStyle(20)
    hSRblinded.SetMarkerSize(0.8)
    hSRblinded.SetMarkerColor(r.kBlack)

    ### Signal histograms
    s_histos = []
    if hname_SI == '': hname_SI = hname_SR
    if treeSI:

        hSIS = treeSI.getLoopStack(inputdir, hname_SI)

        for _i, _h in enumerate(hSIS.GetHists()):
            _h.Scale(lumi/35.87)
            s_histos.append(copy.deepcopy(_h))
            """
            Sover = 0.0
            for n in range(_h.GetNbinsX() + 1):
                if limit and _h.GetBinLowEdge(n) >= limit:
                    Sover += _h.GetBinContent(n)
            print(_h.GetTitle(), Sover, Sover/math.sqrt(Bover+Sover))
            """


    ### Get maximum
    maxValSR = hSRblinded.GetMaximum()
    maxValCR = hCR.GetMaximum()
    maxValSI = 0 if not treeSI else max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSR, maxValCR, maxValSI])

    ### Set Maximum
    if not ylog:
        hCR.SetMaximum(1.3*maxVal)
    else:
        hCR.SetMaximum(100.0*maxVal)
        hCR.SetMinimum(0.1)

    if ymax: hCR.SetMaximum(ymax) 


    ### Canvas object
    plot = Canvas.Canvas('SR_'+hname_SR, 'png', 0.35, 0.5, 0.6, 0.87, 1)
    plot.addHisto(hCR, 'HIST', 'Background (Data-driven)', 'f', '', 1, 0)
    
    ### Add signals:
    if treeSI:
        for i,_h in enumerate(s_histos):
            _h.SetLineWidth(2) # provisional
            masses = eval(_h.GetTitle()[3:])
            print(masses)
            legend = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
            plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    plot.addHisto(hSRblinded, 'P, SAME', 'Data', 'p', '', 1, len(s_histos)+1)
    for line in lines:
        plot.addLine(line, hCR.GetMinimum(), line, hCR.GetMaximum(), r.kBlack)

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)
    plot.addLatex(0.17, 0.85, extralabel, font = 62)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
    plot.saveRatio(1, 0, ylog, luminosity, hSRblinded, hCR, r_ymin = 0.9, r_ymax = 1.1, label="Data/Est.Bkg", outputDir = outdir)

    

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
    #DoubleMuon_list.append(DoubleMuonB)
    DoubleMuon_list.append(DoubleMuonC)
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
    Signals.append('HXX_400_50_4mm')
    Signals.append('HXX_400_50_40mm')
    Signals.append('HXX_400_50_400mm')
    Signals.append('HXX_400_150_400mm')
    Signals.append('HXX_1000_150_10mm')
    Signals.append('HXX_1000_150_100mm')
    Signals.append('HXX_1000_350_35mm')
    Signals.append('HXX_1000_350_350mm')

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




    filename = 'dat/Samples_cern_fillingv2.dat'


    ### Tree SI Tree
    treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Signals, 'SI'), name = 'SI', isdata = 0 )

    ################################
    ######## DoubleEG Plots ########
    ################################
       
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )

    #makeBlindedPlot(lumi = lumi_EG, hname_SR = 'hEEoffZSR_nBSEE', hname_CR = 'hEEoffZCR_nBSEE', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_withCounters', treeSI = treeSI, limit = False, lines = [0], xlabel = '', outtag = '', ymax = 1e13, LLlabel = 'EE', DATAlabel = '', extralabel = 'Signal region') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SR = 'hEEoffZSR_trackIxy', hname_CR = 'hEEoffZCR_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_newElectronVariables', treeSI = treeSI, limit = 6, lines = [6.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'Signal region') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SR = 'hEESROS_trackIxy', hname_CR = 'hEECROS_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_newElectronVariables', treeSI = treeSI, limit = 6, lines = [6.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'Signal region') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SR = 'hEEonZSR_mass', hname_CR = 'hEECROS_mass', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_newElectronVariables', treeSI = treeSI, hname_SI = 'hEESROS_mass',limit = 0.0, lines = [80.0, 100.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_EG, hname_SR = 'hEEdispSR_normalizedChi2', hname_CR = 'hEEdispCR_normalizedChi2', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_chi2test', treeSI = treeSI, hname_SI = 'hEEdispSR_normalizedChi2',limit = 0.0, lines = [10.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )

    #makeBlindedPlot(lumi = lumi_Muon, hname_SR = 'hMMoffZSR_nBSMM', hname_CR = 'hMMoffZCR_nBSMM', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_withCounters', treeSI = treeSI, limit = False, lines = [0], xlabel = '', outtag = '', ymax = 1e13, LLlabel = 'MM', DATAlabel = '', extralabel = 'Signal region') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SR = 'hMMoffZSR_trackIxy', hname_CR = 'hMMoffZCR_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_newElectronVariables', treeSI = treeSI, limit = 6, lines = [6.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'Signal region') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SR = 'hMMSROS_trackIxy', hname_CR = 'hMMCROS_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_newElectronVariables', treeSI = treeSI, limit = 6, lines = [6.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'Signal region') 

    #makeBlindedPlot(lumi = lumi_Muon, hname_SR = 'hMMonZSR_mass', hname_CR = 'hMMCROS_mass', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_newElectronVariables', treeSI = treeSI, hname_SI = 'hMMSROS_mass',limit = 0.0, lines = [80.0, 100.0], xlabel = '', outtag = '', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
