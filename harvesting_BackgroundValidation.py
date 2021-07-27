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




def makePromptBKGPlot(name, lumi, hname_SR, hname_CR, ylog, treeDATA, inputdir, rebin = False, limit = 0.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', xlog = False):


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

    ### Blinding limits:
    if limit:
        for n in range(1, hSR.GetNbinsX()):
            if hSR.GetBinLowEdge(n) > limit: hSR.SetBinContent(n, 0.0)
            #if hCR_.GetBinLowEdge(n) > limit: hCR_.SetBinContent(n, 0.0)

    hCR.SetFillColorAlpha(r.kCyan-6, 0.8) 
    hCR.SetLineColor(r.kCyan+3) 
    hSR.SetMarkerStyle(20)
    hSR.SetMarkerSize(0.8)
    hSR.SetMarkerColor(r.kBlack)
    hCR.GetXaxis().SetTitleSize(0.045)
    hSR.GetXaxis().SetTitleSize(0.045)
    hCR.GetYaxis().SetTitleSize(0.045)
    hSR.GetYaxis().SetTitleSize(0.045)

    ### Get maximum
    maxValSR = hSR.GetMaximum()
    maxValCR = hCR.GetMaximum()
    maxVal = max([maxValSR, maxValCR])

    ### Set Maximum
    if not ylog:
        hCR.SetMaximum(1.3*maxVal)
    else:
        hCR.SetMaximum(100.0*maxVal)
        hCR.SetMinimum(0.1)
 

    ### Canvas object
    plot = Canvas.Canvas('BKGVal_'+name, 'png', 0.51, 0.75, 0.7, 0.83, 1)
    plot.addHisto(hCR, 'HIST', 'Background (Data-driven)', 'f', '', 1, 0)
    plot.addHisto(hSR, 'P, SAME', 'Data', 'p', '', 1, 1)
    plot.addLatex(0.17, 0.8, extralabel, font = 62)

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.74, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.74, '#mu^{+}#mu^{-} channel', font = 42)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/BKGVal_' + outtag + '/'
    plot.saveRatio(1, 0, ylog, luminosity, hSR, hCR, r_ymin = 0.8, r_ymax = 1.2, label="Obs./Pred.", outputDir = outdir, xlog = xlog)

    

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
    DoubleEGB = 'DoubleEG_Run2016B_HIPM'
    DoubleEGC = 'DoubleEG_Run2016C_HIPM'
    DoubleEGD = 'DoubleEG_Run2016D_HIPM'
    DoubleEGE = 'DoubleEG_Run2016E_HIPM'
    DoubleEGF1 = 'DoubleEG_Run2016F_HIPM'
    DoubleEGF2 = 'DoubleEG_Run2016F_noHIPM'
    DoubleEGG = 'DoubleEG_Run2016G_noHIPM'
    DoubleEGH = 'DoubleEG_Run2016H_noHIPM'

    DoubleEG_HIPM = []
    DoubleEG_noHIPM = []
    DoubleEG_HIPM.append(DoubleEGB)
    DoubleEG_HIPM.append(DoubleEGC)
    DoubleEG_HIPM.append(DoubleEGD)
    DoubleEG_HIPM.append(DoubleEGE)
    DoubleEG_HIPM.append(DoubleEGF1)
    DoubleEG_noHIPM.append(DoubleEGF2)
    DoubleEG_noHIPM.append(DoubleEGG)
    DoubleEG_noHIPM.append(DoubleEGH)

    ############# Muon data definition
    DoubleMuonB = 'DoubleMuon_Run2016B_HIPM'
    DoubleMuonC = 'DoubleMuon_Run2016C_HIPM'
    DoubleMuonD = 'DoubleMuon_Run2016D_HIPM'
    DoubleMuonE = 'DoubleMuon_Run2016E_HIPM'
    DoubleMuonF1 = 'DoubleMuon_Run2016F_HIPM'
    DoubleMuonF2 = 'DoubleMuon_Run2016F_noHIPM'
    DoubleMuonG = 'DoubleMuon_Run2016G_noHIPM'
    DoubleMuonH = 'DoubleMuon_Run2016H_noHIPM'

    DoubleMuon_HIPM = []
    DoubleMuon_noHIPM = []
    DoubleMuon_HIPM.append(DoubleMuonB)
    DoubleMuon_HIPM.append(DoubleMuonC)
    DoubleMuon_HIPM.append(DoubleMuonD)
    DoubleMuon_HIPM.append(DoubleMuonE)
    DoubleMuon_HIPM.append(DoubleMuonF1)
    DoubleMuon_noHIPM.append(DoubleMuonF2)
    DoubleMuon_noHIPM.append(DoubleMuonG)
    DoubleMuon_noHIPM.append(DoubleMuonH)


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
    #for dataset in DoubleEG_list: lumi_EG += lumiEG[dataset]


    lumiMuon = {}
    lumiMuon['DoubleMuon_Run2016B'] = lumiB
    lumiMuon['DoubleMuon_Run2016C'] = lumiC
    lumiMuon['DoubleMuon_Run2016D'] = lumiD
    lumiMuon['DoubleMuon_Run2016E'] = lumiE
    lumiMuon['DoubleMuon_Run2016F'] = lumiF
    lumiMuon['DoubleMuon_Run2016G'] = lumiG
    lumiMuon['DoubleMuon_Run2016H'] = lumiH

    lumi_Muon = 0.0
    #for dataset in DoubleMuon_list: lumi_Muon += lumiMuon[dataset]
    



    filename = 'dat/Samples_cern_UltraLegacy.dat'

    ######################################
    ######## Some bin definitions ########
    ######################################
    mass_rebin =  np.concatenate((np.arange(0, 90, 5, float), np.arange(90, 130, 10), np.arange(130, 190, 15), np.arange(190, 260, 40), np.array([300])))
    mass_rebin =  np.concatenate((np.arange(0, 15, 5, float), 
                                  np.arange(15, 65, 10, float), 
                                  np.arange(65, 100, 5, float), 
                                  np.arange(100, 220, 15, float), 
                                  np.arange(220, 400, 20, float)))
#                                  np.arange(15, 65, 10, float), 
#                                  np.arange(100, 400, 15, float)))


    Ixy_reb = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 14.0, 20.0, 30.0, 40.0]
    Ixy_rebin = np.array(Ixy_reb)

    ################################
    ######## DoubleEG Plots ########
    ################################
    treeDATA_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_HIPM, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_noHIPM, 'DATA'), name = 'DATA', isdata = 1 )

    #### Prompt validation

    makePromptBKGPlot(name = 'EEprompt_mass_HIPM', lumi = 19.7, hname_SR = 'hEEpromptSR_mass', hname_CR = 'hEEpromptCR_mass', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False) 
    makePromptBKGPlot(name = 'EEprompt_mass_noHIPM', lumi = 16.15, hname_SR = 'hEEpromptSR_mass', hname_CR = 'hEEpromptCR_mass', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False) 

    #### On-Z validation

    makePromptBKGPlot(name = 'EEOnZ_trackIxy_log_HIPM', lumi = 19.7, hname_SR = 'hEEonZSR_trackIxy_log', hname_CR = 'hEEonZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'On-Z control region', xlog = True) 
    makePromptBKGPlot(name = 'EEOnZ_trackIxy_HIPM', lumi = 19.7, hname_SR = 'hEEonZSR_trackIxy', hname_CR = 'hEEonZCR_trackIxy', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, rebin = Ixy_rebin, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'EEOffZ_trackIxy_log_HIPM', lumi = 19.7, hname_SR = 'hEEoffZSR_trackIxy_log', hname_CR = 'hEEoffZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
    makePromptBKGPlot(name = 'EEOffZ_trackIxy_HIPM', lumi = 19.7, hname_SR = 'hEEoffZSR_trackIxy', hname_CR = 'hEEoffZCR_trackIxy', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False) 
    
    makePromptBKGPlot(name = 'EEOnZ_trackIxy_log_noHIPM', lumi = 16.2, hname_SR = 'hEEonZSR_trackIxy_log', hname_CR = 'hEEonZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'On-Z control region', xlog = True) 
    makePromptBKGPlot(name = 'EEOnZ_trackIxy_noHIPM', lumi = 16.2, hname_SR = 'hEEonZSR_trackIxy', hname_CR = 'hEEonZCR_trackIxy', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, rebin = Ixy_rebin, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'EEOffZ_trackIxy_log_noHIPM', lumi = 16.2, hname_SR = 'hEEoffZSR_trackIxy_log', hname_CR = 'hEEoffZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
    makePromptBKGPlot(name = 'EEOffZ_trackIxy_noHIPM', lumi = 16.2, hname_SR = 'hEEoffZSR_trackIxy', hname_CR = 'hEEoffZCR_trackIxy', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False) 

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    treeDATA_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_HIPM, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_noHIPM, 'DATA'), name = 'DATA', isdata = 1 )

    #### Prompt validation

    makePromptBKGPlot(name = 'MMPrompt_mass_HIPM', lumi = 19.7, hname_SR = 'hMMpromptSR_mass', hname_CR = 'hMMpromptCR_mass', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'Prompt control region', xlog = False) 
    makePromptBKGPlot(name = 'MMPrompt_mass_noHIPM', lumi = 16.2, hname_SR = 'hMMpromptSR_mass', hname_CR = 'hMMpromptCR_mass', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'Prompt control region', xlog = False) 

    
    #### On-Z validation

    makePromptBKGPlot(name = 'MMOnZ_trackIxy_log_HIPM', lumi = 19.7, hname_SR = 'hMMonZSR_trackIxy_log', hname_CR = 'hMMonZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input,  xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = True) 
    makePromptBKGPlot(name = 'MMOnZ_trackIxy_HIPM', lumi = 19.7, hname_SR = 'hMMonZSR_trackIxy', hname_CR = 'hMMonZCR_trackIxy', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input,  rebin = Ixy_rebin, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'MMOnZ_Lxy_HIPM', lumi = 19.7, hname_SR = 'hMMonZSR_Lxy', hname_CR = 'hMMonZCR_Lxy', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'MMOnZ_Ixy_HIPM', lumi = 19.7, hname_SR = 'hMMonZSR_Ixy', hname_CR = 'hMMonZCR_Ixy', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'MMOffZ_trackIxy_log_HIPM', lumi = 19.7, hname_SR = 'hMMoffZSR_trackIxy_log', hname_CR = 'hMMoffZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Off-Z (prompt) control region', extralabel = '', xlog = True) 
    makePromptBKGPlot(name = 'MMOffZ_trackIxy_HIPM', lumi = 19.7, hname_SR = 'hMMoffZSR_trackIxy', hname_CR = 'hMMoffZCR_trackIxy', ylog = True, treeDATA = treeDATA_HIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Off-Z (prompt) control region', extralabel = '', xlog = False) 

    makePromptBKGPlot(name = 'MMOnZ_trackIxy_log_noHIPM', lumi = 16.2, hname_SR = 'hMMonZSR_trackIxy_log', hname_CR = 'hMMonZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input,  xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = True) 
    makePromptBKGPlot(name = 'MMOnZ_trackIxy_noHIPM', lumi = 16.2, hname_SR = 'hMMonZSR_trackIxy', hname_CR = 'hMMonZCR_trackIxy', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input,  rebin = Ixy_rebin, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'MMOnZ_Lxy_noHIPM', lumi = 16.2, hname_SR = 'hMMonZSR_Lxy', hname_CR = 'hMMonZCR_Lxy', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'MMOnZ_Ixy_noHIPM', lumi = 16.2, hname_SR = 'hMMonZSR_Ixy', hname_CR = 'hMMonZCR_Ixy', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'On-Z control region', xlog = False) 
    makePromptBKGPlot(name = 'MMOffZ_trackIxy_log_noHIPM', lumi = 16.2, hname_SR = 'hMMoffZSR_trackIxy_log', hname_CR = 'hMMoffZCR_trackIxy_log', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Off-Z (prompt) control region', extralabel = '', xlog = True) 
    makePromptBKGPlot(name = 'MMOffZ_trackIxy_noHIPM', lumi = 16.2, hname_SR = 'hMMoffZSR_trackIxy', hname_CR = 'hMMoffZCR_trackIxy', ylog = True, treeDATA = treeDATA_noHIPM, inputdir = opts.input, limit = 6.0, xlabel = '', outtag = '', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Off-Z (prompt) control region', extralabel = '', xlog = False) 
