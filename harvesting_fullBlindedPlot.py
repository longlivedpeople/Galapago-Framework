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




def makeBlindedPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, rebin = False, lines = [], xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = ''):


    ### Get histograms
    luminosity = lumi

    print(inputdir, hname_bkg)
    hbkg_ = treeDATA.getLoopTH1F(inputdir, hname_bkg)
    
    ### rebinins:
    if type(rebin) != bool:
        if len(rebin) > 1:
            hbkg = hbkg_.Rebin(len(rebin)-1, hbkg_.GetName() + '_rebined', rebin)
        else:
            hbkg = hbkg_.Rebin(rebin)
    else:
        hbkg = hbkg_.Clone()

    hbkg.SetFillColorAlpha(r.kCyan-6, 0.8) 
    hbkg.SetLineColor(r.kCyan-6) 
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)

    for _i, _h in enumerate(hSIS.GetHists()):
        _h.Scale(lumi/35.87)
        s_histos.append(copy.deepcopy(_h))


    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])

    ### Set Maximum
    if not ylog:
        hbkg.SetMaximum(1.3*maxVal)
    else:
        hbkg.SetMaximum(100.0*maxVal)
        hbkg.SetMinimum(0.0001)

    if ymax: hbkg.SetMaximum(ymax) 

    ### Count background events
    backtotal = 0.0
    for n in range(1, hbkg.GetNbinsX() + 1):
        if hbkg.GetBinLowEdge(n) > 20.0:
            backtotal += hbkg.GetBinContent(n)

    print("background total", backtotal)

    ### Canvas object
    plot = Canvas.Canvas('Blinded_'+hname_bkg, 'png', 0.15, 0.65, 0.6, 0.89, 1)
    plot.addHisto(hbkg, 'HIST', 'Background (Data-driven)', 'f', '', 1, 0)
    
    ### Add signals:
    for i,_h in enumerate(s_histos):
        stotal_6 = 0.0
        stotal_2 = 0.0
        for n in range(1, _h.GetNbinsX() + 1):
            if _h.GetBinLowEdge(n) > 6.0:
                stotal_6 += _h.GetBinContent(n)
            if _h.GetBinLowEdge(n) > 2.0:
                stotal_2 += _h.GetBinContent(n)

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        legend = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, hbkg.GetMaximum(), r.kBlack)

    ### Extralabel
    #plot.addLatex(0.17, 0.8, extralabel, font = 62)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
    plot.save(1, 0, ylog, luminosity, '', outputDir = outdir, xlog = False)

    

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
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C')
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

    """
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
    DoubleEG_list = DoubleEG_HIPM + DoubleEG_noHIPM
    """
    DoubleEG_list = []
    DoubleEG_list.append('DoubleEG_Run2016B_HIPM')
    DoubleEG_list.append('DoubleEG_Run2016C_HIPM')
    DoubleEG_list.append('DoubleEG_Run2016D_HIPM')
    DoubleEG_list.append('DoubleEG_Run2016E_HIPM')
    DoubleEG_list.append('DoubleEG_Run2016F_HIPM')
    DoubleEG_list.append('DoubleEG_Run2016F_noHIPM')
    DoubleEG_list.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG_list.append('DoubleEG_Run2016H_noHIPM')

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
    DoubleMuon_list = DoubleMuon_HIPM + DoubleMuon_noHIPM


    ############# Signal definition
    Signals = []
    Signals.append('HSS_400_50_1_2016')
    Signals.append('HSS_400_50_10_2016')
    Signals.append('HSS_400_50_100_2016')
    Signals.append('HSS_400_50_1000_2016')
    Signals.append('HSS_400_50_10000_2016')


    ############# Luminosity definition
    lumi_Muon = 35.9
    lumi_EG = 35.9



    filename = 'dat/Samples_cern_UltraLegacy.dat'


    ### Tree SI Tree
    treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'signals_2016.dat', Signals, 'SI'), name = 'SI', isdata = 0 )

    ################################
    ######## DoubleEG Plots ########
    ################################
       
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )

    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRI_Ixy', hname_bkg = 'hEEBCRI_Ixy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRII_Ixy', hname_bkg = 'hEEBCRII_Ixy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRI_trackIxy', hname_bkg = 'hEEBCRI_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRII_trackIxy', hname_bkg = 'hEEBCRII_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRI_Lxy', hname_bkg = 'hEEBCRI_Lxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRII_Lxy', hname_bkg = 'hEEBCRII_Lxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRI_mass', hname_bkg = 'hEEBCRI_mass', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEESRII_mass', hname_bkg = 'hEEBCRII_mass', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )

   # makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRI_Ixy', hname_bkg = 'hMMBCRI_Ixy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
   # makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRII_Ixy', hname_bkg = 'hMMBCRII_Ixy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 

    makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRI_trackIxy', hname_bkg = 'hMMBCRI_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRII_trackIxy', hname_bkg = 'hMMBCRII_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 

   # makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRI_Lxy', hname_bkg = 'hMMBCRI_Lxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
   # makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRII_Lxy', hname_bkg = 'hMMBCRII_Lxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 

   # makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRI_mass', hname_bkg = 'hMMBCRI_mass', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
   # makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMSRII_mass', hname_bkg = 'hMMBCRII_mass', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI, xlabel = '', outtag = 'SR_plots', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
