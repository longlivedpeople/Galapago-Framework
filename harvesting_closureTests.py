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




def makeClosureTestMC(lumi, name, hBKG_A, hBKG_B, ylog, tree, inputdir, labelA = 'A', labelB = 'B', xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', rmin = 0.0, rmax = 2.0):

    print(rmin, rmax)

    ### Get histograms
    luminosity = lumi

    #hBKG_A = treeMC.getLoopTH1F(inputdir, Aname)
    #hBKG_B = treeMC.getLoopTH1F(inputdir, Bname)
    

    ### Initialize cumulative histograms
    cumA = hBKG_A.Clone('cumA')
    cumB = hBKG_B.Clone('cumB')
    cumA.Reset()
    cumB.Reset()

    #cumA = r.TH1F('cumA', '', hBKG_A.GetNbinsX(), hBKG_A.GetXaxis().GetXmin(), hBKG_A.GetXaxis().GetXmax())
    #cumB = r.TH1F('cumB', '', hBKG_B.GetNbinsX(), hBKG_B.GetXaxis().GetXmin(), hBKG_B.GetXaxis().GetXmax())
    cumA.SetTitle(';min('+xlabel+');Number of events with '+xlabel+ ' > ('+xlabel+')_{min}')

    ### Set cumulative values
    for n in range(1, hBKG_A.GetNbinsX() + 1):
        valA = 0.0
        errA = 0.0
        valB = 0.0
        errB = 0.0
        for j in range(n, hBKG_A.GetNbinsX() + 1):
            valA = valA + hBKG_A.GetBinContent(j)
            errA = errA + hBKG_A.GetBinError(j)
            valB = valB + hBKG_B.GetBinContent(j)
            errB = errB + hBKG_B.GetBinError(j)
        cumA.SetBinContent(n, valA)
        cumA.SetBinError(n, errA)
        cumB.SetBinContent(n, valB)
        cumB.SetBinError(n, errB)

    ### Histogram tuning 
    cumA.SetMarkerStyle(24)
    cumA.SetMarkerSize(0.8)
    cumA.SetMarkerColor(r.kBlue)
    cumB.SetMarkerStyle(25)
    cumB.SetMarkerSize(0.8)
    cumB.SetMarkerColor(r.kRed)
    hBKG_A.SetMarkerStyle(24)
    hBKG_A.SetMarkerSize(0.8)
    hBKG_A.SetMarkerColor(r.kBlue)
    hBKG_B.SetMarkerStyle(25)
    hBKG_B.SetMarkerSize(0.8)
    hBKG_B.SetMarkerColor(r.kRed)
    
    
    ### Get maximum
    maxValA = cumA.GetMaximum()
    maxValB = cumB.GetMaximum()
    maxVal = max(maxValA, maxValB)

    ### Set Maximum
    if not ylog:
        cumA.SetMaximum(1.3*maxVal)
        cumB.SetMaximum(1.3*maxVal)
    else:
        cumA.SetMaximum(100.0*maxVal)
        cumB.SetMaximum(100.0*maxVal)
        cumA.SetMinimum(0.1)
        cumB.SetMaximum(0.1)
 
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
        hBKG_A.SetMinimum(0.1)
        hBKG_B.SetMaximum(0.1)

    ### Canvas object
    plot = Canvas.Canvas('closuretest_'+name, 'png,pdf', 0.5, 0.75, 0.7, 0.87, 1)
    plot.addHisto(cumA, 'P', labelA, 'p', r.kBlue, 1, 0)
    plot.addHisto(cumB, 'P, SAME', labelB, 'p', r.kRed, 1, 1)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)

    plot.addLatex(0.85, 0.64, 'Tail cumulative of x = #int^{#infty}_{x_{min}} N(x) dx', font = 42, align = 31)
    plot.addLatex(0.17, 0.86, 'Closure test', font = 62, align = 11, size = 0.044)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/ClosureTests_' + outtag + '/'
    plot.saveRatio(1, 0, ylog, luminosity, cumA, cumB, r_ymin = rmin, r_ymax = rmax, label="r_{SR/CR}", outputDir = outdir)

    ### Main comparison
    cplot = Canvas.Canvas('plaincomparison_'+name, 'png,pdf', 0.5, 0.75, 0.7, 0.87, 1)
    cplot.addHisto(hBKG_A, 'P', labelA, 'p', r.kRed, 1, 0)
    cplot.addHisto(hBKG_B, 'P, SAME', labelB, 'p', r.kBlue, 1, 1)
    if LLlabel == 'EE':
        cplot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        cplot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)
    cplot.addLatex(0.17, 0.86, 'Closure test', font = 62, align = 11, size = 0.044)
    cplot.saveRatio(1, 0, ylog, luminosity, hBKG_A, hBKG_B, label="SR/CR", outputDir = outdir)

    

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
    #DoubleEG_list.append(DoubleEGH)

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
    #DoubleMuon_list.append(DoubleMuonG)
    #DoubleMuon_list.append(DoubleMuonH)


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

    lumiTotal = lumiB + lumiC + lumiD + lumiE + lumiF + lumiG + lumiH


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

    ################################
    ######## DoubleEG Plots ########
    ################################
    
    #
    # -- DY Closure
    #
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50', 'DYJetsToLL_M-10to50'], 'MC'), name = 'MC', isdata = 0 )
    hBKG_A = treeMC.getLoopTH1F(opts.input, 'hEESROS_trackIxy')
    hBKG_B = treeMC.getLoopTH1F(opts.input, 'hEECROS_trackIxy')

    makeClosureTestMC(lumi = lumiTotal, name = 'EEDY', hBKG_A = hBKG_A, hBKG_B = hBKG_B, ylog = True, tree = treeMC, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{ee} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{ee} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannelV2', yshift = 0.0, LLlabel = 'EE', extralabel = '')

    hBKG_A = treeMC.getLoopTH1F(opts.input, 'hMMSROS_trackIxy')
    hBKG_B = treeMC.getLoopTH1F(opts.input, 'hMMCROS_trackIxy')

    makeClosureTestMC(lumi = lumiTotal, name = 'MMDY', hBKG_A = hBKG_A, hBKG_B = hBKG_B, ylog = True, tree = treeMC, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{#mu#mu} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{#mu#mu} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannelV2', yshift = 0.0, LLlabel = 'MM', extralabel = '')

    #
    # -- ttbar Closure
    #
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TT'], 'MC'), name = 'MC', isdata = 0 )
    hBKG_A = treeMC.getLoopTH1F(opts.input, 'hEESROS_trackIxy')
    hBKG_B = treeMC.getLoopTH1F(opts.input, 'hEECROS_trackIxy')
    newbin = np.array((0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 14.0, 20.0))
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)


    makeClosureTestMC(lumi = lumiTotal, name = 'TT', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeMC, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{ee} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{ee} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannelV2', yshift = 0.0, LLlabel = 'EE', extralabel = '', rmin = 0.5, rmax = 2.1)

    hBKG_A = treeMC.getLoopTH1F(opts.input, 'hMMSROS_trackIxy')
    hBKG_B = treeMC.getLoopTH1F(opts.input, 'hMMCROS_trackIxy')
    newbin = np.array((0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0))
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)


    makeClosureTestMC(lumi = lumiTotal, name = 'TT', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeMC, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{#mu#mu} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{#mu#mu} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannelV2', yshift = 0.0, LLlabel = 'MM', extralabel = '')

    #
    # -- Diboson Closure
    #
    treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW', 'WZ', 'ZZ'], 'MC'), name = 'MC', isdata = 0 )
    hBKG_A = treeMC.getLoopTH1F(opts.input, 'hEESROS_trackIxy')
    hBKG_B = treeMC.getLoopTH1F(opts.input, 'hEECROS_trackIxy')
    newbin = np.array([0.0, 1.0, 2.0, 4.0, 20.0])
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, hBKG_A.GetName()+'_rebined', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, hBKG_B.GetName()+'_rebined', newbin)

    makeClosureTestMC(lumi = lumiTotal, name = 'Diboson', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeMC, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{ee} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{ee} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannelV2', yshift = 0.0, LLlabel = 'EE', extralabel = '')
    
    hBKG_A = treeMC.getLoopTH1F(opts.input, 'hMMSROS_trackIxy')
    hBKG_B = treeMC.getLoopTH1F(opts.input, 'hMMCROS_trackIxy')
    newbin = np.array([0.0, 1.0, 2.0, 4.0, 20.0])
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, hBKG_A.GetName()+'_rebined', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, hBKG_B.GetName()+'_rebined', newbin)

    makeClosureTestMC(lumi = lumiTotal, name = 'Diboson', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeMC, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{#mu#mu} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{#mu#mu} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannelV2', yshift = 0.0, LLlabel = 'MM', extralabel = '')


    #
    # -- QCD Correlation
    #
    treeEG = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )
    treeMuon = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )

    hBKG_A = treeEG.getLoopTH1F(opts.input, 'hEESRII_trackIxy')
    hBKG_B = treeEG.getLoopTH1F(opts.input, 'hEECRII_trackIxy')
    newbin = np.array([0.0, 1.0, 2.0, 4.0, 7.0, 10.0, 14.0, 20.0])
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, hBKG_A.GetName()+'_rebined', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, hBKG_B.GetName()+'_rebined', newbin)

    makeClosureTestMC(lumi = lumi_EG, name = 'QCDEst', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeEG, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{ee} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{ee} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannelV2', yshift = 0.0, LLlabel = 'EE', extralabel = '')

    hBKG_A = treeMuon.getLoopTH1F(opts.input, 'hMMSRII_trackIxy')
    hBKG_B = treeMuon.getLoopTH1F(opts.input, 'hMMCRII_trackIxy')
    newbin = np.array([0.0, 1.0, 2.0, 4.0, 7.0, 10.0, 14.0, 20.0])
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, hBKG_A.GetName()+'_rebined', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, hBKG_B.GetName()+'_rebined', newbin)

    makeClosureTestMC(lumi = lumi_EG, name = 'QCDEst', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeMuon, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{#mu#mu} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{#mu#mu} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannelV2', yshift = 0.0, LLlabel = 'MM', extralabel = '')

    #
    # -- Wjets Correlation
    #
    hBKG_A = treeEG.getLoopTH1F(opts.input, 'hEESRI_trackIxy')
    hBKG_B = treeEG.getLoopTH1F(opts.input, 'hEECRI_trackIxy')
    newbin = np.array([0.0, 1.0, 2.0, 3.0, 6.0, 10.0, 20.0])
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, hBKG_A.GetName()+'_rebined', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, hBKG_B.GetName()+'_rebined', newbin)

    makeClosureTestMC(lumi = lumi_EG, name = 'WjetsEst', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeEG, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{ee} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{ee} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannelV2', yshift = 0.0, LLlabel = 'EE', extralabel = '', rmin = 0.5, rmax = 1.5)

    hBKG_A = treeMuon.getLoopTH1F(opts.input, 'hMMSRI_trackIxy')
    hBKG_B = treeMuon.getLoopTH1F(opts.input, 'hMMCRI_trackIxy')
    newbin = np.array([0.0, 1.0, 2.0, 3.0, 6.0, 10.0, 20.0])
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, hBKG_A.GetName()+'_rebined', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, hBKG_B.GetName()+'_rebined', newbin)

    makeClosureTestMC(lumi = lumi_Muon, name = 'WjetsEst', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeMuon, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi|_{#mu#mu} > #pi/2', labelA = 'Signal Region: |#Delta#Phi|_{#mu#mu} < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannelV2', yshift = 0.0, LLlabel = 'MM', extralabel = '', rmin = 0.5, rmax = 1.5)




