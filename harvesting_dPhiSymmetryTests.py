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
#from include.Utils import *



def makeSymmetryTest(lumi, hname, ylog, tree, inputdir, label, name, isData, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', rebin = 0, maxY = False):

    ## Get the histogram
    if type(hname) != list:
        hist = tree.getLoopTH1F(inputdir, hname)
    else:
        hist = tree.getLoopTH1F(inputdir, hname[0])
        for n in range(1, len(hname)):
            hist.Add(tree.getLoopTH1F(inputdir, hname[n]))

    nbins = int(hist.GetNbinsX()/2.0)
    xmin = 0.0 # default
    xmax = 3.14/2.0 # default

    ## Define forward and inverted backward histograms
    forward = r.TH1F('Forward', ';;Event yield', nbins, xmin, xmax)
    invBackward = r.TH1F('InvBackward', ';;Event yield', nbins, xmin, xmax)

    for n in range(1, nbins + 1):
        forward.SetBinContent(n, hist.GetBinContent(n))
        invBackward.SetBinContent(n, hist.GetBinContent(2*nbins + 1 - n))
        forward.SetBinError(n, hist.GetBinError(n))
        invBackward.SetBinError(n, hist.GetBinError(2*nbins + 1 - n))


    ## Histogram tunning
    if rebin:
        forward.Rebin(rebin)
        invBackward.Rebin(rebin)

    forward.SetMaximum(1.45*forward.GetMaximum())
    forward.SetMinimum(0.0)
    forward.SetMarkerStyle(24)
    invBackward.SetMarkerStyle(24)

    plot = Canvas.Canvas(name, 'png,pdf', 0.15, 0.72, 0.45, 0.82, 1)
    plot.addHisto(forward, 'P', 'Forward: |#Delta#Phi|, for |#Delta#Phi| #in [0, #pi/2]', 'p', r.kGreen+2, 1, 0)
    plot.addHisto(invBackward, 'P,SAME', 'Backward inverted: |#Delta#Phi| - #pi/2, for |#Delta#Phi| #in [#pi/2, #pi]', 'p', r.kRed+2, 1, 1)
    plot.addLatex(0.18, 0.85, label, font = 62, size = 0.035, align = 11)
    plot.saveRatio(1, 0, 0, '', forward, invBackward, r_ymin = 0.8, r_ymax = 1.2, label = 'Forward / Inv. Backward',  xlog = False, outputDir = WORKPATH + 'SymmetryResults/', maxYnumbers = maxY)

    print('>>>>>>>>>> KOLMOGOROV test for ' + label)
    print('>>>>>>>>>> ' + str(forward.KolmogorovTest(invBackward)))
    print('>>>>>>>>>> Chi2 test for ' + label)
    print('>>>>>>>>>> ' + str(forward.Chi2Test(invBackward, "WWP")))

    return 


def makeSymmetryTest2(lumi, hname, ylog, tree, inputdir, label, name, isData, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', color = r.kBlack):

    ## Get the histogram
    if type(hname) != list:
        hist = tree.getLoopTH1F(inputdir, hname)
    else:
        hist = tree.getLoopTH1F(inputdir, hname[0])
        for n in range(1, len(hname)):
            hist.Add(tree.getLoopTH1F(inputdir, hname[n]))

    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(0.8)

    ## Adjust the maximum
    hist.SetMaximum(2*hist.GetMaximum())
    hist.SetMinimum(0.5*hist.GetMinimum())

    ## Get Mean and Median:
    mean = hist.GetMean()
    print('Mean: ', hist.GetMean())
    y = 0.
    q = np.array([0.])
    prob = np.array([0.5])
    y = hist.GetQuantiles(1, q, prob)
    print('Median: ', q[0])
    print('Reference: ', 3.14/2.0)

    ## Get populations
    nbins = hist.GetNbinsX()
    simbin = int(nbins/2)
    forwardErr = np.array([0.])
    backwardErr = np.array([0.])
    Nforward = hist.IntegralAndError(1, simbin, forwardErr)
    Nbackward = hist.IntegralAndError(simbin + 1, nbins, backwardErr)
    print("N in forward: ", Nforward, forwardErr[0])
    print("N in backward: ", Nbackward, backwardErr[0])

    fullStat = [Nforward, forwardErr[0], Nbackward, backwardErr[0]]

    ## Plot
    plot = Canvas.Canvas(name, 'png,pdf', 0.14, 0.8, 0.3, 0.9, 1)
    plot.addHisto(hist, 'P, SAME', label, 'p', color, 1, 0)
    plot.addLine(mean, hist.GetMinimum(), mean, hist.GetMaximum(), r.kBlue)
    plot.addLatex(0.18, 0.75, 'Mean: #bar{{|#Delta#Phi|}} = {0:.4f}'.format(mean), font = 42, size = 0.03, align = 11, color = r.kBlue)
    plot.addLine(q[0], hist.GetMinimum(), q[0], hist.GetMaximum(), r.kRed)
    plot.addLatex(0.18, 0.71, 'Median: |#Delta#Phi|_{{1/2}} = {0:.4f}'.format(q[0]), font = 42, size = 0.03, align = 11, color = r.kRed)
    plot.addLatex(0.18, 0.66, 'N_{{forward}} = {0} #pm {1}'.format(int(Nforward), int(forwardErr[0])), font = 42, size = 0.03, align = 11, color = r.kBlack)
    plot.addLatex(0.18, 0.62, 'N_{{backward}} = {0} #pm {1}'.format(int(Nbackward), int(backwardErr[0])), font = 42, size = 0.03, align = 11, color = r.kBlack)
    plot.save(1, 0, ylog, '', '', outputDir = WORKPATH + 'SymmetryResults/')

    return fullStat


def computeRatio(NSR, eSR, NCR, eCR):

    # value:
    r = NSR/NCR

    # error:
    error = math.sqrt((1/NCR**4)*(NCR*NCR*eSR*eSR + NSR*NSR*eCR*eCR))

    return r, error


def make2DPlot(h2, name, outdir):

    plot = Canvas.Canvas(name, 'png', 0.17, 0.8, 0.3, 0.9, 1)
    plot.addHisto(h2, 'COLZ', '', '', '', 1, 0)
    plot.save(0, 0, 0, '', '', outputDir = WORKPATH + outdir)


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
    r.gStyle.SetOptFit(0)

    ############# Muon data definition
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
    Signals.append('HXX_1000_350_350mm')
    Signals.append('HXX_1000_350_35mm')
    Signals.append('HXX_1000_150_100mm')
    Signals.append('HXX_1000_150_10mm')
    Signals.append('HXX_400_150_400mm')
    Signals.append('HXX_400_50_400mm')
    Signals.append('HXX_400_50_40mm')
    Signals.append('HXX_400_50_4mm')

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
    

    filename = 'dat/Samples_cern_fillingv2.dat'

    letter = ['B', 'C', 'D', 'E', 'F', 'G', 'H']

    
    #
    # -- QCD Correlation
    #
   
    treeEG = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )
    treeMuon = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )
    makeSymmetryTest(lumi_total, 'hEESS0_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (SS0 Region)', 'DATAEG_hEESS0_dPhi_corr', 1)
    makeSymmetryTest2(lumi_total, 'hEESS0_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (SS0 Region)', 'DATAEG_hEESS0_dPhi_test', 1)
    makeSymmetryTest2(lumi_total, 'hEEOS0disp_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (SS0 Region)', 'DATAEG_hEEOS0disp_dPhi_test', 1, color = r.kBlue)
    makeSymmetryTest(lumi_total, ['hEESSII_dPhi', 'hEEOSII_dPhi'], 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSII + SSII Regions)', 'DATAEG_hEEOSSSII_dPhi_corr', 1)
    EE_QCDstat = makeSymmetryTest2(lumi_total, ['hEESSII_dPhi', 'hEEOSII_dPhi'], 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSII + SSII Regions)', 'DATAEG_hEEOSSSII_dPhi_test', 1)
    makeSymmetryTest(lumi_total, ['hMMSSII_dPhi', 'hMMOSII_dPhi'], 0, treeMuon, WORKPATH + opts.input, 'DoubleMuon_Run2016[B-H] (OSII + SSII Regions)', 'DATAMuon_hMMOSSSII_dPhi_corr', 1)
    MM_QCDstat = makeSymmetryTest2(lumi_total, ['hMMSSII_dPhi', 'hMMOSII_dPhi'], 0, treeMuon, WORKPATH + opts.input, 'DoubleMuon_Run2016[B-H] (OSII + SSII Regions)', 'DATAMuon_hMMOSSSII_dPhi_test', 1)



    #
    # -- W+jets Correlation
    #
    makeSymmetryTest(lumi_total, ['hEESSI_dPhi', 'hEEOSI_dPhi'], 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSI + SSI Regions)', 'DATAEG_hEEOSSSI_dPhi_corr', 1)
    EE_WJetsstat = makeSymmetryTest2(lumi_total, ['hEESSI_dPhi', 'hEEOSI_dPhi'], 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSI + SSI Regions)', 'DATAEG_hEEOSSSI_dPhi_test', 1)
    makeSymmetryTest(lumi_total, ['hMMSSI_dPhi', 'hMMOSI_dPhi'], 0, treeMuon, WORKPATH + opts.input, 'DoubleMuon_Run2016[B-H] (OSI + SSI Regions)', 'DATAMuon_hMMOSSSI_dPhi_corr', 1)
    MM_WJetsstat = makeSymmetryTest2(lumi_total, ['hMMSSI_dPhi', 'hMMOSI_dPhi'], 0, treeMuon, WORKPATH + opts.input, 'DoubleMuon_Run2016[B-H] (OSI + SSI Regions)', 'DATAMuon_hMMOSSSI_dPhi_test', 1)
    
    
    #
    # -- DY Correlation
    #

    treeDY = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50', 'DYJetsToLL_M-10to50'], 'MC'), name = 'MC', isdata = 0 )
    makeSymmetryTest(lumi_total, 'hEEOS0_dPhi', 0, treeDY, WORKPATH + opts.input, 'Monte Carlo: DYJetsToLL_[M-50 + M-10to50] (OS0 Region)', 'DY_hEE_dPhi_corr', 0)
    EE_DYstat = makeSymmetryTest2(lumi_total, 'hEEOS0_dPhi', 0, treeDY, WORKPATH + opts.input, 'Monte Carlo: DYJetsToLL_[M-50 + M-10to50] (OS0 Region)', 'DY_hEE_dPhi_test', 0)
    makeSymmetryTest(lumi_total, 'hMMOS0_dPhi', 0, treeDY, WORKPATH + opts.input, 'Monte Carlo: DYJetsToLL_[M-50 + M-10to50] (OS0 Region)', 'DY_hMM_dPhi_corr', 0)
    MM_DYstat = makeSymmetryTest2(lumi_total, 'hMMOS0_dPhi', 0, treeDY, WORKPATH + opts.input, 'Monte Carlo: DYJetsToLL_[M-50 + M-10to50] (OS0 Region)', 'DY_hMM_dPhi_test', 0)
    hmm = treeDY.getLoopTH2F(WORKPATH + opts.input, 'hMMOS0_DeltaPhi_dPhi') 
    hee = treeDY.getLoopTH2F(WORKPATH + opts.input, 'hEEOS0_DeltaPhi_dPhi') 
    make2DPlot(hmm, 'MM_DY_DeltaPhi_dPhi', '2DPlots')
    make2DPlot(hee, 'EE_DY_DeltaPhi_dPhi', '2DPlots')

   
    #
    # -- ttbar Correlation
    #
    
    treeTT = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TT'], 'MC'), name = 'MC', isdata = 0 )
    makeSymmetryTest(lumi_total, 'hEEOS0_dPhi', 0, treeTT, WORKPATH + opts.input, 'Monte Carlo: TT (OS0 Region)', 'TT_hEE_dPhi_corr', 0)
    EE_TTstat = makeSymmetryTest2(lumi_total, 'hEEOS0_dPhi', 0, treeTT, WORKPATH + opts.input, 'Monte Carlo: TT (OS0 Region)', 'TT_hEE_dPhi_test', 0)
    makeSymmetryTest(lumi_total, 'hMMOS0_dPhi', 0, treeTT, WORKPATH + opts.input, 'Monte Carlo: TT (OS0 Region)', 'TT_hMM_dPhi_corr', 0)
    MM_TTstat = makeSymmetryTest2(lumi_total, 'hMMOS0_dPhi', 0, treeTT, WORKPATH + opts.input, 'Monte Carlo: TT (OS0 Region)', 'TT_hMM_dPhi_test', 0)
    hmm = treeTT.getLoopTH2F(WORKPATH + opts.input, 'hMMOS0_DeltaPhi_dPhi') 
    hee = treeTT.getLoopTH2F(WORKPATH + opts.input, 'hEEOS0_DeltaPhi_dPhi') 
    make2DPlot(hmm, 'MM_TT_DeltaPhi_dPhi', '2DPlots')
    make2DPlot(hee, 'EE_TT_DeltaPhi_dPhi', '2DPlots')

    #
    # -- Diboson Correlation
    #

    treeVV = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW', 'WZ', 'ZZ'], 'MC'), name = 'MC', isdata = 0 )
    makeSymmetryTest(lumi_total, 'hEEOS0_dPhi', 0, treeVV, WORKPATH + opts.input, 'Monte Carlo: WW + WZ + ZZ (OS0 Region)', 'Diboson_hEE_dPhi_corr', 0, rebin = 3)
    EE_VVstat = makeSymmetryTest2(lumi_total, 'hEEOS0_dPhi', 0, treeVV, WORKPATH + opts.input, 'Monte Carlo: WW + WZ + ZZ (OS0 Region)', 'Diboson_hEE_dPhi_test', 0)
    makeSymmetryTest(lumi_total, 'hMMOS0_dPhi', 0, treeVV, WORKPATH + opts.input, 'Monte Carlo: WW + WZ + ZZ (OS0 Region)', 'Diboson_hMM_dPhi_corr', 0, rebin = 3)
    MM_VVstat = makeSymmetryTest2(lumi_total, 'hMMOS0_dPhi', 0, treeVV, WORKPATH + opts.input, 'Monte Carlo: WW + WZ + ZZ (OS0 Region)', 'Diboson_hMM_dPhi_test', 0)
    hmm = treeVV.getLoopTH2F(WORKPATH + opts.input, 'hMMOS0_DeltaPhi_dPhi') 
    hee = treeVV.getLoopTH2F(WORKPATH + opts.input, 'hEEOS0_DeltaPhi_dPhi') 
    make2DPlot(hmm, 'MM_VV_DeltaPhi_dPhi', '2DPlots')
    make2DPlot(hee, 'EE_VV_DeltaPhi_dPhi', '2DPlots')    

    #
    # -- Plot with N(SR)/N(CR) for different background sources
    #    
    """
    leftMargin = 0.13
    rightMargin = 0.10
    topMargin = 0.08
    bottomMargin = 0.16
    r.gStyle.SetPadBottomMargin(bottomMargin)
    r.gStyle.SetPadLeftMargin(leftMargin)
    xWidth = 1.0 - leftMargin - rightMargin
    binWidth = xWidth/5


    ###### -> Electron and muon channel summary
    EE_summary = r.TH1F('EE_Summary', '', 5, 0, 5) 
    EE_summary.GetXaxis().SetLabelSize(0)
    EE_summary.GetYaxis().SetTitle("r_{SR/CR} factor")
    EE_summary.GetYaxis().SetTitleOffset(1.32)
    EE_summary.GetXaxis().SetNdivisions(5)
    EE_summary.SetMarkerColor(r.kBlack)
    EE_summary.SetMarkerStyle(24)
    EE_summary.SetMarkerSize(0.8)
    

    MM_summary = r.TH1F('MM_Summary', '', 5, 0, 5) 
    MM_summary.GetXaxis().SetLabelSize(0)
    MM_summary.GetXaxis().SetNdivisions(5)
    MM_summary.SetMarkerStyle(24)
    MM_summary.SetMarkerColor(r.kRed)
    MM_summary.SetLineColor(r.kRed)
    MM_summary.SetMarkerSize(0.8)

    # Drell-Yan:
    EE_DYr, EE_DYe = computeRatio(EE_DYstat[0], EE_DYstat[1], EE_DYstat[2], EE_DYstat[3])
    EE_summary.SetBinContent(1, EE_DYr)
    EE_summary.SetBinError(1, EE_DYe)
    MM_DYr, MM_DYe = computeRatio(MM_DYstat[0], MM_DYstat[1], MM_DYstat[2], MM_DYstat[3])
    MM_summary.SetBinContent(1, MM_DYr)
    MM_summary.SetBinError(1, MM_DYe)
    # TT:
    EE_TTr, EE_TTe = computeRatio(EE_TTstat[0], EE_TTstat[1], EE_TTstat[2], EE_TTstat[3])
    EE_summary.SetBinContent(2, EE_TTr)
    EE_summary.SetBinError(2, EE_TTe)
    MM_TTr, MM_TTe = computeRatio(MM_TTstat[0], MM_TTstat[1], MM_TTstat[2], MM_TTstat[3])
    MM_summary.SetBinContent(2, MM_TTr)
    MM_summary.SetBinError(2, MM_TTe)
    # VV:
    EE_VVr, EE_VVe = computeRatio(EE_VVstat[0], EE_VVstat[1], EE_VVstat[2], EE_VVstat[3])
    EE_summary.SetBinContent(3, EE_VVr)
    EE_summary.SetBinError(3, EE_VVe)
    MM_VVr, MM_VVe = computeRatio(MM_VVstat[0], MM_VVstat[1], MM_VVstat[2], MM_VVstat[3])
    MM_summary.SetBinContent(3, MM_VVr)
    MM_summary.SetBinError(3, MM_VVe)
    # QCD:
    EE_QCDr, EE_QCDe = computeRatio(EE_QCDstat[0], EE_QCDstat[1], EE_QCDstat[2], EE_QCDstat[3])
    EE_summary.SetBinContent(4, EE_QCDr)
    EE_summary.SetBinError(4, EE_QCDe)
    MM_QCDr, MM_QCDe = computeRatio(MM_QCDstat[0], MM_QCDstat[1], MM_QCDstat[2], MM_QCDstat[3])
    MM_summary.SetBinContent(4, MM_QCDr)
    MM_summary.SetBinError(4, MM_QCDe)
    # WJets 
    EE_WJetsr, EE_WJetse = computeRatio(EE_WJetsstat[0], EE_WJetsstat[1], EE_WJetsstat[2], EE_WJetsstat[3])
    EE_summary.SetBinContent(5, EE_WJetsr)
    EE_summary.SetBinError(5, EE_WJetse)
    MM_WJetsr, MM_WJetse = computeRatio(MM_WJetsstat[0], MM_WJetsstat[1], MM_WJetsstat[2], MM_WJetsstat[3])
    MM_summary.SetBinContent(5, MM_WJetsr)
    MM_summary.SetBinError(5, MM_WJetse)
    

    # Set Limits
    EE_summary.SetMaximum(1.2)
    EE_summary.SetMinimum(0.8)
    MM_summary.SetMaximum(1.2)
    MM_summary.SetMinimum(0.8)

    sumplot = Canvas.Canvas('Ratio', 'png,pdf', 0.17, 0.81, 0.5, 0.89, 1)
    sumplot.addHisto(EE_summary, 'P', 'Baseline dielectron', 'p', '', 1, 0)
    sumplot.addHisto(MM_summary, 'P,SAME', 'Baseline dimuons', 'p', '', 1, 1)
    sumplot.addLatex(leftMargin + (0+0.5)*binWidth, 0.12, 'Drell-Yan', size = 0.033, align = 22)
    sumplot.addLatex(leftMargin + (1+0.5)*binWidth, 0.12, 'TT', size = 0.033, align = 22)
    sumplot.addLatex(leftMargin + (2+0.5)*binWidth, 0.12, 'Diboson', size = 0.033, align = 22)
    sumplot.addLatex(leftMargin + (3+0.5)*binWidth, 0.12, 'QCD', size = 0.033, align = 22)
    sumplot.addLatex(leftMargin + (4+0.5)*binWidth, 0.12, 'W+jets', size = 0.033, align = 22)
    sumplot.addLatex(leftMargin + (1 - leftMargin - rightMargin)/2.0, 0.07, 'Background sources', size = 0.04, align = 22, font = 62)
    sumplot.addLine(0, 1.0, 5, 1.0, r.kRed)
    sumplot.save(1, 0, 0, '', '', outputDir = 'SymmetryResults/')

    """


    #
    # -- cos(alpha) // dPhi correlation in signal models
    #

    """
    for signal in Signals:
        treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, [signal], 'SI'), name = 'SI', isdata = 0 )
        labels = signal[4:-2].split('_')
        mH = int(labels[0])
        mX = int(labels[1])
        ctau = int(labels[2])

        histo = treeSI.getLoopTH2F(WORKPATH + opts.input, 'hMMOS0_cosAlpha_dPhi')
        plot = Canvas.Canvas(signal + 'MM_dPhi_Alpha', 'png', 0.17, 0.8, 0.3, 0.9, 1)
        plot.addHisto(histo, 'COLZ', '', '', '', 1, 0)
        plot.addLine(-1, 3.14/2, 1, 3.14/2, r.kBlack)
        plot.addLatex(0.95, 0.935, 'H#rightarrowXX: m_{{H}} = {0}, m_{{X}} = {1}, c#tau = {2} mm'.format(mH, mX, ctau), font = 42, size = 0.03, align = 31)
        plot.addLatex(0.7, 0.55, 'Control Region', font = 62, size = 0.03, align = 22)
        plot.addLatex(0.7, 0.5, 'Signal Region', font = 62, size = 0.03, align = 22)
        plot.save(0, 0, 0, '', '', outputDir = WORKPATH + 'dPhi_correlations')

        
        
        histo = treeSI.getLoopTH2F(WORKPATH + opts.input, 'hEEOS0_cosAlpha_dPhi')
        plot = Canvas.Canvas(signal + 'EE_dPhi_Alpha', 'png', 0.17, 0.8, 0.3, 0.9, 1)
        plot.addHisto(histo, 'COLZ', '', '', '', 1, 0)
        plot.addLine(-1, 3.14/2, 1, 3.14/2, r.kBlack)
        plot.addLatex(0.95, 0.935, 'H#rightarrowXX: m_{{H}} = {0}, m_{{X}} = {1}, c#tau = {2} mm'.format(mH, mX, ctau), font = 42, size = 0.03, align = 31)
        plot.addLatex(0.7, 0.55, 'Control Region', font = 62, size = 0.03, align = 22)
        plot.addLatex(0.7, 0.5, 'Signal Region', font = 62, size = 0.03, align = 22)
        plot.save(0, 0, 0, '', '', outputDir = WORKPATH + 'dPhi_correlations')
        """


