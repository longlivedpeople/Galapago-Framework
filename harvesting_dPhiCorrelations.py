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




def makeLinearFit(lumi, hname, ylog, tree, inputdir, label, name, isData, xlabel = '', outtag = '', yshift = 0.0, LLlabel = ''):

    ## Get the histogram
    if type(hname) != list:
        hist = tree.getLoopTH1F(inputdir, hname)
    else:
        hist = tree.getLoopTH1F(inputdir, hname[0])
        for n in range(1, len(hname)):
            hist.Add(tree.getLoopTH1F(inputdir, hname[n]))

    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(0.8)

    ## Create the function
    option = 'L'
    #if isData: option = 'L'
    #else: option = 'WL'

    f1 = r.TF1("f1","pol 0",0.0,3.14)
    f1.SetParameters(0, hist.GetMaximum())
    f1.SetLineColor(r.kRed)
    f2 = r.TF1("f2","pol 1",0.0,3.14)
    f2.SetParameters(0, hist.GetMaximum())
    f2.SetLineColor(r.kBlue)
    f3 = r.TF1("f3","pol 2",0.0,3.14)
    f3.SetParameters(0, hist.GetMaximum())
    f3.SetLineColor(r.kGreen+2)
    hist.Fit('f1', option, '', 0.0, 3.14)
    hist.Fit('f2', option + '+', '', 0.0, 3.14)
    hist.Fit('f3', option + '+', '', 0.0, 3.14)
    chi2ndof = f1.GetChisquare()/f1.GetNDF()
    p0_value = f1.GetParameter(0)
    p0_err = f1.GetParError(0)

    ## Adjust the maximum
    hist.SetMaximum(2*hist.GetMaximum())
    hist.SetMinimum(0.5*hist.GetMinimum())

    ## Get fit histo for the pull histogram
    fithisto1 = hist.Clone()
    fithisto1.SetMarkerColor(f1.GetLineColor())
    fithisto1.SetLineColor(f1.GetLineColor())
    fithisto2 = hist.Clone()
    fithisto2.SetMarkerColor(f2.GetLineColor())
    fithisto2.SetLineColor(f2.GetLineColor())
    fithisto3 = hist.Clone()
    fithisto3.SetMarkerColor(f3.GetLineColor())
    fithisto3.SetLineColor(f3.GetLineColor())
    for n in range(0, fithisto1.GetNbinsX()+1):
        fithisto1.SetBinContent(n, p0_value)
        fithisto1.SetBinError(n, p0_err)
        fithisto2.SetBinContent(n, f2.Eval(hist.GetBinCenter(n)))
        fithisto2.SetBinError(n, p0_err)
        fithisto3.SetBinContent(n, f3.Eval(hist.GetBinCenter(n)))
        fithisto3.SetBinError(n, p0_err)


    ## Draw the plot
    plot = Canvas.Canvas(name, 'png', 0.17, 0.7, 0.3, 0.9, 1)
    plot.addHisto(hist, 'P, SAME', label, 'p', r.kBlack, 1, 0)
    plot.addHisto(f1, '', 'Fit: P(0) = p_{0}', 'l', '', 0, 1)
    plot.addHisto(f2, '', 'Fit: P(1) = p_{0} + p_{1}x', 'l', '', 0, 2)
    plot.addHisto(f3, '', 'Fit: P(2) = p_{0} + p_{1}x + p_{2}x^{2}', 'l', '', 0, 3)
    plot.addLatex(0.26, 0.65, '#chi^{{2}}/ndof = {0:.1f}/{1}'.format(f1.GetChisquare(), f1.GetNDF()), font = 42, size = 0.035, align = 21, color = f1.GetLineColor())
    plot.addLatex(0.26, 0.6, 'p_{{0}} = {0:.1f} #pm {1:.1f}'.format(p0_value, p0_err), font = 42, size = 0.035, align = 21, color = f1.GetLineColor())
    plot.addLatex(0.51, 0.65, '#chi^{{2}}/ndof = {0:.1f}/{1}'.format(f2.GetChisquare(), f2.GetNDF()), font = 42, size = 0.035, align = 21, color = f2.GetLineColor())
    plot.addLatex(0.51, 0.6, 'p_{{0}} = {0:.1f} #pm {1:.1f}'.format(f2.GetParameter(0), f2.GetParError(0)), font = 42, size = 0.035, align = 21, color = f2.GetLineColor())
    plot.addLatex(0.51, 0.55, 'p_{{1}} = {0:.1f} #pm {1:.1f}'.format(f2.GetParameter(1), f2.GetParError(1)), font = 42, size = 0.035, align = 21, color = f2.GetLineColor())
    plot.addLatex(0.76, 0.65, '#chi^{{2}}/ndof = {0:.1f}/{1}'.format(f3.GetChisquare(), f3.GetNDF()), font = 42, size = 0.035, align = 21, color = f3.GetLineColor())
    plot.addLatex(0.76, 0.6, 'p_{{0}} = {0:.1f} #pm {1:.1f}'.format(f3.GetParameter(0), f3.GetParError(0)), font = 42, size = 0.035, align = 21, color = f3.GetLineColor())
    plot.addLatex(0.76, 0.55, 'p_{{1}} = {0:.1f} #pm {1:.1f}'.format(f3.GetParameter(1), f3.GetParError(1)), font = 42, size = 0.035, align = 21, color = f3.GetLineColor())
    plot.addLatex(0.76, 0.5, 'p_{{2}} = {0:.1f} #pm {1:.1f}'.format(f3.GetParameter(2), f3.GetParError(2)), font = 42, size = 0.035, align = 21, color = f3.GetLineColor())
    plot.saveRatio(1, 0, 0, '', hist, [fithisto1, fithisto2, fithisto3], r_ymin = 0.9, r_ymax = 1.1, outputDir = WORKPATH + 'dPhi_correlations')


    return 

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
    #treeMuon = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )
    makeLinearFit(lumi_total, 'hEESS0_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (SS0 Region)', 'DATAEG_hEESS0_dPhi_corr', 1)
    makeLinearFit(lumi_total, ['hEESSI_dPhi', 'hEEOSI_dPhi'], 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSI + SSI Regions)', 'DATAEG_hEEOSSSI_dPhi_corr', 1)
    makeLinearFit(lumi_total, ['hEESSII_dPhi', 'hEEOSII_dPhi'], 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSII + SSII Regions)', 'DATAEG_hEEOSSSII_dPhi_corr', 1)
#    makeLinearFit(lumi_total, 'hEEOSI_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSI Region)', 'DATAEG_hEEOSI_dPhi_corr', 1)
#    makeLinearFit(lumi_total, 'hEEOSII_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (OSII Region)', 'DATAEG_hEEOSII_dPhi_corr', 1)
#    makeLinearFit(lumi_total, 'hEESSI_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (SSI Region)', 'DATAEG_hEESSI_dPhi_corr', 1)
#    makeLinearFit(lumi_total, 'hEESSII_dPhi', 0, treeEG, WORKPATH + opts.input, 'DoubleEG_Run2016[B-H] (SSII Region)', 'DATAEG_hEESSII_dPhi_corr', 1)
    #makeLinearFit(lumi_total, 'hMMSS0_dPhi', 0, treeMuon, WORKPATH + opts.input, 'DoubleMuon_Run2016[B-H] (SS0 Region)', 'DATAMuon_hMM_dPhi_corr', 1)

    
    #
    # -- DY Correlation
    #

    treeDY = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50', 'DYJetsToLL_M-10to50'], 'MC'), name = 'MC', isdata = 0 )
    makeLinearFit(lumi_total, 'hEEOS0_dPhi', 0, treeDY, WORKPATH + opts.input, 'Monte Carlo: DYJetsToLL_[M-50 + M-10to50] (OS0 Region)', 'DY_hEE_dPhi_corr', 0)
    #makeLinearFit(lumi_total, 'hMMOS0_dPhi', 0, treeDY, WORKPATH + opts.input, 'Monte Carlo: DYJetsToLL_[M-50 + M-10to50] (OS0 Region)', 'DY_hMM_dPhi_corr', 0)


    #
    # -- ttbar Correlation
    #

    treeTT = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TT'], 'MC'), name = 'MC', isdata = 0 )
    makeLinearFit(lumi_total, 'hEEOS0_dPhi', 0, treeTT, WORKPATH + opts.input, 'Monte Carlo: TT (OS0 Region)', 'TT_hEE_dPhi_corr', 0)
    #makeLinearFit(lumi_total, 'hMMOS0_dPhi', 0, treeTT, WORKPATH + opts.input, 'Monte Carlo: TT (OS0 Region)', 'TT_hMM_dPhi_corr', 0)

    #
    # -- Diboson Correlation
    #

    treeVV = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW', 'WZ', 'ZZ'], 'MC'), name = 'MC', isdata = 0 )
    makeLinearFit(lumi_total, 'hEEOS0_dPhi', 0, treeVV, WORKPATH + opts.input, 'Monte Carlo: WW + WZ + ZZ (OS0 Region)', 'Diboson_hEE_dPhi_corr', 0)
    #makeLinearFit(lumi_total, 'hMMOS0_dPhi', 0, treeVV, WORKPATH + opts.input, 'Monte Carlo: WW + WZ + ZZ (OS0 Region)', 'Diboson_hMM_dPhi_corr', 0)
    
    

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


