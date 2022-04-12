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



def makeAgreementTest(lumi, hname1, hname2, ylog, tree, inputdir, label1, label2, labela, labelb, name, isData, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', rebin = 0, rmin = 0.9, rmax = 1.1, maxY = False):

    ## Get the histogram
    histo1 = tree.getLoopTH1F(inputdir, hname1, doOF = False)
    histo2 = tree.getLoopTH1F(inputdir, hname2, doOF = False)


    ## Histogram tunning
    if rebin:
        histo1.Rebin(rebin)
        histo2.Rebin(rebin)

    histo1.SetMaximum(1.6*histo1.GetMaximum())
    histo1.SetMinimum(0.0)
    histo1.SetMarkerStyle(20)
    histo2.SetMarkerStyle(20)

    plot = Canvas.Canvas(name, 'png,pdf', 0.15, 0.6, 0.45, 0.75, 1, lsize = 0.045)
    plot.addHisto(histo1, 'P', label1, 'p', r.kBlack, 1, 0)
    plot.addHisto(histo2, 'P,SAME', label2, 'p', r.kRed, 1, 1)
    plot.addLatex(0.18, 0.78, labela, font = 42, size = 0.045, align = 11)
    plot.addLatex(0.9, 0.88, labelb, font = 42, size = 0.045, align = 31)
    plot.saveRatio(1, isData, 0, str(lumi), histo1, histo2, r_ymin = rmin, r_ymax = rmax, label = 'Ratio',  xlog = False, outputDir = WORKPATH + 'SymmetryResults/', maxYnumbers = maxY)

    #print('>>>>>>>>>> KOLMOGOROV test for ' + labela)
    #print('>>>>>>>>>> ' + str(histo1.KolmogorovTest(histo2)))
    #print('>>>>>>>>>> Chi2 test for ' + labela)
    #print('>>>>>>>>>> ' + str(histo1.Chi2Test(histo2, "WWP")))

    return 

def makeComparison(lumi, hnameList, ylog, treeList, inputdir, labelList, labela, labelb, name, isData, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', rebin = 0, maxY = False):

    ## Color palette
    galapago_palette = [r.kAzure+10, r.kAzure-4, r.kBlue-3, r.kMagenta-3, r.kMagenta+2]
    point_types = [20, 21, 22, 23, 47]

    ## Get the histogram
    histoList = []
    ymax = 0
    for t,tree in enumerate(treeList):
        histoList.append(tree.getLoopTH1F(inputdir, hnameList[t], doOF = False))
        if rebin:
            histoList[-1].Rebin(rebin)
        histoList[-1].Scale(1./histoList[-1].Integral())
        histoList[-1].SetMarkerSize(1)
        if histoList[-1].GetMaximum() > ymax:
            ymax = histoList[-1].GetMaximum() 
    
    if not maxY:
        histoList[0].SetMaximum(1.6*ymax)
    else:
        histoList[0].SetMaximum(maxY)
    histoList[0].SetMinimum(0.0)

    plot = Canvas.Canvas(name, 'png,pdf', 0.16, 0.7, 0.35, 0.89, 1, lsize = 0.035)
    for h,histo in enumerate(histoList):
        if not h:
            histo.GetXaxis().SetTitleSize(0.045)
            histo.GetYaxis().SetTitleSize(0.045)
            plot.addHisto(histo, 'P', labelList[h], 'p', galapago_palette[h], 1, h, marker = point_types[h])
        else:
            plot.addHisto(histo, 'P,SAME', labelList[h], 'p', galapago_palette[h], 1, h, marker = point_types[h])
    #plot.addLatex(0.18, 0.78, labela, font = 42, size = 0.04, align = 11)
    plot.addLatex(0.9, 0.93, labelb, font = 42, size = 0.035, align = 31)
    plot.save(1, isData, 0, str(lumi), '', outputDir = WORKPATH + 'Distribution/', maxYnumbers = maxY)

    return 

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


    ############# Load data 
    DoubleMuon2016_HIPM = []
    DoubleMuon2016_noHIPM = []
    DoubleMuon2018 = []
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016B_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016C_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016D_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016E_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016F_HIPM')
    DoubleMuon2016_noHIPM.append('DoubleMuon_Run2016F_noHIPM')
    DoubleMuon2016_noHIPM.append('DoubleMuon_Run2016G_noHIPM')
    DoubleMuon2016_noHIPM.append('DoubleMuon_Run2016H_noHIPM')
    DoubleMuon2018.append('DoubleMuon_Run2018A')
    DoubleMuon2018.append('DoubleMuon_Run2018B')
    DoubleMuon2018.append('DoubleMuon_Run2018C')
    DoubleMuon2018.append('DoubleMuon_Run2018D')

    DoubleEG2016_HIPM = []
    DoubleEG2016_noHIPM = []
    DoubleEG2017 = []
    DoubleEG2018 = []
    DoubleEG2016_HIPM.append('DoubleEG_Run2016B_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016C_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016D_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016E_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016F_HIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016F_noHIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016H_noHIPM')
    DoubleEG2017.append('DoubleEG_Run2017B')
    DoubleEG2017.append('DoubleEG_Run2017C')
    DoubleEG2017.append('DoubleEG_Run2017D')
    DoubleEG2017.append('DoubleEG_Run2017E')
    DoubleEG2017.append('DoubleEG_Run2017F')
    DoubleEG2018.append('EGamma_Run2018A')
    DoubleEG2018.append('EGamma_Run2018B')
    DoubleEG2018.append('EGamma_Run2018C')
    DoubleEG2018.append('EGamma_Run2018D')

    ############# Load background simulation
    Backgrounds_preVFP = []
    Backgrounds_preVFP.append('DYJetsToLL_M-50_preVFP')
    Backgrounds_preVFP.append('TTTo2L2Nu_preVFP')
    Backgrounds_preVFP.append('WW_preVFP')
    Backgrounds_preVFP.append('WZ_preVFP')
    Backgrounds_postVFP = []
    Backgrounds_postVFP.append('DYJetsToLL_M-50_postVFP')
    Backgrounds_postVFP.append('TTTo2L2Nu_postVFP')
    Backgrounds_postVFP.append('WW_postVFP')
    Backgrounds_postVFP.append('WZ_postVFP')
    Backgrounds_2017 = []
    Backgrounds_2017.append('DYJetsToLL_M-50_2017')
    Backgrounds_2017.append('TTTo2L2Nu_2017')
    Backgrounds_2017.append('WW_2017')
    Backgrounds_2017.append('WZ_2017')
    Backgrounds_2017.append('ZZ_2017')
    Backgrounds_2018 = []
    Backgrounds_2018.append('DYJetsToLL_M-50_2018')
    Backgrounds_2018.append('TTTo2L2Nu_2018')
    Backgrounds_2018.append('WW_2018')
    Backgrounds_2018.append('WZ_2018')
    Backgrounds_2018.append('ZZ_2018')



    ################ Define the .dat file
    filename = 'dat/Samples_cern_UltraLegacy.dat'

    lumi_total = ""

    #
    # -- Symmetry tests in Drell-Yan simulation
    #

    treeDY_preVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50_preVFP'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeDY_preVFP, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Drell-Yan (Baseline selection)', '2016 UL: preVFP', 'DYpreVFP_hEE_dPhi_corr', 0, rebin = 0, maxY = 4)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeDY_preVFP, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Drell-Yan (Baseline selection)', '2016 UL: preVFP', 'DYpreVFP_hMM_dPhi_corr', 0, rebin = 0, maxY = 4)


    treeDY_postVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50_postVFP'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeDY_postVFP, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Drell-Yan (Baseline selection)', '2016 UL: postVFP', 'DYpostVFP_hEE_dPhi_corr', 0, rebin = 0, maxY = 4)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeDY_postVFP, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Drell-Yan (Baseline selection)', '2016 UL: postVFP', 'DYpostVFP_hMM_dPhi_corr', 0, rebin = 0, maxY = 4)

    treeDY_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50_2017'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeDY_2017, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Drell-Yan (Baseline selection)', '2017 UL', 'DY2017_hEE_dPhi_corr', 0, rebin = 0, maxY = 4)

    treeDY_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-50_2018'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeDY_2018, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Drell-Yan (Baseline selection)', '2018 UL', 'DY2018_hEE_dPhi_corr', 0, rebin = 0, maxY = 4)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeDY_2018, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Drell-Yan (Baseline selection)', '2018 UL', 'DY2018_hMM_dPhi_corr', 0, rebin = 0, maxY = 4)


    #
    # -- Symmetry tests in TTbar simulation
    #

    treeTT_preVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_preVFP'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeTT_preVFP, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 't#bar{t} dilep. (Baseline selection)', '2016 UL: preVFP', 'TTpreVFP_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 0)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeTT_preVFP, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 't#bar{t} dilep. (Baseline selection)', '2016 UL: preVFP', 'TTpreVFP_hMM_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 0)

    treeTT_postVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_postVFP'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeTT_postVFP, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 't#bar{t} dilep. (Baseline selection)', '2016 UL: postVFP', 'TTpostVFP_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 0)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeTT_postVFP, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 't#bar{t} dilep. (Baseline selection)', '2016 UL: postVFP', 'TTpostVFP_hMM_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 0)

    treeTT_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_2017'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeTT_2017, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 't#bar{t} dilep. (Baseline selection)', '2017 UL', 'TT2017_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 0)

    treeTT_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_2018'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeTT_2018, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 't#bar{t} dilep. (Baseline selection)', '2018 UL', 'TT2018_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 4)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeTT_2018, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 't#bar{t} dilep. (Baseline selection)', '2018 UL', 'TT2018_hMM_dPhi_corr', 0, rebin = 0, rmin = 0.9, rmax= 1.1, maxY = 0)
    

    #
    # -- Symmetry tests in diboson simulation
    #

    treeVV_preVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_preVFP', 'WZ_preVFP'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeVV_preVFP, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'WZ, WW (Baseline selection)', '2016 UL: preVFP', 'VVpreVFP_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeVV_preVFP, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'WZ, WW (Baseline selection)', '2016 UL: preVFP', 'VVpreVFP_hMM_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)

    treeVV_postVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_postVFP', 'WZ_postVFP'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeVV_postVFP, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'WZ, WW (Baseline selection)', '2016 UL: postVFP', 'VVpostVFP_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeVV_postVFP, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'WZ, WW (Baseline selection)', '2016 UL: postVFP', 'VVpostVFP_hMM_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)


    treeVV_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_2017', 'WZ_2017', 'ZZ_2017'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeVV_2017, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'ZZ, WZ, WW (Baseline selection)', '2017 UL', 'VV2017_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)

    treeVV_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_2018', 'WZ_2018', 'ZZ_2018'], 'MC'), name = 'MC', isdata = 0 )
    makeAgreementTest(lumi_total, 'hEEsel_dPhi', 'hEEsel_dPhi_inv', 0, treeVV_2018, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'ZZ, WZ, WW (Baseline selection)', '2018 UL', 'VV2018_hEE_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest(lumi_total, 'hMMsel_dPhi', 'hMMsel_dPhi_inv', 0, treeVV_2018, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'ZZ, WZ, WW (Baseline selection)', '2018 UL', 'VV2018_hMM_dPhi_corr', 0, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)


    #
    # -- Symmetry tests in Data QCD control region
    #
    tree_EE2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_HIPM, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_HIPM, 'MC'), name = 'DATA', isdata = 1 )
    tree_EE2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_noHIPM, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_noHIPM, 'MC'), name = 'DATA', isdata = 1 )
    makeAgreementTest('', 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2016_HIPM, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', 'HIP affected runs', 'QCD2016_HIPM_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest('', 'hMMQCD_dPhi', 'hMMQCD_dPhi_inv', 0, tree_MM2016_HIPM, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region', 'HIP affected runs', 'QCD2016_HIPM_hMM_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest('', 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2016_noHIPM, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', 'HIP unaffected runs', 'QCD2016_noHIPM_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest('', 'hMMQCD_dPhi', 'hMMQCD_dPhi_inv', 0, tree_MM2016_noHIPM, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region', 'HIP unaffected runs', 'QCD2016_noHIPM_hMM_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    

    tree_EE2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'MC'), name = 'DATA', isdata = 1 )
    makeAgreementTest(41.5, 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2017, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', '', 'QCD2017_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)


    tree_EE2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2018, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'MC'), name = 'DATA', isdata = 1 )
    makeAgreementTest(59.8, 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2018, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', '', 'QCD2018_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 4)
    makeAgreementTest(59.8, 'hMMQCD_dPhi', 'hMMQCD_dPhi_inv', 0, tree_MM2018, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region', '', 'QCD2018_hMM_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 4)


    #
    # -- Symmetry tests in Data W+jets control region
    #
    makeAgreementTest('', 'hEEWjets_dPhi', 'hEEWjets_dPhi_inv', 0, tree_EE2016_HIPM, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'W+jets Control Region', 'HIP affected runs', 'Wjets2016_HIPM_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest('', 'hMMWjets_dPhi', 'hMMWjets_dPhi_inv', 0, tree_MM2016_HIPM, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Wjets Control Region', 'HIP affected runs', 'Wjets2016_HIPM_hMM_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest('', 'hEEWjets_dPhi', 'hEEWjets_dPhi_inv', 0, tree_EE2016_noHIPM, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'W+jets Control Region', 'HIP unaffected runs', 'Wjets2016_noHIPM_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    makeAgreementTest('', 'hMMWjets_dPhi', 'hMMWjets_dPhi_inv', 0, tree_MM2016_noHIPM, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'W+jets Control Region', 'HIP unaffected runs', 'Wjets2016_noHIPM_hMM_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)
    

    tree_EE2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'MC'), name = 'DATA', isdata = 1 )
    makeAgreementTest(41.5, 'hEEWjets_dPhi', 'hEEWjets_dPhi_inv', 0, tree_EE2017, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'W+jets Control Region', '', 'Wjets2017_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 0)


    tree_EE2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2018, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'MC'), name = 'DATA', isdata = 1 )
    makeAgreementTest(59.8, 'hEEWjets_dPhi', 'hEEWjets_dPhi_inv', 0, tree_EE2018, WORKPATH + opts.input, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'W+jets Control Region', '', 'Wjets2018_hEE_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 4)
    makeAgreementTest(59.8, 'hMMWjets_dPhi', 'hMMWjets_dPhi_inv', 0, tree_MM2018, WORKPATH + opts.input, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'W+jets Control Region', '', 'Wjets2018_hMM_dPhi_corr', 1, rebin = 0, rmin = 0.8, rmax= 1.2, maxY = 4)


    #
    # -- Collinearity simple distributions
    #
    makeComparison('', ['hEEsel_dPhi', 'hEEsel_dPhi', 'hEEsel_dPhi', 'hEEQCD_dPhi', 'hEEWjets_dPhi'], 0, [treeDY_postVFP, treeTT_postVFP, treeVV_postVFP, tree_EE2016_noHIPM, tree_EE2016_noHIPM], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ, simulation', 'QCD control region', 'W+jets control region'], '', '2016 UL: no HIPM / postVFP', 'hdPhi_EE_noHIPMpostVFP', 1)
    makeComparison('', ['hEEsel_dPhi', 'hEEsel_dPhi', 'hEEsel_dPhi', 'hEEQCD_dPhi', 'hEEWjets_dPhi'], 0, [treeDY_preVFP, treeTT_preVFP, treeVV_preVFP, tree_EE2016_HIPM, tree_EE2016_HIPM], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ simulation', 'QCD control region', 'W+jets control region'], '', '2016 UL: HIPM / preVFP', 'hdPhi_EE_HIPMpreVFP', 1)
    makeComparison('', ['hEEsel_dPhi', 'hEEsel_dPhi', 'hEEsel_dPhi', 'hEEQCD_dPhi', 'hEEWjets_dPhi'], 0, [treeDY_2017, treeTT_2017, treeVV_2017, tree_EE2017, tree_EE2017], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ, ZZ simulation', 'QCD control region', 'Wjets control region'], '', '2017 UL', 'hdPhi_EE_2017', 1)
    makeComparison('', ['hEEsel_dPhi', 'hEEsel_dPhi', 'hEEsel_dPhi', 'hEEQCD_dPhi', 'hEEWjets_dPhi'], 0, [treeDY_2018, treeTT_2018, treeVV_2018, tree_EE2018, tree_EE2018], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ, ZZ simulation', 'QCD control region', 'W+jets control region'], '', '2018 UL', 'hdPhi_EE_2018', 1)

    makeComparison('', ['hMMsel_dPhi', 'hMMsel_dPhi', 'hMMsel_dPhi', 'hMMQCD_dPhi', 'hMMWjets_dPhi'], 0, [treeDY_postVFP, treeTT_postVFP, treeVV_postVFP, tree_MM2016_noHIPM, tree_MM2016_noHIPM], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ, simulation', 'QCD control region', 'W+jets control region'], '', '2016 UL: no HIPM / postVFP', 'hdPhi_MM_noHIPMpostVFP', 1)
    makeComparison('', ['hMMsel_dPhi', 'hMMsel_dPhi', 'hMMsel_dPhi', 'hMMQCD_dPhi', 'hMMWjets_dPhi'], 0, [treeDY_preVFP, treeTT_preVFP, treeVV_preVFP, tree_MM2016_HIPM, tree_MM2016_HIPM], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ simulation', 'QCD control region', 'W+jets control region'], '', '2016 UL: HIPM / preVFP', 'hdPhi_MM_HIPMpreVFP', 1)
    makeComparison('', ['hMMsel_dPhi', 'hMMsel_dPhi', 'hMMsel_dPhi', 'hMMQCD_dPhi', 'hMMWjets_dPhi'], 0, [treeDY_2018, treeTT_2018, treeVV_2018, tree_MM2018, tree_MM2018], WORKPATH + opts.input, ['Drell-Yan simulation', 't#bar{t} (dilep.) simulation', 'WW, WZ, ZZ simulation', 'QCD control region', 'W+jets control region'], '', '2018 UL', 'hdPhi_MM_2018', 1)
