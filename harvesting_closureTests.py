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

    ### Get histograms
    luminosity = lumi

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
    cumA.SetMarkerStyle(20)
    cumA.SetMarkerSize(1)
    cumA.SetMarkerColor(r.kBlue)
    cumA.GetYaxis().SetTitleOffset(1.3)
    cumB.SetMarkerStyle(25)
    cumB.SetMarkerSize(1)
    cumB.SetMarkerColor(r.kRed)
    hBKG_A.SetMarkerStyle(20)
    hBKG_A.SetMarkerSize(1)
    hBKG_A.SetMarkerColor(r.kBlue)
    hBKG_B.SetMarkerStyle(25)
    hBKG_B.SetMarkerSize(1)
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
    plot = Canvas.Canvas('closuretest_'+name, 'png,pdf', 0.45, 0.66, 0.65, 0.78, 1)
    plot.addHisto(cumA, 'P', labelA, 'p', r.kBlue, 1, 0)
    plot.addHisto(cumB, 'P, SAME', labelB, 'p', r.kRed, 1, 1)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.75, 'Dielectron vertices', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.75, 'Dimuon vertices', font = 42)

    plot.addLatex(0.85, 0.56, 'Tail cumulative of x = #int^{#infty}_{x_{min}} N(x) dx', font = 42, align = 31)
    plot.addLatex(0.17, 0.81, extralabel, font = 42, align = 11, size = 0.045)
    plot.addLatex(0.9, 0.88, DATAlabel, font = 42, align = 31, size = 0.045)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/ClosureTests_' + outtag + '/'
    plot.saveRatio(1, tree.isData, ylog, luminosity, cumA, cumB, r_ymin = rmin, r_ymax = rmax, label="SR/BCR", outputDir = outdir)

    ### Main comparison
    cplot = Canvas.Canvas('plaincomparison_'+name, 'png,pdf', 0.45, 0.66, 0.65, 0.78, 1)
    cplot.addHisto(hBKG_A, 'P', labelA, 'p', r.kBlue, 1, 0)
    cplot.addHisto(hBKG_B, 'P, SAME', labelB, 'p', r.kRed, 1, 1)
    if LLlabel == 'EE':
        cplot.addLatex(0.17, 0.75, 'Dielectron vertices', font = 42)
    if LLlabel == 'MM':
        cplot.addLatex(0.17, 0.75, 'Dimuon vertices', font = 42)
    cplot.addLatex(0.17, 0.81, extralabel, font = 42, align = 11, size = 0.045)
    cplot.addLatex(0.9, 0.88, DATAlabel, font = 42, align = 31, size = 0.045)
    cplot.saveRatio(1, tree.isData, ylog, luminosity, hBKG_A, hBKG_B, label="SR/BCR", outputDir = outdir)

    

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
    Backgrounds_preVFP.append('DYJetsToLL_M-10to50_preVFP')
    Backgrounds_preVFP.append('TTTo2L2Nu_preVFP')
    Backgrounds_preVFP.append('WJetsToLNu_preVFP')
    Backgrounds_preVFP.append('WW_preVFP')
    Backgrounds_preVFP.append('WZ_preVFP')
    Backgrounds_preVFP.append('ZZ_preVFP')
    Backgrounds_postVFP = []
    Backgrounds_postVFP.append('DYJetsToLL_M-50_postVFP')
    Backgrounds_postVFP.append('DYJetsToLL_M-10to50_postVFP')
    Backgrounds_postVFP.append('TTTo2L2Nu_postVFP')
    Backgrounds_postVFP.append('WW_postVFP')
    Backgrounds_postVFP.append('WZ_postVFP')
    Backgrounds_postVFP.append('ZZ_postVFP')
    Backgrounds_2017 = []
    Backgrounds_2017.append('DYJetsToLL_M-50_2017')
    Backgrounds_2017.append('DYJetsToLL_M-10to50_2017')
    Backgrounds_2017.append('TTTo2L2Nu_2017')
    Backgrounds_2017.append('WJetsToLNu_2017')
    Backgrounds_2017.append('WW_2017')
    Backgrounds_2017.append('WZ_2017')
    Backgrounds_2017.append('ZZ_2017')
    Backgrounds_2018 = []
    Backgrounds_2018.append('DYJetsToLL_M-50_2018')
    Backgrounds_2018.append('DYJetsToLL_M-10to50_2018')
    Backgrounds_2018.append('TTTo2L2Nu_2018')
    Backgrounds_2018.append('WW_2018')
    Backgrounds_2018.append('WZ_2018')
    Backgrounds_2018.append('ZZ_2018')

    ################ Define the .dat file
    filename = 'dat/Samples_cern_UltraLegacy.dat'

    stabin = np.array([0., 1., 2., 3., 4., 6., 8., 10., 15., 20., 30., 40.])

    ################################
    ######## DoubleEG Plots ########
    ################################
    #
    # -- DY Closure
    #
    treeDY_preVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-10to50_preVFP', 'DYJetsToLL_M-50_preVFP'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeDY_preVFP.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeDY_preVFP.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    #newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2016preVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: preVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_preVFP.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeDY_preVFP.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2016preVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: preVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)


    hBKG_A = treeDY_preVFP.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeDY_preVFP.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 40., 21)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_2016preVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: preVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3.0)

    hBKG_A = treeDY_preVFP.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeDY_preVFP.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 40., 21)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_2016preVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: preVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3.0)


    treeDY_postVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-10to50_postVFP', 'DYJetsToLL_M-50_postVFP'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeDY_postVFP.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeDY_postVFP.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2016postVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: postVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_postVFP.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeDY_postVFP.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2016postVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: postVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_postVFP.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeDY_postVFP.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 40., 21)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_2016postVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: postVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3.0)

    hBKG_A = treeDY_postVFP.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeDY_postVFP.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 40., 21)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_2016postVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: postVFP', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3.0)

    treeDY_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-10to50_2017' , 'DYJetsToLL_M-50_2017'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeDY_2017.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeDY_2017.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2017', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2017 UL', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_2017.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeDY_2017.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2017_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2017 UL', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    treeDY_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-10to50_2018', 'DYJetsToLL_M-50_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeDY_2018.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeDY_2018.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2018 UL', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_2018.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeDY_2018.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2018 UL', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_2018.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeDY_2018.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 40., 21)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2018 UL', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3.0)

    hBKG_A = treeDY_2018.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeDY_2018.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 40., 21)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2018 UL', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3.0)


    treeDY_EEFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-10to50_postVFP', 'DYJetsToLL_M-50_postVFP', 'DYJetsToLL_M-10to50_2017' , 'DYJetsToLL_M-50_2017', 'DYJetsToLL_M-10to50_2018', 'DYJetsToLL_M-50_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeDY_EEFull.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeDY_EEFull.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'Run 2', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_EEFull.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeDY_EEFull.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_EE_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'Run 2', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    treeDY_MMFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['DYJetsToLL_M-10to50_preVFP', 'DYJetsToLL_M-50_preVFP', 'DYJetsToLL_M-10to50_postVFP', 'DYJetsToLL_M-50_postVFP', 'DYJetsToLL_M-10to50_2018', 'DYJetsToLL_M-50_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeDY_MMFull.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeDY_MMFull.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Run 2', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)

    hBKG_A = treeDY_MMFull.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeDY_MMFull.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'DY_MM_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeDY_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Run 2', extralabel = 'Drell-Yan (Baseline selection)', rmax = 3)


    #
    # -- ttbar Closure
    #
    treeTT_preVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_preVFP'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeTT_preVFP.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeTT_preVFP.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2016preVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: preVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_preVFP.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeTT_preVFP.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_2016preVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: preVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3.0)

    hBKG_A = treeTT_preVFP.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeTT_preVFP.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2016preVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: preVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_preVFP.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeTT_preVFP.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_2016preVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: preVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3.0)

    treeTT_postVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_postVFP'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeTT_postVFP.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeTT_postVFP.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2016postVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: postVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_postVFP.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeTT_postVFP.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_2016postVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: postVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3.0)

    hBKG_A = treeTT_postVFP.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeTT_postVFP.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2016postVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: postVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_postVFP.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeTT_postVFP.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_2016postVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: postVFP', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3.0)

    treeTT_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_2017'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeTT_2017.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeTT_2017.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2017', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2017 UL', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_2017.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeTT_2017.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2017_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #3pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2017 UL', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    treeTT_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeTT_2018.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeTT_2018.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2018 UL', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_2018.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeTT_2018.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2018 UL', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3.0)

    hBKG_A = treeTT_2018.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeTT_2018.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2018 UL', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_2018.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeTT_2018.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2018 UL', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3.0)

    treeTT_EEFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_postVFP', 'TTTo2L2Nu_2017', 'TTTo2L2Nu_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeTT_EEFull.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeTT_EEFull.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'Run 2', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_EEFull.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeTT_EEFull.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_EE_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'Run 2', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    treeTT_MMFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['TTTo2L2Nu_preVFP', 'TTTo2L2Nu_postVFP', 'TTTo2L2Nu_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeTT_MMFull.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeTT_MMFull.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Run 2', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    hBKG_A = treeTT_MMFull.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeTT_MMFull.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'TT_MM_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeTT_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Run 2', extralabel = 't#bar{t} dilep. (Baseline selection)', rmax = 3)

    #
    # -- Diboson Closure
    #
    treeVV_preVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_preVFP', 'WZ_preVFP'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeVV_preVFP.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeVV_preVFP.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2016preVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: preVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_preVFP.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeVV_preVFP.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_2016preVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: preVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3.0)

    hBKG_A = treeVV_preVFP.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeVV_preVFP.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2016preVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: preVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_preVFP.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeVV_preVFP.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_2016preVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_preVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: preVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3.0)


    treeVV_postVFP = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_postVFP', 'WZ_postVFP'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeVV_postVFP.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeVV_postVFP.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2016postVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: postVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_postVFP.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeVV_postVFP.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_2016postVFP', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: postVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3.0)

    hBKG_A = treeVV_postVFP.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeVV_postVFP.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2016postVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2016 UL: postVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_postVFP.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeVV_postVFP.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_2016postVFP_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_postVFP, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2016 UL: postVFP', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3.0)


    treeVV_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_2017', 'WZ_2017', 'ZZ_2017'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeVV_2017.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeVV_2017.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2017', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2017 UL', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_2017.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeVV_2017.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2017_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2017 UL', extralabel = 'WW, WZ, ZZ (Baseline selection)')


    treeVV_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_2018', 'WZ_2018', 'ZZ_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeVV_2018.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeVV_2018.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2018 UL', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_2018.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeVV_2018.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2018 UL', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3.0)

    hBKG_A = treeVV_2018.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeVV_2018.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '2018 UL', extralabel = 'WW, WZ, ZZ (Baseline selection)')

    hBKG_A = treeVV_2018.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeVV_2018.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '2018 UL', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3.0)

    treeVV_EEFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_postVFP', 'WZ_postVFP', 'ZZ_postVFP', 'WW_2017', 'WZ_2017', 'ZZ_2017', 'WW_2018', 'WZ_2018', 'ZZ_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeVV_EEFull.getLoopTH1F(opts.input, 'hEESR_trackIxy')
    hBKG_B = treeVV_EEFull.getLoopTH1F(opts.input, 'hEEBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'Run 2', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3)

    hBKG_A = treeVV_EEFull.getLoopTH1F(opts.input, 'hEESRTight_trackIxy')
    hBKG_B = treeVV_EEFull.getLoopTH1F(opts.input, 'hEEBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_EE_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'Run 2', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3)

    treeVV_MMFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, ['WW_preVFP', 'WZ_preVFP', 'ZZ_preVFP', 'WW_postVFP', 'WZ_postVFP', 'ZZ_postVFP', 'WW_2018', 'WZ_2018', 'ZZ_2018'], 'MC'), name = 'MC', isdata = 0 )

    hBKG_A = treeVV_MMFull.getLoopTH1F(opts.input, 'hMMSR_trackIxy')
    hBKG_B = treeVV_MMFull.getLoopTH1F(opts.input, 'hMMBCR_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Run 2', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3)

    hBKG_A = treeVV_MMFull.getLoopTH1F(opts.input, 'hMMSRTight_trackIxy')
    hBKG_B = treeVV_MMFull.getLoopTH1F(opts.input, 'hMMBCRTight_trackIxy')
    newbin = np.linspace(0., 30., 16)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'VV_MM_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = treeVV_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'Run 2', extralabel = 'WW, WZ, ZZ (Baseline selection)', rmax = 3)



    """

    #
    # -- QCD and Wjets Correlation
    #
    tree_EE2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_HIPM, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_HIPM, 'MC'), name = 'DATA', isdata = 1 )

    hBKG_A = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEQCDSR_trackIxy')
    hBKG_B = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_EE_2016HIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP affected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMQCDSR_trackIxy')
    hBKG_B = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_MM_2016HIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP affected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEQCDSRTight_trackIxy')
    hBKG_B = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_EE_2016HIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP affected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMQCDSRTight_trackIxy')
    hBKG_B = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_MM_2016HIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP affected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEWjetsSR_trackIxy')
    hBKG_B = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_EE_2016HIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP affected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMWjetsSR_trackIxy')
    hBKG_B = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_MM_2016HIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP affected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEWjetsSRTight_trackIxy')
    hBKG_B = tree_EE2016_HIPM.getLoopTH1F(opts.input, 'hEEWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_EE_2016HIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP affected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMWjetsSRTight_trackIxy')
    hBKG_B = tree_MM2016_HIPM.getLoopTH1F(opts.input, 'hMMWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_MM_2016HIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_HIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP affected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    tree_EE2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_noHIPM, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_noHIPM, 'MC'), name = 'DATA', isdata = 1 )

    hBKG_A = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEQCDSR_trackIxy')
    hBKG_B = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_EE_2016noHIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP unaffected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMQCDSR_trackIxy')
    hBKG_B = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_MM_2016noHIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP unaffected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEQCDSRTight_trackIxy')
    hBKG_B = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_EE_2016noHIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP unaffected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMQCDSRTight_trackIxy')
    hBKG_B = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'QCD_MM_2016noHIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP unaffected runs', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEWjetsSR_trackIxy')
    hBKG_B = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_EE_2016noHIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP unaffected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMWjetsSR_trackIxy')
    hBKG_B = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_MM_2016noHIPM', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP unaffected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEWjetsSRTight_trackIxy')
    hBKG_B = tree_EE2016_noHIPM.getLoopTH1F(opts.input, 'hEEWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_EE_2016noHIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = 'HIP unaffected runs', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMWjetsSRTight_trackIxy')
    hBKG_B = tree_MM2016_noHIPM.getLoopTH1F(opts.input, 'hMMWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '', name = 'Wjets_MM_2016noHIPM_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2016_noHIPM, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = 'HIP unaffected runs', extralabel = 'W+jets Control Region', rmax = 3.0)



    tree_EE2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'MC'), name = 'DATA', isdata = 1 )

    hBKG_A = tree_EE2017.getLoopTH1F(opts.input, 'hEEQCDSR_trackIxy')
    hBKG_B = tree_EE2017.getLoopTH1F(opts.input, 'hEEQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '41.5', name = 'QCD_EE_2017', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2017.getLoopTH1F(opts.input, 'hEEQCDSRTight_trackIxy')
    hBKG_B = tree_EE2017.getLoopTH1F(opts.input, 'hEEQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '41.5', name = 'QCD_EE_2017_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2017.getLoopTH1F(opts.input, 'hEEWjetsSR_trackIxy')
    hBKG_B = tree_EE2017.getLoopTH1F(opts.input, 'hEEWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '41.5', name = 'Wjets_EE_2017', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_EE2017.getLoopTH1F(opts.input, 'hEEWjetsSRTight_trackIxy')
    hBKG_B = tree_EE2017.getLoopTH1F(opts.input, 'hEEWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = '41.5', name = 'Wjets_EE_2017_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2017, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)



    tree_EE2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2018, 'MC'), name = 'DATA', isdata = 1 )
    tree_MM2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'MC'), name = 'DATA', isdata = 1 )

    hBKG_A = tree_EE2018.getLoopTH1F(opts.input, 'hEEQCDSR_trackIxy')
    hBKG_B = tree_EE2018.getLoopTH1F(opts.input, 'hEEQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'QCD_EE_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MM2018.getLoopTH1F(opts.input, 'hMMQCDSR_trackIxy')
    hBKG_B = tree_MM2018.getLoopTH1F(opts.input, 'hMMQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'QCD_MM_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2018.getLoopTH1F(opts.input, 'hEEQCDSRTight_trackIxy')
    hBKG_B = tree_EE2018.getLoopTH1F(opts.input, 'hEEQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'QCD_EE_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MM2018.getLoopTH1F(opts.input, 'hMMQCDSRTight_trackIxy')
    hBKG_B = tree_MM2018.getLoopTH1F(opts.input, 'hMMQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'QCD_MM_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EE2018.getLoopTH1F(opts.input, 'hEEWjetsSR_trackIxy')
    hBKG_B = tree_EE2018.getLoopTH1F(opts.input, 'hEEWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'Wjets_EE_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MM2018.getLoopTH1F(opts.input, 'hMMWjetsSR_trackIxy')
    hBKG_B = tree_MM2018.getLoopTH1F(opts.input, 'hMMWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'Wjets_MM_2018', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_EE2018.getLoopTH1F(opts.input, 'hEEWjetsSRTight_trackIxy')
    hBKG_B = tree_EE2018.getLoopTH1F(opts.input, 'hEEWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'Wjets_EE_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EE2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MM2018.getLoopTH1F(opts.input, 'hMMWjetsSRTight_trackIxy')
    hBKG_B = tree_MM2018.getLoopTH1F(opts.input, 'hMMWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 59.8, name = 'Wjets_MM_2018_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MM2018, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    """

    """
    tree_EEFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_noHIPM + DoubleEG2017 + DoubleEG2018, 'DATA'), name = 'DATA', isdata = 1 )

    hBKG_A = tree_EEFull.getLoopTH1F(opts.input, 'hEEWjetsSR_trackIxy')
    hBKG_B = tree_EEFull.getLoopTH1F(opts.input, 'hEEWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 117, name = 'Wjets_EE_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_EEFull.getLoopTH1F(opts.input, 'hEEWjetsSRTight_trackIxy')
    hBKG_B = tree_EEFull.getLoopTH1F(opts.input, 'hEEWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 117, name = 'Wjets_EE_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)


    hBKG_A = tree_EEFull.getLoopTH1F(opts.input, 'hEEQCDSR_trackIxy')
    hBKG_B = tree_EEFull.getLoopTH1F(opts.input, 'hEEQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 117, name = 'QCD_EE_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_EEFull.getLoopTH1F(opts.input, 'hEEQCDSRTight_trackIxy')
    hBKG_B = tree_EEFull.getLoopTH1F(opts.input, 'hEEQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 117, name = 'QCD_EE_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_EEFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'EEChannel', yshift = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)


    tree_MMFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_noHIPM + DoubleMuon2016_HIPM + DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    hBKG_A = tree_MMFull.getLoopTH1F(opts.input, 'hMMWjetsSR_trackIxy')
    hBKG_B = tree_MMFull.getLoopTH1F(opts.input, 'hMMWjetsBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 96, name = 'Wjets_MM_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)

    hBKG_A = tree_MMFull.getLoopTH1F(opts.input, 'hMMWjetsSRTight_trackIxy')
    hBKG_B = tree_MMFull.getLoopTH1F(opts.input, 'hMMWjetsBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 96, name = 'Wjets_MM_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'W+jets Control Region', rmax = 3.0)


    hBKG_A = tree_MMFull.getLoopTH1F(opts.input, 'hMMQCDSR_trackIxy')
    hBKG_B = tree_MMFull.getLoopTH1F(opts.input, 'hMMQCDBCR_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 96, name = 'QCD_MM_Full', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > #pi/2', labelA = 'Signal Region: |#Delta#Phi| < #pi/2', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    hBKG_A = tree_MMFull.getLoopTH1F(opts.input, 'hMMQCDSRTight_trackIxy')
    hBKG_B = tree_MMFull.getLoopTH1F(opts.input, 'hMMQCDBCRTight_trackIxy')
    newbin = np.linspace(0., 25., 26)
    newbin = stabin
    hBKG_A_rebin = hBKG_A.Rebin(len(newbin)-1, 'hBKG_A_rebin', newbin)
    hBKG_B_rebin = hBKG_B.Rebin(len(newbin)-1, 'hBKG_B_rebin', newbin)
    makeClosureTestMC(lumi = 96, name = 'QCD_MM_Full_Tight', hBKG_A = hBKG_A_rebin, hBKG_B = hBKG_B_rebin, ylog = True, tree = tree_MMFull, inputdir = opts.input, labelB = 'Control Region: |#Delta#Phi| > 3#pi/4', labelA = 'Signal Region: |#Delta#Phi| < #pi/4', xlabel = '|d_{xy}|/#sigma_{d}', outtag = 'MMChannel', yshift = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = 'QCD Control Region', rmax = 3.0)

    """

