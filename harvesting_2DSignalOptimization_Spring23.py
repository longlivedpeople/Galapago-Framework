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

xsecs = {}
xsecs['400'] = 5040. # pb-1

### Define palettes
bkgc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),  
                                            np.array([1.00, 0.00]),  
                                            np.array([1.00, 153/255.]), 
                                            np.array([1.00, 153./255.]), 255);

sigc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),  
                                            np.array([1.00, 204./255.]),  
                                            np.array([1.00, 0.00]), 
                                            np.array([1.00, 0.00]), 255);
bkgpalette_ = []
sigpalette_ = []
for i in range(0, 255):
    bkgpalette_.append(bkgc + i)
    sigpalette_.append(sigc + i)

bkgpalette = np.array(bkgpalette_, dtype=np.int32)
sigpalette = np.array(sigpalette_, dtype=np.int32)

#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###     - xsec to normalize the signal given in fb
###
#####
#####################

def makeBackgroundPlot2D(lumi, hname_bkg, zlog, treeDATA, inputdir, rebin = False, lines = [], xlabel = '', outtag = '', outdir = '', LLlabel = '', extralabel = '', xlog = False, ylog = False):


    ### Get histograms
    luminosity = lumi

    hbkg = treeDATA.getLoopTH2F(inputdir, hname_bkg)
    
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)
    hbkg.SetMinimum(0.1)

    ### Define palette
    nc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),  
                                              np.array([1.00, 0.00]),  
                                              np.array([1.00, 153/255.]), 
                                              np.array([1.00, 153./255.]), 255);
    r.gStyle.SetPalette(255, bkgpalette)

    ### Canvas object
    plot = Canvas.Canvas(outtag+hname_bkg, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, ww = 610, hh = 600, lsize = 0.028)
    plot.addHisto(hbkg, 'COLZ', '', '', '', 1, 0)
    
    for line in lines:
        plot.addLine(line[0], line[1], line[2], line[3], r.kRed, 2)

    ### Extralabel
    plot.addLatex(0.13, 0.93, 'Background (predicted)', font = 42)
    #plot.addLatex(0.4, 0.85, extralabel, font = 42, size = 0.03)
    plot.addLatex(0.17, 0.75, extralabel, font = 42, size = 0.03)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, zlog = zlog, is2d = False, xlog = xlog, isPrivate = True)

    
#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###     - xsec to normalize the signal given in fb
###
#####
#####################

def makeSignalPlot2D(name, lumi, hname_sig, zlog, treeSI, inputdir, rebin = False, lines = [], legend = '', xlabel = '', outtag = '', outdir = '', LLlabel = '', extralabel = '', xlog = False, ylog = False):


    ### Get histograms
    luminosity = lumi

    hsig = treeSI.getLoopTH2F(inputdir, hname_sig)
    
    hsig.GetXaxis().SetTitleSize(0.045)
    hsig.GetYaxis().SetTitleSize(0.045)

    r.gStyle.SetPalette(r.kBird)


    ### Canvas object
    plot = Canvas.Canvas('SIOnly_'+name, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, ww = 610, hh = 600, lsize = 0.028)
    plot.addHisto(hsig, 'COLZ', '', '', '', 1, 0)
    
    ### Extralabel
    plot.addLatex(0.13, 0.93, legend, font = 42, size = 0.032)
    plot.addLatex(0.17, 0.75, extralabel, font = 42, size = 0.03)

    for line in lines:
        plot.addLine(line[0], line[1], line[2], line[3], r.kRed, 2)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, zlog = zlog, xlog = xlog, is2d = False, isPrivate = True)


def countJointYields2D(hname_bkg, hname_sig, treeDATA, treeSI, inputdir, xmins, ymins):

    histo_bkg = treeDATA.getLoopTH2F(inputdir, hname_bkg)
    histo_sig = treeSI.getLoopTH2F(inputdir, hname_sig)
    xbinmax_bkg = histo_bkg.GetNbinsX() + 1 
    ybinmax_bkg = histo_bkg.GetNbinsY() + 1
    xbinmax_sig = histo_sig.GetNbinsX() + 1 
    ybinmax_sig = histo_sig.GetNbinsY() + 1
    xbins_bkg = [histo_bkg.GetXaxis().FindBin(x) for x in xmins]
    ybins_bkg = [histo_bkg.GetYaxis().FindBin(y) for y in ymins]
    xbins_sig = [histo_sig.GetXaxis().FindBin(x) for x in xmins]
    ybins_sig = [histo_sig.GetYaxis().FindBin(y) for y in ymins]

    for i in range(0, len(xmins)):
        for j in range(0, len(ymins)):
            sig_yield = histo_sig.Integral(xbins_sig[i], xbinmax_sig, ybins_sig[j], ybinmax_sig)
            bkg_yield = histo_bkg.Integral(xbins_bkg[i], xbinmax_bkg, ybins_bkg[j], ybinmax_bkg)
            print(">> For (%f, %f):  " % (xmins[i], ymins[j]), "  Bkg: %f" % bkg_yield, "  Sig: %f" % sig_yield)


def countYields2D(hname, tree, inputdir, xedges, yedges):

    histo = tree.getLoopTH2F(inputdir, hname)

    # edges to bins
    xbinmax = histo.GetNbinsX() + 1 
    ybinmax = histo.GetNbinsY() + 1
    xbins = [histo.GetXaxis().FindBin(x) for x in xedges]
    xbins.append(xbinmax)
    ybins = [histo.GetYaxis().FindBin(y) for y in yedges]
    ybins.append(ybinmax)
    for i in range(0, len(xbins)-1):
        for j in range(0, len(ybins)-1):
            print('Bin in [{0},{1},{2},{3}] : '.format(xbins[i], xbins[i+1], ybins[j], ybins[j+1]), histo.Integral(xbins[i], xbins[i+1], ybins[j], ybins[j+1]))
            

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
    r.gStyle.SetPadRightMargin(0.12)

    ############# Dat file
    filename = 'dat/Samples_cern_UltraLegacy.dat'

    ############# EG data definition
    DoubleEG2016 = []
    DoubleEG2016.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG2016.append('DoubleEG_Run2016H_noHIPM')
    DoubleEG2017 = []
    DoubleEG2017.append('DoubleEG_Run2017B')
    DoubleEG2017.append('DoubleEG_Run2017C')
    DoubleEG2017.append('DoubleEG_Run2017D')
    DoubleEG2017.append('DoubleEG_Run2017E')
    DoubleEG2017.append('DoubleEG_Run2017F')
    EGamma2018 = []
    EGamma2018.append('EGamma_Run2018A')
    EGamma2018.append('EGamma_Run2018B')
    EGamma2018.append('EGamma_Run2018C')
    EGamma2018.append('EGamma_Run2018D')


    ############# Muon data definition
    DoubleMuon2016 = []
    DoubleMuon2016.append('DoubleMuon_Run2016B_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016C_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016D_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016E_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016F_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016F_noHIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016G_noHIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016H_noHIPM')
    DoubleMuon2018 = []
    DoubleMuon2018.append('DoubleMuon_Run2018A')
    DoubleMuon2018.append('DoubleMuon_Run2018B')
    DoubleMuon2018.append('DoubleMuon_Run2018C')
    DoubleMuon2018.append('DoubleMuon_Run2018D')

    ############# Luminosity definition
    lumi2016_MM = 35.9 # fb-1
    lumi2016_EE = 16.2 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018_MM = 59.8 # fb-1
    lumi2018_EE = 54.5 # fb-1


    ############# Galapago Tree definitions
    EE_SRIa_lines = []
    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )



    ########################################################
    ######## Background optimization (Lxy/dxy bins) ########
    ########################################################

    www = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/SignalRegionOptimization/Spring23/2DPlots_nLL1_for2MassBins'

    #### -> Electron plots
    EE_SRIa_lines = []
    EE_SRIa_lines.append([3e-2, 20, 1e2, 20])
    EE_SRIa_lines.append([3e-2, 6, 1e2, 6])
    EE_SRIa_lines.append([3e-2, 3, 1e2, 3])
    EE_SRIa_lines.append([3e-2, 3, 3e-2, 1e4])
    makeBackgroundPlot2D(lumi = lumi2016_EE, hname_bkg = 'hEEBCRIa_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, lines = EE_SRIa_lines, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'm_{ee} < 81 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRIa_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, lines = EE_SRIa_lines, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'm_{ee} < 81 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True)
    makeBackgroundPlot2D(lumi = lumi2018_EE, hname_bkg = 'hEEBCRIa_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, lines = EE_SRIa_lines, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'm_{ee} < 81 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 

    EE_SRIb_lines = []
    EE_SRIb_lines.append([2e-2, 15, 1e2, 15])
    EE_SRIb_lines.append([2e-2, 6, 1e2, 6])
    EE_SRIb_lines.append([2e-2, 3, 1e2, 3])
    EE_SRIb_lines.append([2e-2, 3, 2e-2, 1e4])
    makeBackgroundPlot2D(lumi = lumi2016_EE, hname_bkg = 'hEEBCRIb_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, lines = EE_SRIb_lines, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'm_{ee} > 101 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRIb_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, lines = EE_SRIb_lines, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'm_{ee} > 101 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True)
    makeBackgroundPlot2D(lumi = lumi2018_EE, hname_bkg = 'hEEBCRIb_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, lines = EE_SRIb_lines, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'm_{ee} > 101 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 

    #### -> Muon plots
    MM_SRIa_lines = []
    MM_SRIa_lines.append([2e-2, 3, 1e2, 3])
    MM_SRIa_lines.append([2e-2, 9, 1e2, 9])
    MM_SRIa_lines.append([0.2, 3, 0.2, 1e4])
    MM_SRIa_lines.append([2e-2, 3, 2e-2, 1e4])
    makeBackgroundPlot2D(lumi = lumi2018_MM, hname_bkg = 'hMMBCRIa_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, lines = MM_SRIa_lines, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'm_{#mu#mu} < 81 GeV,  N_{#mu#mu} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 
    makeBackgroundPlot2D(lumi = lumi2016_MM, hname_bkg = 'hMMBCRIa_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, lines = MM_SRIa_lines, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'm_{#mu#mu} < 81 GeV,  N_{#mu#mu} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 

    MM_SRIb_lines = []
    MM_SRIb_lines.append([2e-2, 3, 1e2, 3])
    MM_SRIb_lines.append([2e-2, 9, 1e2, 9])
    MM_SRIb_lines.append([7e-2, 3, 7e-2, 1e4])
    MM_SRIb_lines.append([2e-2, 3, 2e-2, 1e4])
    makeBackgroundPlot2D(lumi = lumi2018_MM, hname_bkg = 'hMMBCRIb_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, lines = MM_SRIb_lines, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'm_{#mu#mu} < 101 GeV,  N_{#mu#mu} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 
    makeBackgroundPlot2D(lumi = lumi2016_MM, hname_bkg = 'hMMBCRIb_Lxy_trackIxy_log', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, lines = MM_SRIb_lines, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'm_{#mu#mu} > 101 GeV,  N_{#mu#mu} = 1,  |#Delta#Phi| < #pi/4', outdir = www, ylog = True, xlog = True) 


    ####################################################
    ######## Signal optimization (Lxy/dxy bins) ########
    ####################################################
    Signals = []
    Signals.append('HSS_500_50_1')
    Signals.append('HSS_500_50_10')
    Signals.append('HSS_500_50_100')
    Signals.append('HSS_500_50_1000')
    Signals.append('HSS_500_150_1')
    Signals.append('HSS_500_150_10')
    Signals.append('HSS_500_150_100')
    Signals.append('HSS_500_150_1000')
    #Signals.append('RPV_350_148_1')
    #Signals.append('RPV_350_148_10')
    #Signals.append('RPV_350_148_100')
    #Signals.append('RPV_350_148_1000')

    for signal in Signals:
        label = signal + '_2018'
        values = signal.split('_')
        point = '({0} GeV, {1} GeV, {2} mm)'.format(values[1], values[2], values[3])
        llegend = 'H#rightarrowSS ' + point if 'HSS' in signal else 'RPV ' + point
        treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2018UL_Fall22.dat', [label], 'SI'), name = 'SI', isdata = 0 )


        makeSignalPlot2D(name = 'hMMSRIa_'+label, lumi = lumi2018_MM, hname_sig = 'hMMSRIa_Lxy_trackIxy_log', zlog = True, treeSI = treeSI_2018, inputdir = opts.input, legend = llegend, lines = MM_SRIa_lines, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'm_{#mu#mu} < 81 GeV,  N_{#mu#mu} = 1,  |#Delta#Phi| < #pi/4', outdir = www, xlog = True, ylog = True) 
        makeSignalPlot2D(name = 'hEESRIa_'+label, lumi = lumi2018_EE, hname_sig = 'hEESRIa_Lxy_trackIxy_log', zlog = True, treeSI = treeSI_2018, inputdir = opts.input, legend = llegend, lines = EE_SRIa_lines, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'm_{ee} < 81 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, xlog = True, ylog = True) 

        makeSignalPlot2D(name = 'hMMSRIb_'+label, lumi = lumi2018_MM, hname_sig = 'hMMSRIb_Lxy_trackIxy_log', zlog = True, treeSI = treeSI_2018, inputdir = opts.input, legend = llegend, lines = MM_SRIb_lines, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'm_{#mu#mu} > 101 GeV,  N_{#mu#mu} = 1,  |#Delta#Phi| < #pi/4', outdir = www, xlog = True, ylog = True) 
        makeSignalPlot2D(name = 'hEESRIb_'+label, lumi = lumi2018_EE, hname_sig = 'hEESRIb_Lxy_trackIxy_log', zlog = True, treeSI = treeSI_2018, inputdir = opts.input, legend = llegend, lines = EE_SRIb_lines, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'm_{ee} > 101 GeV,  N_{ee} = 1,  |#Delta#Phi| < #pi/4', outdir = www, xlog = True, ylog = True) 



