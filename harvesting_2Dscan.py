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

def makeBackgroundPlot2D(lumi, hname_bkg, zlog, treeDATA, inputdir, rebin = False, lines = [], xlabel = '', outtag = '', outdir = '', LLlabel = '', extralabel = ''):


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
        plot.addLine(line, hbkg.GetMinimum(), line, hbkg.GetMaximum(), r.kBlack, 2)

    ### Extralabel
    plot.addLatex(0.13, 0.93, 'Background (predicted)', font = 42)
    plot.addLatex(0.5, 0.85, extralabel, font = 42, size = 0.03)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, 0, luminosity, '', outputDir = outdir, zlog = zlog, is2d = True)

    
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

def makeSignalPlot2D(lumi, hname_sig, zlog, treeSI, inputdir, rebin = False, lines = [], legend = '', xlabel = '', outtag = '', outdir = '', LLlabel = '', extralabel = ''):


    ### Get histograms
    luminosity = lumi

    hsig = treeSI.getLoopTH2F(inputdir, hname_sig)
    
    hsig.GetXaxis().SetTitleSize(0.045)
    hsig.GetYaxis().SetTitleSize(0.045)

    r.gStyle.SetPalette(255, sigpalette)


    ### Canvas object
    plot = Canvas.Canvas('SIOnly_'+hname_sig, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, ww = 610, hh = 600, lsize = 0.028)
    plot.addHisto(hsig, 'COLZ', '', '', '', 1, 0)
    
    ### Extralabel
    plot.addLatex(0.13, 0.93, legend, font = 42, size = 0.032)
    plot.addLatex(0.5, 0.85, extralabel, font = 42, size = 0.03)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, 0, luminosity, '', outputDir = outdir, zlog = zlog, is2d = True)


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
    DoubleEG2016.append('DoubleEG_Run2016B_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016C_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016D_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016E_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016F_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016F_noHIPM')
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
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1


    ############# Galapago Tree definitions

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    #####################################################
    ######## Background optimization (mass bins) ########
    #####################################################
    """   
    #### -> Electron plots

    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hEEBCRIM50_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'Mass range: [15, 76] GeV') 
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hEEBCRIM150_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'Mass range: [106, 200] GeV') 
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hEEBCRIM250_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'Mass range: [200, 300] GeV') 
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hEEBCRIM350_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'Mass range: [300, #inf) GeV') 

    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRIM50_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'Mass range: [15, 76] GeV') 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRIM150_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'Mass range: [106, 200] GeV') 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRIM250_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'Mass range: [200, 300] GeV') 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRIM350_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'Mass range: [300, #inf) GeV') 

    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEEBCRIM50_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'Mass range: [15, 76] GeV') 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEEBCRIM150_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'Mass range: [106, 200] GeV') 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEEBCRIM250_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'Mass range: [200, 300] GeV') 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEEBCRIM350_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'Mass range: [300, #inf) GeV') 

    #### -> Muon plots

    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMMBCRIM50_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'Mass range: [15, 76] GeV') 
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMMBCRIM150_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'Mass range: [106, 200] GeV') 
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMMBCRIM250_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'Mass range: [200, 300] GeV') 
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMMBCRIM350_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'Mass range: [300, #inf) GeV') 

    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMMBCRIM50_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'Mass range: [15, 76] GeV') 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMMBCRIM150_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'Mass range: [106, 200] GeV') 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMMBCRIM250_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'Mass range: [200, 300] GeV') 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMMBCRIM350_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'Mass range: [300, #inf) GeV') 


    ### Signal optimization
    treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2018UL.dat', ['HSS_400_50_100_2018'], 'SI'), name = 'SI', isdata = 0 )
    makeSignalPlot2D(lumi = lumi2018, hname_sig = 'hMMSRIM50_Lxy_trackIxy', zlog = True, treeSI = treeSI_2018, inputdir = opts.input, legend = 'H#rightarrowSS (400 GeV, 50 GeV, 100 mm)', xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'Mass range: [15, 76] GeV') 
    makeSignalPlot2D(lumi = lumi2018, hname_sig = 'hEESRIM50_Lxy_trackIxy', zlog = True, treeSI = treeSI_2018, inputdir = opts.input, legend = 'H#rightarrowSS (400 GeV, 50 GeV, 100 mm)', xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'Mass range: [15, 76] GeV') 

    countJointYields2D(hname_bkg = 'hMMBCRIM50_Lxy_trackIxy', hname_sig = 'hMMSRIM50_Lxy_trackIxy', treeDATA = treeDATA_Mu2018, treeSI = treeSI_2018, inputdir = opts.input, xmins = [0.0, 0.1], ymins = [5.0, 10.0, 14.0]) 
    countJointYields2D(hname_bkg = 'hEEBCRIM50_Lxy_trackIxy', hname_sig = 'hEESRIM50_Lxy_trackIxy', treeDATA = treeDATA_EG2018, treeSI = treeSI_2018, inputdir = opts.input, xmins = [0.0, 0.1], ymins = [5.0, 10.0, 14.0]) 
    """


    ########################################################
    ######## Background optimization (Lxy/dxy bins) ########
    ########################################################

    www = '/eos/user/f/fernance/www/LLP/SignalRegion-optimization/2DPlots_nLL1'

    #### -> Electron plots
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hEEBCRI_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'EE', extralabel = 'Off-Z region, N_{DV} = 1', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEEBCRI_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = 'Off-Z region, N_{DV} = 1', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEEBCRI_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = 'Off-Z region, N_{DV} = 1', outdir = www) 

    #### -> Muon plots
    makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMMBCRI_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = 'Off-Z region, N_{DV} = 1', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMMBCRI_Lxy_trackIxy', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = 'Off-Z region, N_{DV} = 1', outdir = www) 








