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

#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###
#####
#####################

def makeBlindYieldPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, xsec = 1.0, rebin = False, lines = [], xlabel = '', outtag = '', ymax = 0.0, LLlabel = ''):


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
    hbkg.SetLineColor(r.kCyan-2) 
    hbkg.GetXaxis().SetLabelSize(0.0)
    hbkg.GetXaxis().SetNdivisions(0)
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)

    for _i, _h in enumerate(hSIS.GetHists()):
        s_histos.append(copy.deepcopy(_h))
        s_histos[-1].Scale(xsec)

    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    minValSI = min([s_histos[i].GetMinimum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])
    minVal = max([minValSI, 0.0])

    ### Set Maximum
    if not ylog:
        hbkg.SetMaximum(1.3*maxVal)
        hbkg.SetMaximum(0.7*minVal)
    else:
        hbkg.SetMaximum(1e4*maxVal)
        hbkg.SetMinimum(0.01*minVal)

    if ymax: hbkg.SetMaximum(ymax) 


    ### Canvas object
    plot = Canvas.Canvas('BlindedYields_'+hname_SI, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, lsize = 0.028)
    plot.addHisto(hbkg, 'HIST', 'Background (predicted)', 'f', '', 1, 0)
    
    ### Add signals:
    for i,_h in enumerate(s_histos):

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        legend = 'H#rightarrowSS (%d GeV,%d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    if LLlabel == 'EE':
        plot.addLatex(0.13, 0.93, 'Electron channel')
    elif LLlabel == 'MM':
        plot.addLatex(0.13, 0.93, 'Muon channel')

    ### X label
    plot.addLatex(0.44, 0.06, hname_SI.split('_')[1] + ' region', size = 0.045)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/BlindYieldsPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, xlog = False)

    

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


    ############# Signal definition
    Signals2016 = []
    Signals2016.append('HSS_400_50_1_2016')
    Signals2016.append('HSS_400_50_10_2016')
    Signals2016.append('HSS_400_50_100_2016')
    Signals2016.append('HSS_400_50_1000_2016')
    Signals2016.append('HSS_400_50_10000_2016')
    Signals2017 = []
    Signals2017.append('HSS_400_50_1_2017')
    Signals2017.append('HSS_400_50_10_2017')
    Signals2017.append('HSS_400_50_100_2017')
    Signals2017.append('HSS_400_50_1000_2017')
    Signals2017.append('HSS_400_50_10000_2017')
    Signals2018 = []
    Signals2018.append('HSS_400_50_1_2018')
    Signals2018.append('HSS_400_50_10_2018')
    Signals2018.append('HSS_400_50_100_2018')
    Signals2018.append('HSS_400_50_1000_2018')
    Signals2018.append('HSS_400_50_10000_2018')


    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1


    ############# Galapago Tree definitions
    treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2016UL.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2017UL.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2018UL.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    ################################
    ######## DoubleEG Plots ########
    ################################
       
    #### -> 2016 plots
    makeBlindYieldPlot(lumi = lumi2016, hname_SI = 'hEE_SRI', hname_bkg = 'hEE_BCRI', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', LLlabel = 'EE') 

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    #### -> 2016 plots

    #makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_mass', hname_bkg = 'hMMBCRI_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 


