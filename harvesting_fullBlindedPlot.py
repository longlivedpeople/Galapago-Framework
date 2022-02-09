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
    hbkg.SetLineColor(r.kCyan-2) 
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)

    for _i, _h in enumerate(hSIS.GetHists()):
        #_h.Scale(lumi/35.87)
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
    plot = Canvas.Canvas('Blinded_'+hname_bkg, 'png,pdf', 0.15, 0.65, 0.6, 0.89, 1)
    plot.addHisto(hbkg, 'HIST', 'Background (Data-driven)', 'f', '', 1, 0)
    
    ### Add signals:
    for i,_h in enumerate(s_histos):

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        legend = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, hbkg.GetMaximum(), r.kBlack, 2)

    ### Extralabel
    #plot.addLatex(0.17, 0.8, extralabel, font = 62)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
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
    treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'signals_2016.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'signals_2017.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'signals_2018.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    ################################
    ######## DoubleEG Plots ########
    ################################
       
    #### -> 2016 plots
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRI_mass', hname_bkg = 'hEEBCRI_mass', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRI_trackIxy', hname_bkg = 'hEEBCRI_trackIxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRII_trackIxy', hname_bkg = 'hEEBCRII_trackIxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRI_trackDxy', hname_bkg = 'hEEBCRI_trackDxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [0.03], xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRII_trackDxy', hname_bkg = 'hEEBCRII_trackDxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRI_Lxy', hname_bkg = 'hEEBCRI_Lxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRII_Lxy', hname_bkg = 'hEEBCRII_Lxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRI_Ixy', hname_bkg = 'hEEBCRI_Ixy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRII_Ixy', hname_bkg = 'hEEBCRII_Ixy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    ### Cut variables in Signal Region
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_nBSEE', hname_bkg = 'hEEBCR_nBSEE', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [15.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [40.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [25.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_normalizedChi2', hname_bkg = 'hEEBCR_normalizedChi2', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, rebin = np.linspace(0.0, 12.0, 13), lines = [10.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    #### -> 2017 plots
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRI_mass', hname_bkg = 'hEEBCRI_mass', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRI_trackIxy', hname_bkg = 'hEEBCRI_trackIxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRII_trackIxy', hname_bkg = 'hEEBCRII_trackIxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRI_trackDxy', hname_bkg = 'hEEBCRI_trackDxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [0.03], xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRII_trackDxy', hname_bkg = 'hEEBCRII_trackDxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRI_Ixy', hname_bkg = 'hEEBCRI_Ixy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRII_Ixy', hname_bkg = 'hEEBCRII_Ixy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRI_Lxy', hname_bkg = 'hEEBCRI_Lxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRII_Lxy', hname_bkg = 'hEEBCRII_Lxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    ### Cut variables in Signal Region
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_nBSEE', hname_bkg = 'hEEBCR_nBSEE', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [15.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [40.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [25.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_normalizedChi2', hname_bkg = 'hEEBCR_normalizedChi2', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, rebin = np.linspace(0.0, 12.0, 13), lines = [10.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    #### -> 2018 plots
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_mass', hname_bkg = 'hEEBCRI_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_trackIxy', hname_bkg = 'hEEBCRI_trackIxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_trackIxy', hname_bkg = 'hEEBCRII_trackIxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_trackDxy', hname_bkg = 'hEEBCRI_trackDxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [0.03], xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_trackDxy', hname_bkg = 'hEEBCRII_trackDxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_Lxy', hname_bkg = 'hEEBCRI_Lxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_Lxy', hname_bkg = 'hEEBCRII_Lxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_Ixy', hname_bkg = 'hEEBCRI_Ixy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_Ixy', hname_bkg = 'hEEBCRII_Ixy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_nBSEE', hname_bkg = 'hEEBCR_nBSEE', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [15.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [40.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [25.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_normalizedChi2', hname_bkg = 'hEEBCR_normalizedChi2', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, rebin = np.linspace(0.0, 12.0, 13), lines = [10.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    #### -> 2016 plots

    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_mass', hname_bkg = 'hMMBCRI_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_trackIxy', hname_bkg = 'hMMBCRI_trackIxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_trackIxy', hname_bkg = 'hMMBCRII_trackIxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_trackDxy', hname_bkg = 'hMMBCRI_trackDxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [0.02], xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_trackDxy', hname_bkg = 'hMMBCRII_trackDxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_Lxy', hname_bkg = 'hMMBCRI_Lxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_Lxy', hname_bkg = 'hMMBCRII_Lxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_Ixy', hname_bkg = 'hMMBCRI_Ixy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_Ixy', hname_bkg = 'hMMBCRII_Ixy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 


    ### Cut variables in Signal Region
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_nBSMM', hname_bkg = 'hMMBCR_nBSMM', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [15.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_leadingPt', hname_bkg = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [30.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_subleadingPt', hname_bkg = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [30.0], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_normalizedChi2', hname_bkg = 'hMMBCR_normalizedChi2', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, rebin = np.linspace(0.0, 12.0, 13), lines = [10.0], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_cosAlpha', hname_bkg = 'hMMBCR_cosAlpha', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [-0.8], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 


    #### -> 2018 plots

    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_mass', hname_bkg = 'hMMBCRI_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_trackIxy', hname_bkg = 'hMMBCRI_trackIxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_trackIxy', hname_bkg = 'hMMBCRII_trackIxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_trackDxy', hname_bkg = 'hMMBCRI_trackDxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [0.02], xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_trackDxy', hname_bkg = 'hMMBCRII_trackDxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_Lxy', hname_bkg = 'hMMBCRI_Lxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_Lxy', hname_bkg = 'hMMBCRII_Lxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_Ixy', hname_bkg = 'hMMBCRI_Ixy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_Ixy', hname_bkg = 'hMMBCRII_Ixy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 


    ### Cut variables in Signal Region
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_nBSMM', hname_bkg = 'hMMBCR_nBSMM', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [15.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_leadingPt', hname_bkg = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [25.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_subleadingPt', hname_bkg = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [25.0], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_normalizedChi2', hname_bkg = 'hMMBCR_normalizedChi2', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, rebin = np.linspace(0.0, 12.0, 13), lines = [10.0], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_cosAlpha', hname_bkg = 'hMMBCR_cosAlpha', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [-0.9], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 


