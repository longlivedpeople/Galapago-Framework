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

def makeBlindedPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, treeBKG = False, treeBKGlabel = '', rebin = False, lines = [], line_ymax = False, xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', xsec = False, xlog = False):


    ### Get histograms
    luminosity = lumi

    ### Get background estimation from data
    hbkg_ = treeDATA.getLoopTH1F(inputdir, hname_bkg)

    ### Get background estimation from simulation
    if treeBKG:
        hbkgsim_ = treeBKG.getLoopTH1F(inputdir, hname_SI)
    
    ### rebinins:
    if type(rebin) != bool:
        if type(rebin) == int:
            hbkg = hbkg_.Rebin(rebin)
            if treeBKG:
                hbkgsim = hbkgsim_.Rebin(rebin)
        else:
            if len(rebin) > 1:
                hbkg = hbkg_.Rebin(len(rebin)-1, hbkg_.GetName() + '_rebined', rebin)
                if treeBKG:
                    hbkgsim = hbkgsim_.Rebin(len(rebin)-1, hbkgsim_.GetName() + '_rebined', rebin)
    else:
        hbkg = hbkg_.Clone()
        if treeBKG:
            hbkgsim = hbkgsim_.Clone()

    ### Set background histos style
    hbkg.SetFillColorAlpha(r.kCyan-6, 0.8) 
    hbkg.SetLineColor(r.kCyan-2) 
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    if treeBKG:
        hbkgsim.SetLineColor(r.kBlack) 
        hbkgsim.SetLineStyle(10)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)

    for _i, _h in enumerate(hSIS.GetHists()):

        if type(rebin) != bool:
            if type(rebin) == int:
                _h2 = _h.Rebin(rebin)
            else:
                if len(rebin) > 1:
                    _h2 = _h.Rebin(len(rebin)-1, hbkg_.GetName() + '_rebined', rebin)
        else:
            _h2 = _h.Clone()

        s_histos.append(copy.deepcopy(_h2))
        s_histos[-1].Scale(xsecs['400'])


    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])

    ### Set Maximum
    if not ylog:
        hbkg.SetMaximum(1.3*maxVal)
    else:
        hbkg.SetMaximum(100.0*maxVal)
        hbkg.SetMinimum(0.1)

    if ymax: hbkg.SetMaximum(ymax) 


    ### Canvas object
    plot = Canvas.Canvas('Blinded_'+hname_bkg, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, lsize = 0.028)
    plot.addHisto(hbkg, 'HIST', 'Background (predicted)', 'f', '', 1, 0)
    
    ### Add signals:
    for i,_h in enumerate(s_histos):

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        #legend = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
        legend = 'H#rightarrowSS (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    if treeBKG:
        plot.addHisto(hbkgsim, 'HIST, SAME', treeBKGlabel, 'l', '', 1, 2 + len(s_histos)) # Signal

    ## Lines
    if not line_ymax:
        line_ymax = hbkg.GetMaximum()
    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, line_ymax, r.kBlack, 2)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, xlog = xlog, maxYnumbers = 4)

    

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
    Signals = []
    #Signals.append('HSS_300_50_1')
    #Signals.append('HSS_300_50_10')
    #Signals.append('HSS_300_50_100')
    #Signals.append('HSS_300_50_1000')
    #Signals.append('HSS_300_50_10000')
    Signals.append('HSS_500_50_1')
    Signals.append('HSS_500_50_10')
    Signals.append('HSS_500_50_100')
    Signals.append('HSS_500_50_1000')
    Signals.append('HSS_500_50_10000')
    """
    Signals.append('HSS_1000_250_1')
    Signals.append('HSS_1000_250_10')
    Signals.append('HSS_1000_250_100')
    Signals.append('HSS_1000_250_1000')
    Signals.append('HSS_1000_250_10000')
    Signals.append('RPV_350_148_1')
    Signals.append('RPV_350_148_10')
    Signals.append('RPV_350_148_100')
    Signals.append('RPV_350_148_1000')
    Signals.append('RPV_350_148_10000')
    Signals.append('RPV_1500_494_1')
    Signals.append('RPV_1500_494_10')
    Signals.append('RPV_1500_494_100')
    Signals.append('RPV_1500_494_1000')
    Signals.append('RPV_1500_494_10000')
    """

    Signals_2016preVFP = [i + '_2016APV' for i in Signals]
    Signals_2016postVFP = [i + '_2016' for i in Signals]
    Signals_2017 = [i + '_2017' for i in Signals]
    Signals_2018 = [i + '_2018' for i in Signals]


    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1


    ############# Galapago Tree definitions
    #treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2016UL_Summer22.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
    #treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2017UL_Summer22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
    #treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2018UL_Summer22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )

    treeSI_2016    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2016UL_Fall22.dat', Signals_2016preVFP + Signals_2016postVFP, 'SI'), name = 'SI', isdata = 0, close = True)
    treeSI_2017    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2017UL_Fall22.dat', Signals_2017, 'SI'), name = 'SI', isdata = 0, close = True)
    treeSI_2018    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2018UL_Fall22.dat', Signals_2018, 'SI'), name = 'SI', isdata = 0, close = True) 


    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )


    """
    mass50_bin = np.linspace(15, 75, 13)
    mass150_bin = np.linspace(105, 200, 20)
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRIM50_mass', hname_bkg = 'hEEBCRIM50_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, rebin = mass50_bin, xlabel = '', outtag = '2018', LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRIM50_mass', hname_bkg = 'hMMBCRIM50_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, rebin = mass50_bin, xlabel = '', outtag = '2018', LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRIM150_mass', hname_bkg = 'hEEBCRIM150_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_1000_150_2018, rebin = mass150_bin, xlabel = '', outtag = '2018', LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRIM150_mass', hname_bkg = 'hMMBCRIM150_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_1000_150_2018, rebin = mass150_bin, xlabel = '', outtag = '2018', LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    """


    ################################
    ######## DoubleEG Plots ########
    ################################
    #### -> 2016 plots
    if True:
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
    if False:
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_nBSEE', hname_bkg = 'hEEBCR_nBSEE', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [15.], line_ymax = 1e7, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [40.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [25.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        #makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESRdisp_normalizedChi2_log', hname_bkg = 'hEEBCRdisp_normalizedChi2_log', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, rebin = 2, lines = [20.], line_ymax = 7e3, xlabel = '', outtag = '2016', ymax = 1e7, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hEESR_trackIxy', hname_bkg = 'hEEBCR_trackIxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [6.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    #### -> 2017 plots
    if True:
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
    if False:
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_nBSEE', hname_bkg = 'hEEBCR_nBSEE', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [15.], line_ymax = 1e7,xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [40.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [25.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        #makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRdisp_normalizedChi2_log', hname_bkg = 'hEEBCRdisp_normalizedChi2_log', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, rebin = 2, lines = [20.], line_ymax = 7e3, xlabel = '', outtag = '2017', ymax = 1e7, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_trackIxy', hname_bkg = 'hEEBCR_trackIxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, lines = [6.], xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    #### -> 2018 plots
    if True:
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_mass', hname_bkg = 'hEEBCRI_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_trackIxy', hname_bkg = 'hEEBCRI_trackIxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_trackIxy', hname_bkg = 'hEEBCRII_trackIxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_trackDxy', hname_bkg = 'hEEBCRI_trackDxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [0.03], xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_trackDxy', hname_bkg = 'hEEBCRII_trackDxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_Lxy', hname_bkg = 'hEEBCRI_Lxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_Lxy', hname_bkg = 'hEEBCRII_Lxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRI_Ixy', hname_bkg = 'hEEBCRI_Ixy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRII_Ixy', hname_bkg = 'hEEBCRII_Ixy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    if False:
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_nBSEE', hname_bkg = 'hEEBCR_nBSEE', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [15.], line_ymax = 1e7, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [40.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [25.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
        #makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESRdisp_normalizedChi2_log', hname_bkg = 'hEEBCRdisp_normalizedChi2_log', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, rebin = 2, lines = [20.], line_ymax = 7e3, xlabel = '', outtag = '2018', ymax = 1e7, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hEESR_trackIxy', hname_bkg = 'hEEBCR_trackIxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [6.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = '') 

    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    #### Signal variables in Signal Region (2016)
    if True:
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_mass', hname_bkg = 'hMMBCRI_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_trackIxy', hname_bkg = 'hMMBCRI_trackIxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_trackIxy', hname_bkg = 'hMMBCRII_trackIxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_trackDxy', hname_bkg = 'hMMBCRI_trackDxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [0.02], xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_trackDxy', hname_bkg = 'hMMBCRII_trackDxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_Lxy', hname_bkg = 'hMMBCRI_Lxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_Lxy', hname_bkg = 'hMMBCRII_Lxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRI_Ixy', hname_bkg = 'hMMBCRI_Ixy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRII_Ixy', hname_bkg = 'hMMBCRII_Ixy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 

    ### Cut variables in Signal Region (2016)
    if False:
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_nBSMM', hname_bkg = 'hMMBCR_nBSMM', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [15.], line_ymax = 1e7, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_leadingPt', hname_bkg = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [30.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_subleadingPt', hname_bkg = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [30.0], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        #makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSRdisp_normalizedChi2_log', hname_bkg = 'hMMBCRdisp_normalizedChi2_log', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, rebin = 2, lines = [20.0], line_ymax = 5e4, xlabel = '', outtag = '2016', ymax = 1e8, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True) 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_cosAlpha', hname_bkg = 'hMMBCR_cosAlpha', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = 'histograms_Summer2022_noCos/', treeSI = treeSI_2016, lines = [-0.8], line_ymax = 1e7, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2016, hname_SI = 'hMMSR_trackIxy', hname_bkg = 'hMMBCR_trackIxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, lines = [6.], xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 


    #### Signal variables in Signal Region (2018)
    if True:
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_mass', hname_bkg = 'hMMBCRI_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_trackIxy', hname_bkg = 'hMMBCRI_trackIxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_trackIxy', hname_bkg = 'hMMBCRII_trackIxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_trackDxy', hname_bkg = 'hMMBCRI_trackDxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [0.02], xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_trackDxy', hname_bkg = 'hMMBCRII_trackDxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_Lxy', hname_bkg = 'hMMBCRI_Lxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_Lxy', hname_bkg = 'hMMBCRII_Lxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRI_Ixy', hname_bkg = 'hMMBCRI_Ixy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRII_Ixy', hname_bkg = 'hMMBCRII_Ixy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 

    ### Cut variables in Signal Region (2018)
    if False:
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_nBSMM', hname_bkg = 'hMMBCR_nBSMM', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [15.], line_ymax = 3e7,xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_mass_log', hname_bkg = 'hMMBCR_mass_log', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [15.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True) 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_leadingPt', hname_bkg = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [25.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_subleadingPt', hname_bkg = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [25.0], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
        #makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSRdisp_normalizedChi2_log', hname_bkg = 'hMMBCRdisp_normalizedChi2_log', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, rebin = 2, lines = [20.0], line_ymax = 5e4, xlabel = '', outtag = '2018', ymax = 1e8, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True) 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_cosAlpha', hname_bkg = 'hMMBCR_cosAlpha', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = 'histograms_Summer2022_noCos/', treeSI = treeSI_2018, lines = [-0.99], line_ymax = 1e7, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False) 
        makeBlindedPlot(lumi = lumi2018, hname_SI = 'hMMSR_trackIxy', hname_bkg = 'hMMBCR_trackIxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, lines = [6.], xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = '') 

