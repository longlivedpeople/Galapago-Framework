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
from include.Utils import XSECS
from include.galapagoStyle import sigpalette, gcolors, dcolors


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

def makeBlindedPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, treeBKG = False, treeBKGlabel = '', rebin = False, lines = [], line_ymax = False, xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', xsec = False, xlog = False, text = False, outpath = '', drawZero = True):


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

    ### Uncertainty bkg
    hunc = hbkg.Clone("uncertainty")
    hunc.Reset()
    hunc.Sumw2()
    sys = 0.1 # Overwritten
    for n in range(1, hbkg.GetNbinsX() + 1):
        hunc.SetBinContent(n, hbkg.GetBinContent(n))
        hunc.SetBinError(n, math.sqrt(sys*hbkg.GetBinContent(n)*sys*hbkg.GetBinContent(n) + hbkg.GetBinError(n)*hbkg.GetBinError(n)))
        if drawZero and  hunc.GetBinContent(n) < 1.:
            hunc.SetBinError(n, 1.8)
            if ylog: hunc.SetBinContent(n, 0.1)
        #print(hunc.GetBinContent(n), hunc.GetBinError(n))


    ### Set background histos style
    hbkg.SetFillColorAlpha(r.kCyan-6, 0.8) 
    hbkg.SetLineColor(r.kCyan-2) 
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)
    hunc.SetFillStyle(3244)
    hunc.SetFillColor(r.kCyan+3)
    hunc.SetLineColor(r.kCyan+3)
    hunc.SetMarkerSize(0)

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
        s_histos[-1].Scale(XSECS[_h2.GetTitle()])


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
    plot = Canvas.Canvas('Blinded_'+hname_bkg, 'png,pdf', 0.15, 0.65, 0.45, 0.89, 1, lsize = 0.03)
    if text:
        plot.addHisto(hbkg, 'HIST, TEXT', 'Background (predicted)', 'f', '', 1, 0)
    else:
        plot.addHisto(hbkg, 'HIST', 'Background (predicted)', 'f', '', 1, 0)
        plot.addHisto(hunc, 'E2, SAME', 'Background uncertainty', 'f', '', 1, 0)
    
    ### Add signals:
    colors = [r.kRed, r.kOrange, r.kGreen+2, r.kBlue, r.kMagenta]
    for i,_h in enumerate(s_histos):

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        if 'HSS' in _h.GetTitle():
            legend = 'H#rightarrowSS (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        else:
            legend = 'RPV (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        if text:
            plot.addHisto(_h, 'HIST TEXT, SAME', legend, 'l', colors[i], 1, i+1) # Signal
        else:
            plot.addHisto(_h, 'HIST, SAME', legend, 'l', colors[i], 1, i+1) # Signal

    if treeBKG:
        plot.addHisto(hbkgsim, 'HIST, SAME', treeBKGlabel, 'l', '', 1, 2 + len(s_histos)) # Signal

    if LLlabel == 'EE':
        plot.addLatex(0.7, 0.86, 'e^{+}e^{-} channel', font = 42, size = 0.035)
    if LLlabel == 'MM':
        plot.addLatex(0.7, 0.86, '#mu^{+}#mu^{-} channel', font = 42, size = 0.035)

    ## Lines
    if not line_ymax:
        line_ymax = hbkg.GetMaximum()
    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, line_ymax, r.kBlack, 2)

    if extralabel:
        plot.addLatex(0.5, 0.6, extralabel, font = 42, align = 22, size = 0.029)

    ### Save it
    #outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
    outdir = outpath + '/SRPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, xlog = xlog, maxYnumbers = False, is2d = True, isPrivate = True)

    

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
    filename = 'dat/Samples_cern_UltraLegacy_Spring23.dat'

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


    ############# Signal definition
    Signals = []
    Signals.append('HSS_300_50_100')
    Signals.append('HSS_400_150_100')
    Signals.append('HSS_500_50_100')
    Signals.append('HSS_800_50_100')
    Signals.append('HSS_1000_250_100')
    #Signals.append('RPV_350_148_100')
    #Signals.append('RPV_1500_494_100')


    Signals_2016preVFP = [i + '_2016APV' for i in Signals]
    Signals_2016postVFP = [i + '_2016' for i in Signals]
    Signals_2017 = [i + '_2017' for i in Signals]
    Signals_2018 = [i + '_2018' for i in Signals]


    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1


    ############# Galapago Tree definitions

    treeSI_2016postVFP    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2016UL_Fall22.dat', Signals_2016postVFP, 'SI'), name = 'SI', isdata = 0, close = True)
    treeSI_2016    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2016UL_Fall22.dat', Signals_2016preVFP + Signals_2016postVFP, 'SI'), name = 'SI', isdata = 0, close = True)
    treeSI_2017    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2017UL_Fall22.dat', Signals_2017, 'SI'), name = 'SI', isdata = 0, close = True)
    treeSI_2018    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2018UL_Fall22.dat', Signals_2018, 'SI'), name = 'SI', isdata = 0, close = True) 


    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    www = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/Vertex-selection/Spring23/'

    nomass_ee     = www + 'histograms_Spring23_VertexPlots_Electrons_Mass'
    chi2_ee       = www + 'histograms_Spring23_VertexPlots_Electrons_Chi2'
    pt_ee         = www + 'histograms_Spring23_VertexPlots_Electrons_Pt'

    nomass_mm     = www + 'histograms_Spring23_VertexPlots_Muons_Mass'
    nocosAlpha_mm = www + 'histograms_Spring23_VertexPlots_Muons_cosAlpha'
    chi2_mm       = www + 'histograms_Spring23_VertexPlots_Muons_Chi2'
    pt_mm         = www + 'histograms_Spring23_VertexPlots_Muons_pt'
    pt_mm         = www + 'histograms_Spring23_VertexPlots_Muons_Pt'


    ################################
    ######## DoubleEG Plots ########
    ################################
    #### -> 2016 plots
    if True:
        makeBlindedPlot(lumi = 16.2, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016, inputdir = nomass_ee, treeSI = treeSI_2016postVFP, lines = [15.], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 16.2, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2016, inputdir = pt_ee, treeSI = treeSI_2016postVFP, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 16.2, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2016, inputdir = pt_ee, treeSI = treeSI_2016postVFP, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 16.2, hname_SI = 'hEESR_normalizedChi2_log', hname_bkg = 'hEEBCR_normalizedChi2_log', ylog = True, treeDATA = treeDATA_EG2016, inputdir = chi2_ee, treeSI = treeSI_2016postVFP, lines = [20.], line_ymax = 1e4, xlabel = '', outtag = '2016', ymax = 1e8, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True, outpath = www) 
        makeBlindedPlot(lumi = 16.2, hname_SI = 'hEESR_trackIxy', hname_bkg = 'hEEBCR_trackIxy', ylog = True, treeDATA = treeDATA_EG2016, inputdir = pt_ee, treeSI = treeSI_2016, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 


        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = nomass_mm, treeSI = treeSI_2016, lines = [15.], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_cosAlpha', hname_bkg = 'hMMBCR_cosAlpha', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = nocosAlpha_mm, treeSI = treeSI_2016, lines = [-0.8], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_leadingPt', hname_bkg = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = pt_mm, treeSI = treeSI_2016, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_subleadingPt', hname_bkg = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = pt_mm, treeSI = treeSI_2016, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_normalizedChi2_log', hname_bkg = 'hMMBCR_normalizedChi2_log', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = chi2_mm, treeSI = treeSI_2016, lines = [20.], line_ymax = 1e4, xlabel = '', outtag = '2016', ymax = 1e8, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True, outpath = www) 
        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_trackIxy', hname_bkg = 'hMMBCR_trackIxy', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = pt_mm, treeSI = treeSI_2016, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 35.9, hname_SI = 'hMMSR_dR', hname_bkg = 'hMMBCR_dR', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = pt_mm, treeSI = treeSI_2016, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2016', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 


    #### -> 2017 plots
    if True:
        makeBlindedPlot(lumi = 41.5, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2017, inputdir = nomass_ee, treeSI = treeSI_2017, lines = [100.], line_ymax = 1e6, xlabel = '', outtag = '2017', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 41.5, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2017, inputdir = pt_ee, treeSI = treeSI_2017, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2017', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 41.5, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2017, inputdir = pt_ee, treeSI = treeSI_2017, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2017', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 41.5, hname_SI = 'hEESR_normalizedChi2_log', hname_bkg = 'hEEBCR_normalizedChi2_log', ylog = True, treeDATA = treeDATA_EG2017, inputdir = chi2_ee, treeSI = treeSI_2017, lines = [20.], line_ymax = 1e4, xlabel = '', outtag = '2017', ymax = 1e8, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True, outpath = www) 
        makeBlindedPlot(lumi = 41.5, hname_SI = 'hEESR_trackIxy', hname_bkg = 'hEEBCR_trackIxy', ylog = True, treeDATA = treeDATA_EG2017, inputdir = pt_ee, treeSI = treeSI_2017, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2017', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 

    #### -> 2018 plots
    if True:
        makeBlindedPlot(lumi = 54.5, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = nomass_ee, treeSI = treeSI_2018, lines = [15.], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 54.5, hname_SI = 'hEESR_leadingEt', hname_bkg = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2018, inputdir = pt_ee, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 54.5, hname_SI = 'hEESR_subleadingEt', hname_bkg = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2018, inputdir = pt_ee, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 54.5, hname_SI = 'hEESR_normalizedChi2_log', hname_bkg = 'hEEBCR_normalizedChi2_log', ylog = True, treeDATA = treeDATA_EG2018, inputdir = chi2_ee, treeSI = treeSI_2018, lines = [20.], line_ymax = 1e4, xlabel = '', outtag = '2018', ymax = 1e8, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True, outpath = www) 
        makeBlindedPlot(lumi = 54.5, hname_SI = 'hEESR_trackIxy', hname_bkg = 'hEEBCR_trackIxy', ylog = True, treeDATA = treeDATA_EG2018, inputdir = pt_ee, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e10, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 


        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = nomass_mm, treeSI = treeSI_2018, lines = [15.], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_cosAlpha', hname_bkg = 'hMMBCR_cosAlpha', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = nocosAlpha_mm, treeSI = treeSI_2018, lines = [-0.9], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_leadingPt', hname_bkg = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = pt_mm, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_subleadingPt', hname_bkg = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = pt_mm, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www, drawZero = False) 
        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_normalizedChi2_log', hname_bkg = 'hMMBCR_normalizedChi2_log', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = chi2_mm, treeSI = treeSI_2018, lines = [20.], line_ymax = 1e4, xlabel = '', outtag = '2018', ymax = 1e8, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True, outpath = www) 
        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_dR', hname_bkg = 'hMMBCR_dR', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = pt_mm, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 
        makeBlindedPlot(lumi = 59.8, hname_SI = 'hMMSR_trackIxy', hname_bkg = 'hMMBCR_trackIxy', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = pt_mm, treeSI = treeSI_2018, lines = [], line_ymax = 1e6, xlabel = '', outtag = '2018', ymax = 1e11, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = False, outpath = www) 



