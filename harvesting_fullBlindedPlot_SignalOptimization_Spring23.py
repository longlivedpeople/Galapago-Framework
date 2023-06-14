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
from include.Utils import XSECS, makeBlindedPlot
from include.galapagoStyle import sigpalette, gcolors, dcolors


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
    lumi2016_EE = 16.2 # fb-1
    lumi2016_MM = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018_MM = 59.8 # fb-1
    lumi2018_EE = 54.5 # fb-1


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

    www = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/SignalRegionOptimization/Spring23/'

    input_ee     = www + 'histograms_Spring23_SROptimization_2MassBins'
    input_mm     = www + 'histograms_Spring23_SROptimization_2MassBins'


    ################################
    ######## DoubleEG Plots ########
    ################################
    #### -> 2016 plots
    if True:
        makeBlindedPlot(lumi = lumi2016_EE, hname_SI = 'hEESR_nEE', hname_bkg = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2016, inputdir = input_ee, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1])
        makeBlindedPlot(lumi = lumi2016_EE, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016, inputdir = input_ee, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1,2])
        makeBlindedPlot(lumi = lumi2016_EE, hname_SI = 'hEESRI_dPhi_scan', hname_bkg = 'hEEBCRI_dPhi_scan', ylog = True, treeDATA = treeDATA_EG2016, inputdir = input_ee, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'N_{ee} = 1, Off-Z region', outpath = www, cutBins = [1,2])

        makeBlindedPlot(lumi = lumi2016_MM, hname_SI = 'hMMSR_nMM', hname_bkg = 'hMMBCR_nMM', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = input_mm, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1])
        makeBlindedPlot(lumi = lumi2016_MM, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = input_mm, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1,2])
        makeBlindedPlot(lumi = lumi2016_MM, hname_SI = 'hMMSRI_dPhi_scan', hname_bkg = 'hMMBCRI_dPhi_scan', ylog = True, treeDATA = treeDATA_Mu2016, inputdir = input_mm, treeSI = treeSI_2016, xlabel = '', outtag = '2016', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = 'N_{#mu#mu} = 1, Off-Z region', outpath = www, cutBins = [1,2])


    #### -> 2017 plots
    if True:
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_nEE', hname_bkg = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2017, inputdir = input_ee, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1])
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2017, inputdir = input_ee, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1,2])
        makeBlindedPlot(lumi = lumi2017, hname_SI = 'hEESRI_dPhi_scan', hname_bkg = 'hEEBCRI_dPhi_scan', ylog = True, treeDATA = treeDATA_EG2017, inputdir = input_ee, treeSI = treeSI_2017, xlabel = '', outtag = '2017', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'N_{ee} = 1, Off-Z region', outpath = www, cutBins = [1,2])


    #### -> 2018 plots
    if True:
        makeBlindedPlot(lumi = lumi2018_EE, hname_SI = 'hEESR_nEE', hname_bkg = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2018, inputdir = input_ee, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1])
        makeBlindedPlot(lumi = lumi2018_EE, hname_SI = 'hEESR_mass', hname_bkg = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2018, inputdir = input_ee, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1,2])
        makeBlindedPlot(lumi = lumi2018_EE, hname_SI = 'hEESRI_dPhi_scan', hname_bkg = 'hEEBCRI_dPhi_scan', ylog = True, treeDATA = treeDATA_EG2018, inputdir = input_ee, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'EE', DATAlabel = '', extralabel = 'N_{ee} = 1, Off-Z region', outpath = www, cutBins = [1,2])
        makeBlindedPlot(lumi = lumi2018_MM, hname_SI = 'hMMSR_nMM', hname_bkg = 'hMMBCR_nMM', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = input_mm, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1])
        makeBlindedPlot(lumi = lumi2018_MM, hname_SI = 'hMMSR_mass', hname_bkg = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = input_mm, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = 'Baseline selection, Sig(d_{0}) > 1', outpath = www, cutBins = [1,2])
        makeBlindedPlot(lumi = lumi2018_MM, hname_SI = 'hMMSRI_dPhi_scan', hname_bkg = 'hMMBCRI_dPhi_scan', ylog = True, treeDATA = treeDATA_Mu2018, inputdir = input_mm, treeSI = treeSI_2018, xlabel = '', outtag = '2018', ymax = 1e12, LLlabel = 'MM', DATAlabel = '', extralabel = 'N_{#mu#mu} = 1, Off-Z region', outpath = www, cutBins = [1,2])



