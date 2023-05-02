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
from include.Utils import XSECS, buildSummaryPlot

################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

def printSignals(signalNames):
    '''
    Function to print the legend of the tables
    '''
    print("_"*80)
    print("Signal Legend")
    print("_"*80)
    for i in range(len(signalNames)):
        print("["+str(i)+"] {:<12}".format(signalNames[i]))


def printTable(title, regions, BKG, signal):
    '''
    Function to print the tables with the yields in LaTeX format
    '''
    nsignals = len(signal[0])
    print("    \\begin{center}")
    print("        {\\footnotesize")
    print("        \\begin{tabular}{ |c|c|"+"c|"*nsignals+" }")
    print("        \\hline")
    title_line = "        \\multicolumn{%d}{|c|}{"%(nsignals+2)+title+"} \\\\"
    print(title_line)
    print("        \\hline")
    header = "        {:<6} & {:<7} "+"& {:<7} "*nsignals+" \\\\"
    print(header.format('Region', 'BKG', *["["+str(i)+"]" for i in range(len(signal[0]))]))
    print("        \\hline")
    format_line = "        {:<6} & {:<7.1f} "+"& {:<7.1f} "*len(signal[0])+" \\\\"
    for n in range(len(regions)):
        print(format_line.format(regions[n], BKG[n], *signal[n]))
    print("        \\hline")
    print("        \\end{tabular}")
    print("        }")
    print("    \\end{center}")



if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--inputMuon', action='store', type=str, dest='inputMuons', default='', help='Target directory')
    parser.add_option('-e', '--inputElectron', action='store', type=str, dest='inputElectrons', default='', help='Target directory')
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
    Signals.append('HSS_500_50_100')
    Signals.append('HSS_1000_250_100')
    #Signals.append('HSS_600_50_100')
    #Signals.append('HSS_1000_350_100')
    Signals.append('RPV_350_148_100')
    Signals.append('RPV_1500_494_100')
    
    Signals_2016preVFP = [i + '_2016APV' for i in Signals]
    Signals_2016postVFP = [i + '_2016' for i in Signals]
    Signals2016 = Signals_2016preVFP + Signals_2016postVFP
    Signals2017 = [i + '_2017' for i in Signals]
    Signals2018 = [i + '_2018' for i in Signals]

    ############# Luminosity definition
    lumiB = 5.79
    lumiC = 2.57
    lumiD = 4.25
    lumiE = 4.01
    lumiF = 2.53 # total 3.10
    lumiF_noHIMP = 0.57
    lumiG = 7.54
    lumiH = 8.61
    lumi =  lumiB + lumiC + lumiD + lumiE + lumiF + lumiG + lumiH# luminosity

    lumi_2016 = 35.9
    lumi_2016_GH = 16.2
    lumi_2017 = 41.5
    lumi_2018_EE = 54.5
    lumi_2018_MM = 59.8

    ##### Cross sections
    #xsecs = np.array([0.107e3, 4.938e3, 2.588e3, 0.67, 10e3])

    ############ Define regions
    regions_mu = ["IaA", "IaB", "IaC", "IaD", "IbA", "IbB", "IbC", "IbD", "II"]
    regions_ee = ["IaA", "IaB", "IaC", "IbA", "IbB", "IbC", "II"]

    ############ Define output
    www = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/SignalYieldsAndTables/No-Parrot/'

    ############ Dielectron plots
    if opts.inputElectrons:
        treeSI_2016_GH = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2016UL_Fall22.dat', Signals_2016postVFP, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2017    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2017UL_Fall22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2018    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2018UL_Fall22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )
        treeSI_Full    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_Full_Spring23.dat', Signals_2016postVFP + Signals2017 + Signals2018, 'SI'), name = 'SI', isdata = 0 )

        treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_EGFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016 + DoubleEG2017 + EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

        #unblinded = ["IaA", "IaB", "IaC", "IbA", "IbB", "IbC", "II"]
        unblinded = []
        BKG_EG_2016, Signal_EG_2016, titles = buildSummaryPlot("EG_2016", treeDATA_EG2016, treeSI_2016_GH, opts.inputElectrons, regions_ee, lumi_2016_GH, LLlabel='EE', sys = 0.1, unblinded = unblinded, outpath = www)
        BKG_EG_2017, Signal_EG_2017, _ = buildSummaryPlot("EG_2017", treeDATA_EG2017, treeSI_2017, opts.inputElectrons, regions_ee, lumi_2017, LLlabel='EE', sys = 0.1, unblinded = unblinded, outpath = www) 
        BKG_EG_2018, Signal_EG_2018, _ = buildSummaryPlot("EG_2018", treeDATA_EG2018, treeSI_2018, opts.inputElectrons, regions_ee, lumi_2018_EE, LLlabel='EE', sys = 0.1, unblinded = unblinded, outpath = www)
        BKG_EG_Full, Signal_EG_Full, _ = buildSummaryPlot("EG_Full", treeDATA_EGFull, treeSI_Full, opts.inputElectrons, regions_ee, 112, LLlabel='EE', sys = 0.1, unblinded = unblinded, outpath = www)


    if opts.inputMuons:
        treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2016UL_Fall22.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2017UL_Fall22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2018UL_Fall22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )
        treeSI_Full    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_Full_Spring23.dat', Signals2016 + Signals2017 + Signals2018, 'SI'), name = 'SI', isdata = 0 )

        treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_MuFull = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016 + DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

        unblinded = ["IaA", "IaB", "IaC", "IaD", "IbA", "IbB", "IbC", "IbD", "II"]
        unblinded = []
        sys_value = 0.0
        BKG_Mu_2016, Signal_Mu_2016, titles = buildSummaryPlot("Mu_2016", treeDATA_Mu2016, treeSI_2016, opts.inputMuons, regions_mu, lumi_2016, LLlabel='MM', sys = 0.1, unblinded = unblinded, outpath = www)
        BKG_Mu_2018, Signal_Mu_2018, _ = buildSummaryPlot("Mu_2018", treeDATA_Mu2018, treeSI_2018, opts.inputMuons, regions_mu, lumi_2018_MM, LLlabel='MM', sys = 0.1, unblinded = unblinded, outpath = www)
        BKG_Mu_Full, Signal_Mu_Full, _ = buildSummaryPlot("Mu_Full", treeDATA_MuFull, treeSI_Full, opts.inputMuons, regions_mu, 96, LLlabel='MM', sys = 0.1, unblinded = unblinded, outpath = www)

