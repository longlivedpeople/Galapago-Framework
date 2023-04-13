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


def buildPlot(title, treeBKG, treeSI, inputdir, regions, luminosity, LLlabel = 'MM'):
    '''
    Function to make the plots
    Output (neccesary arrays to print the tables with 'printTable' and 'printSignals')
    '''
    # Arrays to store yields for table
    BKG_yields = []
    signal_yields = []
    titles = []

    nbins = len(regions)
    
    hbkg = r.TH1F("hbkg","hbkg;;N events", nbins, 0, nbins)
    hbkg.Sumw2()

    ## Fill BKG hist
    for n,reg in enumerate(regions):
        hname = 'h'+LLlabel+'_BCR'+reg
        h_ = treeBKG.getLoopTH1F(inputdir, hname)
        hbkg.GetXaxis().SetBinLabel(n+1,reg)
        hbkg.SetBinContent(n+1,h_.GetBinContent(1))
        hbkg.SetBinError(n+1,h_.GetBinError(1))
        BKG_yields.append(h_.GetBinContent(1))
        
    ### Set background histos style
    hbkg.SetFillColorAlpha(r.kCyan-6, 0.8) 
    hbkg.SetLineColor(r.kCyan-2) 
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)
    
    # Get signal titles
    hname = 'h'+LLlabel+'_SR'+regions[0]
    print(hname)
    h_stack = treeSI.getLoopStack(inputdir, hname)
    s_histos = []
    for i,h_ in enumerate(h_stack):
        s_histos.append(r.TH1F(h_.GetTitle(),h_.GetTitle(),nbins,0,nbins)) # declare histograms
        signal_yields.append([]) # declare array for yields

    ## Fill S_histos
    for n,reg in enumerate(regions):
        hname = 'h'+LLlabel+'_SR'+reg
        h_stack = treeSI.getLoopStack(inputdir, hname)
        for i,h_ in enumerate(h_stack):
            s_histos[i].GetXaxis().SetBinLabel(n+1,reg)
            print(h_.GetTitle(), XSECS[h_.GetTitle()])
            h_.Scale(XSECS[h_.GetTitle()])
            content = h_.GetBinContent(1)
            s_histos[i].SetBinContent(n+1,content)
            s_histos[i].SetBinError(n+1,h_.GetBinError(1))
            signal_yields[i].append(content)
            titles.append(h_.GetTitle())

    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])

    hbkg.SetMaximum(1e3*maxVal)
    
    ### -> Canvas object
    plot = Canvas.Canvas(title, 'png,pdf', 0.15, 0.65, 0.45, 0.89, 1)

    plot.addHisto(hbkg, 'HIST', 'Background (predicted)', 'f', '', 1, 0)
 
    colors = [r.kRed, r.kOrange, r.kGreen+2, r.kBlue, r.kMagenta]
    for i,_h in enumerate(s_histos):
        _h.SetLineWidth(2)
        masses = eval(_h.GetTitle()[3:])
        if 'HSS' in _h.GetTitle():
            legend = 'H#rightarrowSS (%d GeV, %d GeV, %d mm)'%(masses[0], masses[1], masses[2])
        else:
            legend = 'RPV (%d GeV, %d GeV, %d mm)'%(masses[0], masses[1], masses[2])
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', colors[i], 1, i+1) # Signal

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.7, 0.86, 'e^{+}e^{-} channel', font = 42, size = 0.035)
    if LLlabel == 'MM':
        plot.addLatex(0.7, 0.86, '#mu^{+}#mu^{-} channel', font = 42, size = 0.035)

    ### Save it
    outdir = 'Plots_yieldRegions'
    plot.save(1, 1, True, luminosity, '', outputDir = outdir, xlog = False, maxYnumbers = 4, is2d = True, inProgress = True)

    # Return info for table
    return np.asarray(BKG_yields), np.transpose(np.asarray(signal_yields)), titles

################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

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

    ############ Dielectron plots
    if opts.inputElectrons:
        treeSI_2016_GH = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2016UL_Fall22.dat', Signals_2016postVFP, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2017UL_Fall22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2018UL_Fall22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )

        treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

        BKG_EG_2016, Signal_EG_2016, titles = buildPlot("EG_2016", treeDATA_EG2016_GH, treeSI_2016_GH, opts.inputElectrons, regions_ee, lumi_2016_GH, LLlabel='EE')
        BKG_EG_2017, Signal_EG_2017, _ = buildPlot("EG_2017", treeDATA_EG2017, treeSI_2017, opts.inputElectrons, regions_ee, lumi_2017, LLlabel='EE') 
        BKG_EG_2018, Signal_EG_2018, _ = buildPlot("EG_2018", treeDATA_EG2018, treeSI_2018, opts.inputElectrons, regions_ee, lumi_2018_EE, LLlabel='EE')

        original_stdout = sys.stdout # Save a reference to the original standard output
        with open('Plots_yieldRegions/tables_EE.txt', 'w') as f:
            sys.stdout = f
            printSignals(titles)
            printTable("Electron channel 2016", regions_ee, BKG_EG_2016, Signal_EG_2016)
            printTable("Electron channel 2017", regions_ee, BKG_EG_2017, Signal_EG_2017)
            printTable("Electron channel 2018", regions_ee, BKG_EG_2018, Signal_EG_2018)
            sys.stdout = original_stdout

    if opts.inputMuons:
        treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2016UL_Fall22.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2017UL_Fall22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
        treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2018UL_Fall22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )

        treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
        treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

        BKG_Mu_2016, Signal_Mu_2016, titles = buildPlot("Mu_2016", treeDATA_Mu2016, treeSI_2016, opts.inputMuons, regions_mu, lumi_2016, LLlabel='MM')
        BKG_Mu_2018, Signal_Mu_2018, _ = buildPlot("Mu_2018", treeDATA_Mu2018, treeSI_2018, opts.inputMuons, regions_mu, lumi_2018_MM, LLlabel='MM')

        original_stdout = sys.stdout # Save a reference to the original standard output
        with open('Plots_yieldRegions/tables_MuMu.txt', 'w') as f:
            sys.stdout = f
            printSignals(titles)
            printTable("Muon channel 2016", regions_mu, BKG_Mu_2016, Signal_Mu_2016)
            printTable("Muon channel 2018", regions_mu, BKG_Mu_2018, Signal_Mu_2018)
            sys.stdout = original_stdout

    """
    ############# Galapago Tree definitions
    treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2016UL_Fall22.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2016_GH = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2016UL_Fall22.dat', Signals_2016postVFP, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2017UL_Fall22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_2018UL_Fall22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2016_GH = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_GH, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    
    ########### Print plot 
    BKG_Mu_2016, Signal_Mu_2016, titles = buildPlot("Mu_2016", treeDATA_Mu2016, treeSI_2016, opts.inputdir, regions_mu, lumi_preVFP+lumi_postVFP, LLlabel='MM')
    BKG_Mu_2018, Signal_Mu_2018, _ = buildPlot("Mu_2018", treeDATA_Mu2018, treeSI_2018, opts.inputdir, regions_mu, lumi_2018, LLlabel='MM')
    BKG_EG_2016, Signal_EG_2016, _ = buildPlot("EG_2016", treeDATA_EG2016_GH, treeSI_2016_GH, opts.inputdir, regions_ee, lumi_2016_GH, LLlabel='EE')
    BKG_EG_2017, Signal_EG_2017, _ = buildPlot("EG_2017", treeDATA_EG2017, treeSI_2017, opts.inputdir, regions_ee, lumi_2017, LLlabel='EE') 
    BKG_EG_2018, Signal_EG_2018, _ = buildPlot("EG_2018", treeDATA_EG2018, treeSI_2018, opts.inputdir, regions_ee, lumi_2018, LLlabel='EE')
    
    ########### Write table
    original_stdout = sys.stdout # Save a reference to the original standard output
    with open('Plots_yieldRegions/tables_dummy.txt', 'w') as f:
        sys.stdout = f
        printSignals(titles)
        printTable("Muon channel 2016", regions_mu, BKG_Mu_2016, Signal_Mu_2016)
        printTable("Muon channel 2018", regions_mu, BKG_Mu_2018, Signal_Mu_2018)
        printTable("Electron channel 2016", regions_ee, BKG_EG_2016, Signal_EG_2016)
        printTable("Electron channel 2017", regions_ee, BKG_EG_2017, Signal_EG_2017)
        printTable("Electron channel 2018", regions_ee, BKG_EG_2018, Signal_EG_2018)
        sys.stdout = original_stdout
    """
