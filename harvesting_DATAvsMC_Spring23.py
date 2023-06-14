import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os, json
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
from include.galapagoStyle import sigpalette, gcolors, dcolors
from include.Utils import makeSystematicsHist


def makeDataMCPlot(lumi, hname_DATA, hname_MC, ylog, treeDATA, treeMC, inputdir, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '', yshift = 100.0, LLlabel = '', leftlabel = '', xlog = False, outdir = '', sys_errors = None, drawZero = False):

    
    luminosity = lumi

    hMC      = copy.deepcopy(treeMC.getLoopStack(inputdir, hname_MC)) # Stacked with MC
    hMCtotal = copy.deepcopy(treeMC.getLoopTH1F(inputdir, hname_MC))
    hDATA    = copy.deepcopy(treeDATA.getLoopTH1F(inputdir, hname_DATA))

    hBKG     = copy.deepcopy(treeMC.getLoopStack(inputdir, hname_MC)) # Dummy stacked to be filled with contributions
    for _h in hBKG.GetHists():
        hBKG.RecursiveRemove(_h)

    if treeEXTRA and hname_extra:
        hEXTRA = copy.deepcopy(treeDATA.getLoopTH1F(inputdir, hname_extra))
        hEXTRA.SetLineColor(r.kBlack)
        #hEXTRA.SetFillColor(r.kMagenta+1)
        hEXTRA.SetFillColor(r.kCyan-1)
        hEXTRA.SetTitle(label_extra)
        r.SetOwnership(hEXTRA, 0)

    ### Combine BKG + EXTRA contributions:
    nBKG = len(hMC.GetHists()) #+ 1
    #hBKG = hMC

    #hMCtotal.Scale(7.10/59.83)
    if treeEXTRA and hname_extra:
        hBKG.Add(hEXTRA)
        hMCtotal.Add(hEXTRA)
        nBKG += 1

    for _h in hMC.GetHists():
        #_h.Scale(7.10/59.83)
        hBKG.Add(copy.deepcopy(_h))

    ### Systematics
    """
    hsys = None
    if sys_errors is not None: hsys = makeSystematicsHist(sys_errors, hMCtotal)
    """
    if sys_errors is not None: 
        error_values = 1 * np.array(sys_errors)
        sys = np.linalg.norm(error_values)

    ### Uncertainty bkg
    hunc = hMCtotal.Clone("uncertainty")
    hunc.Reset()
    hunc.Sumw2()
    for n in range(1, hMCtotal.GetNbinsX() + 1):
        hunc.SetBinContent(n, hMCtotal.GetBinContent(n))
        hunc.SetBinError(n, math.sqrt(sys*hMCtotal.GetBinContent(n)*sys*hMCtotal.GetBinContent(n) + hMCtotal.GetBinError(n)*hMCtotal.GetBinError(n)))
        if drawZero and  hunc.GetBinContent(n) < 1.:
            hunc.SetBinError(n, 1.8)
            if ylog: hunc.SetBinContent(n, 0.1)


    ### Tune the histograms
    hMCtotal.SetMarkerStyle(20) # Auxiliar to save the ratio correctly 
    hDATA.SetMarkerStyle(20)
    hDATA.SetMarkerSize(0.8)
    hDATA.SetLineWidth(2)
    hunc.SetFillStyle(3244)
    hunc.SetFillColor(r.kCyan+3)
    hunc.SetLineColor(r.kCyan+3)
    hunc.SetMarkerSize(0)
    
    ### Get maximum
    maxValMC = hMCtotal.GetMaximum()
    maxValDATA = hDATA.GetMaximum()
    maxVal = max([maxValMC, maxValDATA])
    
        ### Set Maximum
    if not ylog:
        hBKG.SetMaximum(1.3*maxVal)
        hBKG.SetMinimum(0.0)
        hMCtotal.SetMaximum(1.3*maxVal)
        hMCtotal.SetMinimum(0.0)
        if treeDATA:
            if not yshift:
                hDATA.SetMaximum(1.3*maxVal)
            else:
                hDATA.SetMaximum(yshift*maxVal)
            hDATA.SetMinimum(0.0)
    else:
        if not yshift:
            hBKG.SetMaximum(10.0*maxVal)
        else:
            hBKG.SetMaximum(yshift*maxVal)
        hBKG.SetMinimum(0.1)
        if not yshift:
            hMCtotal.SetMaximum(10.0*maxVal)
        else:
            hMCtotal.SetMaximum(yshift*maxVal)
        hMCtotal.SetMinimum(0.1)
        if treeDATA:
            if not yshift:
                hDATA.SetMaximum(10.0*maxVal)
            else:
                hDATA.SetMaximum(yshift*maxVal)
            hDATA.SetMinimum(0.1)
    
    ### SetOwnership
    r.SetOwnership(hDATA, 0)
    r.SetOwnership(hMC, 0)
    r.SetOwnership(hBKG, 0)
    r.SetOwnership(hMCtotal, 0)
    
    ### -> Canvas object
    plot = Canvas.Canvas(hname_DATA, 'png,pdf', 0.52, 0.55, 0.77, 0.85, 1)
    
    ### Add background:
    plot.addStack(hBKG, 'HIST', 1, 0) # Background
    plot.addHisto(hunc, 'E2, SAME', 'Background uncertainty', 'f', '', 1, nBKG)
    plot.addHisto(hDATA, 'P, SAME', 'Data', 'p', r.kBlack, 1, nBKG+1)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.79, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.79, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if leftlabel:
        plot.addLatex(0.17, 0.74, leftlabel, font = 42)


    ### Save it
    if not outdir:
        outfolder = os.path.dirname(os.path.abspath(__main__.__file__)) + '/PlotsDATAMC_' + outtag + '/'
    else:
        outfolder = outdir + '/PlotsDATAMC_' + outtag + '/'
    plot.saveRatio2(1, 1, ylog, luminosity, hDATA, hMCtotal, r_ymin = 0.0, r_ymax = 2.0, label="Data/BKG", outputDir = outfolder, xlog = xlog, isPrivate = True, sys = sys)
    


def makeMCPlot(lumi, hname_MC, ylog, treeMC, inputdir, treeSI = False, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '', yshift = 100.0, LLlabel = '', leftlabel = '', xlog = False):

    
    luminosity = lumi

    hMC      = copy.deepcopy(treeMC.getLoopStack(inputdir, hname_MC)) # Stacked with MC
    hMCtotal = copy.deepcopy(treeMC.getLoopTH1F(inputdir, hname_MC))

    hBKG     = copy.deepcopy(treeMC.getLoopStack(inputdir, hname_MC)) # Dummy stacked to be filled with contributions
    for _h in hBKG.GetHists():
        hBKG.RecursiveRemove(_h)

    if treeEXTRA and hname_extra:
        hEXTRA = copy.deepcopy(treeEXTRA.getLoopTH1F(inputdir, hname_extra))
        hEXTRA.SetLineColor(r.kBlack)
        hEXTRA.SetFillColor(r.kMagenta+1)
        hEXTRA.SetTitle(label_extra)
        r.SetOwnership(hEXTRA, 0)

    ### Combine BKG + EXTRA contributions:
    nBKG = len(hMC.GetHists()) #+ 1

    if treeEXTRA and hname_extra:
        hBKG.Add(hEXTRA)
        hMCtotal.Add(hEXTRA)
        nBKG += 1

    for _h in hMC.GetHists():
        hBKG.Add(copy.deepcopy(_h))

    ### Tune the histograms
    hMCtotal.SetMarkerStyle(20) # Auxiliar to save the ratio correctly 
    
    ### Signal histograms
    s_histos = []
    #xsecs = [4.938e3, 2.588e3, 0.107e3, 10e3, 0.67] # HARDCODED THIS NEEDS TO CHANGE
    xsecs = [0.107e3, 4.938e3, 2.588e3, 0.67, 10e3] # HARDCODED THIS NEEDS TO CHANGE
    if treeSI:
        hSIS = treeSI.getLoopStack(inputdir, hname_MC)
        for _i, _h in enumerate(hSIS.GetHists()):
            s_histos.append(copy.deepcopy(_h))
            s_histos[-1].Scale(xsecs[_i])


    ### Get maximum
    maxValMC = hMCtotal.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValMC, maxValSI])
    
        ### Set Maximum
    if not ylog:
        hBKG.SetMaximum(1.3*maxVal)
        hBKG.SetMinimum(0.0)
        hMCtotal.SetMaximum(1.3*maxVal)
        hMCtotal.SetMinimum(0.0)
    else:
        if not yshift:
            hBKG.SetMaximum(10.0*maxVal)
        else:
            hBKG.SetMaximum(yshift*maxVal)
        hBKG.SetMinimum(0.1)
        if not yshift:
            hMCtotal.SetMaximum(10.0*maxVal)
        else:
            hMCtotal.SetMaximum(yshift*maxVal)
        hMCtotal.SetMinimum(0.1)
    
    ### SetOwnership
    r.SetOwnership(hMC, 0)
    r.SetOwnership(hBKG, 0)
    r.SetOwnership(hMCtotal, 0)
    
    ### -> Canvas object
    plot = Canvas.Canvas(hname_MC, 'png', 0.35, 0.6, 0.6, 0.85, 1)
    
    ### Add background:
    plot.addStack(hBKG, 'HIST', 1, 0) # Background
    color_order = ['red', 'orange', 'yellow', 'green', 'magenta']
    for i,_h in enumerate(s_histos):
        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        if 'HSS' in _h.GetTitle():
            legend = 'H#rightarrowSS (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        else:
            legend = 'RPV (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', gcolors[color_order[i]], 1, i+nBKG) # Signal    

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.76, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.76, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if leftlabel:
        plot.addLatex(0.17, 0.7, leftlabel, font = 42)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/PlotsMC_' + outtag + '/'
    plot.save(1, 0, ylog, luminosity, '', outputDir = outdir, xlog = xlog, maxYnumbers = 4)



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


    ############# Dat file
    filename = 'dat/Samples_cern_UltraLegacy_thesisColors.dat'
    filename = 'dat/Samples_cern_UltraLegacy_thesisColors_NLO.dat'
    filename = 'dat/Samples_cern_UltraLegacy_Spring23.dat'

    ############# Systematics
    sys = json.loads(open('systematics/systematics.json').read())

    ############# EG data definition
    DoubleEG2016_HIPM = []
    DoubleEG2016_noHIPM = []
    DoubleEG2016_HIPM.append('DoubleEG_Run2016B_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016C_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016D_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016E_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016F_HIPM')
    #DoubleEG2016_noHIPM.append('DoubleEG_Run2016F_noHIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016H_noHIPM')
    DoubleEG2016 = DoubleEG2016_noHIPM + DoubleEG2016_HIPM
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


    ############# Background simulation definition
    Backgrounds_preVFP = []
    Backgrounds_preVFP.append('DYJetsToLL_M-50_preVFP')
    Backgrounds_preVFP.append('DYJetsToLL_M-10to50_preVFP')
    Backgrounds_preVFP.append('TTTo2L2Nu_preVFP')
    #Backgrounds_preVFP.append('WJetsToLNu_preVFP')
    Backgrounds_preVFP.append('WW_preVFP')
    Backgrounds_preVFP.append('WZ_preVFP')
    Backgrounds_preVFP.append('ZZ_preVFP')
    Backgrounds_postVFP = []
    Backgrounds_postVFP.append('DYJetsToLL_M-50_postVFP')
    Backgrounds_postVFP.append('DYJetsToLL_M-10to50_postVFP')
    Backgrounds_postVFP.append('TTTo2L2Nu_postVFP')
    #Backgrounds_postVFP.append('WJetsToLNu_postVFP')
    Backgrounds_postVFP.append('WW_postVFP')
    Backgrounds_postVFP.append('WZ_postVFP')
    Backgrounds_postVFP.append('ZZ_postVFP')
    Backgrounds_2017 = []
    Backgrounds_2017.append('DYJetsToLL_M-50_2017')
    Backgrounds_2017.append('DYJetsToLL_M-10to50_2017')
    Backgrounds_2017.append('TTTo2L2Nu_2017')
    #Backgrounds_2017.append('WJetsToLNu_2017')
    Backgrounds_2017.append('WW_2017')
    Backgrounds_2017.append('WZ_2017')
    Backgrounds_2017.append('ZZ_2017')
    Backgrounds_2018 = []
    Backgrounds_2018.append('DYJetsToLL_M-50_2018')
    Backgrounds_2018.append('DYJetsToLL_M-10to50_2018')
    Backgrounds_2018.append('TTTo2L2Nu_2018')
    #Backgrounds_2018.append('WJetsToLNu_2018')
    Backgrounds_2018.append('WW_2018')
    Backgrounds_2018.append('WZ_2018')
    Backgrounds_2018.append('ZZ_2018')

    ############# Signal definition
    Signals_HSS = []
    Signals_HSS.append('HSS_300_50_100')
    Signals_HSS.append('HSS_500_50_100')
    Signals_HSS.append('HSS_1000_250_100')
    Signals_HSS.append('RPV_350_148_100')
    Signals_HSS.append('RPV_1500_494_100')

    Signals_HSS_2016preVFP = [i + '_2016APV' for i in Signals_HSS]
    Signals_HSS_2016postVFP = [i + '_2016' for i in Signals_HSS]
    Signals_HSS_2017 = [i + '_2017' for i in Signals_HSS]
    Signals_HSS_2018 = [i + '_2018' for i in Signals_HSS]

    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2016_HIPM = 19.7 # fb-1
    lumi2016_noHIPM = 16.2 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 54.5 # fb-1
    lumi2018 = 59.8 # fb-1

    #treeMC_2016     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_postVFP, 'MC'), name = 'MC', isdata = 0, close = True)


    treeHSS_2016    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2016UL_Fall22.dat', Signals_HSS_2016preVFP + Signals_HSS_2016postVFP, 'SI'), name = 'SI', isdata = 0, close = True)
    treeHSS_2017    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2017UL_Fall22.dat', Signals_HSS_2017, 'SI'), name = 'SI', isdata = 0, close = True)
    treeHSS_2018    = Sample.Tree( fileName = helper.selectSamples(WORKPATH + '/dat/CombSignal_2018UL_Fall22.dat', Signals_HSS_2018, 'SI'), name = 'SI', isdata = 0, close = True)



    ############ Plot everything
    if False:
        ### 2016
        makeMCPlot(lumi = lumi2016, hname_MC = 'hEESel_dPhi', ylog = True, treeMC = treeMC_2016, inputdir = opts.input, treeSI = treeHSS_2016, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2016', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2016, hname_MC = 'hMMSel_dPhi', ylog = True, treeMC = treeMC_2016, inputdir = opts.input, treeSI = treeHSS_2016, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2016', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2016, hname_MC = 'hEESel_mass', ylog = True, treeMC = treeMC_2016, inputdir = opts.input, treeSI = treeHSS_2016, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2016', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2016, hname_MC = 'hMMSel_mass', ylog = True, treeMC = treeMC_2016, inputdir = opts.input, treeSI = treeHSS_2016, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2016', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 

        ## 2017
        makeMCPlot(lumi = lumi2017, hname_MC = 'hEESel_dPhi', ylog = True, treeMC = treeMC_2017, inputdir = opts.input, treeSI = treeHSS_2017, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2017', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2017, hname_MC = 'hEESel_mass', ylog = True, treeMC = treeMC_2017, inputdir = opts.input, treeSI = treeHSS_2017, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2017', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 

        ## 2018
        makeMCPlot(lumi = lumi2018, hname_MC = 'hEESel_dPhi', ylog = True, treeMC = treeMC_2018, inputdir = opts.input, treeSI = treeHSS_2018, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2018', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2018, hname_MC = 'hMMSel_dPhi', ylog = True, treeMC = treeMC_2018, inputdir = opts.input, treeSI = treeHSS_2018, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2018', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2018, hname_MC = 'hEESel_mass', ylog = True, treeMC = treeMC_2018, inputdir = opts.input, treeSI = treeHSS_2018, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2018', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 
        makeMCPlot(lumi = lumi2018, hname_MC = 'hMMSel_mass', ylog = True, treeMC = treeMC_2018, inputdir = opts.input, treeSI = treeHSS_2018, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '2018', yshift = 1e4, LLlabel = '', leftlabel = '', xlog = False) 

    www = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/Data-vs-MonteCarlo/Data-vs-MonteCarlo-muons-final-private'

    if True:

        ## 2016
        treeMC_2016     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_preVFP + Backgrounds_postVFP, 'MC'), name = 'MC', isdata = 0, close = True)

        """
        treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1, close = True )
        # electrons
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_nEE', hname_MC = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_mass', hname_MC = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_MET', hname_MC = 'hEEBCR_MET', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_leadingEt', hname_MC = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_subleadingEt', hname_MC = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_eta', hname_MC = 'hEEBCR_eta', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])

        treeMC_postVFP     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_postVFP, 'MC'), name = 'MC', isdata = 0, close = True)
        treeDATA_EG2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_noHIPM, 'DATA'), name = 'DATA', isdata = 1, close = True )
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_nEE', hname_MC = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_mass', hname_MC = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_MET', hname_MC = 'hEEBCR_MET', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 10000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_leadingEt', hname_MC = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_subleadingEt', hname_MC = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_noHIPM, hname_DATA = 'hEEBCR_eta', hname_MC = 'hEEBCR_eta', ylog = True, treeDATA = treeDATA_EG2016_noHIPM, treeMC = treeMC_postVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_noHIPM', yshift = 10000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])

        treeDATA_EG2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_HIPM, 'DATA'), name = 'DATA', isdata = 1, close = True )
        treeMC_preVFP     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_preVFP, 'MC'), name = 'MC', isdata = 0, close = True)
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_nEE', hname_MC = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_mass', hname_MC = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_MET', hname_MC = 'hEEBCR_MET', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 10000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_leadingEt', hname_MC = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_subleadingEt', hname_MC = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        makeDataMCPlot(lumi = lumi2016_HIPM, hname_DATA = 'hEEBCR_eta', hname_MC = 'hEEBCR_eta', ylog = True, treeDATA = treeDATA_EG2016_HIPM, treeMC = treeMC_preVFP, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016_HIPM', yshift = 100000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2016'])
        """

        # muons
        treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1, close = True)
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_nMM', hname_MC = 'hMMBCR_nMM', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_mass', hname_MC = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_MET', hname_MC = 'hMMBCR_MET', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_mass_Z', hname_MC = 'hMMBCR_mass_Z', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 10000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_leadingPt', hname_MC = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_subleadingPt', hname_MC = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])
        makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_eta', hname_MC = 'hMMBCR_eta', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 100000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2016'])

    if False:

        treeMC_2017     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_2017, 'MC'), name = 'MC', isdata = 0, close = True)
        treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1, close = True)


        ## 2017
        makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEEBCR_mass', hname_MC = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2017, treeMC = treeMC_2017, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2017, xlabel = '', outtag = '2017', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2017'])
        makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2017, treeMC = treeMC_2017, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2017, xlabel = '', outtag = '2017', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2017'])
        makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEEBCR_leadingEt', hname_MC = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2017, treeMC = treeMC_2017, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2017, xlabel = '', outtag = '2017', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2017'])
        makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEEBCR_subleadingEt', hname_MC = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2017, treeMC = treeMC_2017, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2017, xlabel = '', outtag = '2017', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2017'])
        makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEEBCR_eta', hname_MC = 'hEEBCR_eta', ylog = True, treeDATA = treeDATA_EG2017, treeMC = treeMC_2017, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2017, xlabel = '', outtag = '2017', yshift = 10000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2017'])


    if True:

        treeMC_2018     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_2018, 'MC'), name = 'MC', isdata = 0, close = True)
        """
        treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1, close = True)

        ## 2018
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_nEE', hname_MC = 'hEEBCR_nEE', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_mass', hname_MC = 'hEEBCR_mass', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 100000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_leadingEt', hname_MC = 'hEEBCR_leadingEt', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_subleadingEt', hname_MC = 'hEEBCR_subleadingEt', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_eta', hname_MC = 'hEEBCR_eta', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 100000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'EE', xlog = False, outdir = www, sys_errors = sys['EE_2018'])
        """

        treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1, close = True)
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMMBCR_mass', hname_MC = 'hMMBCR_mass', ylog = True, treeDATA = treeDATA_Mu2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMMBCR_mass_Z', hname_MC = 'hMMBCR_mass_Z', ylog = True, treeDATA = treeDATA_Mu2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2018, xlabel = '', outtag = '2018', yshift = 10000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMMBCR_leadingPt', hname_MC = 'hMMBCR_leadingPt', ylog = True, treeDATA = treeDATA_Mu2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMMBCR_subleadingPt', hname_MC = 'hMMBCR_subleadingPt', ylog = True, treeDATA = treeDATA_Mu2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2018, xlabel = '', outtag = '2018', yshift = 1000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2018'])
        makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMMBCR_eta', hname_MC = 'hMMBCR_eta', ylog = True, treeDATA = treeDATA_Mu2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = '', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2018, xlabel = '', outtag = '2018', yshift = 100000.0, leftlabel = 'Control Region: |#Delta#Phi| > #pi/2', LLlabel = 'MM', xlog = False, outdir = www, sys_errors = sys['MM_2018'])
