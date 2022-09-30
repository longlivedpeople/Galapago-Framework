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





def makeDataMCPlot(lumi, hname_DATA, hname_MC, ylog, treeDATA, treeMC, inputdir, hname_extra = '', label_extra = '', treeEXTRA = False, xlabel = '', outtag = '', yshift = 100.0, LLlabel = '', leftlabel = '', xlog = False):

    
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
        hEXTRA.SetFillColor(r.kMagenta+1)
        hEXTRA.SetTitle(label_extra)
        r.SetOwnership(hEXTRA, 0)

    ### Combine BKG + EXTRA contributions:
    nBKG = len(hMC.GetHists()) #+ 1
    #hBKG = hMC

    if treeEXTRA and hname_extra:
        hBKG.Add(hEXTRA)
        hMCtotal.Add(hEXTRA)
        nBKG += 1

    for _h in hMC.GetHists():
        hBKG.Add(copy.deepcopy(_h))

    ### Tune the histograms
    hMCtotal.SetMarkerStyle(20) # Auxiliar to save the ratio correctly 
    hDATA.SetMarkerStyle(20)
    hDATA.SetMarkerSize(0.8)
    
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
    plot = Canvas.Canvas(hname_DATA, 'png', 0.6, 0.6, 0.85, 0.84, 1)
    
    ### Add background:
    plot.addStack(hBKG, 'HIST', 1, 0) # Background
    plot.addHisto(hDATA, 'P, SAME', 'Data', 'p', r.kBlack, 1, nBKG)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if leftlabel:
        plot.addLatex(0.17, 0.76, leftlabel, font = 42)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/PlotsDATAMC_' + outtag + '/'
    plot.saveRatio(1, 1, ylog, luminosity, hDATA, hMCtotal, r_ymin = 0.0, r_ymax = 2.0, label="Data/BKG", outputDir = outdir, xlog = xlog)
    

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


    ############# Background simulation definition
    Backgrounds_preVFP = []
    Backgrounds_preVFP.append('DYJetsToLL_M-50_preVFP')
    Backgrounds_preVFP.append('TTTo2L2Nu_preVFP')
    Backgrounds_preVFP.append('WW_preVFP')
    Backgrounds_preVFP.append('WZ_preVFP')
    Backgrounds_postVFP = []
    Backgrounds_postVFP.append('DYJetsToLL_M-50_postVFP')
    Backgrounds_postVFP.append('TTTo2L2Nu_postVFP')
    Backgrounds_postVFP.append('WW_postVFP')
    Backgrounds_postVFP.append('WZ_postVFP')
    Backgrounds_2017 = []
    Backgrounds_2017.append('DYJetsToLL_M-50_2017')
    Backgrounds_2017.append('TTTo2L2Nu_2017')
    Backgrounds_2017.append('WW_2017')
    Backgrounds_2017.append('WZ_2017')
    Backgrounds_2017.append('ZZ_2017')
    Backgrounds_2018 = []
    Backgrounds_2018.append('DYJetsToLL_M-50_2018')
    Backgrounds_2018.append('TTTo2L2Nu_2018')
    Backgrounds_2018.append('WW_2018')
    Backgrounds_2018.append('WZ_2018')
    Backgrounds_2018.append('ZZ_2018')




    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1

    treeMC_2016     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_preVFP + Backgrounds_postVFP, 'MC'), name = 'MC', isdata = 0, close = True)
    treeMC_2017     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_2017, 'MC'), name = 'MC', isdata = 0, close = True)
    treeMC_2018     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_2018, 'MC'), name = 'MC', isdata = 0, close = True)

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1, close = True )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1, close = True)
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1, close = True)

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1, close = True)


    ############ Plot everything
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = 'hEESS_mass_Z', label_extra = 'SS background', treeEXTRA = treeDATA_EG2016, xlabel = '', outtag = '2016', yshift = 1000.0, LLlabel = 'Control Region: |#Delta#Phi| > #pi/2', leftlabel = '', xlog = False)
    makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2017, treeMC = treeMC_2017, inputdir = opts.input, hname_extra = 'hEESS_mass_Z', label_extra = 'SS background', treeEXTRA = treeDATA_EG2017, xlabel = '', outtag = '2017', yshift = 1000.0, LLlabel = 'Control Region: |#Delta#Phi| > #pi/2', leftlabel = '', xlog = False)
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEEBCR_mass_Z', hname_MC = 'hEEBCR_mass_Z', ylog = True, treeDATA = treeDATA_EG2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = 'hEESS_mass_Z', label_extra = 'SS background', treeEXTRA = treeDATA_EG2018, xlabel = '', outtag = '2018', yshift = 1000.0, LLlabel = 'Control Region: |#Delta#Phi| > #pi/2', leftlabel = '', xlog = False)


    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMMBCR_mass_Z', hname_MC = 'hMMBCR_mass_Z', ylog = True, treeDATA = treeDATA_Mu2016, treeMC = treeMC_2016, inputdir = opts.input, hname_extra = 'hMMSS_mass_Z', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2016, xlabel = '', outtag = '2016', yshift = 1000.0, LLlabel = 'Control Region: |#Delta#Phi| > #pi/2', leftlabel = '', xlog = False)
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMMBCR_mass_Z', hname_MC = 'hMMBCR_mass_Z', ylog = True, treeDATA = treeDATA_Mu2018, treeMC = treeMC_2018, inputdir = opts.input, hname_extra = 'hMMSS_mass_Z', label_extra = 'SS background', treeEXTRA = treeDATA_Mu2018, xlabel = '', outtag = '2018', yshift = 1000.0, LLlabel = 'Control Region: |#Delta#Phi| > #pi/2', leftlabel = '', xlog = False)
    
