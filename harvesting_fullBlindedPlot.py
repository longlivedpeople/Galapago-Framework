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




def makeBlindedPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, rebin = False, lines = [], xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = ''):


    ### Get histograms
    luminosity = lumi

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
    hbkg.SetLineColor(r.kCyan-6) 
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)

    for _i, _h in enumerate(hSIS.GetHists()):
        _h.Scale(lumi/35.87)
        s_histos.append(copy.deepcopy(_h))


    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])

    ### Set Maximum
    if not ylog:
        hbkg.SetMaximum(1.3*maxVal)
    else:
        hbkg.SetMaximum(100000.0*maxVal)
        hbkg.SetMinimum(0.0001)

    if ymax: hbkg.SetMaximum(ymax) 

    ### Count background events
    backtotal = 0.0
    for n in range(1, hbkg.GetNbinsX() + 1):
        if hbkg.GetBinLowEdge(n) > 6.0:
            backtotal += hbkg.GetBinContent(n)

    print(backtotal)



    ### Canvas object
    plot = Canvas.Canvas('Blinded_'+hname_bkg, 'png', 0.15, 0.6, 0.6, 0.89, 1)
    plot.addHisto(hbkg, 'HIST', 'Background (Data-driven)', 'f', '', 1, 0)
    
    ### Add signals:
    for i,_h in enumerate(s_histos):
        stotal_6 = 0.0
        stotal_2 = 0.0
        for n in range(1, _h.GetNbinsX() + 1):
            if _h.GetBinLowEdge(n) > 6.0:
                stotal_6 += _h.GetBinContent(n)
            if _h.GetBinLowEdge(n) > 2.0:
                stotal_2 += _h.GetBinContent(n)
                #print('bin', _h.GetBinContent(n))

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        print(masses, stotal_6, stotal_2)
        legend = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, hbkg.GetMaximum(), r.kBlack)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
    plot.save(1, 0, ylog, luminosity, '', outputDir = outdir, xlog = False)

    

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


    ############# EG data definition
    DoubleEGB = 'DoubleEG_Run2016B'
    DoubleEGC = 'DoubleEG_Run2016C'
    DoubleEGD = 'DoubleEG_Run2016D'
    DoubleEGE = 'DoubleEG_Run2016E'
    DoubleEGF = 'DoubleEG_Run2016F'
    DoubleEGG = 'DoubleEG_Run2016G'
    DoubleEGH = 'DoubleEG_Run2016H'

    DoubleEG_list = []
    DoubleEG_list.append(DoubleEGB)
    DoubleEG_list.append(DoubleEGC)
    DoubleEG_list.append(DoubleEGD)
    DoubleEG_list.append(DoubleEGE)
    DoubleEG_list.append(DoubleEGF)
    DoubleEG_list.append(DoubleEGG)
    DoubleEG_list.append(DoubleEGH)

    ############# Muon data definition
    DoubleMuonB = 'DoubleMuon_Run2016B'
    DoubleMuonC = 'DoubleMuon_Run2016C'
    DoubleMuonD = 'DoubleMuon_Run2016D'
    DoubleMuonE = 'DoubleMuon_Run2016E'
    DoubleMuonF = 'DoubleMuon_Run2016F'
    DoubleMuonG = 'DoubleMuon_Run2016G'
    DoubleMuonH = 'DoubleMuon_Run2016H'

    DoubleMuon_list = []
    DoubleMuon_list.append(DoubleMuonB)
    DoubleMuon_list.append(DoubleMuonC)
    DoubleMuon_list.append(DoubleMuonD)
    DoubleMuon_list.append(DoubleMuonE)
    DoubleMuon_list.append(DoubleMuonF)
    DoubleMuon_list.append(DoubleMuonG)
    DoubleMuon_list.append(DoubleMuonH)


    ############# Background definition
    Backgrounds = []
    Backgrounds.append('DYJetsToLL_M-50') 
    Backgrounds.append('DYJetsToLL_M-10to50') 
    Backgrounds.append('WW') 
    Backgrounds.append('WZ') 
    Backgrounds.append('ZZ') 
    Backgrounds.append('TT') 

    ############# Signal definition
    Signals = []
    Signals.append('HXX_400_50_4mm')
    Signals.append('HXX_400_50_40mm')
    Signals.append('HXX_400_50_400mm')
    Signals.append('HXX_400_150_400mm')
    Signals.append('HXX_1000_150_100mm')
    Signals.append('HXX_1000_150_10mm')
    Signals.append('HXX_1000_350_350mm')
    Signals.append('HXX_1000_350_35mm')

    Signals_400_50 = []
    Signals_400_50.append('HXX_400_50_1mm')
    Signals_400_50.append('HXX_400_50_10mm')
    Signals_400_50.append('HXX_400_50_100mm')
    Signals_400_50.append('HXX_400_50_1000mm')


    ############# Luminosity definition
    lumiB = 5.79
    lumiC = 2.57
    lumiD = 4.25
    lumiE = 4.01
    lumiF = 3.10
    lumiG = 7.54
    lumiH = 8.61

    lumiEG = {}
    lumiEG['DoubleEG_Run2016B'] = lumiB
    lumiEG['DoubleEG_Run2016C'] = lumiC
    lumiEG['DoubleEG_Run2016D'] = lumiD
    lumiEG['DoubleEG_Run2016E'] = lumiE
    lumiEG['DoubleEG_Run2016F'] = lumiF
    lumiEG['DoubleEG_Run2016G'] = lumiG
    lumiEG['DoubleEG_Run2016H'] = lumiH

    lumi_EG = 0.0
    for dataset in DoubleEG_list: lumi_EG += lumiEG[dataset]


    lumiMuon = {}
    lumiMuon['DoubleMuon_Run2016B'] = lumiB
    lumiMuon['DoubleMuon_Run2016C'] = lumiC
    lumiMuon['DoubleMuon_Run2016D'] = lumiD
    lumiMuon['DoubleMuon_Run2016E'] = lumiE
    lumiMuon['DoubleMuon_Run2016F'] = lumiF
    lumiMuon['DoubleMuon_Run2016G'] = lumiG
    lumiMuon['DoubleMuon_Run2016H'] = lumiH

    lumi_Muon = 0.0
    for dataset in DoubleMuon_list: lumi_Muon += lumiMuon[dataset]




    #filename = 'dat/Samples_cern_Legacy.dat'
    filename = 'dat/Samples_cern_fillingv2.dat'


    ### Tree SI Tree
    treeSI_400_50_Legacy = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/Samples_cern_Legacy.dat', Signals_400_50, 'SI'), name = 'SI', isdata = 0 )
    treeSI_400_50_Old = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/Samples_cern_fillinglow.dat', Signals, 'SI'), name = 'SI', isdata = 0 )

    ################################
    ######## DoubleEG Plots ########
    ################################
       
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )

    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEoffZSR_nBSEE', hname_bkg = 'hEEoffZCR_nBSEE', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_first', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEoffZSR_ne1_nBSEE', hname_bkg = 'hEEoffZCR_ne1_nBSEE', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_nLL', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEoffZSR_ng1_nBSEE', hname_bkg = 'hEEoffZCR_ng1_nBSEE', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_nLL', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEoffZSRne1_trackIxy', hname_bkg = 'hEEoffZCRne1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Legacy, xlabel = '', outtag = 'nLL-splitting_400_50_Legacy', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEfullZSRng1_trackIxy', hname_bkg = 'hEEfullZCRng1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Legacy, xlabel = '', outtag = 'nLL-splitting_400_50_Legacy', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEoffZSRne1_trackIxy', hname_bkg = 'hEEoffZCRne1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Old, xlabel = '', outtag = 'nLL-splitting_400_50_Old', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_EG, hname_SI = 'hEEfullZSRng1_trackIxy', hname_bkg = 'hEEfullZCRng1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Old, xlabel = '', outtag = 'nLL-splitting_400_50_Old', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '') 


    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )

    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMoffZSR_nBSMM', hname_bkg = 'hMMoffZCR_nBSMM', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_first', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMoffZSR_ne1_nBSMM', hname_bkg = 'hMMoffZCR_ne1_nBSMM', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_nLL', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMoffZSR_ng1_nBSMM', hname_bkg = 'hMMoffZCR_ng1_nBSMM', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_nLL', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMoffZSR_ne1_trackIxy_log', hname_bkg = 'hMMoffZCR_ne1_trackIxy_log', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_nLLnewlog', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMfullZSR_ng1_trackIxy_log', hname_bkg = 'hMMfullZCR_ng1_trackIxy_log', ylog = True, treeDATA = treeDATA, inputdir = 'histograms_nLLnewlog', treeSI = treeSI, xlabel = '', outtag = 'nLL-splitting', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMoffZSRne1_trackIxy', hname_bkg = 'hMMoffZCRne1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Legacy, xlabel = '', outtag = 'nLL-splitting_400_50_Legacy', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMfullZSRng1_trackIxy', hname_bkg = 'hMMfullZCRng1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Legacy, xlabel = '', outtag = 'nLL-splitting_400_50_Legacy', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMoffZSRne1_trackIxy', hname_bkg = 'hMMoffZCRne1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Old, xlabel = '', outtag = 'nLL-splitting_400_50_Old', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
    #makeBlindedPlot(lumi = lumi_Muon, hname_SI = 'hMMfullZSRng1_trackIxy', hname_bkg = 'hMMfullZCRng1_trackIxy', ylog = True, treeDATA = treeDATA, inputdir = opts.input, treeSI = treeSI_400_50_Old, xlabel = '', outtag = 'nLL-splitting_400_50_Old', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '') 
