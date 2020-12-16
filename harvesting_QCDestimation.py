import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager
#from include.Utils import *




def makeQCDEstimationPlot(lumi, hname, ylog, treeDATAs, labels, inputdir, xlabel = '', outtag = '', yshift = 0.0, LLlabel = ''):

    luminosity = lumi
    SShname = hname.split('__')[0] + '_SS'

    #colors = [r.kRed, r.kOrange, r.kYellow, r.kGreen+1, r.kAzure, r.kBlue+2, r.kMagenta+2]
    colors = [r.kYellow, r.kOrange, r.kRed, r.kMagenta+2, r.kBlue+2, r.kAzure, r.kGreen+2]

    hSS = r.THStack('hSS_%s'%(hname), '') # SS stacked 
    for i,tree in enumerate(treeDATAs):
        _h = tree.getLoopTH1F(inputdir, SShname)
        _h.SetFillColorAlpha(colors[i], 0.5)
        _h.SetLineColor(colors[i])
        _h.SetTitle(labels[i].split('_')[1])
        hSS.Add(_h)
        if not i:
            hSS_total = copy.deepcopy(_h)
        else:
            hSS_total.Add(_h)

    ### Axis definition
    """
    can_aux = TCanvas("can_%s"%(hname))
    can_aux.cd()
    hSS.Draw()
    hSS.GetXaxis().SetTitle(hSS_total.GetXaxis().GetTitle())
    hSS.GetYaxis().SetTitle(hSS_total.GetYaxis().GetTitle())
    del can_aux
    """
    hSS.SetTitle(';'+hSS_total.GetXaxis().GetTitle()+';'+hSS_total.GetYaxis().GetTitle())

    
    maxVal = hSS_total.GetMaximum()
    
    ### Set Maximum
    if not ylog:
        hSS.SetMaximum(1.4*maxVal)
        hSS.SetMinimum(0.0)
    else:
        if not yshift:
            hSS.SetMaximum(10.0*maxVal)
        else:
            hSS.SetMaximum(yshift*maxVal)
        hSS.SetMinimum(0.1)

    ### -> Canvas object
    plot = Canvas.Canvas(hname, 'png', 0.65, 0.65, 0.9, 0.9, 1)
    
    ### Add background:
    plot.addStack(hSS, 'HIST', 1, 0) # Background
    
    ### DATAlabel banner:
    plot.addLatex(0.17, 0.85, '2016 Era comparison', font = 62, size = 0.034)
    
    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.81, 'e^{+}e^{-} channel', font = 42, size = 0.034)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.81, '#mu^{+}#mu^{-} channel', font = 42, size = 0.034)

    ### Extralabel:
    plot.addLatex(0.17, 0.76, 'SS Region: |#Delta#Phi| > #pi/2, q_{1}*q_{2} > 0', font = 42, size = 0.03)

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/Plots_' + outtag + '/'
    plot.save(1, 0, ylog, luminosity, '', outputDir = outdir)

    del plot
    return 

################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-m', '--Muoninput', action='store', type=str, dest='Muoninput', default='', help='Target directory')
    parser.add_option('-e', '--EGinput', action='store', type=str, dest='EGinput', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C+')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()


    ############# Muon data definition
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
    #Signals.append('DisplacedSUSY_350_148_173')
    Signals.append('HXX_400_50_400mm')
    #Signals.append('HXX_400_50_40mm')
    #Signals.append('HXX_400_50_4mm')

    ############# Parameter definition
    lumiB = 5.79
    lumiC = 2.57
    lumiD = 4.25
    lumiE = 4.01
    lumiF = 3.10
    lumiG = 7.54
    lumiH = 8.61
    lumi_total = lumiB + lumiC + lumiD + lumiE + lumiF + lumiG + lumiH # luminosity
    lumi_list = [lumiB, lumiC, lumiD, lumiE, lumiF, lumiG, lumiH]
    

    filename = 'dat/Samples_cern_filling.dat'

    letter = ['B', 'C', 'D', 'E', 'F', 'G', 'H']


    ### -> Muon plots:
    MuonTrees = []
    for i in range(0, len(DoubleMuon_list)):
        dataset = DoubleMuon_list[i]
        lumi = lumi_list[i]

        treeEra = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, [dataset], 'DATA'), name = 'DATA', isdata = 1 )
        MuonTrees.append(treeEra)


    makeQCDEstimationPlot(lumi_total, 'hMM_trackIxy', True, MuonTrees, DoubleMuon_list, WORKPATH + opts.Muoninput, outtag = 'QCDEstimation', LLlabel = 'MM')
    makeQCDEstimationPlot(lumi_total, 'hMM_trackDxy', True, MuonTrees, DoubleMuon_list, WORKPATH + opts.Muoninput, outtag = 'QCDEstimation', LLlabel = 'MM')
    makeQCDEstimationPlot(lumi_total, 'hMM_Ixy', True, MuonTrees, DoubleMuon_list, WORKPATH + opts.Muoninput, outtag = 'QCDEstimation', LLlabel = 'MM')
    makeQCDEstimationPlot(lumi_total, 'hMM_Lxy', True, MuonTrees, DoubleMuon_list, WORKPATH + opts.Muoninput, outtag = 'QCDEstimation', LLlabel = 'MM')


    ### -> Electron plots:
    EGTrees = []
    for i in range(0, len(DoubleEG_list)):
        dataset = DoubleEG_list[i]
        lumi = lumi_list[i]

        treeEra = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, [dataset], 'DATA'), name = 'DATA', isdata = 1 )
        EGTrees.append(treeEra)


    makeQCDEstimationPlot(lumi_total, 'hEE_trackIxy', True, EGTrees, DoubleMuon_list, WORKPATH + opts.EGinput, outtag = 'QCDEstimation', LLlabel = 'EE')
    makeQCDEstimationPlot(lumi_total, 'hEE_trackDxy', True, EGTrees, DoubleMuon_list, WORKPATH + opts.EGinput, outtag = 'QCDEstimation', LLlabel = 'EE')
    makeQCDEstimationPlot(lumi_total, 'hEE_Ixy', True, EGTrees, DoubleMuon_list, WORKPATH + opts.EGinput, outtag = 'QCDEstimation', LLlabel = 'EE')
    makeQCDEstimationPlot(lumi_total, 'hEE_Lxy', True, EGTrees, DoubleMuon_list, WORKPATH + opts.EGinput, outtag = 'QCDEstimation', LLlabel = 'EE')


