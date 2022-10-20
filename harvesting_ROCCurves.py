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




def makeROCCurves(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, rebin = False, lines = [], xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', xlog = False):


    ### Get histograms
    luminosity = lumi

    hbkg = treeDATA.getLoopTH1F(inputdir, hname_bkg, doOF = False)
   
    ### Signal histograms
    s_histos = []
    hSIS = treeSI.getLoopStack(inputdir, hname_SI, doOF = False)
    for _i, _h in enumerate(hSIS.GetHists()):
        s_histos.append(copy.deepcopy(_h))

    ## Get the ROC values
    nMax = hbkg.GetNbinsX()
    nBkgMax = hbkg.Integral() 
    nSigMax = [_h.Integral() for _h in s_histos]
    BkgRej = []
    SigEffs = [ [] for _h in hSIS.GetHists() ]
    for n in range(1, hbkg.GetNbinsX() + 1):
        BkgRej.append( hbkg.Integral(n + 1, hbkg.GetNbinsX() + 1) / nBkgMax )
        for _i, _h in enumerate(s_histos):
            SigEffs[_i].append( _h.Integral(1, n) / nSigMax[_i] )

    ## Get the ROC TGraphs
    plot = Canvas.Canvas('ROCs_'+hname_bkg, 'png', 0.15, 0.2, 0.4, 0.4, 1)
    aux_h = r.TH2F('aux', ';Background rejection;Signal efficiency', 1, 0., 1., 1, 0., 1.)
    plot.addHisto(aux_h, 'HIST', '', 'l', _h.GetFillColor(), 1, 6) # Signal
    for _i, _h in enumerate(s_histos):
        tgraph = r.TGraph(len(BkgRej), np.array(BkgRej), np.array(SigEffs[_i])) 
        tgraph.SetLineWidth(2)
        masses = eval(_h.GetTitle()[3:])
        legend = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
        plot.addGraph(tgraph, 'C, SAME', legend, 'l', _h.GetFillColor(), 1, _i) # Signal

    for n in range(1, hbkg.GetNbinsX()):
        print(hbkg.GetBinLowEdge(n), BkgRej[n], SigEffs[2][n], SigEffs[3][n])

    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/ROCPlots_' + outtag + '/'
    plot.save(1, 0, ylog, luminosity, '', outputDir = outdir, xlog = False, is2d = True)

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

    Signals_1000_150_2016 = []
    Signals_1000_150_2016.append('HSS_1000_150_1_2016')
    Signals_1000_150_2016.append('HSS_1000_150_10_2016')
    Signals_1000_150_2016.append('HSS_1000_150_100_2016')
    Signals_1000_150_2016.append('HSS_1000_150_1000_2016')
    Signals_1000_150_2016.append('HSS_1000_150_10000_2016')
    Signals_1000_150_2017 = []
    Signals_1000_150_2017.append('HSS_1000_150_1_2017')
    Signals_1000_150_2017.append('HSS_1000_150_10_2017')
    Signals_1000_150_2017.append('HSS_1000_150_100_2017')
    Signals_1000_150_2017.append('HSS_1000_150_1000_2017')
    Signals_1000_150_2017.append('HSS_1000_150_10000_2017')
    Signals_1000_150_2018 = []
    Signals_1000_150_2018.append('HSS_1000_150_1_2018')
    Signals_1000_150_2018.append('HSS_1000_150_10_2018')
    Signals_1000_150_2018.append('HSS_1000_150_100_2018')
    Signals_1000_150_2018.append('HSS_1000_150_1000_2018')
    Signals_1000_150_2018.append('HSS_1000_150_10000_2018')


    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1



    filename = 'dat/Samples_cern_UltraLegacy.dat'


    ### Tree SI Tree
    treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2016UL_Summer22.dat', Signals2016, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2017UL_Summer22.dat', Signals2017, 'SI'), name = 'SI', isdata = 0 )
    treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2018UL_Summer22.dat', Signals2018, 'SI'), name = 'SI', isdata = 0 )


    ##################################
    ######## DoubleMuon Plots ########
    ##################################
    
    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    makeROCCurves(lumi = lumi2016, hname_SI = 'hMMSRdisp_normalizedChi2_log', hname_bkg = 'hMMBCRdisp_normalizedChi2_log', ylog = False, treeDATA = treeDATA_Mu2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = 'ROC2016', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True) 
    makeROCCurves(lumi = lumi2018, hname_SI = 'hMMSRdisp_normalizedChi2_log', hname_bkg = 'hMMBCRdisp_normalizedChi2_log', ylog = False, treeDATA = treeDATA_Mu2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = 'ROC2018', ymax = 0.0, LLlabel = 'MM', DATAlabel = '', extralabel = '', xlog = True) 

    ################################
    ######## DoubleEG Plots ########
    ################################

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    makeROCCurves(lumi = lumi2016, hname_SI = 'hEESRdisp_normalizedChi2_log', hname_bkg = 'hEEBCRdisp_normalizedChi2_log', ylog = False, treeDATA = treeDATA_EG2016, inputdir = opts.input, treeSI = treeSI_2016, xlabel = '', outtag = 'ROC2016', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
    makeROCCurves(lumi = lumi2017, hname_SI = 'hEESRdisp_normalizedChi2_log', hname_bkg = 'hEEBCRdisp_normalizedChi2_log', ylog = False, treeDATA = treeDATA_EG2017, inputdir = opts.input, treeSI = treeSI_2017, xlabel = '', outtag = 'ROC2017', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
    makeROCCurves(lumi = lumi2018, hname_SI = 'hEESRdisp_normalizedChi2_log', hname_bkg = 'hEEBCRdisp_normalizedChi2_log', ylog = False, treeDATA = treeDATA_EG2018, inputdir = opts.input, treeSI = treeSI_2018, xlabel = '', outtag = 'ROC2018', ymax = 0.0, LLlabel = 'EE', DATAlabel = '', extralabel = '', xlog = True) 
