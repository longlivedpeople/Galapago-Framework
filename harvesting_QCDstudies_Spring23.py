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
from include.Utils import makeAgreementTest, makeBackgroundValidationPlot


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
    r.gStyle.SetOptFit(0)


    ############# Load data 
    DoubleMuon2016_HIPM = []
    DoubleMuon2016_noHIPM = []
    DoubleMuon2018 = []
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016B_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016C_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016D_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016E_HIPM')
    DoubleMuon2016_HIPM.append('DoubleMuon_Run2016F_HIPM')
    DoubleMuon2016_noHIPM.append('DoubleMuon_Run2016F_noHIPM')
    DoubleMuon2016_noHIPM.append('DoubleMuon_Run2016G_noHIPM')
    DoubleMuon2016_noHIPM.append('DoubleMuon_Run2016H_noHIPM')
    DoubleMuon2018.append('DoubleMuon_Run2018A')
    DoubleMuon2018.append('DoubleMuon_Run2018B')
    DoubleMuon2018.append('DoubleMuon_Run2018C')
    DoubleMuon2018.append('DoubleMuon_Run2018D')

    DoubleEG2016_HIPM = []
    DoubleEG2016_noHIPM = []
    DoubleEG2017 = []
    DoubleEG2018 = []
    DoubleEG2016_HIPM.append('DoubleEG_Run2016B_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016C_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016D_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016E_HIPM')
    DoubleEG2016_HIPM.append('DoubleEG_Run2016F_HIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016F_noHIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG2016_noHIPM.append('DoubleEG_Run2016H_noHIPM')
    DoubleEG2017.append('DoubleEG_Run2017B')
    DoubleEG2017.append('DoubleEG_Run2017C')
    DoubleEG2017.append('DoubleEG_Run2017D')
    DoubleEG2017.append('DoubleEG_Run2017E')
    DoubleEG2017.append('DoubleEG_Run2017F')
    DoubleEG2018.append('EGamma_Run2018A')
    DoubleEG2018.append('EGamma_Run2018B')
    DoubleEG2018.append('EGamma_Run2018C')
    DoubleEG2018.append('EGamma_Run2018D')

    ############# Load background simulation
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



    ################ Define the .dat file
    filename = 'dat/Samples_cern_UltraLegacy_Spring23.dat'

    lumi_total = ""

    rebin1 = np.linspace(0, 3.14/2., 16)

    inputDir = 'histograms_Spring23_FullQCD_upperCut'
    www = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/Background-studies/QCD-studies-Full_upperCut/'

    plotMuons = True
    plotElectrons = False


    #
    # -- Muon plots
    #
    if plotMuons:

        tree_MM2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_HIPM, 'MC'), name = 'DATA', isdata = 1, close = True)
        tree_MM2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_noHIPM, 'MC'), name = 'DATA', isdata = 1, close = True)
        tree_MM2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016_HIPM + DoubleMuon2016_noHIPM, 'MC'), name = 'DATA', isdata = 1, close = True)
        tree_MM2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'MC'), name = 'DATA', isdata = 1, close = True)


        # QCD Control Regioni (2016) #

        makeAgreementTest('', 'hMMQCD_dPhi', 'hMMQCD_dPhi_inv', 0, tree_MM2016_HIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region', 'HIP affected runs', 'QCD2016_HIPM_hMM_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest('', 'hMMQCDDisp_dPhi', 'hMMQCDDisp_dPhi_inv', 0, tree_MM2016_HIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP affected runs', 'QCD2016_HIPM_hMM_dPhi_Disp_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)
        makeAgreementTest('', 'hMMQCD_dPhi', 'hMMQCD_dPhi_inv', 0, tree_MM2016_noHIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region', 'HIP unaffected runs', 'QCD2016_noHIPM_hMM_dPhi_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest('', 'hMMQCDDisp_dPhi', 'hMMQCDDisp_dPhi_inv', 0, tree_MM2016_noHIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP unaffected runs', 'QCD2016_noHIPM_hMM_dPhi_Disp_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)

        makeBackgroundValidationPlot('2016_hMMQCDSRDisp_mass', 35.9, 'hMMQCDSRDisp_mass', 'hMMQCDBCRDisp_mass', 0, tree_MM2016, inputDir, False, 0.0, '', www, 0.0, 'MM', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSRDisp_cosAlpha', 35.9, 'hMMQCDSRDisp_cosAlpha', 'hMMQCDBCRDisp_cosAlpha', 0, tree_MM2016, inputDir, False, 0.0, '', www, 0.0, 'MM', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSRDisp_normalizedChi2', 35.9, 'hMMQCDSRDisp_normalizedChi2', 'hMMQCDBCRDisp_normalizedChi2', 0, tree_MM2016, inputDir, False, 0.0, '', www, 0.0, 'MM', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSR_Lxy', 35.9, 'hMMQCDSR_Lxy', 'hMMQCDBCR_Lxy', 1, tree_MM2016, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSR_trackIxy', 35.9, 'hMMQCDSR_trackIxy', 'hMMQCDBCR_trackIxy', 1, tree_MM2016, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSSSR_Lxy', 35.9, 'hMMQCDSSSR_Lxy', 'hMMQCDSSBCR_Lxy', 1, tree_MM2016, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSSSR_trackIxy', 35.9, 'hMMQCDSSSR_trackIxy', 'hMMQCDSSBCR_trackIxy', 1, tree_MM2016, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDOSSR_Lxy', 35.9, 'hMMQCDOSSR_Lxy', 'hMMQCDOSBCR_Lxy', 1, tree_MM2016, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDOSSR_trackIxy', 35.9, 'hMMQCDOSSR_trackIxy', 'hMMQCDOSBCR_trackIxy', 1, tree_MM2016, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSR_eta', 35.9, 'hMMQCDSR_eta', 'hMMQCDBCR_eta', 1, tree_MM2016, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSRDisp_dR', 35.9, 'hMMQCDSRDisp_dR', 'hMMQCDBCRDisp_dR', 1, tree_MM2016, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMQCDSRDisp_ptDiff', 35.9, 'hMMQCDSRDisp_ptDiff', 'hMMQCDBCRDisp_ptDiff', 1, tree_MM2016, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)

        # QCD Control Regioni (2018) #

        makeAgreementTest(59.8, 'hMMQCD_dPhi', 'hMMQCD_dPhi_inv', 0, tree_MM2018, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region', '', 'QCD2018_hMM_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 4, outpath = www)
        makeAgreementTest(59.8, 'hMMQCDDisp_dPhi', 'hMMQCDDisp_dPhi_inv', 0, tree_MM2018, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', '', 'QCD2018_hMM_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)

        makeBackgroundValidationPlot('2018_hMMQCDSRDisp_mass', 59.8, 'hMMQCDSRDisp_mass', 'hMMQCDBCRDisp_mass', 0, tree_MM2018, inputDir, False, 0.0, '', www, 0.0, 'MM', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSRDisp_cosAlpha', 59.8, 'hMMQCDSRDisp_cosAlpha', 'hMMQCDBCRDisp_cosAlpha', 0, tree_MM2018, inputDir, False, 0.0, '', www, 0.0, 'MM', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSRDisp_normalizedChi2', 59.8, 'hMMQCDSRDisp_normalizedChi2', 'hMMQCDBCRDisp_normalizedChi2', 0, tree_MM2018, inputDir, False, 0.0, '', www, 0.0, 'MM', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSR_Lxy', 59.8, 'hMMQCDSR_Lxy', 'hMMQCDBCR_Lxy', 1, tree_MM2018, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSR_trackIxy', 59.8, 'hMMQCDSR_trackIxy', 'hMMQCDBCR_trackIxy', 1, tree_MM2018, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSSSR_Lxy', 59.8, 'hMMQCDSSSR_Lxy', 'hMMQCDSSBCR_Lxy', 1, tree_MM2018, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSSSR_trackIxy', 59.8, 'hMMQCDSSSR_trackIxy', 'hMMQCDSSBCR_trackIxy', 1, tree_MM2018, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDOSSR_Lxy', 59.8, 'hMMQCDOSSR_Lxy', 'hMMQCDOSBCR_Lxy', 1, tree_MM2018, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDOSSR_trackIxy', 59.8, 'hMMQCDOSSR_trackIxy', 'hMMQCDOSBCR_trackIxy', 1, tree_MM2018, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSR_eta', 59.8, 'hMMQCDSR_eta', 'hMMQCDBCR_eta', 1, tree_MM2018, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSRDisp_dR', 59.8, 'hMMQCDSRDisp_dR', 'hMMQCDBCRDisp_dR', 1, tree_MM2018, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMQCDSRDisp_ptDiff', 59.8, 'hMMQCDSRDisp_ptDiff', 'hMMQCDBCRDisp_ptDiff', 1, tree_MM2018, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'QCD Control Region', False, 0.1)

        # SS Control Regioni (2016) #

        makeAgreementTest('', 'hMMSS_dPhi', 'hMMSS_dPhi_inv', 0, tree_MM2016_HIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'SS Control Region', 'HIP affected runs', 'SS2016_HIPM_hMM_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.4, rmax= 1.6, maxY = 0, outpath = www)
        makeAgreementTest('', 'hMMSSDisp_dPhi', 'hMMSSDisp_dPhi_inv', 0, tree_MM2016_HIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP affected runs', 'SS2016_HIPM_hMM_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.0, rmax= 2.0, maxY = 0, outpath = www)
        makeAgreementTest('', 'hMMSS_dPhi', 'hMMSS_dPhi_inv', 0, tree_MM2016_noHIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Same-sign Control Region', 'HIP unaffected runs', 'SS2016_noHIPM_hMM_dPhi_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.4, rmax= 1.6, maxY = 0, outpath = www)
        makeAgreementTest('', 'hMMSSDisp_dPhi', 'hMMSSDisp_dPhi_inv', 0, tree_MM2016_noHIPM, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Same-sign Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP unaffected runs', 'SS2016_noHIPM_hMM_dPhi_Disp_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.0, rmax= 2.0, maxY = 0, outpath = www)


        makeBackgroundValidationPlot('2016_hMMSSSRDisp_mass', 35.9, 'hMMSSSRDisp_mass', 'hMMSSBCRDisp_mass', 0, tree_MM2016, inputDir, False, 0.0, '', www, 0.0, 'MM', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSRDisp_cosAlpha', 35.9, 'hMMSSSRDisp_cosAlpha', 'hMMSSBCRDisp_cosAlpha', 0, tree_MM2016, inputDir, False, 0.0, '', www, 0.0, 'MM', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSRDisp_normalizedChi2', 35.9, 'hMMSSSRDisp_normalizedChi2', 'hMMSSBCRDisp_normalizedChi2', 0, tree_MM2016, inputDir, False, 0.0, '', www, 0.0, 'MM', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSR_Lxy', 35.9, 'hMMSSSR_Lxy', 'hMMSSBCR_Lxy', 1, tree_MM2016, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSR_trackIxy', 35.9, 'hMMSSSR_trackIxy', 'hMMSSBCR_trackIxy', 1, tree_MM2016, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSR_eta', 35.9, 'hMMSSSR_eta', 'hMMSSBCR_eta', 1, tree_MM2016, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSRDisp_dR', 35.9, 'hMMSSSRDisp_dR', 'hMMSSBCRDisp_dR', 1, tree_MM2016, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2016_hMMSSSRDisp_ptDiff', 35.9, 'hMMSSSRDisp_ptDiff', 'hMMSSBCRDisp_ptDiff', 1, tree_MM2016, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)


        # SS Control Region (2018) #

        makeAgreementTest(59.8, 'hMMSS_dPhi', 'hMMSS_dPhi_inv', 0, tree_MM2018, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Same-sign Control Region', '', 'SS2018_hMM_dPhi_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 4, outpath = www)
        makeAgreementTest(59.8, 'hMMSSDisp_dPhi', 'hMMSSDisp_dPhi_inv', 0, tree_MM2018, inputDir, '|#Delta#Phi(#mu^{+}, #mu^{-})|', '#pi - |#Delta#Phi(#mu^{+}, #mu^{-})|', 'Same-sign Control Region, |d_{0}|/#sigma_{d}  > 2.5', '', 'SS2018_hMM_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.6, rmax= 1.4, maxY = 0, outpath = www)

        makeBackgroundValidationPlot('2018_hMMSSSRDisp_mass', 59.8, 'hMMSSSRDisp_mass', 'hMMSSBCRDisp_mass', 0, tree_MM2018, inputDir, False, 0.0, '', www, 0.0, 'MM', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSRDisp_cosAlpha', 59.8, 'hMMSSSRDisp_cosAlpha', 'hMMSSBCRDisp_cosAlpha', 0, tree_MM2018, inputDir, False, 0.0, '', www, 0.0, 'MM', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSRDisp_normalizedChi2', 59.8, 'hMMSSSRDisp_normalizedChi2', 'hMMSSBCRDisp_normalizedChi2', 0, tree_MM2018, inputDir, False, 0.0, '', www, 0.0, 'MM', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSR_Lxy', 59.8, 'hMMSSSR_Lxy', 'hMMSSBCR_Lxy', 1, tree_MM2018, inputDir, 5, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSR_trackIxy', 59.8, 'hMMSSSR_trackIxy', 'hMMSSBCR_trackIxy', 1, tree_MM2018, inputDir, 2, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSR_eta', 59.8, 'hMMSSSR_eta', 'hMMSSBCR_eta', 1, tree_MM2018, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSRDisp_dR', 59.8, 'hMMSSSRDisp_dR', 'hMMSSBCRDisp_dR', 1, tree_MM2018, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)
        makeBackgroundValidationPlot('2018_hMMSSSRDisp_ptDiff', 59.8, 'hMMSSSRDisp_ptDiff', 'hMMSSBCRDisp_ptDiff', 1, tree_MM2018, inputDir, 1, 0.0, '', www, 100.0, 'MM', 'SS Control Region', False, 0.1)


    #
    # -- Electron plots
    #
    if plotElectrons:

        tree_EE2016_HIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_HIPM, 'MC'), name = 'DATA', isdata = 1, close = True)
        tree_EE2016_noHIPM = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016_noHIPM, 'MC'), name = 'DATA', isdata = 1, close = True)
        tree_EE2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'MC'), name = 'DATA', isdata = 1, close = True)
        tree_EE2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2018, 'MC'), name = 'DATA', isdata = 1, close = True)

        # QCD Control Region (2016) #
        makeAgreementTest('', 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2016_HIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', 'HIP affected runs', 'QCD2016_HIPM_hEE_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest('', 'hEEQCDDisp_dPhi', 'hEEQCDDisp_dPhi_inv', 0, tree_EE2016_HIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP affected runs', 'QCD2016_HIPM_hEE_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)
        makeAgreementTest('', 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2016_noHIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', 'HIP unaffected runs', 'QCD2016_noHIPM_hEE_dPhi_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest('', 'hEEQCDDisp_dPhi', 'hEEQCDDisp_dPhi_inv', 0, tree_EE2016_noHIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP unaffected runs', 'QCD2016_noHIPM_hEE_dPhi_Disp_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)
    
        # QCD Control Region (2017) #

        makeAgreementTest(41.5, 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2017, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', '', 'QCD2017_hEE_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest(41.5, 'hEEQCDDisp_dPhi', 'hEEQCDDisp_dPhi_inv', 0, tree_EE2017, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', '', 'QCD2017_hEE_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)

        # QCD Control Region (2018) #

        makeAgreementTest(54.5, 'hEEQCD_dPhi', 'hEEQCD_dPhi_inv', 0, tree_EE2018, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region', '', 'QCD2018_hEE_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 4, outpath = www)
        makeAgreementTest(54.5, 'hEEQCDDisp_dPhi', 'hEEQCDDisp_dPhi_inv', 0, tree_EE2018, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', '', 'QCD2018_hEE_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.5, rmax= 1.5, maxY = 0, outpath = www)

        makeBackgroundValidationPlot('2018_hEEQCDSRDisp_mass', 54.5, 'hEEQCDSRDisp_mass', 'hEEQCDBCRDisp_mass', 0, tree_EE2018, inputDir, False, 0.0, '', www, 0.0, 'EE', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hEEQCDSRDisp_normalizedChi2', 54.5, 'hEEQCDSRDisp_normalizedChi2', 'hEEQCDBCRDisp_normalizedChi2', 0, tree_EE2018, inputDir, False, 0.0, '', www, 0.0, 'EE', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hEEQCDSR_Lxy', 54.5, 'hEEQCDSR_Lxy', 'hEEQCDBCR_Lxy', 1, tree_EE2018, inputDir, 5, 0.0, '', www, 100.0, 'EE', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hEEQCDSR_trackIxy', 54.5, 'hEEQCDSR_trackIxy', 'hEEQCDBCR_trackIxy', 1, tree_EE2018, inputDir, 2, 0.0, '', www, 100.0, 'EE', 'QCD Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)


        # SS Control Region (2016) #

        makeAgreementTest('', 'hEESS_dPhi', 'hEESS_dPhi_inv', 0, tree_EE2016_HIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region', 'HIP affected runs', 'SS2016_HIPM_hEE_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest('', 'hEESSDisp_dPhi', 'hEESSDisp_dPhi_inv', 0, tree_EE2016_HIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP affected runs', 'SS2016_HIPM_hEE_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.6, rmax= 1.4, maxY = 0, outpath = www)

        makeAgreementTest('', 'hEESS_dPhi', 'hEESS_dPhi_inv', 0, tree_EE2016_noHIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region', 'HIP unaffected runs', 'SS2016_noHIPM_hEE_dPhi_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest('', 'hEESSDisp_dPhi', 'hEESSDisp_dPhi_inv', 0, tree_EE2016_noHIPM, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region, |d_{0}|/#sigma_{d}  > 2.5', 'HIP unaffected runs', 'SS2016_noHIPM_hEE_dPhi_Disp_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.6, rmax= 1.4, maxY = 0, outpath = www)
    
        # SS Control Region (2017) #

        makeAgreementTest(41.5, 'hEESS_dPhi', 'hEESS_dPhi_inv', 0, tree_EE2017, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region', '', 'SS2017_hEE_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 0, outpath = www)
        makeAgreementTest(41.5, 'hEESSDisp_dPhi', 'hEESSDisp_dPhi_inv', 0, tree_EE2017, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region, |d_{0}|/#sigma_{d}  > 2.5', '', 'SS2017_hEE_dPhi_Disp_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.6, rmax= 1.4, maxY = 0, outpath = www)

        # SS Control Region (2018) #

        makeAgreementTest(59.8, 'hEESS_dPhi', 'hEESS_dPhi_inv', 0, tree_EE2018, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region', '', 'SS2018_hEE_dPhi_corr', 1, sys = 0.1, ranges = [0, 1.56], rebin = 3, rmin = 0.8, rmax= 1.2, maxY = 4, outpath = www)
        makeAgreementTest(59.8, 'hEESSDisp_dPhi', 'hEESSDisp_dPhi_inv', 0, tree_EE2018, inputDir, '|#Delta#Phi(e^{+}, e^{-})|', '#pi - |#Delta#Phi(e^{+}, e^{-})|', 'Same-sign Control Region, |d_{0}|/#sigma_{d}  > 2.5', '', 'SS2018_hEE_dPhi_Disp_corr', 1, sys = 0.1,  ranges = [0, 1.56], rebin = 3, rmin = 0.6, rmax= 1.4, maxY = 0, outpath = www)

        makeBackgroundValidationPlot('2018_hEESSSRDisp_mass', 54.5, 'hEESSSRDisp_mass', 'hEESSBCRDisp_mass', 0, tree_EE2018, inputDir, False, 0.0, '', www, 0.0, 'EE', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hEESSSRDisp_normalizedChi2', 54.5, 'hEESSSRDisp_normalizedChi2', 'hEESSBCRDisp_normalizedChi2', 0, tree_EE2018, inputDir, False, 0.0, '', www, 0.0, 'EE', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hEESSSR_Lxy', 54.5, 'hEESSSR_Lxy', 'hEESSBCR_Lxy', 1, tree_EE2018, inputDir, 5, 0.0, '', www, 100.0, 'EE', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)
        makeBackgroundValidationPlot('2018_hEESSSR_trackIxy', 54.5, 'hEESSSR_trackIxy', 'hEESSBCR_trackIxy', 1, tree_EE2018, inputDir, 2, 0.0, '', www, 100.0, 'EE', 'SS Control Region, |d_{0}|/#sigma_{d}  > 2.5', False, 0.1)




