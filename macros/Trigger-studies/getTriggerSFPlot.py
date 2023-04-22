import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
import math, sys, optparse, array, copy, os
import gc, inspect
import numpy as np



################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]:
    WORKPATH += level
    WORKPATH += '/'

GALAPAGOPATH = ''
for d in WORKPATH.split('/'):
    GALAPAGOPATH += d
    GALAPAGOPATH += '/'
    if d == 'Galapago-Framework': break

sys.path.insert(0, GALAPAGOPATH)

import include.Canvas as Canvas

print(WORKPATH, WORKPATH)
print(GALAPAGOPATH, GALAPAGOPATH)

import include.Canvas as Canvas

##################################### FUNCTION DEFINITION ########################################

def getObject(filename, key):

    _f = r.TFile(filename)
    _h = _f.Get(key)
    _hcopy = copy.deepcopy(_h)
    _f.Close()

    return _hcopy


def rebinAxis(eff, axis):

    totals = eff.GetTotalHistogram()
    passed = eff.GetPassedHistogram()
    totals_rebin = totals.Rebin(len(axis)-1, totals.GetName()+'_rebined', axis)
    passed_rebin = passed.Rebin(len(axis)-1, passed.GetName()+'_rebined', axis)
    neweff = r.TEfficiency(passed_rebin, totals_rebin)

    c1 = r.TCanvas()
    neweff.Draw('AP')
#    c1.SaveAs('Pruebita.png')

    return(neweff)


def getHistoFromEff(eff):

    histo = eff.GetTotalHistogram().Clone()
    histo.Reset()

    for i in range(1, histo.GetNbinsX() + 1):
        error = max([eff.GetEfficiencyErrorLow(i), eff.GetEfficiencyErrorUp(i)])
        histo.SetBinContent(i, eff.GetEfficiency(i))
        histo.SetBinError(i, error)

    histo.SetMaximum(1.2)
    histo.SetMinimum(0.0)

    return histo


if __name__ == "__main__":


    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    print('WORKPATH: ' + WORKPATH)
    r.setTDRStyle()

    ### Get plots (2016)
    Efficiency_HLT_Full_pt_DYJetsToLL_M50 = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/FullComparisons/PhotonTrigger-SFs_Mass15_DYOptimized/TH1F_photontrigger_2016.root', 'Efficiency_HLT_Full_pt_DYJetsToLL_M-50')
    Efficiency_HLT_Full_pt_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/FullComparisons/PhotonTrigger-SFs_Mass15_DYOptimized/TH1F_photontrigger_2016.root', 'Efficiency_HLT_Full_pt_DATA')

    efficiencyDY = Efficiency_HLT_Full_pt_DYJetsToLL_M50.GetPassedHistogram()
    totalDY = Efficiency_HLT_Full_pt_DYJetsToLL_M50.GetTotalHistogram()
    efficiencyDY.Divide(totalDY)

    efficiencyData = Efficiency_HLT_Full_pt_DATA.GetPassedHistogram()
    totalData = Efficiency_HLT_Full_pt_DATA.GetTotalHistogram()
    efficiencyData.Divide(totalData)

    EE_SF2016 = efficiencyData.Clone('SF_subleadingEt_2016')
    EE_SF2016.Divide(efficiencyDY)     

    Efficiency_HLT_Full_pt_TT = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/MuonTrigger-SFs_AllPaths_2016Cuts/TH1F_muontrigger_2016.root', 'Efficiency_HLT_Full_pt_DYJetsToLL_M-50')
    Efficiency_HLT_Full_pt_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/MuonTrigger-SFs_AllPaths_2016Cuts/TH1F_muontrigger_2016.root', 'Efficiency_HLT_Full_pt_DATA')

    efficiencyTT = Efficiency_HLT_Full_pt_TT.GetPassedHistogram()
    totalTT = Efficiency_HLT_Full_pt_TT.GetTotalHistogram()
    efficiencyTT.Divide(totalTT)

    efficiencyData = Efficiency_HLT_Full_pt_DATA.GetPassedHistogram()
    totalData = Efficiency_HLT_Full_pt_DATA.GetTotalHistogram()
    efficiencyData.Divide(totalData)

    MM_SF2016 = efficiencyData.Clone('SF_subleadingPt_2016')
    MM_SF2016.Divide(efficiencyTT)     

    Efficiency_HLT_Full_pt_TT = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/MuonTrigger-SFs_AllPaths_2016Cuts/TH1F_muontrigger_2016APV.root', 'Efficiency_HLT_Full_pt_DYJetsToLL_M-50')
    Efficiency_HLT_Full_pt_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/MuonTrigger-SFs_AllPaths_2016Cuts/TH1F_muontrigger_2016APV.root', 'Efficiency_HLT_Full_pt_DATA')

    efficiencyTT = Efficiency_HLT_Full_pt_TT.GetPassedHistogram()
    totalTT = Efficiency_HLT_Full_pt_TT.GetTotalHistogram()
    efficiencyTT.Divide(totalTT)

    efficiencyData = Efficiency_HLT_Full_pt_DATA.GetPassedHistogram()
    totalData = Efficiency_HLT_Full_pt_DATA.GetTotalHistogram()
    efficiencyData.Divide(totalData)

    MM_SF2016APV = efficiencyData.Clone('SF_subleadingPt_2016APV')
    MM_SF2016APV.Divide(efficiencyTT)     

    ### Get plots (2017)
    Efficiency_HLT_Full_pt_Sim = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/FullComparisons/PhotonTrigger-SFs_2017Full_m100/TH1F_photontrigger_2017.root', 'Efficiency_HLT_Full_pt_TTTo2L2Nu')
    Efficiency_HLT_Full_pt_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/FullComparisons/PhotonTrigger-SFs_2017Full_m100/TH1F_photontrigger_2017.root', 'Efficiency_HLT_Full_pt_DATA')

    efficiencySim = Efficiency_HLT_Full_pt_Sim.GetPassedHistogram()
    totalSim = Efficiency_HLT_Full_pt_Sim.GetTotalHistogram()
    efficiencySim.Divide(totalSim)

    efficiencyData = Efficiency_HLT_Full_pt_DATA.GetPassedHistogram()
    totalData = Efficiency_HLT_Full_pt_DATA.GetTotalHistogram()
    efficiencyData.Divide(totalData)

    EE_SF2017 = efficiencyData.Clone('SF_subleadingEt_2017')
    EE_SF2017.Divide(efficiencySim)     


    ### Get plots (2018)
    Efficiency_HLT_Full_pt_DYJetsToLL_M50 = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/FullComparisons/PhotonTrigger-SFs_2018Full_DYOptimized/TH1F_photontrigger_2018.root', 'Efficiency_HLT_Full_pt_DYJetsToLL_M-50')
    Efficiency_HLT_Full_pt_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/FullComparisons/PhotonTrigger-SFs_2018Full_DYOptimized/TH1F_photontrigger_2018.root', 'Efficiency_HLT_Full_pt_DATA')

    efficiencyDY = Efficiency_HLT_Full_pt_DYJetsToLL_M50.GetPassedHistogram()
    totalDY = Efficiency_HLT_Full_pt_DYJetsToLL_M50.GetTotalHistogram()
    efficiencyDY.Divide(totalDY)

    efficiencyData = Efficiency_HLT_Full_pt_DATA.GetPassedHistogram()
    totalData = Efficiency_HLT_Full_pt_DATA.GetTotalHistogram()
    efficiencyData.Divide(totalData)

    EE_SF2018 = efficiencyData.Clone('SF_subleadingEt_2018')
    EE_SF2018.Divide(efficiencyDY)     

    Efficiency_HLT_Full_pt_TT = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/MuonTrigger-SFs_AllPaths_2018Cuts/TH1F_muontrigger_2018.root', 'Efficiency_HLT_Full_pt_DYJetsToLL_M-50')
    Efficiency_HLT_Full_pt_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/MuonTrigger-SFs_AllPaths_2018Cuts/TH1F_muontrigger_2018.root', 'Efficiency_HLT_Full_pt_DATA')

    efficiencyTT = Efficiency_HLT_Full_pt_TT.GetPassedHistogram()
    totalTT = Efficiency_HLT_Full_pt_TT.GetTotalHistogram()
    efficiencyTT.Divide(totalTT)

    efficiencyData = Efficiency_HLT_Full_pt_DATA.GetPassedHistogram()
    totalData = Efficiency_HLT_Full_pt_DATA.GetTotalHistogram()
    efficiencyData.Divide(totalData)

    MM_SF2018 = efficiencyData.Clone('SF_subleadingPt_2018')
    MM_SF2018.Divide(efficiencyTT)     

    ### Save plots
    file_ = r.TFile("ScaleFactors_Trigger_Spring23.root", "RECREATE")
    EE_SF2016.Write()
    EE_SF2017.Write()
    EE_SF2018.Write()
    MM_SF2016.Write()
    MM_SF2016APV.Write()
    MM_SF2018.Write()
    file_.Close()






