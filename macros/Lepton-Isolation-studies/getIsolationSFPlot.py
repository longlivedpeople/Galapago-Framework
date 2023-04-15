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


if __name__ == "__main__":


    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    print('WORKPATH: ' + WORKPATH)
    r.setTDRStyle()

    ### Get plots (2016)
    h2016_EE_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_EE_Electron_et_total_DATA')
    h2016_EE_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_EE_Electron_et_passedIso_DATA')
    h2016_EE_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_EE_Electron_et_total_MC')
    h2016_EE_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_EE_Electron_et_passedIso_MC')
    h2016_DMDM_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_DMDM_DGM_pt_total_DATA')
    h2016_DMDM_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_DMDM_DGM_pt_passedIso_DATA')
    h2016_DMDM_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_DMDM_DGM_pt_total_MC')
    h2016_DMDM_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016/TH1F_Isolation_2016.root', 'TH1F_DMDM_DGM_pt_passedIso_MC')
 
    ### Get plots (2016APV)
    h2016APV_EE_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_EE_Electron_et_total_DATA')
    h2016APV_EE_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_EE_Electron_et_passedIso_DATA')
    h2016APV_EE_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_EE_Electron_et_total_MC')
    h2016APV_EE_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_EE_Electron_et_passedIso_MC')
    h2016APV_DMDM_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_DMDM_DGM_pt_total_DATA')
    h2016APV_DMDM_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_DMDM_DGM_pt_passedIso_DATA')
    h2016APV_DMDM_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_DMDM_DGM_pt_total_MC')
    h2016APV_DMDM_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2016APV/TH1F_Isolation_2016APV.root', 'TH1F_DMDM_DGM_pt_passedIso_MC')
 
    ### Get plots (2017)
    h2017_EE_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2017/TH1F_Isolation_2017.root', 'TH1F_EE_Electron_et_total_DATA')
    h2017_EE_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2017/TH1F_Isolation_2017.root', 'TH1F_EE_Electron_et_passedIso_DATA')
    h2017_EE_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2017/TH1F_Isolation_2017.root', 'TH1F_EE_Electron_et_total_MC')
    h2017_EE_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/Isolation-SFs_2017/TH1F_Isolation_2017.root', 'TH1F_EE_Electron_et_passedIso_MC')
 
    ### Get plots (2018)
    h2018_EE_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_EE_Electron_et_total_DATA')
    h2018_EE_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_EE_Electron_et_passedIso_DATA')
    h2018_EE_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_EE_Electron_et_total_MC')
    h2018_EE_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_EE_Electron_et_passedIso_MC')
    h2018_DMDM_total_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_DMDM_DGM_pt_total_DATA')
    h2018_DMDM_passed_DATA = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_DMDM_DGM_pt_passedIso_DATA')
    h2018_DMDM_total_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_DMDM_DGM_pt_total_MC')
    h2018_DMDM_passed_MC = getObject('/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/Isolation-SFs_2018/TH1F_Isolation_2018.root', 'TH1F_DMDM_DGM_pt_passedIso_MC')
 

    ### Get SF plot (2016)
    h2016_EE_eff_DATA = h2016_EE_passed_DATA.Clone()
    h2016_EE_eff_DATA.Divide(h2016_EE_total_DATA)
    h2016_EE_eff_MC = h2016_EE_passed_MC.Clone()
    h2016_EE_eff_MC.Divide(h2016_EE_total_MC)
    h2016_EE_SF = h2016_EE_eff_DATA.Clone('SF_Iso_EE_2016')
    h2016_EE_SF.Divide(h2016_EE_eff_MC)
    h2016_DMDM_eff_DATA = h2016_DMDM_passed_DATA.Clone()
    h2016_DMDM_eff_DATA.Divide(h2016_DMDM_total_DATA)
    h2016_DMDM_eff_MC = h2016_DMDM_passed_MC.Clone()
    h2016_DMDM_eff_MC.Divide(h2016_DMDM_total_MC)
    h2016_DMDM_SF = h2016_DMDM_eff_DATA.Clone('SF_Iso_DMDM_2016')
    h2016_DMDM_SF.Divide(h2016_DMDM_eff_MC)

    ### Get SF plot (2016APV)
    h2016APV_EE_eff_DATA = h2016APV_EE_passed_DATA.Clone()
    h2016APV_EE_eff_DATA.Divide(h2016APV_EE_total_DATA)
    h2016APV_EE_eff_MC = h2016APV_EE_passed_MC.Clone()
    h2016APV_EE_eff_MC.Divide(h2016APV_EE_total_MC)
    h2016APV_EE_SF = h2016APV_EE_eff_DATA.Clone('SF_Iso_EE_2016APV')
    h2016APV_EE_SF.Divide(h2016APV_EE_eff_MC)
    h2016APV_DMDM_eff_DATA = h2016APV_DMDM_passed_DATA.Clone()
    h2016APV_DMDM_eff_DATA.Divide(h2016APV_DMDM_total_DATA)
    h2016APV_DMDM_eff_MC = h2016APV_DMDM_passed_MC.Clone()
    h2016APV_DMDM_eff_MC.Divide(h2016APV_DMDM_total_MC)
    h2016APV_DMDM_SF = h2016APV_DMDM_eff_DATA.Clone('SF_Iso_DMDM_2016APV')
    h2016APV_DMDM_SF.Divide(h2016APV_DMDM_eff_MC)

    ### Get SF plot (2017)
    h2017_EE_eff_DATA = h2017_EE_passed_DATA.Clone()
    h2017_EE_eff_DATA.Divide(h2017_EE_total_DATA)
    h2017_EE_eff_MC = h2017_EE_passed_MC.Clone()
    h2017_EE_eff_MC.Divide(h2017_EE_total_MC)
    h2017_EE_SF = h2017_EE_eff_DATA.Clone('SF_Iso_EE_2017')
    h2017_EE_SF.Divide(h2017_EE_eff_MC)

    ### Get SF plot (2018)
    h2018_EE_eff_DATA = h2018_EE_passed_DATA.Clone()
    h2018_EE_eff_DATA.Divide(h2018_EE_total_DATA)
    h2018_EE_eff_MC = h2018_EE_passed_MC.Clone()
    h2018_EE_eff_MC.Divide(h2018_EE_total_MC)
    h2018_EE_SF = h2018_EE_eff_DATA.Clone('SF_Iso_EE_2018')
    h2018_EE_SF.Divide(h2018_EE_eff_MC)
    h2018_DMDM_eff_DATA = h2018_DMDM_passed_DATA.Clone()
    h2018_DMDM_eff_DATA.Divide(h2018_DMDM_total_DATA)
    h2018_DMDM_eff_MC = h2018_DMDM_passed_MC.Clone()
    h2018_DMDM_eff_MC.Divide(h2018_DMDM_total_MC)
    h2018_DMDM_SF = h2018_DMDM_eff_DATA.Clone('SF_Iso_DMDM_2018')
    h2018_DMDM_SF.Divide(h2018_DMDM_eff_MC)


    outputFile = TFile('Results_sf.root', 'RECREATE') 
    h2016_EE_SF.Write()
    h2016APV_EE_SF.Write()
    h2017_EE_SF.Write()
    h2018_EE_SF.Write()
    h2016_DMDM_SF.Write()
    h2016APV_DMDM_SF.Write()
    h2018_DMDM_SF.Write()
    outputFile.Close()




