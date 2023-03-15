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


def lorentzianPeak(x, par):

    # par[0]
    term1 = 0.5*par[0]*par[1]/3.14
    term2 = (x[0]-par[2])*(x[0]-par[2]) + .25*par[1]*par[1]

    return term1/term2

def gaussian(x, par):
    
    return par[0]*r.TMath.Gaus(x[0], par[1], par[2], 0)


def signalplusbackground2(x, par):

    background = par[0] + par[1]*x[0] + par[2]*x[0]*x[0]
    signal = par[3]*r.TMath.Gaus(x[0], par[4], par[5], 0)

    return background + signal


def RooCMSShape(x, par):

  # par[0]: alpha
  # par[1]: beta
  # par[2]: gamma
  # par[3]: peak

    erf = r.TMath.Erfc((par[0] - x[0]) * par[1])
    u = (x[0] - par[3])*par[2]

    if(u < -70):
        u = 1e20
    elif( u>70 ):
        u = 0
    else: 
        u = math.exp(-u) #exponential decay
    return erf*u;

def exponential(x, par):

    arg = 0;
    if (par[2]):
        arg = (x[0] - par[1])/par[2]
    return par[0]*r.TMath.Exp(-0.5*arg*arg)*x[0]*x[0]

def crystalBall(x, par):

    return r.Math.crystalball_function(x[0], par[0], par[1], par[2], par[3])

def backgroundPDFError(x, par):

    if par[1] < 1e-10:
        term1 = 0.0
    else:
        term1 = math.exp(-x[0]/par[1])
    if par[3] < 1e-10:
        term2 = 0.
    else:
        term2 = (1+r.TMath.Erfc((x[0]-par[2])/ par[3]))

    return par[0]*term1*term2

def polynomialBkg(x, par):
    background = par[0] + par[1]*x[0] + par[2]*x[0]*x[0]
    return background

def CBPlusPDF(x, par):

    if par[1] < 1e-10:
        term1 = 0.0
    else:
        term1 = math.exp(-x[0]/par[1])
    if par[3] < 1e-10:
        term2 = 0.
    else:
        term2 = (1+r.TMath.Erfc((x[0]-par[2])/ par[3]))

    background =  par[0]*term1*term2
    signal = par[4]*r.Math.crystalball_function(x[0], par[5], par[6], par[7], par[8])

    return signal + background

def CBPlusPolynomial(x, par):

    background = par[0] + par[1]*x[0] + par[2]*x[0]*x[0]
    signal = par[3]*r.Math.crystalball_function(x[0], par[4], par[5], par[6], par[7])
    return background + signal

def CBPlusRooCMSShape(x, par):

  # par[0]: alpha
  # par[1]: beta
  # par[2]: gamma
  # par[3]: peak

    erf = r.TMath.Erfc((par[0] - x[0]) * par[1])
    u = (x[0] - par[3])*par[2]

    if(u < -70):
        u = 1e20
    elif( u>70 ):
        u = 0
    else: 
        u = math.exp(-u) #exponential decay

    background = erf*u;
    signal = par[4]*r.Math.crystalball_function(x[0], par[5], par[6], par[7], par[8])

    return background + signal

def CBPlusExp(x, par):

    arg = 0;
    if (par[2]):
        arg = (x[0] - par[1])/par[2]
    
    background = par[0]*r.TMath.Exp(-0.5*arg*arg)*x[0]*x[0]
    signal = par[3]*r.Math.crystalball_function(x[0], par[4], par[5], par[6], par[7])
    return background + signal

def CBPlusExpSimp(x, par):


    """
    arg = 0;
    if (par[2]):
        arg = (x[0] - par[1])/par[2]
    
    background = par[0]*r.TMath.Exp(-0.5*arg*arg)*x[0]*x[0]
    """
    # par[0]
    # par[1] : exp rate
    background = par[0]*math.exp(-0.5*x[0]/par[1])

    #signal = par[3]*r.TMath.Gaus(x[0], par[4], par[5], 0)
    # par[2]
    # par[3, 4, 5, 6] : alpha, n sigma, mu
    signal = par[2]*r.Math.crystalball_function(x[0], par[3], par[4], par[5], par[6])

    return background + signal

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

    ###########################
    ####   Parser object   ####
    ###########################
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    (opts, args) = parser.parse_args()

    ##############################
    ####   Some definitions   ####
    ##############################
    ptbin = np.concatenate((np.arange(0, 90, 5, float), np.arange(90, 130, 10), np.arange(130, 190, 15), np.arange(190, 260, 40), np.array([300])))

    #####################################
    ####   Construct TEfficiencies   ####
    #####################################
    
    ### Pt (IFCA)

    TH1F_probeIsoTrack_pt_total_MC = getObject('plots_condor/DYJetsToLL_M50_dispEle_2018.root', 'TH1F_probeIsoTrack_pt_total_DATA')
    TH1F_probeIsoTrack_pt_passed_MC = getObject('plots_condor/DYJetsToLL_M50_dispEle_2018.root', 'TH1F_probeIsoTrack_pt_passed_DATA')
    TEfficiency_probeIsoTrack_pt_MC = r.TEfficiency(TH1F_probeIsoTrack_pt_passed_MC, TH1F_probeIsoTrack_pt_total_MC)
    TH1F_probeIsoTrack_eta_total_MC = getObject('plots_condor/DYJetsToLL_M50_dispEle_2018.root', 'TH1F_probeIsoTrack_eta_total_DATA')
    TH1F_probeIsoTrack_eta_passed_MC = getObject('plots_condor/DYJetsToLL_M50_dispEle_2018.root', 'TH1F_probeIsoTrack_eta_passed_DATA')
    TEfficiency_probeIsoTrack_eta_MC = r.TEfficiency(TH1F_probeIsoTrack_eta_passed_MC, TH1F_probeIsoTrack_eta_total_MC)

    print(TH1F_probeIsoTrack_pt_total_MC.Integral(), TH1F_probeIsoTrack_pt_passed_MC.Integral())

    TH1F_probeIsoTrack_pt_total_DATA = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_probeIsoTrack_pt_total_DATA')
    TH1F_probeIsoTrack_pt_passed_DATA = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_probeIsoTrack_pt_passed_DATA')
    TEfficiency_probeIsoTrack_pt_DATA = r.TEfficiency(TH1F_probeIsoTrack_pt_passed_DATA, TH1F_probeIsoTrack_pt_total_DATA)
    TH1F_probeIsoTrack_eta_total_DATA = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_probeIsoTrack_eta_total_DATA')
    TH1F_probeIsoTrack_eta_passed_DATA = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_probeIsoTrack_eta_passed_DATA')
    TEfficiency_probeIsoTrack_eta_DATA = r.TEfficiency(TH1F_probeIsoTrack_eta_passed_DATA, TH1F_probeIsoTrack_eta_total_DATA)

    EFF_pt = Canvas.Canvas("TEfficiency_probeIsoTrack_pt", 'png', 0.65, 0.77, 0.85, 0.9, 1) 
    EFF_pt.addRate(TEfficiency_probeIsoTrack_pt_MC, 'AP', 'Simulation', 'p', r.kBlue, True, 0, marker = 24)
    EFF_pt.addRate(TEfficiency_probeIsoTrack_pt_DATA, 'P,SAME', 'DATA', 'p', r.kBlack, True, 0, marker = 24)
    EFF_pt.addLatex(0.9, 0.93, 'On-Z Control Region (Tag-and-probe selection)', size = 0.032, align = 31)
    EFF_pt.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_TrackSC/', inProgress = True)
    
    EFF_eta = Canvas.Canvas("TEfficiency_probeIsoTrack_eta", 'png', 0.65, 0.77, 0.85, 0.9, 1) 
    EFF_eta.addRate(TEfficiency_probeIsoTrack_eta_MC, 'AP', 'Simulation', 'p', r.kBlue, True, 0, marker = 24)
    EFF_eta.addRate(TEfficiency_probeIsoTrack_eta_DATA, 'P,SAME', 'Data', 'p', r.kBlack, True, 0, marker = 24)
    EFF_eta.addLatex(0.9, 0.93, 'On-Z Control Region (Tag-and-probe selection)', size = 0.032, align = 31)
    EFF_eta.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_TrackSC/', inProgress = True)

    ### Mass distribution of failed and passed
    TH1F_tnp_mass_passed = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_tnp_mass_passed_DATA')
    #TH1F_tnp_mass_passed_C = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_tnp_mass_passed_DATA')
    #TH1F_tnp_mass_passed_D = getObject('plots_condor/EGammaRun2018D.root', 'TH1F_tnp_mass_passed_DATA')
    TH1F_tnp_mass_failed = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_tnp_mass_failed_DATA')
    #TH1F_tnp_mass_failed_C = getObject('plots_condor/EGammaRun2018C.root', 'TH1F_tnp_mass_failed_DATA')
    #TH1F_tnp_mass_failed_D = getObject('plots_condor/EGammaRun2018D.root', 'TH1F_tnp_mass_failed_DATA')
    #TH1F_tnp_mass_passed.Add(TH1F_tnp_mass_passed_C)
    #TH1F_tnp_mass_passed.Add(TH1F_tnp_mass_passed_D)
    #TH1F_tnp_mass_failed.Add(TH1F_tnp_mass_failed_C)
    #TH1F_tnp_mass_failed.Add(TH1F_tnp_mass_failed_D)
    TH1F_tnp_mass_failed.SetLineColor(r.kGray+2)
    TH1F_tnp_mass_failed.SetFillColor(r.kGray)
    TH1F_tnp_mass_passed.SetLineColor(r.kCyan-2)
    TH1F_tnp_mass_passed.SetFillColor(r.kCyan-6)
    TH1F_tnp_mass_passed.SetTitle('Passing probes')
    TH1F_tnp_mass_failed.SetTitle('Failing probes')
    Stack_mass = r.THStack("aux_stack", ";Tnp mass (GeV);Tnp pairs")
    Stack_mass.Add(TH1F_tnp_mass_failed) 
    Stack_mass.Add(TH1F_tnp_mass_passed) 
    StackPlot = Canvas.Canvas("StackPlot", 'png', 0.6, 0.77, 0.8, 0.9, 1)
    StackPlot.addStack(Stack_mass, 'HIST', 1, 0)
    StackPlot.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_TrackSC/', inProgress = True)

