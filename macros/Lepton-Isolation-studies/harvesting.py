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

    ###########################
    ####   Parser object   ####
    ###########################
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    (opts, args) = parser.parse_args()

    ##############################
    ####   Some definitions   ####
    ##############################

    #####################################
    ####   Construct TEfficiencies   ####
    #####################################
    
    hist_bkg = getObject('Results/th1f.root', 'hist_E_relTrkIso_log')
    hist_bkg.Scale(1./hist_bkg.Integral())
    hist_bkg.SetLineColor(r.kGray+1)
    hist_bkg.SetFillColor(r.kGray)
    hist_bkg.GetXaxis().SetTitle('Rel. Tracker isolation')
    hist_bkg.GetYaxis().SetTitleOffset(1.5)
    hist_bkg.GetYaxis().SetTitle('Normalized disp. Electron yield')
    hist_bkg.SetMaximum(2*hist_bkg.GetMaximum())
    hist_125_30_1 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_125_30_1_2016')
    hist_125_30_1.Scale(1./hist_125_30_1.Integral())
    hist_125_30_1.SetLineColor(r.kGreen+2)
    hist_125_30_1.SetLineWidth(2)
    hist_125_30_10 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_125_30_10_2016')
    hist_125_30_10.Scale(1./hist_125_30_10.Integral())
    hist_125_30_10.SetLineColor(r.kGreen+3)
    hist_125_30_10.SetLineWidth(2)
    hist_125_30_100 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_125_30_100_2016')
    hist_125_30_100.Scale(1./hist_125_30_100.Integral())
    hist_125_30_100.SetLineColor(r.kGreen-7)
    hist_125_30_100.SetLineWidth(2)
    hist_400_50_1 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_400_50_1_2016')
    hist_400_50_1.Scale(1./hist_400_50_1.Integral())
    hist_400_50_1.SetLineColor(r.kRed+1)
    hist_400_50_1.SetLineWidth(2)
    hist_400_50_10 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_400_50_10_2016')
    hist_400_50_10.Scale(1./hist_400_50_10.Integral())
    hist_400_50_10.SetLineColor(r.kRed-4)
    hist_400_50_10.SetLineWidth(2)
    hist_400_50_100 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_400_50_100_2016')
    hist_400_50_100.Scale(1./hist_400_50_100.Integral())
    hist_400_50_100.SetLineColor(r.kRed-7)
    hist_400_50_100.SetLineWidth(2)
    hist_1000_150_1 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_1000_150_1_2016')
    hist_1000_150_1.Scale(1./hist_1000_150_1.Integral())
    hist_1000_150_1.SetLineColor(r.kBlue+1)
    hist_1000_150_1.SetLineWidth(2)
    hist_1000_150_10 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_1000_150_10_2016')
    hist_1000_150_10.Scale(1./hist_1000_150_10.Integral())
    hist_1000_150_10.SetLineColor(r.kBlue-4)
    hist_1000_150_10.SetLineWidth(2)
    hist_1000_150_100 = getObject('Results/th1f_signals.root', 'hist_E_relTrkIso_log_HSS_1000_150_100_2016')
    hist_1000_150_100.Scale(1./hist_1000_150_100.Integral())
    hist_1000_150_100.SetLineColor(r.kBlue-7)
    hist_1000_150_100.SetLineWidth(2)
    plot = Canvas.Canvas('hist_electron_relTrkIso', 'png,pdf', 0.4, 0.55, 0.8, 0.89, 1, lsize=0.025)
    plot.addHisto(hist_bkg, 'HIST', 'Background: Hadronic WZ, WW decays', 'f', '', 1, 0)
    plot.addHisto(hist_125_30_1, 'HIST,SAME', '125 GeV, 30 GeV, 1 mm', 'l', '', 1, 1)
    plot.addHisto(hist_125_30_10, 'HIST,SAME', '125 GeV, 30 GeV, 10 mm', 'l', '', 1, 2)
    plot.addHisto(hist_125_30_100, 'HIST,SAME', '125 GeV, 30 GeV, 100 mm', 'l', '', 1, 3)
    plot.addHisto(hist_400_50_1, 'HIST,SAME', '400 GeV, 50 GeV, 1 mm', 'l', '', 1, 4)
    plot.addHisto(hist_400_50_10, 'HIST,SAME', '400 GeV, 50 GeV, 10 mm', 'l', '', 1, 5)
    plot.addHisto(hist_400_50_100, 'HIST,SAME', '400 GeV, 50 GeV, 100 mm', 'l', '', 1, 6)
    plot.addHisto(hist_1000_150_1, 'HIST,SAME', '1000 GeV, 150 GeV, 1 mm', 'l', '', 1, 7)
    plot.addHisto(hist_1000_150_10, 'HIST,SAME', '1000 GeV, 150 GeV, 10 mm', 'l', '', 1, 8)
    plot.addHisto(hist_1000_150_100, 'HIST,SAME', '1000 GeV, 150 GeV, 100 mm', 'l', '', 1, 9)
    plot.addLatex(0.9, 0.93, '2016 UL', size = 0.03, align = 31)
    plot.addLatex(0.75, 0.5, 'disp. Electrons:', size = 0.03, align = 22)
    plot.addLatex(0.75, 0.44, 'E_{T} > 25 GeV', size = 0.03, align = 22)
    plot.addLatex(0.75, 0.4, '|#eta| < 2', size = 0.03, align = 22)
    plot.addLine(0.1, hist_bkg.GetMinimum(), 0.1, 0.09, r.kBlack, 2)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested/', isPrivate = True, xlog = True, is2d = True)


    hist_bkg = getObject('Results/th1f.root', 'hist_E_relPFIso_log')
    hist_bkg.Scale(1./hist_bkg.Integral())
    hist_bkg.SetLineColor(r.kGray+1)
    hist_bkg.SetFillColor(r.kGray)
    hist_bkg.GetXaxis().SetTitle('Rel. PF isolation')
    hist_bkg.GetYaxis().SetTitleOffset(1.5)
    hist_bkg.GetYaxis().SetTitle('disp. Electron yield')
    hist_bkg.SetMaximum(2*hist_bkg.GetMaximum())
    hist_125_30_1 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_125_30_1_2016')
    hist_125_30_1.Scale(1./hist_125_30_1.Integral())
    hist_125_30_1.SetLineColor(r.kGreen+2)
    hist_125_30_1.SetLineWidth(2)
    hist_125_30_10 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_125_30_10_2016')
    hist_125_30_10.Scale(1./hist_125_30_10.Integral())
    hist_125_30_10.SetLineColor(r.kGreen+3)
    hist_125_30_10.SetLineWidth(2)
    hist_125_30_100 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_125_30_100_2016')
    hist_125_30_100.Scale(1./hist_125_30_100.Integral())
    hist_125_30_100.SetLineColor(r.kGreen-7)
    hist_125_30_100.SetLineWidth(2)
    hist_400_50_1 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_400_50_1_2016')
    hist_400_50_1.Scale(1./hist_400_50_1.Integral())
    hist_400_50_1.SetLineColor(r.kRed+1)
    hist_400_50_1.SetLineWidth(2)
    hist_400_50_10 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_400_50_10_2016')
    hist_400_50_10.Scale(1./hist_400_50_10.Integral())
    hist_400_50_10.SetLineColor(r.kRed-4)
    hist_400_50_10.SetLineWidth(2)
    hist_400_50_100 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_400_50_100_2016')
    hist_400_50_100.Scale(1./hist_400_50_100.Integral())
    hist_400_50_100.SetLineColor(r.kRed-7)
    hist_400_50_100.SetLineWidth(2)
    hist_1000_150_1 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_1000_150_1_2016')
    hist_1000_150_1.Scale(1./hist_1000_150_1.Integral())
    hist_1000_150_1.SetLineColor(r.kBlue+1)
    hist_1000_150_1.SetLineWidth(2)
    hist_1000_150_10 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_1000_150_10_2016')
    hist_1000_150_10.Scale(1./hist_1000_150_10.Integral())
    hist_1000_150_10.SetLineColor(r.kBlue-4)
    hist_1000_150_10.SetLineWidth(2)
    hist_1000_150_100 = getObject('Results/th1f_signals.root', 'hist_E_relPFIso_log_HSS_1000_150_100_2016')
    hist_1000_150_100.Scale(1./hist_1000_150_100.Integral())
    hist_1000_150_100.SetLineColor(r.kBlue-7)
    hist_1000_150_100.SetLineWidth(2)
    plot = Canvas.Canvas('hist_electron_relPFIso', 'png,pdf', 0.4, 0.55, 0.8, 0.89, 1, lsize=0.025)
    plot.addHisto(hist_bkg, 'HIST', 'Background: Hadronic WZ, WW decays', 'f', '', 1, 0)
    plot.addHisto(hist_125_30_1, 'HIST,SAME', '125 GeV, 30 GeV, 1 mm', 'l', '', 1, 1)
    plot.addHisto(hist_125_30_10, 'HIST,SAME', '125 GeV, 30 GeV, 10 mm', 'l', '', 1, 2)
    plot.addHisto(hist_125_30_100, 'HIST,SAME', '125 GeV, 30 GeV, 100 mm', 'l', '', 1, 3)
    plot.addHisto(hist_400_50_1, 'HIST,SAME', '400 GeV, 50 GeV, 1 mm', 'l', '', 1, 4)
    plot.addHisto(hist_400_50_10, 'HIST,SAME', '400 GeV, 50 GeV, 10 mm', 'l', '', 1, 5)
    plot.addHisto(hist_400_50_100, 'HIST,SAME', '400 GeV, 50 GeV, 100 mm', 'l', '', 1, 6)
    plot.addHisto(hist_1000_150_1, 'HIST,SAME', '1000 GeV, 150 GeV, 1 mm', 'l', '', 1, 7)
    plot.addHisto(hist_1000_150_10, 'HIST,SAME', '1000 GeV, 150 GeV, 10 mm', 'l', '', 1, 8)
    plot.addHisto(hist_1000_150_100, 'HIST,SAME', '1000 GeV, 150 GeV, 100 mm', 'l', '', 1, 9)
    plot.addLatex(0.9, 0.93, '2016 UL', size = 0.03, align = 31)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested/', isPrivate = False, xlog = True, is2d = True)


    hist_bkg = getObject('Results/th1f.root', 'hist_M_relPFIso_log')
    hist_bkg.Scale(1./hist_bkg.Integral())
    hist_bkg.SetLineColor(r.kGray+1)
    hist_bkg.SetFillColor(r.kGray)
    hist_bkg.GetXaxis().SetTitle('Rel. PF isolation')
    hist_bkg.GetYaxis().SetTitleOffset(1.5)
    hist_bkg.GetYaxis().SetTitle('Normalized dGlobal muon yield')
    hist_bkg.SetMaximum(2*hist_bkg.GetMaximum())
    hist_125_30_1 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_125_30_1_2016')
    hist_125_30_1.Scale(1./hist_125_30_1.Integral())
    hist_125_30_1.SetLineColor(r.kGreen+2)
    hist_125_30_1.SetLineWidth(2)
    hist_125_30_10 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_125_30_10_2016')
    hist_125_30_10.Scale(1./hist_125_30_10.Integral())
    hist_125_30_10.SetLineColor(r.kGreen+3)
    hist_125_30_10.SetLineWidth(2)
    hist_125_30_100 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_125_30_100_2016')
    hist_125_30_100.Scale(1./hist_125_30_100.Integral())
    hist_125_30_100.SetLineColor(r.kGreen-7)
    hist_125_30_100.SetLineWidth(2)
    hist_400_50_1 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_400_50_1_2016')
    hist_400_50_1.Scale(1./hist_400_50_1.Integral())
    hist_400_50_1.SetLineColor(r.kRed+1)
    hist_400_50_1.SetLineWidth(2)
    hist_400_50_10 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_400_50_10_2016')
    hist_400_50_10.Scale(1./hist_400_50_10.Integral())
    hist_400_50_10.SetLineColor(r.kRed-4)
    hist_400_50_10.SetLineWidth(2)
    hist_400_50_100 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_400_50_100_2016')
    hist_400_50_100.Scale(1./hist_400_50_100.Integral())
    hist_400_50_100.SetLineColor(r.kRed-7)
    hist_400_50_100.SetLineWidth(2)
    hist_1000_150_1 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_1000_150_1_2016')
    hist_1000_150_1.Scale(1./hist_1000_150_1.Integral())
    hist_1000_150_1.SetLineColor(r.kBlue+1)
    hist_1000_150_1.SetLineWidth(2)
    hist_1000_150_10 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_1000_150_10_2016')
    hist_1000_150_10.Scale(1./hist_1000_150_10.Integral())
    hist_1000_150_10.SetLineColor(r.kBlue-4)
    hist_1000_150_10.SetLineWidth(2)
    hist_1000_150_100 = getObject('Results/th1f_signals.root', 'hist_M_relPFIso_log_HSS_1000_150_100_2016')
    hist_1000_150_100.Scale(1./hist_1000_150_100.Integral())
    hist_1000_150_100.SetLineColor(r.kBlue-7)
    hist_1000_150_100.SetLineWidth(2)
    plot = Canvas.Canvas('hist_muon_relPFIso', 'png,pdf', 0.4, 0.55, 0.8, 0.89, 1, lsize=0.025)
    plot.addHisto(hist_bkg, 'HIST', 'Background: Hadronic WZ, WW decays', 'f', '', 1, 0)
    plot.addHisto(hist_125_30_1, 'HIST,SAME', '125 GeV, 30 GeV, 1 mm', 'l', '', 1, 1)
    plot.addHisto(hist_125_30_10, 'HIST,SAME', '125 GeV, 30 GeV, 10 mm', 'l', '', 1, 2)
    plot.addHisto(hist_125_30_100, 'HIST,SAME', '125 GeV, 30 GeV, 100 mm', 'l', '', 1, 3)
    plot.addHisto(hist_400_50_1, 'HIST,SAME', '400 GeV, 50 GeV, 1 mm', 'l', '', 1, 4)
    plot.addHisto(hist_400_50_10, 'HIST,SAME', '400 GeV, 50 GeV, 10 mm', 'l', '', 1, 5)
    plot.addHisto(hist_400_50_100, 'HIST,SAME', '400 GeV, 50 GeV, 100 mm', 'l', '', 1, 6)
    plot.addHisto(hist_1000_150_1, 'HIST,SAME', '1000 GeV, 150 GeV, 1 mm', 'l', '', 1, 7)
    plot.addHisto(hist_1000_150_10, 'HIST,SAME', '1000 GeV, 150 GeV, 10 mm', 'l', '', 1, 8)
    plot.addHisto(hist_1000_150_100, 'HIST,SAME', '1000 GeV, 150 GeV, 100 mm', 'l', '', 1, 9)
    plot.addLatex(0.9, 0.93, '2016 UL', size = 0.03, align = 31)
    plot.addLatex(0.75, 0.5, 'dGlobal muons:', size = 0.03, align = 22)
    plot.addLatex(0.75, 0.44, 'p_{T} > 30 GeV', size = 0.03, align = 22)
    plot.addLatex(0.75, 0.4, '|#eta| < 2', size = 0.03, align = 22)
    plot.addLine(0.2, hist_bkg.GetMinimum(), 0.2, 0.11, r.kBlack, 2)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested/', isPrivate = True, xlog = True, is2d = True)





