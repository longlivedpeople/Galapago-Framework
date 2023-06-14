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
    
    hSR  = getObject('Results/th1f.root', 'hSR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90')
    hBCR = getObject('Results/th1f.root', 'hBCR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90')
    hSR.Sumw2()
    hBCR.Sumw2()

    plot = Canvas.Canvas('hist_fullPath', 'png,pdf', 0.17, 0.55, 0.45, 0.89, 1, lsize=0.028)
    plot.addHisto(hSR, 'P', 'SR (fullPath)', 'p', r.kBlue, 1, 0)
    plot.addHisto(hBCR, 'P,SAME', 'BCR (fullPath)', 'p', r.kRed, 1, 1)
    plot.addLatex(0.9, 0.93, '2016 UL', size = 0.03, align = 31)
    plot.addLatex(0.75, 0.5, 'disp. Electrons:', size = 0.03, align = 22)
    plot.addLatex(0.75, 0.44, 'E_{T} > 25 GeV', size = 0.03, align = 22)
    plot.addLatex(0.75, 0.4, '|#eta| < 2', size = 0.03, align = 22)
    plot.saveRatio(1, 1, 1, '', hSR, hBCR, r_ymin = 0.5, r_ymax = 1.5, outputDir = WORKPATH + 'harvested/')

    h_fullPath = []
    h_m90 = []
    h_m55 = []
    h_p70 = []

    h_fullPath.append(getObject('Results/th1f.root', 'hSR_fullPath'))
    h_fullPath.append(getObject('Results/th1f.root', 'hBCR_fullPath'))
    h_m90.append(getObject('Results/th1f.root', 'hSR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90'))
    h_m90.append(getObject('Results/th1f.root', 'hBCR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90'))
    h_m55.append(getObject('Results/th1f.root', 'hSR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55'))
    h_m55.append(getObject('Results/th1f.root', 'hBCR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55'))
    h_p70.append(getObject('Results/th1f.root', 'hSR_fullPath'))
    h_p70.append(getObject('Results/th1f.root', 'hBCR_fullPath'))

    h_all = [h_fullPath, h_m90, h_m55, h_p70]

    ### Plot
    canvas = r.TCanvas("c1", "", 600, 600)
    canvas.cd()
    pad1 = r.TPad("pad1", "pad1", 0, 0.3, 1, 0.99)
    pad1.SetBottomMargin(0.015)
    pad1.SetTopMargin(0.13)
    pad1.Draw()
    #pad1.SetLogx(1)
    #r.gStyle.SetOptStat(0)
    pad2 = r.TPad("pad2", "pad2", 0, 0.01, 1, 0.3)
    pad2.SetTopMargin(0.05);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    #pad2.SetLogx(1)

    h_fullPath[0].SetMarkerColor(r.kBlack)
    h_fullPath[1].SetMarkerColor(r.kBlack)
    h_fullPath[0].SetLineColor(r.kBlack)
    h_fullPath[1].SetLineColor(r.kBlack)
    h_fullPath[0].SetMarkerStyle(20)
    h_fullPath[1].SetMarkerStyle(24)
    h_fullPath[0].SetMarkerSize(1)
    h_fullPath[1].SetMarkerSize(1)
    h_fullPath[0].Sumw2()
    h_fullPath[1].Sumw2()
    ratio_fullPath = h_fullPath[0].Clone()
    ratio_fullPath.Divide(h_fullPath[1])

    h_m90[0].SetMarkerColor(r.TColor.GetColor('#048ba8'))
    h_m90[1].SetMarkerColor(r.TColor.GetColor('#048ba8'))
    h_m90[0].SetLineColor(r.TColor.GetColor('#048ba8'))
    h_m90[1].SetLineColor(r.TColor.GetColor('#048ba8'))
    h_m90[0].SetMarkerStyle(21)
    h_m90[1].SetMarkerStyle(25)
    h_m90[0].SetMarkerSize(1)
    h_m90[1].SetMarkerSize(1)
    h_m90[0].Sumw2()
    h_m90[1].Sumw2()
    ratio_m90 = h_m90[0].Clone()
    ratio_m90.Divide(h_m90[1])
    ratio_m90.SetMarkerStyle(20)

    h_m55[0].SetMarkerColor(r.kAzure+2)
    h_m55[1].SetMarkerColor(r.kAzure+2)
    h_m55[0].SetLineColor(r.kAzure+2)
    h_m55[1].SetLineColor(r.kAzure+2)
    h_m55[0].SetMarkerStyle(22)
    h_m55[1].SetMarkerStyle(26)
    h_m55[0].SetMarkerSize(1)
    h_m55[1].SetMarkerSize(1)
    h_m55[0].Sumw2()
    h_m55[1].Sumw2()
    ratio_m55 = h_m55[0].Clone()
    ratio_m55.Divide(h_m55[1])
    ratio_m55.SetMarkerStyle(20)

    h_p70[0].SetMarkerColor(r.TColor.GetColor('#f18f01'))
    h_p70[1].SetMarkerColor(r.TColor.GetColor('#f18f01'))
    h_p70[0].SetLineColor(r.TColor.GetColor('#f18f01'))
    h_p70[1].SetLineColor(r.TColor.GetColor('#f18f01'))
    h_p70[0].SetMarkerStyle(23)
    h_p70[1].SetMarkerStyle(32)
    h_p70[0].SetMarkerSize(1)
    h_p70[1].SetMarkerSize(1)
    h_p70[0].Sumw2()
    h_p70[1].Sumw2()
    ratio_p70 = h_p70[0].Clone()
    ratio_p70.Divide(h_p70[1])
    ratio_p70.SetMarkerStyle(20)

    pad1.cd()
    pad1.SetLogy(1)
    h_fullPath[0].SetMaximum(1e10)
    h_fullPath[0].SetMinimum(1)
    h_fullPath[0].GetXaxis().SetLabelSize(0)
    h_fullPath[0].Draw('P E0')
    h_fullPath[1].Draw('P E0, SAME')
    h_m90[0].Draw('P E0,SAME')
    h_m90[1].Draw('P E0,SAME')
    #h_m55[0].Draw('P E0 E1,SAME')
    #h_m55[1].Draw('P E0 E1,SAME')
    h_p70[0].Draw('P E0,SAME')
    h_p70[1].Draw('P E0,SAME')

    legend = r.TLegend(0.18, 0.50, 0.3, 0.82)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0) # added (Celia)
    legend.AddEntry(h_fullPath[0], 'Full path (Forward)', 'p')
    legend.AddEntry(h_fullPath[1], 'Full path (Backward)', 'p')
    legend.AddEntry(h_m90[0], 'HLT_Diphoton30_22_*_Mass90 (Forward)', 'p')
    legend.AddEntry(h_m90[1], 'HLT_Diphoton30_22_*_Mass90 (Backward)', 'p')
    legend.AddEntry(h_p70[0], 'HLT_DoublePhoton70 (Forward)', 'p')
    legend.AddEntry(h_p70[1], 'HLT_DoublePhoton70 (Backward)', 'p')
    legend.Draw()

    latex = r.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(r.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.13, 0.88, "#bf{Private work}")
    latexb = r.TLatex()
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextColor(r.kBlack);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(31);
    latexb.SetTextSize(0.045);
    latexb.DrawLatex(0.6, 0.88, "#it{CMS data/simulation}")

    latexc = r.TLatex()
    latexc.SetNDC();
    latexc.SetTextAngle(0);
    latexc.SetTextColor(r.kBlack);
    latexc.SetTextFont(42);
    latexc.SetTextAlign(31);
    latexc.SetTextSize(0.05);
    latexc.DrawLatex(0.90, 0.88, '41.5 fb^{-1} (13 TeV)')

    

    pad2.cd()
    ratio_fullPath.SetMaximum(1.4)
    ratio_fullPath.SetMinimum(0.6)
    ratio_fullPath.GetXaxis().SetTitleSize(0.11)
    ratio_fullPath.GetXaxis().SetLabelSize(0.10)
    ratio_fullPath.GetYaxis().SetLabelSize(0.10)
    ratio_fullPath.GetYaxis().SetTitleOffset(0.43)
    ratio_fullPath.GetYaxis().SetTitleSize(0.12)
    ratio_fullPath.GetYaxis().SetTitle('Obs./Pred.')
    ratio_fullPath.GetYaxis().SetNdivisions(4)
    ratio_fullPath.SetLineWidth(2)
    ratio_m90.SetLineWidth(2)
    ratio_p70.SetLineWidth(2)
    ratio_fullPath.Draw('P E0')
    ratio_m90.Draw('P E0,SAME')
    ratio_p70.Draw('P E0,SAME')

    canvas.Print('EE_mass_2017_triggersplit.png')
    canvas.Print('EE_mass_2017_triggersplit.pdf')




