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


def makeAgreementTest(lumi, histo1, histo2, ylog, label1, label2, labela, labelb, name, isData, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', rebin = 0, rmin = 0.9, rmax = 1.1, maxY = False):


    ## Histogram tunning
    if rebin:
        histo1.Rebin(rebin)
        histo2.Rebin(rebin)

    histo1.SetTitle('; Collinearity |#Delta#Phi|; Dimuon pairs')


    histo1.SetMaximum(1.45*histo1.GetMaximum())
    histo1.SetMinimum(0.0)
    histo1.SetMarkerStyle(20)
    histo2.SetMarkerStyle(24)

    plot = Canvas.Canvas(name, 'png,pdf', 0.15, 0.66, 0.45, 0.78, 1)
    plot.addHisto(histo1, 'P', label1, 'p', r.kBlack, 1, 0)
    plot.addHisto(histo2, 'P,SAME', label2, 'p', r.kBlue, 1, 1)
    plot.addLatex(0.16, 0.8, labela, font = 42, size = 0.037, align = 11)
    plot.addLatex(0.88, 0.88, labelb, font = 42, size = 0.045, align = 31)
    plot.saveRatio(1, 0, 0, '', histo1, histo2, r_ymin = rmin, r_ymax = rmax, label = 'Ratio',  xlog = False, outputDir = WORKPATH + 'SymmetryResults/', maxYnumbers = maxY)

    print('>>>>>>>>>> KOLMOGOROV test for ' + labela)
    print('>>>>>>>>>> ' + str(histo1.KolmogorovTest(histo2)))
    print('>>>>>>>>>> Chi2 test for ' + labela)
    print('>>>>>>>>>> ' + str(histo1.Chi2Test(histo2, "WWP")))

    return




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
    

    hist_1d_r_chi = getObject('Results/th1f_2018ABCD.root', 'hist_1d_r_chi')
    hist_1d_r_char  = getObject('Results/th1f_2018ABCD.root', 'hist_1d_r_char')
    hist_1d_r_cos = getObject('Results/th1f_2018ABCD.root', 'hist_1d_r_cos')
    hist_1d_r_id = getObject('Results/th1f_2018ABCD.root', 'hist_1d_r_id')
    hist_1d_r_chi.SetTitle(';Vertex radius r [cm]; Dimuon vertex yield')
    hist_1d_r_chi.SetMaximum(100.0*hist_1d_r_chi.GetMaximum())
    hist_1d_r_chi.SetMinimum(0.1)
    hist_1d_r_chi.Sumw2()
    hist_1d_r_chi.SetLineColor(r.kBlack)
    hist_1d_r_char.SetLineColor(r.kBlack)
    hist_1d_r_cos.SetLineColor(r.kBlack)
    hist_1d_r_id.SetLineColor(r.kBlack)
    hist_1d_r_chi.SetFillColor(r.kMagenta-6)
    hist_1d_r_char.SetFillColor(r.kBlue-6)
    hist_1d_r_cos.SetFillColor(r.kCyan-6)
    hist_1d_r_id.SetFillColor(r.kGray)

    plot = Canvas.Canvas('hist_dr', 'png,pdf', 0.14, 0.75, 0.5, 0.89, 1)
    plot.addHisto(hist_1d_r_chi, 'HIST', 'Vertices passing kinematic cuts and #chi^{2}/ndf < 7.5', 'f', '', 1, 0)
    plot.addHisto(hist_1d_r_char, 'HIST,SAME', 'Vertices with OS dGlobals', 'f', '', 1, 1)
    plot.addHisto(hist_1d_r_cos, 'HIST,SAME', 'Vertices passing dR, mass, isolation and cos(#alpha) cuts', 'f', '', 1, 2)
    plot.addHisto(hist_1d_r_id, 'HIST,SAME', 'Vertices with dGlobals passing ID', 'f', '', 1, 3)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(1, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_2d_BPregion = getObject('Results/th1f_2018ABCD.root', 'hist_2d_BPregion')
    hist_2d_BPregion.SetTitle(';x [cm]; y [cm]')
    plot = Canvas.Canvas('hist_BPregion', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_BPregion, 'COLZ', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_2d_20region = getObject('Results/th1f_2018ABCD.root', 'hist_2d_20region')
    hist_2d_20region.SetTitle(';x [cm]; y [cm]')
    plot = Canvas.Canvas('hist_20region', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_20region, 'COLZ', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)



    ## SS pairs studies

    hist_2d_20region_SS = getObject('Results/th1f_2018ABCD.root', 'hist_2d_20region_SS')
    hist_2d_20region_SS.SetTitle(';x [cm]; y [cm]')
    plot = Canvas.Canvas('hist_20region_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_20region_SS, 'COLZ', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_r_SS = getObject('Results/th1f_2018ABCD.root', 'hist_1d_r_SS')
    hist_1d_r_SS.SetFillColor(r.kRed-9)
    hist_1d_r_SS.SetLineColor(r.kRed-6)
    hist_1d_r_SS.SetTitle('; Vertex radius r [cm]; SS Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_r_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_r_SS, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_dr_SS = getObject('Results/th1f_2018ABCD.root', 'hist_1d_dr_SS')
    hist_1d_dr_SS.SetFillColor(r.kRed-9)
    hist_1d_dr_SS.SetLineColor(r.kRed-6)
    hist_1d_dr_SS.SetTitle('; Dimuon #Delta R; SS Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_dr_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_dr_SS, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False, xlog = True)

    hist_1d_mass_SS = getObject('Results/th1f_2018ABCD.root', 'hist_1d_mass_SS')
    hist_1d_mass_SS.SetFillColor(r.kRed-9)
    hist_1d_mass_SS.SetLineColor(r.kRed-6)
    hist_1d_mass_SS.SetTitle('; Dimuon invariant mass [GeV]; SS Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_mass_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_mass_SS, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_iso_SS = getObject('Results/th1f_2018ABCD.root', 'hist_1d_iso_SS')
    hist_1d_iso_SS.SetFillColor(r.kRed-9)
    hist_1d_iso_SS.SetLineColor(r.kRed-6)
    hist_1d_iso_SS.SetTitle('; SS dGlobal isolation; SS dGlobal yield')
    plot = Canvas.Canvas('hist_1d_iso_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_iso_SS, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False, xlog = True)


    hist_1d_dPhi_SS = getObject('Results/th1f_2018ABCD.root', 'hist_1d_dPhi_SS')
    hist_1d_dPhi_SS.SetFillColor(r.kRed-9)
    hist_1d_dPhi_SS.SetLineColor(r.kRed-6)
    hist_1d_dPhi_SS.SetTitle('; Dimuon invariant dPhi [GeV]; SS Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_dPhi_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_dPhi_SS, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)


    hist_1d_dpt_SS = getObject('Results/th1f_2018ABCD.root', 'hist_1d_dpt_SS')
    hist_1d_dpt_SS.SetFillColor(r.kRed-9)
    hist_1d_dpt_SS.SetLineColor(r.kRed-6)
    hist_1d_dpt_SS.SetTitle('; Dimuon invariant dpt [GeV]; SS Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_dpt_SS', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_dpt_SS, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    
    ## SS pairs studies in the ring

    hist_2d_20region_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_2d_20region_SS_ring')
    hist_2d_20region_SS_ring.SetTitle(';x [cm]; y [cm]')
    plot = Canvas.Canvas('hist_20region_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_20region_SS_ring, 'COLZ', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_r_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_1d_r_SS_ring')
    hist_1d_r_SS_ring.SetFillColor(r.kRed-9)
    hist_1d_r_SS_ring.SetLineColor(r.kRed-6)
    hist_1d_r_SS_ring.SetTitle('; Vertex radius r [cm]; SS_ring Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_r_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_r_SS_ring, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_dr_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_1d_dr_SS_ring')
    hist_1d_dr_SS_ring.SetFillColor(r.kRed-9)
    hist_1d_dr_SS_ring.SetLineColor(r.kRed-6)
    hist_1d_dr_SS_ring.SetTitle('; Dimuon #Delta R; SS_ring Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_dr_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_dr_SS_ring, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False, xlog = True)

    hist_1d_mass_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_1d_mass_SS_ring')
    hist_1d_mass_SS_ring.SetFillColor(r.kRed-9)
    hist_1d_mass_SS_ring.SetLineColor(r.kRed-6)
    hist_1d_mass_SS_ring.SetTitle('; Dimuon invariant mass [GeV]; SS_ring Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_mass_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_mass_SS_ring, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_iso_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_1d_iso_SS_ring')
    hist_1d_iso_SS_ring.SetFillColor(r.kRed-9)
    hist_1d_iso_SS_ring.SetLineColor(r.kRed-6)
    hist_1d_iso_SS_ring.SetTitle('; SS_ring dGlobal isolation; SS_ring dGlobal yield')
    plot = Canvas.Canvas('hist_1d_iso_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_iso_SS_ring, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False, xlog = True)


    hist_1d_dPhi_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_1d_dPhi_SS_ring')
    hist_1d_dPhi_SS_ring.SetFillColor(r.kRed-9)
    hist_1d_dPhi_SS_ring.SetLineColor(r.kRed-6)
    hist_1d_dPhi_SS_ring.SetTitle('; Dimuon invariant dPhi [GeV]; SS_ring Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_dPhi_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_dPhi_SS_ring, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)


    hist_1d_dpt_SS_ring = getObject('Results/th1f_2018ABCD.root', 'hist_1d_dpt_SS_ring')
    hist_1d_dpt_SS_ring.SetFillColor(r.kRed-9)
    hist_1d_dpt_SS_ring.SetLineColor(r.kRed-6)
    hist_1d_dpt_SS_ring.SetTitle('; Dimuon invariant dpt [GeV]; SS_ring Dimuon vertex yield')
    plot = Canvas.Canvas('hist_1d_dpt_SS_ring', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_dpt_SS_ring, 'HIST', '', 'f', '', 1, 0)
    plot.addLatex(0.9, 0.93, '2018 Data (6.0 fb)', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    ### QCD
    hist_2d_BPregion_QCD = getObject('Results/th1f_MM_QCD.root', 'hist_2d_BPregion')
    hist_2d_BPregion_QCD.SetTitle(';x [cm]; y [cm]')
    hist_2d_BPregion_QCD.SetMarkerSize(0.4)
    plot = Canvas.Canvas('hist_BPregion_QCD_MM', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_BPregion_QCD, 'P', '', 'f', '', 1, 0, marker = 8)
    plot.addEllipse(0.171, -0.176, 2.210, 2.210, r.kRed)
    plot.addEllipse(0.0, 0.0, 1.5, 1.5, r.kGray, r.kGray)
    plot.addLatex(0.9, 0.93, 'QCD Control Region    2018 (59.8 fb^{-1})', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_2d_20region_QCD = getObject('Results/th1f_MM_QCD.root', 'hist_2d_20region')
    hist_2d_20region_QCD.SetTitle(';x [cm]; y [cm]')
    hist_2d_20region_QCD.SetMarkerSize(0.4)
    plot = Canvas.Canvas('hist_20region_QCD_MM', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_20region_QCD, 'P', '', 'f', '', 1, 0, marker = 8)
    plot.addEllipse(0.0, 0.0, 1.5, 1.5, r.kGray, r.kGray)
    plot.addLatex(0.9, 0.93, 'QCD Control Region    2018 (59.8 fb^{-1})', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_r_QCD = getObject('Results/th1f_MM_QCD.root', 'hist_1d_r_all')
    hist_1d_r_QCD.SetTitle(';Decay distance R [cm]; Counts')
    hist_1d_r_QCD.SetFillColor(r.kRed-9)
    hist_1d_r_QCD.SetLineColor(r.kRed-6)
    plot = Canvas.Canvas('hist_r_QCD_MM', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_r_QCD, 'HIST', '', 'f', '', 1, 0, marker = 8)
    plot.addLatex(0.9, 0.93, 'QCD Control Region    2018 (59.8 fb^{-1})', size = 0.026, align = 31)
    plot.save(0, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_2d_BPregion_QCD = getObject('Results/th1f_EE_QCD.root', 'hist_2d_BPregion')
    hist_2d_BPregion_QCD.SetTitle(';x [cm]; y [cm]')
    hist_2d_BPregion_QCD.SetMarkerSize(0.4)
    plot = Canvas.Canvas('hist_BPregion_QCD_EE', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_BPregion_QCD, 'COLZ', '', 'f', '', 1, 0, marker = 8)
    #plot.addEllipse(0.171, -0.176, 2.210, 2.210, r.kRed)
    plot.addEllipse(0.0, 0.0, 1.5, 1.5, r.kGray, r.kGray)
    plot.addLatex(0.9, 0.93, 'QCD Control Region    2018 (59.8 fb^{-1})', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_2d_20region_QCD = getObject('Results/th1f_EE_QCD.root', 'hist_2d_20region')
    hist_2d_20region_QCD.SetTitle(';x [cm]; y [cm]')
    hist_2d_20region_QCD.SetMarkerSize(0.4)
    plot = Canvas.Canvas('hist_20region_QCD_EE', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_2d_20region_QCD, 'COLZ', '', 'f', '', 1, 0, marker = 8)
    plot.addEllipse(0.0, 0.0, 1.5, 1.5, r.kGray, r.kGray)
    plot.addLatex(0.9, 0.93, 'QCD Control Region    2018 (59.8 fb^{-1})', size = 0.026, align = 31)
    plot.save(0, 1, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

    hist_1d_r_QCD = getObject('Results/th1f_EE_QCD.root', 'hist_1d_r_all')
    hist_1d_r_QCD.SetTitle(';Decay distance R [cm]; Counts')
    hist_1d_r_QCD.SetFillColor(r.kRed-9)
    hist_1d_r_QCD.SetLineColor(r.kRed-6)
    plot = Canvas.Canvas('hist_r_QCD_EE', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_1d_r_QCD, 'HIST', '', 'f', '', 1, 0, marker = 8)
    plot.addLatex(0.9, 0.93, 'QCD Control Region    2018 (59.8 fb^{-1})', size = 0.026, align = 31)
    plot.save(0, 1, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)

