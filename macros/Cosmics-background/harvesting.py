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
    

    hist_dPhi = getObject('Results/th1f_all.root', 'hist_dPhi')
    hist_dPhi.Rebin(4)
    hist_dPhi.SetTitle(';Dimuon collinearity |#Delta#Phi|; Dimuon vertex yield')
    hist_dPhi.SetMaximum(1.4*hist_dPhi.GetMaximum())
    hist_dPhi.SetMinimum(0.0)
    hist_dPhi.Sumw2()
    plot = Canvas.Canvas('hist_cosmics_dPhi', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_dPhi, 'P', '', 'p', r.kBlue, 1, 0)
    plot.addLatex(0.9, 0.93, 'Cosmic data 2018', size = 0.026, align = 31)
    plot.addLatex(0.85, 0.85, '/NoBPTX/Run2018B-12Nov2019_UL2018-v1/AOD', size = 0.026, align = 31)
    plot.addLatex(0.85, 0.8, 'HLT_L2Mu10_NoVertex_NoBPTX3BX', size = 0.026, align = 31)
    plot.save(0, 0, 0, '', '', outputDir = WORKPATH + 'harvested/', is2d = True, isPrivate = True, maxYnumbers = 3)


    hist_alpha = getObject('Results/th1f_all.root', 'hist_Alpha')
    hist_alpha.SetTitle(';Dimuon 3D angle #alpha; Dimuon vertex yield')
    hist_alpha.SetMaximum(10.0*hist_alpha.GetMaximum())
    hist_alpha.SetMinimum(0.1)
    hist_alpha.Sumw2()
    plot = Canvas.Canvas('hist_cosmics_alpha', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_alpha, 'P', '', 'p', r.kBlue, 1, 0)
    plot.addLatex(0.9, 0.93, 'Cosmic data 2018', size = 0.026, align = 31)
    plot.addLatex(0.85, 0.85, '/NoBPTX/Run2018B-12Nov2019_UL2018-v1/AOD', size = 0.026, align = 31)
    plot.addLatex(0.85, 0.8, 'HLT_L2Mu10_NoVertex_NoBPTX3BX', size = 0.026, align = 31)
    plot.save(0, 0, 1, '', '', outputDir = WORKPATH + 'harvested/', is2d = True, isPrivate = True)


    hist_cosalpha = getObject('Results/th1f_all.root', 'hist_cosAlpha')
    print(hist_cosalpha.GetEntries())
    hist_cosalpha.SetTitle(';Dimuon cos(#alpha); Dimuon vertex yield')
    hist_cosalpha.SetMaximum(1e7)
    hist_cosalpha.SetMinimum(0.1)
    hist_cosalpha.Sumw2()
    plot = Canvas.Canvas('hist_cosmics_cosalpha', 'png,pdf', 0.3, 0.75, 0.89, 0.8, 1)
    plot.addHisto(hist_cosalpha, 'P', '', 'p', r.kBlue, 1, 0)
    plot.addLine(-0.8, 0.1, -0.8, 1e5,r.kRed)
    plot.addLine(-0.9, 0.1, -0.9, 1e5,r.kRed)
    plot.addLatex(0.9, 0.93, 'Cosmic data 2018', size = 0.026, align = 31)
    plot.addLatex(0.85, 0.85, '/NoBPTX/Run2018B-12Nov2019_UL2018-v1/AOD', size = 0.026, align = 31)
    plot.addLatex(0.85, 0.8, 'HLT_L2Mu10_NoVertex_NoBPTX3BX', size = 0.026, align = 31)
    plot.save(0, 0, 1, '', '', outputDir = WORKPATH + 'harvested/', is2d = True, isPrivate = True)


