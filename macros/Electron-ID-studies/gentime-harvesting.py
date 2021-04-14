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

def makeOFHisto(h):

    """
    Function to make overflow histograms
       Parameter: Histogram
       Return:    Same histogram with an additional bin with events out of x axis
    """

    nbin = h.GetNbinsX()
    bw = h.GetBinWidth(1)
    xmin = h.GetXaxis().GetBinLowEdge(1)
    xmax = h.GetXaxis().GetBinUpEdge(nbin)
    _h = r.TH1F(h.GetName() + '_OF', '', nbin + 1, xmin, xmax + bw)

    for _bin in range(1, nbin + 2):
        _h.SetBinContent(_bin, h.GetBinContent(_bin))
        _h.SetBinError(_bin, h.GetBinError(_bin))

    _h.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
    _h.GetYaxis().SetTitle(h.GetYaxis().GetTitle())

    return _h



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

    ########################
    ####   Histograms   ####
    ########################

    h400_50_4 = getObject('GenResults/th1fHXX_400_50_4mm.root', 'hist_Lxy')
    h400_50_40 = getObject('GenResults/th1fHXX_400_50_40mm.root', 'hist_Lxy')
    h400_50_400 = getObject('GenResults/th1fHXX_400_50_400mm.root', 'hist_Lxy')
    h400_150_400 = getObject('GenResults/th1fHXX_400_150_400mm.root', 'hist_Lxy')
    h1000_150_10 = getObject('GenResults/th1fHXX_1000_150_10mm.root', 'hist_Lxy')
    h1000_150_100 = getObject('GenResults/th1fHXX_1000_150_100mm.root', 'hist_Lxy')
    h1000_350_35 = getObject('GenResults/th1fHXX_1000_350_35mm.root', 'hist_Lxy')
    h1000_350_350 = getObject('GenResults/th1fHXX_1000_350_350mm.root', 'hist_Lxy')

    h400_50_4.Scale(1.0/h400_50_4.GetEntries())
    h400_50_40.Scale(1.0/h400_50_40.GetEntries())
    h400_50_400.Scale(1.0/h400_50_400.GetEntries())
    h400_150_400.Scale(1.0/h400_150_400.GetEntries())
    h1000_150_10.Scale(1.0/h1000_150_10.GetEntries())
    h1000_150_100.Scale(1.0/h1000_150_100.GetEntries())
    h1000_350_35.Scale(1.0/h1000_350_35.GetEntries())
    h1000_350_350.Scale(1.0/h1000_350_350.GetEntries())

    h400_50_4.SetMaximum(100.0)

    HIST_ptRes = Canvas.Canvas("HIST_genLxy", 'png', 0.17, 0.66, 0.55, 0.9, 1) 
    HIST_ptRes.addHisto(h400_50_4, 'P', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kGreen+2, True, 0, marker = 20)
    HIST_ptRes.addHisto(h400_50_40, 'P,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kSpring, True, 0, marker = 20)
    HIST_ptRes.addHisto(h400_50_400, 'P,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kAzure, True, 0, marker = 20)
    HIST_ptRes.addHisto(h400_150_400, 'P,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kBlue+2, True, 0, marker = 20)
    HIST_ptRes.addHisto(h1000_150_10, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 150 GeV, c#tau = 10 mm', 'p', r.kMagenta+2, True, 0, marker = 20)
    HIST_ptRes.addHisto(h1000_150_100, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 150 GeV, c#tau = 100 mm', 'p', r.kRed, True, 0, marker = 20)
    HIST_ptRes.addHisto(h1000_350_35, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 350 GeV, c#tau = 35 mm', 'p', r.kOrange, True, 0, marker = 20)
    HIST_ptRes.addHisto(h1000_350_350, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 350 GeV, c#tau = 350 mm', 'p', r.kYellow, True, 0, marker = 20)
    HIST_ptRes.save(1, 0, 1, '', '', outputDir = WORKPATH + 'harvested_time_'+opts.tag+'/', inProgress = False)


    h400_50_4 = getObject('GenResults/th1fHXX_400_50_4mm.root', 'hist_t')
    h400_50_40 = getObject('GenResults/th1fHXX_400_50_40mm.root', 'hist_t')
    h400_50_400 = getObject('GenResults/th1fHXX_400_50_400mm.root', 'hist_t')
    h400_150_400 = getObject('GenResults/th1fHXX_400_150_400mm.root', 'hist_t')
    h1000_150_10 = getObject('GenResults/th1fHXX_1000_150_10mm.root', 'hist_t')
    h1000_150_100 = getObject('GenResults/th1fHXX_1000_150_100mm.root', 'hist_t')
    h1000_350_35 = getObject('GenResults/th1fHXX_1000_350_35mm.root', 'hist_t')
    h1000_350_350 = getObject('GenResults/th1fHXX_1000_350_350mm.root', 'hist_t')

    h400_50_4.Scale(1.0/h400_50_4.GetEntries())
    h400_50_40.Scale(1.0/h400_50_40.GetEntries())
    h400_50_400.Scale(1.0/h400_50_400.GetEntries())
    h400_150_400.Scale(1.0/h400_150_400.GetEntries())
    h1000_150_10.Scale(1.0/h1000_150_10.GetEntries())
    h1000_150_100.Scale(1.0/h1000_150_100.GetEntries())
    h1000_350_35.Scale(1.0/h1000_350_35.GetEntries())
    h1000_350_350.Scale(1.0/h1000_350_350.GetEntries())

    print(h400_50_4.Integral())
    #h400_50_4.SetMaximum(10.0)
    h400_50_4.SetMaximum(1.5*h400_50_4.GetMaximum())
    h400_50_4.SetMaximum(100.0)


    P_ptRes = Canvas.Canvas("HIST_gent_in", 'png', 0.17, 0.66, 0.55, 0.9, 1) 
    P_ptRes.addHisto(h400_50_4, 'P', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kGreen+2, True, 0, marker = 20)
    P_ptRes.addHisto(h400_50_40, 'P,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kSpring, True, 0, marker = 20)
    P_ptRes.addHisto(h400_50_400, 'P,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kAzure, True, 0, marker = 20)
    P_ptRes.addHisto(h400_150_400, 'P,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kBlue+2, True, 0, marker = 20)
    P_ptRes.addHisto(h1000_150_10, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 150 GeV, c#tau = 10 mm', 'p', r.kMagenta+2, True, 0, marker = 20)
    P_ptRes.addHisto(h1000_150_100, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 150 GeV, c#tau = 100 mm', 'p', r.kRed, True, 0, marker = 20)
    P_ptRes.addHisto(h1000_350_35, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 350 GeV, c#tau = 35 mm', 'p', r.kOrange, True, 0, marker = 20)
    P_ptRes.addHisto(h1000_350_350, 'P,SAME', 'H#rightarrowXX: m_{H} = 1000 GeV, m_{X} = 350 GeV, c#tau = 350 mm', 'p', r.kYellow, True, 0, marker = 20)
    P_ptRes.save(1, 0, 1, '', '', outputDir = WORKPATH + 'harvested_time_'+opts.tag+'/', inProgress = False)
