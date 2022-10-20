import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
import math, sys, optparse, array, copy, os
import gc, inspect
import numpy as np

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
import include.Sample as Sample
import include.helper as helper
import include.CutManager as CutManager

from include.galapagoStyle import sigpalette, gcolors, dcolors


##################################### FUNCTION DEFINITION ########################################

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

def doFakeStack(matched, unmatched, doOF = True, fakedown = False, ymin = False, ymax = False):

    matched = makeOFHisto(matched)
    unmatched = makeOFHisto(unmatched)

    matched.SetLineColor(r.kBlack)
    unmatched.SetLineColor(r.kBlack)
    unmatched.SetFillColorAlpha(r.kRed+1, 0.7)
    matched.SetFillColorAlpha(r.kGreen+2, 0.7)
    matched.SetTitle('True displaced global')
    unmatched.SetTitle('Fake displaced global')


    Stack = r.THStack("aux_stack", ";"+matched.GetXaxis().GetTitle()+";DG yield")

    if fakedown:
        Stack.Add(matched)    
        Stack.Add(unmatched)    
    else:
        Stack.Add(unmatched)    
        Stack.Add(matched)    

    return Stack


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
    r.setTDRStyle()

    ###########################
    ####   Parser object   ####
    ###########################
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    (opts, args) = parser.parse_args()

    #### No cuts applied on gen muons


    filenames = {}
    filenames['1000_150_1'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-1_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220521_082831/0000/output_1.root'
    filenames['1000_150_10'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-10_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220521_082826/0000/output_1.root'
    filenames['1000_150_100'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-100_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220521_082820/0000/output_1.root'
    filenames['1000_150_1000'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-1000_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220521_081715/0000/output_1.root'
    filenames['1000_150_10000'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-10000_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220521_081710/0000/output_1.root'

    ### Efficiency vs pt
    
    dglEff_pt = getObject(filenames['1000_150_1'], 'dglEff_pt')
    glEff_pt = getObject(filenames['1000_150_1'], 'glEff_pt')
    tkEff_pt = getObject(filenames['1000_150_1'], 'tkEff_pt')

    plot_ = Canvas.Canvas('recoEff_pt', 'png', 0.16, 0.76, 0.35, 0.89, 1) 
    plot_.addRate(dglEff_pt, 'AP', 'Displaced Global efficiency', 'p', r.kBlack, 1, 0, marker = 20)
    plot_.addRate(glEff_pt, 'P,SAME', 'Global efficiency', 'p', r.kBlue, 1, 1, marker = 24)
    plot_.addRate(tkEff_pt, 'P,SAME', 'Tracker efficiency', 'p', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, 'H#rightarrowSS (2018 UL)', size = 0.035, align = 31, font = 42)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/')

    dglEff_Lxy = getObject(filenames['1000_150_1'], 'dglEff_Lxy')
    dglEff_Lxy.Add(getObject(filenames['1000_150_10'], 'dglEff_Lxy'))
    dglEff_Lxy.Add(getObject(filenames['1000_150_100'], 'dglEff_Lxy'))
    dglEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'dglEff_Lxy'))
    dglEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'dglEff_Lxy'))
    glEff_Lxy = getObject(filenames['1000_150_1'], 'glEff_Lxy')
    glEff_Lxy.Add(getObject(filenames['1000_150_10'], 'glEff_Lxy'))
    glEff_Lxy.Add(getObject(filenames['1000_150_100'], 'glEff_Lxy'))
    glEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'glEff_Lxy'))
    glEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'glEff_Lxy'))
    tkEff_Lxy = getObject(filenames['1000_150_1'], 'tkEff_Lxy')
    tkEff_Lxy.Add(getObject(filenames['1000_150_10'], 'tkEff_Lxy'))
    tkEff_Lxy.Add(getObject(filenames['1000_150_100'], 'tkEff_Lxy'))
    tkEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'tkEff_Lxy'))
    tkEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'tkEff_Lxy'))

    plot_ = Canvas.Canvas('recoEff_Lxy', 'png', 0.16, 0.76, 0.35, 0.89, 1) 
    plot_.addRate(dglEff_Lxy, 'AP', 'Displaced Global efficiency', 'p', r.kBlack, 1, 0, marker = 20)
    plot_.addRate(glEff_Lxy, 'P,SAME', 'Global efficiency', 'p', r.kBlue, 1, 1, marker = 24)
    plot_.addRate(tkEff_Lxy, 'P,SAME', 'Tracker efficiency', 'p', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, 'H#rightarrowSS (2018 UL)', size = 0.035, align = 31, font = 42)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/')



    ### Resolution plot (DGL)

    dgl_ptRes_1 = getObject(filenames['1000_150_1'], 'dgl_ptRes')
    dgl_ptRes_10 = getObject(filenames['1000_150_10'], 'dgl_ptRes')
    dgl_ptRes_100 = getObject(filenames['1000_150_100'], 'dgl_ptRes')
    dgl_ptRes_1000 = getObject(filenames['1000_150_1000'], 'dgl_ptRes')
    dgl_ptRes_1.Scale(1./dgl_ptRes_1.Integral())
    dgl_ptRes_10.Scale(1./dgl_ptRes_10.Integral())
    dgl_ptRes_100.Scale(1./dgl_ptRes_100.Integral())
    dgl_ptRes_1000.Scale(1./dgl_ptRes_1000.Integral())

    plot_ = Canvas.Canvas('dgl_ptRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(dgl_ptRes_1, 'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], 1, 0)
    plot_.addHisto(dgl_ptRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], 1, 1)
    plot_.addHisto(dgl_ptRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], 1, 2)
    plot_.addHisto(dgl_ptRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/')

    ### Resolution plot (TK)

    tk_ptRes_1 = getObject(filenames['1000_150_1'], 'tk_ptRes')
    tk_ptRes_10 = getObject(filenames['1000_150_10'], 'tk_ptRes')
    tk_ptRes_100 = getObject(filenames['1000_150_100'], 'tk_ptRes')
    tk_ptRes_1000 = getObject(filenames['1000_150_1000'], 'tk_ptRes')
    tk_ptRes_10000 = getObject(filenames['1000_150_10000'], 'tk_ptRes')

    plot_ = Canvas.Canvas('tk_ptRes', 'png', 0.16, 0.76, 0.35, 0.89, 1) 
    plot_.addHisto(tk_ptRes_1, 'HIST', '1', 'p', r.kBlack, 1, 0)
    plot_.addHisto(tk_ptRes_10, 'HIST,SAME', '1', 'p', r.kRed, 1, 1)
    plot_.addHisto(tk_ptRes_100, 'HIST,SAME', '1', 'p', r.kBlue, 1, 2)
    plot_.addHisto(tk_ptRes_1000, 'HIST,SAME', '1', 'p', r.kGreen, 1, 3)
    plot_.addLatex(0.9, 0.93, 'H#rightarrowSS (2018 UL)', size = 0.035, align = 31, font = 42)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/')

