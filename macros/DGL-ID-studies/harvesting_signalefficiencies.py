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

from include.galapagoStyle import sigpalette, gcolors, dcolors, acolors



################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]:
    WORKPATH += level
    WORKPATH += '/'

print('runningfile: ' + runningfile)

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

def getObject(filename, key):

    _f = r.TFile(filename)
    _h = _f.Get(key)
    _hcopy = copy.deepcopy(_h)
    r.SetOwnership(_hcopy, 0)
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


def copyHisto(h):

    """
    Function to make overflow histograms
       Parameter: Histogram
       Return:    Same histogram with an additional bin with events out of x axis
    """

    nbin = h.GetNbinsX()
    bw = h.GetBinWidth(1)
    xmin = h.GetXaxis().GetBinLowEdge(1)
    xmax = h.GetXaxis().GetBinUpEdge(nbin)
    _h = r.TH1F(h.GetName() + '_copy', '', nbin, xmin, xmax)

    for _bin in range(1, nbin + 1):
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

    #### No cuts applied on gen muons


    filenames = {}
    filenames['1000_150_1'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-1_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220526_110419/0000/output_1.root'
    filenames['1000_150_10'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-10_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220526_110411/0000/output_1.root'
    filenames['1000_150_100'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-100_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220526_105806/0000/output_1.root'
    filenames['1000_150_1000'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-1000_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220526_105758/0000/output_1.root'
    filenames['1000_150_10000'] = '/eos/user/f/fernance/MuonEff_HSS/ggH_HToSSTo4l_MH-1000_MS-150_ctauS-10000_TuneCP5_13TeV-powheg-pythia8/LLP_RECOIDresults/220526_105650/0000/output_1.root'


    ### Efficiency vs pt
    
    dglEff_pt = getObject(filenames['1000_150_1'], 'dglEff_pt')
    glEff_pt = getObject(filenames['1000_150_1'], 'glEff_pt')
    tkEff_pt = getObject(filenames['1000_150_1'], 'tkEff_pt')

    plot_ = Canvas.Canvas('recoEff_pt_1mm', 'png,pdf', 0.35, 0.26, 0.7, 0.39, 1) 
    plot_.addHisto(dglEff_pt, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_pt, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_pt, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.34, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.34, 0.81, 'M_{H} = 1000 GeV, M_{S} = 150 GeV, c#tau = 1 mm', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

    dglEff_pt = getObject(filenames['1000_150_10'], 'dglEff_pt')
    glEff_pt = getObject(filenames['1000_150_10'], 'glEff_pt')
    tkEff_pt = getObject(filenames['1000_150_10'], 'tkEff_pt')

    plot_ = Canvas.Canvas('recoEff_pt_10mm', 'png,pdf', 0.35, 0.26, 0.7, 0.39, 1)
    plot_.addHisto(dglEff_pt, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_pt, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_pt, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.34, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.34, 0.81, 'M_{H} = 1000 GeV, M_{S} = 150 GeV, c#tau = 10 mm', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

    dglEff_pt = getObject(filenames['1000_150_100'], 'dglEff_pt')
    glEff_pt = getObject(filenames['1000_150_100'], 'glEff_pt')
    tkEff_pt = getObject(filenames['1000_150_100'], 'tkEff_pt')

    plot_ = Canvas.Canvas('recoEff_pt_100mm', 'png,pdf', 0.35, 0.26, 0.7, 0.39, 1)
    plot_.addHisto(dglEff_pt, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_pt, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_pt, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.34, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.34, 0.81, 'M_{H} = 1000 GeV, M_{S} = 150 GeV, c#tau = 100 mm', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)


    ### Efficiency vs eta
    
    dglEff_eta = getObject(filenames['1000_150_1'], 'dglEff_eta')
    glEff_eta = getObject(filenames['1000_150_1'], 'glEff_eta')
    tkEff_eta = getObject(filenames['1000_150_1'], 'tkEff_eta')

    plot_ = Canvas.Canvas('recoEff_eta_1mm', 'png,pdf', 0.35, 0.26, 0.7, 0.39, 1) 
    plot_.addHisto(dglEff_eta, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_eta, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_eta, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.34, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.34, 0.81, 'M_{H} = 1000 GeV, M_{S} = 150 GeV, c#tau = 1 mm', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

    dglEff_eta = getObject(filenames['1000_150_10'], 'dglEff_eta')
    glEff_eta = getObject(filenames['1000_150_10'], 'glEff_eta')
    tkEff_eta = getObject(filenames['1000_150_10'], 'tkEff_eta')

    plot_ = Canvas.Canvas('recoEff_eta_10mm', 'png,pdf', 0.35, 0.26, 0.7, 0.39, 1)
    plot_.addHisto(dglEff_eta, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_eta, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_eta, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.34, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.34, 0.81, 'M_{H} = 1000 GeV, M_{S} = 150 GeV, c#tau = 10 mm', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

    dglEff_eta = getObject(filenames['1000_150_100'], 'dglEff_eta')
    glEff_eta = getObject(filenames['1000_150_100'], 'glEff_eta')
    tkEff_eta = getObject(filenames['1000_150_100'], 'tkEff_eta')

    plot_ = Canvas.Canvas('recoEff_eta_100mm', 'png,pdf', 0.35, 0.26, 0.7, 0.39, 1)
    plot_.addHisto(dglEff_eta, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_eta, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_eta, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.34, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.34, 0.81, 'M_{H} = 1000 GeV, M_{S} = 150 GeV, c#tau = 100 mm', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)


    ### Resolution plot (DGL)

    dgl_ptRes_1 = copyHisto(getObject(filenames['1000_150_1'], 'dgl_ptRes'))
    dgl_ptRes_10 = copyHisto(getObject(filenames['1000_150_10'], 'dgl_ptRes'))
    dgl_ptRes_100 = copyHisto(getObject(filenames['1000_150_100'], 'dgl_ptRes'))
    dgl_ptRes_1000 = copyHisto(getObject(filenames['1000_150_1000'], 'dgl_ptRes'))
    dgl_ptRes_1.Scale(1./dgl_ptRes_1.Integral())
    dgl_ptRes_10.Scale(1./dgl_ptRes_10.Integral())
    dgl_ptRes_100.Scale(1./dgl_ptRes_100.Integral())
    dgl_ptRes_1000.Scale(1./dgl_ptRes_1000.Integral())
    dgl_ptRes_1.SetLineWidth(2)
    dgl_ptRes_10.SetLineWidth(2)
    dgl_ptRes_100.SetLineWidth(2)
    dgl_ptRes_1000.SetLineWidth(2)
    plot_ = Canvas.Canvas('dgl_ptRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(dgl_ptRes_1, 'HIST', 'c#tau = 1 mm', 'l', acolors['1'], 1, 0)
    plot_.addHisto(dgl_ptRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', acolors['2'], 1, 1)
    plot_.addHisto(dgl_ptRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', acolors['3'], 1, 2)
    plot_.addHisto(dgl_ptRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', acolors['4'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.14, 0.93, 'Displaced Global Muons', size = 0.035, align = 11)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.0, ymax = 0.15, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True)

    ### Resolution plot (TK)

    tk_ptRes_1 = copyHisto(getObject(filenames['1000_150_1'], 'tk_ptRes'))
    tk_ptRes_10 = copyHisto(getObject(filenames['1000_150_10'], 'tk_ptRes'))
    tk_ptRes_100 = copyHisto(getObject(filenames['1000_150_100'], 'tk_ptRes'))
    tk_ptRes_1000 = copyHisto(getObject(filenames['1000_150_1000'], 'tk_ptRes'))
    tk_ptRes_1.Scale(1./tk_ptRes_1.Integral())
    tk_ptRes_10.Scale(1./tk_ptRes_10.Integral())
    tk_ptRes_100.Scale(1./tk_ptRes_100.Integral())
    tk_ptRes_1000.Scale(1./tk_ptRes_1000.Integral())
    tk_ptRes_1.SetMaximum(0.15)
    tk_ptRes_1.SetLineWidth(2)
    tk_ptRes_10.SetLineWidth(2)
    tk_ptRes_100.SetLineWidth(2)
    tk_ptRes_1000.SetLineWidth(2)
    plot_ = Canvas.Canvas('tk_ptRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(tk_ptRes_1, 'HIST', 'c#tau = 1 mm', 'l', acolors['1'], 1, 0)
    plot_.addHisto(tk_ptRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', acolors['2'], 1, 1)
    plot_.addHisto(tk_ptRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', acolors['3'], 1, 2)
    plot_.addHisto(tk_ptRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', acolors['4'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.14, 0.93, 'Tracker Muons', size = 0.035, align = 11)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.0, ymax = 0.15, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True)

    ### Resolution plot (GL)

    gl_ptRes_1 = copyHisto(getObject(filenames['1000_150_1'], 'gl_ptRes'))
    gl_ptRes_10 = copyHisto(getObject(filenames['1000_150_10'], 'gl_ptRes'))
    gl_ptRes_100 = copyHisto(getObject(filenames['1000_150_100'], 'gl_ptRes'))
    gl_ptRes_1000 = copyHisto(getObject(filenames['1000_150_1000'], 'gl_ptRes'))
    gl_ptRes_1.Scale(1./gl_ptRes_1.Integral())
    gl_ptRes_10.Scale(1./gl_ptRes_10.Integral())
    gl_ptRes_100.Scale(1./gl_ptRes_100.Integral())
    gl_ptRes_1000.Scale(1./gl_ptRes_1000.Integral())
    gl_ptRes_1.SetLineWidth(2)
    gl_ptRes_10.SetLineWidth(2)
    gl_ptRes_100.SetLineWidth(2)
    gl_ptRes_1000.SetLineWidth(2)
    plot_ = Canvas.Canvas('gl_ptRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(gl_ptRes_1, 'HIST', 'c#tau = 1 mm', 'l', acolors['1'], 1, 0)
    plot_.addHisto(gl_ptRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', acolors['2'], 1, 1)
    plot_.addHisto(gl_ptRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', acolors['3'], 1, 2)
    plot_.addHisto(gl_ptRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', acolors['4'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.14, 0.93, 'Global Muons', size = 0.035, align = 11)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.0, ymax = 0.15, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True)

    ### Resolution plot (DGL)

    dgl_dxyRes_1 = getObject(filenames['1000_150_1'], 'dgl_dxyRes').Clone()
    dgl_dxyRes_10 = getObject(filenames['1000_150_10'], 'dgl_dxyRes').Clone()
    dgl_dxyRes_100 = getObject(filenames['1000_150_100'], 'dgl_dxyRes').Clone()
    dgl_dxyRes_1000 = getObject(filenames['1000_150_1000'], 'dgl_dxyRes').Clone()
    dgl_dxyRes_1.Scale(1./dgl_dxyRes_1.Integral())
    dgl_dxyRes_10.Scale(1./dgl_dxyRes_10.Integral())
    dgl_dxyRes_100.Scale(1./dgl_dxyRes_100.Integral())
    dgl_dxyRes_1000.Scale(1./dgl_dxyRes_1000.Integral())

    plot_ = Canvas.Canvas('dgl_dxyRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(dgl_dxyRes_1, 'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], 1, 0)
    plot_.addHisto(dgl_dxyRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], 1, 1)
    plot_.addHisto(dgl_dxyRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], 1, 2)
    plot_.addHisto(dgl_dxyRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.14, 0.93, 'Displaced Global Muons', size = 0.035, align = 11)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.0, ymax = 0.7, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True)

    ### Resolution plot (TK)

    tk_dxyRes_1 = getObject(filenames['1000_150_1'], 'tk_dxyRes')
    tk_dxyRes_10 = getObject(filenames['1000_150_10'], 'tk_dxyRes')
    tk_dxyRes_100 = getObject(filenames['1000_150_100'], 'tk_dxyRes')
    tk_dxyRes_1000 = getObject(filenames['1000_150_1000'], 'tk_dxyRes')
    tk_dxyRes_1.Scale(1./tk_dxyRes_1.Integral())
    tk_dxyRes_10.Scale(1./tk_dxyRes_10.Integral())
    tk_dxyRes_100.Scale(1./tk_dxyRes_100.Integral())
    tk_dxyRes_1000.Scale(1./tk_dxyRes_1000.Integral())

    plot_ = Canvas.Canvas('tk_dxyRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(tk_dxyRes_1, 'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], 1, 0)
    plot_.addHisto(tk_dxyRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], 1, 1)
    plot_.addHisto(tk_dxyRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], 1, 2)
    plot_.addHisto(tk_dxyRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.14, 0.93, 'Tracker Muons', size = 0.035, align = 11)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.0, ymax = 0.7, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True)

    ### Resolution plot (GL)

    gl_dxyRes_1 = getObject(filenames['1000_150_1'], 'gl_dxyRes')
    gl_dxyRes_10 = getObject(filenames['1000_150_10'], 'gl_dxyRes')
    gl_dxyRes_100 = getObject(filenames['1000_150_100'], 'gl_dxyRes')
    gl_dxyRes_1000 = getObject(filenames['1000_150_1000'], 'gl_dxyRes')
    gl_dxyRes_1.Scale(1./gl_dxyRes_1.Integral())
    gl_dxyRes_10.Scale(1./gl_dxyRes_10.Integral())
    gl_dxyRes_100.Scale(1./gl_dxyRes_100.Integral())
    gl_dxyRes_1000.Scale(1./gl_dxyRes_1000.Integral())

    plot_ = Canvas.Canvas('gl_dxyRes', 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1) 
    plot_.addHisto(gl_dxyRes_1, 'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], 1, 0)
    plot_.addHisto(gl_dxyRes_10, 'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], 1, 1)
    plot_.addHisto(gl_dxyRes_100, 'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], 1, 2)
    plot_.addHisto(gl_dxyRes_1000, 'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], 1, 3)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.14, 0.93, 'Global Muons', size = 0.035, align = 11)
    plot_.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    plot_.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.0, ymax = 0.7, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True)

    ### Effiiciency vs Lxy

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

    plot_ = Canvas.Canvas('recoEff_Lxy', 'png,pdf', 0.4, 0.76, 0.75, 0.89, 1) 
    plot_.addHisto(dglEff_Lxy, 'AP', 'Displaced Global Muon efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(glEff_Lxy, 'P,SAME', 'Global Muon efficiency', 'pl', r.kBlue, 1, 1, marker = 24)
    plot_.addHisto(tkEff_Lxy, 'P,SAME', 'Tracker Muon efficiency', 'pl', r.kRed, 1, 2, marker = 25)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.2, 0.27, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.2, 0.23, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    plot_.addLatex(0.2, 0.19, '(all lifetimes combined)', size = 0.03)
    plot_.save(1, 0, 0, '','', outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

    ### ID efficiency (vs Lxy)
    idEff_Lxy = getObject(filenames['1000_150_1'], 'idEff_Lxy')
    idEff_Lxy.Add(getObject(filenames['1000_150_10'], 'idEff_Lxy'))
    idEff_Lxy.Add(getObject(filenames['1000_150_100'], 'idEff_Lxy'))
    idEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'idEff_Lxy'))
    idEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'idEff_Lxy'))
    muHitEff_Lxy = getObject(filenames['1000_150_1'], 'muHitEff_Lxy')
    muHitEff_Lxy.Add(getObject(filenames['1000_150_10'], 'muHitEff_Lxy'))
    muHitEff_Lxy.Add(getObject(filenames['1000_150_100'], 'muHitEff_Lxy'))
    muHitEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'muHitEff_Lxy'))
    muHitEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'muHitEff_Lxy'))
    tkHitEff_Lxy = getObject(filenames['1000_150_1'], 'tkHitEff_Lxy')
    tkHitEff_Lxy.Add(getObject(filenames['1000_150_10'], 'tkHitEff_Lxy'))
    tkHitEff_Lxy.Add(getObject(filenames['1000_150_100'], 'tkHitEff_Lxy'))
    tkHitEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'tkHitEff_Lxy'))
    tkHitEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'tkHitEff_Lxy'))
    ptResEff_Lxy = getObject(filenames['1000_150_1'], 'ptResEff_Lxy')
    ptResEff_Lxy.Add(getObject(filenames['1000_150_10'], 'ptResEff_Lxy'))
    ptResEff_Lxy.Add(getObject(filenames['1000_150_100'], 'ptResEff_Lxy'))
    ptResEff_Lxy.Add(getObject(filenames['1000_150_1000'], 'ptResEff_Lxy'))
    ptResEff_Lxy.Add(getObject(filenames['1000_150_10000'], 'ptResEff_Lxy'))

    plot_ = Canvas.Canvas('idEff_Lxy', 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1) 
    plot_.addHisto(idEff_Lxy, 'AP', 'Full ID efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(muHitEff_Lxy, 'P,SAME', 'Muon hit selection', 'pl', gcolors['red'], 1, 0, marker = 24)
    plot_.addHisto(tkHitEff_Lxy, 'P,SAME', 'Track hit selection', 'pl', gcolors['blue'], 1, 0, marker = 24)
    plot_.addHisto(ptResEff_Lxy, 'P,SAME', 'p_{T} resolution selection', 'pl', gcolors['green'], 1, 0, marker = 24)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    plot_.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.7, ymax = 1.1, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)


    ### ID efficiency (vs pT)
    idEff_pt = getObject(filenames['1000_150_1'], 'idEff_pt')
    idEff_pt.Add(getObject(filenames['1000_150_10'], 'idEff_pt'))
    idEff_pt.Add(getObject(filenames['1000_150_100'], 'idEff_pt'))
    idEff_pt.Add(getObject(filenames['1000_150_1000'], 'idEff_pt'))
    idEff_pt.Add(getObject(filenames['1000_150_10000'], 'idEff_pt'))
    muHitEff_pt = getObject(filenames['1000_150_1'], 'muHitEff_pt')
    muHitEff_pt.Add(getObject(filenames['1000_150_10'], 'muHitEff_pt'))
    muHitEff_pt.Add(getObject(filenames['1000_150_100'], 'muHitEff_pt'))
    muHitEff_pt.Add(getObject(filenames['1000_150_1000'], 'muHitEff_pt'))
    muHitEff_pt.Add(getObject(filenames['1000_150_10000'], 'muHitEff_pt'))
    tkHitEff_pt = getObject(filenames['1000_150_1'], 'tkHitEff_pt')
    tkHitEff_pt.Add(getObject(filenames['1000_150_10'], 'tkHitEff_pt'))
    tkHitEff_pt.Add(getObject(filenames['1000_150_100'], 'tkHitEff_pt'))
    tkHitEff_pt.Add(getObject(filenames['1000_150_1000'], 'tkHitEff_pt'))
    tkHitEff_pt.Add(getObject(filenames['1000_150_10000'], 'tkHitEff_pt'))
    ptResEff_pt = getObject(filenames['1000_150_1'], 'ptResEff_pt')
    ptResEff_pt.Add(getObject(filenames['1000_150_10'], 'ptResEff_pt'))
    ptResEff_pt.Add(getObject(filenames['1000_150_100'], 'ptResEff_pt'))
    ptResEff_pt.Add(getObject(filenames['1000_150_1000'], 'ptResEff_pt'))
    ptResEff_pt.Add(getObject(filenames['1000_150_10000'], 'ptResEff_pt'))

    plot_ = Canvas.Canvas('idEff_pt', 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1) 
    plot_.addHisto(idEff_pt, 'AP', 'Full ID efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(muHitEff_pt, 'P,SAME', 'Muon hit selection', 'pl', gcolors['red'], 1, 0, marker = 24)
    plot_.addHisto(tkHitEff_pt, 'P,SAME', 'Track hit selection', 'pl', gcolors['blue'], 1, 0, marker = 24)
    plot_.addHisto(ptResEff_pt, 'P,SAME', 'p_{T} resolution selection', 'pl', gcolors['green'], 1, 0, marker = 24)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    plot_.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.7, ymax = 1.1, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

    ### ID efficiency (vs eta)
    idEff_eta = getObject(filenames['1000_150_1'], 'idEff_eta')
    idEff_eta.Add(getObject(filenames['1000_150_10'], 'idEff_eta'))
    idEff_eta.Add(getObject(filenames['1000_150_100'], 'idEff_eta'))
    idEff_eta.Add(getObject(filenames['1000_150_1000'], 'idEff_eta'))
    idEff_eta.Add(getObject(filenames['1000_150_10000'], 'idEff_eta'))
    muHitEff_eta = getObject(filenames['1000_150_1'], 'muHitEff_eta')
    muHitEff_eta.Add(getObject(filenames['1000_150_10'], 'muHitEff_eta'))
    muHitEff_eta.Add(getObject(filenames['1000_150_100'], 'muHitEff_eta'))
    muHitEff_eta.Add(getObject(filenames['1000_150_1000'], 'muHitEff_eta'))
    muHitEff_eta.Add(getObject(filenames['1000_150_10000'], 'muHitEff_eta'))
    tkHitEff_eta = getObject(filenames['1000_150_1'], 'tkHitEff_eta')
    tkHitEff_eta.Add(getObject(filenames['1000_150_10'], 'tkHitEff_eta'))
    tkHitEff_eta.Add(getObject(filenames['1000_150_100'], 'tkHitEff_eta'))
    tkHitEff_eta.Add(getObject(filenames['1000_150_1000'], 'tkHitEff_eta'))
    tkHitEff_eta.Add(getObject(filenames['1000_150_10000'], 'tkHitEff_eta'))
    ptResEff_eta = getObject(filenames['1000_150_1'], 'ptResEff_eta')
    ptResEff_eta.Add(getObject(filenames['1000_150_10'], 'ptResEff_eta'))
    ptResEff_eta.Add(getObject(filenames['1000_150_100'], 'ptResEff_eta'))
    ptResEff_eta.Add(getObject(filenames['1000_150_1000'], 'ptResEff_eta'))
    ptResEff_eta.Add(getObject(filenames['1000_150_10000'], 'ptResEff_eta'))

    plot_ = Canvas.Canvas('idEff_eta', 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1) 
    plot_.addHisto(idEff_eta, 'AP', 'Full ID efficiency', 'pl', r.kBlack, 1, 0, marker = 20)
    plot_.addHisto(muHitEff_eta, 'P,SAME', 'Muon hit selection', 'pl', gcolors['red'], 1, 0, marker = 24)
    plot_.addHisto(tkHitEff_eta, 'P,SAME', 'Track hit selection', 'pl', gcolors['blue'], 1, 0, marker = 24)
    plot_.addHisto(ptResEff_eta, 'P,SAME', 'p_{T} resolution selection', 'pl', gcolors['green'], 1, 0, marker = 24)
    plot_.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    plot_.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    plot_.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    plot_.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    plot_.save(1, 0, 0, '','', ymin = 0.7, ymax = 1.1, outputDir = WORKPATH + 'output_signalefficiencies/', isPrivate = True, is2d = True)

