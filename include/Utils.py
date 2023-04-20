##################################################################################################
#                           ____       _                                                         #
#                          / ___| __ _| | __ _ _ __   __ _  __ _  ___                            #
#                         | |  _ / _` | |/ _` | '_ \ / _` |/ _` |/ _ \                           #
#                         | |_| | (_| | | (_| | |_) | (_| | (_| | (_) |                          # 
#                     _____\____|\__,_|_|\__,_| .__/ \__,_|\__, |\___/ _                         #  
#                    |  ___| __ __ _ _ __ ___ |_|____      |___/  _ __| | _                      #_
#                    | |_ | '__/ _` | '_ ` _ \ / _ \ \ /\ / / _ \| '__| |/ /                     #
#                    |  _|| | | (_| | | | | | |  __/\ V  V / (_) | |  |   <                      # 
#                    |_|  |_|  \__,_|_| |_| |_|\___| \_/\_/ \___/|_|  |_|\_\                     #
#                                                                                                # 
##################################################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time

import Sample as Sample
import Launcher as Launcher
import helper as helper
import Canvas as Canvas
import CutManager as CutManager

#### Cross sections of the models
XSECS = {}

XSECS['HSS(1000,450,10000)'] = 0.107e3
XSECS['HSS(1000,450,1000)'] = 0.107e3
XSECS['HSS(1000,450,100)'] = 0.107e3
XSECS['HSS(1000,450,10)'] = 0.107e3
XSECS['HSS(1000,450,1)'] = 0.107e3
XSECS['HSS(1000,350,10000)'] = 0.107e3
XSECS['HSS(1000,350,1000)'] = 0.107e3
XSECS['HSS(1000,350,100)'] = 0.107e3
XSECS['HSS(1000,350,10)'] = 0.107e3
XSECS['HSS(1000,350,1)'] = 0.107e3
XSECS['HSS(1000,250,10000)'] = 0.107e3
XSECS['HSS(1000,250,1000)'] = 0.107e3
XSECS['HSS(1000,250,100)'] = 0.107e3
XSECS['HSS(1000,250,10)'] = 0.107e3
XSECS['HSS(1000,250,1)'] = 0.107e3
XSECS['HSS(1000,150,10000)'] = 0.107e3
XSECS['HSS(1000,150,1000)'] = 0.107e3
XSECS['HSS(1000,150,100)'] = 0.107e3
XSECS['HSS(1000,150,10)'] = 0.107e3
XSECS['HSS(1000,150,1)'] = 0.107e3
XSECS['HSS(1000,50,10000)'] = 0.107e3
XSECS['HSS(1000,50,1000)'] = 0.107e3
XSECS['HSS(1000,50,100)'] = 0.107e3
XSECS['HSS(1000,50,10)'] = 0.107e3
XSECS['HSS(1000,50,1)'] = 0.107e3

XSECS['HSS(800,350,10000)'] = 0.3025e3
XSECS['HSS(800,350,1000)'] = 0.3025e3
XSECS['HSS(800,350,100)'] = 0.3025e3
XSECS['HSS(800,350,10)'] = 0.3025e3
XSECS['HSS(800,350,1)'] = 0.3025e3
XSECS['HSS(800,250,10000)'] = 0.3025e3
XSECS['HSS(800,250,1000)'] = 0.3025e3
XSECS['HSS(800,250,100)'] = 0.3025e3
XSECS['HSS(800,250,10)'] = 0.3025e3
XSECS['HSS(800,250,1)'] = 0.3025e3
XSECS['HSS(800,150,10000)'] = 0.3025e3
XSECS['HSS(800,150,1000)'] = 0.3025e3
XSECS['HSS(800,150,100)'] = 0.3025e3
XSECS['HSS(800,150,10)'] = 0.3025e3
XSECS['HSS(800,150,1)'] = 0.3025e3
XSECS['HSS(800,50,10000)'] = 0.3025e3
XSECS['HSS(800,50,1000)'] = 0.3025e3
XSECS['HSS(800,50,100)'] = 0.3025e3
XSECS['HSS(800,50,10)'] = 0.3025e3
XSECS['HSS(800,50,1)'] = 0.3025e3

XSECS['HSS(600,250,10000)'] = 1.17846e3
XSECS['HSS(600,250,1000)'] = 1.17846e3
XSECS['HSS(600,250,100)'] = 1.17846e3
XSECS['HSS(600,250,10)'] = 1.17846e3
XSECS['HSS(600,250,1)'] = 1.17846e3
XSECS['HSS(600,150,10000)'] = 1.17846e3
XSECS['HSS(600,150,1000)'] = 1.17846e3
XSECS['HSS(600,150,100)'] = 1.17846e3
XSECS['HSS(600,150,10)'] = 1.17846e3
XSECS['HSS(600,150,1)'] = 1.17846e3
XSECS['HSS(600,50,10000)'] = 1.17846e3
XSECS['HSS(600,50,1000)'] = 1.17846e3
XSECS['HSS(600,50,100)'] = 1.17846e3
XSECS['HSS(600,50,10)'] = 1.17846e3
XSECS['HSS(600,50,1)'] = 1.17846e3

XSECS['HSS(500,150,10000)'] = 2.5884e3
XSECS['HSS(500,150,1000)'] = 2.5884e3
XSECS['HSS(500,150,100)'] = 2.5884e3
XSECS['HSS(500,150,10)'] = 2.5884e3
XSECS['HSS(500,150,1)'] = 2.5884e3
XSECS['HSS(500,50,10000)'] = 2.5884e3
XSECS['HSS(500,50,1000)'] = 2.5884e3
XSECS['HSS(500,50,100)'] = 2.5884e3
XSECS['HSS(500,50,10)'] = 2.5884e3
XSECS['HSS(500,50,1)'] = 2.5884e3

XSECS['HSS(400,150,10000)'] = 5.0403e3
XSECS['HSS(400,150,1000)'] = 5.0403e3
XSECS['HSS(400,150,100)'] = 5.0403e3
XSECS['HSS(400,150,10)'] = 5.0403e3
XSECS['HSS(400,150,1)'] = 5.0403e3
XSECS['HSS(400,50,10000)'] = 5.0403e3
XSECS['HSS(400,50,1000)'] = 5.0403e3
XSECS['HSS(400,50,100)'] = 5.0403e3
XSECS['HSS(400,50,10)'] = 5.0403e3
XSECS['HSS(400,50,1)'] = 5.0403e3

XSECS['HSS(300,50,10000)']  = 4.93829e3
XSECS['HSS(300,50,1000)']  = 4.93829e3
XSECS['HSS(300,50,100)']  = 4.93829e3
XSECS['HSS(300,50,10)']  = 4.93829e3
XSECS['HSS(300,50,1)']  = 4.93829e3

XSECS['HSS(125,50,10000)']  = 21.459897e3
XSECS['HSS(125,50,1000)']  = 21.459897e3
XSECS['HSS(125,50,100)']  = 21.459897e3
XSECS['HSS(125,50,10)']  = 21.459897e3
XSECS['HSS(125,50,1)']  = 21.459897e3

XSECS['RPV(1500,494,10000)'] = 0.67
XSECS['RPV(1500,494,1000)'] = 0.67
XSECS['RPV(1500,494,100)'] = 0.67
XSECS['RPV(1500,494,10)'] = 0.67
XSECS['RPV(1500,494,1)'] = 0.67

XSECS['RPV(350,148,10000)'] = 10e3
XSECS['RPV(350,148,1000)'] = 10e3
XSECS['RPV(350,148,100)'] = 10e3
XSECS['RPV(350,148,10)'] = 10e3
XSECS['RPV(350,148,1)'] = 10e3



def makeClosureTest(lumi, name, hBKG_A, hBKG_B, ylog, tree, inputdir, labelA = 'A', labelB = 'B', xlabel = '', outpath = '', sys = 0.1, yshift = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', rmin = 0.0, rmax = 2.0):

    ### Get histograms
    luminosity = lumi

    ### Initialize cumulative histograms
    cumA = hBKG_A.Clone('cumA')
    cumB = hBKG_B.Clone('cumB')
    cumA.Reset()
    cumB.Reset()

    #cumA = r.TH1F('cumA', '', hBKG_A.GetNbinsX(), hBKG_A.GetXaxis().GetXmin(), hBKG_A.GetXaxis().GetXmax())
    #cumB = r.TH1F('cumB', '', hBKG_B.GetNbinsX(), hBKG_B.GetXaxis().GetXmin(), hBKG_B.GetXaxis().GetXmax())
    cumA.SetTitle(';min('+xlabel+');Number of events with '+xlabel+ ' > ('+xlabel+')_{min}')

    ### Set cumulative values
    for n in range(1, hBKG_A.GetNbinsX() + 1):
        valA = 0.0
        errA = 0.0
        valB = 0.0
        errB = 0.0
        for j in range(n, hBKG_A.GetNbinsX() + 1):
            valA = valA + hBKG_A.GetBinContent(j)
            errA = errA + hBKG_A.GetBinError(j)
            valB = valB + hBKG_B.GetBinContent(j)
            errB = errB + hBKG_B.GetBinError(j)
        cumA.SetBinContent(n, valA)
        cumA.SetBinError(n, errA)
        cumB.SetBinContent(n, valB)
        cumB.SetBinError(n, errB)

    ### Histogram tuning 
    cumA.SetMarkerStyle(20)
    cumA.SetMarkerSize(1)
    cumA.SetMarkerColor(r.kBlue)
    cumA.SetLineWidth(2)
    cumA.GetYaxis().SetTitleOffset(1.3)
    cumB.SetMarkerStyle(25)
    cumB.SetMarkerSize(1)
    cumB.SetLineWidth(2)
    cumB.SetMarkerColor(r.kRed)
    hBKG_A.SetMarkerStyle(20)
    hBKG_A.SetMarkerSize(1)
    hBKG_A.SetLineWidth(2)
    hBKG_A.SetMarkerColor(r.kBlue)
    hBKG_B.SetMarkerStyle(25)
    hBKG_B.SetMarkerSize(1)
    hBKG_B.SetLineWidth(2)
    hBKG_B.SetMarkerColor(r.kRed)


    ### Get maximum
    maxValA = cumA.GetMaximum()
    maxValB = cumB.GetMaximum()
    maxVal = max(maxValA, maxValB)

    ### Set Maximum
    if not ylog:
        cumA.SetMaximum(1.3*maxVal)
        cumB.SetMaximum(1.3*maxVal)
    else:
        cumA.SetMaximum(100.0*maxVal)
        cumB.SetMaximum(100.0*maxVal)
        cumA.SetMinimum(0.1)
        cumB.SetMaximum(0.1)

    ### Get maximum
    maxValA = hBKG_A.GetMaximum()
    maxValB = hBKG_B.GetMaximum()
    maxVal = max(maxValA, maxValB)

    ### Set Maximum
    if not ylog:
        hBKG_A.SetMaximum(1.3*maxVal)
        hBKG_B.SetMaximum(1.3*maxVal)
    else:
        hBKG_A.SetMaximum(100.0*maxVal)
        hBKG_B.SetMaximum(100.0*maxVal)
        hBKG_A.SetMinimum(0.1)
        hBKG_B.SetMaximum(0.1)

    ## Systematic bar
    hsys = False
    if sys:
        hsys = hBKG_A.Clone("sys_up")
        hsys.Reset()
        for n in range(1, hBKG_A.GetNbinsX() +1):
            hsys.SetBinContent(n, 1.0 )
            hsys.SetBinError(n, sys)

    ### Canvas object
    plot = Canvas.Canvas('closuretest_'+name, 'png,pdf', 0.45, 0.66, 0.65, 0.78, 1)
    plot.addHisto(cumA, 'P', labelA, 'p', r.kBlue, 1, 0)
    plot.addHisto(cumB, 'P, SAME', labelB, 'p', r.kRed, 1, 1)

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.75, 'Dielectron vertices', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.75, 'Dimuon vertices', font = 42)

    plot.addLatex(0.85, 0.56, 'Tail cumulative of x = #int^{#infty}_{x_{min}} N(x) dx', font = 42, align = 31)
    plot.addLatex(0.17, 0.81, extralabel, font = 42, align = 11, size = 0.045)
    plot.addLatex(0.9, 0.88, DATAlabel, font = 42, align = 31, size = 0.045)

    ### Save it
    #outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/ClosureTests_' + outtag + '/'
    outdir = outpath + '/'
    plot.saveRatio(1, tree.isData, ylog, luminosity, cumA, cumB, r_ymin = rmin, r_ymax = rmax, label="SR/BCR", hsys = hsys, outputDir = outdir, inProgress = True)

    ### Main comparison
    cplot = Canvas.Canvas('plaincomparison_'+name, 'png,pdf', 0.45, 0.66, 0.65, 0.78, 1)
    cplot.addHisto(hBKG_A, 'P', labelA, 'p', r.kBlue, 1, 0)
    cplot.addHisto(hBKG_B, 'P, SAME', labelB, 'p', r.kRed, 1, 1)
    if LLlabel == 'EE':
        cplot.addLatex(0.17, 0.75, 'Dielectron vertices', font = 42)
    if LLlabel == 'MM':
        cplot.addLatex(0.17, 0.75, 'Dimuon vertices', font = 42)
    cplot.addLatex(0.17, 0.81, extralabel, font = 42, align = 11, size = 0.045)
    cplot.addLatex(0.9, 0.88, DATAlabel, font = 42, align = 31, size = 0.045)
    cplot.saveRatio(1, tree.isData, ylog, luminosity, hBKG_A, hBKG_B, label="SR/BCR", hsys = hsys, outputDir = outdir, inProgress = True)




#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###     - xsec to normalize the signal given in fb
###
#####
#####################


def makeBlindedPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, treeBKG = False, treeBKGlabel = '', rebin = False, lines = [], line_ymax = False, xlabel = '', outtag = '', ymax = 0.0, LLlabel = '', DATAlabel = '', extralabel = '', xsec = False, xlog = False, text = False, outpath = ''):


    ### Get histograms
    luminosity = lumi

    ### Get background estimation from data
    hbkg_ = treeDATA.getLoopTH1F(inputdir, hname_bkg)

    ### Get background estimation from simulation
    if treeBKG:
        hbkgsim_ = treeBKG.getLoopTH1F(inputdir, hname_SI)

    ### rebinins:
    if type(rebin) != bool:
        if type(rebin) == int:
            hbkg = hbkg_.Rebin(rebin)
            if treeBKG:
                hbkgsim = hbkgsim_.Rebin(rebin)
        else:
            if len(rebin) > 1:
                hbkg = hbkg_.Rebin(len(rebin)-1, hbkg_.GetName() + '_rebined', rebin)
                if treeBKG:
                    hbkgsim = hbkgsim_.Rebin(len(rebin)-1, hbkgsim_.GetName() + '_rebined', rebin)
    else:
        hbkg = hbkg_.Clone()
        if treeBKG:
            hbkgsim = hbkgsim_.Clone()

    ### Set background histos style
    hbkg.SetFillColorAlpha(r.kCyan-6, 0.8)
    hbkg.SetLineColor(r.kCyan-2)
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    if treeBKG:
        hbkgsim.SetLineColor(r.kBlack)
        hbkgsim.SetLineStyle(10)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)
    xsecs = [0.107e3, 4.938e3, 2.588e3, 0.67, 10e3] # HARDCODED THIS NEEDS TO CHANGE

    for _i, _h in enumerate(hSIS.GetHists()):

        if type(rebin) != bool:
            if type(rebin) == int:
                _h2 = _h.Rebin(rebin)
            else:
                if len(rebin) > 1:
                    _h2 = _h.Rebin(len(rebin)-1, hbkg_.GetName() + '_rebined', rebin)
        else:
            _h2 = _h.Clone()

        s_histos.append(copy.deepcopy(_h2))
        print(_h2.GetName(), xsecs[_i])
        s_histos[-1].Scale(XSECS[_h2.GetTitle()])


    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])

    ### Set Maximum
    if not ylog:
        hbkg.SetMaximum(1.3*maxVal)
    else:
        hbkg.SetMaximum(100.0*maxVal)
        hbkg.SetMinimum(0.1)

    if ymax: hbkg.SetMaximum(ymax)


    ### Canvas object
    plot = Canvas.Canvas('Blinded_'+hname_bkg, 'png,pdf', 0.15, 0.65, 0.45, 0.89, 1, lsize = 0.03)
    if text:
        plot.addHisto(hbkg, 'HIST, TEXT', 'Background (predicted)', 'f', '', 1, 0)
    else:
        plot.addHisto(hbkg, 'HIST', 'Background (predicted)', 'f', '', 1, 0)

    ### Add signals:
    colors = [r.kRed, r.kOrange, r.kGreen+2, r.kBlue, r.kMagenta]
    for i,_h in enumerate(s_histos):

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        if 'HSS' in _h.GetTitle():
            legend = 'H#rightarrowSS (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        else:
            legend = 'RPV (%d GeV, %d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        if text:
            plot.addHisto(_h, 'HIST TEXT, SAME', legend, 'l', colors[i], 1, i+1) # Signal
        else:
            plot.addHisto(_h, 'HIST, SAME', legend, 'l', colors[i], 1, i+1) # Signal

    if treeBKG:
        plot.addHisto(hbkgsim, 'HIST, SAME', treeBKGlabel, 'l', '', 1, 2 + len(s_histos)) # Signal

    if LLlabel == 'EE':
        plot.addLatex(0.7, 0.86, 'e^{+}e^{-} channel', font = 42, size = 0.035)
    if LLlabel == 'MM':
        plot.addLatex(0.7, 0.86, '#mu^{+}#mu^{-} channel', font = 42, size = 0.035)

    ## Lines
    if not line_ymax:
        line_ymax = hbkg.GetMaximum()
    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, line_ymax, r.kBlack, 2)

    if extralabel:
        plot.addLatex(0.5, 0.6, extralabel, font = 42, align = 22, size = 0.035)

    ### Save it
    #outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/SRPlots_' + outtag + '/'
    outdir = outpath + '/SRPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, xlog = xlog, maxYnumbers = 4, is2d = True, inProgress = True)




########################################################################
####   Function to compare observed data and predicted data
#
#        * Originally developed for background validation tests
#


def makeBackgroundValidationPlot(name, lumi, hname_SR, hname_CR, ylog, treeDATA, inputdir, rebin = False, limit = 0.0, xlabel = '', outpath = False, yshift = 0.0, LLlabel = '', extralabel = '', xlog = False, sys = 0.0):


    ### Get histograms
    luminosity = lumi

    hSR_ = treeDATA.getLoopTH1F(inputdir, hname_SR)
    hCR_ = treeDATA.getLoopTH1F(inputdir, hname_CR)

    ### rebinins:
    if type(rebin) != bool and type(rebin) != int:
        if len(rebin) > 1:
            hSR = hSR_.Rebin(len(rebin)-1, hSR_.GetName() + '_rebined', rebin)
            hCR = hCR_.Rebin(len(rebin)-1, hCR_.GetName() + '_rebined', rebin)
    if type(rebin) == int:
        hSR = hSR_.Rebin(rebin)
        hCR = hCR_.Rebin(rebin)
    else:
        hSR = hSR_.Clone()
        hCR = hCR_.Clone()

    ### Blinding limits:
    if limit:
        for n in range(1, hSR.GetNbinsX()):
            if hSR.GetBinLowEdge(n) > limit: hSR.SetBinContent(n, 0.0)
            #if hCR_.GetBinLowEdge(n) > limit: hCR_.SetBinContent(n, 0.0)

    hCR.SetFillColorAlpha(r.kCyan-6, 0.8)
    #hCR.SetLineColor(r.kCyan-6) 
    hCR.SetLineWidth(0)
    hSR.SetMarkerStyle(20)
    hSR.SetMarkerSize(0.8)
    hSR.SetMarkerColor(r.kBlack)
    hSR.SetLineColor(r.kBlack)
    hSR.SetLineWidth(2)
    hCR.GetXaxis().SetTitleSize(0.045)
    hSR.GetXaxis().SetTitleSize(0.045)
    hCR.GetYaxis().SetTitleSize(0.045)
    hSR.GetYaxis().SetTitleSize(0.045)

    ### Get maximum
    maxValSR = hSR.GetMaximum()
    maxValCR = hCR.GetMaximum()
    maxVal = max([maxValSR, maxValCR])

    ### Set Maximum
    if not ylog:
        hCR.SetMaximum(1.3*maxVal)
    else:
        hCR.SetMaximum(100.0*maxVal)
        hCR.SetMinimum(0.1)

    ## Systematic bar
    hsys = False
    if sys:
        hsys = hSR.Clone("sys_up")
        hsys.Reset()
        for n in range(1, hSR.GetNbinsX() +1):
            hsys.SetBinContent(n, 1.0 )
            hsys.SetBinError(n, sys)

    ### Canvas object
    plot = Canvas.Canvas('BKGVal_'+name, 'png,pdf', 0.51, 0.65, 0.7, 0.77, 1)
    plot.addHisto(hCR, 'HIST', 'Background (Data-driven)', 'f', '', 1, 0)
    plot.addHisto(hSR, 'P, SAME', 'Data', 'p', '', 1, 1)
    plot.addLatex(0.17, 0.8, extralabel, font = 42)

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.74, 'Dielectron channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.74, 'Dimuon channel', font = 42)


    ### Save it
    if not outpath:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/BKGValidation/'
    else:
        outdir = outpath
    plot.saveRatio(1, 1, ylog, luminosity, hSR, hCR, r_ymin = 0.0, r_ymax = 2.0, label="Obs./Pred.", hsys = hsys, outputDir = outdir, xlog = xlog, inProgress = True)




########################################################################
####   Function to compare two histograms in the same Galapago Tree
#
#        * Originally developed for symmetry tests
#

def makeAgreementTest(lumi, hname1, hname2, ylog, tree, inputdir, label1, label2, labela, labelb, name, isData, sys = False, xlabel = '', outtag = '', yshift = 0.0, LLlabel = '', rebin = 0, ranges = [], rmin = 0.8, rmax = 1.2, maxY = False, outpath = False):

    ## Get the histogram
    histo1 = tree.getLoopTH1F(inputdir, hname1, doOF = False)
    histo2 = tree.getLoopTH1F(inputdir, hname2, doOF = False)

    ## Histogram tunning
    if rebin:
        histo1.Rebin(rebin)
        histo2.Rebin(rebin)

    if len(ranges) > 0:
        histo1.SetAxisRange(ranges[0], ranges[1], "X")
        histo2.SetAxisRange(ranges[0], ranges[1], "X")
        r_xmin = ranges[0]
        r_xmax = ranges[1]
    else:
        r_xmin = 0
        r_xmax = 0

    histo1.SetMaximum(1.6*histo1.GetMaximum())
    histo1.SetMinimum(0.0)
    histo1.SetMarkerStyle(20)
    histo2.SetMarkerStyle(20)

    ## Systematic bar
    hsys = False
    if sys:
        hsys = histo1.Clone("sys_up")
        hsys.Reset()
        for n in range(1, histo1.GetNbinsX() +1):
            hsys.SetBinContent(n, 1.0 )
            hsys.SetBinError(n, sys)

    if outpath:
        outdir = outpath
    else:
        outdir = WORKPATH + 'SymmetryResults/'

    plot = Canvas.Canvas(name, 'png,pdf', 0.15, 0.6, 0.45, 0.75, 1, lsize = 0.045)
    plot.addHisto(histo1, 'P', label1, 'p', r.kBlack, 1, 0)
    plot.addHisto(histo2, 'P,SAME', label2, 'p', r.kRed, 1, 1)
    plot.addLatex(0.18, 0.78, labela, font = 42, size = 0.045, align = 11)
    plot.addLatex(0.9, 0.88, labelb, font = 42, size = 0.045, align = 31)
    plot.saveRatio(1, isData, 0, str(lumi), histo1, histo2, r_ymin = rmin, r_ymax = rmax, r_xmin = r_xmin, r_xmax = r_xmax, label = 'Ratio', hsys = hsys, xlog = False, outputDir = outdir, maxYnumbers = maxY)

    return


##################################################################################################
#################
####   Function to launch each task to a queue by using Launcher.py
#
#        * It creates a new script and submission file each time that the 
#          function is called with arg queue set to True
#

def launchToQueue(funname, queue, name, outtag):

    launcher = Launcher.Launcher(script = os.path.dirname(os.path.abspath(__main__.__file__)) +'/'+ __main__.__file__, env = queue, ID = name, output = outtag)

    # create the new script file:
    _file = open(launcher.script, 'r')
    for _l in _file.readlines():
        if funname and 'name' and name in _l:
            raworder = _l
            break
    _file.close()

    order = ''
    for i,_e in enumerate(raworder.split(',')):
        if 'queue' in _e:
            order += _e.replace('True', 'False')
        else:
            order += _e
        if (i != len(raworder.split(',')) - 1):
            order += ','

    launcher.addOrder(order)
    launcher.launch()
    time.sleep(1.0)


##################################################################################################
#################
####   Function to Make a 2D histogram for a same tree 'tree' with 2 variables
#
#        * [Format]: var = 'y:x'
#

def make2DPlot(queue, lumi, var, name, nbinx, xmin, xmax, nbiny, ymin, ymax, xlabel, ylabel, tree, cuts, outtag = '', options = 'colz', normed = False):
    

    if queue:

        launchToQueue('make2DPlot', queue, name, outtag)

    else:

        h = tree.getTH2F(lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cuts, options, xlabel, ylabel)

        ### TH2F tunning
        h.GetZaxis().SetLabelSize(0.03)

        plot = Canvas.Canvas('hist_'+name, 'png', 0.34, 0.68, 0.9, 0.9, 2)
        plot.addHisto(h, options, '', '', False, 1, 0)

        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/' 
        plot.save(0,0, 0, lumi, '', outputDir = outdir, zlog = True)
        

##################################################################################################
#################
####   Function to Make a single 1D histogram for ONE tree 'tree'
#
#        * [Format]: var = 'x'
#

def make1DPlot(queue, lumi, var, name, nbinx, xmin, xmax, xlabel, tree, cuts, ylog = False, outtag = '', split = False, normed = False):
    

    if queue:

        launchToQueue('make1DPlot', queue, name, outtag)

    else:

        h = tree.getTH1F(lumi, name, var, nbinx, xmin, xmax, cuts, '', xlabel)
        h.SetFillColor(r.kAzure)
        h.SetLineColor(r.kBlack)

        plot = Canvas.Canvas('hist_'+name, 'png', 0.34, 0.68, 0.9, 0.9, 2)
        plot.addHisto(h, 'hist', '', 'f', '', 1, 0)

        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/1DPlots_' + outtag + '/' 
        plot.save(0,1, ylog, lumi, '', outputDir = outdir)


def makePlot(queue, lumi, var, name, nbin, xmin, xmax, xlabel, logx, treeMC, cuts, outtag = '', treeSI = False, treeDATA = False, LLlabel = '', normed = False):

    if queue:

        launchToQueue('makePlot', queue, name, outtag)

    else:

        hMCS = treeMC.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, cuts, "", xlabel)
        luminosity = lumi
        backgroundCounter = 0
    
        ### MC total histogram
        f = 0
        for _i, _h in enumerate(hMCS.GetHists()):
            if not f: hMC = copy.deepcopy(_h)
            else: hMC.Add(_h, 1.)
            f = 1
            backgroundCounter+=1
    
        hMC.SetMarkerStyle(20) # Auxiliar to save the ratio correctly 
    
        ### Signal histograms
        if treeSI:
    
            hSIS = treeSI.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, cuts, "", xlabel)
    
            s_histos = []
            for _i, _h in enumerate(hSIS.GetHists()):
                s_histos.append(copy.deepcopy(_h))
    
    
        # Normalization
        if normed:
            hMC.Scale(1.0/hMC.Integral())
            for _h in s_histos: _h.Scale(1.0/_h.Integral())
            luminosity = 1.0
    
        ### Get maximum
        maxValMC = hMC.GetMaximum()
        maxValSI = 0 if not treeSI else max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
        maxValDATA = 0 if not treeDATA else hDATA.GetMaximum()
        maxVal = max([maxValMC, maxValSI, maxValDATA])
    
        ### Set Maximum
        if not logx:
            hMCS.SetMaximum(1.3*maxVal)
            hMCS.SetMinimum(0.0)
            hMC.SetMaximum(1.3*maxVal)
            hMC.SetMinimum(0.0)
            if treeSI:
                for _h in s_histos: 
                    _h.SetMaximum(1.3*maxVal)
                    _h.SetMinimum(0.0)
            if treeDATA:
                hDATA.SetMaximum(1.3*maxVal)
                hDATA.SetMinimum(0.0)
        else:
            hMCS.SetMaximum(10.0*maxVal)
            hMCS.SetMinimum(0.1)
            hMC.SetMaximum(10.0*maxVal)
            hMC.SetMinimum(0.1)
            if treeSI:
                for _h in s_histos: 
                    _h.SetMaximum(10.0*maxVal)
                    _h.SetMinimum(0.1)
            if treeDATA:
                hDATA.SetMaximum(10.0*maxVal)
                hDATA.SetMinimum(0.1)
    
    
        ### Canvas object
        if treeDATA:
            plot = Canvas.Canvas('hist_'+name, 'png', 0.34, 0.68, 0.9, 0.9, 2)
        else:
            plot = Canvas.Canvas('hist_'+name, 'png', 0.55, 0.7, 0.9, 0.9, 1)
    
        if normed: plot.addHisto(hMC, 'HIST', '', 'l', r.kBlue, 1, 0) # Background
        else: plot.addStack(hMCS, 'HIST', 1, 0) # Background
    
        if treeSI:
            for i,_h in enumerate(s_histos):
                _h.SetLineWidth(2) # provisional
                plot.addHisto(_h, 'HIST, SAME', _h.GetTitle(), 'l', _h.GetFillColor(), 1, i+backgroundCounter) # Signal
    
        if treeDATA:
            plot.addHisto(hDATA, 'P, SAME', '', 'p', r.kBlack, 1, 4)
    
    
       ### Dilepton banner
        if LLlabel == 'EE':
            plot.addLatex(0.17, 0.86, 'e^{+}e^{-} channel')
        if LLlabel == 'MM':
            plot.addLatex(0.17, 0.86, '#mu^{+}#mu^{-} channel')
    
        ### Save it
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/Plots_' + outtag + '/'
        if treeDATA:
            plot.saveRatio(1, 0, logx, luminosity, hDATA, hMC, label="Data/MC", outputDir = outdir)
        else:
            plot.save(1, 0, logx, luminosity, '', outputDir = outdir)


        del plot
        return 


#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###     - xsec to normalize the signal given in fb
###
#####
#####################

def makeSignalPlot2D(name, lumi, hname_sig, zlog, treeSI, inputdir, rebin = False, lines = [], legend = '', xlabel = '', outtag = '', outdir = '', LLlabel = '', extralabel = '', xlog = False, ylog = False, text = False):


    ### Get histograms
    luminosity = lumi

    hsig = treeSI.getLoopTH2F(inputdir, hname_sig)

    hsig.GetXaxis().SetTitleSize(0.045)
    hsig.GetYaxis().SetTitleSize(0.045)
    
    hsig.GetYaxis().SetTitle('Dielectron regions')
    hsig.GetYaxis().SetTitleOffset(1.3)
    r.gStyle.SetPalette(r.kBird)


    ### Canvas object
    plot = Canvas.Canvas('SIOnly_'+name, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, ww = 650, hh = 600, lsize = 0.028)
    if text:
        plot.addHisto(hsig, 'COLZ, TEXT', '', '', '', 1, 0)
    else:
        plot.addHisto(hsig, 'COLZ', '', '', '', 1, 0)

    ### Extralabel
    plot.addLatex(0.13, 0.93, legend, font = 42, size = 0.032)
    plot.addLatex(0.4, 0.85, extralabel, font = 42, size = 0.03)

    for line in lines:
        plot.addLine(line[0], line[1], line[2], line[3], r.kRed, 2)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, zlog = zlog, xlog = xlog, is2d = False)



##################################################################################################
#################
####   Function to do the background validation closure test in Monte Carlo
####   with charge-symmetric background contributions estimated from data
#
#        * Input: MC and DATA trees, cuts 1 and 2
#

def makeBkgClosureTestInMC(queue, lumi, var, name, nbin, xmin, xmax, xlabel, ylog, treeMC, treeDATA, cuts1, cuts2, outtag = '', normed = False, yshift = 0.0):

    if queue:

        launchToQueue('makeBkgClosureTestInMC', queue, name, outtag)

    else:

        luminosity = lumi

        ### Get charged-symmetric background contributions from data:
        SScut = ''
        if 'EEBase' in var:
            SScut = 'IsoTrackSel_charge[ElectronCandidate_isotrackIdx[EEBase_idxA[EEBase_maxIxy]]]*IsoTrackSel_charge[ElectronCandidate_isotrackIdx[EEBase_idxB[EEBase_maxIxy]]] > 0'
        else:
            SScut = 'DGM_charge[DMDMBase_idxA[DMDMBase_maxIxy]]*DGM_charge[DMDMBase_idxB[DMDMBase_maxIxy]] > 0' 
        OScut = ''
        if 'EEBase' in var:
            OScut = 'IsoTrackSel_charge[ElectronCandidate_isotrackIdx[EEBase_idxA[EEBase_maxIxy]]]*IsoTrackSel_charge[ElectronCandidate_isotrackIdx[EEBase_idxB[EEBase_maxIxy]]] < 0'
        else:
            OScut = 'DGM_charge[DMDMBase_idxA[DMDMBase_maxIxy]]*DGM_charge[DMDMBase_idxB[DMDMBase_maxIxy]] < 0' 
        

        ### Cut definition:
        cutManager = CutManager.CutManager()
        SSRegion1 = cutManager.AddListB([cuts1, SScut]) # Same-sign in region 1
        OSRegion1 = cutManager.AddListB([cuts1, OScut]) # Opposite-sign in region 1
        SSRegion2 = cutManager.AddListB([cuts2, SScut]) # Same-sign in region 2
        OSRegion2 = cutManager.AddListB([cuts2, OScut]) # Opposite-sign in region 2

        ### Histogram definition
        hSS1 = treeDATA.getTH1F(lumi, 'hSS1_%s'%(name), var, nbin, xmin, xmax, SSRegion1, '', xlabel)
        hOS1 = treeMC.getTH1F(lumi, "hOS1_%s"%(name), var, nbin, xmin, xmax, OSRegion1, '', xlabel)
        hSS2 = treeDATA.getTH1F(lumi, 'hSS2_%s'%(name), var, nbin, xmin, xmax, SSRegion2, '', xlabel)
        hOS2 = treeMC.getTH1F(lumi, "hOS2_%s"%(name), var, nbin, xmin, xmax, OSRegion2, '', xlabel)

        ### Final histograms by adding QCD contribution
        hOS1.Add(hSS1)
        hOS2.Add(hSS2)

        ### Get maximum
        maxVal = max([hOS1.GetMaximum(), hOS2.GetMaximum()])
    
        ### Set Maximum
        if not ylog:
            hOS1.SetMaximum(1.3*maxVal)
            hOS2.SetMaximum(1.3*maxVal)
            hOS1.SetMinimum(0)
            hOS2.SetMinimum(0)
        else:
            hOS1.SetMaximum(10.0*maxVal)
            hOS2.SetMaximum(10.0*maxVal)
            hOS1.SetMinimum(0.1)
            hOS2.SetMinimum(0.1)
    
        ### Histogram tunning:
        hOS1.SetMarkerSize(1)
        hOS2.SetMarkerSize(1)
        hOS1.SetMarkerStyle(24)
        hOS2.SetMarkerStyle(24)

        ### Create canvas
        plot = Canvas.Canvas('hist_'+name, 'png', 0.3, 0.79, 0.8, 0.9, 1) 
        plot.addHisto(hOS1, 'P', 'Background estimation in SR (|#Delta#Phi| < #pi/2)', 'p', r.kBlue, 1, 0)
        plot.addHisto(hOS2, 'P, SAME', 'Background estimation in CR (|#Delta#Phi| > #pi/2)', 'p', r.kRed, 1, 1)
    
        ### Dilepton banner
        if 'EE' in var:
            plot.addLatex(0.65, 0.65, 'e^{+}e^{-} channel')
        if 'DMDM' in var:
            plot.addLatex(0.65, 0.65, '#mu^{+}#mu^{-} channel')
        

        ### Save it
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/ClosureTests_' + outtag + '/'
        plot.saveRatio(1, 0, ylog, luminosity, hOS1, hOS2, label="SR/CR", outputDir = outdir)

        del plot
        return 


##################################################################################################
#################
####   Function to compute the sensitivity histogram
#
#        * Input: Background [TH1F], Signal [TH1F]
#        * Returns: Sensitivity [TH1F], Background cumulative [TH1F], Signal cumulative [TH1F]
#

def computeSensitivity(hBkg, hSig, Scolor = False):

    Nbins = hBkg.GetNbinsX()

    # Histogram booking:
    cumBkg = r.TH1F('cumBkg', '', Nbins, hBkg.GetXaxis().GetXmin(), hBkg.GetXaxis().GetXmax()) # Background cumulative
    cumSig = r.TH1F('cumSig', '', Nbins, hSig.GetXaxis().GetXmin(), hSig.GetXaxis().GetXmax()) # Signal cumulative
    Sensitivity = r.TH1F('sens', '', Nbins, hSig.GetXaxis().GetXmin(), hSig.GetXaxis().GetXmax()) # Sensitivity

    # Scan over axis:
    for n in range(1, Nbins + 1):

        Svalue = hSig.Integral(n, Nbins)
        Bvalue = hBkg.Integral(n, Nbins)

        if (Svalue + Bvalue == 0): 
            sens = 0.0 # Avoid math error
        else:
            sens = Svalue/math.sqrt(Svalue + Bvalue)

        cumSig.SetBinContent(n, Svalue)
        cumBkg.SetBinContent(n, Bvalue)
        Sensitivity.SetBinContent(n, sens)

    # Histogram details:
    if Scolor:
        cumSig.SetLineColor(Scolor)
        Sensitivity.SetLineColor(Scolor)

    cumSig.SetLineWidth(2)
    cumBkg.SetLineWidth(2)
    Sensitivity.SetLineWidth(2)

    # Set Y axis labels:
    cumSig.GetYaxis().SetTitle('Yield')
    cumBkg.GetYaxis().SetTitle('Yield')
    Sensitivity.GetYaxis().SetTitle('S/#sqrt{S + B}')
    # Set X labels:
    cumSig.GetXaxis().SetTitle('min. ' + hSig.GetXaxis().GetTitle())
    cumBkg.GetXaxis().SetTitle('min. ' + hBkg.GetXaxis().GetTitle())
    Sensitivity.GetXaxis().SetTitle('min, ' + hSig.GetXaxis().GetTitle())

    return Sensitivity, cumBkg, cumSig
    

def makeSensitivity(queue, lumi, var, name, nbin, xmin, xmax, xlabel, logx, treeMC, treeSI, cuts, outtag = '', treeDATA = False, LLlabel = ''):

    if queue:

        launchToQueue('makeSensitivity', queue, name, outtag)

    else:
        ### Get background histogram
        hMC = treeMC.getTH1F(lumi, "hMC_%s"%(name), var, nbin, xmin, xmax, cuts, "", xlabel)
        print("GetEntries: ", hMC.GetEntries())
        hSIS = treeSI.getStack(lumi, "hMCS_%s"%(name), var, nbin, xmin, xmax, cuts, "", xlabel)
        luminosity = lumi

        ### Get signal histogram
        s_histos = []
        for _i, _h in enumerate(hSIS.GetHists()):
            s_histos.append(copy.deepcopy(_h))

        ### Get significance values
        plot = Canvas.Canvas('sensitivity_'+name, 'png', 0.5, 0.7, 0.9, 0.9, 1)
        significances = []
        signalYields = []

        for _i, _h in enumerate(s_histos):

            if (_i == 0):
                bh = r.TH1F('gr', '', hMC.GetNbinsX(), hMC.GetXaxis().GetXmin(), hMC.GetXaxis().GetXmax())

            sh = r.TH1F('sh', '', hMC.GetNbinsX(), hMC.GetXaxis().GetXmin(), hMC.GetXaxis().GetXmax())
            gh = r.TH1F('gh', '', hMC.GetNbinsX(), hMC.GetXaxis().GetXmin(), hMC.GetXaxis().GetXmax())
            
            for n in range(1, hMC.GetNbinsX() + 1):
                
                Svalue = _h.Integral(n, hMC.GetNbinsX())
                Bvalue = hMC.Integral(n, hMC.GetNbinsX())
                print('Signal: ', Svalue)
                print('Background: ', Bvalue)

                if (Svalue + Bvalue == 0): sig = 0
                else: sig = Svalue/math.sqrt(Svalue + Bvalue)
                print("Sensitivity: ", sig)

                gh.SetBinContent(n, sig)
                sh.SetBinContent(n, Svalue)
                if (_i == 0): bh.SetBinContent(n, Bvalue)

            gh.SetTitle(_h.GetTitle())
            gh.SetLineWidth(2)
            gh.GetYaxis().SetTitle('S/#sqrt{S + B}')
            gh.GetXaxis().SetTitle(xlabel)
            gh.SetLineColor(_h.GetFillColor())
            sh.SetTitle(_h.GetTitle())
            sh.SetLineWidth(2)
            sh.GetYaxis().SetTitle('S/#sqrt{S + B}')
            sh.GetXaxis().SetTitle(xlabel)
            sh.SetLineColor(_h.GetFillColor())
            significances.append(copy.deepcopy(gh))
            signalYields.append(copy.deepcopy(sh))

        ### Compute the maximum significance
        s_max = 0
        for s in significances:
            if s.GetMaximum() > s_max: s_max = s.GetMaximum()

        ### Compute the maximum sample yield
        y_max = 0
        for s in signalYields:
            if s.GetMaximum() > s_max: y_max = s.GetMaximum()
        if bh.GetMaximum() > y_max: y_max = bh.GetMaximum()

        ### Plot significances
        for _i,s in enumerate(significances):
            if _i == 0:
                if logx: 
                    s.SetMaximum(10.0*s_max)
                    s.SetMinimum(0.1)
                else: 
                    s.SetMaximum(1.4*s_max)
                    s.SetMinimum(0.0)

                plot.addHisto(s, 'l', s.GetTitle(), 'l', s.GetLineColor(), 1, _i)
            else:
                plot.addHisto(s, 'l, same', s.GetTitle(), 'l', s.GetLineColor(), 1, _i)

        ### Dilepton banner
        if LLlabel == 'EE':
            plot.addLatex(0.17, 0.86, 'e^{+}e^{-} channel')
        if LLlabel == 'MM':
            plot.addLatex(0.17, 0.86, '#mu^{+}#mu^{-} channel')

        ### Save it
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/sensitivity_' + outtag + '/'
        plot.save(1, 0, logx, luminosity, '', outputDir = outdir)

        ### yield plot
        yieldplot = Canvas.Canvas('yields_'+name, 'png', 0.5, 0.7, 0.9, 0.9, 1)
        bh.SetMaximum(10.0*y_max)
        bh.SetLineWidth(2)
        yieldplot.addHisto(bh, 'l', 'Background', 'l', r.kBlack, 1, 0)
        for _i,sy in enumerate(signalYields):
            yieldplot.addHisto(sy, 'l, same', sy.GetTitle(), 'l', sy.GetLineColor(), 1, _i +1)

        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/yields_' + outtag + '/' # Redefinition
        yieldplot.save(1, 0, logx, luminosity, '', outputDir = outdir)

        return

'''
Function to build the histogram of systematic errors
Arguments:
   - sys_errors: list with all the systematics to apply, e.g. [0.02, 0.1, 0.015]
'''
def makeSystematicsHist(sys_errors, hMC):
    hsys = hMC.Clone()
    for i in range(hsys.GetNbinsX() + 1):
        # compute MC systematic error
        error_values = 1 * np.array(sys_errors)
        band_error = np.linalg.norm(error_values)
        # Fill histogram
        hsys.SetBinContent(i,1)
        hsys.SetBinError(i,band_error)
    return hsys


