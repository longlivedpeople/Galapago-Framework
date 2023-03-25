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


