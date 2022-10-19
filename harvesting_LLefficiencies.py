import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time
import shutil

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager
#from include.Utils import *

colors = {}
colors['125'] = {'1' : r.kGreen-9, '10' : r.kGreen-4, '100' : r.kGreen+2}
colors['400'] = {'1' : r.kRed-9, '10' : r.kRed-4, '100' : r.kRed+2}
colors['1000'] = {'1' : r.kBlue-9, '10' : r.kBlue-4, '100' : r.kBlue+2}


def makeEffPlotStack(lumi, name, hname, ylog, tree, inputdir, labels, outtag = '', LLlabel = '', era = 2016):

    luminosity = lumi

    _hBKG = tree.getLoopStack(inputdir, hname)
    hBKGtotal = tree.getLoopTH1F(inputdir, hname)

    ## Get maximum events
    nmax = hBKGtotal.GetBinContent(1)

 
    ## Geometric pad parameters
    leftMargin = 0.1
    topMargin = 0.13
    bottomMargin = 0.17
    rightMargin = 0.13
    r.gStyle.SetPadLeftMargin(leftMargin)   
    r.gStyle.SetPadTopMargin(topMargin)   
    r.gStyle.SetPadBottomMargin(bottomMargin)   
    r.gStyle.SetHistTopMargin(0.)
    r.gStyle.SetPadRightMargin(rightMargin)   
    binw = (1.0 - leftMargin - rightMargin)/hBKGtotal.GetNbinsX()
    x_labels = []
    for n in range(0, hBKGtotal.GetNbinsX()):
        x_labels.append(leftMargin + binw*(n + 0.5))


    if len(labels) != len(x_labels): return


    ## Init legend
    legend = r.TLegend(0.7, 0.7, 1.0 - rightMargin, 1.0 - topMargin)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(1)

    ## Declare normalize stack
    hBKG = r.THStack()
    hBKG.SetTitle(';; Selection efficiency')
    for _hist in _hBKG.GetHists():
        hist = _hist.Clone()
        for n in range(1, hist.GetNbinsX() +1 ): hist.SetBinContent(n, _hist.GetBinContent(n)/nmax)
        hBKG.Add(hist)
        legend.AddEntry(hist, hist.GetTitle(), 'f')


    ## Draw
    canvas = TCanvas('canvas', '', 700, 500)
    if ylog: canvas.SetLogy(1)
    canvas.cd()

    hBKG.SetMaximum(1.0)
    hBKG.SetMinimum(0.0)
    hBKG.Draw('')
    hBKG.GetXaxis().SetLabelSize(0)
    hBKG.GetXaxis().SetNdivisions(len(labels))
    hBKG.GetYaxis().SetTitleOffset(1.05)
    hBKG.SetMaximum(1.0)
    hBKG.SetMinimum(0.0)
    hBKG.Draw('hist')

    x_ticks = []
    for n in range(0, len(labels)):
        x_ticks.append(r.TLatex(x_labels[n], bottomMargin - 0.01, labels[n]))
        x_ticks[-1].SetNDC()
        x_ticks[-1].SetTextAlign(13)
        x_ticks[-1].SetTextSize(0.03)
        x_ticks[-1].SetTextFont(42)
        x_ticks[-1].SetTextAngle(-25)
        x_ticks[-1].Draw('same')
        
    ## Add legend
    legend.Draw()


    ## Make the histogram bonito
    CMS = r.TLatex(leftMargin, 1.0 - topMargin + 0.01, '#bf{CMS}')
    CMS.SetTextSize(0.075)
    CMS.SetTextAlign(11)
    CMS.SetTextFont(42)
    CMS.SetNDC()
    Simulation = r.TLatex(leftMargin + 0.11, 1.0 - topMargin + 0.01, '#it{Simulation}')
    Simulation.SetTextSize(0.05)
    Simulation.SetTextAlign(11)
    Simulation.SetTextFont(42)
    Simulation.SetNDC()
    if LLlabel == 'MM':
        year = r.TLatex(1.0 - rightMargin, 1.0 - topMargin + 0.01, 'Dimuons ({0})'.format(str(era)))
    elif LLlabel == 'EE':
        year = r.TLatex(1.0 - rightMargin, 1.0 - topMargin + 0.01, 'Dielectrons ({0})'.format(str(era))) 
    else:
        year = r.TLatex(1.0 - rightMargin, 1.0 - topMargin + 0.01, '({0})'.format(str(era)))
    year.SetTextSize(0.045)
    year.SetTextAlign(31)
    year.SetTextFont(42)
    year.SetNDC()

    CMS.Draw('Same')
    Simulation.Draw('Same')
    year.Draw('same')

    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/LLefficiencies_' + outtag + '/'
    if not os.path.exists(outdir): os.makedirs(outdir)
    canvas.SaveAs(outdir + name+'.png')
    canvas.SaveAs(outdir + name+'.pdf')



def makeEffPlotJoint(lumi, name, hname, ylog, tree, inputdir, labels, outtag = '', LLlabel = '', era = 2016, isSignal = True):

    luminosity = lumi

    _hSI = tree.getLoopStack(inputdir, hname)

    ## Geometric pad parameters
    leftMargin = 0.1
    topMargin = 0.13
    bottomMargin = 0.17
    rightMargin = 0.13
    r.gStyle.SetPadLeftMargin(leftMargin)   
    r.gStyle.SetPadTopMargin(topMargin)   
    r.gStyle.SetPadBottomMargin(bottomMargin)   
    r.gStyle.SetHistTopMargin(0.)
    r.gStyle.SetPadRightMargin(rightMargin)   
    binw = (1.0 - leftMargin - rightMargin)/len(labels)
    x_labels = []
    for n in range(0, len(labels)):
        x_labels.append(leftMargin + binw*(n + 0.5))

    if len(labels) != len(x_labels): return

    ## Init legend
    legend = r.TLegend(0.12, 0.18, 0.3, 0.4)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(1)

    ## Declare and draw normalized signal histos
    canvas = TCanvas('canvasSI', '', 700, 500)
    if ylog: canvas.SetLogy(1)
    canvas.cd()

    histos = []
    for _i, _hist in enumerate(_hSI.GetHists()):
        nmax = _hist.GetBinContent(1)
        histos.append(_hist.Clone())
        histos[-1].SetTitle(';; Selection efficiency')
        histos[-1].GetXaxis().SetLabelSize(0)
        histos[-1].GetYaxis().SetLabelSize(0.04)
        histos[-1].GetYaxis().SetTitleSize(0.05)
        histos[-1].GetXaxis().SetNdivisions(len(labels))
        histos[-1].GetYaxis().SetTitleOffset(1.05)
        histos[-1].SetMaximum(1.0)
        histos[-1].SetMinimum(0.0)
        histos[-1].SetFillColor(0)
        histos[-1].SetLineColor(_hist.GetFillColor())
        histos[-1].SetLineWidth(2)
        for n in range(1, histos[-1].GetNbinsX() +1 ): histos[-1].SetBinContent(n, _hist.GetBinContent(n)/nmax)
        if not _i: 
            histos[-1].Draw("hist")
        else:
            histos[-1].Draw("hist, same")

        if isSignal:
            masses = eval(_hist.GetTitle()[3:])
            legendtxt = 'm_{H} = '+str(masses[0])+' GeV, m_{X} = '+str(masses[1])+' GeV, c#tau = '+str(masses[2])+' mm'
            histos[-1].SetLineColor(colors[str(masses[0])][str(masses[2])])
        else: 
            legendtxt = _hist.GetTitle()

        legend.AddEntry(histos[-1], legendtxt, 'l')


    x_ticks = []
    for n in range(0, len(labels)):
        x_ticks.append(r.TLatex(x_labels[n], bottomMargin - 0.01, labels[n]))
        x_ticks[-1].SetNDC()
        x_ticks[-1].SetTextAlign(13)
        x_ticks[-1].SetTextSize(0.03)
        x_ticks[-1].SetTextFont(42)
        x_ticks[-1].SetTextAngle(-25)
        x_ticks[-1].Draw('same')
        
    ## Add legend
    legend.Draw()


    ## Make the histogram bonito
    CMS = r.TLatex(leftMargin, 1.0 - topMargin + 0.01, '#bf{CMS}')
    CMS.SetTextSize(0.075)
    CMS.SetTextAlign(11)
    CMS.SetTextFont(42)
    CMS.SetNDC()
    Simulation = r.TLatex(leftMargin + 0.11, 1.0 - topMargin + 0.01, '#it{Simulation}')
    Simulation.SetTextSize(0.05)
    Simulation.SetTextAlign(11)
    Simulation.SetTextFont(42)
    Simulation.SetNDC()
    if LLlabel == 'MM':
        year = r.TLatex(1.0 - rightMargin, 1.0 - topMargin + 0.01, 'Dimuons ({0})'.format(str(era)))
    elif LLlabel == 'EE':
        year = r.TLatex(1.0 - rightMargin, 1.0 - topMargin + 0.01, 'Dielectrons ({0})'.format(str(era))) 
    else:
        year = r.TLatex(1.0 - rightMargin, 1.0 - topMargin + 0.01, '({0})'.format(str(era)))
    year.SetTextSize(0.045)
    year.SetTextAlign(31)
    year.SetTextFont(42)
    year.SetNDC()

    CMS.Draw('Same')
    Simulation.Draw('Same')
    year.Draw('same')

    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/LLefficiencies_' + outtag + '/'
    if not os.path.exists(outdir): os.makedirs(outdir)
    canvas.SaveAs(outdir + name+'.png')
    canvas.SaveAs(outdir + name+'.pdf')


    

################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C+')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()

    ############# Background definition
    Backgrounds = []
    Backgrounds.append('DYJetsToLL_M-50') 
    Backgrounds.append('DYJetsToLL_M-10to50') 
    Backgrounds.append('WW') 
    Backgrounds.append('WZ') 
    Backgrounds.append('ZZ') 
    Backgrounds.append('TT') 

    ############# Signal definition
    Signals_2016 = []
    Signals_2016.append('HSS_400_50_1_2016')
    Signals_2016.append('HSS_400_50_10_2016')
    Signals_2016.append('HSS_400_50_100_2016')
    Signals_2016.append('HSS_1000_150_1_2016')
    Signals_2016.append('HSS_1000_150_10_2016')
    Signals_2016.append('HSS_1000_150_100_2016')
    Signals_2017 = []
    Signals_2017.append('HSS_400_50_1_2017')
    Signals_2017.append('HSS_400_50_10_2017')
    Signals_2017.append('HSS_400_50_100_2017')
    Signals_2017.append('HSS_1000_150_1_2017')
    Signals_2017.append('HSS_1000_150_10_2017')
    Signals_2017.append('HSS_1000_150_100_2017')
    Signals_2018 = []
    Signals_2018.append('HSS_400_50_1_2018')
    Signals_2018.append('HSS_400_50_10_2018')
    Signals_2018.append('HSS_400_50_100_2018')
    Signals_2018.append('HSS_1000_150_1_2018')
    Signals_2018.append('HSS_1000_150_10_2018')
    Signals_2018.append('HSS_1000_150_100_2018')

    filename = 'dat/Samples_cern_UltraLegacy.dat'
    #treeMC = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'MC'), name = 'MC', isdata = 0 )
    treeSI_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2016UL_Summer22.dat', Signals_2016, 'MC'), name = 'SI', isdata = 0 )
    treeSI_2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2017UL_Summer22.dat', Signals_2017, 'MC'), name = 'SI', isdata = 0 )
    treeSI_2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + 'dat/signals_2018UL_Summer22.dat', Signals_2018, 'MC'), name = 'SI', isdata = 0 )

    MMlabels = []
    MMlabels.append('All preselected vertices')
    MMlabels.append('Opposite charge')
    MMlabels.append('Min. p_{T}')
    MMlabels.append('Max. |#eta|')
    MMlabels.append('Isolation')
    MMlabels.append('Vertex quality')
    MMlabels.append('Cosmic muon rejection')
    MMlabels.append('Min. mass')
    MMlabels.append('|#Delta#Phi| < #pi/2')

    #makeEffPlotJoint(lumi_Muon, 'BKG_MMefficiency_2016', 'hMM_cutEfficiency', False, treeMC, WORKPATH + opts.input, MMlabels, outtag = 'LLefficiencies', LLlabel = 'MM', isSignal = False)
    makeEffPlotJoint(59.8, 'SI_MMefficiency_2016', 'hMM_efficiency', False, treeSI_2016, WORKPATH + opts.input, MMlabels, outtag = 'LLefficiencies', LLlabel = 'MM', era = 2016)
    makeEffPlotJoint(59.8, 'SI_MMefficiency_2018', 'hMM_efficiency', False, treeSI_2018, WORKPATH + opts.input, MMlabels, outtag = 'LLefficiencies', LLlabel = 'MM', era = 2018)


    EElabels = []
    EElabels.append('All preselected vertices')
    EElabels.append('Opposite charge')
    EElabels.append('Leading E_{T}')
    EElabels.append('Subleading E_{T}')
    EElabels.append('Max. |#eta|')
    EElabels.append('Isolation')
    EElabels.append('Vertex quality')
    EElabels.append('Min. mass')
    EElabels.append('|#Delta#Phi| < #pi/2')

    #makeEffPlotJoint(lumi_EG, 'BKG_EEefficiency_2016', 'hEE_cutEfficiency', False, treeMC, WORKPATH + opts.input, EElabels, outtag = 'LLefficiencies', LLlabel = 'EE', isSignal = False)
    makeEffPlotJoint(59.8, 'SI_EEefficiency_2016', 'hEE_efficiency', False, treeSI_2016, WORKPATH + opts.input, EElabels, outtag = 'LLefficiencies', LLlabel = 'EE', era = 2016)
    makeEffPlotJoint(59.8, 'SI_EEefficiency_2017', 'hEE_efficiency', False, treeSI_2017, WORKPATH + opts.input, EElabels, outtag = 'LLefficiencies', LLlabel = 'EE', era = 2017)
    makeEffPlotJoint(59.8, 'SI_EEefficiency_2018', 'hEE_efficiency', False, treeSI_2018, WORKPATH + opts.input, EElabels, outtag = 'LLefficiencies', LLlabel = 'EE', era = 2018)

