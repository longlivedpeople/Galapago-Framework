import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import operator
import time
import shutil

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager
#from include.Utils import *



def makeCutStability(name, hname, cut, ylog, tree, inputdir, outtag = '', LLlabel = '', era = 2016):


    _hSI = tree.getLoopStack(inputdir, hname)

    models = []

    for _i, _hist in enumerate(_hSI.GetHists()):

        sample = [] # eff; mH; mX; ctau;

        after = _hist.GetBinContent(cut)
        before = _hist.GetBinContent(cut - 1)
        sample.append(float(after/before))

        params = eval(_hist.GetTitle()[3:])
        sample.append(params[0])
        sample.append(params[1])
        sample.append(params[2])

        models.append(sample)

    ### Sort the models by mH > mX > ctau
    sorted_models = sorted(models, key=operator.itemgetter(1, 2, 3))

    print(sorted_models)

    ### Stability histogram
    shisto = r.TH1F("shisto", "", len(sorted_models), 0, len(sorted_models))
    shisto.GetXaxis().SetLabelSize(0)
    shisto.GetYaxis().SetTitle("Cut stability")
    shisto.GetXaxis().SetNdivisions(len(sorted_models))
    shisto.SetFillColorAlpha(r.kAzure-4, 0.4)
    shisto.SetMarkerColor(r.kRed+2)
    shisto.SetMarkerStyle(24)
    shisto.SetLineWidth(1)
    shisto.SetLineColor(r.kBlack)

    for _m, model in enumerate(sorted_models):
        shisto.SetBinContent(_m +1, model[0])

    

    ## Geometric pad parameters
    leftMargin = 0.13
    topMargin = 0.13
    bottomMargin = 0.15
    rightMargin = 0.13
    r.gStyle.SetPadLeftMargin(leftMargin)   
    r.gStyle.SetPadTopMargin(topMargin)   
    r.gStyle.SetPadBottomMargin(bottomMargin)   
    r.gStyle.SetHistTopMargin(0.)
    r.gStyle.SetPadRightMargin(rightMargin)   
    binw = (1.0 - leftMargin - rightMargin)/len(sorted_models)
    x_labels = []
    for n in range(0, len(sorted_models)):
        x_labels.append(leftMargin + binw*(n + 0.5))


    ## Declare and draw normalized signal histos
    canvas = TCanvas('canvasSI', '', 600, 500)
    if ylog: canvas.SetLogy(1)
    canvas.cd()
    

    ## Draw histogram
    y_max = 1.00
    y_min = 0.0
    shisto.SetMaximum(y_max)
    shisto.SetMinimum(y_min)
    shisto.Draw('hist')


    ## Draw x axis
    mH_level = 0.04
    mX_level = 0.08
    ctau_level = 0.12

    mH_heading = r.TLatex(0.01, mH_level, 'M_{H} (GeV):')
    mH_heading.SetNDC()
    mH_heading.SetTextAlign(13)
    mH_heading.SetTextSize(0.03)
    mH_heading.SetTextFont(42)
    mH_heading.Draw('same')
  
    mX_heading = r.TLatex(0.01, mX_level, 'M_{X} (GeV):')
    mX_heading.SetNDC()
    mX_heading.SetTextAlign(13)
    mX_heading.SetTextSize(0.03)
    mX_heading.SetTextFont(42)
    mX_heading.Draw('same')


    ctau_heading = r.TLatex(0.01, ctau_level, 'c#tau (mm):')
    ctau_heading.SetNDC()
    ctau_heading.SetTextAlign(13)
    ctau_heading.SetTextSize(0.03)
    ctau_heading.SetTextFont(42)
    ctau_heading.Draw('same')

    limit = r.TLine(0, 1.0, len(sorted_models), 1.0)
    limit.SetLineStyle(8)
    limit.Draw('same')


    mH_ticks = []
    mX_ticks = []
    ctau_ticks = []
    for n in range(0, len(x_labels)):
        mH_ticks.append(r.TLatex(x_labels[n], mH_level, str(sorted_models[n][1])))
        mH_ticks[-1].SetNDC()
        mH_ticks[-1].SetTextAlign(23)
        mH_ticks[-1].SetTextSize(0.03)
        mH_ticks[-1].SetTextFont(42)
        mH_ticks[-1].Draw('same')

        mX_ticks.append(r.TLatex(x_labels[n], mX_level, str(sorted_models[n][2])))
        mX_ticks[-1].SetNDC()
        mX_ticks[-1].SetTextAlign(23)
        mX_ticks[-1].SetTextSize(0.03)
        mX_ticks[-1].SetTextFont(42)
        mX_ticks[-1].Draw('same')

        ctau_ticks.append(r.TLatex(x_labels[n], ctau_level, str(sorted_models[n][3])))
        ctau_ticks[-1].SetNDC()
        ctau_ticks[-1].SetTextAlign(23)
        ctau_ticks[-1].SetTextSize(0.03)
        ctau_ticks[-1].SetTextFont(42)
        ctau_ticks[-1].Draw('same')


    ## Make the histogram bonito
    CMS = r.TLatex(leftMargin, 1.0 - topMargin + 0.01, '#bf{CMS}')
    CMS.SetTextSize(0.075)
    CMS.SetTextAlign(11)
    CMS.SetTextFont(42)
    CMS.SetNDC()
    Simulation = r.TLatex(leftMargin + 0.14, 1.0 - topMargin + 0.01, '#it{Simulation}')
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
    print('hello')
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
    parser.add_option('-e', '--EGinput', action='store', type=str, dest='EGinput', default='', help='Target directory')
    parser.add_option('-m', '--Muoninput', action='store', type=str, dest='Muoninput', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C+')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()
    #r.gStyle.SetErrorX(1)


    ############# Signal definition
    Signals = []
    Signals.append('HXX_1000_350_350mm')
    Signals.append('HXX_1000_350_35mm')
    Signals.append('HXX_1000_150_100mm')
    Signals.append('HXX_1000_150_10mm')
    Signals.append('HXX_400_150_400mm')
    Signals.append('HXX_400_50_400mm')
    Signals.append('HXX_400_50_40mm')
    Signals.append('HXX_400_50_4mm')

    ############# Luminosity definition
    lumiB = 5.79
    lumiC = 2.57
    lumiD = 4.25
    lumiE = 4.01
    lumiF = 3.10
    lumiG = 7.54
    lumiH = 8.61

    lumi_EG = lumiB + lumiC + lumiD + lumiE + lumiF + lumiG + lumiH
    lumi_Muon = lumiB + lumiC + lumiD + lumiE + lumiF + lumiG + lumiH


    filename = 'dat/Samples_cern_fillingv2.dat'
    treeSI = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Signals, 'MC'), name = 'SI', isdata = 0 )

    makeCutStability('SI_MMdPhiStability_2016', 'hMM_cutEfficiency', 12, False, treeSI, WORKPATH + opts.Muoninput, outtag = 'LLefficiencies', LLlabel = 'MM')
    makeCutStability('SI_EEdPhiStability_2016', 'hEE_cutEfficiency', 11, False, treeSI, WORKPATH + opts.Muoninput, outtag = 'LLefficiencies', LLlabel = 'EE')


