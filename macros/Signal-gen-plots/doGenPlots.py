import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
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

from include.galapagoStyle import sigpalette

################################################################################################################

if __name__ == "__main__":

    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    #########################
    ####   Load sample   ####
    #########################
    year = '2018'
    HSS_Signals = []
    HSS_Signals.append('HSS_125_30_1_'+year)
    HSS_Signals.append('HSS_125_30_10_'+year)
    HSS_Signals.append('HSS_125_30_100_'+year)
    HSS_Signals.append('HSS_400_50_1_'+year)
    HSS_Signals.append('HSS_400_50_10_'+year)
    HSS_Signals.append('HSS_400_50_100_'+year)
    HSS_Signals.append('HSS_1000_150_1_'+year)
    HSS_Signals.append('HSS_1000_150_10_'+year)
    HSS_Signals.append('HSS_1000_150_100_'+year)


    #########################
    ####   Init plots    ####
    #########################
    plot = {}
    for signal in HSS_Signals:
        plot[signal + '_muon_R']     = copy.deepcopy(r.TH1F(signal + '_muon_R', ";Generated LLP decay radius (cm);Number of generated muons", 100, 0, 100)) # in cm    
        plot[signal + '_electron_R'] = copy.deepcopy(r.TH1F(signal + '_electron_R', ";Generated LLP decay radius (cm);Number of generated electrons", 100, 0, 100)) # in cm    
        plot[signal + '_muon_dxy']     = copy.deepcopy(r.TH1F(signal + '_muon_dxy', ";Generated muon |d_{0}| (cm);Number of generated muons", 50, 0, 100)) # in cm    
        plot[signal + '_electron_dxy'] = copy.deepcopy(r.TH1F(signal + '_electron_dxy', ";Generated electron |d_{0}| (cm);Number of generated electrons", 50, 0, 100)) # in cm    
        plot[signal + '_muon_pt']     = copy.deepcopy(r.TH1F(signal + '_muon_pt', ";Generated muon p_{T} (GeV);Number of generated muons", 100, 0, 300)) # in cm    
        plot[signal + '_electron_pt'] = copy.deepcopy(r.TH1F(signal + '_electron_pt', ";Generated electron p_{T} (GeV);Number of generated electrons", 100, 0, 300)) # in cm    
        plot[signal + '_muon_vx_vy']     = copy.deepcopy(r.TH2F(signal + '_muon_vx_vy', ";Generated LLP vertex v_{x} (cm);Generated LLP vertex v_{y} (cm)", 100, -50, 50, 100, -50, 50)) # in cm    
        plot[signal + '_electron_vx_vy']     = copy.deepcopy(r.TH2F(signal + '_electron_vx_vy', ";Generated LLP vertex v_{x} (cm);Generated LLP vertex v_{y} (cm)", 100, -50, 50, 100, -50, 50)) # in cm    

        r.SetOwnership(plot[signal + '_muon_R'], 0)
        r.SetOwnership(plot[signal + '_electron_R'], 0)
        r.SetOwnership(plot[signal + '_muon_dxy'], 0)
        r.SetOwnership(plot[signal + '_electron_dxy'], 0)
        r.SetOwnership(plot[signal + '_muon_vx_vy'], 0)
        r.SetOwnership(plot[signal + '_electron_vx_vy'], 0)


        tree = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_'+year+'UL.dat', [signal], 'MC'), name = year, isdata = 0 )

        print("Processing " + signal)

        for b in tree.blocks: # assumed to be 1
            for s in b.samples: # assumed to be 1
                for t in s.ttrees:
                    for e,ev in enumerate(t):

                        for n in range(0, ev.nHardProcessParticle): 

                            lxy = math.sqrt(ev.HardProcessParticle_vx[n]*ev.HardProcessParticle_vx[n] + ev.HardProcessParticle_vy[n]*ev.HardProcessParticle_vy[n])
                            dxy = abs(-ev.HardProcessParticle_vx[n]*math.sin(ev.HardProcessParticle_phi[n]) + ev.HardProcessParticle_vy[n]*math.cos(ev.HardProcessParticle_phi[n]))
                            if abs(ev.HardProcessParticle_pdgId[n]) == 11:
                                plot[signal + '_electron_R'].Fill(lxy)
                                plot[signal + '_electron_dxy'].Fill(dxy)
                                plot[signal + '_electron_pt'].Fill(ev.HardProcessParticle_pt[n])
                                plot[signal + '_electron_vx_vy'].Fill(ev.HardProcessParticle_vx[n], ev.HardProcessParticle_vy[n])
                            if abs(ev.HardProcessParticle_pdgId[n]) == 13:
                                plot[signal + '_muon_R'].Fill(lxy)
                                plot[signal + '_muon_dxy'].Fill(dxy)
                                plot[signal + '_muon_pt'].Fill(ev.HardProcessParticle_pt[n])
                                plot[signal + '_muon_vx_vy'].Fill(ev.HardProcessParticle_vx[n], ev.HardProcessParticle_vy[n])


    ### Scale the plots
    for key in plot.keys():
        plot[key].Scale(plot[key].GetEntries() / plot[key].Integral())


    ### H->SS decay radius plots

    plot['HSS_400_50_1_2018_muon_dxy'].SetLineWidth(2)
    plot['HSS_400_50_10_2018_muon_dxy'].SetLineWidth(2)
    plot['HSS_400_50_100_2018_muon_dxy'].SetLineWidth(2)
    plot['HSS_1000_150_1_2018_muon_dxy'].SetLineWidth(2)
    plot['HSS_1000_150_10_2018_muon_dxy'].SetLineWidth(2)
    plot['HSS_1000_150_100_2018_muon_dxy'].SetLineWidth(2)
    plot['HSS_400_50_1_2018_muon_dxy'].SetMinimum(0.1)
    plot['HSS_400_50_10_2018_muon_dxy'].SetMinimum(0.1)
    plot['HSS_400_50_100_2018_muon_dxy'].SetMinimum(0.1)
    plot['HSS_1000_150_1_2018_muon_dxy'].SetMinimum(0.1)
    plot['HSS_1000_150_10_2018_muon_dxy'].SetMinimum(0.1)
    plot['HSS_1000_150_100_2018_muon_dxy'].SetMinimum(0.1)
    plot['HSS_400_50_1_2018_muon_dxy'].SetMaximum(1e6)
    plot['HSS_400_50_10_2018_muon_dxy'].SetMaximum(1e6)
    plot['HSS_400_50_100_2018_muon_dxy'].SetMaximum(1e6)
    canvas = Canvas.Canvas("HSS_400_50_muon_dxy", 'png,pdf', 0.33, 0.65, 0.63, 0.9, 1)
    canvas.addHisto(plot['HSS_400_50_1_2018_muon_dxy'],'HIST', 'H#rightarrowSS (400 GeV, 50 GeV, 1 mm)', 'l', r.kRed-9, True, 0, doOF = True)
    canvas.addHisto(plot['HSS_400_50_10_2018_muon_dxy'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 10 mm)', 'l', r.kRed-4, True, 1, doOF = True)
    canvas.addHisto(plot['HSS_400_50_100_2018_muon_dxy'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 100 mm)', 'l', r.kRed+2, True, 2, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_1_2018_muon_dxy'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 1 mm)', 'l', r.kBlue-9, True, 3, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_10_2018_muon_dxy'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 10 mm)', 'l', r.kBlue-4, True, 4, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_100_2018_muon_dxy'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 100 mm)', 'l', r.kBlue+2, True, 5, doOF = True)
    canvas.save(1, 0, 1, '', '', outputDir = WORKPATH + 'plots/', inProgress = False)

    plot['HSS_400_50_1_2018_electron_dxy'].SetLineWidth(2)
    plot['HSS_400_50_10_2018_electron_dxy'].SetLineWidth(2)
    plot['HSS_400_50_100_2018_electron_dxy'].SetLineWidth(2)
    plot['HSS_1000_150_1_2018_electron_dxy'].SetLineWidth(2)
    plot['HSS_1000_150_10_2018_electron_dxy'].SetLineWidth(2)
    plot['HSS_1000_150_100_2018_electron_dxy'].SetLineWidth(2)
    plot['HSS_400_50_1_2018_electron_dxy'].SetMinimum(0.1)
    plot['HSS_400_50_10_2018_electron_dxy'].SetMinimum(0.1)
    plot['HSS_400_50_100_2018_electron_dxy'].SetMinimum(0.1)
    plot['HSS_1000_150_1_2018_electron_dxy'].SetMinimum(0.1)
    plot['HSS_1000_150_10_2018_electron_dxy'].SetMinimum(0.1)
    plot['HSS_1000_150_100_2018_electron_dxy'].SetMinimum(0.1)
    plot['HSS_400_50_1_2018_electron_dxy'].SetMaximum(1e6)
    plot['HSS_400_50_10_2018_electron_dxy'].SetMaximum(1e6)
    plot['HSS_400_50_100_2018_electron_dxy'].SetMaximum(1e6)
    canvas = Canvas.Canvas("HSS_400_50_electron_dxy", 'png,pdf', 0.33, 0.65, 0.63, 0.9, 1)
    canvas.addHisto(plot['HSS_400_50_1_2018_electron_dxy'],'HIST', 'H#rightarrowSS (400 GeV, 50 GeV, 1 mm)', 'l', r.kRed-9, True, 0, doOF = True)
    canvas.addHisto(plot['HSS_400_50_10_2018_electron_dxy'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 10 mm)', 'l', r.kRed-4, True, 1, doOF = True)
    canvas.addHisto(plot['HSS_400_50_100_2018_electron_dxy'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 100 mm)', 'l', r.kRed+2, True, 2, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_1_2018_electron_dxy'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 1 mm)', 'l', r.kBlue-9, True, 3, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_10_2018_electron_dxy'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 10 mm)', 'l', r.kBlue-4, True, 4, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_100_2018_electron_dxy'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 100 mm)', 'l', r.kBlue+2, True, 5, doOF = True)
    canvas.save(1, 0, 1, '', '', outputDir = WORKPATH + 'plots/', inProgress = False)

    plot['HSS_400_50_10_2018_electron_pt'].SetLineWidth(2)
    plot['HSS_400_50_100_2018_electron_pt'].SetLineWidth(2)
    plot['HSS_125_30_10_2018_electron_pt'].SetLineWidth(2)
    plot['HSS_125_30_100_2018_electron_pt'].SetLineWidth(2)
    plot['HSS_1000_150_10_2018_electron_pt'].SetLineWidth(2)
    plot['HSS_1000_150_100_2018_electron_pt'].SetLineWidth(2)
    plot['HSS_125_30_10_2018_electron_pt'].SetMaximum(1e4)
    canvas = Canvas.Canvas("HSS_400_50_electron_pt", 'png,pdf', 0.33, 0.65, 0.63, 0.9, 1)
    canvas.addHisto(plot['HSS_125_30_10_2018_electron_pt'],'HIST', 'H#rightarrowSS (125 GeV, 30 GeV, 10 mm)', 'l', r.kGreen-4, True, 0, doOF = True)
    canvas.addHisto(plot['HSS_125_30_100_2018_electron_pt'],'HIST,SAME', 'H#rightarrowSS (125 GeV, 30 GeV, 100 mm)', 'l', r.kGreen+2, True, 1, doOF = True)
    canvas.addHisto(plot['HSS_400_50_10_2018_electron_pt'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 10 mm)', 'l', r.kRed-4, True, 2, doOF = True)
    canvas.addHisto(plot['HSS_400_50_100_2018_electron_pt'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 100 mm)', 'l', r.kRed+2, True, 3, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_10_2018_electron_pt'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 10 mm)', 'l', r.kBlue-4, True, 4, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_100_2018_electron_pt'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 100 mm)', 'l', r.kBlue+2, True, 5, doOF = True)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, maxYnumbers = 3)

    plot['HSS_400_50_10_2018_muon_pt'].SetLineWidth(2)
    plot['HSS_400_50_100_2018_muon_pt'].SetLineWidth(2)
    plot['HSS_125_30_10_2018_muon_pt'].SetLineWidth(2)
    plot['HSS_125_30_100_2018_muon_pt'].SetLineWidth(2)
    plot['HSS_1000_150_10_2018_muon_pt'].SetLineWidth(2)
    plot['HSS_1000_150_100_2018_muon_pt'].SetLineWidth(2)
    plot['HSS_125_30_10_2018_muon_pt'].SetMaximum(1e4)
    canvas = Canvas.Canvas("HSS_400_50_muon_pt", 'png,pdf', 0.33, 0.65, 0.63, 0.9, 1)
    canvas.addHisto(plot['HSS_125_30_10_2018_muon_pt'],'HIST', 'H#rightarrowSS (125 GeV, 30 GeV, 10 mm)', 'l', r.kGreen-4, True, 0, doOF = True)
    canvas.addHisto(plot['HSS_125_30_100_2018_muon_pt'],'HIST,SAME', 'H#rightarrowSS (125 GeV, 30 GeV, 100 mm)', 'l', r.kGreen+2, True, 1, doOF = True)
    canvas.addHisto(plot['HSS_400_50_10_2018_muon_pt'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 10 mm)', 'l', r.kRed-4, True, 2, doOF = True)
    canvas.addHisto(plot['HSS_400_50_100_2018_muon_pt'],'HIST,SAME', 'H#rightarrowSS (400 GeV, 50 GeV, 100 mm)', 'l', r.kRed+2, True, 3, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_10_2018_muon_pt'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 10 mm)', 'l', r.kBlue-4, True, 4, doOF = True)
    canvas.addHisto(plot['HSS_1000_150_100_2018_muon_pt'],'HIST,SAME', 'H#rightarrowSS (1000 GeV, 150 GeV, 100 mm)', 'l', r.kBlue+2, True, 5, doOF = True)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, maxYnumbers = 3)

    ### H->SS 2D plots
    r.gStyle.SetPalette(255, sigpalette)
    r.gStyle.SetPadRightMargin(0.12)

    #plot['HSS_400_50_100_2018_muon_vx_vy'].Scale(plot['HSS_400_50_100_2018_muon_vx_vy'].GetEntries()/plot['HSS_400_50_100_2018_muon_vx_vy'].Integral())
    canvas = Canvas.Canvas("HSS_400_50_100_muon_vx_vy", 'png,pdf', 0.35, 0.7, 0.75, 0.9, 1,  ww = 610, hh = 600)
    canvas.addHisto(plot['HSS_400_50_100_2018_muon_vx_vy'],'COLZ', '', '', r.kBlue+3, True, 0)
    canvas.addLatex(0.14, 0.93, '2018 UL      H#rightarrowSS#rightarrow#mu#mu+ll      (400 GeV, 50 GeV, 100 mm)', font = 42, size = 0.032)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, zlog = True)

    #plot['HSS_400_50_10_2018_muon_vx_vy'].Scale(plot['HSS_400_50_10_2018_muon_vx_vy'].GetEntries()/plot['HSS_400_50_10_2018_muon_vx_vy'].Integral())
    canvas = Canvas.Canvas("HSS_400_50_10_muon_vx_vy", 'png,pdf', 0.35, 0.7, 0.75, 0.9, 1,  ww = 610, hh = 600)
    canvas.addHisto(plot['HSS_400_50_10_2018_muon_vx_vy'],'COLZ', '', '', r.kBlue+3, True, 0)
    canvas.addLatex(0.14, 0.93, '2018 UL      H#rightarrowSS#rightarrow#mu#mu+ll      (400 GeV, 50 GeV, 10 mm)', font = 42, size = 0.032)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, zlog = True)

    #plot['HSS_400_50_100_2018_electron_vx_vy'].Scale(plot['HSS_400_50_100_2018_electron_vx_vy'].GetEntries()/plot['HSS_400_50_100_2018_electron_vx_vy'].Integral())
    canvas = Canvas.Canvas("HSS_400_50_100_electron_vx_vy", 'png,pdf', 0.35, 0.7, 0.75, 0.9, 1,  ww = 610, hh = 600)
    canvas.addHisto(plot['HSS_400_50_100_2018_electron_vx_vy'],'COLZ', '', '', r.kBlue+3, True, 0)
    canvas.addLatex(0.14, 0.93, '2018 UL      H#rightarrowSS#rightarrow#mu#mu+ll      (400 GeV, 50 GeV, 100 mm)', font = 42, size = 0.032)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, zlog = True)

    #plot['HSS_400_50_10_2018_electron_vx_vy'].Scale(plot['HSS_400_50_10_2018_electron_vx_vy'].GetEntries()/plot['HSS_400_50_10_2018_electron_vx_vy'].Integral())
    canvas = Canvas.Canvas("HSS_400_50_10_electron_vx_vy", 'png,pdf', 0.35, 0.7, 0.75, 0.9, 1,  ww = 610, hh = 600)
    canvas.addHisto(plot['HSS_400_50_10_2018_electron_vx_vy'],'COLZ', '', '', r.kBlue+3, True, 0)
    canvas.addLatex(0.14, 0.93, '2018 UL      H#rightarrowSS#rightarrow#mu#mu+ll      (400 GeV, 50 GeV, 10 mm)', font = 42, size = 0.032)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, zlog = True)
