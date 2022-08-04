import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3, TLorentzVector, TMath
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

from include.galapagoStyle import sigpalette, gcolors, dcolors

def combinePlots(plotList):

    plot = plotList[0]
    for i in range(1, len(plotList)):
        plot.Add(plotList[i])

    return plot

def getSFPlot(data, mc):

    total_data = data.GetTotalHistogram()
    h_sf = total_data.Clone('sf')
    h_sfErr = total_data.Clone('sfErr')
    h_sf.Reset()
    h_sfErr.Reset()

    for i in range(1, total_data.GetNbinsX() + 1):
        for j in range(1, total_data.GetNbinsY() + 1):
            ijbin = data.GetGlobalBin(i, j)
            dataEff = data.GetEfficiency(ijbin)
            dataErrUp = data.GetEfficiencyErrorUp(ijbin)
            dataErrDown = data.GetEfficiencyErrorLow(ijbin)
            dataErr = max([dataErrUp, dataErrDown])
            mcEff = mc.GetEfficiency(ijbin)
            mcErrUp = mc.GetEfficiencyErrorUp(ijbin)
            mcErrDown = mc.GetEfficiencyErrorLow(ijbin)
            mcErr = max([mcErrUp, mcErrDown])
            if dataEff != 0 and mcEff != 0:
                sf = dataEff/mcEff
                sfErr = sf*((dataErr / dataEff)**2 + (mcErr / mcEff)**2)**0.5
            else:
                sf = 0.0
                sfErr = 0.0
            h_sf.SetBinContent(i, j, sf)
            h_sfErr.SetBinContent(i, j, sfErr)

    #h_sf.GetZaxis().SetTitle("Scale factor")
    #h_sfErr.GetZaxis().SetTitle("Scale factor uncertainty (stat)")

    return h_sf, h_sfErr


################################################################################################################

if __name__ == "__main__":

    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    r.gStyle.SetPalette(r.kBird)
    r.gStyle.SetPaintTextFormat("0.2f")


    #########################
    ####   Load sample   ####
    #########################
    year = '2016'

    Datasets_2016 = []
    Datasets_2016.append('MET_Muon_Run2016B_HIPM')
    Datasets_2016.append('MET_Muon_Run2016C_HIPM')
    Datasets_2016.append('MET_Muon_Run2016D_HIPM')
    Datasets_2016.append('MET_Muon_Run2016E_HIPM')
    Datasets_2016.append('MET_Muon_Run2016F_HIPM')
    Datasets_2016.append('MET_Muon_Run2016F_noHIPM')
    Datasets_2016.append('MET_Muon_Run2016G_noHIPM')
    Datasets_2016.append('MET_Muon_Run2016H_noHIPM')

    MC_2016 = ['TTTo2L2Nu_postVFP']

    if year == '2016':

        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2016, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2016, 'DATA'), name = year, isdata = 1 )

    #########################
    ####   Init plots    ####
    #########################

    pt1_bin = np.array([30., 60., 80., 100., 125., 150., 200., 250.])
    pt2_bin = np.array([30., 40., 60., 80., 100., 125., 150., 200.])

    plot = {}
    plot['Efficiency_HLT_Full_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_2d_DATA', ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['Efficiency_HLT_Full_2d_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_2d_MC', ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['Efficiency_HLT_Full_1d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_2d_DATA', ";Subleading muon p_{T} (GeV) ;Efficiency", len(pt2_bin)-1, pt2_bin))
    plot['Efficiency_HLT_Full_1d_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_2d_MC', ";Subleading muon p_{T} (GeV) ;Efficiency", len(pt2_bin)-1, pt2_bin))

    for p in plot.keys():
        r.SetOwnership(plot[p], 0)

    ###############################
    ####   Loop over events    ####
    ###############################

    ### Data:
    for b in treeDATA.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                for e,ev in enumerate(t):

                     if not (ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_MET200 or ev.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight or ev.HLT_PFMET170_HBHECleaned or ev.HLT_PFMET300 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight):
                         continue

                     pt_values = []
                     indexes = []
                     for i in range(0, ev.nDGM):
                         if abs(ev.DGM_eta[i]) > 2.:
                             continue
                         pt_values.append(ev.DGM_pt[i])
                         indexes.append(i)

                     if len(pt_values) < 2:
                         continue
                     pt_values.sort(reverse = True)


                     mu1 = r.TLorentzVector()
                     mu2 = r.TLorentzVector()
                     mu1.SetPtEtaPhiM(ev.DGM_pt[indexes[0]], ev.DGM_eta[indexes[0]], ev.DGM_phi[indexes[0]], 105e-3)
                     mu2.SetPtEtaPhiM(ev.DGM_pt[indexes[1]], ev.DGM_eta[indexes[1]], ev.DGM_phi[indexes[1]], 105e-3)

                     if pt_values[0] < 30 or pt_values[1] < 30:
                         continue
                     if (mu1+mu2).M() < 10.:
                         continue
                     if mu1.Angle(mu2.Vect()) > 2.5:
                         continue
                     if mu1.DeltaR(mu2) < 0.1:
                         continue


                     #print(ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, pt_values[0], pt_values[1], (mu1+mu2).M(), mu1.DeltaR(mu2), mu1.Angle(mu2.Vect()))

                     plot['Efficiency_HLT_Full_2d_DATA'].Fill(ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, pt_values[1], pt_values[0])
                     plot['Efficiency_HLT_Full_1d_DATA'].Fill(ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, pt_values[1])
                     

    ### MC:
    num = 0
    for b in treeMC.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                for e,ev in enumerate(t):

                     if not (ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_MET200 or ev.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight or ev.HLT_PFMET170_HBHECleaned or ev.HLT_PFMET300 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight):
                         continue

                     pt_values = []
                     indexes = []
                     for i in range(0, ev.nDGM):
                         if abs(ev.DGM_eta[i]) > 2.:
                             continue
                         pt_values.append(ev.DGM_pt[i])
                         indexes.append(i)

                     if len(pt_values) < 2:
                         continue
                     pt_values.sort(reverse = True)


                     mu1 = r.TLorentzVector()
                     mu2 = r.TLorentzVector()
                     mu1.SetPtEtaPhiM(ev.DGM_pt[indexes[0]], ev.DGM_eta[indexes[0]], ev.DGM_phi[indexes[0]], 105e-3)
                     mu2.SetPtEtaPhiM(ev.DGM_pt[indexes[1]], ev.DGM_eta[indexes[1]], ev.DGM_phi[indexes[1]], 105e-3)

                     if pt_values[0] < 30 or pt_values[1] < 30:
                         continue
                     if (mu1+mu2).M() < 10.:
                         continue
                     if mu1.Angle(mu2.Vect()) > 2.5:
                         continue
                     if mu1.DeltaR(mu2) < 0.1:
                         continue


                     num += 1
                     if num > 100000:
                          break

                     plot['Efficiency_HLT_Full_2d_MC'].Fill(ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, pt_values[1], pt_values[0])
                     plot['Efficiency_HLT_Full_1d_MC'].Fill(ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, pt_values[1])


    ##################################################################################################

    ### Efficiency vs Lxy
    # (vertex efficiency on top of reco)

    canvas = Canvas.Canvas("MuonTrigger_2017Data_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(plot['Efficiency_HLT_Full_2d_DATA'],'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(plot['Efficiency_HLT_Full_2d_MC'],'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.saveRatio(1, 1, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False)


    r.gStyle.SetPadRightMargin(0.19)

    canvas = Canvas.Canvas("MuonTrigger_2017Data_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_Full_2d_DATA'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("MuonTrigger_2017MC_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_Full_2d_MC'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_2d_DATA'], plot['Efficiency_HLT_Full_2d_MC'])
    canvas = Canvas.Canvas("MuonTrigger_2017SF_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("MuonTrigger_2017SFErr_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    r.gStyle.SetPadRightMargin(0.1)


