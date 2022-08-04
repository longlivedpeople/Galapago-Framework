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
    year = '2017'

    Datasets_2017 = []
    Datasets_2017.append('MET_EG_Run2017B')
    Datasets_2017.append('MET_EG_Run2017C')
    Datasets_2017.append('MET_EG_Run2017D')
    Datasets_2017.append('MET_EG_Run2017E')
    Datasets_2017.append('MET_EG_Run2017F')

    MC_2017 = ['TTTo2L2Nu_2017']

    if year == '2017':

        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2017, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2017, 'DATA'), name = year, isdata = 1 )

    #########################
    ####   Init plots    ####
    #########################

    et1_bin = np.array([40., 60., 80., 100., 125., 150., 200., 250.])
    et2_bin = np.array([25., 40., 60., 80., 100., 125., 150., 200.])

    plot = {}
    plot['Efficiency_HLT_Full_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_2d_DATA', ";Subleading electron E_{T};Leading electron E_{T}", len(et2_bin)-1, et2_bin, len(et1_bin)-1, et1_bin))
    plot['Efficiency_HLT_Full_2d_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_2d_MC', ";Subleading electron E_{T};Leading electron E_{T}", len(et2_bin)-1, et2_bin, len(et1_bin)-1, et1_bin))
    plot['Efficiency_HLT_path1_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_path1_2d_DATA', ";Subleading electron E_{T};Leading electron E_{T}", len(et2_bin)-1, et2_bin, len(et1_bin)-1, et1_bin))
    plot['Efficiency_HLT_path2_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_path2_2d_DATA', ";Subleading electron E_{T};Leading electron E_{T}", len(et2_bin)-1, et2_bin, len(et1_bin)-1, et1_bin))
    plot['Efficiency_HLT_path3_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_path3_2d_DATA', ";Subleading electron E_{T};Leading electron E_{T}", len(et2_bin)-1, et2_bin, len(et1_bin)-1, et1_bin))

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

                     if not (ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight):
                         continue

                     et_values = []
                     for i in range(0, ev.nElectronCandidate):
                         if abs(ev.ElectronCandidate_eta[i]) > 2.:
                             continue
                         et_values.append(ev.ElectronCandidate_et[i])
                     if len(et_values) < 2:
                         continue
                     et_values.sort(reverse = True)

                     plot['Efficiency_HLT_Full_2d_DATA'].Fill(ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 or ev.HLT_DoublePhoton70, et_values[1], et_values[0])
                     plot['Efficiency_HLT_path1_2d_DATA'].Fill(ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, et_values[1], et_values[0])
                     plot['Efficiency_HLT_path2_2d_DATA'].Fill(ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55, et_values[1], et_values[0])
                     plot['Efficiency_HLT_path3_2d_DATA'].Fill(ev.HLT_DoublePhoton70, et_values[1], et_values[0])
                     

    ### MC:
    num = 0
    for b in treeMC.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                for e,ev in enumerate(t):

                     if not (ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight):
                         continue

                     et_values = []
                     for i in range(0, ev.nElectronCandidate):
                         if abs(ev.ElectronCandidate_eta[i]) > 2.:
                             continue
                         et_values.append(ev.ElectronCandidate_et[i])
                     if len(et_values) < 2:
                         continue
                     et_values.sort(reverse = True)

                     num += 1
                     if num > 5000:
                          break

                     plot['Efficiency_HLT_Full_2d_MC'].Fill(ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 or ev.HLT_DoublePhoton70, et_values[1], et_values[0])


    ##################################################################################################

    ### Efficiency vs Lxy
    # (vertex efficiency on top of reco)
    r.gStyle.SetPadRightMargin(0.19)

    canvas = Canvas.Canvas("PhotonTrigger_2017Data_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_Full_2d_DATA'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_2017Data_path1_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_path1_2d_DATA'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_2017Data_path2_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_path2_2d_DATA'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_2017Data_path3_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_path3_2d_DATA'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_2017MC_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(plot['Efficiency_HLT_Full_2d_MC'],'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_2d_DATA'], plot['Efficiency_HLT_Full_2d_MC'])

    canvas = Canvas.Canvas("PhotonTrigger_2017SF_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("PhotonTrigger_2017SFErr_full_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, '2017 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    r.gStyle.SetPadRightMargin(0.1)


