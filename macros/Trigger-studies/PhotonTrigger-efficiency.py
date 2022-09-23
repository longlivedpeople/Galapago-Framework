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

EOSPATH = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/'

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

def getHistoFromEff(eff):

    histo = eff.GetTotalHistogram().Clone()
    histo.Reset()

    for i in range(1, histo.GetNbinsX() + 1):
        error = max([eff.GetEfficiencyErrorLow(i), eff.GetEfficiencyErrorUp(i)])
        histo.SetBinContent(i, eff.GetEfficiency(i))
        histo.SetBinError(i, error)

    histo.SetMaximum(1.2) 
    histo.SetMinimum(0.0) 

    return histo


def getSFPlot(data, mc):

    total_data = data.GetTotalHistogram()
    h_sf = total_data.Clone(total_data.GetName() + '_sf')
    h_sfErr = total_data.Clone(total_data.GetName() + '_sfErr')
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



def passedMETTrigger(ev, year):

    passed = False

    if year == '2016':
        passed = ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_MET200 or ev.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight or ev.HLT_PFMET170_HBHECleaned or ev.HLT_PFMET300 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
    elif year == '2017':
        passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
    elif year == '2018':
        passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMET200_HBHE_BeamHaloCleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight

    return passed


def passedPhotonTrigger(ev, year):

    passed = False

    if year == '2016':
        passed = ev.HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 or ev.HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55
    elif year == '2017':
        passed = ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 or ev.HLT_DoublePhoton70
    elif year == '2018':
        passed = ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto or ev.HLT_DoublePhoton70

    return passed


################################################################################################################

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-e', '--era', action='store', type=str, dest='era', default='', help='Era')
    (opts, args) = parser.parse_args()


    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    r.gStyle.SetPalette(r.kBird)
    r.gStyle.SetPaintTextFormat("0.2f")


    ## Set the era and year
    era = opts.era
    if era != '2016APV' and era != '2016' and era != '2017' and era != '2018':
       print('The code is gonna crash -> solve it with exception')

    year = ''
    if '2016' in era:
        year = '2016'
    else:
        year = era

    print(">> Running trigger efficiency measurement for")
    print("     era: " + era)
    print("     year: " + year)


    #########################
    ####   Load sample   ####
    #########################

    Datasets_2016APV = []
    Datasets_2016APV.append('MET_EG_Run2016B_HIPM')
    Datasets_2016APV.append('MET_EG_Run2016C_HIPM')
    Datasets_2016APV.append('MET_EG_Run2016D_HIPM')
    Datasets_2016APV.append('MET_EG_Run2016E_HIPM')
    Datasets_2016APV.append('MET_EG_Run2016F_HIPM')

    Datasets_2016 = []
    Datasets_2016.append('MET_EG_Run2016F_noHIPM')
    Datasets_2016.append('MET_EG_Run2016G_noHIPM')
    Datasets_2016.append('MET_EG_Run2016H_noHIPM')

    Datasets_2017 = []
    Datasets_2017.append('MET_EG_Run2017B')
    Datasets_2017.append('MET_EG_Run2017C')
    Datasets_2017.append('MET_EG_Run2017D')
    Datasets_2017.append('MET_EG_Run2017E')
    Datasets_2017.append('MET_EG_Run2017F')

    Datasets_2018 = []
    Datasets_2018.append('MET_EG_Run2018A')
    Datasets_2018.append('MET_EG_Run2018B')
    Datasets_2018.append('MET_EG_Run2018C')
    Datasets_2018.append('MET_EG_Run2018D')



    MC_2016APV = ['TTTo2L2Nu_preVFP']
    MC_2016 = ['TTTo2L2Nu_postVFP']
    MC_2017 = ['TTTo2L2Nu_2017']
    MC_2018 = ['TTTo2L2Nu_2018']


    if era == '2016APV':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2016APV, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2016APV, 'DATA'), name = year, isdata = 1 )
    elif era == '2016':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2016, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2016, 'DATA'), name = year, isdata = 1 )
    elif era == '2017':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2017, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2017, 'DATA'), name = year, isdata = 1 )
    elif era == '2018':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2018, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2018, 'DATA'), name = year, isdata = 1 )

    #########################
    ####   Init plots    ####
    #########################

    pt1_bin = np.array([30., 60., 80., 100., 125., 150., 200., 250.])
    pt2_bin = np.array([30., 40., 60., 80., 100., 125., 200.])
    eta2_bin = np.array([-2.0, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.0])
    eta2_bin = np.array([-2.0, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.0])

    plot = {}
    plot['Efficiency_HLT_Full_pt_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_2d_DATA', ";Subleading electron E_{T};Leading electron E_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['Efficiency_HLT_Full_pt_2d_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_2d_MC', ";Subleading electron E_{T};Leading electron E_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['Efficiency_HLT_Full_eta_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_eta_2d_DATA', ";Subleading electron #eta;Leading electron #eta", len(eta2_bin)-1, eta2_bin, len(eta2_bin)-1, eta2_bin))
    plot['Efficiency_HLT_Full_eta_2d_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_eta_2d_MC', ";Subleading electron #eta;Leading electron #eta", len(eta2_bin)-1, eta2_bin, len(eta2_bin)-1, eta2_bin))
    plot['Efficiency_HLT_Full_pt_eta_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_eta_2d_DATA', ";Subleading electron E_{T};Subeading electron #eta", len(pt2_bin)-1, pt2_bin, len(eta2_bin)-1, eta2_bin))
    plot['Efficiency_HLT_Full_pt_eta_2d_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_eta_2d_MC', ";Subleading electron E_{T};Subleading electron #eta", len(pt2_bin)-1, pt2_bin, len(eta2_bin)-1, eta2_bin))
    plot['Efficiency_HLT_Full_pt_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_DATA', ";Subleading electron E_{T} (GeV) ;Efficiency", len(pt2_bin)-1, pt2_bin))
    plot['Efficiency_HLT_Full_pt_MC'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MC', ";Subleading electron E_{T} (GeV) ;Efficiency", len(pt2_bin)-1, pt2_bin))

    for p in plot.keys():
        r.SetOwnership(plot[p], 0)


    #########################
    ####   Cut values    ####
    #########################
    minPt = 30 if year == '2016' else 25
    minCosAlpha = -0.8 if year == '2016' else -0.99


    ###############################
    ####   Loop over events    ####
    ###############################

    ### Data:
    for b in treeDATA.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     if not (ev.nEE > 0):
                         continue

                     if not passedMETTrigger(ev, year):
                         continue

                     et_values = []
                     eta_values = []
                     indexes = []
                     for i in range(0, ev.nElectronCandidate):
                         if abs(ev.ElectronCandidate_eta[i]) > 2.:
                             continue
                         if ev.ElectronCandidate_relTrkiso[i] > 0.1:
                             continue
                         et_values.append(ev.ElectronCandidate_et[i])
                         eta_values.append(ev.ElectronCandidate_eta[i])
                         indexes.append(i)

                     if len(et_values) < 2:
                         continue

                     et_eta = zip(et_values, eta_values)
                     et_eta.sort(reverse = True)
                     et_ord = [x[0] for x in et_eta]
                     eta_ord = [x[1] for x in et_eta]


                     el1 = r.TLorentzVector()
                     el2 = r.TLorentzVector()
                     el1.SetPtEtaPhiM(ev.ElectronCandidate_et[indexes[0]], ev.ElectronCandidate_eta[indexes[0]], ev.ElectronCandidate_phi[indexes[0]], 105e-3)
                     el2.SetPtEtaPhiM(ev.ElectronCandidate_et[indexes[1]], ev.ElectronCandidate_eta[indexes[1]], ev.ElectronCandidate_phi[indexes[1]], 5e-6)

                     if et_values[0] < 40 or et_values[1] < 25:
                         continue
                     if (el1+el2).M() < 15.:
                         continue


                     plot['Efficiency_HLT_Full_pt_2d_DATA'].Fill(passedPhotonTrigger(ev, year), et_ord[1], et_ord[0])
                     plot['Efficiency_HLT_Full_eta_2d_DATA'].Fill(passedPhotonTrigger(ev, year), eta_ord[1], eta_ord[0])
                     plot['Efficiency_HLT_Full_pt_eta_2d_DATA'].Fill(passedPhotonTrigger(ev, year), et_ord[1], eta_ord[1])
                     plot['Efficiency_HLT_Full_pt_DATA'].Fill(passedPhotonTrigger(ev, year), et_ord[1])
                     
                     num += 1
                     if num > 1000:
                         num = 0
                         #break

    ### MC:
    for b in treeMC.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     if not (ev.nEE > 0):
                         continue

                     if not passedMETTrigger(ev, year):
                         continue

                     et_values = []
                     eta_values = []
                     indexes = []
                     for i in range(0, ev.nElectronCandidate):
                         if abs(ev.ElectronCandidate_eta[i]) > 2.:
                             continue
                         if ev.ElectronCandidate_relTrkiso[i] > 0.1:
                             continue
                         et_values.append(ev.ElectronCandidate_et[i])
                         eta_values.append(ev.ElectronCandidate_eta[i])
                         indexes.append(i)

                     if len(et_values) < 2:
                         continue

                     et_eta = zip(et_values, eta_values)
                     et_eta.sort(reverse = True)
                     et_ord = [x[0] for x in et_eta]
                     eta_ord = [x[1] for x in et_eta]


                     el1 = r.TLorentzVector()
                     el2 = r.TLorentzVector()
                     el1.SetPtEtaPhiM(ev.ElectronCandidate_et[indexes[0]], ev.ElectronCandidate_eta[indexes[0]], ev.ElectronCandidate_phi[indexes[0]], 105e-3)
                     el2.SetPtEtaPhiM(ev.ElectronCandidate_et[indexes[1]], ev.ElectronCandidate_eta[indexes[1]], ev.ElectronCandidate_phi[indexes[1]], 5e-6)

                     if et_values[0] < 40 or et_values[1] < 25:
                         continue
                     if (el1+el2).M() < 15.:
                         continue


                     plot['Efficiency_HLT_Full_pt_2d_MC'].Fill(passedPhotonTrigger(ev, year), et_ord[1], et_ord[0])
                     plot['Efficiency_HLT_Full_eta_2d_MC'].Fill(passedPhotonTrigger(ev, year), eta_ord[1], eta_ord[0])
                     plot['Efficiency_HLT_Full_pt_eta_2d_MC'].Fill(passedPhotonTrigger(ev, year), et_ord[1], eta_ord[1])
                     plot['Efficiency_HLT_Full_pt_MC'].Fill(passedPhotonTrigger(ev, year), et_ord[1])

                     num += 1
                     if num > 1000:
                          num = 0
                          #break

    ##################################################################################################

    outputFile = TFile(EOSPATH + 'PhotonTrigger-SFs/TH1F_photontrigger_'+era+'.root', 'RECREATE')
    for key in plot.keys():
        plot[key].Write()


    ### Efficiency 

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_Eff_full_pt_1D", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hMC_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MC'])
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hMC_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.addLatex(0.9, 0.88, era, size = 0.035, align = 31)
    canvas.saveRatio(1, 1, 0, '', hdata = hdata_, hMC = hMC_, r_ymin = 0.7, r_ymax = 1.0, label = 'Scale factor',outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False)

    r.gStyle.SetPadRightMargin(0.19)

    ### pt dependence

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_Data_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_2d_DATA'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_MC_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_2d_MC'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_pt_2d_DATA'], plot['Efficiency_HLT_Full_pt_2d_MC'])
    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_SF_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_SFErr_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    SF.Write()
    SFErr.Write()


    ### eta dependence

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_Data_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_eta_2d_DATA'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_MC_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_eta_2d_MC'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_eta_2d_DATA'], plot['Efficiency_HLT_Full_eta_2d_MC'])
    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_SF_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_SFErr_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    SF.Write()
    SFErr.Write()

    ### Full subleading muon dependence

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_Data_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_eta_2d_DATA'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_MC_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_eta_2d_MC'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_pt_eta_2d_DATA'], plot['Efficiency_HLT_Full_pt_eta_2d_MC'])
    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_SF_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("PhotonTrigger_"+era+"_SFErr_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'PhotonTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    SF.Write()
    SFErr.Write()

    r.gStyle.SetPadRightMargin(0.1)

    outputFile.Close()
