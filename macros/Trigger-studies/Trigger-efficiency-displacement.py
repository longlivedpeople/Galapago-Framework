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


def passedMuonTrigger(ev, year):

    passed = False

    if year == '2016':
        passed = ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10
    elif year == '2018':
        passed = ev.HLT_DoubleL2Mu23NoVtx_2Cha or ev.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed

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

    #year = '2018'

    HSS_Signals_2016 = []
    HSS_Signals_2016.append('HSS_125_30_1_2016')
    HSS_Signals_2016.append('HSS_125_30_10_2016')
    HSS_Signals_2016.append('HSS_125_30_100_2016')
    HSS_Signals_2016.append('HSS_125_30_1000_2016')
    HSS_Signals_2016.append('HSS_125_30_10000_2016')
    HSS_Signals_2016.append('HSS_400_50_1_2016')
    HSS_Signals_2016.append('HSS_400_50_10_2016')
    HSS_Signals_2016.append('HSS_400_50_100_2016')
    HSS_Signals_2016.append('HSS_400_50_1000_2016')
    HSS_Signals_2016.append('HSS_400_50_10000_2016')
    HSS_Signals_2016.append('HSS_1000_150_1_2016')
    HSS_Signals_2016.append('HSS_1000_150_10_2016')
    HSS_Signals_2016.append('HSS_1000_150_100_2016')
    HSS_Signals_2016.append('HSS_1000_150_1000_2016')
    HSS_Signals_2016.append('HSS_1000_150_10000_2016')

    HSS_Signals_2017 = []
    HSS_Signals_2017.append('HSS_125_30_1_2017')
    HSS_Signals_2017.append('HSS_125_30_10_2017')
    HSS_Signals_2017.append('HSS_125_30_100_2017')
    HSS_Signals_2017.append('HSS_125_30_1000_2017')
    HSS_Signals_2017.append('HSS_125_30_10000_2017')
    HSS_Signals_2017.append('HSS_400_50_1_2017')
    HSS_Signals_2017.append('HSS_400_50_10_2017')
    HSS_Signals_2017.append('HSS_400_50_100_2017')
    HSS_Signals_2017.append('HSS_400_50_1000_2017')
    HSS_Signals_2017.append('HSS_400_50_10000_2017')
    HSS_Signals_2017.append('HSS_1000_150_1_2017')
    HSS_Signals_2017.append('HSS_1000_150_10_2017')
    HSS_Signals_2017.append('HSS_1000_150_100_2017')
    HSS_Signals_2017.append('HSS_1000_150_1000_2017')
    HSS_Signals_2017.append('HSS_1000_150_10000_2017')

    HSS_Signals_2018 = []
    HSS_Signals_2018.append('HSS_125_30_1_2018')
    HSS_Signals_2018.append('HSS_125_30_10_2018')
    HSS_Signals_2018.append('HSS_125_30_100_2018')
    HSS_Signals_2018.append('HSS_125_30_1000_2018')
    HSS_Signals_2018.append('HSS_125_30_1000_2018')
    HSS_Signals_2018.append('HSS_400_50_1_2018')
    HSS_Signals_2018.append('HSS_400_50_10_2018')
    HSS_Signals_2018.append('HSS_400_50_100_2018')
    HSS_Signals_2018.append('HSS_400_50_1000_2018')
    HSS_Signals_2018.append('HSS_400_50_10000_2018')
    HSS_Signals_2018.append('HSS_1000_150_1_2018')
    HSS_Signals_2018.append('HSS_1000_150_10_2018')
    HSS_Signals_2018.append('HSS_1000_150_100_2018')
    HSS_Signals_2018.append('HSS_1000_150_1000_2018')
    HSS_Signals_2018.append('HSS_1000_150_10000_2018')

    sys_1000_150 = {}
    sys_1000_150['muon2016'] = 0.02
    sys_1000_150['muon2018'] = 0.03 
    sys_1000_150['photon2016'] = 0.01 
    sys_1000_150['photon2017'] = 0.01
    sys_1000_150['photon2018'] = 0.01
    sys_400_50 = {}
    sys_400_50['muon2016'] = 0.05
    sys_400_50['muon2018'] = 0.03 
    sys_400_50['photon2016'] = 0.05 
    sys_400_50['photon2017'] = 0.06
    sys_400_50['photon2018'] = 0.04
    sys_125_30 = {}
    sys_125_30['muon2016'] = 0.05
    sys_125_30['muon2018'] = 0.03 
    sys_125_30['photon2016'] = 0.15 
    sys_125_30['photon2017'] = 0.20
    sys_125_30['photon2018'] = 0.20
    


    #########################
    ####   Init plots    ####
    #########################
    dxy_bin = np.array([0.0, 0.5, 1., 2., 4.,6., 8., 10., 12., 14.])


    plot = {}
    if era == '2016':
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_2016UL_Summer22.dat', HSS_Signals_2016, 'DATA'), name = year, isdata = 1 )
        for signal in HSS_Signals_2016:
            key = signal.split('_')[1] + '_' + signal.split('_')[2]
            if key not in plot.keys():
                plot['Efficiency_muon_2d_dxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_muon_2d_dxy_' + key, ";Subleading muon |d_{dxy}| cm;Leading muon |d_{xy}| cm", len(dxy_bin )-1, dxy_bin, len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_photon_2d_dxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_photon_2d_dxy_' + key, ";Subleading electron |d_{dxy}|;Leading electron |d_{xy}| cm", len(dxy_bin )-1, dxy_bin, len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_muon_maxdxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_muon_maxdxy_' + key, "; Max. muon |d_{xy}| cm; Effiency", len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_photon_maxdxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_photon_maxdxy_' + key, ";Max. electron |d_{dxy}| cm;Efficiency", len(dxy_bin )-1, dxy_bin))
    elif era == '2017':
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_2017UL_Summer22.dat', HSS_Signals_2017, 'DATA'), name = year, isdata = 1 )
        for signal in HSS_Signals_2017:
            key = signal.split('_')[1] + '_' + signal.split('_')[2]
            if key not in plot.keys():
                plot['Efficiency_muon_2d_dxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_muon_2d_dxy_' + key, ";Subleading muon |d_{dxy}| cm;Leading muon |d_{xy}| cm", len(dxy_bin )-1, dxy_bin, len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_photon_2d_dxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_photon_2d_dxy_' + key, ";Subleading electron |d_{dxy}|;Leading electron |d_{xy}| cm", len(dxy_bin )-1, dxy_bin, len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_muon_maxdxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_muon_maxdxy_' + key, "; Max. muon |d_{xy}| cm; Effiency", len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_photon_maxdxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_photon_maxdxy_' + key, ";Max. electron |d_{dxy}| cm;Efficiency", len(dxy_bin )-1, dxy_bin))
    elif era == '2018':
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_2018UL_Summer22.dat', HSS_Signals_2018, 'DATA'), name = year, isdata = 1 )
        for signal in HSS_Signals_2018:
            key = signal.split('_')[1] + '_' + signal.split('_')[2]
            if key not in plot.keys():
                plot['Efficiency_muon_2d_dxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_muon_2d_dxy_' + key, ";Subleading muon |d_{dxy}| cm;Leading muon |d_{xy}| cm", len(dxy_bin )-1, dxy_bin, len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_photon_2d_dxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_photon_2d_dxy_' + key, ";Subleading electron |d_{dxy}|;Leading electron |d_{xy}| cm", len(dxy_bin )-1, dxy_bin, len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_muon_maxdxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_muon_maxdxy_' + key, "; Max. muon |d_{xy}| cm; Effiency", len(dxy_bin )-1, dxy_bin))
                plot['Efficiency_photon_maxdxy_' + key] = copy.deepcopy(r.TEfficiency('Efficiency_photon_maxdxy_' + key, ";Max. electron |d_{dxy}| cm;Efficiency", len(dxy_bin )-1, dxy_bin))


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

    ### MC:
    for b in treeMC.blocks:
        for s in b.samples: 
            key = s.name.split('_')[1] + '_' + s.name.split('_')[2]
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     ### Mu mu triggers:

                     if (ev.nDMDM > 0):

                         pt_values = []
                         eta_values = []
                         dxy_values = []
                         indexes = []
                         for i in range(0, ev.nDGM):
                             if abs(ev.DGM_eta[i]) > 2.:
                                 continue
                             if ev.DGM_relPFiso[i] > 0.2:
                                 continue
                             pt_values.append(ev.DGM_pt[i])
                             eta_values.append(ev.DGM_eta[i])
                             dxy_values.append(abs(ev.DGM_dxy_PV[i]))
                             indexes.append(i)

                         if len(pt_values) > 1:

                             pt_eta = zip(pt_values, eta_values)
                             pt_dxy = zip(pt_values, dxy_values)
                             pt_eta.sort(reverse = True)
                             pt_dxy.sort(reverse = True)
                             pt_ord = [x[0] for x in pt_eta]
                             eta_ord = [x[1] for x in pt_eta]
                             dxy_ord = [x[1] for x in pt_dxy]
                             dxy_ord.sort(reverse = True)

                             mu1 = r.TLorentzVector()
                             mu2 = r.TLorentzVector()
                             mu1.SetPtEtaPhiM(ev.DGM_pt[indexes[0]], ev.DGM_eta[indexes[0]], ev.DGM_phi[indexes[0]], 105e-3)
                             mu2.SetPtEtaPhiM(ev.DGM_pt[indexes[1]], ev.DGM_eta[indexes[1]], ev.DGM_phi[indexes[1]], 105e-3)

                             if pt_values[0] > minPt and pt_values[1] > minPt and (mu1+mu2).M() > 15. and math.cos(mu1.Angle(mu2.Vect())) > minCosAlpha and mu1.DeltaR(mu2) > 0.1:

                                 plot['Efficiency_muon_2d_dxy_' + key].Fill(passedMuonTrigger(ev, year), dxy_ord[1], dxy_ord[0])
                                 plot['Efficiency_muon_maxdxy_' + key].Fill(passedMuonTrigger(ev, year), dxy_ord[0])

                     if (ev.nEE > 0):

                         et_values = []
                         eta_values = []
                         dxy_values = []
                         indexes = []
                         for i in range(0, ev.nElectronCandidate):
                             if abs(ev.ElectronCandidate_eta[i]) > 2.:
                                 continue
                             if ev.ElectronCandidate_relTrkiso[i] > 0.1:
                                 continue
                             et_values.append(ev.ElectronCandidate_et[i])
                             eta_values.append(ev.ElectronCandidate_eta[i])
                             dxy_values.append(abs(ev.ElectronCandidate_dxy_PV[i]))
                             indexes.append(i)

                         if len(et_values) > 1:

                             et_eta = zip(et_values, eta_values)
                             et_dxy = zip(et_values, dxy_values)
                             et_eta.sort(reverse = True)
                             et_dxy.sort(reverse = True)
                             et_ord = [x[0] for x in et_eta]
                             eta_ord = [x[1] for x in et_eta]
                             dxy_ord = [x[1] for x in et_dxy]
                             dxy_ord.sort(reverse = True)


                             el1 = r.TLorentzVector()
                             el2 = r.TLorentzVector()
                             el1.SetPtEtaPhiM(ev.ElectronCandidate_et[indexes[0]], ev.ElectronCandidate_eta[indexes[0]], ev.ElectronCandidate_phi[indexes[0]], 105e-3)
                             el2.SetPtEtaPhiM(ev.ElectronCandidate_et[indexes[1]], ev.ElectronCandidate_eta[indexes[1]], ev.ElectronCandidate_phi[indexes[1]], 5e-6)

                             if et_values[0] > 40 and et_values[1] > 25 and (el1+el2).M() > 15.:
                                plot['Efficiency_photon_2d_dxy_' + key].Fill(passedPhotonTrigger(ev, year), dxy_ord[1], dxy_ord[0])
                                plot['Efficiency_photon_maxdxy_' + key].Fill(passedPhotonTrigger(ev, year), dxy_ord[0])




                     num += 1
                     if num > 1000:
                          num = 0
                     #     break



    ##################################################################################################

    outputFile = TFile(EOSPATH + 'SignalTrigger-SFs/TH1F_signal_'+era+'.root', 'RECREATE')
    for key in plot.keys():
        plot[key].Write()
        plot[key].Write()
        plot[key].GetTotalHistogram().Write()


    ### Efficiency 



    r.gStyle.SetPadRightMargin(0.19)

    ### pt dependence
    canvas = Canvas.Canvas("Efficiency_muon_2d_dxy_125_30_" + era, 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_muon_2d_dxy_125_30'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("Efficiency_photon_2d_dxy_125_30_" + era, 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_photon_2d_dxy_125_30'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("Efficiency_muon_2d_dxy_400_50_" + era, 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_muon_2d_dxy_400_50'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("Efficiency_photon_2d_dxy_400_50_" + era, 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_photon_2d_dxy_400_50'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("Efficiency_muon_2d_dxy_1000_150_" + era, 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_muon_2d_dxy_1000_150'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("Efficiency_photon_2d_dxy_1000_150_" + era, 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_photon_2d_dxy_1000_150'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False, is2d = True, labelz = 'Efficiency')


    r.gStyle.SetPadRightMargin(0.1)

    ## 1d histograms
    if year != '2017':

        canvas = Canvas.Canvas("MuonTrigger_maxdxy_" + era, 'png,pdf', 0.36, 0.72, 0.76, 0.89, 1)
        h_1000_150_muon = getHistoFromEff(plot['Efficiency_muon_maxdxy_1000_150'])
        e_1000_150_muon = h_1000_150_muon.GetBinContent(1)
        h_400_50_muon = getHistoFromEff(plot['Efficiency_muon_maxdxy_400_50'])
        e_400_50_muon = h_400_50_muon.GetBinContent(1)
        h_125_30_muon = getHistoFromEff(plot['Efficiency_muon_maxdxy_125_30'])
        e_125_30_muon = h_125_30_muon.GetBinContent(1)
        canvas.addHisto(h_1000_150_muon,'P', 'm_{H} = 1000 GeV, m_{S} = 150 GeV', 'pl', r.kBlack, True, 0, marker = 20)
        canvas.addHisto(h_400_50_muon,'P,SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV', 'pl', r.kBlue, True, 0, marker = 20)
        canvas.addHisto(h_125_30_muon,'P,SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV', 'pl', r.kRed, True, 0, marker = 20)
        canvas.addBand(0, e_1000_150_muon*(1-sys_1000_150['muon'+era]), 14, e_1000_150_muon*(1+sys_1000_150['muon'+era]), r.kBlack, 0.2)
        canvas.addLine(0, e_1000_150_muon, 14, e_1000_150_muon, r.kBlack)
        canvas.addBand(0, e_400_50_muon*(1-sys_400_50['muon'+era]), 14, e_400_50_muon*(1+sys_400_50['muon'+era]), r.kBlue, 0.2)
        canvas.addLine(0, e_400_50_muon, 14, e_400_50_muon, r.kBlue)
        canvas.addBand(0, e_125_30_muon*(1-sys_125_30['muon'+era]), 14, e_125_30_muon*(1+sys_125_30['muon'+era]), r.kRed, 0.2)
        canvas.addLine(0, e_125_30_muon, 14, e_125_30_muon, r.kRed)
        canvas.addLatex(0.9, 0.93, era, size = 0.035, align = 31)
        canvas.save(1, 0, 0, '', '', ymin=0.4, ymax=1.2, outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False)

    ## 1d histograms
    canvas = Canvas.Canvas("PhotonTrigger_maxdxy_" + era, 'png,pdf', 0.36, 0.72, 0.76, 0.89, 1)
    h_1000_150_photon = getHistoFromEff(plot['Efficiency_photon_maxdxy_1000_150'])
    e_1000_150_photon = h_1000_150_photon.GetBinContent(1)
    h_400_50_photon = getHistoFromEff(plot['Efficiency_photon_maxdxy_400_50'])
    e_400_50_photon = h_400_50_photon.GetBinContent(1)
    h_125_30_photon = getHistoFromEff(plot['Efficiency_photon_maxdxy_125_30'])
    e_125_30_photon = h_125_30_photon.GetBinContent(1)
    canvas.addHisto(h_1000_150_photon,'P', 'm_{H} = 1000 GeV, m_{S} = 150 GeV', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(h_400_50_photon,'P,SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV', 'pl', r.kBlue, True, 0, marker = 20)
    canvas.addHisto(h_125_30_photon,'P,SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV', 'pl', r.kRed, True, 0, marker = 20)
    canvas.addBand(0, e_1000_150_photon*(1-sys_1000_150['photon'+era]), 14, e_1000_150_photon*(1+sys_1000_150['photon'+era]), r.kBlack, 0.2)
    canvas.addLine(0, e_1000_150_photon, 14, e_1000_150_photon, r.kBlack)
    canvas.addBand(0, e_400_50_photon*(1-sys_400_50['photon'+era]), 14, e_400_50_photon*(1+sys_400_50['photon'+era]), r.kBlue, 0.2)
    canvas.addLine(0, e_400_50_photon, 14, e_400_50_photon, r.kBlue)
    canvas.addBand(0, e_125_30_photon*(1-sys_125_30['photon'+era]), 14, e_125_30_photon*(1+sys_125_30['photon'+era]), r.kRed, 0.2)
    canvas.addLine(0, e_125_30_photon, 14, e_125_30_photon, r.kRed)
    canvas.addLatex(0.9, 0.93, era, size = 0.035, align = 31)
    canvas.save(1, 0, 0, '', '', ymin=0.4, ymax=1.2, outputDir = EOSPATH + 'SignalTrigger-SFs/', inProgress = False)

    outputFile.Close()



