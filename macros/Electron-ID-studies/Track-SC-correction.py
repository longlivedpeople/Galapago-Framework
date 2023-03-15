import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3, TLorentzVector, TMath
import math, sys, optparse, array, copy, os, json
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

EOSPATH = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/SCcorrection-studies/'

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

def createSysPlot(hist, sys):
    hsys = hist.Clone(hist.GetName() + '_sys')
    hsys.Reset()
    for i in range(1, hsys.GetNbinsX() + 1):
        hsys.SetBinContent(i, hist.GetBinContent(i))
        hsys.SetBinError(i, hist.GetBinContent(i)*sys)
    hsys.SetMarkerSize(0)
    hsys.SetFillColorAlpha(r.kBlue, 0.3)
    hsys.SetLineColorAlpha(r.kBlue, 0.3)
    return hsys

def passedMETTrigger(ev, year):

    passed = False

    if year == '2016':
        #passed = ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_MET200 or ev.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight or ev.HLT_PFMET170_HBHECleaned or ev.HLT_PFMET300 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        passed = ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_MET200 or ev.HLT_PFMET170_HBHECleaned or ev.HLT_PFMET300 
    elif year == '2017':
        passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
    elif year == '2018':
        #passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMET200_HBHE_BeamHaloCleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMET200_HBHE_BeamHaloCleaned

    return passed


def passedMuonTrigger(ev, year):

    passed = False

    if year == '2016':
        passed = ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10
    elif year == '2018':
        passed = ev.HLT_DoubleL2Mu23NoVtx_2Cha or ev.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed

    return passed


def combineEfficiency(effs, weights):

    final_total = effs[0].GetTotalHistogram().Clone(effs[0].GetName()+'_w')
    final_passed = effs[0].GetPassedHistogram().Clone(effs[0].GetName()+'_w')
    final_total.Sumw2()
    final_passed.Sumw2()
    final_total.Scale(weights[0])
    final_passed.Scale(weights[0])
    
    for i in range(1, len(effs) - 1):
        total = effs[i].GetTotalHistogram().Clone(effs[i].GetName()+'_w')
        passed = effs[i].GetPassedHistogram().Clone(effs[i].GetName()+'_w')
        total.Sumw2()
        passed.Sumw2()
        total.Scale(weights[i])
        passed.Scale(weights[i])
        final_total.Add(total)
        final_passed.Add(passed)

    final_eff = final_passed.Clone(final_passed.GetName() + '_combined')
    final_eff.Divide(final_total)
        
    return final_eff



################################################################################################################

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Era')
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


    lumi = {}
    lumi['2016APV'] = 19.7
    lumi['2016'] = 16.2
    lumi['2017'] = 41.5
    lumi['2018'] = 59.7

    #########################
    ####   Init plots    ####
    #########################

    pt_bin = np.array([30., 60., 125., 200.])
    eta_bin = np.array([-2., -1.56, -1.4, -0.9, 0.9, 1.4, 1.56, 2.])
    abseta_bin = np.array([0., 0.9, 1.42, 1.56, 2.])

    plot = {}
    plot['nPairs'] = copy.deepcopy(r.TH1F('nPairs', ";nPairs; Events", 6, 0, 6))
    plot['TH1F_probeIsoTrack_pt_total_DATA'] = copy.deepcopy(r.TH1F('TH1F_probeIsoTrack_pt_total_DATA', ";Isolated track probes;Isolated track p_{T} (GeV)", len(pt_bin)-1, pt_bin))
    plot['TH1F_probeIsoTrack_pt_passed_DATA'] = copy.deepcopy(r.TH1F('TH1F_probeIsoTrack_pt_passed_DATA', ";Isolated track probes;Isolated track p_{T} (GeV)", len(pt_bin)-1, pt_bin))
    plot['TH1F_probeIsoTrack_eta_total_DATA'] = copy.deepcopy(r.TH1F('TH1F_probeIsoTrack_eta_total_DATA', ";Isolated track probes;Isolated track #eta", len(eta_bin)-1, eta_bin))
    plot['TH1F_probeIsoTrack_eta_passed_DATA'] = copy.deepcopy(r.TH1F('TH1F_probeIsoTrack_eta_passed_DATA', ";Isolated track probes;Isolated track #eta", len(eta_bin)-1, eta_bin))
    plot['TH1F_probeIsoTrack_abseta_total_DATA'] = copy.deepcopy(r.TH1F('TH1F_probeIsoTrack_abseta_total_DATA', ";Isolated track probes;Isolated track #abseta", len(abseta_bin)-1, abseta_bin))
    plot['TH1F_probeIsoTrack_abseta_passed_DATA'] = copy.deepcopy(r.TH1F('TH1F_probeIsoTrack_abseta_passed_DATA', ";Isolated track probes;Isolated track #abseta", len(abseta_bin)-1, abseta_bin))
    plot['TH1F_tnp_mass_DATA'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_DATA', ";Pairs;Electron-track pair mass (GeV)", 50, 70, 120))
    plot['TH1F_tnp_mass_passed_eta1'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_eta1', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_eta1'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_eta1', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_eta2'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_eta2', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_eta2'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_eta2', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_eta3'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_eta3', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_eta3'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_eta3', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_eta4'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_eta4', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_eta4'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_eta4', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_pt1'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_pt1', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_pt1'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_pt1', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_pt2'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_pt2', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_pt2'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_pt2', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_pt3'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_pt3', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_pt3'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_pt3', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_pt4'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_pt4', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_pt4'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_pt4', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_passed_DATA'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_passed_DATA', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))
    plot['TH1F_tnp_mass_failed_DATA'] = copy.deepcopy(r.TH1F('TH1F_tnp_mass_failed_DATA', ";Pairs;Electron-track pair mass (GeV)", 60, 60, 120))

    for p in plot.keys():
        r.SetOwnership(plot[p], 0)

    ###############################
    ####   Loop over events    ####
    ###############################
    file_ = r.TFile(opts.input)
    tree_ = file_.Get("Events")
    if True:
         if True:
            if True:
                num = 0
                for e,ev in enumerate(tree_):

                     tnp_pairs = [] ## List of indexes [e, t] (electron ; track)
                     tnp_mass = [] ## List of mass

                     #### Electron channel
                     if (ev.HLT_Ele27_WPTight_Gsf):

                         for e in range(0, ev.nElectronCandidate):
                             #print(ev.ElectronCandidate_pt[e], ev.ElectronCandidate_et[e], abs(ev.ElectronCandidate_eta[e]), abs(ev.ElectronCandidate_dxy_PV[e]), abs(ev.IsoTrackSel_dz[ev.ElectronCandidate_isotrackIdx[e]]))
                             if ev.ElectronCandidate_pt[e] < 15: continue
                             if ev.ElectronCandidate_et[e] < 40: continue
                             if ev.ElectronCandidate_relTrkiso[e] > 0.1: continue
                             if abs(ev.ElectronCandidate_eta[e]) > 2: continue
                             if abs(ev.ElectronCandidate_dxy_PV[e]) > 0.2: continue
                             if abs(ev.IsoTrackSel_dz[ev.ElectronCandidate_isotrackIdx[e]]) > 0.5: continue
                             if ev.ElectronCandidate_dR[e] > 0.1: continue
                             if ev.PhotonSel_hadronicOverEm[ev.ElectronCandidate_photonIdx[e]] > 0.05: continue
                             if ev.PhotonSel_full5x5_sigmaIetaIeta[ev.ElectronCandidate_photonIdx[e]] > 0.0425 and abs(ev.ElectronCandidate_eta[e]) > 1.56: continue
                             if ev.PhotonSel_full5x5_sigmaIetaIeta[ev.ElectronCandidate_photonIdx[e]] > 0.0112 and abs(ev.ElectronCandidate_eta[e]) < 1.42: continue

                             for t in range(0, ev.nIsoTrack):
                                 #print(e, t, ev.IsoTrackSel_pt[t], ev.IsoTrackSel_eta[t], ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[e]]*ev.IsoTrackSel_charge[t])
                                 if t == ev.ElectronCandidate_isotrackIdx[e]: continue
                                 if ev.IsoTrackSel_pt[t] < 15: continue
                                 if abs(ev.IsoTrackSel_eta[t]) > 2: continue
                                 if ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[e]]*ev.IsoTrackSel_charge[t] > 0: continue

                                 tnp_pairs.append([e, t])

                                 tprobe = r.TLorentzVector()
                                 tprobe.SetPtEtaPhiM(ev.IsoTrackSel_pt[t], ev.IsoTrackSel_eta[t], ev.IsoTrackSel_phi[t], 0.0)
                                 ttag = r.TLorentzVector()
                                 ttag.SetPtEtaPhiM(ev.ElectronCandidate_et[e], ev.ElectronCandidate_eta[e], ev.ElectronCandidate_phi[e], 0.0)
                                 mass = (tprobe + ttag).M()                           
                                 tnp_mass.append(mass)
                                 #print(mass)

                                 plot['TH1F_tnp_mass_DATA'].Fill(mass)
                                 if mass > 81 and mass < 101:
                                     plot['TH1F_probeIsoTrack_pt_total_DATA'].Fill(ev.IsoTrackSel_pt[t])
                                     plot['TH1F_probeIsoTrack_eta_total_DATA'].Fill(ev.IsoTrackSel_eta[t])
                                     plot['TH1F_probeIsoTrack_abseta_total_DATA'].Fill(abs(ev.IsoTrackSel_eta[t]))
                                     if t in ev.ElectronCandidate_isotrackIdx:
                                         plot['TH1F_probeIsoTrack_pt_passed_DATA'].Fill(ev.IsoTrackSel_pt[t])
                                         plot['TH1F_probeIsoTrack_abseta_passed_DATA'].Fill(abs(ev.IsoTrackSel_eta[t]))
                                 
                         if (len(tnp_pairs) == 1):
                             if (tnp_mass[0] > 60 and tnp_mass[0] < 120):
                                 abseta_t = abs(ev.IsoTrackSel_eta[tnp_pairs[0][1]])
                                 pt_t = abs(ev.IsoTrackSel_pt[tnp_pairs[0][1]])
                                 if tnp_pairs[0][1] in ev.ElectronCandidate_isotrackIdx:
                                     plot['TH1F_tnp_mass_passed_DATA'].Fill(tnp_mass[0])
                                     ### abseta bins:
                                     if abseta_t < 0.9:
                                         plot['TH1F_tnp_mass_passed_eta1'].Fill(tnp_mass[0])
                                     elif abseta_t > 0.9 and abseta_t < 1.42:
                                         plot['TH1F_tnp_mass_passed_eta2'].Fill(tnp_mass[0])
                                     elif abseta_t > 1.42 and abseta_t < 1.56:
                                         plot['TH1F_tnp_mass_passed_eta3'].Fill(tnp_mass[0])
                                     elif abseta_t > 1.56 and abseta_t < 2.0:
                                         plot['TH1F_tnp_mass_passed_eta4'].Fill(tnp_mass[0])
                                     ### pt bins:
                                     if pt_t > 25 and pt_t < 60:
                                         plot['TH1F_tnp_mass_passed_pt1'].Fill(tnp_mass[0])
                                     elif pt_t > 60 and pt_t < 120:
                                         plot['TH1F_tnp_mass_passed_pt2'].Fill(tnp_mass[0])
                                     elif pt_t > 120 and pt_t < 200:
                                         plot['TH1F_tnp_mass_passed_pt3'].Fill(tnp_mass[0])
                                     elif pt_t > 200:
                                         plot['TH1F_tnp_mass_passed_pt4'].Fill(tnp_mass[0])
                                 else:
                                     plot['TH1F_tnp_mass_failed_DATA'].Fill(tnp_mass[0])
                                     ### abseta bins:
                                     if abseta_t < 0.9:
                                         plot['TH1F_tnp_mass_failed_eta1'].Fill(tnp_mass[0])
                                     elif abseta_t > 0.9 and abseta_t < 1.42:
                                         plot['TH1F_tnp_mass_failed_eta2'].Fill(tnp_mass[0])
                                     elif abseta_t > 1.42 and abseta_t < 1.56:
                                         plot['TH1F_tnp_mass_failed_eta3'].Fill(tnp_mass[0])
                                     elif abseta_t > 1.56 and abseta_t < 2.0:
                                         plot['TH1F_tnp_mass_failed_eta4'].Fill(tnp_mass[0])
                                     ### pt bins:
                                     if pt_t > 25 and pt_t < 60:
                                         plot['TH1F_tnp_mass_failed_pt1'].Fill(tnp_mass[0])
                                     elif pt_t > 60 and pt_t < 120:
                                         plot['TH1F_tnp_mass_failed_pt2'].Fill(tnp_mass[0])
                                     elif pt_t > 120 and pt_t < 200:
                                         plot['TH1F_tnp_mass_failed_pt3'].Fill(tnp_mass[0])
                                     elif pt_t > 200:
                                         plot['TH1F_tnp_mass_failed_pt4'].Fill(tnp_mass[0])


                             plot['nPairs'].Fill(len(tnp_pairs))


                     num += 1
                     if num > 100:
                         num = 0
                         #break


    ##################################################################################################


    ## Save
    inputname = opts.input.split('/')[-1]
    print(inputname)
    outputFile = TFile('/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/UL-analysis/CMSSW_10_6_20/src/MyAnalysis/Galapago-Framework/macros/Electron-ID-studies/plots_condor/'+inputname, 'RECREATE')
    for key in plot.keys():
        print('>> ' + key + ' saved!')
        plot[key].Write()

    outputFile.Close()


