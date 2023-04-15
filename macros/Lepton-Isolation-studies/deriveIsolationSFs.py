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

EOSPATH = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/Isolation-studies/WithHLTLoose/'

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

    ############# Muon data definition
    DoubleMuon2016APV = []
    DoubleMuon2016 = []
    DoubleMuon2018 = []
    DoubleMuon2016APV.append('DoubleMuon_Run2016B_HIPM')
    DoubleMuon2016APV.append('DoubleMuon_Run2016C_HIPM')
    DoubleMuon2016APV.append('DoubleMuon_Run2016D_HIPM')
    DoubleMuon2016APV.append('DoubleMuon_Run2016E_HIPM')
    DoubleMuon2016APV.append('DoubleMuon_Run2016F_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016F_noHIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016G_noHIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016H_noHIPM')
    DoubleMuon2018.append('DoubleMuon_Run2018A')
    DoubleMuon2018.append('DoubleMuon_Run2018B')
    DoubleMuon2018.append('DoubleMuon_Run2018C')
    DoubleMuon2018.append('DoubleMuon_Run2018D')

    DoubleEG2016APV = []
    DoubleEG2016 = []
    DoubleEG2017 = []
    DoubleEG2018 = []
    DoubleEG2016APV.append('DoubleEG_Run2016B_HIPM')
    DoubleEG2016APV.append('DoubleEG_Run2016C_HIPM')
    DoubleEG2016APV.append('DoubleEG_Run2016D_HIPM')
    DoubleEG2016APV.append('DoubleEG_Run2016E_HIPM')
    #DoubleEG2016APV.append('DoubleEG_Run2016F_HIPM')
    #DoubleEG2016.append('DoubleEG_Run2016F_noHIPM')
    DoubleEG2016.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG2016.append('DoubleEG_Run2016H_noHIPM')
    DoubleEG2017.append('DoubleEG_Run2017B')
    DoubleEG2017.append('DoubleEG_Run2017C')
    DoubleEG2017.append('DoubleEG_Run2017D')
    DoubleEG2017.append('DoubleEG_Run2017E')
    DoubleEG2017.append('DoubleEG_Run2017F')
    DoubleEG2018.append('EGamma_Run2018A')
    DoubleEG2018.append('EGamma_Run2018B')
    DoubleEG2018.append('EGamma_Run2018C')
    DoubleEG2018.append('EGamma_Run2018D')

    MC_2016APV = ['DYJetsToLL_M-50_preVFP']
    MC_2016 = ['DYJetsToLL_M-50_postVFP']
    MC_2017 = ['DYJetsToLL_M-50_2017']
    MC_2018 = ['DYJetsToLL_M-50_2018']

    if era == '2016APV':
        treeMuonDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Isolation.dat', DoubleMuon2016APV, 'DATA'), name = year, isdata = 1 )
        treeEGDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Isolation.dat', DoubleEG2016APV, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Isolation.dat', MC_2016APV, 'DATA'), name = year, isdata = 1 )
    elif era == '2016':
        treeMuonDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', DoubleMuon2016, 'DATA'), name = year, isdata = 1 )
        treeEGDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', DoubleEG2016, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', MC_2016, 'DATA'), name = year, isdata = 1 )
    elif era == '2017':
        treeEGDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', DoubleEG2017, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', MC_2017, 'DATA'), name = year, isdata = 1 )
    elif era == '2018':
        treeMuonDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', DoubleMuon2018, 'DATA'), name = year, isdata = 1 )
        treeEGDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', DoubleEG2018, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', MC_2018, 'DATA'), name = year, isdata = 1 )

    lumi = {}
    lumi['2016APV'] = 19.7
    lumi['2016'] = 16.2
    lumi['2017'] = 41.5
    lumi['2018'] = 59.7

    #########################
    ####   Init plots    ####
    #########################

    pt1_bin = np.array([0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 80., 100.])
    pt2_bin = np.array([30., 40., 60., 80., 100., 125., 200.])

    plot = {}
    plot['TH1F_DMDM_DGM_pt_total_DATA'] = copy.deepcopy(r.TH1F('TH1F_DMDM_DGM_pt_total_DATA', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_DMDM_DGM_pt_passedIso_DATA'] = copy.deepcopy(r.TH1F('TH1F_DMDM_DGM_pt_passedIso_DATA', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_DMDM_DGM_pt_total_MC'] = copy.deepcopy(r.TH1F('TH1F_DMDM_DGM_pt_total_MC', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_DMDM_DGM_pt_passedIso_MC'] = copy.deepcopy(r.TH1F('TH1F_DMDM_DGM_pt_passedIso_MC', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_EE_Electron_et_total_DATA'] = copy.deepcopy(r.TH1F('TH1F_EE_Electron_et_total_DATA', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_EE_Electron_et_passedIso_DATA'] = copy.deepcopy(r.TH1F('TH1F_EE_Electron_et_passedIso_DATA', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_EE_Electron_et_total_MC'] = copy.deepcopy(r.TH1F('TH1F_EE_Electron_et_total_MC', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))
    plot['TH1F_EE_Electron_et_passedIso_MC'] = copy.deepcopy(r.TH1F('TH1F_EE_Electron_et_passedIso_MC', ";Events;Muon p_{T}", len(pt1_bin)-1, pt1_bin))

    for p in plot.keys():
        r.SetOwnership(plot[p], 0)


    ############################
    ####   Config values    ####
    ############################
    cm = CutManager.CutManager()
    configfile = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/UL-analysis/CMSSW_10_6_20/src/MyAnalysis/Galapago-Framework/macros/Lepton-Isolation-studies/config_Spring23_Isolation.json'
    config = {}
    with open(configfile) as f:
        config = json.load(f)

    if era=='2016' or era=='2016APV':
        mumu_path      = cm.ORList(config['triggerPaths']['muons']['2016'], 'ev.')
        ee_path        = cm.ORList(config['triggerPaths']['electrons']['2016'], 'ev.')
        mumu_selection = cm.AddList(config['selection']['muons']['2016'])
        ee_selection   = cm.AddList(config['selection']['electrons']['2016'])
    elif era=='2017':
        ee_path        = cm.ORList(config['triggerPaths']['electrons']['2017'], 'ev.')
        ee_selection   = cm.AddList(config['selection']['electrons']['2017'])
    elif era=='2018':
        mumu_path      = cm.ORList(config['triggerPaths']['muons']['2018'], 'ev.')
        ee_path        = cm.ORList(config['triggerPaths']['electrons']['2018'], 'ev.')
        mumu_selection = cm.AddList(config['selection']['muons']['2018'])
        ee_selection   = cm.AddList(config['selection']['electrons']['2018'])


    ###############################
    ####   Loop over events    ####
    ###############################

    ### Data:
    if era != '2017':
        for b in treeMuonDATA.blocks:
            for s in b.samples: 
                for t in s.ttrees:
                    num = 0
                    for e,ev in enumerate(t):

                         #### Muon channel
                         if (eval(mumu_path) and ev.nDMDM == 1):

                             # Index for baseline selection application
                             imm = 0
                             # Basic selection to isolate On-Z region
                             if not eval(mumu_selection):
                                 continue
                             if abs(ev.DMDM_mass[0] - 91) > 10:
                                 continue
                             if ev.DGM_charge[ev.DMDM_idxA[0]]*ev.DGM_charge[ev.DMDM_idxB[0]] > 0:
                                 continue

                             plot['TH1F_DMDM_DGM_pt_total_DATA'].Fill(ev.DGM_pt[ev.DMDM_idxA[0]])
                             plot['TH1F_DMDM_DGM_pt_total_DATA'].Fill(ev.DGM_pt[ev.DMDM_idxB[0]])
                             if ev.DGM_relPFiso[ev.DMDM_idxA[0]] < 0.2:
                                 plot['TH1F_DMDM_DGM_pt_passedIso_DATA'].Fill(ev.DGM_pt[ev.DMDM_idxA[0]])
                             if ev.DGM_relPFiso[ev.DMDM_idxB[0]] < 0.2:
                                 plot['TH1F_DMDM_DGM_pt_passedIso_DATA'].Fill(ev.DGM_pt[ev.DMDM_idxB[0]])


                         num += 1
                         if num > 30:
                             num = 0
                             break

    for b in treeEGDATA.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     #### Electron channel
                     if (eval(ee_path) and ev.nEE == 1):

                         # Index for baseline selection application
                         iee = 0
                         # Basic selection to isolate On-Z region
                         if not eval(ee_selection):
                             continue
                         if abs(ev.EE_mass[0] - 91) > 10:
                             continue
                         if ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[0]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[0]]] > 0:
                             continue

                         if ev.EE_relisoA[0] < 0.1: 
                             plot['TH1F_EE_Electron_et_total_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]], ev.wPU)
                             if ev.EE_relisoB[0] < 0.1:
                                 plot['TH1F_EE_Electron_et_passedIso_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]], ev.wPU)
                         if ev.EE_relisoB[0] < 0.1: 
                             plot['TH1F_EE_Electron_et_total_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]], ev.wPU)
                             if ev.EE_relisoA[0] < 0.1:
                                 plot['TH1F_EE_Electron_et_passedIso_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]], ev.wPU)
                         """
                         plot['TH1F_EE_Electron_et_total_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]])
                         plot['TH1F_EE_Electron_et_total_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]])
                         if ev.EE_relisoA[0] < 0.1:
                             plot['TH1F_EE_Electron_et_passedIso_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]])
                         if ev.EE_relisoB[0] < 0.1:
                             plot['TH1F_EE_Electron_et_passedIso_DATA'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]])
                         """

                     num += 1
                     if num > 10000000:
                         num = 0
                         break


    ### MC:
    weights = {}
    for b in treeMC.blocks:
        for s in b.samples: 
            if era == '2016APV': 
                key = s.name.replace('_preVFP', '')
            elif era == '2016': 
                key = s.name.replace('_postVFP', '')
            elif era == '2017': 
                key = s.name.replace('_2017', '')
            elif era == '2018': 
                key = s.name.replace('_2018', '')
            weights[key] = s.lumWeight
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     #### Muon channel
                     if era != '2017':
                         if (eval(mumu_path) and ev.nDMDM == 1):

                             # Index for baseline selection application
                             imm = 0
                             # Basic selection to isolate On-Z region
                             if not eval(mumu_selection):
                                 continue
                             if abs(ev.DMDM_mass[0] - 91) > 10:
                                 continue
                             if ev.DGM_charge[ev.DMDM_idxA[0]]*ev.DGM_charge[ev.DMDM_idxB[0]] > 0:
                                 continue

                             plot['TH1F_DMDM_DGM_pt_total_MC'].Fill(ev.DGM_pt[ev.DMDM_idxA[0]], ev.wPU)
                             plot['TH1F_DMDM_DGM_pt_total_MC'].Fill(ev.DGM_pt[ev.DMDM_idxB[0]], ev.wPU)
                             if ev.DGM_relPFiso[ev.DMDM_idxA[0]] < 0.2:
                                 plot['TH1F_DMDM_DGM_pt_passedIso_MC'].Fill(ev.DGM_pt[ev.DMDM_idxA[0]], ev.wPU)
                             if ev.DGM_relPFiso[ev.DMDM_idxB[0]] < 0.2:
                                 plot['TH1F_DMDM_DGM_pt_passedIso_MC'].Fill(ev.DGM_pt[ev.DMDM_idxB[0]], ev.wPU)

                     #### Electron channel
                     if (eval(ee_path) and ev.nEE == 1):
                     #if (ev.nEE == 1):

                         # Index for baseline selection application
                         iee = 0
                         # Basic selection to isolate On-Z region
                         if not eval(ee_selection):
                             continue
                         if abs(ev.EE_mass[0] - 91) > 10:
                             continue
                         if ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[0]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[0]]] > 0:
                             continue

                         if ev.EE_relisoA[0] < 0.1: 
                             plot['TH1F_EE_Electron_et_total_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]], ev.wPU)
                             if ev.EE_relisoB[0] < 0.1:
                                 plot['TH1F_EE_Electron_et_passedIso_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]], ev.wPU)
                         if ev.EE_relisoB[0] < 0.1: 
                             plot['TH1F_EE_Electron_et_total_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]], ev.wPU)
                             if ev.EE_relisoA[0] < 0.1:
                                 plot['TH1F_EE_Electron_et_passedIso_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]], ev.wPU)

                         """
                         plot['TH1F_EE_Electron_et_total_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]], ev.wPU)
                         plot['TH1F_EE_Electron_et_total_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]], ev.wPU)
                         if ev.EE_relisoA[0] < 0.1:
                             plot['TH1F_EE_Electron_et_passedIso_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxA[0]], ev.wPU)
                         if ev.EE_relisoB[0] < 0.1:
                             plot['TH1F_EE_Electron_et_passedIso_MC'].Fill(ev.ElectronCandidate_et[ev.EE_idxB[0]], ev.wPU)
                         """

                     num += 1
                     if num >200000:
                          num = 0
                          break


    ##################################################################################################

    ## Build efficiencies
    plot['Efficiency_DMDM_DGM_pt_MC']     = r.TEfficiency(plot['TH1F_DMDM_DGM_pt_passedIso_MC'], plot['TH1F_DMDM_DGM_pt_total_MC'])
    plot['Efficiency_DMDM_DGM_pt_MC'].SetTitle(';Reconstructed muon p_{T} (GeV); Isolation cut efficiency')
    plot['Efficiency_DMDM_DGM_pt_DATA']     = r.TEfficiency(plot['TH1F_DMDM_DGM_pt_passedIso_DATA'], plot['TH1F_DMDM_DGM_pt_total_DATA'])
    plot['Efficiency_DMDM_DGM_pt_DATA'].SetTitle(';Reconstructed muon p_{T} (GeV); Isolation cut efficiency')
    plot['Efficiency_EE_Electron_et_MC']     = r.TEfficiency(plot['TH1F_EE_Electron_et_passedIso_MC'], plot['TH1F_EE_Electron_et_total_MC'])
    plot['Efficiency_EE_Electron_et_MC'].SetTitle(';Reconstructed electron E_{T} (GeV); Isolation cut efficiency')
    plot['Efficiency_EE_Electron_et_DATA']     = r.TEfficiency(plot['TH1F_EE_Electron_et_passedIso_DATA'], plot['TH1F_EE_Electron_et_total_DATA'])
    plot['Efficiency_EE_Electron_et_DATA'].SetTitle(';Reconstructed electron E_{T} (GeV); Isolation cut efficiency')



    ## Save
    outputFile = TFile(EOSPATH + 'Isolation-SFs_'+era+'/TH1F_Isolation_'+era+'.root', 'RECREATE')
    for key in plot.keys():
        print('>> ' + key + ' saved!')
        plot[key].Write()

    canvas = Canvas.Canvas("MuonIsolation_"+era+"_Eff_pt_1D", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_DMDM_DGM_pt_DATA'])
    hMC_ = getHistoFromEff(plot['Efficiency_DMDM_DGM_pt_MC'])
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hMC_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.saveRatio(1, 1, 0, lumi[era], hdata = hdata_, hMC = hMC_, r_ymin = 0.96, r_ymax = 1.04, label = 'Scale factor',outputDir = EOSPATH + 'Isolation-SFs_'+era, inProgress = True)

    SF_ = hdata_.Clone(hdata_.GetName() + '_SF')
    SF_.Divide(hMC_)
    SF_.Write()

    canvas = Canvas.Canvas("ElectronIsolation_"+era+"_Eff_pt_1D", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_EE_Electron_et_DATA'])
    hMC_ = getHistoFromEff(plot['Efficiency_EE_Electron_et_MC'])
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hMC_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.saveRatio(1, 1, 0, lumi[era], hdata = hdata_, hMC = hMC_, r_ymin = 0.92, r_ymax = 1.08, label = 'Scale factor',outputDir = EOSPATH + 'Isolation-SFs_'+era, inProgress = True)

    SF_ = hdata_.Clone(hdata_.GetName() + '_SF')
    SF_.Divide(hMC_)
    SF_.Write()


    outputFile.Close()


