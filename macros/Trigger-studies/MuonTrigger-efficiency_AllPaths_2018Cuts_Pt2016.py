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

EOSPATH = '/eos/user/f/fernance/www/DisplacedLeptons-analysis/Trigger-studies/Spring23/'

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
        passed = ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_MET200 or ev.HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight or ev.HLT_PFMET170_HBHECleaned or ev.HLT_PFMET300 or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        #passed = ev.HLT_PFMET120_PFMHT90_IDTight or ev.HLT_PFMET120_PFMHT100_IDTight or ev.HLT_PFMET120_PFMHT110_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight 
    elif year == '2017':
        passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
    elif year == '2018':
        passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ev.HLT_CaloMET350_HBHECleaned or ev.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or ev.HLT_PFMET250_HBHECleaned or ev.HLT_PFMET200_HBHE_BeamHaloCleaned or ev.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        #passed = ev.HLT_PFMET120_PFMHT120_IDTight or ev.HLT_PFMET120_PFMHT120_IDTight_PFHT60

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

    Datasets_2016APV = []
    Datasets_2016APV.append('MET_Muon_Run2016B_HIPM')
    Datasets_2016APV.append('MET_Muon_Run2016C_HIPM')
    Datasets_2016APV.append('MET_Muon_Run2016D_HIPM')
    Datasets_2016APV.append('MET_Muon_Run2016E_HIPM')
    Datasets_2016APV.append('MET_Muon_Run2016F_HIPM')

    Datasets_2016 = []
    Datasets_2016.append('MET_Muon_Run2016F_noHIPM')
    Datasets_2016.append('MET_Muon_Run2016G_noHIPM')
    Datasets_2016.append('MET_Muon_Run2016H_noHIPM')

    Datasets_2018 = []
    Datasets_2018.append('MET_Muon_Run2018A')
    Datasets_2018.append('MET_Muon_Run2018B')
    Datasets_2018.append('MET_Muon_Run2018C')
    Datasets_2018.append('MET_Muon_Run2018D')



    MC_2016APV = ['TTTo2L2Nu_preVFP', 'DYJetsToLL_M-50_preVFP']
    MC_2016 = ['TTTo2L2Nu_postVFP', 'DYJetsToLL_M-50_postVFP']
    MC_2018 = ['TTTo2L2Nu_2018', 'DYJetsToLL_M-50_2018']

    if era == '2016APV':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2016APV, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', MC_2016APV, 'DATA'), name = year, isdata = 1 )
    elif era == '2016':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2016, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', MC_2016 , 'DATA'), name = year, isdata = 1 )
    elif era == '2018':
        treeDATA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_MET.dat', Datasets_2018, 'DATA'), name = year, isdata = 1 )
        treeMC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy_Spring23.dat', MC_2018, 'DATA'), name = year, isdata = 1 )


    #########################
    ####   Init plots    ####
    #########################

    pt1_bin = np.array([25., 60., 100., 200.])
    pt2_bin = np.array([25., 60., 100., 200.])
    eta2_bin = np.array([-2.0, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.0])
    eta2_bin = np.array([-2.0, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.0])
    mass_bin = np.array([15., 60., 100., 200.])

    plot = {}
    plot['PassMET_pt2_pt1_DATA'] = copy.deepcopy(r.TH2F('PassMET_pt2_pt1_DATA', ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['PassTRG_pt2_pt1_DATA'] = copy.deepcopy(r.TH2F('PassTRG_pt2_pt1_DATA', ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['Efficiency_HLT_Full_pt_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_2d_DATA', ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
    plot['Efficiency_HLT_Full_eta_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_eta_2d_DATA', ";Subleading muon #eta;Leading muon #eta", len(eta2_bin)-1, eta2_bin, len(eta2_bin)-1, eta2_bin))
    plot['Efficiency_HLT_Full_pt_eta_2d_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_eta_2d_DATA', ";Subleading muon p_{T};Subeading muon #eta", len(pt2_bin)-1, pt2_bin, len(eta2_bin)-1, eta2_bin))
    plot['Efficiency_HLT_Full_pt_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_DATA', ";Subleading muon p_{T} (GeV) ;Efficiency", len(pt2_bin)-1, pt2_bin))
    plot['Efficiency_HLT_Full_mass_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_mass_DATA', ";Dimuon mass m_{#mu#mu} (GeV) ;Efficiency", len(mass_bin)-1, mass_bin))
    plot['Efficiency_HLT_Full_pt_MET40_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MET40_DATA', ";Subleading muon p_{T}; Efficiency", len(pt2_bin)-1, pt2_bin))
    plot['Efficiency_HLT_Full_pt_MET60_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MET60_DATA', ";Subleading muon p_{T};Efficiency", len(pt2_bin)-1, pt2_bin))
    plot['Efficiency_HLT_Full_pt_MET80_DATA'] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MET80_DATA', ";Subleading muon p_{T};Efficiency", len(pt2_bin)-1, pt2_bin))

    for key in ['TTTo2L2Nu', 'DYJetsToLL_M-50']:
        plot['PassMET_pt2_pt1_'+key] = copy.deepcopy(r.TH2F('PassMET_pt2_pt1_'+key, ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
        plot['PassTRG_pt2_pt1_'+key] = copy.deepcopy(r.TH2F('PassTRG_pt2_pt1_'+key, ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
        plot['Efficiency_HLT_Full_pt_2d_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_2d_'+key, ";Subleading muon p_{T};Leading muon p_{T}", len(pt2_bin)-1, pt2_bin, len(pt1_bin)-1, pt1_bin))
        plot['Efficiency_HLT_Full_eta_2d_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_eta_2d_+key', ";Subleading muon #eta;Leading muon #eta", len(eta2_bin)-1, eta2_bin, len(eta2_bin)-1, eta2_bin))
        plot['Efficiency_HLT_Full_pt_eta_2d_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_eta_2d_'+key, ";Subleading muon p_{T};Subleading muon #eta", len(pt2_bin)-1, pt2_bin, len(eta2_bin)-1, eta2_bin))
        plot['Efficiency_HLT_Full_pt_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_'+key, ";Subleading muon p_{T} (GeV) ;Efficiency", len(pt2_bin)-1, pt2_bin))
        plot['Efficiency_HLT_Full_mass_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_mass_'+key, ";Dimuon mass m_{#mu#mu} (GeV) ;Efficiency", len(mass_bin)-1, mass_bin))
        plot['Efficiency_HLT_Full_pt_MET40_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MET40_'+key, ";Subleading muon p_{T};Efficiency", len(pt2_bin)-1, pt2_bin))
        plot['Efficiency_HLT_Full_pt_MET60_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MET60_'+key, ";Subleading muon p_{T};Efficiency", len(pt2_bin)-1, pt2_bin))
        plot['Efficiency_HLT_Full_pt_MET80_'+key] = copy.deepcopy(r.TEfficiency('Efficiency_HLT_Full_pt_MET80_'+key, ";Subleading muon p_{T};Efficiency", len(pt2_bin)-1, pt2_bin))
        





    for p in plot.keys():
        r.SetOwnership(plot[p], 0)


    #########################
    ####   Cut values    ####
    #########################
    #minPt = 30 if year == '2016' else 30
    #minCosAlpha = -0.8 if year == '2016' else -0.8
    minPt = 30
    minCosAlpha = -0.9 


    ###############################
    ####   Loop over events    ####
    ###############################

    ### Data:
    for b in treeDATA.blocks:
        for s in b.samples: 
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     if not (ev.nDMDM > 0):
                         continue

                     if not passedMETTrigger(ev, year):
                         continue

                     pt_values = []
                     eta_values = []
                     indexes = []
                     for i in range(0, ev.nDGM):
                         if abs(ev.DGM_eta[i]) > 2.:
                             continue
                         if ev.DGM_relPFiso[i] > 0.2:
                             continue
                         pt_values.append(ev.DGM_pt[i])
                         eta_values.append(ev.DGM_eta[i])
                         indexes.append(i)

                     if len(pt_values) < 2:
                         continue

                     pt_eta = zip(pt_values, eta_values)
                     pt_idx = zip(pt_values, indexes)
                     pt_eta.sort(reverse = True)
                     pt_idx.sort(reverse = True)
                     pt_ord = [x[0] for x in pt_eta]
                     eta_ord = [x[1] for x in pt_eta]
                     idx_ord = [x[1] for x in pt_idx]


                     mu1 = r.TLorentzVector()
                     mu2 = r.TLorentzVector()
                     mu1.SetPtEtaPhiM(ev.DGM_pt[idx_ord[0]], ev.DGM_eta[idx_ord[0]], ev.DGM_phi[idx_ord[0]], 105e-3)
                     mu2.SetPtEtaPhiM(ev.DGM_pt[idx_ord[1]], ev.DGM_eta[idx_ord[1]], ev.DGM_phi[idx_ord[1]], 105e-3)

                     if pt_values[0] < minPt or pt_values[1] < minPt:
                         continue
                     if (mu1+mu2).M() < 15.:
                         continue
                     if math.cos(mu1.Angle(mu2.Vect())) < minCosAlpha:
                         continue
                     if mu1.DeltaR(mu2) < 0.1:
                         continue


                     plot['Efficiency_HLT_Full_pt_2d_DATA'].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])
                     plot['Efficiency_HLT_Full_eta_2d_DATA'].Fill(passedMuonTrigger(ev, year), eta_ord[1], eta_ord[0])
                     plot['Efficiency_HLT_Full_pt_eta_2d_DATA'].Fill(passedMuonTrigger(ev, year), pt_ord[1], eta_ord[1])
                     plot['Efficiency_HLT_Full_pt_DATA'].Fill(passedMuonTrigger(ev, year), pt_ord[1])
                     plot['Efficiency_HLT_Full_mass_DATA'].Fill(passedMuonTrigger(ev, year), (mu1+mu2).M())

                     plot['PassMET_pt2_pt1_DATA'].Fill(pt_ord[1], pt_ord[0])
                     if passedMuonTrigger(ev, year):
                         plot['PassTRG_pt2_pt1_DATA'].Fill(pt_ord[1], pt_ord[0])
                     
                     if ev.MET_pt > 40:
                         plot['Efficiency_HLT_Full_pt_MET40_DATA'].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])
                     if ev.MET_pt > 60:
                         plot['Efficiency_HLT_Full_pt_MET60_DATA'].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])
                     if ev.MET_pt > 80:
                         plot['Efficiency_HLT_Full_pt_MET80_DATA'].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])

                     num += 1
                     if num > 10:
                         num = 0
                         #break

    ### MC:
    weights = {}
    for b in treeMC.blocks:
        for s in b.samples: 
            if era == '2016APV': 
                key = s.name.replace('_preVFP', '')
            elif era == '2016': 
                if '_postVFP' in s.name: key = s.name.replace('_postVFP', '')
                if '_preVFP' in s.name: key = s.name.replace('_preVFP', '')
            elif era == '2018': 
                key = s.name.replace('_2018', '')
            weights[key] = s.lumWeight
            for t in s.ttrees:
                num = 0
                for e,ev in enumerate(t):

                     if not (ev.nDMDM > 0):
                         continue

                     if not passedMETTrigger(ev, year):
                         continue

                     pt_values = []
                     eta_values = []
                     indexes = []
                     for i in range(0, ev.nDGM):
                         if abs(ev.DGM_eta[i]) > 2.:
                             continue
                         if ev.DGM_relPFiso[i] > 0.2:
                             continue
                         pt_values.append(ev.DGM_pt[i])
                         eta_values.append(ev.DGM_eta[i])
                         indexes.append(i)

                     if len(pt_values) < 2:
                         continue

                     pt_eta = zip(pt_values, eta_values)
                     pt_idx = zip(pt_values, indexes)
                     pt_eta.sort(reverse = True)
                     pt_idx.sort(reverse = True)
                     pt_ord = [x[0] for x in pt_eta]
                     eta_ord = [x[1] for x in pt_eta]
                     idx_ord = [x[1] for x in pt_idx]


                     mu1 = r.TLorentzVector()
                     mu2 = r.TLorentzVector()
                     mu1.SetPtEtaPhiM(ev.DGM_pt[idx_ord[0]], ev.DGM_eta[idx_ord[0]], ev.DGM_phi[idx_ord[0]], 105e-3)
                     mu2.SetPtEtaPhiM(ev.DGM_pt[idx_ord[1]], ev.DGM_eta[idx_ord[1]], ev.DGM_phi[idx_ord[1]], 105e-3)

                     if pt_values[0] < minPt or pt_values[1] < minPt:
                         continue
                     if (mu1+mu2).M() < 15.:
                         continue
                     if math.cos(mu1.Angle(mu2.Vect())) < minCosAlpha:
                         continue
                     if mu1.DeltaR(mu2) < 0.1:
                         continue

                     num += 1
                     if num > 10:
                          num = 0
                          #break

                     plot['Efficiency_HLT_Full_pt_2d_' + key].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])
                     plot['Efficiency_HLT_Full_eta_2d_' + key].Fill(passedMuonTrigger(ev, year), eta_ord[1], eta_ord[0])
                     plot['Efficiency_HLT_Full_pt_eta_2d_' + key].Fill(passedMuonTrigger(ev, year), pt_ord[1], eta_ord[1])
                     plot['Efficiency_HLT_Full_pt_' + key].Fill(passedMuonTrigger(ev, year), pt_ord[1])
                     plot['Efficiency_HLT_Full_mass_' + key].Fill(passedMuonTrigger(ev, year), (mu1+mu2).M())

                     plot['PassMET_pt2_pt1_' + key].Fill(pt_ord[1], pt_ord[0])
                     if passedMuonTrigger(ev, year):
                         plot['PassTRG_pt2_pt1_' + key].Fill(pt_ord[1], pt_ord[0])

	             if ev.MET_pt > 40:
                         plot['Efficiency_HLT_Full_pt_MET40_' + key].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])
                     if ev.MET_pt > 60:
                         plot['Efficiency_HLT_Full_pt_MET60_' + key].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])
                     if ev.MET_pt > 80:
                         plot['Efficiency_HLT_Full_pt_MET80_' + key].Fill(passedMuonTrigger(ev, year), pt_ord[1], pt_ord[0])


    ##################################################################################################
    ## Weight the histograms 
    """
    for key in plot.keys():
        if 'DATA' in key:
            continue
        if 'DYJetsToLL_M-50' in key:
            plot[key].SetWeight(weights['DYJetsToLL_M-50'])
        if 'TTTo2L2Nu' in key:
            plot[key].SetWeight(weights['TTTo2L2Nu'])
    """       
 

    ##################################################################################################
    ## Plot

    outputFile = TFile(EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/TH1F_muontrigger_'+era+'.root', 'RECREATE')
    for key in plot.keys():
        plot[key].Write()

    plot['Efficiency_HLT_Full_pt_2d_MC']     = plot['Efficiency_HLT_Full_pt_2d_TTTo2L2Nu'] 
    plot['Efficiency_HLT_Full_eta_2d_MC']    = plot['Efficiency_HLT_Full_eta_2d_TTTo2L2Nu']
    plot['Efficiency_HLT_Full_pt_eta_2d_MC'] = plot['Efficiency_HLT_Full_pt_eta_2d_TTTo2L2Nu'] 
    plot['Efficiency_HLT_Full_pt_MC']        = combineEfficiency([plot['Efficiency_HLT_Full_pt_DYJetsToLL_M-50'], plot['Efficiency_HLT_Full_eta_2d_TTTo2L2Nu']], [weights['DYJetsToLL_M-50'], weights['TTTo2L2Nu']])
    plot['Efficiency_HLT_Full_pt_MET40_MC'] = combineEfficiency([plot['Efficiency_HLT_Full_pt_MET40_DYJetsToLL_M-50'], plot['Efficiency_HLT_Full_pt_MET40_TTTo2L2Nu']], [weights['DYJetsToLL_M-50'], weights['TTTo2L2Nu']])
    plot['Efficiency_HLT_Full_pt_MET60_MC'] = combineEfficiency([plot['Efficiency_HLT_Full_pt_MET60_DYJetsToLL_M-50'], plot['Efficiency_HLT_Full_pt_MET60_TTTo2L2Nu']], [weights['DYJetsToLL_M-50'], weights['TTTo2L2Nu']])
    plot['Efficiency_HLT_Full_pt_MET80_MC'] = combineEfficiency([plot['Efficiency_HLT_Full_pt_MET80_DYJetsToLL_M-50'], plot['Efficiency_HLT_Full_pt_MET80_TTTo2L2Nu']], [weights['DYJetsToLL_M-50'], weights['TTTo2L2Nu']])

    ### Scale factor plot
    plot['Efficiency_pt2_pt1_DATA'] = plot['PassTRG_pt2_pt1_DATA'].Clone('Efficiency_pt2_pt1_DATA')
    plot['Efficiency_pt2_pt1_DATA'].Divide(plot['PassMET_pt2_pt1_DATA'])
    plot['Efficiency_pt2_pt1_TTTo2L2Nu'] = plot['PassTRG_pt2_pt1_TTTo2L2Nu'].Clone('Efficiency_pt2_pt1_TTTo2L2Nu')
    plot['Efficiency_pt2_pt1_TTTo2L2Nu'].Divide(plot['PassMET_pt2_pt1_TTTo2L2Nu'])
    plot['ScaleFactor_pt2_pt1'] = plot['Efficiency_pt2_pt1_DATA'].Clone('ScaleFactor_pt2_pt') 
    plot['ScaleFactor_pt2_pt1'].Divide(plot['Efficiency_pt2_pt1_TTTo2L2Nu'])

    plot['Efficiency_pt2_pt1_DATA'].Write()
    plot['Efficiency_pt2_pt1_TTTo2L2Nu'].Write()
    plot['ScaleFactor_pt2_pt1'].Write()


    ### Efficiency (1d)

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Eff_full_pt_1D", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hMC_ = plot['Efficiency_HLT_Full_pt_MC']
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hMC_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.addLatex(0.9, 0.88, era, size = 0.035, align = 31)
    canvas.saveRatio(1, 1, 0, '', hdata = hdata_, hMC = hMC_, r_ymin = 0.7, r_ymax = 1.0, label = 'Scale factor',outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Eff_full_pt_1D_DYJetsToLL_M-50", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hDYJetsToLL_M50_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DYJetsToLL_M-50'])
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hDYJetsToLL_M50_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.addLatex(0.9, 0.88, era, size = 0.035, align = 31)
    canvas.saveRatio(1, 1, 0, '', hdata = hdata_, hMC = hDYJetsToLL_M50_, r_ymin = 0.7, r_ymax = 1.0, label = 'Scale factor',outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Eff_full_pt_1D_TTTo2L2Nu", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hTTTo2L2Nu_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_TTTo2L2Nu'])
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hTTTo2L2Nu_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.addLatex(0.9, 0.88, era, size = 0.035, align = 31)
    canvas.saveRatio(1, 1, 0, '', hdata = hdata_, hMC = hTTTo2L2Nu_, r_ymin = 0.7, r_ymax = 1.0, label = 'Scale factor',outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Eff_full_mass_1D_TTTo2L2Nu", 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    hdata_ = getHistoFromEff(plot['Efficiency_HLT_Full_mass_DATA'])
    hTTTo2L2Nu_ = getHistoFromEff(plot['Efficiency_HLT_Full_mass_TTTo2L2Nu'])
    canvas.addHisto(hdata_,'P', 'Data', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hTTTo2L2Nu_,'P,SAME', 'Simulation', 'pl', r.kBlue, True, 0, marker = 25)
    canvas.addLatex(0.9, 0.88, era, size = 0.035, align = 31)
    canvas.saveRatio(1, 1, 0, '', hdata = hdata_, hMC = hTTTo2L2Nu_, r_ymin = 0.7, r_ymax = 1.0, label = 'Scale factor',outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)

    ### Sys variations

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SFvar_full_pt_1D_MC", 'png,pdf', 0.46, 0.72, 0.86, 0.89, 1)
    hSF_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hMC_ = plot['Efficiency_HLT_Full_pt_MC']
    hSF_MET40 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET40_DATA'])
    hMC_MET40 = plot['Efficiency_HLT_Full_pt_MET40_MC']
    hSF_MET60 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET60_DATA'])
    hMC_MET60 = plot['Efficiency_HLT_Full_pt_MET60_MC']
    hSF_MET80 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET80_DATA'])
    hMC_MET80 = plot['Efficiency_HLT_Full_pt_MET80_MC']
    hSF_.Divide(hMC_)
    hSF_MET40.Divide(hMC_MET40)
    hSF_MET60.Divide(hMC_MET60)
    hSF_MET80.Divide(hMC_MET80)
    hSF_.GetYaxis().SetTitle('Scale factor')
    hsys_ = createSysPlot(hSF_, 0.03)
    canvas.addHisto(hsys_,'E2', '3% syst.', 'f', '', True, 1)
    canvas.addHisto(hSF_,'P,SAME', 'Scale factor', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hSF_MET80,'P,SAME', 'MET > 80 GeV', 'pl', r.kBlue, True, 2, marker = 26)
    canvas.addHisto(hSF_MET60,'P,SAME', 'MET > 60 GeV', 'pl', r.kBlue, True, 3, marker = 25)
    canvas.addHisto(hSF_MET40,'P,SAME', 'MET > 40 GeV', 'pl', r.kBlue, True, 4, marker = 32)
    canvas.addLatex(0.9, 0.93, era, size = 0.035, align = 31)
    canvas.save(1, 1, 0, '', '', ymin=0.65, ymax=1.2, outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SFvar_full_pt_1D_DYJetsToLL_M-50", 'png,pdf', 0.46, 0.72, 0.86, 0.89, 1)
    hSF_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hMC_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DYJetsToLL_M-50'])
    hSF_MET40 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET40_DATA'])
    hMC_MET40 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET40_DYJetsToLL_M-50'])
    hSF_MET60 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET60_DATA'])
    hMC_MET60 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET60_DYJetsToLL_M-50'])
    hSF_MET80 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET80_DATA'])
    hMC_MET80 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET80_DYJetsToLL_M-50'])
    hSF_.Divide(hMC_)
    hSF_MET40.Divide(hMC_MET40)
    hSF_MET60.Divide(hMC_MET60)
    hSF_MET80.Divide(hMC_MET80)
    hSF_.GetYaxis().SetTitle('Scale factor')
    hsys_ = createSysPlot(hSF_, 0.03)
    canvas.addHisto(hsys_,'E2', '3% syst.', 'f', '', True, 1)
    canvas.addHisto(hSF_,'P,SAME', 'Scale factor', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hSF_MET80,'P,SAME', 'MET > 80 GeV', 'pl', r.kBlue, True, 2, marker = 26)
    canvas.addHisto(hSF_MET60,'P,SAME', 'MET > 60 GeV', 'pl', r.kBlue, True, 3, marker = 25)
    canvas.addHisto(hSF_MET40,'P,SAME', 'MET > 40 GeV', 'pl', r.kBlue, True, 4, marker = 32)
    canvas.addLatex(0.9, 0.93, era, size = 0.035, align = 31)
    canvas.save(1, 1, 0, '', '', ymin=0.65, ymax=1.2, outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SFvar_full_pt_1D_TTTo2L2Nu", 'png,pdf', 0.46, 0.72, 0.86, 0.89, 1)
    hSF_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_DATA'])
    hMC_ = getHistoFromEff(plot['Efficiency_HLT_Full_pt_TTTo2L2Nu'])
    hSF_MET40 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET40_DATA'])
    hMC_MET40 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET40_TTTo2L2Nu'])
    hSF_MET60 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET60_DATA'])
    hMC_MET60 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET60_TTTo2L2Nu'])
    hSF_MET80 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET80_DATA'])
    hMC_MET80 = getHistoFromEff(plot['Efficiency_HLT_Full_pt_MET80_TTTo2L2Nu'])
    hSF_.Divide(hMC_)
    hSF_MET40.Divide(hMC_MET40)
    hSF_MET60.Divide(hMC_MET60)
    hSF_MET80.Divide(hMC_MET80)
    hSF_.GetYaxis().SetTitle('Scale factor')
    hsys_ = createSysPlot(hSF_, 0.03)
    canvas.addHisto(hsys_,'E2', '3% syst.', 'f', '', True, 1)
    canvas.addHisto(hSF_,'P,SAME', 'Scale factor', 'pl', r.kBlack, True, 0, marker = 20)
    canvas.addHisto(hSF_MET80,'P,SAME', 'MET > 80 GeV', 'pl', r.kBlue, True, 2, marker = 26)
    canvas.addHisto(hSF_MET60,'P,SAME', 'MET > 60 GeV', 'pl', r.kBlue, True, 3, marker = 25)
    canvas.addHisto(hSF_MET40,'P,SAME', 'MET > 40 GeV', 'pl', r.kBlue, True, 4, marker = 32)
    canvas.addLatex(0.9, 0.93, era, size = 0.035, align = 31)
    canvas.save(1, 1, 0, '', '', ymin=0.65, ymax=1.2, outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False)


    ### Efficiency (2d)

    r.gStyle.SetPadRightMargin(0.19)

    ### pt dependence

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Data_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_2d_DATA'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_MC_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_2d_MC'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_pt_2d_DATA'], plot['Efficiency_HLT_Full_pt_2d_MC'])
    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SF_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SFErr_full_pt_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    SF.Write()
    SFErr.Write()


    ### eta dependence

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Data_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_eta_2d_DATA'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_MC_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_eta_2d_MC'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_eta_2d_DATA'], plot['Efficiency_HLT_Full_eta_2d_MC'])
    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SF_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SFErr_full_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    SF.Write()
    SFErr.Write()


    ### Full subleading muon dependence

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_Data_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_eta_2d_DATA'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Efficiency')

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_MC_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.add2DRate(plot['Efficiency_HLT_Full_pt_eta_2d_MC'],'COLZ,TEXT', 0.0, 1.0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Efficiency')

    SF, SFErr = getSFPlot(plot['Efficiency_HLT_Full_pt_eta_2d_DATA'], plot['Efficiency_HLT_Full_pt_eta_2d_MC'])
    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SF_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SF,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Scale factor')

    canvas = Canvas.Canvas("MuonTrigger_"+era+"_SFErr_full_pt_eta_2D", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(SFErr,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.8, 0.93, era, size = 0.035, align = 31)
    canvas.save(0, 1, 0, '', '', outputDir = EOSPATH + 'MuonTrigger-SFs_AllPaths_2018Cuts_Pt2016/', inProgress = False, is2d = True, labelz = 'Scale factor uncertainty (stat)')

    SF.Write()
    SFErr.Write()

    r.gStyle.SetPadRightMargin(0.1)

    outputFile.Close()



