
import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
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

import include.Canvas as Canvas
import include.Sample as Sample
import include.helper as helper
import include.CutManager as CutManager

#print(WORKPATH, WORKPATH)
#print(GALAPAGOPATH, GALAPAGOPATH)


if __name__ == "__main__":


    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ###########################
    ####   Parser object   ####
    ###########################
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    #parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    parser.add_option('-e', '--era', action='store', type=str, dest='era', default='', help='2018 era')
    #parser.add_option('-n', '--nmax', action='store', type=int, dest='nmax', default=0, help='Path to file')
    (opts, args) = parser.parse_args()

    era = opts.era

    ##################################
    ####   Variable declaration   ####
    ##################################

    dr_logbin = np.logspace(-3, 1, 101)
    iso_logbin = np.logspace(-2, 2, 101)

    #### -----------------
    #### ---- Histograms
    #### -----------------

    bin_iso_log = np.logspace(-3, 2)
    histograms = {}

    #########################
    ####   Load sample   ####
    #########################
    Signals = []
    Signals.append('HSS_125_50_1')
    Signals.append('HSS_125_50_10')
    Signals.append('HSS_125_50_100')
    Signals.append('HSS_125_50_1000')
    Signals.append('HSS_300_50_1')
    Signals.append('HSS_300_50_10')
    Signals.append('HSS_300_50_100')
    Signals.append('HSS_300_50_1000')
    Signals.append('HSS_500_50_1')
    Signals.append('HSS_500_50_10')
    Signals.append('HSS_500_50_100')
    Signals.append('HSS_500_50_1000')
    Signals.append('HSS_800_50_1')
    Signals.append('HSS_800_50_10')
    Signals.append('HSS_800_50_100')
    Signals.append('HSS_800_50_1000')
    Signals.append('HSS_1000_350_1')
    Signals.append('HSS_1000_350_10')
    Signals.append('HSS_1000_350_100')
    Signals.append('HSS_1000_350_1000')

    Signals = [x + '_' + era for x in Signals]


    if era=='2016' or era=='2016APV':
        treeSI = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + '/dat/CombSignal_2016UL_Fall22.dat', Signals, 'SI'), name = 'SI', isdata = 0)
    elif era=='2017':
        treeSI = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + '/dat/CombSignal_2017UL_Fall22.dat', Signals, 'SI'), name = 'SI', isdata = 0)
    elif era=='2018':
        treeSI = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + '/dat/CombSignal_2018UL_Fall22.dat', Signals, 'SI'), name = 'SI', isdata = 0)


    file_sf_iso = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/Isolation_ScaleFactors_Spring23.root")
    if '2016APV' == era:
        sf_mm_iso = file_sf_iso.Get("SF_Iso_DMDM_2016APV")
        sf_ee_iso = file_sf_iso.Get("SF_Iso_EE_2016APV")
    elif '2016' == era:
        sf_mm_iso = file_sf_iso.Get("SF_Iso_DMDM_2016")
        sf_ee_iso = file_sf_iso.Get("SF_Iso_EE_2016")
    elif '2017' == era:
        sf_ee_iso = file_sf_iso.Get("SF_Iso_EE_2017")
    elif '2018' == era:
        sf_mm_iso = file_sf_iso.Get("SF_Iso_DMDM_2018")
        sf_ee_iso = file_sf_iso.Get("SF_Iso_EE_2018")




    ###################################
    ####   Loop over tree events   ####
    ###################################
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


    for b in treeSI.blocks:
        for s in b.samples:
            histograms[s.name + '_MM_corrected']   = r.TH1F(s.name + '_MM_corrected', "", 1, 0, 1)
            histograms[s.name + '_MM_uncorrected'] = r.TH1F(s.name + '_MM_uncorrected', "", 1, 0, 1)
            histograms[s.name + '_EE_corrected']   = r.TH1F(s.name + '_EE_corrected', "", 1, 0, 1)
            histograms[s.name + '_EE_uncorrected'] = r.TH1F(s.name + '_EE_uncorrected', "", 1, 0, 1)
            for t in s.ttrees:

                print('New tree with:', t.GetEntries())

                for e,ev in enumerate(t):

                    for j in range(0, ev.nEE):
                        iee = j
                        if not eval(ee_selection):
                            continue
                        if not ev.EE_relisoA[iee] < 0.1 and ev.EE_relisoB[iee] < 0.1: continue

                        bxA = sf_ee_iso.GetXaxis().FindBin(abs(ev.ElectronCandidate_et[ev.EE_idxA[iee]]))
                        bxB = sf_ee_iso.GetXaxis().FindBin(abs(ev.ElectronCandidate_et[ev.EE_idxB[iee]]))
                        if bxA < 7: bxA = 7
                        if bxB < 7: bxB = 7
                        sf_A = sf_ee_iso.GetBinContent(bxA) 
                        sf_B = sf_ee_iso.GetBinContent(bxB) 
                        sf_ee = sf_A*sf_B

                        histograms[s.name + '_EE_corrected'].Fill(0, sf_ee)
                        histograms[s.name + '_EE_uncorrected'].Fill(0)


                    for j in range(0, ev.nDMDM):
                        imm = j
                        if not eval(mumu_selection):
                            continue
                        if not ev.DGM_relPFiso[ev.DMDM_idxA[imm]] < 0.2 and ev.DGM_relPFiso[ev.DMDM_idxB[imm]] < 0.2: continue

                        sf_A = sf_mm_iso.GetBinContent(sf_mm_iso.GetXaxis().FindBin(abs(ev.DGM_pt[ev.DMDM_idxA[imm]]))) 
                        sf_B = sf_mm_iso.GetBinContent(sf_mm_iso.GetXaxis().FindBin(abs(ev.DGM_pt[ev.DMDM_idxB[imm]]))) 
                        sf_mm = sf_A*sf_B

                        histograms[s.name + '_MM_corrected'].Fill(0, sf_mm)
                        histograms[s.name + '_MM_uncorrected'].Fill(0)


    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f_systematics'+era+'.root', 'RECREATE')




    #### Plot final plot
    labels = []
    labels.append('HSS_125_50')
    labels.append('HSS_300_50')
    labels.append('HSS_500_50')
    labels.append('HSS_800_50')
    labels.append('HSS_1000_350')
    EE_counter = {}
    MM_counter = {}
    for label in labels:
        EE_counter[label] = [0, 0] # corrected : uncorrected
        MM_counter[label] = [0, 0] # corrected : uncorrected

    for key in histograms.keys():
        if not 'EE' in key: continue
        for key2 in EE_counter:
            if key2 in key:
                if 'uncorrected' in key:
                    EE_counter[key2][1] = EE_counter[key2][1] + histograms[key].GetBinContent(1)
                else:
                    EE_counter[key2][0] = EE_counter[key2][0] + histograms[key].GetBinContent(1)

    for key in histograms.keys():
        if not 'MM' in key: continue
        for key2 in MM_counter:
            if key2 in key:
                if 'uncorrected' in key:
                    MM_counter[key2][1] = MM_counter[key2][1] + histograms[key].GetBinContent(1)
                else:
                    MM_counter[key2][0] = MM_counter[key2][0] + histograms[key].GetBinContent(1)

    EE_correctedhisto = r.TH1F('EE_correctedhisto', '', len(labels), 0, len(labels))
    EE_uncorrectedhisto = r.TH1F('EE_uncorrectedhisto', '', len(labels), 0, len(labels))
    EE_diffhisto = r.TH1F('EE_diffhisto', '', len(labels), 0, len(labels))
    MM_correctedhisto = r.TH1F('MM_correctedhisto', '', len(labels), 0, len(labels))
    MM_uncorrectedhisto = r.TH1F('MM_uncorrectedhisto', '', len(labels), 0, len(labels))
    MM_diffhisto = r.TH1F('MM_diffhisto', '', len(labels), 0, len(labels))


    for l,label in enumerate(labels):
         EE_correctedhisto.SetBinContent(l+1, EE_counter[label][0])
         EE_uncorrectedhisto.SetBinContent(l+1, EE_counter[label][1])
         EE_diffhisto.SetBinContent(l+1, EE_counter[label][0]/EE_counter[label][1])
         if era != '2017':
             MM_correctedhisto.SetBinContent(l+1, MM_counter[label][0])
             MM_uncorrectedhisto.SetBinContent(l+1, MM_counter[label][1])
             MM_diffhisto.SetBinContent(l+1, MM_counter[label][0]/MM_counter[label][1])

    EE_correctedhisto.Write()
    EE_uncorrectedhisto.Write()
    EE_diffhisto.Write()
    MM_correctedhisto.Write()
    MM_uncorrectedhisto.Write()
    MM_diffhisto.Write()
    


    #### Write everything to use later:
    for key in histograms.keys():
        histograms[key].Write()
    outputFile.Close()




