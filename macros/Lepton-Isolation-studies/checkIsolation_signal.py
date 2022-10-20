
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
    Signals2016 = []
    Signals2016.append('HSS_125_30_1_2016')
    Signals2016.append('HSS_125_30_10_2016')
    Signals2016.append('HSS_125_30_100_2016')
    Signals2016.append('HSS_125_30_1000_2016')
    #Signals2016.append('HSS_125_50_10000_2016')
    Signals2016.append('HSS_400_50_1_2016')
    Signals2016.append('HSS_400_50_10_2016')
    Signals2016.append('HSS_400_50_100_2016')
    Signals2016.append('HSS_400_50_1000_2016')
    #Signals2016.append('HSS_400_50_10000_2016')
    Signals2016.append('HSS_1000_150_1_2016')
    Signals2016.append('HSS_1000_150_10_2016')
    Signals2016.append('HSS_1000_150_100_2016')
    Signals2016.append('HSS_1000_150_1000_2016')
    #Signals2016.append('HSS_1000_350_10000_2016')

    treeSI = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_2016UL_Summer22.dat', Signals2016, 'MC'), name = '2016', isdata = 0 )


    ###################################
    ####   Loop over tree events   ####
    ###################################
    cm = CutManager.CutManager()

    for b in treeSI.blocks:
        for s in b.samples:
            histograms['hist_E_relTrkIso_' + s.name] = r.TH1F("hist_E_relTrkIso_" + s.name, "", 50, 0, 1)
            histograms['hist_E_relTrkIso_log_' + s.name] = r.TH1F("hist_E_relTrkIso_log_" + s.name, "", len(bin_iso_log) -1, bin_iso_log)
            histograms['hist_E_relPFIso_' + s.name] = r.TH1F("hist_E_relPFIso_" + s.name, "", 50, 0, 1)
            histograms['hist_E_relPFIso_log_' + s.name] = r.TH1F("hist_E_relPFIso_log_" + s.name, "", len(bin_iso_log) -1, bin_iso_log)
            histograms['hist_M_relPFIso_' + s.name] = r.TH1F("hist_M_relPFIso_" + s.name, "", 50, 0, 1)
            histograms['hist_M_relPFIso_log_' + s.name] = r.TH1F("hist_M_relPFIso_log_" + s.name, "", len(bin_iso_log) -1, bin_iso_log)
            for t in s.ttrees:

                print('New tree with:', t.GetEntries())

                for e,ev in enumerate(t):

                    #if e > 100000: break

                    for j in range(0, ev.nElectronCandidate):
                        if abs(ev.ElectronCandidate_eta[j]) > 2: continue
                        if abs(ev.ElectronCandidate_et[j]) < 25: continue
                        histograms['hist_E_relTrkIso_' + s.name].Fill(ev.ElectronCandidate_relTrkiso[j])
                        histograms['hist_E_relTrkIso_log_' + s.name].Fill(ev.ElectronCandidate_relTrkiso[j])
                        histograms['hist_E_relPFIso_' + s.name].Fill(ev.ElectronCandidate_relPFiso[j])
                        histograms['hist_E_relPFIso_log_' + s.name].Fill(ev.ElectronCandidate_relPFiso[j])

                        for j in range(0, ev.nDGM):
                            if abs(ev.DGM_eta[j]) > 2: continue
                            if abs(ev.DGM_pt[j]) < 30: continue
                            histograms['hist_M_relPFIso_' + s.name].Fill(ev.DGM_relPFiso[j])
                            histograms['hist_M_relPFIso_log_' + s.name].Fill(ev.DGM_relPFiso[j])


    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f_signals.root', 'RECREATE')


    #### Write everything to use later:
    for key in histograms.keys():
        histograms[key].Write()
    outputFile.Close()



