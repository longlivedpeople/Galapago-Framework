
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
    hist_E_relTrkIso     = r.TH1F("hist_E_relTrkIso", "", 50, 0, 1)
    hist_E_relTrkIso_log = r.TH1F("hist_E_relTrkIso_log", "", len(bin_iso_log) -1, bin_iso_log)
    hist_E_relPFIso      = r.TH1F("hist_E_relPFIso", "", 50, 0, 1)
    hist_E_relPFIso_log  = r.TH1F("hist_E_relPFIso_log", "", len(bin_iso_log) -1, bin_iso_log)
    hist_M_relPFIso      = r.TH1F("hist_M_relPFIso", "", 50, 0, 1)
    hist_M_relPFIso_log  = r.TH1F("hist_M_relPFIso_log", "", len(bin_iso_log) -1, bin_iso_log)

    #########################
    ####   Load sample   ####
    #########################

    treeBKG = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['WW_postVFP', 'WZ_postVFP'], 'MC'), name = '2016', isdata = 0 )


    ###################################
    ####   Loop over tree events   ####
    ###################################
    cm = CutManager.CutManager()

    for b in treeBKG.blocks:
        for s in b.samples:
            for t in s.ttrees:

                print('New tree with:', t.GetEntries())

                for e,ev in enumerate(t):

                    #if e > 500000: break
                    #if not eval(cm.epath2016): continue

                    # Count number of hard process leptons:
                    hasLepton = False
                    for j in range(0, ev.nHardProcessParticle):
                        if (abs(ev.HardProcessParticle_pdgId[j]) == 11) or (abs(ev.HardProcessParticle_pdgId[j]) == 13) or (abs(ev.HardProcessParticle_pdgId[j]) == 15):
                            hasLepton = True
                            break
                    if hasLepton: continue

                    for j in range(0, ev.nElectronCandidate):
                        if abs(ev.ElectronCandidate_eta[j]) > 2: continue
                        if abs(ev.ElectronCandidate_et[j]) < 25: continue
                        hist_E_relTrkIso.Fill(ev.ElectronCandidate_relTrkiso[j])
                        hist_E_relTrkIso_log.Fill(ev.ElectronCandidate_relTrkiso[j])
                        hist_E_relPFIso.Fill(ev.ElectronCandidate_relPFiso[j])
                        hist_E_relPFIso_log.Fill(ev.ElectronCandidate_relPFiso[j])

                    for j in range(0, ev.nDGM):
                        if abs(ev.DGM_eta[j]) > 2: continue
                        if abs(ev.DGM_pt[j]) < 30: continue
                        hist_M_relPFIso.Fill(ev.DGM_relPFiso[j])
                        hist_M_relPFIso_log.Fill(ev.DGM_relPFiso[j])


    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f.root', 'RECREATE')


    #### Write everything to use later:
    hist_E_relTrkIso.Write()
    hist_E_relTrkIso_log.Write()
    hist_E_relPFIso.Write()
    hist_E_relPFIso_log.Write()
    hist_M_relPFIso.Write()
    hist_M_relPFIso_log.Write()
    outputFile.Close()



