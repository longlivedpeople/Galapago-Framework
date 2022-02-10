
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
    parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    parser.add_option('-f', '--filename', action='store', type=str, dest='filename', default='', help='Path to file')
    parser.add_option('-n', '--nmax', action='store', type=int, dest='nmax', default=0, help='Path to file')
    (opts, args) = parser.parse_args()


    ##################################
    ####   Variable declaration   ####
    ##################################

    #### -----------------
    #### ---- Histograms
    #### -----------------
    hist_counts = r.TH1F("hist_counts", "", 4, 0, 4)
    hist_cosAlpha = r.TH1F("hist_cosAlpha", "", 100, -1.0, 1.0)
    hist_alpha = r.TH1F("hist_Alpha", "", 100, 0.0, 3.15)
    hist_deltaPhi = r.TH1F("hist_deltaPhi", "", 100, 0.0, 3.15)
    hist_dPhi = r.TH1F("hist_dPhi", "", 100, 0.0, 3.15)
    hist_normChi2 = r.TH1F("hist_normChi2", "", 100, 0.0, 30.0)

    
    #########################
    ####   Load sample   ####
    #########################
    """
    _dirName = opts.filename
    _tree = r.TChain('Events')
    for _file in os.listdir(_dirName):
        if '.root' not in _file: continue
        _tree.Add(_dirName + _file)
    """
    _file = r.TFile(opts.filename)
    _tree = _file.Get('Events')
    print("TTree with " + str(_tree.GetEntries()) + " entries")


    ###################################
    ####   Loop over tree events   ####
    ###################################
    cm = CutManager.CutManager()

    for i, ev in enumerate(_tree):

        if opts.nmax and i > opts.nmax: break
        hist_counts.Fill(0)
        if not ev.HLT_L2Mu10_NoVertex_NoBPTX3BX: continue
        hist_counts.Fill(1)
        if ev.nDMDM < 1: continue
        hist_counts.Fill(2)

        qual = False
        for j in range(0, ev.nDMDM):

            imm = j # index to handle DMDM pair
            if not ev.DMDM_mass[imm] > 15: continue
            if not ev.DMDM_normalizedChi2[imm] < 10: continue
            if not ev.DMDM_dR[imm] > 0.2: continue
            if not ev.DGM_pt[ev.DMDM_idxA[imm]] > 25 or not ev.DGM_pt[ev.DMDM_idxB[imm]] > 25: continue
            if not abs(ev.DGM_eta[ev.DMDM_idxA[imm]]) < 2.0 or not abs(ev.DGM_eta[ev.DMDM_idxB[imm]]) < 2.0: continue
            if not ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] < 0: continue
            if not eval(cm.MM_iso2l): continue
            
            #ID:
            if not ev.DGM_ptError[ev.DMDM_idxB[imm]]/ev.DGM_pt[ev.DMDM_idxB[imm]] < 0.3: continue
            if not ev.DGM_ptError[ev.DMDM_idxA[imm]]/ev.DGM_pt[ev.DMDM_idxA[imm]] < 0.3: continue
            if not ev.DGM_normChi2[ev.DMDM_idxA[imm]] < 7.5: continue
            if not ev.DGM_normChi2[ev.DMDM_idxB[imm]] < 7.5: continue
            if not ev.DGM_muonHits[ev.DMDM_idxA[imm]] > 11: continue  
            if not ev.DGM_muonHits[ev.DMDM_idxB[imm]] > 11: continue  
            if not ev.DGM_outerTrackerHits[ev.DMDM_idxA[imm]] > 8: continue  
            if not ev.DGM_outerTrackerHits[ev.DMDM_idxB[imm]] > 8: continue  

            qual = True

            #if not eval(cm.MM_ID): continue
            #if not eval(cm.MM_cosAlpha0p8): continue

            dgl1 = r.TVector3()
            dgl2 = r.TVector3()
            dgl1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[j]], ev.DGM_eta[ev.DMDM_idxA[j]], ev.DGM_phi[ev.DMDM_idxA[j]])
            dgl2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[j]], ev.DGM_eta[ev.DMDM_idxB[j]], ev.DGM_phi[ev.DMDM_idxB[j]])

            #print("Eta: ", dgl1.Eta(), dgl2.Eta(), "Theta: " ,dgl1.Theta(), dgl2.Theta())

            hist_cosAlpha.Fill(ev.DMDM_cosAlpha[j])
            hist_alpha.Fill(abs(dgl1.Angle(dgl2)))
            hist_deltaPhi.Fill(abs(dgl1.DeltaPhi(dgl2)))
            hist_normChi2.Fill(ev.DMDM_normalizedChi2[j])
            hist_dPhi.Fill(abs(ev.DMDM_dPhi[j]))
             
        if qual:
            hist_counts.Fill(3)

    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    hist_counts.Write()
    hist_cosAlpha.Write()
    hist_deltaPhi.Write()
    hist_alpha.Write()
    hist_normChi2.Write()
    hist_dPhi.Write()
    outputFile.Close()



