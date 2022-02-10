
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

print(WORKPATH, WORKPATH)
print(GALAPAGOPATH, GALAPAGOPATH)


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
    (opts, args) = parser.parse_args()


    ##################################
    ####   Variable declaration   ####
    ##################################
    MAX_DELTAR = 0.2
    logbin_iso = np.logspace(-3, 2, 100)


    #### -----------------
    #### ---- Histograms
    #### -----------------
    hEE_vxvy = r.TH2F("EE_vxvy", ";Vertex x coordinate (cm); Vertex y coordinate (cm)", 100, -1, 1, 100, -1, 1)
    hEE_vxIvyI = r.TH2F("EE_vxIvyI", ";Vertex x_{PV} coordinate (cm); Vertex y_{PV} coordinate (cm)", 100, -0.01, 0.01, 100, -0.01, 0.01)
    hEE_vxIIvyII = r.TH2F("EE_vxIIvyII", ";Vertex x_{ee} coordinate (cm); Vertex y_{ee} coordinate (cm)", 100, -0.01, 0.01, 100, -0.01, 0.01)
    hEE_vx1vy1 = r.TH2F("EE_vx1vy1", ";Vertex x_{e1} coordinate (cm); Vertex y_{e1} coordinate (cm)", 100, -0.01, 0.01, 100, -0.01, 0.01)
    hEE_vx2vy2 = r.TH2F("EE_vx2vy2", ";Vertex x_{e2} coordinate (cm); Vertex y_{e2} coordinate (cm)", 100, -0.01, 0.01, 100, -0.01, 0.01)
    hEE_vxII = r.TH1F("EE_vxII", ";Vertex x'' coordinate (cm); EE yield", 100, -0.01, 0.01)
    hEE_vyII = r.TH1F("EE_vyII", ";Vertex y'' coordinate (cm); EE yield", 100, -0.01, 0.01)
    hEE_isodPhi = r.TH1F("EE_isodPhi", ";Collinearity |#Delta#Phi|; EE yield", 50, 0.0, 3.14)
    hEE_trackdPhi = r.TH1F("EE_trackdPhi", ";Collinearity |#Delta#Phi|; EE yield", 50, 0.0, 3.14)
    hEE_llangle = r.TH1F("EE_llangle", ";Lepton alpha angle; EE yield", 50, 0.0, 3.14)

    hMM_vxvy = r.TH2F("MM_vxvy", ";Vertex x coordinate (cm); Vertex y coordinate (cm)", 100, -1, 1, 100, -1, 1)
    hMM_vxIvyI = r.TH2F("MM_vxIvyI", ";Vertex x_{PV} coordinate (cm); Vertex y_{PV} coordinate (cm)", 100, -0.01, 0.01, 100, -0.01, 0.01)
    hMM_vxIIvyII = r.TH2F("MM_vxIIvyII", ";Vertex x_{#mu#mu} coordinate (cm); Vertex y_{#mu#mu} coordinate (cm)", 100, -0.004, 0.004, 100, -0.004, 0.004)
    hMM_vx1vy1 = r.TH2F("MM_vx1vy1", ";Vertex x_{#mu1} coordinate (cm); Vertex y_{#mu1} coordinate (cm)", 100, -0.004, 0.004, 100, -0.004, 0.004)
    hMM_vx2vy2 = r.TH2F("MM_vx2vy2", ";Vertex x_{#mu2} coordinate (cm); Vertex y_{#mu2} coordinate (cm)", 100, -0.004, 0.004, 100, -0.004, 0.004)
    hMM_vxII = r.TH1F("MM_vxII", ";Vertex x'' coordinate (cm); MM yield", 100, -0.004, 0.004)
    hMM_vyII = r.TH1F("MM_vyII", ";Vertex y'' coordinate (cm); MM yield", 100, -0.004, 0.004)
    hMM_dPhi = r.TH1F("MM_dPhi", ";Collinearity |#Delta#Phi|; MM yield", 50, 0.0, 3.14)
    hMM_llangle = r.TH1F("MM_llangle", ";Lepton alpha angle; MM yield", 50, 0.0, 3.14)

    #########################
    ####   Load sample   ####
    #########################
    _sampleNames = (opts.filename).split(',')
    _tree = r.TChain('Events')
    for _name in _sampleNames:
        _tree.Add(_name)
    print("TTree with " + str(_tree.GetEntries()) + " entries")


    ###################################
    ####   Loop over tree events   ####
    ###################################
    cm = CutManager.CutManager()
    for i,ev in enumerate(_tree):

        #if i > 1000000: break

        #### -------------------
        #### ----
        #### ---- Muon channel
        #### ----
        #### -------------------

        if eval(cm.mupath2018) and ev.nDMDM > 0:

            #### -----------------
            #### ---- Get DMDM
            #### -----------------
	    for j in range(0, ev.nDMDM):
                imm = j
                if not eval(cm.MM_BS2018): continue
                if not eval(cm.MM_iso2l): continue
                if not eval(cm.MM_OS): continue

                #### -----------------------------
                #### ---- Change of coordinates
                #### -----------------------------

                l1 = TVector3()
                l2 = TVector3()

                if ev.DGM_pt[ev.DMDM_idxA[imm]] > ev.DGM_pt[ev.DMDM_idxB[imm]]:
                    l1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[imm]], ev.DGM_eta[ev.DMDM_idxA[imm]], ev.DGM_phi[ev.DMDM_idxA[imm]])
                    l2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[imm]], ev.DGM_eta[ev.DMDM_idxB[imm]], ev.DGM_phi[ev.DMDM_idxB[imm]])
                else:
                    l2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[imm]], ev.DGM_eta[ev.DMDM_idxA[imm]], ev.DGM_phi[ev.DMDM_idxA[imm]])
                    l1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[imm]], ev.DGM_eta[ev.DMDM_idxB[imm]], ev.DGM_phi[ev.DMDM_idxB[imm]])

                l1l2 = l1 + l2
                phill = -l1l2.Phi()
                alpha = abs(l1.DeltaPhi(l2))

                vx = ev.DMDM_vx[imm]
                vy = ev.DMDM_vy[imm]

                vxI = vx - ev.PV_vx
                vyI = vy - ev.PV_vy

                vec = TVector3(vx - ev.PV_vx, vy - ev.PV_vy, 0.0)
                dPhi = abs(l1l2.DeltaPhi(vec))

                #print(dPhi, ev.DMDM_dPhi[imm])

                # l1 + l2 oriented:
                vxII = vxI*math.cos(phill) - vyI*math.sin(phill)
                vyII = vxI*math.sin(phill) + vyI*math.cos(phill)

                # l1 oriented:
                phi1 = l1.Phi()
                vx1 = vxI*math.cos(phi1) + vyI*math.sin(phi1)
                vy1 = -vxI*math.sin(phi1) + vyI*math.cos(phi1)
            
                # l2 oriented:
                phi2 = l2.Phi()
                vx2 = vxI*math.cos(phi2) + vyI*math.sin(phi2)
                vy2 = -vxI*math.sin(phi2) + vyI*math.cos(phi2)



                 #### -------------------------
                 #### ---- Filling histograms
                 #### -------------------------
                if abs(phill) < 3.0:
                    hMM_vxvy.Fill(vx, vy)
                    hMM_vxIvyI.Fill(vxI, vyI)
                    hMM_vxIIvyII.Fill(vxII, vyII)
                    hMM_vxII.Fill(vxII)
                    hMM_vyII.Fill(vyII)
                    hMM_vx1vy1.Fill(vx1, vy1)
                    hMM_vx2vy2.Fill(vx2, vy2)
                    hMM_dPhi.Fill(dPhi)
                    hMM_llangle.Fill(alpha)

        #### -----------------------
        #### ----
        #### ---- Electron channel
        #### ----
        #### -----------------------

        if not eval(cm.mupath2018) and eval(cm.epath2018) and ev.nEE > 0:

            #### -----------------
            #### ---- Get EE
            #### -----------------
            for j in range(0, ev.nEE):
                iee = j
                if not eval(cm.EE_BS2018): continue
                if not eval(cm.EE_iso2l): continue
                if not eval(cm.EE_OS): continue


                #### -----------------------------
                #### ---- Change of coordinates
                #### -----------------------------

                l1 = TVector3()
                l2 = TVector3()

                if ev.ElectronCandidate_pt[ev.EE_idxA[iee]] > ev.ElectronCandidate_pt[ev.EE_idxB[iee]]:
                    l1.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxA[iee]], ev.ElectronCandidate_eta[ev.EE_idxA[iee]], ev.ElectronCandidate_phi[ev.EE_idxA[iee]])
                    l2.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxB[iee]], ev.ElectronCandidate_eta[ev.EE_idxB[iee]], ev.ElectronCandidate_phi[ev.EE_idxB[iee]])
                else:
                    l2.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxA[iee]], ev.ElectronCandidate_eta[ev.EE_idxA[iee]], ev.ElectronCandidate_phi[ev.EE_idxA[iee]])
                    l1.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxB[iee]], ev.ElectronCandidate_eta[ev.EE_idxB[iee]], ev.ElectronCandidate_phi[ev.EE_idxB[iee]])

                l1l2 = l1 + l2
                phill = -l1l2.Phi()
                alpha = abs(l1.DeltaPhi(l2))

                vx = ev.EE_vx[iee]
                vy = ev.EE_vy[iee]

                vec = TVector3(vx - ev.PV_vx, vy - ev.PV_vy, 0.0)
                dPhi = abs(l1l2.DeltaPhi(vec))

                vxI = vx - ev.PV_vx
                vyI = vy - ev.PV_vy

                # l1l2 oriented:
                vxII = vxI*math.cos(phill) - vyI*math.sin(phill)
                vyII = vxI*math.sin(phill) + vyI*math.cos(phill)

                # l1 oriented:
                phi1 = l1.Phi()
                vx1 = vxI*math.cos(phi1) + vyI*math.sin(phi1)
                vy1 = -vxI*math.sin(phi1) + vyI*math.cos(phi1)
            
                # l2 oriented:
                phi2 = l2.Phi()
                vx2 = vxI*math.cos(phi2) + vyI*math.sin(phi2)
                vy2 = -vxI*math.sin(phi2) + vyI*math.cos(phi2)



                #### -------------------------
                #### ---- Filling histograms
                #### -------------------------
                hEE_vxvy.Fill(vx, vy)
                hEE_vxIvyI.Fill(vxI, vyI)
                hEE_vx1vy1.Fill(vx1, vy1)
                hEE_vx2vy2.Fill(vx2, vy2)
                hEE_vxIIvyII.Fill(vxII, vyII)
                hEE_vxII.Fill(vxII)
                hEE_vyII.Fill(vyII)
                hEE_isodPhi.Fill(dPhi)
                hEE_trackdPhi.Fill(abs(ev.EE_dPhi[iee]))
                hEE_llangle.Fill(alpha)


    if not os.path.exists(WORKPATH + 'Vertex-results/'): os.makedirs(WORKPATH + 'Vertex-results/')
    outputFile = TFile(WORKPATH + 'Vertex-results/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    hEE_vxvy.Write()
    hEE_vxIvyI.Write()
    hEE_vxIIvyII.Write()
    hEE_vx1vy1.Write()
    hEE_vx2vy2.Write()
    hEE_vxII.Write()
    hEE_vyII.Write()
    hEE_isodPhi.Write()
    hEE_trackdPhi.Write()
    hEE_llangle.Write()

    hMM_vxvy.Write()
    hMM_vxIvyI.Write()
    hMM_vxIIvyII.Write()
    hMM_vx1vy1.Write()
    hMM_vx2vy2.Write()
    hMM_vxII.Write()
    hMM_vyII.Write()
    hMM_dPhi.Write()
    hMM_llangle.Write()

    outputFile.Close()


