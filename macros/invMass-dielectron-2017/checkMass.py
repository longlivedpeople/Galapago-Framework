
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
    hSR_fullPath = r.TH1F("hSR_fullPath", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hBCR_fullPath = r.TH1F("hBCR_fullPath", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hSR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 = r.TH1F("hSR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hBCR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 = r.TH1F("hBCR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hSR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 = r.TH1F("hSR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hBCR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 = r.TH1F("hBCR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hSR_HLT_DoublePhoton70 = r.TH1F("hSR_HLT_DoublePhoton70", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)
    hBCR_HLT_DoublePhoton70 = r.TH1F("hBCR_HLT_DoublePhoton70", ";Dielectron invariant mass m_{ee} (GeV); Events / 5.00 units", 20, 15, 115)

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

        #if i > 10000: break

        #### -------------------
        
        #### ---- Muon channel
        #### ----
        #### -------------------

        if ev.nEE > 0:

            #### -----------------
            #### ---- Get EE
            #### -----------------
            ee_maxIxy = -99
            maxIxy = -1
            nBSEE = 0

	    for j in range(0, ev.nEE):
                iee = j
                if not eval(cm.EE_sel): continue
                if not eval(cm.EE_etanoBE): continue
                nBSEE+=1
                if ev.EE_trackIxy_PV[j] > maxIxy:
                    maxIxy = ev.EE_trackIxy_PV[j]
                    ee_maxIxy = j

            if ee_maxIxy < 0 or nBSEE < 1: continue

            iee = ee_maxIxy

            if eval(cm.AddList([cm.EE_Ixy6prompt, cm.EE_dPhiforward, cm.EE_nBSEEe1])):

                if ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 or ev.HLT_DoublePhoton70:
                     hSR_fullPath.Fill(ev.EE_mass[iee])

                if ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90:
                     hSR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90.Fill(ev.EE_mass[iee])

                if ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55:
                     hSR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55.Fill(ev.EE_mass[iee])

                if ev.HLT_DoublePhoton70:
                     hSR_HLT_DoublePhoton70.Fill(ev.EE_mass[iee])


            if eval(cm.AddList([cm.EE_Ixy6prompt, cm.EE_dPhibackward, cm.EE_nBSEEe1])):

                if ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 or ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55 or ev.HLT_DoublePhoton70:
                     hBCR_fullPath.Fill(ev.EE_mass[iee])

                if ev.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90:
                     hBCR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90.Fill(ev.EE_mass[iee])

                if ev.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55:
                     hBCR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55.Fill(ev.EE_mass[iee])

                if ev.HLT_DoublePhoton70:
                     hBCR_HLT_DoublePhoton70.Fill(ev.EE_mass[iee])


    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    hSR_fullPath.Write()
    hSR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90.Write()
    hSR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55.Write()
    hSR_HLT_DoublePhoton70.Write()
    hBCR_fullPath.Write()
    hBCR_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90.Write()
    hBCR_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55.Write()
    hBCR_HLT_DoublePhoton70.Write()
    outputFile.Close()


