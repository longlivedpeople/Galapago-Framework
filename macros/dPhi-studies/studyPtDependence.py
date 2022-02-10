
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
    hMM_dPhi_sym = r.TH1F("MM_dPhi_sym", ";Collinearity |#Delta#Phi|; Dimuon yield", 20, 0.0, 3.14)
    hMM_dPhi_asym = r.TH1F("MM_dPhi_asym", ";Collinearity |#Delta#Phi|; Dimuon yield", 20, 0.0, 3.14)
    hEE_dPhi_sym = r.TH1F("EE_dPhi_sym", ";Collinearity |#Delta#Phi|; Dielectron yield", 20, 0.0, 3.14)
    hEE_dPhi_asym = r.TH1F("EE_dPhi_asym", ";Collinearity |#Delta#Phi|; Dielectron yield", 20, 0.0, 3.14)

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

        if i > 1000000: break

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

                if ev.DMDM_leadingPt[imm] > 40 and ev.DMDM_subleadingPt[imm] > 25:
                    hMM_dPhi_asym.Fill(abs(ev.DMDM_dPhi[imm]))

                if ev.DMDM_leadingPt[imm] > 40 and ev.DMDM_subleadingPt[imm] > 40:
                    hMM_dPhi_sym.Fill(abs(ev.DMDM_dPhi[imm]))


        #### -----------------------
        #### ----
        #### ---- Electron channel
        #### ----
        #### -----------------------

        if eval(cm.epath2018) and ev.nEE > 0:

            #### -----------------
            #### ---- Get EE
            #### -----------------
            for j in range(0, ev.nEE):
                iee = j
                if not eval(cm.EE_BS2018): continue
                if not eval(cm.EE_iso2l): continue
                if not eval(cm.EE_OS): continue

                hEE_dPhi_asym.Fill(abs(ev.EE_dPhi[iee]))

                if ev.EE_leadingEt[iee] > 40 and ev.EE_subleadingEt[iee] > 40:
                    hEE_dPhi_sym.Fill(abs(ev.EE_dPhi[iee]))


    hMM_dPhi_sym.SetMinimum(0.0)
    hMM_dPhi_sym.SetMaximum(1.5*hMM_dPhi_asym.GetMaximum())
    plot = Canvas.Canvas('MM_dPhi_ptcutsymmetry', 'png,pdf', 0.16, 0.8, 0.35, 0.89, 1, lsize = 0.035)
    plot.addHisto(hMM_dPhi_sym, 'P', 'Symmetric cuts (p_{T} > 40 GeV)', 'p', r.kBlue-3, 1, 0, marker = 20)
    plot.addHisto(hMM_dPhi_asym, 'P,SAME', 'Asymmetric cuts (p_{T} > 40, 25 GeV)', 'p', r.kMagenta-3, 1, 1, marker = 25)
    plot.addLatex(0.9, 0.93, 'Drell-Yan (2018 UL)', font = 42, size = 0.035, align = 31)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + '/', maxYnumbers = 3)

    
    hEE_dPhi_sym.SetMinimum(0.0)
    hEE_dPhi_sym.SetMaximum(1.5*hEE_dPhi_asym.GetMaximum())
    plot = Canvas.Canvas('EE_dPhi_ptcutsymmetry', 'png,pdf', 0.16, 0.8, 0.35, 0.89, 1, lsize = 0.035)
    plot.addHisto(hEE_dPhi_sym, 'P', 'Symmetric cuts (p_{T} > 40 GeV)', 'p', r.kBlue-3, 1, 0, marker = 20)
    plot.addHisto(hEE_dPhi_asym, 'P,SAME', 'Asymmetric cuts (p_{T} > 40, 25 GeV)', 'p', r.kMagenta-3, 1, 1, marker = 25)
    plot.addLatex(0.9, 0.93, 'Drell-Yan (2018 UL)', font = 42, size = 0.035, align = 31)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + '/', maxYnumbers = 3)


