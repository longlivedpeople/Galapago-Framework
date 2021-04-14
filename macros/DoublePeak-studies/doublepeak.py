
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
    iso_bin = np.logspace(-3, 1, 60)
    iso_bin_lin = np.linspace(0.1, 1.2, 50)
    et_bin = np.linspace(20, 100, 50)

    #### -----------------
    #### ---- Histograms
    #### -----------------
    hist_dR = r.TH1F("hist_dR", ";dR; Generated electron yield", 40, 0, 1.0)
    hist_dPhi = r.TH1F("hist_dPhi", ";#Delta#Phi(X, electron); Generated electron yield", 40, -3.14, 3.14)
    hist_relPFisowithE = r.TH1F("hist_relPFisowithE", ";PFiso / SC p_{T}; Displaced electron yield", len(iso_bin)-1, iso_bin)
    hist_neutralsumwithE = r.TH1F("hist_neutralsumwithE", ";NeutralSum / SC p_{T}; Displaced electron yield", len(iso_bin)-1, iso_bin)
    hist_neutralsum_et = r.TH2F("hist_neutralsum_Et", ";NeutralSum / SC p_{T}; SC p_{T}", len(iso_bin)-1, iso_bin, len(et_bin)-1, et_bin)
    hist_neutralsum_pt = r.TH2F("hist_neutralsum_pt", ";NeutralSum / SC p_{T}; isoltrack p_{T}", len(iso_bin)-1, iso_bin, len(et_bin)-1, et_bin)
    hist_neutralsum_ptdet = r.TH2F("hist_neutralsum_ptdet", ";NeutralSum / SC p_{T}; isotrack p_{T} / SC p_{T}", len(iso_bin)-1, iso_bin, len(np.linspace(0, 1.5, 50)) - 1, np.linspace(0, 1.5, 50))
    hist_neutralsum_ptdet_lin = r.TH2F("hist_neutralsum_ptdet_lin", ";NeutralSum / SC p_{T}; isotrack p_{T} / SC p_{T}", len(iso_bin_lin)-1, iso_bin_lin, len(np.linspace(0, 1.5, 50)) - 1, np.linspace(0, 1.5, 50))
    hist_neutralsum_ptmet = r.TH2F("hist_neutralsum_ptmet", ";NeutralSum / SC p_{T}; isotrack p_{T} - SC p_{T}", len(iso_bin)-1, iso_bin, len(et_bin)-1, et_bin)
    hist_neutralsum_dPhi = r.TH2F("hist_neutralsum_dPhi", ";NeutralSum / SC p_{T}; #Delta#Phi(X, electron)", len(iso_bin)-1, iso_bin, len(np.linspace(0, 1, 20)) - 1, np.linspace(0, 1, 20))
   
    hist_2d = r.TH2F("hist_NeutralSum_nhits", ";neutralsum; nhits", len(iso_bin)-1, iso_bin, len(np.linspace(8, 25, 25))-1, np.linspace(8, 25, 25))
    prof_neutralsum_dPhi = r.TProfile("prof_NeutralSum_dPhi", ";NeutralSum / SC p_{T}; #Delta#Phi(X boson, displaced electron)", len(iso_bin_lin)-1, iso_bin_lin, 0.0, 0.7)
    prof_neutralsum_dR = r.TProfile("prof_NeutralSum_dR", ";NeutralSum / SC p_{T}; #DeltaR(X boson, displaced electron)", len(iso_bin_lin)-1, iso_bin_lin, 0.0, 1.0)

    
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
    _tree = _file.Get("Events")

    print("TTree with " + str(_tree.GetEntries()) + " entries")


    ###################################
    ####   Loop over tree events   ####
    ###################################

    for i in range(0, _tree.GetEntries()):

        _tree.GetEntry(i)
        if opts.nmax and i > opts.nmax: break

        # count generated electrons
        iEle = []
        for n in range(0, _tree.nGenLepton):
            pdgId = _tree.GenLeptonSel_pdgId[n]
            if abs(pdgId) != 11: continue
            iEle.append(n)
            
        #if len(iEle) != 2: continue
        if len(iEle) < 2: continue

        # compute dR between generated electrons
        l1 = r.TVector3()
        l2 = r.TVector3()
        i1 = iEle[0]
        i2 = iEle[1]
        l1.SetPtEtaPhi(_tree.GenLeptonSel_pt[i1], _tree.GenLeptonSel_eta[i1], _tree.GenLeptonSel_phi[i1])
        l2.SetPtEtaPhi(_tree.GenLeptonSel_pt[i2], _tree.GenLeptonSel_eta[i2], _tree.GenLeptonSel_phi[i2])
        eedR = l1.DeltaR(l2)
        hist_dR.Fill(eedR)
        l12 = l1 + l2

        Lxy = math.sqrt(_tree.GenLeptonSel_vx[0]*_tree.GenLeptonSel_vx[0] + _tree.GenLeptonSel_vy[0]*_tree.GenLeptonSel_vy[0]) 

        if Lxy < 5.0: continue

        for n in range(0, _tree.nElectronCandidate):
            hist_relPFisowithE.Fill(_tree.ElectronCandidate_pt[n]/_tree.ElectronCandidate_et[n]*_tree.ElectronCandidate_relPFiso[n])
            hist_neutralsumwithE.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n])
            hist_neutralsum_et.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], _tree.ElectronCandidate_et[n])
            hist_neutralsum_pt.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], _tree.ElectronCandidate_pt[n])
            hist_neutralsum_ptdet.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], _tree.ElectronCandidate_pt[n]/_tree.ElectronCandidate_et[n])
            hist_neutralsum_ptdet_lin.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], _tree.ElectronCandidate_pt[n]/_tree.ElectronCandidate_et[n])
            hist_neutralsum_ptmet.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], _tree.ElectronCandidate_pt[n] - _tree.ElectronCandidate_et[n])
            hist_2d.Fill(_tree.ElectronCandidate_relPFiso[n]*_tree.ElectronCandidate_pt[n]/_tree.ElectronCandidate_et[n], _tree.IsoTrackSel_numberOfValidTrackerHits[_tree.ElectronCandidate_isotrackIdx[n]])
            le = r.TVector3()
            le.SetPtEtaPhi(_tree.ElectronCandidate_pt[n], _tree.ElectronCandidate_eta[n], _tree.ElectronCandidate_phi[n])
            dPhi = l12.DeltaPhi(le)
            hist_dPhi.Fill(dPhi)
            dPhi = abs(dPhi)
            dR = l12.DeltaR(le)
            hist_neutralsum_dPhi.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], dPhi)
            prof_neutralsum_dPhi.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], dPhi)
            prof_neutralsum_dR.Fill(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n], dR)
            #print(_tree.ElectronCandidate_PFiso_NeutralSumOut[n]/_tree.ElectronCandidate_et[n])



    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    hist_2d.Write()
    outputFile.Close()

    c1 = TCanvas("c1", "", 600, 500)
    c1.SetLogx(1)

    hist_relPFisowithE.Draw("hist")
    c1.SaveAs("relPFisowithE.png")

    c1.Clear()
    hist_neutralsum_et.Draw("colz")
    c1.SaveAs("neutralsum_et.png")

    c1.Clear()
    hist_neutralsumwithE.Draw("hist")
    c1.SaveAs("neutralsumwithE.png")

    c1.Clear()
    hist_neutralsum_pt.Draw("colz")
    c1.SaveAs("neutralsum_pt.png")

    c1.Clear()
    hist_neutralsum_ptdet.Draw("colz")
    c1.SaveAs("neutralsum_ptdet.png")

    c1.Clear()
    c1.SetLogx(0)
    hist_neutralsum_ptdet_lin.Draw("colz")
    c1.SaveAs("neutralsum_ptdet_lin.png")

    c1.Clear()
    hist_dPhi.Draw("hist")
    c1.SaveAs("dPhi.png")

    c1.Clear()
    c1.SetLogx(1)
    hist_neutralsum_ptmet.Draw("colz")
    c1.SaveAs("neutralsum_ptmet.png")

    c1.Clear()
    hist_neutralsum_dPhi.Draw("colz")
    c1.SaveAs("neutralsum_dPhi.png")

    c1.Clear()
    c1.SetLogx(0)
    prof_neutralsum_dPhi.Draw("P")
    c1.SaveAs("prof_neutralsum_dPhi.png")

    c1.Clear()
    c1.SetLogx(0)
    prof_neutralsum_dR.Draw("P")
    c1.SaveAs("prof_neutralsum_dR.png")

    """
    c1.Clear()
    c1.SetLogx(0)
    hist_dR.Draw("hist")
    c1.SaveAs("dR.png")
    """
