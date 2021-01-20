
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
    parser.add_option('-n', '--nmax', action='store', type=int, dest='nmax', default=0, help='Path to file')
    (opts, args) = parser.parse_args()


    ##################################
    ####   Variable declaration   ####
    ##################################
    MAX_DELTAR = 0.2
    logbin_dxy = np.logspace(-3, 3, 70)
    logbin_Lxy = np.logspace(-3, 3, 70)
    pt_bin = np.linspace(0.0, 300.0, 61)
    eta_bin = np.linspace(-2.4, 2.4, 61)
    dxy_bin = np.linspace(0.0, 10, 101)
    Lxy_bin = np.linspace(0.0, 70, 101)

    #### -----------------
    #### ---- Histograms
    #### -----------------
    eff_IFCA_pt = r.TEfficiency("eff_IFCA_pt", ";Generated electron p_{T} (GeV);Efficiency", len(pt_bin)-1, pt_bin)
    eff_IFCA_eta = r.TEfficiency("eff_IFCA_eta", ";Generated electron #eta;Efficiency", len(eta_bin)-1, eta_bin)
    eff_IFCA_dxy = r.TEfficiency("eff_IFCA_dxy", ";Generated electron |d_{xy}| (cm);Efficiency", len(dxy_bin)-1, dxy_bin)
    eff_IFCA_dxy_log = r.TEfficiency("eff_IFCA_dxy_log", ";Generated electron |d_{xy}| (cm);Efficiency", len(logbin_dxy)-1, logbin_dxy)
    eff_IFCA_Lxy = r.TEfficiency("eff_IFCA_Lxy", ";Generated electron L_{xy} (cm);Efficiency", len(Lxy_bin)-1, Lxy_bin)
    eff_IFCA_Lxy_log = r.TEfficiency("eff_IFCA_Lxy_log", ";Generated electron L_{xy} (cm);Efficiency", len(logbin_Lxy)-1, logbin_Lxy)


    hist_IFCA_dxyRes = r.TH1F("hist_IFCA_dxyRes", ";(d_{0}^{reco} - d_{0}^{gen})/d_{0}^{gen};Electron yield", 40, -1, 1)
    hist_IFCA_ptRes = r.TH1F("hist_IFCA_ptRes", ";(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen};Electron yield", 41, -0.2, 0.2)

    fake_IFCA_pt = r.TEfficiency("fake_IFCA_pt", ";Reconstructed electron p_{T} (GeV);Fake fraction", len(pt_bin)-1, pt_bin)
    fake_IFCA_eta = r.TEfficiency("fake_IFCA_eta", ";Reconstructed electron #eta;Fake fraction", len(eta_bin)-1, eta_bin)
    fake_IFCA_dxy = r.TEfficiency("fake_IFCA_dxy", ";Reconstructed electron |d_{xy}| (cm);Fake fraction", len(dxy_bin)-1, dxy_bin)
    fake_IFCA_dxy_log = r.TEfficiency("fake_IFCA_dxy_log", ";Reconstructed electron |d_{xy}| (cm);Fake fraction", len(logbin_dxy)-1, logbin_dxy)

    eff_CMS_pt = r.TEfficiency("eff_CMS_pt", ";Generated electron p_{T} (GeV);Efficiency", len(pt_bin)-1, pt_bin)
    eff_CMS_eta = r.TEfficiency("eff_CMS_eta", ";Generated electron #eta;Efficiency", len(eta_bin)-1, eta_bin)
    eff_CMS_dxy = r.TEfficiency("eff_CMS_dxy", ";Generated electron |d_{xy}| (cm);Efficiency", len(dxy_bin)-1, dxy_bin)
    eff_CMS_dxy_log = r.TEfficiency("eff_CMS_dxy_log", ";Generated electron |d_{xy}| (cm);Efficiency", len(logbin_dxy)-1, logbin_dxy)
    eff_CMS_Lxy = r.TEfficiency("eff_CMS_Lxy", ";Generated electron L_{xy} (cm);Efficiency", len(Lxy_bin)-1, Lxy_bin)
    eff_CMS_Lxy_log = r.TEfficiency("eff_CMS_Lxy_log", ";Generated electron L_{xy} (cm);Efficiency", len(logbin_Lxy)-1, logbin_Lxy)


    hist_CMS_dxyRes = r.TH1F("hist_CMS_dxyRes", ";(d_{0}^{reco} - d_{0}^{gen})/d_{0}^{gen};Electron yield", 40, -1, 1)
    hist_CMS_ptRes = r.TH1F("hist_CMS_ptRes", ";(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen};Electron yield", 41, -0.2, 0.2)

    qeff_IFCA_pt = r.TEfficiency("qeff_IFCA_pt", ";Generated electron p_{T} (GeV); Charge assignment efficiency", len(pt_bin)-1, pt_bin)
    qeff_IFCA_Lxy = r.TEfficiency("qeff_IFCA_Lxy", ";Generated electron p_{T} (GeV); Charge assignment efficiency", len(Lxy_bin)-1, Lxy_bin)


    #########################
    ####   Load sample   ####
    #########################

    _dirName = opts.filename
    _tree = r.TChain('Events')
    for _file in os.listdir(_dirName):
        if '.root' not in _file: continue
        _tree.Add(_dirName + _file)

    print("TTree with " + str(_tree.GetEntries()) + " entries")

    ###################################
    ####   Loop over tree events   ####
    ###################################

    for i in range(0, _tree.GetEntries()):

        _tree.GetEntry(i)
        if opts.nmax and i > opts.nmax: break


        #### -----------------------
        #### ----
        #### ---- Electron channel
        #### ----
        #### -----------------------

        if _tree.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 == 1:

            #### -----------------------------------
            #### ---- Efficiencies and resolutions
            #### -----------------------------------

            for j in range(0, _tree.nGenLepton):

                if not _tree.GenLeptonSel_isPromptFinalState[j]: continue
                if abs(_tree.GenLeptonSel_pdgId[j]) != 11: continue
                if (abs(_tree.GenLeptonSel_eta[j]) > 1.4442 and abs(_tree.GenLeptonSel_eta[j]) < 1.566) or abs(_tree.GenLeptonSel_eta[j]) > 2.4: continue #2.4 antes
                if _tree.GenLeptonSel_pt[j] < 28: continue

                Lxy = math.sqrt((_tree.GenLeptonSel_vx[j])**2 + (_tree.GenLeptonSel_vy[j])**2)
                pt = _tree.GenLeptonSel_pt[j]
                et = _tree.GenLeptonSel_et[j]
                eta = _tree.GenLeptonSel_eta[j]
                phi = _tree.GenLeptonSel_phi[j]
                dxy = abs(_tree.GenLeptonSel_dxy[j])

                l = TVector3()
                l.SetPtEtaPhi(pt, eta, phi)

                ll = TVector3()
                partner = False
                for jj in range(0, _tree.nGenLepton):
                    if jj == j: continue
                    if _tree.GenLeptonSel_vx[j] == _tree.GenLeptonSel_vx[jj] and _tree.GenLeptonSel_vy[j] == _tree.GenLeptonSel_vy[jj]:
                        ll.SetPtEtaPhi(_tree.GenLeptonSel_pt[jj], _tree.GenLeptonSel_eta[jj], _tree.GenLeptonSel_phi[jj])
                        partner = True
                        break

                if not partner: continue
                lldR = l.DeltaR(ll)


                ####################
                ## IFCA electrons ##
                ####################
                deltaR = 9999.0
                index = -9
                for k in range(0, _tree.nElectronCandidate):
                    if _tree.ElectronCandidate_pt[k] < 28: continue
                    if abs(_tree.ElectronCandidate_eta[k]) > 1.4442 and abs(_tree.ElectronCandidate_eta[k]) < 1.566: continue
                    re = TVector3()
                    re.SetPtEtaPhi(_tree.ElectronCandidate_pt[k], _tree.ElectronCandidate_eta[k], _tree.ElectronCandidate_phi[k])
                    if re.DeltaR(l) < deltaR:
                        deltaR = re.DeltaR(l)
                        index = k

                if deltaR < MAX_DELTAR:
                    eff_IFCA_pt.Fill(True, pt)
                    eff_IFCA_eta.Fill(True, eta)
                    eff_IFCA_dxy.Fill(True, dxy)
                    eff_IFCA_dxy_log.Fill(True, dxy)
                    eff_IFCA_Lxy.Fill(True, Lxy)
                    eff_IFCA_Lxy_log.Fill(True, Lxy)
                else:
                    eff_IFCA_pt.Fill(False, pt)
                    eff_IFCA_eta.Fill(False, eta)
                    eff_IFCA_dxy.Fill(False, dxy)
                    eff_IFCA_dxy_log.Fill(False, dxy)
                    eff_IFCA_Lxy.Fill(False, Lxy)
                    eff_IFCA_Lxy_log.Fill(False, Lxy)

                 
                #
                # -- Resolutions
                # 
                if deltaR < MAX_DELTAR:
                    ptRes = -(pt - _tree.ElectronCandidate_pt[index])/pt
                    dxyRes = -(dxy - abs(_tree.ElectronCandidate_dxy[index]))/abs(dxy)

                    hist_IFCA_ptRes.Fill(ptRes)
                    hist_IFCA_dxyRes.Fill(dxyRes)

                    if _tree.IsoTrackSel_charge[_tree.ElectronCandidate_isotrackIdx[index]]*_tree.GenLeptonSel_pdgId[j] < 0:
                        qeff_IFCA_pt.Fill(True, pt)
                        qeff_IFCA_Lxy.Fill(True, Lxy)
                    else:
                        qeff_IFCA_pt.Fill(False, pt)
                        qeff_IFCA_Lxy.Fill(False, Lxy)


                ###################
                ## CMS electrons ##
                ###################
                deltaR = 9999.0
                index = -9
                for k in range(0, _tree.nElectron):
                    if _tree.ElectronSel_pt[k] < 28: continue
                    if abs(_tree.ElectronSel_eta[k]) > 1.4442 and abs(_tree.ElectronSel_eta[k]) < 1.566: continue
                    if not _tree.ElectronSel_isLoose[k]: continue

                    re = TVector3()
                    re.SetPtEtaPhi(_tree.ElectronSel_pt[k], _tree.ElectronSel_eta[k], _tree.ElectronSel_phi[k])
                    if re.DeltaR(l) < deltaR:
                        deltaR = re.DeltaR(l)
                        index = k

                if deltaR < MAX_DELTAR:
                    eff_CMS_pt.Fill(True, pt)
                    eff_CMS_eta.Fill(True, eta)
                    eff_CMS_dxy.Fill(True, dxy)
                    eff_CMS_dxy_log.Fill(True, dxy)
                    eff_CMS_Lxy.Fill(True, Lxy)
                    eff_CMS_Lxy_log.Fill(True, Lxy)
                else:
                    eff_CMS_pt.Fill(False, pt)
                    eff_CMS_eta.Fill(False, eta)
                    eff_CMS_dxy.Fill(False, dxy)
                    eff_CMS_dxy_log.Fill(False, dxy)
                    eff_CMS_Lxy.Fill(False, Lxy)
                    eff_CMS_Lxy_log.Fill(False, Lxy)

                 
                #
                # -- Resolutions
                # 
                if deltaR < MAX_DELTAR:
                    ptRes = -(pt - _tree.ElectronSel_pt[index])/pt
                    dxyRes = -(dxy - abs(_tree.ElectronSel_dB[index]))/abs(dxy)

                    hist_CMS_ptRes.Fill(ptRes)
                    hist_CMS_dxyRes.Fill(dxyRes)


            #### -----------------
            #### ---- Fake rates
            #### -----------------

            for j in range(0, _tree.nElectronCandidate):

                pt = _tree.ElectronCandidate_pt[j]
                et = _tree.ElectronCandidate_et[j]
                eta = _tree.ElectronCandidate_eta[j]
                phi = _tree.ElectronCandidate_phi[j]
                dxy = abs(_tree.ElectronCandidate_dxy[j])

                l = TVector3()
                l.SetPtEtaPhi(pt, eta, phi)

                ####################
                ## Gen electrons ##
                ####################
                deltaR = 9999.0
                index = -9
                for k in range(0, _tree.nGenLepton):
                    if not _tree.GenLeptonSel_isPromptFinalState[k]: continue
                    if abs(_tree.GenLeptonSel_pdgId[k]) != 11: continue
                    re = TVector3()
                    re.SetPtEtaPhi(_tree.GenLeptonSel_pt[k], _tree.GenLeptonSel_eta[k], _tree.GenLeptonSel_phi[k])
                    if re.DeltaR(l) < deltaR:
                        deltaR = re.DeltaR(l)
                        index = k

                if deltaR < MAX_DELTAR:
                    fake_IFCA_pt.Fill(False, pt)
                    fake_IFCA_eta.Fill(False, eta)
                    fake_IFCA_dxy.Fill(False, dxy)
                    fake_IFCA_dxy_log.Fill(False, dxy)
                else:
                    fake_IFCA_pt.Fill(True, pt)
                    fake_IFCA_eta.Fill(True, eta)
                    fake_IFCA_dxy.Fill(True, dxy)
                    fake_IFCA_dxy_log.Fill(True, dxy)


    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    eff_IFCA_pt.Write()
    eff_IFCA_eta.Write()
    eff_IFCA_dxy.Write()
    eff_IFCA_dxy_log.Write()
    eff_IFCA_Lxy.Write()
    eff_IFCA_Lxy_log.Write()
    hist_IFCA_ptRes.Write()
    hist_IFCA_dxyRes.Write()
    fake_IFCA_pt.Write()
    fake_IFCA_eta.Write()
    fake_IFCA_dxy.Write()
    fake_IFCA_dxy_log.Write()
    qeff_IFCA_pt.Write()
    qeff_IFCA_Lxy.Write()
    eff_CMS_pt.Write()
    eff_CMS_eta.Write()
    eff_CMS_dxy.Write()
    eff_CMS_dxy_log.Write()
    eff_CMS_Lxy.Write()
    eff_CMS_Lxy_log.Write()
    hist_CMS_ptRes.Write()
    hist_CMS_dxyRes.Write()

    outputFile.Close()


