import math
import os
import json
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager
from include.processHandler import processHandler
import numpy as np


class plotHandler(processHandler):

    def __init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, config, year, raw):

        processHandler.__init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, config, year, raw)
        self.initHistograms()



    #### ------------------------
    #### --
    #### ---- declareHistograms
    #### --
    #### ------------------------

    def initHistograms(self):
        
        bins_mass     = np.array([15., 65., 115., 165., 215., 265., 315., 365., 415., 465.])
        bins_trackIxy = np.array([1., 3., 6., 10., 15., 25., 35., 60.])
        #### ---------------
        #### ---- EE plots
        #### ---------------
        self.hEE_nEE            = {}
        self.hEE_nBSEE          = {}
        self.hEE_nPU            = {}
        self.hEE_MET            = {}
        self.hEE_MET_phi        = {}
        self.hEE_dR             = {}
        self.hEE_dPhi           = {}
        self.hEE_dPhi_scan      = {}
        self.hEE_dPhi_inv       = {}
        self.hEE_mass           = {}
        self.hEE_mass_Z         = {}
        self.hEE_mass_log       = {}
        self.hEE_mass_scan      = {}
        self.hEE_cosAlpha       = {}
        self.hEE_trackIxy       = {}
        self.hEE_trackIxy_log   = {}
        self.hEE_trackDxy       = {}
        self.hEE_sigmaD         = {}
        self.hEE_Lxy            = {}
        self.hEE_Ixy            = {}
        self.hEE_leadingEt      = {}
        self.hEE_subleadingEt   = {}
        self.hEE_eta            = {}
        self.hEE_normalizedChi2 = {}
        self.hEE_normalizedChi2_log = {}
        self.hEE_vr             = {}
        self.hEE_vr_low         = {}
        self.hEE_vr_high        = {}
        self.hEE_vx_vy          = {}
        self.hEE_phi_eta        = {}
        self.hEE_mass_trackIxy  = {}
        self.hEE_mass_Lxy       = {}
        self.hEE_Lxy_trackIxy   = {}
        self.hEE_Lxy_trackIxy_log   = {}
        self.hEE_trackDxy_trackIxy_log   = {}
        for reg in self.dielectronRegions:
            region = reg[0]
            self.hEE_nEE[region]            = r.TH1F('hEE{0}_nEE'.format(region) + self.sufix, ';Number of preselected dielectron vertices;', 6, 0, 6)
            self.hEE_nBSEE[region]          = r.TH1F('hEE{0}_nBSEE'.format(region) + self.sufix, ';Number of EE candidates;', 4, 0, 4)
            self.hEE_nPU[region]            = r.TH1F('hEE{0}_nPU'.format(region) + self.sufix, ';Number of true primary vertices;', 40, 0, 80)
            self.hEE_MET[region]            = r.TH1F('hEE{0}_MET'.format(region) + self.sufix, '; Missing E_{T} (GeV);', 50, 0, 200)
            self.hEE_MET_phi[region]        = r.TH1F('hEE{0}_MET_phi'.format(region) + self.sufix, '; Missing E_{T} #phi;', 50, 0, 2*3.14)
            self.hEE_dR[region]             = r.TH1F('hEE{0}_dR'.format(region) + self.sufix, ';Dielectron #Delta#R;', 30, 0, 3)
            self.hEE_dPhi[region]           = r.TH1F('hEE{0}_dPhi'.format(region) + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 30, 0, 3.14)
            self.hEE_dPhi_scan[region]      = r.TH1F('hEE{0}_dPhi_scan'.format(region) + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 4, 0, 3.14)
            self.hEE_dPhi_inv[region]       = r.TH1F('hEE{0}_dPhi_inv'.format(region) + self.sufix, ';Dielectron inverted collinearity #pi - |#Delta#Phi|;', 30, 0, 3.14)
            self.hEE_mass[region]           = r.TH1F('hEE{0}_mass'.format(region) + self.sufix, ';Dielectron invariant mass m_{ee} (GeV);', 50, 0, 300)
            self.hEE_mass_Z[region]         = r.TH1F('hEE{0}_mass_Z'.format(region) + self.sufix, ';Dielectron invariant mass m_{ee} (GeV);', 50, 75, 105)
            self.hEE_mass_log[region]       = r.TH1F('hEE{0}_mass_log'.format(region) + self.sufix, ';Dielectron invariant mass m_{ee} (GeV);', len(np.logspace(-1, 2, 100))-1, np.logspace(-1, 2, 100))
            self.hEE_mass_scan[region]      = r.TH1F('hEE{0}_mass_scan'.format(region) + self.sufix, ';Dielectron invariant mass m_{ee} (GeV);', 8, np.array([0., 21, 41., 61., 81., 101., 121., 141., 161.]))
            self.hEE_cosAlpha[region]       = r.TH1F('hEE{0}_cosAlpha'.format(region) + self.sufix, ';Dielectron cos(#alpha_{ee});', 44, -1.1, 1.1)
            self.hEE_trackIxy[region]       = r.TH1F('hEE{0}_trackIxy'.format(region) + self.sufix, ';Dielectron |d_{0}|/#sigma_{d};', 40, 0, 40)
            self.hEE_trackIxy_log[region]   = r.TH1F('hEE{0}_trackIxy_log'.format(region) + self.sufix, ';Dielectron |d_{0}|/#sigma_{d};', len(np.logspace(-3, 3, 50))-1, np.logspace(-3, 3, 50))
            self.hEE_trackDxy[region]       = r.TH1F('hEE{0}_trackDxy'.format(region) + self.sufix, ';Dielectron |d_{0}| (cm);', 60, 0, 60)
            self.hEE_sigmaD[region]         = r.TH1F('hEE{0}_sigmaD'.format(region) + self.sufix, ';Dielectron #sigma_{d} (cm);', 100, 0, 0.01)
            self.hEE_Lxy[region]            = r.TH1F('hEE{0}_Lxy'.format(region) + self.sufix, ';Dielectron vertex |L_{xy}| (cm);', 50, 0, 100)
            self.hEE_Ixy[region]            = r.TH1F('hEE{0}_Ixy'.format(region) + self.sufix, ';Dielectron vertex |L_{xy}|/#sigma_{L};', 40, 0, 40)
            self.hEE_leadingEt[region]      = r.TH1F('hEE{0}_leadingEt'.format(region) + self.sufix, ';Dielectron leading E_{T};', 40, 0, 200)
            self.hEE_subleadingEt[region]   = r.TH1F('hEE{0}_subleadingEt'.format(region) + self.sufix, ';Dielectron subleading E_{T};', 40, 0, 200)
            self.hEE_eta[region]            = r.TH1F('hEE{0}_eta'.format(region) + self.sufix, ';Electron #eta;', 24, -2.4, 2.4)
            self.hEE_normalizedChi2[region] = r.TH1F('hEE{0}_normalizedChi2'.format(region) + self.sufix, ';Dielectron vertex fit #chi^{2}/ndof;', 50, 0, 50)
            self.hEE_normalizedChi2_log[region] = r.TH1F('hEE{0}_normalizedChi2_log'.format(region) + self.sufix, ';Dielectron vertex fit #chi^{2}/ndof;', 60, np.logspace(-3, 4, 61))
            self.hEE_vr[region]             = r.TH1F('hEE{0}_vr'.format(region) + self.sufix, ';Dielectron vertex v_{r} (cm);', 150, 0, 25)
            self.hEE_vr_low[region]         = r.TH1F('hEE{0}_vr_low'.format(region) + self.sufix, ';Dielectron vertex v_{r} (cm);', 100, 0, 4)
            self.hEE_vr_high[region]        = r.TH1F('hEE{0}_vr_high'.format(region) + self.sufix, ';Dielectron vertex v_{r} (cm);', 100, 4, 8)
            self.hEE_vx_vy[region]          = r.TH2F('hEE{0}_vx_vy'.format(region) + self.sufix, ';Dielectron vertex v_{x} (cm) ; Dielectron vertex v_{y} (cm)', 150, -10.0, 10.0, 150, -10.0, 10.0)
            self.hEE_phi_eta[region]        = r.TH2F('hEE{0}_phi_eta'.format(region) + self.sufix, ';Electron #Phi; Electron #eta', 150, -3.2, 3.2, 150, -2.4, 2.4)
            self.hEE_mass_Lxy[region]       = r.TH2F('hEE{0}_mass_Lxy'.format(region) + self.sufix, ';Dielectron invariant mass m_{ee} (GeV); Dielectron vertex |L_{xy}| (cm)', 7, np.array([0., 15., 75., 105., 200., 300., 400., 500.]), 50, np.linspace(0, 5, 51))
            self.hEE_Lxy_trackIxy[region]   = r.TH2F('hEE{0}_Lxy_trackIxy'.format(region) + self.sufix, ';Dielectron vertex |L_{xy}| (cm);Dielectron |d_{0}|/#sigma_{d}', 40, np.linspace(0, 2, 41), 40, np.linspace(0, 40, 41))
            self.hEE_Lxy_trackIxy_log[region]   = r.TH2F('hEE{0}_Lxy_trackIxy_log'.format(region) + self.sufix, ';Dielectron vertex |L_{xy}| (cm);Dielectron |d_{0}|/#sigma_{d}', 80, np.logspace(-2, 2, 81), 80, np.logspace(0, 4, 81))
            self.hEE_trackDxy_trackIxy_log[region]   = r.TH2F('hEE{0}_trackDxy_trackIxy_log'.format(region) + self.sufix, ';Dielectron vertex |d_{0}| (cm);Dielectron |d_{0}|/#sigma_{d}', 80, np.logspace(-2, 2, 81), 80, np.logspace(0, 4, 81))
            self.hEE_mass_trackIxy[region]  = r.TH2F('hEE{0}_mass_trackIxy'.format(region) + self.sufix, ';Dielectron invariant mass (GeV) ; Dielectron |d_{0}|/#sigma_{d}', len(bins_mass)-1, bins_mass, len(bins_trackIxy)-1, bins_trackIxy)

        #### ---------------
        #### ---- MM plots
        #### ---------------
        self.hMM_nMM            = {}
        self.hMM_nBSMM          = {}
        self.hMM_nPU            = {}
        self.hMM_MET            = {}
        self.hMM_MET_phi        = {}
        self.hMM_ptDiff         = {}
        self.hMM_dR             = {}
        self.hMM_dPhi           = {}
        self.hMM_dPhi_inv       = {}
        self.hMM_dPhi_scan      = {}
        self.hMM_partdPhi1      = {}
        self.hMM_partdPhi2      = {}
        self.hMM_mass           = {}
        self.hMM_mass_Z         = {}
        self.hMM_mass_log       = {}
        self.hMM_mass_scan      = {}
        self.hMM_cosAlpha       = {}
        self.hMM_trackIxy       = {}
        self.hMM_trackIxy_log   = {}
        self.hMM_trackDxy       = {}
        self.hMM_sigmaD         = {}
        self.hMM_Lxy            = {}
        self.hMM_Ixy            = {}
        self.hMM_leadingPt      = {}
        self.hMM_subleadingPt   = {}
        self.hMM_eta            = {}
        self.hMM_normalizedChi2 = {}
        self.hMM_normalizedChi2_log = {}
        self.hMM_vr             = {}
        self.hMM_vr_low         = {}
        self.hMM_vr_high        = {}
        self.hMM_vx_vy          = {}
        self.hMM_phi_eta        = {}
        self.hMM_mass_trackIxy  = {}
        self.hMM_mass_Lxy       = {}
        self.hMM_Lxy_trackIxy   = {}
        self.hMM_Lxy_trackIxy_log   = {}
        self.hMM_trackDxy_trackIxy_log   = {}
        for reg in self.dimuonRegions:
            region = reg[0]
            self.hMM_nMM[region]            = r.TH1F('hMM{0}_nMM'.format(region) + self.sufix, ';Number of preselected dimuon vertices;', 6, 0, 6)
            self.hMM_nBSMM[region]          = r.TH1F('hMM{0}_nBSMM'.format(region) + self.sufix, ';Number of MM candidates;', 4, 0, 4)
            self.hMM_nPU[region]            = r.TH1F('hMM{0}_nPU'.format(region) + self.sufix, ';Number of true primary vertices;', 40, 0, 80)
            self.hMM_MET[region]            = r.TH1F('hMM{0}_MET'.format(region) + self.sufix, '; Missing E_{T} (GeV);', 50, 0, 200)
            self.hMM_MET_phi[region]        = r.TH1F('hMM{0}_MET_phi'.format(region) + self.sufix, '; Missing E_{T} #phi;', 50, 0, 2*3.14)
            self.hMM_ptDiff[region]         = r.TH1F('hMM{0}_ptDiff'.format(region) + self.sufix, ';Dimuon |#Deltap_{T}| GeV;', 30, 0, 30)
            self.hMM_dR[region]             = r.TH1F('hMM{0}_dR'.format(region) + self.sufix, ';Dimuon #Delta#R;', 30, 0, 3)
            self.hMM_dPhi[region]           = r.TH1F('hMM{0}_dPhi'.format(region) + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 30, 0, 3.14)
            self.hMM_dPhi_inv[region]       = r.TH1F('hMM{0}_dPhi_inv'.format(region) + self.sufix, ';Dimuon inverted collinearity #pi - |#Delta#Phi|;', 30, 0, 3.14)
            self.hMM_dPhi_scan[region]      = r.TH1F('hMM{0}_dPhi_scan'.format(region) + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 4, 0, 3.14)
            self.hMM_partdPhi1[region]      = r.TH1F('hMM{0}_partdPhi1'.format(region) + self.sufix, ';|#Delta#Phi(DV, #mu_{1})|;', 30, 0, 3.14)
            self.hMM_partdPhi2[region]      = r.TH1F('hMM{0}_partdPhi2'.format(region) + self.sufix, ';|#Delta#Phi(DV, #mu_{2})|;', 30, 0, 3.14)
            self.hMM_mass[region]           = r.TH1F('hMM{0}_mass'.format(region) + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', 50, 0, 300)
            self.hMM_mass_Z[region]         = r.TH1F('hMM{0}_mass_Z'.format(region) + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', 50, 80, 105)
            self.hMM_mass_log[region]       = r.TH1F('hMM{0}_mass_log'.format(region) + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', len(np.logspace(-1, 2, 100))-1, np.logspace(-1, 2, 100))
            self.hMM_mass_scan[region]      = r.TH1F('hMM{0}_mass_scan'.format(region) + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', 8, np.array([0., 21, 41., 61., 81., 101., 121., 141., 161.]))
            self.hMM_cosAlpha[region]       = r.TH1F('hMM{0}_cosAlpha'.format(region) + self.sufix, ';Dimuon cos(#alpha_{#mu#mu});', 44, -1.1, 1.1)
            self.hMM_trackIxy[region]       = r.TH1F('hMM{0}_trackIxy'.format(region) + self.sufix, ';Dimuon |d_{0}|/#sigma_{d};', 40, 0, 40)
            self.hMM_trackIxy_log[region]   = r.TH1F('hMM{0}_trackIxy_log'.format(region) + self.sufix, ';Dimuon |d_{0}|/#sigma_{d};', len(np.logspace(-3, 3, 50))-1, np.logspace(-3, 3, 50))
            self.hMM_trackDxy[region]       = r.TH1F('hMM{0}_trackDxy'.format(region) + self.sufix, ';Dimuon |d_{0}| (cm);', 60, 0, 60.0)
            self.hMM_sigmaD[region]         = r.TH1F('hMM{0}_sigmaD'.format(region) + self.sufix, ';Dimuon #sigma_{d} (cm);', 100, 0, 0.01)
            self.hMM_Lxy[region]            = r.TH1F('hMM{0}_Lxy'.format(region) + self.sufix, ';Dimuon vertex |L_{xy}| (cm);', 50, 0, 100)
            self.hMM_Ixy[region]            = r.TH1F('hMM{0}_Ixy'.format(region) + self.sufix, ';Dimuon vertex |L_{xy}|/#sigma_{L};', 40, 0, 40)
            self.hMM_leadingPt[region]      = r.TH1F('hMM{0}_leadingPt'.format(region) + self.sufix, ';Dimuon leading p_{T};', 40, 0, 200)
            self.hMM_subleadingPt[region]   = r.TH1F('hMM{0}_subleadingPt'.format(region) + self.sufix, ';Dimuon subleading p_{T};', 40, 0, 200)
            self.hMM_eta[region]            = r.TH1F('hMM{0}_eta'.format(region) + self.sufix, ';Muon #eta;', 24, -2.4, 2.4)
            self.hMM_normalizedChi2[region] = r.TH1F('hMM{0}_normalizedChi2'.format(region) + self.sufix, ';Dimuon vertex fit #chi^{2}/ndof;', 50, 0, 50)
            self.hMM_normalizedChi2_log[region] = r.TH1F('hMM{0}_normalizedChi2_log'.format(region) + self.sufix, ';Dimuon vertex fit #chi^{2}/ndof;', 60, np.logspace(-3, 4, 61))
            self.hMM_vr[region]             = r.TH1F('hMM{0}_vr'.format(region) + self.sufix, ';Dimuon vertex v_{r} (cm);', 150, 0, 25)
            self.hMM_vr_low[region]         = r.TH1F('hMM{0}_vr_low'.format(region) + self.sufix, ';Dimuon vertex v_{r} (cm);', 100, 0, 4)
            self.hMM_vr_high[region]        = r.TH1F('hMM{0}_vr_high'.format(region) + self.sufix, ';Dimuon vertex v_{r} (cm);', 100, 4, 8)
            self.hMM_vx_vy[region]          = r.TH2F('hMM{0}_vx_vy'.format(region) + self.sufix, ';Dimuon vertex v_{x} (cm) ; Dimuon vertex v_{y} (cm)', 150, -10.0, 10.0, 150, -10.0, 10.0)
            self.hMM_phi_eta[region]        = r.TH2F('hMM{0}_phi_eta'.format(region) + self.sufix, ';Muon #Phi; Muon #eta', 150, -3.2, 3.2, 150, -2.4, 2.4)
            self.hMM_mass_Lxy[region]       = r.TH2F('hMM{0}_mass_Lxy'.format(region) + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV); Dimuon vertex |L_{xy}| (cm)', 7, np.array([0., 15., 75., 105., 200., 300., 400., 500.]), 50, np.linspace(0, 5, 51))
            self.hMM_Lxy_trackIxy[region]   = r.TH2F('hMM{0}_Lxy_trackIxy'.format(region) + self.sufix, ';Dimuon vertex |L_{xy}| (cm);Dimuon |d_{0}|/#sigma_{d}', 40, np.linspace(0, 2, 41), 40, np.linspace(0, 40, 41))
            self.hMM_Lxy_trackIxy_log[region]   = r.TH2F('hMM{0}_Lxy_trackIxy_log'.format(region) + self.sufix, ';Dimuon vertex |L_{xy}| (cm);Dimuon |d_{0}|/#sigma_{d}', 80, np.logspace(-2, 2, 81), 80, np.logspace(0, 4, 81))
            self.hMM_trackDxy_trackIxy_log[region]   = r.TH2F('hMM{0}_trackDxy_trackIxy_log'.format(region) + self.sufix, ';Dimuon vertex |d_{0}| (cm);Dimuon |d_{0}|/#sigma_{d}', 80, np.logspace(-2, 2, 81), 80, np.logspace(0, 4, 81))
            self.hMM_mass_trackIxy[region]  = r.TH2F('hMM{0}_mass_trackIxy'.format(region) + self.sufix, ';Dimuon invariant mass (GeV) ; Dimuon |d_{0}|/#sigma_{d}', len(bins_mass)-1, bins_mass, len(bins_trackIxy)-1, bins_trackIxy)


        #### --------------------------------
        #### ---- Apply sumw2 to histograms
        #### --------------------------------

        for attr, value in self.__dict__.items():
            if attr[0] == 'h':
                if type(value) == dict:
                    for key in value.keys():
                        value[key].Sumw2()
                else:
                    value.Sumw2()

    #### ----------------------------
    #### --
    #### ---- Core event processing
    #### --   ( Overriden to fill histograms)
    #### ----------------------------

    def processEvent(self, ev):

        # Weight definition:
        weight = self.getWeight(ev, self.isdata)


        # Test PV acceptance:
        if not ev.PV_passAcceptance: return


        # Dimuon processing
        try:
            mm_maxIxy = -99
            mm_maxIxy, nBSMM = self.processDimuons(ev)
            passMuonTrigger = eval(self.mumu_path)
            if not mm_maxIxy < 0 and passMuonTrigger:
                imm = mm_maxIxy
                for region in self.dimuonRegions:
                    if eval(region[1]):
                        self.fillDimuonHistograms(ev, weight, region[0], mm_maxIxy, nBSMM)
        except AttributeError:
            print('There are some collections missed in this file: Dimuon histograms will be empty')


        # Dielectron processing
        try:
            ee_maxIxy = -99
            ee_maxIxy, nBSEE = self.processDielectrons(ev)
            passElectronTrigger = eval(self.ee_path)
            #passElectronTrigger = eval(self.ee_path) and not eval(self.mumu_path)
            if not ee_maxIxy < 0 and passElectronTrigger:
                iee = ee_maxIxy
                for region in self.dielectronRegions:
                    if eval(region[1]):
                        self.fillDielectronHistograms(ev, weight, region[0], ee_maxIxy, nBSEE)
        except AttributeError:
            print('There are some collections missed in this file: Dielectron histograms will be empty')



    #### ---------------------
    #### --
    #### ---- DiMuon filling
    #### --
    #### ---------------------

    def fillDimuonHistograms(self, ev, weight, region, mm_maxIxy, nBSMM):

        ## Compute SF
        sf = 1.0
        if not self.isdata:
            sf = self.getDimuonSF(ev, mm_maxIxy)

        self.hMM_nMM[region].Fill(ev.nDMDM, weight*sf)
        self.hMM_nPU[region].Fill(ev.nPV, weight*sf)
        self.hMM_MET[region].Fill(ev.MET_pt, weight*sf)
        self.hMM_MET_phi[region].Fill(ev.MET_phi, weight*sf)
        self.hMM_nBSMM[region].Fill(nBSMM, weight*sf)
        self.hMM_dR[region].Fill(abs(ev.DMDM_dR[mm_maxIxy]), weight*sf)
        self.hMM_dPhi[region].Fill(abs(ev.DMDM_dPhi[mm_maxIxy]), weight*sf)
        self.hMM_dPhi_scan[region].Fill(abs(ev.DMDM_dPhi[mm_maxIxy]), weight*sf)
        self.hMM_dPhi_inv[region].Fill(3.14 - abs(ev.DMDM_dPhi[mm_maxIxy]), weight*sf)
        self.hMM_mass[region].Fill(ev.DMDM_mass[mm_maxIxy], weight*sf)
        self.hMM_mass_Z[region].Fill(ev.DMDM_mass[mm_maxIxy], weight*sf)
        self.hMM_mass_log[region].Fill(ev.DMDM_mass[mm_maxIxy], weight*sf)
        self.hMM_mass_scan[region].Fill(ev.DMDM_mass[mm_maxIxy], weight*sf)
        self.hMM_trackIxy[region].Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight*sf)
        self.hMM_trackIxy_log[region].Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight*sf)
        self.hMM_trackDxy[region].Fill(abs(ev.DMDM_trackDxy_PV[mm_maxIxy]), weight*sf)
        self.hMM_Lxy[region].Fill(abs(ev.DMDM_Lxy_PV[mm_maxIxy]), weight*sf)
        self.hMM_Ixy[region].Fill(abs(ev.DMDM_Ixy_PV[mm_maxIxy]), weight*sf)
        self.hMM_normalizedChi2[region].Fill(ev.DMDM_normalizedChi2[mm_maxIxy], weight*sf)
        self.hMM_normalizedChi2_log[region].Fill(ev.DMDM_normalizedChi2[mm_maxIxy], weight*sf)
        self.hMM_cosAlpha[region].Fill(ev.DMDM_cosAlpha[mm_maxIxy], weight*sf)
        self.hMM_ptDiff[region].Fill(ev.DMDM_leadingPt[mm_maxIxy] - ev.DMDM_subleadingPt[mm_maxIxy], weight*sf)
        self.hMM_leadingPt[region].Fill(ev.DMDM_leadingPt[mm_maxIxy], weight*sf)
        self.hMM_subleadingPt[region].Fill(ev.DMDM_subleadingPt[mm_maxIxy], weight*sf)
        self.hMM_eta[region].Fill(ev.DGM_eta[ev.DMDM_idxA[mm_maxIxy]], weight*sf)
        self.hMM_eta[region].Fill(ev.DGM_eta[ev.DMDM_idxB[mm_maxIxy]], weight*sf)
        self.hMM_vr[region].Fill(np.sqrt(ev.DMDM_vx[mm_maxIxy]**2+ev.DMDM_vy[mm_maxIxy]**2), weight*sf)
        self.hMM_vr_low[region].Fill(np.sqrt(ev.DMDM_vx[mm_maxIxy]**2+ev.DMDM_vy[mm_maxIxy]**2), weight*sf)
        self.hMM_vr_high[region].Fill(np.sqrt(ev.DMDM_vx[mm_maxIxy]**2+ev.DMDM_vy[mm_maxIxy]**2), weight*sf)
        self.hMM_vx_vy[region].Fill(ev.DMDM_vx[mm_maxIxy], ev.DMDM_vy[mm_maxIxy], weight*sf)
        self.hMM_phi_eta[region].Fill(ev.DGM_phi[ev.DMDM_idxA[mm_maxIxy]], ev.DGM_eta[ev.DMDM_idxA[mm_maxIxy]], weight*sf)
        self.hMM_mass_trackIxy[region].Fill(ev.DMDM_mass[mm_maxIxy], ev.DMDM_trackIxy_PV[mm_maxIxy], weight*sf)
        self.hMM_mass_Lxy[region].Fill(ev.DMDM_mass[mm_maxIxy], abs(ev.DMDM_Lxy_PV[mm_maxIxy]), weight*sf)
        self.hMM_Lxy_trackIxy[region].Fill(abs(ev.DMDM_Lxy_PV[mm_maxIxy]), ev.DMDM_trackIxy_PV[mm_maxIxy], weight*sf)
        self.hMM_Lxy_trackIxy_log[region].Fill(abs(ev.DMDM_Lxy_PV[mm_maxIxy]), ev.DMDM_trackIxy_PV[mm_maxIxy], weight*sf)
        self.hMM_trackDxy_trackIxy_log[region].Fill(abs(ev.DMDM_trackDxy_PV[mm_maxIxy]), ev.DMDM_trackIxy_PV[mm_maxIxy], weight*sf)

        mu1 = r.TVector3()
        mu2 = r.TVector3()
        mmv = r.TVector3()
        mu1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[mm_maxIxy]], ev.DGM_eta[ev.DMDM_idxA[mm_maxIxy]], ev.DGM_phi[ev.DMDM_idxA[mm_maxIxy]])
        mu2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[mm_maxIxy]], ev.DGM_eta[ev.DMDM_idxB[mm_maxIxy]], ev.DGM_phi[ev.DMDM_idxB[mm_maxIxy]])
        mmv.SetXYZ(ev.DMDM_vx[mm_maxIxy], ev.DMDM_vy[mm_maxIxy], ev.DMDM_vz[mm_maxIxy])

        self.hMM_partdPhi1[region].Fill(mu1.DeltaPhi(mmv), weight*sf)
        self.hMM_partdPhi2[region].Fill(mu2.DeltaPhi(mmv), weight*sf)

    #### -------------------------
    #### --
    #### ---- DiElectron filling
    #### --
    #### -------------------------

    def fillDielectronHistograms(self, ev, weight, region, ee_maxIxy, nBSEE):

        ## Compute SF
        sf = 1.0
        if not self.isdata:
            sf = self.getDielectronSF(ev, ee_maxIxy)

        self.hEE_nEE[region].Fill(ev.nEE, weight*sf)
        self.hEE_nPU[region].Fill(ev.nPV, weight*sf)
        self.hEE_MET[region].Fill(ev.MET_pt, weight*sf)
        self.hEE_MET_phi[region].Fill(ev.MET_phi, weight*sf)
        self.hEE_nBSEE[region].Fill(nBSEE, weight*sf)
        self.hEE_dR[region].Fill(abs(ev.EE_dR[ee_maxIxy]), weight*sf)
        self.hEE_dPhi[region].Fill(abs(ev.EE_dPhi[ee_maxIxy]), weight*sf)
        self.hEE_dPhi_scan[region].Fill(abs(ev.EE_dPhi[ee_maxIxy]), weight*sf)
        self.hEE_dPhi_inv[region].Fill(3.14 - abs(ev.EE_dPhi[ee_maxIxy]), weight*sf)
        self.hEE_mass[region].Fill(ev.EE_mass[ee_maxIxy], weight*sf)
        self.hEE_mass_Z[region].Fill(ev.EE_mass[ee_maxIxy], weight*sf)
        self.hEE_mass_log[region].Fill(ev.EE_mass[ee_maxIxy], weight*sf)
        self.hEE_mass_scan[region].Fill(ev.EE_mass[ee_maxIxy], weight*sf)
        self.hEE_trackDxy[region].Fill(abs(ev.EE_trackDxy_PV[ee_maxIxy]), weight*sf)
        self.hEE_trackIxy[region].Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight*sf)
        self.hEE_trackIxy_log[region].Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight*sf)
        self.hEE_Lxy[region].Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), weight*sf)
        self.hEE_Ixy[region].Fill(abs(ev.EE_Ixy_PV[ee_maxIxy]), weight*sf)
        self.hEE_leadingEt[region].Fill(ev.EE_leadingEt[ee_maxIxy], weight*sf)
        self.hEE_subleadingEt[region].Fill(ev.EE_subleadingEt[ee_maxIxy], weight*sf)
        self.hEE_eta[region].Fill(ev.ElectronCandidate_eta[ev.EE_idxA[ee_maxIxy]], weight*sf)
        self.hEE_eta[region].Fill(ev.ElectronCandidate_eta[ev.EE_idxB[ee_maxIxy]], weight*sf)
        self.hEE_normalizedChi2[region].Fill(ev.EE_normalizedChi2[ee_maxIxy], weight*sf)
        self.hEE_normalizedChi2_log[region].Fill(ev.EE_normalizedChi2[ee_maxIxy], weight*sf)
        self.hEE_vr[region].Fill(np.sqrt(ev.EE_vx[ee_maxIxy]**2+ev.EE_vy[ee_maxIxy]**2), weight*sf)
        self.hEE_vr_low[region].Fill(np.sqrt(ev.EE_vx[ee_maxIxy]**2+ev.EE_vy[ee_maxIxy]**2), weight*sf)
        self.hEE_vr_high[region].Fill(np.sqrt(ev.EE_vx[ee_maxIxy]**2+ev.EE_vy[ee_maxIxy]**2), weight*sf)
        self.hEE_vx_vy[region].Fill(ev.EE_vx[ee_maxIxy], ev.EE_vy[ee_maxIxy], weight*sf)
        self.hEE_phi_eta[region].Fill(ev.ElectronCandidate_phi[ev.EE_idxA[ee_maxIxy]], ev.ElectronCandidate_eta[ev.EE_idxA[ee_maxIxy]], weight*sf)
        self.hEE_mass_trackIxy[region].Fill(ev.EE_mass[ee_maxIxy], ev.EE_trackIxy_PV[ee_maxIxy], weight*sf)
        self.hEE_mass_Lxy[region].Fill(ev.EE_mass[ee_maxIxy], abs(ev.EE_Lxy_PV[ee_maxIxy]), weight*sf)
        self.hEE_Lxy_trackIxy[region].Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), ev.EE_trackIxy_PV[ee_maxIxy], weight*sf)
        self.hEE_Lxy_trackIxy_log[region].Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), ev.EE_trackIxy_PV[ee_maxIxy], weight*sf)
        self.hEE_trackDxy_trackIxy_log[region].Fill(abs(ev.EE_trackDxy_PV[ee_maxIxy]), ev.EE_trackIxy_PV[ee_maxIxy], weight*sf)


