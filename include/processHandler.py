import math
import os
import json
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager
import numpy as np


class processHandler:

    def __init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, configfile, year = '2016', raw = False):

        ### outdir: Path where the root file will be stored
        ### treename: Name of the (Galapago) Tree (MC; DATA; SI)
        ### blockname: Name of the (Galapago) Block
        ### samplename: Name of the sample
        ### samplenumber: Number of the sample file/ttree

        self.outdir = outdir
        if self.outdir[-1] != '/': self.outdir = self.outdir + '/' 
        self.filename = self.outdir + '{0}__{1}__{2}__{3}.root'.format(treename, blockname, samplename, str(samplenumber))
        self.treename = treename
        self.blockname = blockname
        self.samplename = samplename
        self.samplenumber = samplenumber
        self.lumiweight = lumiweight
        self.isdata = isdata
        self.year = year
        self.cm = CutManager.CutManager()
        self.sufix = '__{0}__{1}__{2}__{3}'.format(treename, blockname, samplename, str(samplenumber))
        self.raw = raw

        ## Init configuration
        self.config = {}
        with open(configfile) as f:
            self.config = json.load(f)  


        ## Init histogram efficiencies
        self.declareEfficiencies()

        #### ------------------------------
        #### ---- Set year configuration
        #### ------------------------------
        self.mumu_path      = 'False'
        self.mumu_selection = 'False'
        self.ee_path        = 'False'
        self.ee_selection   = 'False'
        if year=='2016' or year=='2016APV':
            self.mumu_path      = self.cm.ORList(self.config['triggerPaths']['muons']['2016'], 'ev.')
            self.ee_path        = self.cm.ORList(self.config['triggerPaths']['electrons']['2016'], 'ev.')
            self.mumu_selection = self.cm.AddList(self.config['selection']['muons']['2016'])
            self.ee_selection   = self.cm.AddList(self.config['selection']['electrons']['2016'])
        elif year=='2017':
            self.ee_path        = self.cm.ORList(self.config['triggerPaths']['electrons']['2017'], 'ev.')
            self.ee_selection   = self.cm.AddList(self.config['selection']['electrons']['2017'])
        elif year=='2018':
            self.mumu_path      = self.cm.ORList(self.config['triggerPaths']['muons']['2018'], 'ev.')
            self.ee_path        = self.cm.ORList(self.config['triggerPaths']['electrons']['2018'], 'ev.')
            self.mumu_selection = self.cm.AddList(self.config['selection']['muons']['2018'])
            self.ee_selection   = self.cm.AddList(self.config['selection']['electrons']['2018'])

        #### ------------------------------
        #### ---- Regions initialization
        #### ------------------------------
    
        self.dielectronRegions = [] # region name : region cut
        for region in self.config["regions"]["electrons"].keys():
           self.dielectronRegions.append([region, self.cm.AddList(self.config["regions"]["electrons"][region])])
        
        self.dimuonRegions = [] # region name : region cut
        for region in self.config["regions"]["muons"].keys():
           self.dimuonRegions.append([region, self.cm.AddList(self.config["regions"]["muons"][region])])

        #### -------------------------------
        #### ---- Histogram initialization
        #### -------------------------------

        self.hEE_yield = r.TH1F('hEE_yield' + self.sufix, ';Passed photon trigger;', 1, 0, 1)
        self.hMM_yield = r.TH1F('hMM_yield' + self.sufix, ';Passed muon trigger;', 1, 0, 1)


        #### --------------------------------
        #### ---- Load SFs
        #### --------------------------------

        ### Reco + ID SFs
        if '2016' in year:
            self.file_sf_ee_reco = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/Electron_ScaleFactors_2016_Fall22.root")
            self.sf_ee_reco = self.file_sf_ee_reco.Get("NUM_genTracksDown_DEN_tagsInTime_absdxy_2d_absdz_2d")
            self.file_sf_mm_reco = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/Muon_ScaleFactors_2016_Fall22.root")
            self.sf_mm_reco = self.file_sf_mm_reco.Get("NUM_dGlobalsUp_DEN_dGlobalsDown_absdxy_2d_absdz_2d")
            self.file_sf_mm_id = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/MuonID_ScaleFactors_2016_Fall22.root")
            self.sf_mm_id = self.file_sf_mm_id.Get("NUM_dGlobalID_DEN_dGlobalsUp_absdxy_2d_absdz_2d")
        elif '2017' in year:
            self.file_sf_ee_reco = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/Electron_ScaleFactors_2017_Fall22.root")
            self.sf_ee_reco = self.file_sf_ee_reco.Get("NUM_genTracksDown_DEN_tagsInTime_absdxy_2d_absdz_2d")
        elif '2018' in year:
            self.file_sf_ee_reco = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/Electron_ScaleFactors_2018_Fall22.root")
            self.sf_ee_reco = self.file_sf_ee_reco.Get("NUM_genTracksDown_DEN_tagsInTime_absdxy_2d_absdz_2d")
            self.file_sf_mm_reco = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/Muon_ScaleFactors_2018_Fall22.root")
            self.sf_mm_reco = self.file_sf_mm_reco.Get("NUM_dGlobalsUp_DEN_dGlobalsDown_absdxy_2d_absdz_2d")
            self.file_sf_mm_id = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/MuonID_ScaleFactors_2018_Fall22.root")
            self.sf_mm_id = self.file_sf_mm_id.Get("NUM_dGlobalID_DEN_dGlobalsUp_absdxy_2d_absdz_2d")

        ### Trigger
        self.file_sf_ee_trg = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/PhotonTrigger_ScaleFactors_"+year+".root")
        self.sf_ee_trg = self.file_sf_ee_trg.Get("ScaleFactor_pt2_pt")
        if year != '2017':
            self.file_sf_mm_trg = r.TFile("/eos/user/f/fernance/LLP_Analysis/calibration/MuonTrigger_ScaleFactors_"+year+".root")
            self.sf_mm_trg = self.file_sf_mm_trg.Get("ScaleFactor_pt2_pt")


    #### -------------------------------------------
    #### --
    #### ---- Get event weights and Scale Factors 
    #### --   
    #### -------------------------------------------

    def getWeight(self, ev, isData):

        if not self.isdata:
            weight = self.lumiweight*ev.wPU*ev.genWeight/abs(ev.genWeight)
        else: 
            weight = 1.

        return weight


    def getDimuonSF(self, ev, idx):

        sf = 1.0 # default

        if self.year == '2017' or self.raw:
            return sf

        ### Reco 
        ## The corrections to account for early muons are applied by hand
        ## (Numbers were extracted when measuring scale factors)
        corr = {}
        corr['2016'] = 0.976 ## Derived
        corr['2016APV'] = 0.976 ## Derived
        corr['2018'] = 1.097 ## Derived
        bx1 = self.sf_mm_reco.GetXaxis().FindBin(abs(ev.DGM_dxy_PV[ev.DMDM_idxA[idx]]))
        bx2 = self.sf_mm_reco.GetXaxis().FindBin(abs(ev.DGM_dxy_PV[ev.DMDM_idxB[idx]]))
        by1 = self.sf_mm_reco.GetYaxis().FindBin(abs(ev.DGM_dz[ev.DMDM_idxA[idx]]))
        by2 = self.sf_mm_reco.GetYaxis().FindBin(abs(ev.DGM_dz[ev.DMDM_idxB[idx]]))
        sf = sf * self.sf_mm_reco.GetBinContent(bx1, by1)*self.sf_mm_reco.GetBinContent(bx2, by2) * (corr[self.year])**2

        ### ID
        bx1 = self.sf_mm_id.GetXaxis().FindBin(abs(ev.DGM_dxy_PV[ev.DMDM_idxA[idx]]))
        bx2 = self.sf_mm_id.GetXaxis().FindBin(abs(ev.DGM_dxy_PV[ev.DMDM_idxB[idx]]))
        by1 = self.sf_mm_id.GetYaxis().FindBin(abs(ev.DGM_dz[ev.DMDM_idxA[idx]]))
        by2 = self.sf_mm_id.GetYaxis().FindBin(abs(ev.DGM_dz[ev.DMDM_idxB[idx]]))
        sf = sf * self.sf_mm_id.GetBinContent(bx1, by1)*self.sf_mm_id.GetBinContent(bx2, by2)


        ### Trigger
        bx = self.sf_mm_trg.GetXaxis().FindBin(ev.DMDM_subleadingPt[idx])
        by = self.sf_mm_trg.GetYaxis().FindBin(ev.DMDM_leadingPt[idx])
        sf = sf * self.sf_mm_trg.GetBinContent(bx, by)

        return sf


    def getDielectronSF(self, ev, idx):

        sf = 1.0 # default

        if self.raw:
            return 1.0

        ### Reco + ID
        bx1 = self.sf_ee_reco.GetXaxis().FindBin(abs(ev.ElectronCandidate_dxy_PV[ev.EE_idxA[idx]]))
        bx2 = self.sf_ee_reco.GetXaxis().FindBin(abs(ev.ElectronCandidate_dxy_PV[ev.EE_idxB[idx]]))
        by1 = self.sf_ee_reco.GetYaxis().FindBin(abs(ev.IsoTrackSel_dz[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[idx]]]))
        by2 = self.sf_ee_reco.GetYaxis().FindBin(abs(ev.IsoTrackSel_dz[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[idx]]]))
        sf = sf * self.sf_ee_reco.GetBinContent(bx1, by1)*self.sf_ee_reco.GetBinContent(bx2, by2)

        ### Trigger
        bx = self.sf_ee_trg.GetXaxis().FindBin(ev.EE_subleadingEt[idx])
        by = self.sf_ee_trg.GetYaxis().FindBin(ev.EE_leadingEt[idx])
        sf = sf * self.sf_ee_trg.GetBinContent(bx, by)

        return sf

    #### ----------------------------
    #### --
    #### ---- Core event processing
    #### --   ( Overriden in child classes)
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

            if not mm_maxIxy < 0:

                imm = mm_maxIxy
                self.hMM_yield.Fill(0, weight)

        except AttributeError:
            print('There are some collections missed in this file: Dimuon histograms will be empty')
        
        
        # Dielectron processing
        try:

            ee_maxIxy = -99
            ee_maxIxy, nBSEE = self.processDielectrons(ev)

            if not ee_maxIxy < 0:

                iee = ee_maxIxy
                self.hEE_yield.Fill(0, weight)

        except AttributeError:
            print('There are some collections missed in this file: Dielectron histograms will be empty')
        



       
    #### ------------------------
    #### --
    #### ---- DiMuon Processing
    #### --
    #### ------------------------

    def processDimuons(self, ev):

        if not eval(self.mumu_path):
            return -99, -99
        if ev.nDMDM < 1:
            return -99, -99

        #### ----------------------------
        #### ---- Select a MM candidate
        #### ----------------------------

        ### -> Find the MM pair with maximum Ixy that passes the cuts
        mm_maxIxy = -99
        maxIxy = -1
        nBSMM = 0

        for i in range(0, ev.nDMDM):

            imm = i
            if not eval(self.mumu_selection): continue

            if eval(self.cm.MM_iso2l) and eval(self.cm.MM_OS): nBSMM+=1

            if ev.DMDM_trackIxy_PV[i] > maxIxy:
                maxIxy = ev.DMDM_trackIxy_PV[i]
                mm_maxIxy = i

        # return the index of the selected MM and the number of BS MM's:
        return mm_maxIxy, nBSMM 


    #### ----------------------------
    #### --
    #### ---- DiElectron Processing
    #### --
    #### ----------------------------

    def processDielectrons(self, ev):


        if not eval(self.ee_path) or eval(self.mumu_path):
            return -99, -99
        if ev.nEE < 1:
            return -99, -99

        ### -> Find the EE pair with maximum Ixy that passes the cuts
        ee_maxIxy = -99
        maxIxy = -1
        nBSEE = 0
        for i in range(0, ev.nEE):
         
            iee = i
            if not eval(self.ee_selection): continue
           
            if eval(self.cm.EE_iso2l) and eval(self.cm.EE_OS): nBSEE+=1

            if ev.EE_trackIxy_PV[i] > maxIxy:
                maxIxy = ev.EE_trackIxy_PV[i]
                ee_maxIxy = i

        # return the index of the selected MM and the number of BS MM's:
        return ee_maxIxy, nBSEE 



    #### -----------------------
    #### --
    #### ---- fillEfficiencies
    #### --
    #### -----------------------
    def countLLs(self, ev):

        if not self.isdata:
            weight = self.lumiweight*ev.wPU*ev.genWeight/abs(ev.genWeight)
            weightnoPU = self.lumiweight*ev.genWeight/abs(ev.genWeight)
        else:
            weight = 1


        ### EE Channel
        if eval(self.ee_path) and not eval(self.mumu_path):

            self.hEE_cutEfficiency.Fill(0, weight) 

            if not ev.PV_passAcceptance: return
            self.hEE_cutEfficiency.Fill(1, weight)

            if ev.nEE < 1: return
            self.hEE_cutEfficiency.Fill(2, weight)

            pass_etmin = False
            pass_etamax = False
            pass_normChi2max = False
            pass_relisomax = False
            pass_massmin = False
            pass_oscharge = False
            pass_offZ = False
            pass_dPhiSR = False
            pass_Ixymin = False

            for i in range(0, ev.nEE):

                iee = i

                if not eval(self.cm.EE_et1_40): continue
                if not eval(self.cm.EE_et2_25): continue
                pass_etmin = True

                if not eval(self.cm.EE_eta2): continue
                if not eval(self.cm.EE_etanoBE): continue
                pass_etamax = True

                if not eval(self.cm.EE_mass15): continue
                pass_massmin = True

                if not eval(self.cm.EE_iso2l): continue
                pass_relisomax = True

                if not eval(self.cm.EE_normChi2_10): continue
                pass_normChi2max = True 

                if not eval(self.cm.EE_OS): continue
                pass_oscharge = True

                if not eval(self.cm.EE_dPhiforward): continue
                pass_dPhiSR = True

                if not eval(self.cm.EE_OffZ): continue
                pass_offZ = True

                if not eval(self.cm.EE_Ixy6high): continue
                pass_Ixymin = True

            if pass_etmin: self.hEE_cutEfficiency.Fill(3, weight)
            if pass_etamax: self.hEE_cutEfficiency.Fill(4, weight)
            if pass_massmin: self.hEE_cutEfficiency.Fill(5, weight)
            if pass_relisomax: self.hEE_cutEfficiency.Fill(6, weight)
            if pass_normChi2max: self.hEE_cutEfficiency.Fill(7, weight)
            if pass_oscharge: self.hEE_cutEfficiency.Fill(8, weight)
            if pass_dPhiSR: self.hEE_cutEfficiency.Fill(9, weight)
            if pass_offZ: self.hEE_cutEfficiency.Fill(10, weight)
            if pass_Ixymin: self.hEE_cutEfficiency.Fill(11, weight)


        ### MM Channel
        if eval(self.mumu_path):

            self.hMM_cutEfficiency.Fill(0, weight) 

            if not ev.PV_passAcceptance: return
            self.hMM_cutEfficiency.Fill(1, weight)

            if ev.nDMDM < 1: return

            pass_ID = False
            pass_ptmin = False
            pass_etamax = False
            pass_normChi2max = False
            pass_dRmin = False
            pass_relisomax = False
            pass_cosAlphamin = False
            pass_massmin = False
            pass_oscharge = False
            pass_dPhiSR = False
            pass_offZ = False
            pass_Ixymin = False

            for i in range(0, ev.nDMDM):

                ## Index setting
                imm = i

                if not eval(self.cm.MM_ID): continue
                pass_ID = True

                if self.year == '2016' and not eval(self.cm.MM_pt30): continue
                if self.year == '2018' and not eval(self.cm.MM_pt25): continue
                pass_ptmin = True

                if not eval(self.cm.MM_eta2): continue
                pass_etamax = True

                if not eval(self.cm.MM_dR0p2): continue
                pass_dRmin = True

                if not eval(self.cm.MM_mass15): continue
                pass_massmin = True

                if not eval(self.cm.MM_iso2l): continue
                pass_relisomax = True

                if not eval(self.cm.MM_normChi2_10): continue
                pass_normChi2max = True 

                if self.year == '2016' and not eval(self.cm.MM_cosAlpha0p8): continue
                if self.year == '2018' and not eval(self.cm.MM_cosAlpha0p9): continue
                pass_cosAlphamin = True

                if not eval(self.cm.MM_OS): continue
                pass_oscharge = True

                if not eval(self.cm.MM_dPhiforward): continue
                pass_dPhiSR = True

                if not eval(self.cm.MM_OffZ): continue
                pass_offZ = True

                if not eval(self.cm.MM_Ixy6high): continue
                pass_Ixymin = True

            if pass_ID: self.hMM_cutEfficiency.Fill(2, weight)
            if pass_ptmin: self.hMM_cutEfficiency.Fill(3, weight)
            if pass_etamax: self.hMM_cutEfficiency.Fill(4, weight)
            if pass_dRmin: self.hMM_cutEfficiency.Fill(5, weight)
            if pass_massmin: self.hMM_cutEfficiency.Fill(6, weight)
            if pass_relisomax: self.hMM_cutEfficiency.Fill(7, weight)
            if pass_normChi2max: self.hMM_cutEfficiency.Fill(8, weight)
            if pass_cosAlphamin: self.hMM_cutEfficiency.Fill(9, weight)
            if pass_oscharge: self.hMM_cutEfficiency.Fill(10, weight)
            if pass_dPhiSR: self.hMM_cutEfficiency.Fill(11, weight)
            if pass_offZ: self.hMM_cutEfficiency.Fill(12, weight)
            if pass_Ixymin: self.hMM_cutEfficiency.Fill(13, weight)



    #### --------------------------
    #### --
    #### ---- declare Efficiencies
    #### --
    #### --------------------------

    def declareEfficiencies(self):

        self.hEE_cutEfficiency = r.TH1F('hEE_cutEfficiency' + self.sufix, ';; Events with EE / Category', 12, 0, 12)
        self.hMM_cutEfficiency = r.TH1F('hMM_cutEfficiency' + self.sufix, ';; Events with MM / Category', 14, 0, 14)



    #### --------------------------
    #### --
    #### ---- Write output
    #### --
    #### --------------------------


    def Write(self):
    
        while True:
            d = os.path.dirname(self.filename)
            try:
                if not os.path.exists(d): os.makedirs(d)
                break
            except OSError as e:
                if e.errno != os.errno.EEXIST:
                    raise
                print("Sleeping...")
                time.sleep(1.0)
                pass

        output = r.TFile(self.filename, 'RECREATE')

        if self.hEE_cutEfficiency.GetEntries() > 0:
            self.hEE_cutEfficiency.Write()
            self.hMM_cutEfficiency.Write()
        else:
            for attr, value in self.__dict__.items():
                if attr[0] == 'h':
                    if type(value) == dict:
                        for key in value.keys():
                            value[key].Write()
                    else:
                        value.Write()

        
        output.Close()

