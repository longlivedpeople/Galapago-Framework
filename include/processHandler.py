import math
import os
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager
import numpy as np


class processHandler:

    def __init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, year = '2016'):

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
        self.cm = CutManager.CutManager()
        self.sufix = '__{0}__{1}__{2}__{3}'.format(treename, blockname, samplename, str(samplenumber))



        ## Init histogram efficiencies
        self.declareEfficiencies()

        #### ------------------------------
        #### ---- Set year configuration
        #### ------------------------------
        self.mumu_path      = 'False'
        self.mumu_selection = 'False'
        self.ee_path        = 'False'
        self.ee_selection   = 'False'
        if year=='2016':
            self.mumu_path      = self.cm.mupath2016
            self.ee_path        = self.cm.epath2016
            self.mumu_selection = self.cm.MM_BS2016
            self.ee_selection   = self.cm.EE_BS2016
        elif year=='2017':
            self.ee_path        = self.cm.epath2017
            self.ee_selection   = self.cm.EE_BS2017
        elif year=='2018':
            self.mumu_path      = self.cm.mupath2018
            self.ee_path        = self.cm.epath2018
            self.mumu_selection = self.cm.MM_BS2018
            self.ee_selection   = self.cm.EE_BS2018

        #### ------------------------------
        #### ---- Regions initialization
        #### ------------------------------
    
        self.dielectronRegions = [] # region name : region cut
        self.declareDielectronRegions()
        
        self.dimuonRegions = [] # region name : region cut
        self.declareDimuonRegions()


        #### -------------------------------
        #### ---- Histogram initialization
        #### -------------------------------
        self.hEE_yield = r.TH1F('hEE_yield' + self.sufix, ';Passed photon trigger;', 1, 0, 1)
        self.hEE_nBSEE          = {}
        self.hEE_nPU            = {}
        self.hEE_dPhi           = {}
        self.hEE_dPhi_inv       = {}
        self.hEE_mass           = {}
        self.hEE_cosAlpha       = {}
        self.hEE_trackIxy       = {}
        self.hEE_trackIxy_log   = {}
        self.hEE_trackDxy       = {}
        self.hEE_Lxy            = {}
        self.hEE_Ixy            = {}
        self.hEE_leadingEt      = {}
        self.hEE_subleadingEt   = {}
        self.hEE_normalizedChi2 = {}
        self.hEE_vx_vy          = {}
        for region in self.dielectronRegions:
            region_name = region[0]
            self.declareDielectronHistograms(region_name)

        self.hMM_yield = r.TH1F('hMM_yield' + self.sufix, ';Passed muon trigger;', 1, 0, 1)
        self.hMM_nBSMM          = {}
        self.hMM_nPU            = {}
        self.hMM_dPhi           = {}
        self.hMM_dPhi_inv       = {}
        self.hMM_mass           = {}
        self.hMM_cosAlpha       = {}
        self.hMM_trackIxy       = {}
        self.hMM_trackIxy_log   = {}
        self.hMM_trackDxy       = {}
        self.hMM_Lxy            = {}
        self.hMM_Ixy            = {}
        self.hMM_leadingPt      = {}
        self.hMM_subleadingPt   = {}
        self.hMM_normalizedChi2 = {}
        self.hMM_vx_vy          = {}
        for region in self.dimuonRegions:
            region_name = region[0]
            self.declareDimuonHistograms(region_name)



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



    #### ------------------------
    #### --
    #### ---- declareHistograms
    #### --
    #### ------------------------

    def declareDimuonHistograms(self, region):

        binlog_Ixy = np.logspace(-3, 3, 50)

        #### ---------------
        #### ---- MM plots
        #### ---------------

        self.hMM_nBSMM[region]          = r.TH1F('hMM{0}_nBSMM'.format(region) + self.sufix, ';Number of MM candidates;', 4, 0, 4)
        self.hMM_nPU[region]            = r.TH1F('hMM{0}_nPU'.format(region) + self.sufix, ';Number of true primary vertices;', 40, 0, 80)
        self.hMM_dPhi[region]           = r.TH1F('hMM{0}_dPhi'.format(region) + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 30, 0, 3.14)
        self.hMM_dPhi_inv[region]       = r.TH1F('hMM{0}_dPhi_inv'.format(region) + self.sufix, ';Dimuon inverted collinearity #pi - |#Delta#Phi|;', 30, 0, 3.14)
        self.hMM_mass[region]           = r.TH1F('hMM{0}_mass'.format(region) + self.sufix, ';Dimuon invariant mass m_{{#mu#mu}} (GeV);', 80, 0, 400)
        self.hMM_cosAlpha[region]       = r.TH1F('hMM{0}_cosAlpha'.format(region) + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}});', 22, -1.1, 1.1)
        self.hMM_trackIxy[region]       = r.TH1F('hMM{0}_trackIxy'.format(region) + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', 40, 0, 40)
        self.hMM_trackIxy_log[region]   = r.TH1F('hMM{0}_trackIxy_log'.format(region) + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', len(binlog_Ixy)-1, binlog_Ixy)
        self.hMM_trackDxy[region]       = r.TH1F('hMM{0}_trackDxy'.format(region) + self.sufix, ';Dimuon |d_{{0}}| (cm);', 30, 0, 0.5)
        self.hMM_Lxy[region]            = r.TH1F('hMM{0}_Lxy'.format(region) + self.sufix, ';Dimuon vertex |L_{{xy}}| (cm);', 40, 0, 10) 
        self.hMM_Ixy[region]            = r.TH1F('hMM{0}_Ixy'.format(region) + self.sufix, ';Dimuon vertex |L_{{xy}}|/#sigma_{{L}};', 40, 0, 40)
        self.hMM_leadingPt[region]      = r.TH1F('hMM{0}_leadingPt'.format(region) + self.sufix, ';Dimuon leading p_{{T}};', 30, 0, 300)
        self.hMM_subleadingPt[region]   = r.TH1F('hMM{0}_subleadingPt'.format(region) + self.sufix, ';Dimuon subleading p_{{T}};', 30, 0, 300)
        self.hMM_normalizedChi2[region] = r.TH1F('hMM{0}_normalizedChi2'.format(region) + self.sufix, ';Dimuon vertex fit #chi^{{2}}/ndof;', 50, 0, 50)      
        self.hMM_vx_vy[region]          = r.TH2F('hMM{0}_vx_vy'.format(region) + self.sufix, ';Dimuon vertex v_{{x}} (cm) ; Dimuon vertex v_{{y}} (cm)', 200, -4.0, 4.0, 200, -4.0, 4.0)



    def declareDielectronHistograms(self, region):

        binlog_Ixy = np.logspace(-3, 3, 50)

        #### ---------------
        #### ---- EE plots
        #### ---------------

        self.hEE_nBSEE[region]          = r.TH1F('hEE{0}_nBSEE'.format(region) + self.sufix, ';Number of EE candidates;', 4, 0, 4)
        self.hEE_nPU[region]            = r.TH1F('hEE{0}_nPU'.format(region) + self.sufix, ';Number of true primary vertices;', 40, 0, 80)
        self.hEE_dPhi[region]           = r.TH1F('hEE{0}_dPhi'.format(region) + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 30, 0, 3.14)
        self.hEE_dPhi_inv[region]       = r.TH1F('hEE{0}_dPhi_inv'.format(region) + self.sufix, ';Dielectron inverted collinearity #pi - |#Delta#Phi|;', 30, 0, 3.14)
        self.hEE_mass[region]           = r.TH1F('hEE{0}_mass'.format(region) + self.sufix, ';Dielectron invariant mass m_{{ee}} (GeV);', 80, 0, 400) 
        self.hEE_cosAlpha[region]       = r.TH1F('hEE{0}_cosAlpha'.format(region) + self.sufix, ';Dielectron cos(#alpha_{{ee}});', 22, -1.1, 1.1)
        self.hEE_trackIxy[region]       = r.TH1F('hEE{0}_trackIxy'.format(region) + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', 40, 0, 40)
        self.hEE_trackIxy_log[region]   = r.TH1F('hEE{0}_trackIxy_log'.format(region) + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', len(binlog_Ixy)-1, binlog_Ixy)
        self.hEE_trackDxy[region]       = r.TH1F('hEE{0}_trackDxy'.format(region) + self.sufix, ';Dielectron |d_{{0}}| (cm);', 30, 0, 0.5)
        self.hEE_Lxy[region]            = r.TH1F('hEE{0}_Lxy'.format(region) + self.sufix, ';Dielectron vertex |L_{{xy}}| (cm);', 40, 0, 10)
        self.hEE_Ixy[region]            = r.TH1F('hEE{0}_Ixy'.format(region) + self.sufix, ';Dielectron vertex |L_{{xy}}|/#sigma_{{L}};', 40, 0, 40) 
        self.hEE_leadingEt[region]      = r.TH1F('hEE{0}_leadingEt'.format(region) + self.sufix, ';Dielectron leading E_{{T}};', 30, 0, 300)
        self.hEE_subleadingEt[region]   = r.TH1F('hEE{0}_subleadingEt'.format(region) + self.sufix, ';Dielectron subleading E_{{T}};', 30, 0, 300)
        self.hEE_normalizedChi2[region] = r.TH1F('hEE{0}_normalizedChi2'.format(region) + self.sufix, ';Dielectron vertex fit #chi^{{2}}/ndof;', 50, 0, 50)
        self.hEE_vx_vy[region]          = r.TH2F('hEE{0}_vx_vy'.format(region) + self.sufix, ';Dielectron vertex v_{{x}} (cm) ; Dielectron vertex v_{{y}} (cm)', 200, -4.0, 4.0, 200, -4.0, 4.0)




    #### ------------------------
    #### --
    #### ---- declareRegions
    #### --
    #### ------------------------

    def declareDielectronRegions(self):

        #### ---------------------------
        #### ---- EE Region definition
        #### ---------------------------

        EE_onZCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OnZ, self.cm.EE_dPhibackward])
        EE_onZSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OnZ, self.cm.EE_dPhiforward])
        EE_promptSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_Ixy6prompt, self.cm.EE_dPhiforward])
        EE_promptCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_Ixy6prompt, self.cm.EE_dPhibackward])


        #self.dielectronRegions.append(['sel', self.cm.EE_sel])
        #self.dielectronRegions.append(['Wjets', self.cm.EE_Wjets])
        #self.dielectronRegions.append(['QCD', self.cm.EE_QCD])
        #self.dielectronRegions.append(['promptCR', EE_promptCRcut])
        #self.dielectronRegions.append(['promptSR', EE_promptSRcut])
        #self.dielectronRegions.append(['onZCR', EE_onZCRcut])
        #self.dielectronRegions.append(['onZSR', EE_onZSRcut])
        #self.dielectronRegions.append(['SR', self.cm.EE_SR])
        self.dielectronRegions.append(['SRI', self.cm.EE_SRI])
        self.dielectronRegions.append(['SRII', self.cm.EE_SRII])
        #self.dielectronRegions.append(['BCR', self.cm.EE_BCR])
        self.dielectronRegions.append(['BCRI', self.cm.EE_BCRI])
        self.dielectronRegions.append(['BCRII', self.cm.EE_BCRII])


    def declareDimuonRegions(self):


        #### ---------------------------
        #### ---- MM Region definition
        #### ---------------------------

        MM_onZCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OnZ, self.cm.MM_dPhibackward])
        MM_onZSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OnZ, self.cm.MM_dPhiforward])
        MM_promptSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_Ixy6prompt, self.cm.MM_dPhiforward])
        MM_promptCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_Ixy6prompt, self.cm.MM_dPhibackward])

        #self.dimuonRegions.append(['sel', self.cm.MM_sel])
        #self.dimuonRegions.append(['Wjets', self.cm.MM_Wjets])
        #self.dimuonRegions.append(['QCD', self.cm.MM_QCD])
        #self.dimuonRegions.append(['promptCR', MM_promptCRcut])
        #self.dimuonRegions.append(['promptSR', MM_promptSRcut])
        #self.dimuonRegions.append(['onZCR', MM_onZCRcut])
        #self.dimuonRegions.append(['onZSR', MM_onZSRcut])
        #self.dimuonRegions.append(['SR', self.cm.MM_SR])
        self.dimuonRegions.append(['SRI', self.cm.MM_SRI])
        self.dimuonRegions.append(['SRII', self.cm.MM_SRII])
        #self.dimuonRegions.append(['BCR', self.cm.MM_BCR])
        self.dimuonRegions.append(['BCRI', self.cm.MM_BCRI])
        self.dimuonRegions.append(['BCRII', self.cm.MM_BCRII])




    #### ----------------------------
    #### --
    #### ---- Core event processing
    #### --
    #### ----------------------------

    def processEvent(self, ev):


        # Weight definition:
        if not self.isdata:
            weight = self.lumiweight*ev.wPU*ev.genWeight/abs(ev.genWeight)
            weightnoPU = self.lumiweight*ev.genWeight/abs(ev.genWeight)
        else: 
            weight = 1


        # Test PV acceptance:
        if not ev.PV_passAcceptance: return

        # Count
        if eval(self.mumu_path):
            self.hMM_yield.Fill(0, weight)
        if eval(self.ee_path) and not eval(self.mumu_path):
            self.hEE_yield.Fill(0, weight)

        # Dimuon processing
        
        try:

            mm_maxIxy = -99
            mm_maxIxy, nBSMM = self.processDimuons(ev)

            if not mm_maxIxy < 0:

                imm = mm_maxIxy
                for region in self.dimuonRegions:
                    if eval(region[1]):
                        self.fillDimuons(ev, weight, region[0], mm_maxIxy, nBSMM)

        except AttributeError:
            print('There are some collections missed in this file: Dimuon histograms will be empty')
        
        
        # Dielectron processing
        try:

            ee_maxIxy = -99
            ee_maxIxy, nBSEE = self.processDielectrons(ev)

            if not ee_maxIxy < 0:

                iee = ee_maxIxy
                for region in self.dielectronRegions:
                    if eval(region[1]):
                        self.fillDielectrons(ev, weight, region[0], ee_maxIxy, nBSEE)

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



    #### ---------------------
    #### --
    #### ---- DiMuon filling
    #### --
    #### ---------------------

    def fillDimuons(self, ev, weight, region, mm_maxIxy, nBSMM):

        self.hMM_nPU[region].Fill(ev.nPV, weight)
        self.hMM_nBSMM[region].Fill(nBSMM, weight) 
        self.hMM_dPhi[region].Fill(abs(ev.DMDM_dPhi[mm_maxIxy]), weight)
        self.hMM_dPhi_inv[region].Fill(3.14 - abs(ev.DMDM_dPhi[mm_maxIxy]), weight)
        self.hMM_mass[region].Fill(ev.DMDM_mass[mm_maxIxy], weight)
        self.hMM_trackIxy[region].Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight)
        self.hMM_trackIxy_log[region].Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight)
        self.hMM_trackDxy[region].Fill(ev.DMDM_trackDxy_PV[mm_maxIxy], weight)
        self.hMM_Lxy[region].Fill(abs(ev.DMDM_Lxy_PV[mm_maxIxy]), weight)
        self.hMM_Ixy[region].Fill(abs(ev.DMDM_Ixy_PV[mm_maxIxy]), weight)
        self.hMM_normalizedChi2[region].Fill(ev.DMDM_normalizedChi2[mm_maxIxy], weight)
        self.hMM_cosAlpha[region].Fill(ev.DMDM_cosAlpha[mm_maxIxy], weight)
        self.hMM_leadingPt[region].Fill(ev.DMDM_leadingPt[mm_maxIxy], weight) 
        self.hMM_vx_vy[region].Fill(ev.DMDM_vx[mm_maxIxy], ev.DMDM_vy[mm_maxIxy], weight)




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


    #### -------------------------
    #### --
    #### ---- DiElectron filling
    #### --
    #### -------------------------

    def fillDielectrons(self, ev, weight, region, ee_maxIxy, nBSEE):

        self.hEE_nPU[region].Fill(ev.nPV, weight)
        self.hEE_nBSEE[region].Fill(nBSEE, weight)  
        self.hEE_dPhi[region].Fill(abs(ev.EE_dPhi[ee_maxIxy]), weight)
        self.hEE_dPhi_inv[region].Fill(3.14 - abs(ev.EE_dPhi[ee_maxIxy]), weight)
        self.hEE_mass[region].Fill(ev.EE_mass[ee_maxIxy], weight) 
        self.hEE_trackDxy[region].Fill(ev.EE_trackDxy_PV[ee_maxIxy], weight)
        self.hEE_trackIxy[region].Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight)
        self.hEE_trackIxy_log[region].Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight)
        self.hEE_Lxy[region].Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), weight) 
        self.hEE_Ixy[region].Fill(abs(ev.EE_Ixy_PV[ee_maxIxy]), weight)
        self.hEE_leadingEt[region].Fill(ev.EE_leadingEt[ee_maxIxy], weight)
        self.hEE_subleadingEt[region].Fill(ev.EE_subleadingEt[ee_maxIxy], weight)
        self.hEE_normalizedChi2[region].Fill(ev.EE_normalizedChi2[ee_maxIxy], weight)
        self.hEE_vx_vy[region].Fill(ev.EE_vx[ee_maxIxy], ev.EE_vy[ee_maxIxy], weight)



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
        if not ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 and ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15:

            self.hEE_cutEfficiency.Fill(0, weight) 

            if not ev.PV_passAcceptance: return
            self.hEE_cutEfficiency.Fill(1, weight)

            if ev.nEE < 1: return
            self.hEE_cutEfficiency.Fill(2, weight)

            pass_ptmin = False
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

                if not eval(self.cm.EE_leadingPt45): continue
                if not eval(self.cm.EE_subleadingPt28): continue
                pass_ptmin = True

                if not eval(self.cm.EE_leadingEt45): continue
                if not eval(self.cm.EE_subleadingEt28): continue
                pass_etmin = True

                if not eval(self.cm.EE_eta2): continue
                if not eval(self.cm.EE_etanoBE): continue
                pass_etamax = True

                if not eval(self.cm.EE_normChi2_10): continue
                pass_normChi2max = True 

                if not eval(self.cm.EE_iso2l): continue
                pass_relisomax = True

                if not eval(self.cm.EE_mass15): continue
                pass_massmin = True

                if not eval(self.cm.EE_OS): continue
                pass_oscharge = True

                if not eval(self.cm.EE_dPhiforward): continue
                pass_dPhiSR = True

                if not eval(self.cm.EE_OffZ): continue
                pass_offZ = True

                if not eval(self.cm.EE_Ixy6high): continue
                pass_Ixymin = True

            if pass_ptmin: self.hEE_cutEfficiency.Fill(3, weight)
            if pass_etmin: self.hEE_cutEfficiency.Fill(4, weight)
            if pass_etamax: self.hEE_cutEfficiency.Fill(5, weight)
            if pass_normChi2max: self.hEE_cutEfficiency.Fill(6, weight)
            if pass_relisomax: self.hEE_cutEfficiency.Fill(7, weight)
            if pass_massmin: self.hEE_cutEfficiency.Fill(8, weight)
            if pass_oscharge: self.hEE_cutEfficiency.Fill(9, weight)
            if pass_dPhiSR: self.hEE_cutEfficiency.Fill(10, weight)
            if pass_offZ: self.hEE_cutEfficiency.Fill(11, weight)
            if pass_Ixymin: self.hEE_cutEfficiency.Fill(12, weight)


        ### MM Channel
        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10:

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

                if not eval(self.cm.MM_pt31): continue
                pass_ptmin = True

                if not eval(self.cm.MM_eta2): continue
                pass_etamax = True

                if not eval(self.cm.MM_normChi2_10): continue
                pass_normChi2max = True 

                if not eval(self.cm.MM_dR0p2): continue
                pass_dRmin = True

                if not eval(self.cm.MM_iso2l): continue
                pass_relisomax = True

                if not eval(self.cm.MM_cosAlpha0p8): continue
                pass_cosAlphamin = True

                if not eval(self.cm.MM_mass15): continue
                pass_massmin = True

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
            if pass_normChi2max: self.hMM_cutEfficiency.Fill(5, weight)
            if pass_dRmin: self.hMM_cutEfficiency.Fill(6, weight)
            if pass_relisomax: self.hMM_cutEfficiency.Fill(7, weight)
            if pass_cosAlphamin: self.hMM_cutEfficiency.Fill(8, weight)
            if pass_massmin: self.hMM_cutEfficiency.Fill(9, weight)
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

        self.hEE_cutEfficiency = r.TH1F('hEE_cutEfficiency' + self.sufix, ';; Events with EE / Category', 13, 0, 13)
        self.hMM_cutEfficiency = r.TH1F('hMM_cutEfficiency' + self.sufix, ';; Events with MM / Category', 14, 0, 14)






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

