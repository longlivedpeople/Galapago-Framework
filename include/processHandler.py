import math
import os
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager
import numpy as np


class processHandler:

    def __init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata):

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
        #### ---- Regions initialization
        #### ------------------------------
    
        self.dielectronRegions = [] # region name : region cut
        self.declareDielectronRegions()
        
        self.dimuonRegions = [] # region name : region cut
        self.declareDimuonRegions()


        #### -------------------------------
        #### ---- Histogram initialization
        #### -------------------------------
        exec("self.hMM_yield = r.TH1F('hMM_yield' + self.sufix, ';Passed muon trigger;', 1, 0, 1)")
        exec("self.hEE_yield = r.TH1F('hEE_yield' + self.sufix, ';Passed photon trigger;', 1, 0, 1)")

        for region in self.dielectronRegions:
            region_name = region[0]
            self.declareDielectronHistograms(region_name)
        for region in self.dimuonRegions:
            region_name = region[0]
            self.declareDimuonHistograms(region_name)



        #### --------------------------------
        #### ---- Apply sumw2 to histograms
        #### --------------------------------

        for attr, value in self.__dict__.items():
            if attr[0] == 'h': value.Sumw2()



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

        exec("self.hMM{0}_nBSMM = r.TH1F('hMM{0}_nBSMM' + self.sufix, ';Number of MM candidates;', 4, 0, 4)".format(region))
        exec("self.hMM{0}_nPU = r.TH1F('hMM{0}_nPU' + self.sufix, ';Number of true primary vertices;', 40, 0, 80)".format(region))
        exec("self.hMM{0}_dPhi = r.TH1F('hMM{0}_dPhi' + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 30, 0, 3.14)".format(region))
        exec("self.hMM{0}_dPhi_inv = r.TH1F('hMM{0}_dPhi_inv' + self.sufix, ';Dimuon inverted collinearity #pi - |#Delta#Phi|;', 30, 0, 3.14)".format(region))
        exec("self.hMM{0}_mass = r.TH1F('hMM{0}_mass' + self.sufix, ';Dimuon invariant mass m_{{#mu#mu}} (GeV);', 80, 0, 400)".format(region))
        exec("self.hMM{0}_cosAlpha = r.TH1F('hMM{0}_cosAlpha' + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}});', 22, -1.1, 1.1)".format(region))
        exec("self.hMM{0}_trackIxy = r.TH1F('hMM{0}_trackIxy' + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', 40, 0, 40)".format(region))
        exec("self.hMM{0}_trackIxy_log = r.TH1F('hMM{0}_trackIxy_log' + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', len(binlog_Ixy)-1, binlog_Ixy)".format(region))
        exec("self.hMM{0}_trackDxy = r.TH1F('hMM{0}_trackDxy' + self.sufix, ';Dimuon |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hMM{0}_Lxy = r.TH1F('hMM{0}_Lxy' + self.sufix, ';Dimuon vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hMM{0}_Ixy = r.TH1F('hMM{0}_Ixy' + self.sufix, ';Dimuon vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_leadingPt = r.TH1F('hMM{0}_leadingPt' + self.sufix, ';Dimuon leading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hMM{0}_subleadingPt = r.TH1F('hMM{0}_subleadingPt' + self.sufix, ';Dimuon subleading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hMM{0}_normalizedChi2 = r.TH1F('hMM{0}_normalizedChi2' + self.sufix, ';Dimuon vertex fit #chi^{{2}}/ndof;', 50, 0, 50)".format(region))
        exec("self.hMM{0}_cosAlpha_dPhi = r.TH2F('hMM{0}_cosAlpha_dPhi' + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}}) ; Dimuon collinearity |#Delta#Phi|', 20, -1.0, 1.0, 20, 0, 3.14)".format(region))
        exec("self.hMM{0}_vx_vy = r.TH2F('hMM{0}_vx_vy' + self.sufix, ';Dimuon vertex v_{{x}} (cm) ; Dimuon vertex v_{{y}} (cm)', 200, -4.0, 4.0, 200, -4.0, 4.0)".format(region))
        exec("self.hMM{0}_DeltaPhi_dPhi = r.TH2F('hMM{0}_DeltaPhi_dPhi' + self.sufix, ';Dimuon #Delta#Phi(ll) ; Dimuon collinearity |#Delta#Phi|', 20, 0.0, 3.14, 20, 0, 3.14)".format(region))


    def declareDielectronHistograms(self, region):

        binlog_Ixy = np.logspace(-3, 3, 50)

        #### ---------------
        #### ---- EE plots
        #### ---------------

        exec("self.hEE{0}_nBSEE = r.TH1F('hEE{0}_nBSEE' + self.sufix, ';Number of EE candidates;', 4, 0, 4)".format(region))
        exec("self.hEE{0}_nPU = r.TH1F('hEE{0}_nPU' + self.sufix, ';Number of true primary vertices;', 40, 0, 80)".format(region))
        exec("self.hEE{0}_dPhi = r.TH1F('hEE{0}_dPhi' + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 30, 0, 3.14)".format(region))
        exec("self.hEE{0}_dPhi_inv = r.TH1F('hEE{0}_dPhi_inv' + self.sufix, ';Dielectron inverted collinearity #pi - |#Delta#Phi|;', 30, 0, 3.14)".format(region))
        exec("self.hEE{0}_mass = r.TH1F('hEE{0}_mass' + self.sufix, ';Dielectron invariant mass m_{{ee}} (GeV);', 80, 0, 400)".format(region))
        exec("self.hEE{0}_cosAlpha = r.TH1F('hEE{0}_cosAlpha' + self.sufix, ';Dielectron cos(#alpha_{{ee}});', 22, -1.1, 1.1)".format(region))
        exec("self.hEE{0}_trackIxy = r.TH1F('hEE{0}_trackIxy' + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', 40, 0, 40)".format(region))
        exec("self.hEE{0}_trackIxy_log = r.TH1F('hEE{0}_trackIxy_log' + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', len(binlog_Ixy)-1, binlog_Ixy)".format(region))
        exec("self.hEE{0}_trackDxy = r.TH1F('hEE{0}_trackDxy' + self.sufix, ';Dielectron |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hEE{0}_Lxy = r.TH1F('hEE{0}_Lxy' + self.sufix, ';Dielectron vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hEE{0}_Ixy = r.TH1F('hEE{0}_Ixy' + self.sufix, ';Dielectron vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_leadingPt = r.TH1F('hEE{0}_leadingPt' + self.sufix, ';Dielectron leading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_subleadingPt = r.TH1F('hEE{0}_subleadingPt' + self.sufix, ';Dielectron subleading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_leadingEt = r.TH1F('hEE{0}_leadingEt' + self.sufix, ';Dielectron leading E_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_subleadingEt = r.TH1F('hEE{0}_subleadingEt' + self.sufix, ';Dielectron subleading E_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_normalizedChi2 = r.TH1F('hEE{0}_normalizedChi2' + self.sufix, ';Dielectron vertex fit #chi^{{2}}/ndof;', 50, 0, 50)".format(region))
        exec("self.hEE{0}_cosAlpha_dPhi = r.TH2F('hEE{0}_cosAlpha_dPhi' + self.sufix, ';Dielectron cos(#alpha_{{ee}}) ; Dielectron collinearity |#Delta#Phi|', 20, -1.0, 1.0, 20, 0, 3.14)".format(region))
        exec("self.hEE{0}_vx_vy = r.TH2F('hEE{0}_vx_vy' + self.sufix, ';Dielectron vertex v_{{x}} (cm) ; Dielectron vertex v_{{y}} (cm)', 200, -4.0, 4.0, 200, -4.0, 4.0)".format(region))
        exec("self.hEE{0}_Lxy_dPhi = r.TH2F('hEE{0}_Lxy_dPhi' + self.sufix, ';Dielectron |L_{{xy}}| (cm) ; Dielectron collinearity |#Delta#Phi|', 20, 0.0, 1.0, 20, 0, 3.14)".format(region))
        exec("self.hEE{0}_DeltaPhi_dPhi = r.TH2F('hEE{0}_DeltaPhi_dPhi' + self.sufix, ';Dielectron #Delta#Phi(ll) ; Dielectron collinearity |#Delta#Phi|', 20, 0.0, 3.14, 20, 0, 3.14)".format(region))




    #### ------------------------
    #### --
    #### ---- declareRegions
    #### --
    #### ------------------------

    def declareDielectronRegions(self):

        #### ---------------------------
        #### ---- EE Region definition
        #### ---------------------------

        EE_BScut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l])
        EE_OSSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_dPhiforward])
        EE_OSCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_dPhibackward])
        EE_offZCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OffZ, self.cm.EE_dPhibackward])
        EE_offZSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OffZ, self.cm.EE_dPhiforward])
        EE_onZCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OnZ, self.cm.EE_dPhibackward])
        EE_onZSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OnZ, self.cm.EE_dPhiforward])
        EE_promptSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_Ixy6prompt, self.cm.EE_dPhiforward])
        EE_promptCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_Ixy6prompt, self.cm.EE_dPhibackward])
        EE_dispSRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_Ixy6high, self.cm.EE_dPhiforward])
        EE_dispCRcut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_Ixy6high, self.cm.EE_dPhibackward])
        EE_offZCR_ne1_cut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OffZ, self.cm.EE_dPhibackward, self.cm.EE_nBSEEe1])
        EE_offZSR_ne1_cut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OffZ, self.cm.EE_dPhiforward, self.cm.EE_nBSEEe1])
        EE_offZCR_ng1_cut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OffZ, self.cm.EE_dPhibackward, self.cm.EE_nBSEEg1])
        EE_offZSR_ng1_cut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_OffZ, self.cm.EE_dPhiforward, self.cm.EE_nBSEEg1])
        EE_fullZCR_ng1_cut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_dPhibackward, self.cm.EE_nBSEEg1])
        EE_fullZSR_ng1_cut = self.cm.AddList([self.cm.EE_OS, self.cm.EE_iso2l, self.cm.EE_dPhiforward, self.cm.EE_nBSEEg1])


        self.dielectronRegions.append(['BS', EE_BScut])
        self.dielectronRegions.append(['CROS', EE_OSCRcut])
        self.dielectronRegions.append(['SROS', EE_OSSRcut])
        self.dielectronRegions.append(['promptCR', EE_promptCRcut])
        self.dielectronRegions.append(['promptSR', EE_promptSRcut])
        #self.dielectronRegions.append(['dispCR', EE_dispCRcut])
        #self.dielectronRegions.append(['dispSR', EE_dispSRcut])
        self.dielectronRegions.append(['onZCR', EE_onZCRcut])
        self.dielectronRegions.append(['onZSR', EE_onZSRcut])
        self.dielectronRegions.append(['offZCR', EE_offZCRcut])
        self.dielectronRegions.append(['offZSR', EE_offZSRcut])
        #self.dielectronRegions.append(['offZCRne1', EE_offZCR_ne1_cut])
        #self.dielectronRegions.append(['offZSRne1', EE_offZSR_ne1_cut])
        #self.dielectronRegions.append(['offZCRng1', EE_offZCR_ng1_cut])
        #self.dielectronRegions.append(['offZSRng1', EE_offZSR_ng1_cut])
        #self.dielectronRegions.append(['fullZCRng1', EE_fullZCR_ng1_cut])
        #self.dielectronRegions.append(['fullZSRng1', EE_fullZSR_ng1_cut])


    def declareDimuonRegions(self):


        #### ---------------------------
        #### ---- MM Region definition
        #### ---------------------------

        MM_BScut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l])
        MM_OSSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_dPhiforward])
        MM_OSCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_dPhibackward])
        MM_offZCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OffZ, self.cm.MM_dPhibackward])
        MM_offZSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OffZ, self.cm.MM_dPhiforward])
        MM_onZCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OnZ, self.cm.MM_dPhibackward])
        MM_onZSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OnZ, self.cm.MM_dPhiforward])
        MM_promptSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_Ixy6prompt, self.cm.MM_dPhiforward])
        MM_promptCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_Ixy6prompt, self.cm.MM_dPhibackward])
        MM_dispSRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_Ixy6high, self.cm.MM_dPhiforward])
        MM_dispCRcut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_Ixy6high, self.cm.MM_dPhibackward])
        MM_offZCR_ne1_cut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OffZ, self.cm.MM_dPhibackward, self.cm.MM_nBSMMe1])
        MM_offZSR_ne1_cut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OffZ, self.cm.MM_dPhiforward, self.cm.MM_nBSMMe1])
        MM_offZCR_ng1_cut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OffZ, self.cm.MM_dPhibackward, self.cm.MM_nBSMMg1])
        MM_offZSR_ng1_cut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_OffZ, self.cm.MM_dPhiforward, self.cm.MM_nBSMMg1])
        MM_fullZCR_ng1_cut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_dPhibackward, self.cm.MM_nBSMMg1])
        MM_fullZSR_ng1_cut = self.cm.AddList([self.cm.MM_OS, self.cm.MM_iso2l, self.cm.MM_dPhiforward, self.cm.MM_nBSMMg1])


        self.dimuonRegions.append(['BS', MM_BScut])
        self.dimuonRegions.append(['CROS', MM_OSCRcut])
        self.dimuonRegions.append(['SROS', MM_OSSRcut])
        self.dimuonRegions.append(['promptCR', MM_promptCRcut])
        self.dimuonRegions.append(['promptSR', MM_promptSRcut])
        #self.dimuonRegions.append(['dispCR', MM_dispCRcut])
        #self.dimuonRegions.append(['dispSR', MM_dispSRcut])
        self.dimuonRegions.append(['onZCR', MM_onZCRcut])
        self.dimuonRegions.append(['onZSR', MM_onZSRcut])
        self.dimuonRegions.append(['offZCR', MM_offZCRcut])
        self.dimuonRegions.append(['offZSR', MM_offZSRcut])
        #self.dimuonRegions.append(['offZCRne1', MM_offZCR_ne1_cut])
        #self.dimuonRegions.append(['offZSRne1', MM_offZSR_ne1_cut])
        #self.dimuonRegions.append(['offZCRng1', MM_offZCR_ng1_cut])
        #self.dimuonRegions.append(['offZSRng1', MM_offZSR_ng1_cut])
        #self.dimuonRegions.append(['fullZCRng1', MM_fullZCR_ng1_cut])
        #self.dimuonRegions.append(['fullZSRng1', MM_fullZSR_ng1_cut])




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
        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10:
            self.hMM_yield.Fill(0, weight)
        if not ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 and ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15:
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
            print('There is some collections missed in this file: Dimuon histograms will be empty')
        
        
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
            print('There is some collections missed in this file: Dielectron histograms will be empty')
        



       
    #### ------------------------
    #### --
    #### ---- DiMuon Processing
    #### --
    #### ------------------------

    def processDimuons(self, ev):


        #### ----------------------------
        #### ---- Select a MM candidate
        #### ----------------------------

        ### -> Events just passing the trigger and having a valid DG pair
        if not ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 or ev.nDMDM == 0: return -99, -99

        ### -> Find the MM pair with maximum Ixy that passes the cuts
        mm_maxIxy = -99
        maxIxy = -1
        nBSMM = 0

        for i in range(0, ev.nDMDM):

            imm = i
            if not eval(self.cm.MM_ID): continue
            if not eval(self.cm.MM_pt31): continue
            if not eval(self.cm.MM_eta2): continue
            if not eval(self.cm.MM_cosAlpha0p8): continue
            if not eval(self.cm.MM_mass15): continue
            if not eval(self.cm.MM_normChi2_5): continue
            if not eval(self.cm.MM_dR0p2): continue 

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

        exec("self.hMM{0}_nPU.Fill(ev.nPV, weight)".format(region))
        exec("self.hMM{0}_nBSMM.Fill(nBSMM, weight)".format(region))
        exec("self.hMM{0}_dPhi.Fill(abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_dPhi_inv.Fill(3.14 - abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_mass.Fill(ev.DMDM_mass[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackIxy.Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackIxy_log.Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackDxy.Fill(ev.DMDM_trackDxy_PV[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_Lxy.Fill(abs(ev.DMDM_Lxy_PV[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_normalizedChi2.Fill(ev.DMDM_normalizedChi2[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_Ixy.Fill(abs(ev.DMDM_Ixy_PV[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_cosAlpha.Fill(ev.DMDM_cosAlpha[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_leadingPt.Fill(ev.DMDM_leadingPt[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_cosAlpha_dPhi.Fill(ev.DMDM_cosAlpha[mm_maxIxy], abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_vx_vy.Fill(ev.DMDM_vx[mm_maxIxy], ev.DMDM_vy[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_DeltaPhi_dPhi.Fill(ev.DMDM_lldPhi[mm_maxIxy], abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))



    #### ----------------------------
    #### --
    #### ---- DiElectron Processing
    #### --
    #### ----------------------------

    def processDielectrons(self, ev):

        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 or not ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 or ev.nEE == 0: return -99, -99

        ### -> Find the EE pair with maximum Ixy that passes the cuts
        ee_maxIxy = -99
        maxIxy = -1
        nBSEE = 0
        for i in range(0, ev.nEE):
         
            iee = i
            if not eval(self.cm.EE_eta2): continue
            if not eval(self.cm.EE_leadingPt45): continue
            if not eval(self.cm.EE_subleadingPt28): continue
            if not eval(self.cm.EE_leadingEt45): continue
            if not eval(self.cm.EE_subleadingEt28): continue
            if not eval(self.cm.EE_normChi2_7): continue
            if not eval(self.cm.EE_mass15): continue
           
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

        exec("self.hEE{0}_nPU.Fill(ev.nPV, weight)".format(region))
        exec("self.hEE{0}_nBSEE.Fill(nBSEE, weight)".format(region))
        exec("self.hEE{0}_dPhi.Fill(abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_dPhi_inv.Fill(3.14 - abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_mass.Fill(ev.EE_mass[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackDxy.Fill(ev.EE_trackDxy_PV[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackIxy.Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackIxy_log.Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_Lxy.Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_Ixy.Fill(abs(ev.EE_Ixy_PV[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_leadingPt.Fill(ev.EE_leadingPt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_subleadingPt.Fill(ev.EE_subleadingPt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_leadingEt.Fill(ev.EE_leadingEt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_subleadingEt.Fill(ev.EE_subleadingEt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_normalizedChi2.Fill(ev.EE_normalizedChi2[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_cosAlpha_dPhi.Fill(ev.EE_cosAlpha[ee_maxIxy], abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_vx_vy.Fill(ev.EE_vx[ee_maxIxy], ev.EE_vy[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_Lxy_dPhi.Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_DeltaPhi_dPhi.Fill(ev.EE_lldPhi[ee_maxIxy], abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))


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
                if attr[0] == 'h': value.Write()

        
        output.Close()

