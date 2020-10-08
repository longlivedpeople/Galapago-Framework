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
        self.cutManager = CutManager.CutManager()
        self.sufix = '__{0}__{1}__{2}__{3}'.format(treename, blockname, samplename, str(samplenumber))


        #### ---------------------------
        #### ---- EE Region definition
        #### ---------------------------

        EE_OScut = 'ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxA[ee_maxIxy]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxB[ee_maxIxy]]] < 0'
        EE_SScut = 'ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxA[ee_maxIxy]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxB[ee_maxIxy]]] > 0'
        EE_IsoAcut = '(ev.EEBase_relisoA[ee_maxIxy] < 0.1)'
        EE_IsoBcut = '(ev.EEBase_relisoB[ee_maxIxy] < 0.1)'
        EE_SS0cut = self.cutManager.AddList([EE_SScut, EE_IsoAcut, EE_IsoBcut])
        EE_OS0cut = self.cutManager.AddList([EE_OScut, EE_IsoAcut, EE_IsoBcut])
        EE_SSIcut = self.cutManager.OR( self.cutManager.AddList([EE_SScut, 'not '+EE_IsoAcut, EE_IsoBcut]), self.cutManager.AddList([EE_SScut, EE_IsoAcut, 'not '+EE_IsoBcut]) )
        EE_OSIcut = self.cutManager.OR( self.cutManager.AddList([EE_OScut, 'not '+EE_IsoAcut, EE_IsoBcut]), self.cutManager.AddList([EE_OScut, EE_IsoAcut, 'not '+EE_IsoBcut]) )
        EE_SSIIcut = self.cutManager.AddList([EE_SScut, 'not '+EE_IsoAcut, 'not '+EE_IsoBcut])
        EE_OSIIcut = self.cutManager.AddList([EE_OScut, 'not '+EE_IsoAcut, 'not '+EE_IsoBcut])

        self.dielectronRegions = [] # region name : region cut
        self.dielectronRegions.append(['BaseLine', '1'])
        self.dielectronRegions.append(['SS0', EE_SS0cut])
        self.dielectronRegions.append(['OS0', EE_OS0cut])
        #self.dielectronRegions.append(['SSI', EE_SSIcut])
        #self.dielectronRegions.append(['OSI', EE_OSIcut])
        #self.dielectronRegions.append(['SSII', EE_SSIIcut])
        #self.dielectronRegions.append(['OSII', EE_OSIIcut])



        #### ---------------------------
        #### ---- MM Region definition
        #### ---------------------------

        MM_OScut = 'ev.DGM_charge[ev.DMDMBase_idxA[mm_maxIxy]]*ev.DGM_charge[ev.DMDMBase_idxB[mm_maxIxy]] < 0'
        MM_SScut = 'ev.DGM_charge[ev.DMDMBase_idxA[mm_maxIxy]]*ev.DGM_charge[ev.DMDMBase_idxB[mm_maxIxy]] > 0'
        MM_IsoAcut = '(ev.DGM_relPFiso[ev.DMDMBase_idxA[mm_maxIxy]] < 0.1)'
        MM_IsoBcut = '(ev.DGM_relPFiso[ev.DMDMBase_idxB[mm_maxIxy]] < 0.1)'
        MM_SS0cut = self.cutManager.AddList([MM_SScut, MM_IsoAcut, MM_IsoBcut])
        MM_OS0cut = self.cutManager.AddList([MM_OScut, MM_IsoAcut, MM_IsoBcut])
        MM_SSIcut = self.cutManager.OR( self.cutManager.AddList([MM_SScut, 'not '+MM_IsoAcut, MM_IsoBcut]), self.cutManager.AddList([MM_SScut, MM_IsoAcut, 'not '+MM_IsoBcut]) )
        MM_OSIcut = self.cutManager.OR( self.cutManager.AddList([MM_OScut, 'not '+MM_IsoAcut, MM_IsoBcut]), self.cutManager.AddList([MM_OScut, MM_IsoAcut, 'not '+MM_IsoBcut]) )
        MM_SSIIcut = self.cutManager.AddList([MM_SScut, 'not '+MM_IsoAcut, 'not '+MM_IsoBcut])
        MM_OSIIcut = self.cutManager.AddList([MM_OScut, 'not '+MM_IsoAcut, 'not '+MM_IsoBcut])

        self.dimuonRegions = [] # region name : region cut
        self.dimuonRegions.append(['BaseLine', '1'])
        self.dimuonRegions.append(['SS0', MM_SS0cut])
        self.dimuonRegions.append(['OS0', MM_OS0cut])
        #self.dimuonRegions.append(['SSI', MM_SSIcut])
        #self.dimuonRegions.append(['OSI', MM_OSIcut])
        #self.dimuonRegions.append(['SSII', MM_SSIIcut])
        #self.dimuonRegions.append(['OSII', MM_OSIIcut])



        #### -------------------------------
        #### ---- Histogram initialization
        #### -------------------------------
        for region in self.dimuonRegions:
            region_name = region[0]
            self.declareHistograms(region_name)



        #### --------------------------------
        #### ---- Apply sumw2 to histograms
        #### --------------------------------

        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Sumw2()



    #### ----------------------------
    #### --
    #### ---- Core event processing
    #### --
    #### ----------------------------

    def processEvent(self, ev):

        if not self.isdata:
            weight = self.lumiweight*ev.wPU*ev.genWeight/abs(ev.genWeight)
            weightnoPU = self.lumiweight*ev.genWeight/abs(ev.genWeight)
        else: 
            weight = 1


        #### N tracks cut
        nTracks = ev.RefittedPV_nPFTrack + ev.RefittedPV_nLostTrack + ev.RefittedPV_nExcludedTrack
        if nTracks < 4: return

        #self.processDimuons(ev, weight)
        #self.processDielectrons(ev, weight)
        
        try:
            self.processDimuons(ev, weight)
        except AttributeError:
            print('There is some collections missed in this file: Dimuon histograms will be empty')

        try:
            self.processDielectrons(ev, weight)
        except AttributeError:
            print('There is some collections missed in this file: Dielectron histograms will be empty')
        



       
    #### ------------------------
    #### --
    #### ---- DiMuon Processing
    #### --
    #### ------------------------

    def processDimuons(self, ev, weight):


        #### ----------------------------
        #### ---- Select a MM candidate
        #### ----------------------------

        ### -> Events just passing the trigger and having a valid DG pair
        if not ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 or ev.nDMDMBase == 0: return

        ### -> Find the MM pair with maximum Ixy that passes the cuts
        mm_maxIxy = -99
        maxIxy = -1
        nGoodMM = 0
        for i in range(0, ev.nDMDMBase):
            if abs(ev.DGM_eta[ev.DMDMBase_idxA[i]]) > 2: continue
            if abs(ev.DGM_eta[ev.DMDMBase_idxB[i]]) > 2: continue
            if ev.DGM_pt[ev.DMDMBase_idxA[i]] < 31: continue
            if ev.DGM_pt[ev.DMDMBase_idxB[i]] < 31: continue
            if ev.DGM_chi2[ev.DMDMBase_idxA[i]]/ev.DGM_ndof[ev.DMDMBase_idxA[i]] > 10: continue
            if ev.DGM_chi2[ev.DMDMBase_idxB[i]]/ev.DGM_ndof[ev.DMDMBase_idxB[i]] > 10: continue
            if ev.DGM_numberOfValidHits[ev.DMDMBase_idxA[i]] < 22: continue
            if ev.DGM_numberOfValidHits[ev.DMDMBase_idxB[i]] < 22: continue
            if ev.DMDMBase_cosAlpha[i] < -0.80: continue
            #if abs(90 - ev.DMDMBase_mass[i]) < 10: continue

            nGoodMM+=1

            if ev.DMDMBase_trackIxy[i] > maxIxy:
                maxIxy = ev.DMDMBase_trackIxy[i]
                mm_maxIxy = i

        # If not candidate is encountered, we skip the event
        if not nGoodMM: return 


        #### ---------------------
        #### ---- Region filling
        #### ---------------------

        for region in self.dimuonRegions:
            if eval(region[1]):
                self.fillDimuons(ev, weight, region[0], mm_maxIxy, nGoodMM)




    #### ---------------------
    #### --
    #### ---- DiMuon filling
    #### --
    #### ---------------------

    def fillDimuons(self, ev, weight, region, mm_maxIxy, nGoodMM):

        exec("self.hMM{0}_nGoodMM.Fill(nGoodMM, weight)".format(region))
        exec("self.hMM{0}_dPhi.Fill(abs(ev.DMDMBase_dPhi[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_mass.Fill(ev.DMDMBase_mass[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackIxy.Fill(ev.DMDMBase_trackIxy[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackDxy.Fill(ev.DMDMBase_trackDxy[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_Lxy.Fill(abs(ev.DMDMBase_Lxy[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_Ixy.Fill(abs(ev.DMDMBase_Ixy[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_cosAlpha.Fill(ev.DMDMBase_cosAlpha[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_leadingPt.Fill(ev.DMDMBase_leadingPt[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_dPhi.Fill(ev.DMDMBase_dPhi[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_cosAlpha_dPhi.Fill(ev.DMDMBase_cosAlpha[mm_maxIxy], abs(ev.DMDMBase_dPhi[mm_maxIxy]), weight)".format(region))

        for i in range(0, ev.nDGM):

            exec("self.hM{0}_pt.Fill(ev.DGM_pt[i], weight)".format(region))
            exec("self.hM{0}_eta.Fill(ev.DGM_eta[i], weight)".format(region))
            exec("self.hM{0}_dxy.Fill(ev.DGM_dxy[i], weight)".format(region))
            exec("self.hM{0}_dxyError.Fill(ev.DGM_dxyError[i], weight)".format(region))
            exec("self.hM{0}_charge.Fill(ev.DGM_charge[i], weight)".format(region))
            exec("self.hM{0}_numberOfValidHits.Fill(ev.DGM_numberOfValidHits[i], weight)".format(region))
            exec("self.hM{0}_normChi2.Fill(ev.DGM_chi2[i]/ev.DGM_ndof[i], weight)".format(region))
            exec("self.hM{0}_reliso.Fill(ev.DGM_relPFiso[i], weight)".format(region))



    #### ----------------------------
    #### --
    #### ---- DiElectron Processing
    #### --
    #### ----------------------------

    def processDielectrons(self, ev, weight):

        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 or not ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 or ev.nEEBase == 0: return

        ### -> Find the EE pair with maximum Ixy that passes the cuts
        ee_maxIxy = -99
        maxIxy = -1
        nGoodEE = 0
        for i in range(0, ev.nEEBase):
            if abs(ev.ElectronCandidate_eta[ev.EEBase_idxA[i]]) > 2.0: continue
            if abs(ev.ElectronCandidate_eta[ev.EEBase_idxB[i]]) > 2.0: continue
            if ev.EEBase_leadingPt[i] < 48: continue
            if ev.EEBase_leadingEt[i] < 48: continue
            if ev.EEBase_subleadingPt[i] < 28: continue
            if ev.EEBase_subleadingEt[i] < 28: continue
            #if abs(90 - ev.EEBase_mass[i]) < 10: continue

            nGoodEE+=1

            if ev.EEBase_trackIxy[i] > maxIxy:
                maxIxy = ev.EEBase_trackIxy[i]
                ee_maxIxy = i


        if not nGoodEE: return ## If no candidate is encountered, the event is not filled


        #### ---------------------
        #### ---- Region filling
        #### ---------------------

        for region in self.dielectronRegions:
            if eval(region[1]):
                self.fillDielectrons(ev, weight, region[0], ee_maxIxy, nGoodEE)




    #### -------------------------
    #### --
    #### ---- DiElectron filling
    #### --
    #### -------------------------

    def fillDielectrons(self, ev, weight, region, ee_maxIxy, nGoodEE):


        exec("self.hEE{0}_nGoodEE.Fill(nGoodEE, weight)".format(region))
        exec("self.hEE{0}_dPhi.Fill(abs(ev.EEBase_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_mass.Fill(ev.EEBase_mass[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackDxy.Fill(ev.EEBase_trackDxy[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackIxy.Fill(ev.EEBase_trackIxy[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_Lxy.Fill(ev.EEBase_Lxy[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_Ixy.Fill(ev.EEBase_Ixy[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_leadingPt.Fill(ev.EEBase_leadingPt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_subleadingPt.Fill(ev.EEBase_subleadingPt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_leadingEt.Fill(ev.EEBase_leadingEt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_subleadingEt.Fill(ev.EEBase_subleadingEt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_normalizedChi2.Fill(ev.EEBase_normalizedChi2[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_cosAlpha_dPhi.Fill(ev.EEBase_cosAlpha[ee_maxIxy], abs(ev.EEBase_dPhi[ee_maxIxy]), weight)".format(region))

        for i in range(0, ev.nElectronCandidate):

            exec("self.hE{0}_eta.Fill(ev.ElectronCandidate_eta[i], weight)".format(region))
            exec("self.hE{0}_pt.Fill(ev.ElectronCandidate_pt[i], weight)".format(region))
            exec("self.hE{0}_dxy.Fill(ev.ElectronCandidate_dxy[i], weight)".format(region))
            exec("self.hE{0}_dxyError.Fill(ev.ElectronCandidate_dxyError[i], weight)".format(region))
            exec("self.hE{0}_charge.Fill(ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[i]], weight)".format(region))
            exec("self.hE{0}_reliso.Fill(ev.IsoTrackSel_pfIsolationDR03[ev.ElectronCandidate_isotrackIdx[i]], weight)".format(region))





    #### ------------------------
    #### --
    #### ---- declareHistograms
    #### --
    #### ------------------------

    def declareHistograms(self, region):

        #### ---------------
        #### ---- MM plots
        #### ---------------
        exec("self.hMM{0}_nGoodMM = r.TH1F('hMM{0}_nGoodMM' + self.sufix, ';Number of MM candidates;', 0, 0, 4)".format(region))
        exec("self.hMM{0}_nPU = r.TH1F('hMM{0}_nPU' + self.sufix, ';Number of true primary vertices;', 40, 0, 80)".format(region))
        exec("self.hMM{0}_dPhi = r.TH1F('hMM{0}_dPhi' + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 29, 0, 3.14)".format(region))
        exec("self.hMM{0}_mass = r.TH1F('hMM{0}_mass' + self.sufix, ';Dimuon invariant mass m_{{#mu#mu}} (GeV);', 55, 0, 500)".format(region))
        exec("self.hMM{0}_cosAlpha = r.TH1F('hMM{0}_cosAlpha' + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}});', 22, -1.1, 1.1)".format(region))
        exec("self.hMM{0}_trackIxy = r.TH1F('hMM{0}_trackIxy' + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_trackDxy = r.TH1F('hMM{0}_trackDxy' + self.sufix, ';Dimuon |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hMM{0}_Lxy = r.TH1F('hMM{0}_Lxy' + self.sufix, ';Dimuon vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hMM{0}_Ixy = r.TH1F('hMM{0}_Ixy' + self.sufix, ';Dimuon vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_leadingPt = r.TH1F('hMM{0}_leadingPt' + self.sufix, ';Dimuon leading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hMM{0}_subleadingPt = r.TH1F('hMM{0}_subleadingPt' + self.sufix, ';Dimuon subleading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hMM{0}_normalizedChi2 = r.TH1F('hMM{0}_normalizedChi2' + self.sufix, ';Dimuon vertex fit #chi^{{2}}/ndof;', 20, 0, 10)".format(region))
        exec("self.hMM{0}_cosAlpha_dPhi = r.TH2F('hMM{0}_cosAlpha_dPhi' + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}}) ; Dimuon collinearity |#Delta#Phi|', 20, -1.0, 1.0, 20, 0, 3.14)".format(region))


        #### ----------------------
        #### ---- Single DG plots
        #### ----------------------

        exec("self.hM{0}_pt = r.TH1F('hM{0}_pt' + self.sufix, ';DG p_{{T}};', 50, 0, 200)".format(region))
        exec("self.hM{0}_eta = r.TH1F('hM{0}_eta' + self.sufix, ';DG #eta;', 27, -2.6, 2.6)".format(region))
        exec("self.hM{0}_dxy = r.TH1F('hM{0}_dxy' + self.sufix, ';DG transverse impact parameter |d_{{xy}}|;', 30, 0.0, 0.5)".format(region))
        exec("self.hM{0}_dxyError = r.TH1F('hM{0}_dxyError' + self.sufix, ';DG #sigma_{{d}};', 70, 0.0, 0.01)".format(region))
        exec("self.hM{0}_charge = r.TH1F('hM{0}_charge' + self.sufix, ';DG charge q;', 5, -2.5, 2.5)".format(region))
        exec("self.hM{0}_numberOfValidHits = r.TH1F('hM{0}_numberOfValidHits' + self.sufix, ';DG N_{{Hits}};', 80, 0.0, 80.0)".format(region))
        exec("self.hM{0}_normChi2 = r.TH1F('hM{0}_normChi2' + self.sufix, ';DG #chi^{{2}}/ndof;', 30, 0.0, 15.0)".format(region))
        exec("self.hM{0}_reliso = r.TH1F('hM{0}_reliso' + self.sufix, ';DG RelIso;', 15, 0.0, 0.15)".format(region))

        #### ---------------
        #### ---- EE plots
        #### ---------------

        exec("self.hEE{0}_nGoodEE = r.TH1F('hEE{0}_nGoodEE' + self.sufix, ';Number of EE candidates;', 0, 0, 4)".format(region))
        exec("self.hEE{0}_nPU = r.TH1F('hEE{0}_nPU' + self.sufix, ';Number of true primary vertices;', 40, 0, 80)".format(region))
        exec("self.hEE{0}_dPhi = r.TH1F('hEE{0}_dPhi' + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 29, 0, 3.14)".format(region))
        exec("self.hEE{0}_mass = r.TH1F('hEE{0}_mass' + self.sufix, ';Dielectron invariant mass m_{{ee}} (GeV);', 55, 0, 500)".format(region))
        exec("self.hEE{0}_cosAlpha = r.TH1F('hEE{0}_cosAlpha' + self.sufix, ';Dielectron cos(#alpha_{{ee}});', 22, -1.1, 1.1)".format(region))
        exec("self.hEE{0}_trackIxy = r.TH1F('hEE{0}_trackIxy' + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_trackDxy = r.TH1F('hEE{0}_trackDxy' + self.sufix, ';Dielectron |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hEE{0}_Lxy = r.TH1F('hEE{0}_Lxy' + self.sufix, ';Dielectron vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hEE{0}_Ixy = r.TH1F('hEE{0}_Ixy' + self.sufix, ';Dielectron vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_leadingPt = r.TH1F('hEE{0}_leadingPt' + self.sufix, ';Dielectron leading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_subleadingPt = r.TH1F('hEE{0}_subleadingPt' + self.sufix, ';Dielectron subleading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_leadingEt = r.TH1F('hEE{0}_leadingEt' + self.sufix, ';Dielectron leading E_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_subleadingEt = r.TH1F('hEE{0}_subleadingEt' + self.sufix, ';Dielectron subleading E_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_normalizedChi2 = r.TH1F('hEE{0}_normalizedChi2' + self.sufix, ';Dielectron vertex fit #chi^{{2}}/ndof;', 20, 0, 10)".format(region))
        exec("self.hEE{0}_cosAlpha_dPhi = r.TH2F('hMM{0}_cosAlpha_dPhi' + self.sufix, ';Dielectron cos(#alpha_{{ee}}) ; Dielectron collinearity |#Delta#Phi|', 20, -1.0, 1.0, 20, 0, 3.14)".format(region))

        #### ---------------------
        #### ---- Single E plots
        #### ---------------------

        exec("self.hE{0}_pt = r.TH1F('hE{0}_pt' + self.sufix, ';Electron p_{{T}};', 50, 0, 200)".format(region))
        exec("self.hE{0}_eta = r.TH1F('hE{0}_eta' + self.sufix, ';Electron #eta;', 27, -2.6, 2.6)".format(region))
        exec("self.hE{0}_dxy = r.TH1F('hE{0}_dxy' + self.sufix, ';Electron transverse impact parameter |d_{{xy}}|;', 30, 0.0, 0.5)".format(region))
        exec("self.hE{0}_dxyError = r.TH1F('hE{0}_dxyError' + self.sufix, ';Electron #sigma_{{d}};', 70, 0.0, 0.01)".format(region))
        exec("self.hE{0}_charge = r.TH1F('hE{0}_charge' + self.sufix, ';Electron charge q;', 5, -2.5, 2.5)".format(region))
        exec("self.hE{0}_reliso = r.TH1F('hE{0}_reliso' + self.sufix, ';Electron RelIso;', 15, 0.0, 0.15)".format(region))



    def Write(self):
    
        while True:
            d = os.path.dirname(self.filename)
            try:
                if not os.path.exists(d): os.makedirs(d)
                break
            except OSError, e:
                if e.errno != os.errno.EEXIST:
                    raise
                print("Sleeping...")
                time.sleep(1.0)
                pass

        output = r.TFile(self.filename, 'RECREATE')

        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Write()

        
        output.Close()


