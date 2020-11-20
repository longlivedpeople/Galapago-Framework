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

        EE_OScut = 'ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[ee_maxIxy]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[ee_maxIxy]]] < 0'
        EE_SScut = 'ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[ee_maxIxy]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[ee_maxIxy]]] > 0'
        EE_IsoAcut = '(ev.IsoTrackSel_pfIsolationDR03[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[ee_maxIxy]]]/ev.IsoTrackSel_pt[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[ee_maxIxy]]] < 0.2)'
        EE_IsoBcut = '(ev.IsoTrackSel_pfIsolationDR03[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[ee_maxIxy]]]/ev.IsoTrackSel_pt[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[ee_maxIxy]]] < 0.2)'
        EE_SRdPhicut = '(abs(ev.EE_dPhi[ee_maxIxy]) < 3.14/2)'
        EE_CRdPhicut = '(abs(ev.EE_dPhi[ee_maxIxy]) > 3.14/2)'
        EE_masked = '(abs(90 - ev.EE_mass[ee_maxIxy]) > 10)'
        EE_disp = '(abs(ev.EE_trackIxy_PV[ee_maxIxy]) > 4)'
        EE_SS0cut = self.cutManager.AddList([EE_SScut, EE_IsoAcut, EE_IsoBcut])
        EE_OS0cut = self.cutManager.AddList([EE_OScut, EE_IsoAcut, EE_IsoBcut])
        EE_SSIcut = self.cutManager.OR( self.cutManager.AddList([EE_SScut, 'not '+EE_IsoAcut, EE_IsoBcut]), self.cutManager.AddList([EE_SScut, EE_IsoAcut, 'not '+EE_IsoBcut]) )
        EE_OSIcut = self.cutManager.OR( self.cutManager.AddList([EE_OScut, 'not '+EE_IsoAcut, EE_IsoBcut]), self.cutManager.AddList([EE_OScut, EE_IsoAcut, 'not '+EE_IsoBcut]) )
        EE_Icut = self.cutManager.OR(self.cutManager.AddList(['not '+EE_IsoAcut, EE_IsoBcut]), self.cutManager.AddList([EE_IsoAcut, 'not '+EE_IsoBcut]) )
        EE_SSIIcut = self.cutManager.AddList([EE_SScut, 'not '+EE_IsoAcut, 'not '+EE_IsoBcut])
        EE_OSIIcut = self.cutManager.AddList([EE_OScut, 'not '+EE_IsoAcut, 'not '+EE_IsoBcut])
        EE_IIcut = self.cutManager.AddList(['not '+EE_IsoAcut, 'not '+EE_IsoBcut])
        EE_OSCRcut = self.cutManager.AddList([EE_OScut, EE_IsoAcut, EE_IsoBcut, EE_CRdPhicut])
        EE_SSCRcut = self.cutManager.AddList([EE_SScut, EE_IsoAcut, EE_IsoBcut, EE_CRdPhicut])
        EE_OSSRcut = self.cutManager.AddList([EE_OScut, EE_IsoAcut, EE_IsoBcut, EE_SRdPhicut])
        EE_SRIcut = self.cutManager.AddList([EE_SRdPhicut, EE_Icut])
        EE_CRIcut = self.cutManager.AddList([EE_CRdPhicut, EE_Icut])
        EE_SRIIcut = self.cutManager.AddList([EE_SRdPhicut, EE_IIcut])
        EE_CRIIcut = self.cutManager.AddList([EE_CRdPhicut, EE_IIcut])
        EE_OS0Dispcut = self.cutManager.AddList([EE_OScut, EE_IsoAcut, EE_IsoBcut, EE_disp])
        EE_SSCRmaskedcut = self.cutManager.AddList([EE_SScut, EE_IsoAcut, EE_IsoBcut, EE_CRdPhicut, EE_masked])
        EE_OSSRmaskedcut = self.cutManager.AddList([EE_OScut, EE_IsoAcut, EE_IsoBcut, EE_SRdPhicut, EE_masked])

        self.dielectronRegions = [] # region name : region cut
        self.dielectronRegions.append(['BaseLine', '1'])
        self.dielectronRegions.append(['SS0', EE_SS0cut])
        self.dielectronRegions.append(['OS0', EE_OS0cut])
        self.dielectronRegions.append(['SSI', EE_SSIcut])
        self.dielectronRegions.append(['OSI', EE_OSIcut])
        self.dielectronRegions.append(['SSII', EE_SSIIcut])
        self.dielectronRegions.append(['OSII', EE_OSIIcut])
        self.dielectronRegions.append(['CROS', EE_OSCRcut])
        self.dielectronRegions.append(['CRSS', EE_SSCRcut])
        self.dielectronRegions.append(['SROS', EE_OSSRcut])
        self.dielectronRegions.append(['SRI', EE_SRIcut])
        self.dielectronRegions.append(['CRI', EE_CRIcut])
        self.dielectronRegions.append(['SRII', EE_SRIIcut])
        self.dielectronRegions.append(['CRII', EE_CRIIcut])
        self.dielectronRegions.append(['OS0disp', EE_OS0Dispcut])



        #### ---------------------------
        #### ---- MM Region definition
        #### ---------------------------

        MM_OScut = 'ev.DGM_charge[ev.DMDM_idxA[mm_maxIxy]]*ev.DGM_charge[ev.DMDM_idxB[mm_maxIxy]] < 0'
        MM_SScut = 'ev.DGM_charge[ev.DMDM_idxA[mm_maxIxy]]*ev.DGM_charge[ev.DMDM_idxB[mm_maxIxy]] > 0'
        MM_IsoAcut = '(ev.DGM_relPFiso[ev.DMDM_idxA[mm_maxIxy]] < 0.2)'
        MM_IsoBcut = '(ev.DGM_relPFiso[ev.DMDM_idxB[mm_maxIxy]] < 0.2)'
        MM_SRdPhicut = '(abs(ev.DMDM_dPhi[mm_maxIxy]) < 3.14/2)'
        MM_CRdPhicut = '(abs(ev.DMDM_dPhi[mm_maxIxy]) > 3.14/2)'
        MM_disp = '(abs(ev.DMDM_trackIxy_PV[mm_maxIxy]) > 4)'
        MM_SS0cut = self.cutManager.AddList([MM_SScut, MM_IsoAcut, MM_IsoBcut])
        MM_OS0cut = self.cutManager.AddList([MM_OScut, MM_IsoAcut, MM_IsoBcut])
        MM_SSIcut = self.cutManager.OR( self.cutManager.AddList([MM_SScut, 'not '+MM_IsoAcut, MM_IsoBcut]), self.cutManager.AddList([MM_SScut, MM_IsoAcut, 'not '+MM_IsoBcut]) )
        MM_OSIcut = self.cutManager.OR( self.cutManager.AddList([MM_OScut, 'not '+MM_IsoAcut, MM_IsoBcut]), self.cutManager.AddList([MM_OScut, MM_IsoAcut, 'not '+MM_IsoBcut]) )
        MM_Icut = self.cutManager.OR(self.cutManager.AddList(['not '+MM_IsoAcut, MM_IsoBcut]), self.cutManager.AddList([MM_IsoAcut, 'not '+MM_IsoBcut]))
        MM_SSIIcut = self.cutManager.AddList([MM_SScut, 'not '+MM_IsoAcut, 'not '+MM_IsoBcut])
        MM_OSIIcut = self.cutManager.AddList([MM_OScut, 'not '+MM_IsoAcut, 'not '+MM_IsoBcut])
        MM_IIcut = self.cutManager.AddList(['not '+MM_IsoAcut, 'not '+MM_IsoBcut])
        MM_OSCRcut = self.cutManager.AddList([MM_OScut, MM_IsoAcut, MM_IsoBcut, MM_CRdPhicut])
        MM_SSCRcut = self.cutManager.AddList([MM_SScut, MM_IsoAcut, MM_IsoBcut, MM_CRdPhicut])
        MM_OSSRcut = self.cutManager.AddList([MM_OScut, MM_IsoAcut, MM_IsoBcut, MM_SRdPhicut])
        MM_SRIcut = self.cutManager.AddList([MM_SRdPhicut, MM_Icut])
        MM_CRIcut = self.cutManager.AddList([MM_CRdPhicut, MM_Icut])
        MM_SRIIcut = self.cutManager.AddList([MM_SRdPhicut, MM_IIcut])
        MM_CRIIcut = self.cutManager.AddList([MM_CRdPhicut, MM_IIcut])
        MM_OS0Dispcut = self.cutManager.AddList([MM_OScut, MM_IsoAcut, MM_IsoBcut, MM_disp])

        self.dimuonRegions = [] # region name : region cut
        self.dimuonRegions.append(['BaseLine', '1'])
        self.dimuonRegions.append(['SS0', MM_SS0cut])
        self.dimuonRegions.append(['OS0', MM_OS0cut])
        self.dimuonRegions.append(['SSI', MM_SSIcut])
        self.dimuonRegions.append(['OSI', MM_OSIcut])
        self.dimuonRegions.append(['SSII', MM_SSIIcut])
        self.dimuonRegions.append(['OSII', MM_OSIIcut])
        self.dimuonRegions.append(['CROS', MM_OSCRcut])
        self.dimuonRegions.append(['CRSS', MM_SSCRcut])
        self.dimuonRegions.append(['SROS', MM_OSSRcut])
        self.dimuonRegions.append(['SRI', MM_SRIcut])
        self.dimuonRegions.append(['CRI', MM_CRIcut])
        self.dimuonRegions.append(['SRII', MM_SRIIcut])
        self.dimuonRegions.append(['CRII', MM_CRIIcut])
        self.dimuonRegions.append(['OS0disp', MM_OS0Dispcut])


        ## Init histogram efficiencies
        self.declareEfficiencies()


        #### -------------------------------
        #### ---- Histogram initialization
        #### -------------------------------
        for region in self.dielectronRegions:
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
        #nTracks = ev.RefittedPV_nPFTrack + ev.RefittedPV_nLostTrack + ev.RefittedPV_nExcludedTrack
        #if nTracks < 4: return

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
        if not ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 or ev.nDMDM == 0: return

        ### -> Find the MM pair with maximum Ixy that passes the cuts
        mm_maxIxy = -99
        maxIxy = -1
        nGoodMM = 0
        for i in range(0, ev.nDMDM):
            # Displaced Global ID cuts:
            if abs(ev.DGM_eta[ev.DMDM_idxA[i]]) > 2: continue
            if abs(ev.DGM_eta[ev.DMDM_idxB[i]]) > 2: continue
            if ev.DGM_pt[ev.DMDM_idxA[i]] < 31: continue
            if ev.DGM_pt[ev.DMDM_idxB[i]] < 31: continue
            if ev.DGM_ptError[ev.DMDM_idxA[i]]/ev.DGM_pt[ev.DMDM_idxA[i]] > 0.3: continue
            if ev.DGM_ptError[ev.DMDM_idxB[i]]/ev.DGM_pt[ev.DMDM_idxB[i]] > 0.3: continue
            if ev.DGM_ndof[ev.DMDM_idxA[i]] < 0.00001: continue
            if ev.DGM_ndof[ev.DMDM_idxB[i]] < 0.00001: continue
            if ev.DGM_chi2[ev.DMDM_idxA[i]]/ev.DGM_ndof[ev.DMDM_idxA[i]] > 10: continue
            if ev.DGM_chi2[ev.DMDM_idxB[i]]/ev.DGM_ndof[ev.DMDM_idxB[i]] > 10: continue
            if ev.DGM_numberOfValidHits[ev.DMDM_idxA[i]] < 22: continue
            if ev.DGM_numberOfValidHits[ev.DMDM_idxB[i]] < 22: continue
            # Selection cuts
            if ev.DMDM_cosAlpha[i] < -0.80: continue
            if ev.DMDM_normalizedChi2[i] > 10: continue
            if ev.DMDM_mass[i] < 15: continue

            l1 = TVector3()
            l2 = TVector3()
            l1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[i]], ev.DGM_eta[ev.DMDM_idxA[i]], ev.DGM_phi[ev.DMDM_idxA[i]])
            l2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[i]], ev.DGM_eta[ev.DMDM_idxB[i]], ev.DGM_phi[ev.DMDM_idxB[i]])
            if abs(l1.DeltaR(l2)) < 0.2: continue


            nGoodMM+=1

            if ev.DMDM_trackIxy_PV[i] > maxIxy:
                maxIxy = ev.DMDM_trackIxy_PV[i]
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

        l1 = TVector3()
        l2 = TVector3()
        l1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[mm_maxIxy]], ev.DGM_eta[ev.DMDM_idxA[mm_maxIxy]], ev.DGM_phi[ev.DMDM_idxA[mm_maxIxy]])
        l2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[mm_maxIxy]], ev.DGM_eta[ev.DMDM_idxB[mm_maxIxy]], ev.DGM_phi[ev.DMDM_idxB[mm_maxIxy]])
        deltaPhi = abs(l1.DeltaPhi(l2))

        exec("self.hMM{0}_nGoodMM.Fill(nGoodMM, weight)".format(region))
        exec("self.hMM{0}_dPhi.Fill(abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_mass.Fill(ev.DMDM_mass[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackIxy.Fill(ev.DMDM_trackIxy_PV[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackDxy.Fill(ev.DMDM_trackDxy_PV[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_Lxy.Fill(abs(ev.DMDM_Lxy_PV[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_Ixy.Fill(abs(ev.DMDM_Ixy_PV[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_trackDxy_BS.Fill(ev.DMDM_trackDxy_BS[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_trackIxy_BS.Fill(ev.DMDM_trackIxy_BS[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_Lxy_BS.Fill(abs(ev.DMDM_Lxy_BS[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_Ixy_BS.Fill(abs(ev.DMDM_Ixy_BS[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_cosAlpha.Fill(ev.DMDM_cosAlpha[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_leadingPt.Fill(ev.DMDM_leadingPt[mm_maxIxy], weight)".format(region))
        exec("self.hMM{0}_cosAlpha_dPhi.Fill(ev.DMDM_cosAlpha[mm_maxIxy], abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))
        exec("self.hMM{0}_DeltaPhi_dPhi.Fill(deltaPhi, abs(ev.DMDM_dPhi[mm_maxIxy]), weight)".format(region))



    #### ----------------------------
    #### --
    #### ---- DiElectron Processing
    #### --
    #### ----------------------------

    def processDielectrons(self, ev, weight):

        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 or not ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 or ev.nEE == 0: return

        ### -> Find the EE pair with maximum Ixy that passes the cuts
        ee_maxIxy = -99
        maxIxy = -1
        nGoodEE = 0
        for i in range(0, ev.nEE):
            if abs(ev.ElectronCandidate_eta[ev.EE_idxA[i]]) > 2.0: continue
            if abs(ev.ElectronCandidate_eta[ev.EE_idxB[i]]) > 2.0: continue
            if ev.EE_leadingPt[i] < 45: continue
            if ev.EE_leadingEt[i] < 45: continue
            if ev.EE_subleadingPt[i] < 28: continue
            if ev.EE_subleadingEt[i] < 28: continue
            if ev.EE_mass[i] < 15: continue
            #if ev.EE_relisoA[i] < 0.2: continue
            #if ev.EE_relisoB[i] < 0.2: continue
            if ev.EE_normalizedChi2[i] > 10: continue

            nGoodEE+=1

            if ev.EE_trackIxy_PV[i] > maxIxy:
                maxIxy = ev.EE_trackIxy_PV[i]
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

        l1 = TVector3()
        l2 = TVector3()
        l1.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxA[ee_maxIxy]], ev.ElectronCandidate_eta[ev.EE_idxA[ee_maxIxy]], ev.ElectronCandidate_phi[ev.EE_idxA[ee_maxIxy]])
        l2.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxB[ee_maxIxy]], ev.ElectronCandidate_eta[ev.EE_idxB[ee_maxIxy]], ev.ElectronCandidate_phi[ev.EE_idxB[ee_maxIxy]])
        deltaPhi = abs(l1.DeltaPhi(l2))


        exec("self.hEE{0}_nGoodEE.Fill(nGoodEE, weight)".format(region))
        exec("self.hEE{0}_dPhi.Fill(abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_mass.Fill(ev.EE_mass[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackDxy.Fill(ev.EE_trackDxy_PV[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackIxy.Fill(ev.EE_trackIxy_PV[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_Lxy.Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_Ixy.Fill(abs(ev.EE_Ixy_PV[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_trackDxy_BS.Fill(ev.EE_trackDxy_BS[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_trackIxy_BS.Fill(ev.EE_trackIxy_BS[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_Lxy_BS.Fill(abs(ev.EE_Lxy_BS[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_Ixy_BS.Fill(abs(ev.EE_Ixy_BS[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_leadingPt.Fill(ev.EE_leadingPt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_subleadingPt.Fill(ev.EE_subleadingPt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_leadingEt.Fill(ev.EE_leadingEt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_subleadingEt.Fill(ev.EE_subleadingEt[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_normalizedChi2.Fill(ev.EE_normalizedChi2[ee_maxIxy], weight)".format(region))
        exec("self.hEE{0}_cosAlpha_dPhi.Fill(ev.EE_cosAlpha[ee_maxIxy], abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_Lxy_dPhi.Fill(abs(ev.EE_Lxy_PV[ee_maxIxy]), abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))
        exec("self.hEE{0}_DeltaPhi_dPhi.Fill(deltaPhi, abs(ev.EE_dPhi[ee_maxIxy]), weight)".format(region))


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
            pass_massmin = False
            pass_relisomax = False
            pass_normChi2max = False
            pass_oscharge = False
            pass_dPhiSR = False
            pass_Ixymin = False

            for i in range(0, ev.nEE):

                if ev.EE_leadingPt[i] < 45: continue
                if ev.EE_subleadingPt[i] < 28: continue
                pass_ptmin = True

                if ev.EE_leadingEt[i] < 45: continue
                if ev.EE_subleadingEt[i] < 28: continue
                pass_etmin = True

                if abs(ev.ElectronCandidate_eta[ev.EE_idxA[i]]) > 2.0: continue
                if abs(ev.ElectronCandidate_eta[ev.EE_idxB[i]]) > 2.0: continue
                pass_etamax = True

                if ev.EE_mass[i] < 15: continue
                pass_massmin = True

                if (ev.IsoTrackSel_pfIsolationDR03[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[i]]]/ev.IsoTrackSel_pt[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[i]]] > 0.2): continue
                if (ev.IsoTrackSel_pfIsolationDR03[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[i]]]/ev.IsoTrackSel_pt[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[i]]] > 0.2): continue
                pass_relisomax = True

                if ev.EE_normalizedChi2[i] > 10: continue
                pass_normChi2max = True 

                if ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[i]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[i]]] > 0: continue
                pass_oscharge = True

                if abs(ev.EE_dPhi[i]) > 3.14/2.0: continue
                pass_dPhiSR = True

                if abs(ev.EE_trackIxy_PV[i]) < 6: continue
                pass_Ixymin = True

            if pass_ptmin: self.hEE_cutEfficiency.Fill(3, weight)
            if pass_etmin: self.hEE_cutEfficiency.Fill(4, weight)
            if pass_etamax: self.hEE_cutEfficiency.Fill(5, weight)
            if pass_massmin: self.hEE_cutEfficiency.Fill(6, weight)
            if pass_relisomax: self.hEE_cutEfficiency.Fill(7, weight)
            if pass_normChi2max: self.hEE_cutEfficiency.Fill(8, weight)
            if pass_oscharge: self.hEE_cutEfficiency.Fill(9, weight)
            if pass_dPhiSR: self.hEE_cutEfficiency.Fill(10, weight)
            if pass_Ixymin: self.hEE_cutEfficiency.Fill(11, weight)


        ### MM Channel
        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10:

            self.hMM_cutEfficiency.Fill(0, weight) 

            if not ev.PV_passAcceptance: return
            self.hMM_cutEfficiency.Fill(1, weight)

            if ev.nDMDM < 1: return
            self.hMM_cutEfficiency.Fill(2, weight)

            pass_ptmin = False
            pass_etamax = False
            pass_dRmin = False
            pass_massmin = False
            pass_relisomax = False
            pass_cosAlphamin = False
            pass_normChi2max = False
            pass_oscharge = False
            pass_dPhiSR = False
            pass_Ixymin = False

            for i in range(0, ev.nDMDM):

                if ev.DMDM_leadingPt[i] < 31: continue
                if ev.DMDM_subleadingPt[i] < 31: continue
                pass_ptmin = True

                if abs(ev.DGM_eta[ev.DMDM_idxA[i]]) > 2: continue
                if abs(ev.DGM_eta[ev.DMDM_idxB[i]]) > 2: continue
                pass_etamax = True

                l1 = TVector3()
                l2 = TVector3()
                l1.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[i]], ev.DGM_eta[ev.DMDM_idxA[i]], ev.DGM_phi[ev.DMDM_idxA[i]])
                l2.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[i]], ev.DGM_eta[ev.DMDM_idxB[i]], ev.DGM_phi[ev.DMDM_idxB[i]])
                if abs(l1.DeltaR(l2)) < 0.2: continue
                pass_dRmin = True

                if ev.DMDM_mass[i] < 15: continue
                pass_massmin = True

                if ev.DGM_relPFiso[ev.DMDM_idxA[i]] > 0.2: continue
                if ev.DGM_relPFiso[ev.DMDM_idxB[i]] > 0.2: continue
                pass_relisomax = True

                if ev.DMDM_cosAlpha[i] < -0.80: continue
                pass_cosAlphamin = True

                if ev.DMDM_normalizedChi2[i] > 10: continue
                pass_normChi2max = True 

                if ev.DGM_charge[ev.DMDM_idxA[i]]*ev.DGM_charge[ev.DMDM_idxB[i]] > 0: continue
                pass_oscharge = True

                if abs(ev.DMDM_dPhi[i]) > 3.14/2.0: continue
                pass_dPhiSR = True

                if abs(ev.DMDM_trackIxy_PV[i]) < 6: continue
                pass_Ixymin = True

            if pass_ptmin: self.hMM_cutEfficiency.Fill(3, weight)
            if pass_etamax: self.hMM_cutEfficiency.Fill(4, weight)
            if pass_dRmin: self.hMM_cutEfficiency.Fill(5, weight)
            if pass_massmin: self.hMM_cutEfficiency.Fill(6, weight)
            if pass_relisomax: self.hMM_cutEfficiency.Fill(7, weight)
            if pass_cosAlphamin: self.hMM_cutEfficiency.Fill(8, weight)
            if pass_normChi2max: self.hMM_cutEfficiency.Fill(9, weight)
            if pass_oscharge: self.hMM_cutEfficiency.Fill(10, weight)
            if pass_dPhiSR: self.hMM_cutEfficiency.Fill(11, weight)
            if pass_Ixymin: self.hMM_cutEfficiency.Fill(12, weight)




    #### --------------------------
    #### --
    #### ---- declare Efficiencies
    #### --
    #### --------------------------

    def declareEfficiencies(self):

        self.hEE_cutEfficiency = r.TH1F('hEE_cutEfficiency' + self.sufix, ';; Events with EE / Category', 12, 0, 12)
        self.hMM_cutEfficiency = r.TH1F('hMM_cutEfficiency' + self.sufix, ';; Events with MM / Category', 13, 0, 13)



    #### ------------------------
    #### --
    #### ---- declareHistograms
    #### --
    #### ------------------------

    def declareHistograms(self, region):

        #### ---------------
        #### ---- MM plots
        #### ---------------
        exec("self.hMM{0}_nGoodMM = r.TH1F('hMM{0}_nGoodMM' + self.sufix, ';Number of MM candidates;', 4, 0, 4)".format(region))
        exec("self.hMM{0}_nPU = r.TH1F('hMM{0}_nPU' + self.sufix, ';Number of true primary vertices;', 40, 0, 80)".format(region))
        exec("self.hMM{0}_dPhi = r.TH1F('hMM{0}_dPhi' + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 30, 0, 3.14)".format(region))
        exec("self.hMM{0}_mass = r.TH1F('hMM{0}_mass' + self.sufix, ';Dimuon invariant mass m_{{#mu#mu}} (GeV);', 55, 0, 500)".format(region))
        exec("self.hMM{0}_cosAlpha = r.TH1F('hMM{0}_cosAlpha' + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}});', 22, -1.1, 1.1)".format(region))
        exec("self.hMM{0}_trackIxy = r.TH1F('hMM{0}_trackIxy' + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_trackDxy = r.TH1F('hMM{0}_trackDxy' + self.sufix, ';Dimuon |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hMM{0}_trackIxy_BS = r.TH1F('hMM{0}_trackIxy_BS' + self.sufix, ';Dimuon |d_{{0}}|/#sigma_{{d}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_trackDxy_BS = r.TH1F('hMM{0}_trackDxy_BS' + self.sufix, ';Dimuon |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hMM{0}_Lxy = r.TH1F('hMM{0}_Lxy' + self.sufix, ';Dimuon vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hMM{0}_Ixy = r.TH1F('hMM{0}_Ixy' + self.sufix, ';Dimuon vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_Lxy_BS = r.TH1F('hMM{0}_Lxy_BS' + self.sufix, ';Dimuon vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hMM{0}_Ixy_BS = r.TH1F('hMM{0}_Ixy_BS' + self.sufix, ';Dimuon vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hMM{0}_leadingPt = r.TH1F('hMM{0}_leadingPt' + self.sufix, ';Dimuon leading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hMM{0}_subleadingPt = r.TH1F('hMM{0}_subleadingPt' + self.sufix, ';Dimuon subleading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hMM{0}_normalizedChi2 = r.TH1F('hMM{0}_normalizedChi2' + self.sufix, ';Dimuon vertex fit #chi^{{2}}/ndof;', 20, 0, 10)".format(region))
        exec("self.hMM{0}_cosAlpha_dPhi = r.TH2F('hMM{0}_cosAlpha_dPhi' + self.sufix, ';Dimuon cos(#alpha_{{#mu#mu}}) ; Dimuon collinearity |#Delta#Phi|', 20, -1.0, 1.0, 20, 0, 3.14)".format(region))
        exec("self.hMM{0}_DeltaPhi_dPhi = r.TH2F('hMM{0}_DeltaPhi_dPhi' + self.sufix, ';Dimuon #Delta#Phi(ll) ; Dimuon collinearity |#Delta#Phi|', 20, 0.0, 3.14, 20, 0, 3.14)".format(region))


        #### ---------------
        #### ---- EE plots
        #### ---------------

        exec("self.hEE{0}_nGoodEE = r.TH1F('hEE{0}_nGoodEE' + self.sufix, ';Number of EE candidates;', 4, 0, 4)".format(region))
        exec("self.hEE{0}_nPU = r.TH1F('hEE{0}_nPU' + self.sufix, ';Number of true primary vertices;', 40, 0, 80)".format(region))
        exec("self.hEE{0}_dPhi = r.TH1F('hEE{0}_dPhi' + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 30, 0, 3.14)".format(region))
        exec("self.hEE{0}_mass = r.TH1F('hEE{0}_mass' + self.sufix, ';Dielectron invariant mass m_{{ee}} (GeV);', 55, 0, 500)".format(region))
        exec("self.hEE{0}_cosAlpha = r.TH1F('hEE{0}_cosAlpha' + self.sufix, ';Dielectron cos(#alpha_{{ee}});', 22, -1.1, 1.1)".format(region))
        exec("self.hEE{0}_trackIxy = r.TH1F('hEE{0}_trackIxy' + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_trackDxy = r.TH1F('hEE{0}_trackDxy' + self.sufix, ';Dielectron |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hEE{0}_trackIxy_BS = r.TH1F('hEE{0}_trackIxy_BS' + self.sufix, ';Dielectron |d_{{0}}|/#sigma_{{d}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_trackDxy_BS = r.TH1F('hEE{0}_trackDxy_BS' + self.sufix, ';Dielectron |d_{{0}}| (cm);', 30, 0, 0.5)".format(region))
        exec("self.hEE{0}_Lxy = r.TH1F('hEE{0}_Lxy' + self.sufix, ';Dielectron vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hEE{0}_Ixy = r.TH1F('hEE{0}_Ixy' + self.sufix, ';Dielectron vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_Lxy_BS = r.TH1F('hEE{0}_Lxy_BS' + self.sufix, ';Dielectron vertex |L_{{xy}}| (cm);', 20, 0, 10)".format(region))
        exec("self.hEE{0}_Ixy_BS = r.TH1F('hEE{0}_Ixy_BS' + self.sufix, ';Dielectron vertex |L_{{xy}}|/#sigma_{{L}};', 20, 0, 20)".format(region))
        exec("self.hEE{0}_leadingPt = r.TH1F('hEE{0}_leadingPt' + self.sufix, ';Dielectron leading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_subleadingPt = r.TH1F('hEE{0}_subleadingPt' + self.sufix, ';Dielectron subleading p_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_leadingEt = r.TH1F('hEE{0}_leadingEt' + self.sufix, ';Dielectron leading E_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_subleadingEt = r.TH1F('hEE{0}_subleadingEt' + self.sufix, ';Dielectron subleading E_{{T}};', 30, 0, 300)".format(region))
        exec("self.hEE{0}_normalizedChi2 = r.TH1F('hEE{0}_normalizedChi2' + self.sufix, ';Dielectron vertex fit #chi^{{2}}/ndof;', 20, 0, 10)".format(region))
        exec("self.hEE{0}_cosAlpha_dPhi = r.TH2F('hEE{0}_cosAlpha_dPhi' + self.sufix, ';Dielectron cos(#alpha_{{ee}}) ; Dielectron collinearity |#Delta#Phi|', 20, -1.0, 1.0, 20, 0, 3.14)".format(region))
        exec("self.hEE{0}_Lxy_dPhi = r.TH2F('hEE{0}_Lxy_dPhi' + self.sufix, ';Dielectron |L_{{xy}}| (cm) ; Dielectron collinearity |#Delta#Phi|', 20, 0.0, 1.0, 20, 0, 3.14)".format(region))
        exec("self.hEE{0}_DeltaPhi_dPhi = r.TH2F('hEE{0}_DeltaPhi_dPhi' + self.sufix, ';Dielectron #Delta#Phi(ll) ; Dielectron collinearity |#Delta#Phi|', 20, 0.0, 3.14, 20, 0, 3.14)".format(region))



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

        if self.hEE_cutEfficiency.GetEntries() > 0:
            self.hEE_cutEfficiency.Write()
            self.hMM_cutEfficiency.Write()
        else:
            for attr, value in self.__dict__.iteritems():
                if attr[0] == 'h': value.Write()

        
        output.Close()


