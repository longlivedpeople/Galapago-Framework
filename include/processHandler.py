import math
import os
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager



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

        self.hcounts = r.TH1F('hcounts'+self.sufix, '', 1, 0, 1)

        ###########################
        ###  DiMuon Histograms  ###
        ###########################
        # -> Dimuon candidate
        self.hMM_dPhi = r.TH1F('hMM_dPhi' + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 20, -3.3, 3.3)
        self.hMM_mass = r.TH1F('hMM_mass' + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', 35, 0, 200)
        self.hMM_cosAlpha = r.TH1F('hMM_cosAlpha' + self.sufix, ';Dimuon cos(#alpha_{#mu#mu});', 21, -1.1, 1.1) 
        self.hMM_trackIxy = r.TH1F('hMM_trackIxy' + self.sufix, ';Dimuon |d_{0}|/#sigma_{d};', 20, 0, 20)
        self.hMM_trackDxy = r.TH1F('hMM_trackDxy' + self.sufix, ';Dimuon |d_{0}| (cm);', 20, 0, 5)
        self.hMM_Lxy = r.TH1F('hMM_Lxy' + self.sufix, ';Dimuon vertex |L_{xy}| (cm);', 20, 0, 10)
        self.hMM_Ixy = r.TH1F('hMM_Ixy' + self.sufix, ';Dimuon vertex |L_{xy}|/#sigma_{L};', 20, 0, 20)
        self.hMM_leadingPt = r.TH1F('hMM_leadingPt' + self.sufix, ';Dimuon leading p_{T};', 30, 0, 300)
        self.hMM_subleadingPt = r.TH1F('hMM_subleadingPt' + self.sufix, ';Dimuon subleading p_{T};', 30, 0, 300)
        self.hMM_normalizedChi2 = r.TH1F('hMM_normalizedChi2' + self.sufix, ';Dimuon vertex fit #chi^{2}/ndof;', 20, 0, 10)
        # -> Dimuon members
        self.hM_eta = r.TH1F('hM_eta' + self.sufix, ';DG #eta;', 27, -2.6, 2.6)
        self.hM_dxy = r.TH1F('hM_dxy' + self.sufix, ';DG transverse impact parameter d_{xy};', 60, -3.0, 3.0)
        self.hM_dxyError = r.TH1F('hM_dxyError' + self.sufix, ';DG #sigma_{d};', 70, 0.0, 0.01)
        self.hM_charge = r.TH1F('hM_charge' + self.sufix, ';DG charge q;', 5, -2.5, 2.5)
        self.hM_numberOfValidHits = r.TH1F('hM_numberOfValidHits' + self.sufix, ';DG N_{Hits};', 80, 0.0, 80.0)
        self.hM_normChi2 = r.TH1F('hM_normChi2' + self.sufix, ';DG #chi^{2}/ndof;', 30, 0.0, 15.0)
        


        ###############################
        ###  DiElectron Histograms  ###
        ###############################
        self.hEE_dPhi = r.TH1F('hEE_dPhi' + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 20, -3.3, 3.3)
        self.hEE_mass = r.TH1F('hEE_mass' + self.sufix, ';Dielectron invariant mass m_{ee} (GeV);', 35, 0, 200)
        self.hEE_trackIxy = r.TH1F('hEE_trackIxy' + self.sufix, ';Dielectron |d_{0}|/#sigma_{d};', 20, 0, 20)
        self.hEE_trackDxy = r.TH1F('hEE_trackDxy' + self.sufix, ';Dielectron |d_{0}| (cm);', 20, 0, 5)
        self.hEE_Lxy = r.TH1F('hEE_Lxy' + self.sufix, ';Dielectron vertex |L_{xy}| (cm);', 20, 0, 10)
        self.hEE_Ixy = r.TH1F('hEE_Ixy' + self.sufix, ';Dielectron vertex |L_{xy}|/#sigma_{L};', 20, 0, 20)
        self.hEE_leadingPt = r.TH1F('hEE_leadingPt' + self.sufix, ';Dielectron leading p_{T};', 30, 0, 300)
        self.hEE_subleadingPt = r.TH1F('hEE_subleadingPt' + self.sufix, ';Dielectron subleading p_{T};', 30, 0, 300)
        self.hEE_normalizedChi2 = r.TH1F('hEE_normalizedChi2' + self.sufix, ';Dielectron vertex fit #chi^{2}/ndof;', 30, 0, 15)
        # -> Dielectron members
        self.hE_eta = r.TH1F('hE_eta' + self.sufix, ';Electron #eta;', 27, -2.6, 2.6)
        self.hE_dxy = r.TH1F('hE_dxy' + self.sufix, ';Electron transverse impact parameter d_{xy};', 60, -3.0, 3.0)
        self.hE_dxyError = r.TH1F('hE_dxyError' + self.sufix, ';Electron #sigma_{d};', 70, 0.0, 0.01)
        self.hE_charge = r.TH1F('hE_charge' + self.sufix, ';Electron charge q;', 5, -2.5, 2.5)

        #######################
        ###  SS Histograms  ### (Data only)
        #######################
        if self.isdata:

            self.hMM_dPhi_SS = r.TH1F('hMM_dPhi_SS' + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 20, -3.3, 3.3)
            self.hMM_mass_SS = r.TH1F('hMM_mass_SS' + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', 35, 0, 200)
            self.hMM_cosAlpha_SS = r.TH1F('hMM_cosAlpha_SS' + self.sufix, ';Dimuon cos(#alpha_{#mu#mu});', 21, -1.1, 1.1) 
            self.hMM_trackIxy_SS = r.TH1F('hMM_trackIxy_SS' + self.sufix, ';Dimuon |d_{0}|/#sigma_{d};', 20, 0, 20)
            self.hMM_trackDxy_SS = r.TH1F('hMM_trackDxy_SS' + self.sufix, ';Dimuon |d_{0}| (cm);', 20, 0, 5)
            self.hMM_Lxy_SS = r.TH1F('hMM_Lxy_SS' + self.sufix, ';Dimuon vertex |L_{xy}| (cm);', 20, 0, 10)
            self.hMM_Ixy_SS = r.TH1F('hMM_Ixy_SS' + self.sufix, ';Dimuon vertex |L_{xy}|/#sigma_{L};', 20, 0, 20)
            self.hMM_leadingPt_SS = r.TH1F('hMM_leadingPt_SS' + self.sufix, ';Dimuon leading p_{T};', 30, 0, 300)
            self.hMM_subleadingPt_SS = r.TH1F('hMM_subleadingPt_SS' + self.sufix, ';Dimuon subleading p_{T};', 30, 0, 300)
            self.hMM_normalizedChi2_SS = r.TH1F('hMM_normalizedChi2_SS' + self.sufix, ';Dimuon vertex fit #chi^{2}/ndof;', 20, 0, 10)

            self.hM_eta_SS = r.TH1F('hM_eta_SS' + self.sufix, ';DG #eta;', 27, -2.6, 2.6)
            self.hM_dxy_SS = r.TH1F('hM_dxy_SS' + self.sufix, ';DG transverse impact parameter d_{xy};', 60, -3.0, 3.0)
            self.hM_dxyError_SS = r.TH1F('hM_dxyError_SS' + self.sufix, ';DG #sigma_{d};', 70, 0.0, 0.01)
            self.hM_charge_SS = r.TH1F('hM_charge_SS' + self.sufix, ';DG charge q;', 5, -2.5, 2.5)
            self.hM_numberOfValidHits_SS = r.TH1F('hM_numberOfValidHits_SS' + self.sufix, ';DG N_{Hits};', 80, 0.0, 80.0)
            self.hM_normChi2_SS = r.TH1F('hM_normChi2_SS' + self.sufix, ';DG #chi^{2}/ndof;', 30, 0.0, 15.0)

            self.hEE_dPhi_SS = r.TH1F('hEE_dPhi_SS' + self.sufix, ';Dielectron collinearity |#Delta#Phi|;', 20, -3.3, 3.3)
            self.hEE_mass_SS = r.TH1F('hEE_mass_SS' + self.sufix, ';Dielectron invariant mass m_{ee} (GeV);', 35, 0, 200)
            self.hEE_trackIxy_SS = r.TH1F('hEE_trackIxy_SS' + self.sufix, ';Dielectron |d_{0}|/#sigma_{d};', 20, 0, 20)
            self.hEE_trackDxy_SS = r.TH1F('hEE_trackDxy_SS' + self.sufix, ';Dielectron |d_{0}| (cm);', 20, 0, 5)
            self.hEE_Lxy_SS = r.TH1F('hEE_Lxy_SS' + self.sufix, ';Dielectron vertex |L_{xy}| (cm);', 20, 0, 10)
            self.hEE_Ixy_SS = r.TH1F('hEE_Ixy_SS' + self.sufix, ';Dielectron vertex |L_{xy}|/#sigma_{L};', 20, 0, 20)
            self.hEE_leadingPt_SS = r.TH1F('hEE_leadingPt_SS' + self.sufix, ';Dielectron leading p_{T};', 30, 0, 300)
            self.hEE_subleadingPt_SS = r.TH1F('hEE_subleadingPt_SS' + self.sufix, ';Dielectron subleading p_{T};', 30, 0, 300)
            self.hEE_normalizedChi2_SS = r.TH1F('hEE_normalizedChi2_SS' + self.sufix, ';Dielectron vertex fit #chi^{2}/ndof;', 30, 0, 15)

            self.hE_eta_SS = r.TH1F('hE_eta_SS' + self.sufix, ';Electron #eta;', 27, -2.6, 2.6)
            self.hE_dxy_SS = r.TH1F('hE_dxy_SS' + self.sufix, ';Electron transverse impact parameter d_{xy};', 60, -3.0, 3.0)
            self.hE_dxyError_SS = r.TH1F('hE_dxyError_SS' + self.sufix, ';Electron #sigma_{d};', 70, 0.0, 0.01)
            self.hE_charge_SS = r.TH1F('hE_charge_SS' + self.sufix, ';Electron charge q;', 5, -2.5, 2.5)

        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Sumw2()



    def processEvent(self, ev):

        if not self.isdata:
            weight = self.lumiweight*ev.wPU
        else: 
            weight = 1

        #### Generation variables

        self.hcounts.Fill(0, weight)


        try:
            self.processDimuons(ev, weight)
        except AttributeError:
            print('There is some collections missed in this file: Dimuon histograms will be empty')

        try:
            self.processDielectrons(ev, weight)
        except AttributeError:
            print('There is some collections missed in this file: Dielectron histograms will be empty')


    def processDimuons(self, ev, weight):

        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 and ev.nDMDMBase > 0:

            if eval(self.cutManager.LoopMM_OScharge) and eval(self.cutManager.LoopMMCR_dPhi):

                self.hMM_dPhi.Fill(ev.DMDMBase_dPhi[ev.DMDMBase_maxIxy], weight)
                self.hMM_mass.Fill(ev.DMDMBase_mass[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackIxy.Fill(ev.DMDMBase_trackIxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackDxy.Fill(ev.DMDMBase_trackDxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Lxy.Fill(ev.DMDMBase_Lxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Ixy.Fill(ev.DMDMBase_Ixy[ev.DMDMBase_maxIxy], weight)
                self.hMM_cosAlpha.Fill(ev.DMDMBase_cosAlpha[ev.DMDMBase_maxIxy], weight)
                self.hMM_leadingPt.Fill(ev.DMDMBase_leadingPt[ev.DMDMBase_maxIxy], weight)
                self.hMM_subleadingPt.Fill(ev.DMDMBase_subleadingPt[ev.DMDMBase_maxIxy], weight)
                self.hMM_normalizedChi2.Fill(ev.DMDMBase_normalizedChi2[ev.DMDMBase_maxIxy], weight)

                self.hM_eta.Fill(ev.DGM_eta[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_eta.Fill(ev.DGM_eta[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxy.Fill(ev.DGM_dxy[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxy.Fill(ev.DGM_dxy[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxyError.Fill(ev.DGM_dxyError[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxyError.Fill(ev.DGM_dxyError[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_charge.Fill(ev.DGM_charge[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_charge.Fill(ev.DGM_charge[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_numberOfValidHits.Fill(ev.DGM_numberOfValidHits[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_numberOfValidHits.Fill(ev.DGM_numberOfValidHits[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_normChi2.Fill(ev.DGM_chi2[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]]/ev.DGM_ndof[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_normChi2.Fill(ev.DGM_chi2[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]]/ev.DGM_ndof[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)


            if self.isdata and eval(self.cutManager.LoopMM_SScharge) and eval(self.cutManager.LoopMMCR_dPhi):

                self.hMM_dPhi_SS.Fill(ev.DMDMBase_dPhi[ev.DMDMBase_maxIxy], weight)
                self.hMM_mass_SS.Fill(ev.DMDMBase_mass[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackIxy_SS.Fill(ev.DMDMBase_trackIxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackDxy_SS.Fill(ev.DMDMBase_trackDxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Lxy_SS.Fill(ev.DMDMBase_Lxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Ixy_SS.Fill(ev.DMDMBase_Ixy[ev.DMDMBase_maxIxy], weight)
                self.hMM_cosAlpha_SS.Fill(ev.DMDMBase_cosAlpha[ev.DMDMBase_maxIxy], weight)
                self.hMM_leadingPt_SS.Fill(ev.DMDMBase_leadingPt[ev.DMDMBase_maxIxy], weight)
                self.hMM_subleadingPt_SS.Fill(ev.DMDMBase_subleadingPt[ev.DMDMBase_maxIxy], weight)
                self.hMM_normalizedChi2_SS.Fill(ev.DMDMBase_normalizedChi2[ev.DMDMBase_maxIxy], weight)

                self.hM_eta_SS.Fill(ev.DGM_eta[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_eta_SS.Fill(ev.DGM_eta[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxy_SS.Fill(ev.DGM_dxy[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxy_SS.Fill(ev.DGM_dxy[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxyError_SS.Fill(ev.DGM_dxyError[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_dxyError_SS.Fill(ev.DGM_dxyError[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_charge_SS.Fill(ev.DGM_charge[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_charge_SS.Fill(ev.DGM_charge[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_numberOfValidHits_SS.Fill(ev.DGM_numberOfValidHits[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_numberOfValidHits_SS.Fill(ev.DGM_numberOfValidHits[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)
                self.hM_normChi2_SS.Fill(ev.DGM_chi2[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]]/ev.DGM_ndof[ev.DMDMBase_idxA[ev.DMDMBase_maxIxy]], weight)
                self.hM_normChi2_SS.Fill(ev.DGM_chi2[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]]/ev.DGM_ndof[ev.DMDMBase_idxB[ev.DMDMBase_maxIxy]], weight)

    def processDielectrons(self, ev, weight):

        if not ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 and ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 and ev.nEEBase > 0:


            if eval(self.cutManager.LoopEE_OScharge) and eval(self.cutManager.LoopEECR_dPhi):

                self.hEE_dPhi.Fill(ev.EEBase_dPhi[ev.EEBase_maxIxy], weight)
                self.hEE_mass.Fill(ev.EEBase_mass[ev.EEBase_maxIxy], weight)
                self.hEE_trackIxy.Fill(ev.EEBase_trackIxy[ev.EEBase_maxIxy], weight)
                self.hEE_trackDxy.Fill(ev.EEBase_trackDxy[ev.EEBase_maxIxy], weight)
                self.hEE_Lxy.Fill(ev.EEBase_Lxy[ev.EEBase_maxIxy], weight)
                self.hEE_Ixy.Fill(ev.EEBase_Ixy[ev.EEBase_maxIxy], weight)
                self.hEE_leadingPt.Fill(ev.EEBase_leadingPt[ev.EEBase_maxIxy], weight)
                self.hEE_subleadingPt.Fill(ev.EEBase_subleadingPt[ev.EEBase_maxIxy], weight)
                self.hEE_normalizedChi2.Fill(ev.EEBase_normalizedChi2[ev.EEBase_maxIxy], weight)

                self.hE_eta.Fill(ev.ElectronCandidate_eta[ev.EEBase_idxA[ev.EEBase_maxIxy]])
                self.hE_eta.Fill(ev.ElectronCandidate_eta[ev.EEBase_idxB[ev.EEBase_maxIxy]])
                self.hE_dxy.Fill(ev.ElectronCandidate_dxy[ev.EEBase_idxA[ev.EEBase_maxIxy]])
                self.hE_dxy.Fill(ev.ElectronCandidate_dxy[ev.EEBase_idxB[ev.EEBase_maxIxy]])
                self.hE_dxyError.Fill(ev.ElectronCandidate_dxyError[ev.EEBase_idxA[ev.EEBase_maxIxy]])
                self.hE_dxyError.Fill(ev.ElectronCandidate_dxyError[ev.EEBase_idxB[ev.EEBase_maxIxy]])
                self.hE_charge.Fill(ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxA[ev.EEBase_maxIxy]]])
                self.hE_charge.Fill(ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxB[ev.EEBase_maxIxy]]])

            if self.isdata and eval(self.cutManager.LoopEE_SScharge) and eval(self.cutManager.LoopEECR_dPhi):

                self.hEE_dPhi_SS.Fill(ev.EEBase_dPhi[ev.EEBase_maxIxy], weight)
                self.hEE_mass_SS.Fill(ev.EEBase_mass[ev.EEBase_maxIxy], weight)
                self.hEE_trackIxy_SS.Fill(ev.EEBase_trackIxy[ev.EEBase_maxIxy], weight)
                self.hEE_trackDxy_SS.Fill(ev.EEBase_trackDxy[ev.EEBase_maxIxy], weight)
                self.hEE_Lxy_SS.Fill(ev.EEBase_Lxy[ev.EEBase_maxIxy], weight)
                self.hEE_Ixy_SS.Fill(ev.EEBase_Ixy[ev.EEBase_maxIxy], weight)
                self.hEE_leadingPt_SS.Fill(ev.EEBase_leadingPt[ev.EEBase_maxIxy], weight)
                self.hEE_subleadingPt_SS.Fill(ev.EEBase_subleadingPt[ev.EEBase_maxIxy], weight)
                self.hEE_normalizedChi2_SS.Fill(ev.EEBase_normalizedChi2[ev.EEBase_maxIxy], weight)

                self.hE_eta_SS.Fill(ev.ElectronCandidate_eta[ev.EEBase_idxA[ev.EEBase_maxIxy]])
                self.hE_eta_SS.Fill(ev.ElectronCandidate_eta[ev.EEBase_idxB[ev.EEBase_maxIxy]])
                self.hE_dxy_SS.Fill(ev.ElectronCandidate_dxy[ev.EEBase_idxA[ev.EEBase_maxIxy]])
                self.hE_dxy_SS.Fill(ev.ElectronCandidate_dxy[ev.EEBase_idxB[ev.EEBase_maxIxy]])
                self.hE_dxyError_SS.Fill(ev.ElectronCandidate_dxyError[ev.EEBase_idxA[ev.EEBase_maxIxy]])
                self.hE_dxyError_SS.Fill(ev.ElectronCandidate_dxyError[ev.EEBase_idxB[ev.EEBase_maxIxy]])
                self.hE_charge_SS.Fill(ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxA[ev.EEBase_maxIxy]]])
                self.hE_charge_SS.Fill(ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxB[ev.EEBase_maxIxy]]])


        """ CUT EXAMPLE
        if eval(self.cutManager.twoElectrons):
            self.hEE.Fill(ev.nEE)
        if eval(self.cutManager.twoMuons):
            self.hMM.Fill(ev.nMM)
        if eval(self.cutManager.Add(self.cutManager.twoElectrons, self.cutManager.twoMuons)):
            self.hLL.Fill(ev.nEE)
            self.hLL.Fill(ev.nMM)
        """



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


