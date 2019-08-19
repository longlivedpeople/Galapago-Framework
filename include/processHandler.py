import math
import ROOT as r
from ROOT import TVector3
import include.CutManager as CutManager



class processHandler:

    def __init__(self, fileName, treename, blockname, samplename):
        self.fileName = fileName
        self.treename = treename
        self.blockname = blockname
        self.samplename = samplename
        self.cutManager = CutManager.CutManager()


        # Generation histograms
        self.hgenMM_cosAlpha = r.TH1F('hgenMM_cosAlpha_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 22, -1.1, 1.1)

        # Single electron histograms
        self.hE1_pt = r.TH1F('hE1_pt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hE2_pt = r.TH1F('hE2_pt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hE1_et = r.TH1F('hE1_et_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hE2_et = r.TH1F('hE2_et_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hE1_eta = r.TH1F('hE1_eta_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, -3, 3)
        self.hE2_eta = r.TH1F('hE2_eta_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, -3, 3)
        self.hE1_reliso = r.TH1F('hE1_reliso_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 1)
        self.hE2_reliso = r.TH1F('hE2_reliso_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 1)

        # Single muon histograms
        self.hM1_pt = r.TH1F('hM1_pt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hM2_pt = r.TH1F('hM2_pt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hM1_triggerPt = r.TH1F('hM1_triggerPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hM2_triggerPt = r.TH1F('hM2_triggerPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 100)
        self.hM1_eta = r.TH1F('hM1_eta_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, -3, 3)
        self.hM2_eta = r.TH1F('hM2_eta_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, -3, 3)
        self.hM1_reliso = r.TH1F('hM1_reliso_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 1)
        self.hM2_reliso = r.TH1F('hM2_reliso_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 1)

        # Dilepton numbers 
        self.hnLL = r.TH1F('hnLL_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)
        self.hnEE = r.TH1F('hnEE_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)
        self.hnMM = r.TH1F('hnMM_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)

        # Dielectron selection (Signal Region SR)
        self.hSR_EEsel_minIxy = r.TH1F('hSR_EEsel_minIxy_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 20)
        self.hSR_EEsel_invMass = r.TH1F('hSR_EEsel_invMass_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 500)
        self.hSR_EEsel_Chi2 = r.TH1F('hSR_EEsel_Chi2_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 20, 0, 20)
        self.hSR_EEsel_dPhi = r.TH1F('hSR_EEsel_dPhi_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, math.pi)
        self.hSR_EEsel_leadingPt = r.TH1F('hSR_EEsel_leadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 1000)
        self.hSR_EEsel_subleadingPt = r.TH1F('hSR_EEsel_subleadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 700)
        
        # Dimuon selection (Signal Region SR)
        self.hSR_MMsel_minIxy = r.TH1F('hSR_MMsel_minIxy_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 20)
        self.hSR_MMsel_invMass = r.TH1F('hSR_MMsel_invMass_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 500)
        self.hSR_MMsel_Chi2 = r.TH1F('hSR_MMsel_Chi2_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 20, 0, 20)
        self.hSR_MMsel_dPhi = r.TH1F('hSR_MMsel_dPhi_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, math.pi)
        self.hSR_MMsel_leadingPt = r.TH1F('hSR_MMsel_leadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 1000)
        self.hSR_MMsel_subleadingPt = r.TH1F('hSR_MMsel_subleadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 700)
        self.hSR_MMsel_cosAlpha = r.TH1F('hSR_MMsel_cosAlpha_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 22, -1.1, 1.1)

        # Dielectron selection (Control Region CR)
        self.hCR_EEsel_minIxy = r.TH1F('hCR_EEsel_minIxy_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 20)
        self.hCR_EEsel_invMass = r.TH1F('hCR_EEsel_invMass_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 500)
        self.hCR_EEsel_Chi2 = r.TH1F('hCR_EEsel_Chi2_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 20, 0, 20)
        self.hCR_EEsel_dPhi = r.TH1F('hCR_EEsel_dPhi_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, math.pi)
        self.hCR_EEsel_leadingPt = r.TH1F('hCR_EEsel_leadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 1000)
        self.hCR_EEsel_subleadingPt = r.TH1F('hCR_EEsel_subleadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 700)
        
        # Dimuon selection (Control Region CR)
        self.hCR_MMsel_minIxy = r.TH1F('hCR_MMsel_minIxy_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 20)
        self.hCR_MMsel_invMass = r.TH1F('hCR_MMsel_invMass_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 500)
        self.hCR_MMsel_Chi2 = r.TH1F('hCR_MMsel_Chi2_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 20, 0, 20)
        self.hCR_MMsel_dPhi = r.TH1F('hCR_MMsel_dPhi_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, math.pi)
        self.hCR_MMsel_leadingPt = r.TH1F('hCR_MMsel_leadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 1000)
        self.hCR_MMsel_subleadingPt = r.TH1F('hCR_MMsel_subleadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 700)
        self.hCR_MMsel_cosAlpha = r.TH1F('hCR_MMsel_cosAlpha_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 22, -1.1, 1.1)


        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Sumw2()



    def processEvent(self, ev, weight):

        #### Generation variables

        if (ev.nGenLepton == 4):

            ## Pair indentification
            gPair = [[], []] 

            for i in range(0, ev.nGenLepton):
                gPair[ev.GenLeptonSel_motherIdx[i]].append(i)


            ## cos(alpha) studies
            for p in range(0, len(gPair)):
                if abs(ev.GenLeptonSel_pdgId[gPair[p][0]]) == 13:
                    m0 = TVector3()
                    m1 = TVector3()
                    m0.SetPtEtaPhi(ev.GenLeptonSel_pt[gPair[p][0]], ev.GenLeptonSel_eta[gPair[p][0]], ev.GenLeptonSel_phi[gPair[p][0]])
                    m1.SetPtEtaPhi(ev.GenLeptonSel_pt[gPair[p][1]], ev.GenLeptonSel_eta[gPair[p][1]], ev.GenLeptonSel_phi[gPair[p][1]])
 
                    self.hgenMM_cosAlpha.Fill(math.cos(m0.Angle(m1)), weight)



        #### Dielectron channel
        if ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8 and ev.nElectronCandidate > 1:            

            # Single electron histograms
            ei = list(range(0, ev.nElectronCandidate))
            ei = sorted(ei, reverse=True, key=lambda x: ev.ElectronCandidate_pt[x]) # sort electrons by p
            self.hE1_pt.Fill(ev.ElectronCandidate_pt[ei[0]])
            self.hE2_pt.Fill(ev.ElectronCandidate_pt[ei[1]])
            self.hE1_et.Fill(ev.ElectronCandidate_et[ei[0]])
            self.hE2_et.Fill(ev.ElectronCandidate_et[ei[1]])
            self.hE1_eta.Fill(ev.ElectronCandidate_eta[ei[0]])
            self.hE2_eta.Fill(ev.ElectronCandidate_eta[ei[1]])
            self.hE1_reliso.Fill(ev.IsoTrackSel_relPfIsolationDR03[ev.ElectronCandidate_isotrackIdx[ei[0]]])
            self.hE2_reliso.Fill(ev.IsoTrackSel_relPfIsolationDR03[ev.ElectronCandidate_isotrackIdx[ei[1]]])



            
            # Dielectron candidates

            self.hnEE.Fill(ev.nEE, weight)

            maxIxy = -1
            maxIxyi = -99

            for ee in range(0, ev.nEE):

                if abs(ev.ElectronCandidate_eta[ev.EE_idxA[ee]]) > 1.4442: continue
                if abs(ev.ElectronCandidate_eta[ev.EE_idxB[ee]]) > 1.4442: continue
                if ev.EE_relisoA[ee] > 0.1 or ev.EE_relisoB[ee] > 0.1: continue
                if ev.EE_normalizedChi2[ee] > 10: continue
                if ev.EE_leadingPt[ee] < 41: continue
                if ev.EE_subleadingPt[ee] < 24: continue
                if ev.EE_leadingEt[ee] < 45: continue
                if ev.EE_subleadingEt[ee] < 28: continue
                if ev.EE_invMass[ee] < 15: continue
                #if abs(ev.EE_dPhi[ee]) > math.pi/2.0: continue

                IxyA = ev.ElectronCandidate_dxySignificance[ev.EE_idxA[ee]]
                IxyB = ev.ElectronCandidate_dxySignificance[ev.EE_idxB[ee]]

                mIxy = IxyA if IxyA < IxyB else IxyB


                if mIxy > maxIxy:
                    maxIxy = mIxy
                    maxIxyi = ee

            if maxIxyi > -1: 

                if abs(ev.EE_dPhi[maxIxyi]) < math.pi/2.0:

                    self.hSR_EEsel_minIxy.Fill(maxIxy, weight)
                    self.hSR_EEsel_invMass.Fill(ev.EE_invMass[maxIxyi], weight)
                    self.hSR_EEsel_dPhi.Fill(ev.EE_dPhi[maxIxyi], weight)
                    self.hSR_EEsel_leadingPt.Fill(ev.EE_leadingPt[maxIxyi], weight)
                    self.hSR_EEsel_subleadingPt.Fill(ev.EE_subleadingPt[maxIxyi], weight)
                    self.hSR_EEsel_Chi2.Fill(ev.EE_normalizedChi2[maxIxyi], weight)

                if abs(ev.EE_dPhi[maxIxyi]) > math.pi/2.0:

                    self.hCR_EEsel_minIxy.Fill(maxIxy, weight)
                    self.hCR_EEsel_invMass.Fill(ev.EE_invMass[maxIxyi], weight)
                    self.hCR_EEsel_dPhi.Fill(ev.EE_dPhi[maxIxyi], weight)
                    self.hCR_EEsel_leadingPt.Fill(ev.EE_leadingPt[maxIxyi], weight)
                    self.hCR_EEsel_subleadingPt.Fill(ev.EE_subleadingPt[maxIxyi], weight)
                    self.hCR_EEsel_Chi2.Fill(ev.EE_normalizedChi2[maxIxyi], weight)

        #### Dimuon channel
        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6 and not ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8 and ev.nMuonCandidate > 1:

            # Single dimuon histograms
            mi = list(range(0, ev.nMuonCandidate))
            mi = sorted(mi, reverse=True, key=lambda x: ev.MuonCandidate_pt[x]) # sort muons by p
            self.hM1_pt.Fill(ev.MuonCandidate_pt[mi[0]])
            self.hM2_pt.Fill(ev.MuonCandidate_pt[mi[1]])
            self.hM1_triggerPt.Fill(ev.MuonCandidate_triggerPt[mi[0]])
            self.hM2_triggerPt.Fill(ev.MuonCandidate_triggerPt[mi[1]])
            self.hM1_eta.Fill(ev.MuonCandidate_eta[mi[0]])
            self.hM2_eta.Fill(ev.MuonCandidate_eta[mi[1]])
            self.hM1_reliso.Fill(ev.IsoTrackSel_relPfIsolationDR03[ev.MuonCandidate_isotrackIdx[mi[0]]])
            self.hM2_reliso.Fill(ev.IsoTrackSel_relPfIsolationDR03[ev.MuonCandidate_isotrackIdx[mi[1]]])

            
            # Dimuon histograms
            self.hnMM.Fill(ev.nMM, weight)

            maxIxy = -1
            maxIxyi = -99

            for mm in range(0, ev.nMM):

                if abs(ev.MuonCandidate_eta[ev.MM_idxA[mm]]) > 2: continue
                if abs(ev.MuonCandidate_eta[ev.MM_idxB[mm]]) > 2: continue
                if ev.MM_relisoA[mm] > 0.1 or ev.MM_relisoB[mm] > 0.1: continue
                if ev.MM_normalizedChi2[mm] > 5: continue
                if ev.MM_leadingPt[mm] < 31: continue
                if ev.MM_subleadingPt[mm] < 31: continue
                if ev.MM_invMass[mm] < 15: continue
                if ev.MM_cosAlpha[mm] < -0.79: continue
                if ev.IsoTrackSel_charge[ev.MuonCandidate_isotrackIdx[ev.MM_idxA[mm]]]*ev.IsoTrackSel_charge[ev.MuonCandidate_isotrackIdx[ev.MM_idxB[mm]]] > 1: continue
                
                va = r.TVector3()
                vb = r.TVector3()
                va.SetPtEtaPhi(ev.MuonCandidate_pt[ev.MM_idxA[mm]], ev.MuonCandidate_eta[ev.MM_idxA[mm]], ev.MuonCandidate_phi[ev.MM_idxA[mm]])
                vb.SetPtEtaPhi(ev.MuonCandidate_pt[ev.MM_idxB[mm]], ev.MuonCandidate_eta[ev.MM_idxB[mm]], ev.MuonCandidate_phi[ev.MM_idxB[mm]])
                if va.DeltaR(vb) < 0.2: continue

                IxyA = ev.MuonCandidate_dxySignificance[ev.MM_idxA[mm]]
                IxyB = ev.MuonCandidate_dxySignificance[ev.MM_idxB[mm]]

                mIxy = IxyA if IxyA < IxyB else IxyB


                if mIxy > maxIxy:
                    maxIxy = mIxy
                    maxIxyi = mm

            if maxIxyi > -1: 

                if abs(ev.MM_dPhi[maxIxyi]) < math.pi/2.0: 

                    self.hSR_MMsel_minIxy.Fill(maxIxy, weight)
                    self.hSR_MMsel_invMass.Fill(ev.MM_invMass[maxIxyi], weight)
                    self.hSR_MMsel_dPhi.Fill(ev.MM_dPhi[maxIxyi], weight)
                    self.hSR_MMsel_leadingPt.Fill(ev.MM_leadingPt[maxIxyi], weight)
                    self.hSR_MMsel_subleadingPt.Fill(ev.MM_subleadingPt[maxIxyi], weight)
                    self.hSR_MMsel_Chi2.Fill(ev.MM_normalizedChi2[maxIxyi], weight)
                    self.hSR_MMsel_cosAlpha.Fill(ev.MM_cosAlpha[maxIxyi], weight)
        
                if abs(ev.MM_dPhi[maxIxyi]) > math.pi/2.0: 

                    self.hCR_MMsel_minIxy.Fill(maxIxy, weight)
                    self.hCR_MMsel_invMass.Fill(ev.MM_invMass[maxIxyi], weight)
                    self.hCR_MMsel_dPhi.Fill(ev.MM_dPhi[maxIxyi], weight)
                    self.hCR_MMsel_leadingPt.Fill(ev.MM_leadingPt[maxIxyi], weight)
                    self.hCR_MMsel_subleadingPt.Fill(ev.MM_subleadingPt[maxIxyi], weight)
                    self.hCR_MMsel_Chi2.Fill(ev.MM_normalizedChi2[maxIxyi], weight)
                    self.hCR_MMsel_cosAlpha.Fill(ev.MM_cosAlpha[maxIxyi], weight)

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
    
        output = r.TFile(self.fileName, 'UPDATE')

        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Write()

        
        output.Close()


