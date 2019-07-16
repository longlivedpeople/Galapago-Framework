import math
import ROOT as r
import include.CutManager as CutManager



class processHandler:

    def __init__(self, fileName, treename, blockname, samplename):
        self.fileName = fileName
        self.treename = treename
        self.blockname = blockname
        self.samplename = samplename
        self.cutManager = CutManager.CutManager()
 
        self.hnLL = r.TH1F('hnLL_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)
        self.hnEE = r.TH1F('hnEE_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)
        self.hnMM = r.TH1F('hnMM_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)
        self.hEEsel_minIxy = r.TH1F('hEEsel_minIxy_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 20)
        self.hEEsel_invMass = r.TH1F('hEEsel_invMass_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 50, 0, 500)
        self.hEEsel_Chi2 = r.TH1F('hEEsel_Chi2_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 20, 0, 20)
        self.hEEsel_dPhi = r.TH1F('hEEsel_dPhi_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, math.pi)
        self.hEEsel_leadingPt = r.TH1F('hEEsel_leadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 1000)
        self.hEEsel_subleadingPt = r.TH1F('hEEsel_subleadingPt_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 40, 0, 700)
        

        self.hMMsel_cosAlpha = r.TH1F('hMMsel_cosAlpha_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 22, -1.1, 1.1)

        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Sumw2()



    def processEvent(self, ev, weight):

        #### Dielectron channel
        if ev.Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v8:

            self.hnEE.Fill(ev.nEE, weight)

            maxIxy = -1
            maxIxyi = -99

            for ee in range(0, ev.nEE):

                if ev.EE_relisoA[ee] > 0.1 or ev.EE_relisoB[ee] > 0.1: continue
                if ev.EE_normalizedChi2[ee] > 10: continue
                if ev.EE_leadingPt[ee] < 36: continue
                if ev.EE_subleadingPt[ee] < 21: continue
                if ev.EE_leadingEt[ee] < 40: continue
                if ev.EE_subleadingEt[ee] < 25: continue
                if ev.EE_invMass[ee] < 15: continue
                if abs(ev.EE_dPhi[ee]) > math.pi/2.0: continue

                IxyA = ev.IsoTrackSel_dxySignificance[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[ee]]]
                IxyB = ev.IsoTrackSel_dxySignificance[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[ee]]]

                mIxy = IxyA if IxyA < IxyB else IxyB


                if mIxy > maxIxy:
                    maxIxy = mIxy
                    maxIxyi = ee

            if maxIxyi > -1: 

                self.hEEsel_minIxy.Fill(maxIxy, weight)
                self.hEEsel_invMass.Fill(ev.EE_invMass[maxIxyi], weight)
                self.hEEsel_dPhi.Fill(ev.EE_dPhi[maxIxyi], weight)
                self.hEEsel_leadingPt.Fill(ev.EE_leadingPt[maxIxyi], weight)
                self.hEEsel_subleadingPt.Fill(ev.EE_subleadingPt[maxIxyi], weight)
                self.hEEsel_Chi2.Fill(ev.EE_normalizedChi2[maxIxyi], weight)

        
        #### Dimuon channel
        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v6:

            self.hnMM.Fill(ev.nMM, weight)




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

        self.hnLL.Write()
        self.hnEE.Write()
        self.hnMM.Write()
        self.hEEsel_minIxy.Write()
        self.hEEsel_invMass.Write()
        self.hEEsel_Chi2.Write()
        self.hEEsel_dPhi.Write()
        self.hEEsel_leadingPt.Write()
        self.hEEsel_subleadingPt.Write()

        output.Close()


