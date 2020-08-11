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

        self.hcounts = r.TH1F('hcounts__{0}__{1}__{2}'.format(treename, blockname, samplename, str(samplenumber)), '', 1, 0, 1)

        ###########################
        ###  DiMuon Histograms  ###
        ###########################
        self.hMM_dPhi = r.TH1F('hMM_dPhi' + self.sufix, ';Dimuon collinearity |#Delta#Phi|;', 20, -3.3, 3.3)
        self.hMM_mass = r.TH1F('hMM_mass' + self.sufix, ';Dimuon invariant mass m_{#mu#mu} (GeV);', 35, 0, 200)
        self.hMM_cosAlpha = r.TH1F('hMM_cosAlpha' + self.sufix, ';Dimuon cos(#alpha_{#mu#mu});', 21, -1.1, 1.1) 
        self.hMM_trackIxy = r.TH1F('hMM_trackIxy' + self.sufix, ';Dimuon |d_{0}|/#sigma_{d};', 20, 0, 20)
        self.hMM_trackDxy = r.TH1F('hMM_trackDxy' + self.sufix, ';Dimuon |d_{0}| (cm);', 20, 0, 5)
        self.hMM_Lxy = r.TH1F('hMM_Lxy' + self.sufix, ';Dimuon vertex |L_{xy}| (cm);', 20, 0, 10)
        self.hMM_Ixy = r.TH1F('hMM_Ixy' + self.sufix, ';Dimuon vertex |L_{xy}|/#sigma_{L};', 20, 0, 20)

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



        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Sumw2()



    def processEvent(self, ev):

        if not self.isdata:
            weight = self.lumiweight*ev.wPU
        else: 
            weight = 1

        #### Generation variables

        self.hcounts.Fill(0, weight)


        ###############################
        ####   Dimuon histograms   ####
        ###############################

        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 and ev.nDMDMBase > 0:

            if eval(self.cutManager.LoopMM_OScharge):

                self.hMM_dPhi.Fill(ev.DMDMBase_dPhi[ev.DMDMBase_maxIxy], weight)
                self.hMM_mass.Fill(ev.DMDMBase_mass[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackIxy.Fill(ev.DMDMBase_trackIxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackDxy.Fill(ev.DMDMBase_trackDxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Lxy.Fill(ev.DMDMBase_Lxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Ixy.Fill(ev.DMDMBase_Ixy[ev.DMDMBase_maxIxy], weight)
                self.hMM_cosAlpha.Fill(ev.DMDMBase_cosAlpha[ev.DMDMBase_maxIxy], weight)



            if self.isdata and eval(self.cutManager.LoopMM_SScharge):

                self.hMM_dPhi_SS.Fill(ev.DMDMBase_dPhi[ev.DMDMBase_maxIxy], weight)
                self.hMM_mass_SS.Fill(ev.DMDMBase_mass[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackIxy_SS.Fill(ev.DMDMBase_trackIxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_trackDxy_SS.Fill(ev.DMDMBase_trackDxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Lxy_SS.Fill(ev.DMDMBase_Lxy[ev.DMDMBase_maxIxy], weight)
                self.hMM_Ixy_SS.Fill(ev.DMDMBase_Ixy[ev.DMDMBase_maxIxy], weight)
                self.hMM_cosAlpha_SS.Fill(ev.DMDMBase_cosAlpha[ev.DMDMBase_maxIxy], weight)

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


