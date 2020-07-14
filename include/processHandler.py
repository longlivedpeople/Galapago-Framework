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
        self.filename = self.outdir + '{0}_{1}_{2}_{3}.root'.format(treename, blockname, samplename, str(samplenumber))
        self.treename = treename
        self.blockname = blockname
        self.samplename = samplename
        self.samplenumber = samplenumber
        self.lumiweight = lumiweight
        self.isdata = isdata
        self.cutManager = CutManager.CutManager()

        self.hcounts = r.TH1F('hcounts_{0}_{1}_{2}'.format(treename, blockname, samplename, str(samplenumber)), '', 1, 0, 1)

        ###########################
        ###  DiMuon Histograms  ###
        ###########################
        self.hMM_dPhi = r.TH1F('hMM_dPhi_{0}_{1}_{2}_{3}'.format(treename, blockname, samplename, str(samplenumber)), '', 20, -3.3, 3.3)


        for attr, value in self.__dict__.iteritems():
            if attr[0] == 'h': value.Sumw2()



    def processEvent(self, ev, lumiweight, isData):

        if not isData:
            weight = lumiweight*ev.wPU
        else: 
            weight = 1

        #### Generation variables

        self.hcounts.Fill(0, weight)

        """ CUT EXAMPLE
        if eval(self.cutManager.twoElectrons):
            self.hEE.Fill(ev.nEE)
        if eval(self.cutManager.twoMuons):
            self.hMM.Fill(ev.nMM)
        if eval(self.cutManager.Add(self.cutManager.twoElectrons, self.cutManager.twoMuons)):
            self.hLL.Fill(ev.nEE)
            self.hLL.Fill(ev.nMM)
        """
    def processDimuons(self, ev):

        if not self.isdata:
            weight = self.lumiweight*ev.wPU
        else: 
            weight = 1

        if ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 and ev.nDMDMBase > 0:

            self.hMM_dPhi.Fill(ev.DMDMBase_dPhi[ev.DMDMBase_maxIxy], weight)


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


