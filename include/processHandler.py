
import ROOT as r




class processHandler:

    def __init__(self, fileName, treename, blockname, samplename):
        self.fileName = fileName
        self.treename = treename
        self.blockname = blockname
        self.samplename = samplename

        self.hEE = r.TH1F('hEE_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)
        self.hMM = r.TH1F('hMM_{0}_{1}_{2}'.format(treename, blockname, samplename), '', 5, 0, 5)


    def processEvent(self, ev):

        self.hEE.Fill(ev.nEE)
        self.hMM.Fill(ev.nMM)


    def Write(self):
    
        output = r.TFile(self.fileName, 'UPDATE')
        self.hEE.Write()
        self.hMM.Write()
        output.Close()


