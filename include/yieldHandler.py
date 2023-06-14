import math
import os
import json
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager
from include.processHandler import processHandler
import numpy as np


class yieldHandler(processHandler):

    def __init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, config, year, raw):

        processHandler.__init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, config, year, raw)
        self.initHistograms()



    #### ------------------------
    #### --
    #### ---- declareHistograms
    #### --
    #### ------------------------

    def initHistograms(self):

        #### ---------------
        #### ---- EE plots
        #### ---------------
        self.hEE_yields = {}
        self.hEE_wPU = {}
        self.hEE_wSF = {}
        for reg in self.dielectronRegions:
            region = reg[0]
            self.hEE_yields[region] = r.TH1F('hEE_{0}'.format(region) + self.sufix, ';;Dielectron vertex yield', 1, 0, 1)
            self.hEE_wPU[region] = r.TH1F('hEE_{0}_wPU'.format(region) + self.sufix, ';;Dielectron w_{PU}', 100, 0, 2)
            self.hEE_wSF[region] = r.TH1F('hEE_{0}_wSF'.format(region) + self.sufix, ';;Dielectron w_{SF}', 100, 0, 10)

        #### ---------------
        #### ---- MM plots
        #### ---------------
        self.hMM_yields  = {}
        self.hMM_wPU = {}
        self.hMM_wSF = {}
        for reg in self.dimuonRegions:
            region = reg[0]
            self.hMM_yields[region] = r.TH1F('hMM_{0}'.format(region) + self.sufix, ';;Dimuon vertex yield', 1, 0, 1)
            self.hMM_wPU[region] = r.TH1F('hMM_{0}_wPU'.format(region) + self.sufix, ';;Dimuon w_{PU}', 100, 0, 2)
            self.hMM_wSF[region] = r.TH1F('hMM_{0}_wSF'.format(region) + self.sufix, ';;Dimuon w_{SF}', 100, 0, 10)

        #### --------------------
        #### ---- Summary plots
        #### --------------------
        self.hEE_summary = r.TH2F('hEE_summary'+ self.sufix, ';Dimuon regions;Dielectron regions', len(self.dimuonRegions) + 1, 0, len(self.dimuonRegions) + 1, len(self.dielectronRegions) + 1, 0, len(self.dielectronRegions) + 1)
        self.hMM_summary = r.TH2F('hMM_summary'+ self.sufix, ';Dimuon regions;Dielectron regions', len(self.dimuonRegions) + 1, 0, len(self.dimuonRegions) + 1, len(self.dielectronRegions) + 1, 0, len(self.dielectronRegions) + 1)
        self.hEEMM_summary = r.TH2F('hEEMM_summary'+ self.sufix, ';Dimuon regions;Dielectron regions', len(self.dimuonRegions) + 1, 0, len(self.dimuonRegions) + 1, len(self.dielectronRegions) + 1, 0, len(self.dielectronRegions) + 1)
        for re,region in enumerate(self.dimuonRegions):
            self.hEE_summary.GetXaxis().SetBinLabel(re+2, region[0])
            self.hMM_summary.GetXaxis().SetBinLabel(re+2, region[0])
            self.hEEMM_summary.GetXaxis().SetBinLabel(re+2, region[0])
        for re,region in enumerate(self.dielectronRegions):
            self.hEE_summary.GetYaxis().SetBinLabel(re+2, region[0])
            self.hMM_summary.GetYaxis().SetBinLabel(re+2, region[0])
            self.hEEMM_summary.GetYaxis().SetBinLabel(re+2, region[0])


        #### --------------------
        #### ---- Records
        #### --------------------
        self.record = True
        self.recordpath = '/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/UL-analysis/CMSSW_10_6_20/src/MyAnalysis/FastPR-Galapago/Galapago-Framework' + '/records/'
        if self.record:
            for re,region in enumerate(self.dimuonRegions):
                with open(self.recordpath + "DMDM_" + self.samplename + "_" + region[0] + "_" + self.samplenumber +".txt", "w") as f: f.close()
            for re,region in enumerate(self.dielectronRegions):
                with open(self.recordpath + "EE_" + self.samplename + "_" + region[0] + "_" + self.samplenumber + ".txt", "w") as f: f.close()

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

    #### ----------------------------
    #### --
    #### ---- Core event processing
    #### --   ( Overriden to fill histograms)
    #### ----------------------------

    def processEvent(self, ev):

        passMuonTrigger = eval(self.mumu_path)
        passElectronTrigger = eval(self.ee_path)

        if (not passMuonTrigger) and (not passElectronTrigger):
            return

        # Weight definition:
        weight = self.getWeight(ev, self.isdata)


        # Test PV acceptance:
        if not ev.PV_passAcceptance: return


        # Select MM and EE candidates
        mm_maxIxy = -99
        try:
            mm_maxIxy, nBSMM = self.processDimuons(ev)
        except AttributeError:
            print('There are some collections missed in this file: Dimuon histograms will be empty')

        ee_maxIxy = -99
        try:
            ee_maxIxy, nBSEE = self.processDielectrons(ev)
        except AttributeError:
            print('There are some collections missed in this file: Dielectron histograms will be empty')


        # Fill Dimuon channel histograms

        if not mm_maxIxy < 0 and passMuonTrigger:

            imm = mm_maxIxy

            ## Compute SF
            sf = 1.0
            if not self.isdata:
                sf = self.getDimuonSF(ev, mm_maxIxy)

            for region in self.dimuonRegions:
                if eval(region[1]):
                    self.hMM_yields[region[0]].Fill(0, weight*sf)
                    self.hMM_wPU[region[0]].Fill(ev.wPU)
                    self.hMM_wSF[region[0]].Fill(sf)
                    if self.record:
                        with open(self.recordpath + "DMDM_" + self.samplename + "_" + region[0] + "_" + self.samplenumber + ".txt", "a") as f: 
                            f.write("New: " + str(ev.DMDM_Lxy_PV[imm]) + '\t' + str(ev.DMDM_trackIxy_PV[imm]) + '\t' + str(ev.DMDM_leadingPt[imm]) + '\t' + str(ev.DMDM_subleadingPt[imm]) + '\t' + str(ev.DMDM_dR[imm]) + '\t' + str(ev.DMDM_mass[imm]) + '\n')
                            f.close()


        # Fill Dielectron channel histograms
        #passElectronTrigger = eval(self.ee_path) and not eval(self.mumu_path)
        ### Note:
        # Orthogonal: if not ee_maxIxy < 0 and passElectronTrigger and not (not mm_maxIxy < 0 and passMuonTrigger):
        # With overlap: if not ee_maxIxy < 0 and passElectronTrigger:
        #if not ee_maxIxy < 0 and passElectronTrigger and not (not mm_maxIxy < 0 and passMuonTrigger):
        if not ee_maxIxy < 0 and passElectronTrigger:

            iee = ee_maxIxy

            ## Compute SF
            sf = 1.0

            if not self.isdata:
                sf = self.getDielectronSF(ev, ee_maxIxy)

            for region in self.dielectronRegions:
                if eval(region[1]):
                    self.hEE_yields[region[0]].Fill(0, weight*sf)
                    self.hEE_wPU[region[0]].Fill(ev.wPU)
                    self.hEE_wSF[region[0]].Fill(sf)
                    if self.record:
                        with open(self.recordpath + "EE_" + self.samplename + "_" + region[0] + "_" + self.samplenumber + ".txt", "a") as f: 
                            f.write(str(ev.Event_event) + '\t' + str(ev.EE_Lxy_PV[iee]) + '\t' + str(ev.EE_trackIxy_PV[iee]) + '\t' + str(ev.EE_leadingEt[iee]) + '\t' + str(ev.EE_subleadingEt[iee]) + '\t' + str(ev.EE_dR[iee]) + '\t' + str(ev.EE_mass[iee]) +'\n')
                            for i in range(ev.nEE):
                                f.write(str(ev.EE_vx[i]) + '\t' + str(ev.EE_vy[i]) + '\t' + str(ev.ElectronCandidate_phi[ev.EE_idxA[i]]) + '\t' + str(ev.ElectronCandidate_phi[ev.EE_idxB[i]])+'\n')
                            f.close()

        ### Fill summary plot
        # Note: Bin label and bin content are shifted by one i.e. bin 1 will contain value 0
        ee_bin = 0
        mm_bin = 0
        for er,ee_region in enumerate(self.dielectronRegions):
            if ee_maxIxy < 0:
                ee_bin = 0
            else:
                iee = ee_maxIxy
                if eval(ee_region[1]):
                    ee_bin = er + 1
            for mr,mm_region in enumerate(self.dimuonRegions):
                if mm_maxIxy < 0:
                    mm_bin = 0
                else:
                    imm = mm_maxIxy
                    if eval(mm_region[1]):
                        mm_bin = mr + 1

        if eval(self.ee_path):
            self.hEE_summary.Fill(mm_bin, ee_bin)
        if eval(self.mumu_path):
            self.hMM_summary.Fill(mm_bin, ee_bin)
        if eval(self.ee_path) and eval(self.mumu_path):
            self.hEEMM_summary.Fill(mm_bin, ee_bin)





