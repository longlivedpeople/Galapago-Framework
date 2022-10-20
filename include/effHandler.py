import math
import os
import json
import ROOT as r
from ROOT import TVector3, TLorentzVector
import include.CutManager as CutManager
from include.processHandler import processHandler
import numpy as np


class effHandler(processHandler):

    def __init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, config, year):

        processHandler.__init__(self, outdir, treename, blockname, samplename, samplenumber, lumiweight, isdata, config, year)


        ### Get the vertex selection from configuration file to define
        ### the sequential cuts that are evaluated

        self.mumu_sequence = []
        self.ee_sequence   = []
        if year=='2016':
            self.mumu_sequence = self.config['selection']['muons']['2016']
            self.ee_sequence   = self.config['selection']['electrons']['2016']
        elif year=='2017':
            self.ee_sequence   = self.config['selection']['electrons']['2017']
        elif year=='2018':
            self.mumu_sequence = self.config['selection']['muons']['2018']
            self.ee_sequence   = self.config['selection']['electrons']['2018']


        ### Create histograms for sequential cut efficiency

        nCutsEE = len(self.ee_sequence)
        self.hEE_efficiency = r.TH1F('hEE_efficiency' + self.sufix, ';;Efficiency', nCutsEE, 0, nCutsEE)
        self.hEE_efficiency.Sumw2()

        nCutsMM = len(self.mumu_sequence)
        self.hMM_efficiency = r.TH1F('hMM_efficiency' + self.sufix, ';;Efficiency', nCutsMM, 0, nCutsMM)
        self.hMM_efficiency.Sumw2()


    #### ----------------------------
    #### --
    #### ---- Core event processing
    #### --   ( Overriden to fill histograms)
    #### ----------------------------

    def processEvent(self, ev):

        # Weight definition:
        weight = self.getWeight(ev, self.isdata)


        # Test PV acceptance:
        if not ev.PV_passAcceptance: return


        # Dimuon processing
        try:

            if eval(self.mumu_path) and ev.nDMDM > 0:
                haveMM = [False]*len(self.mumu_sequence)
                for i in range(0, ev.nDMDM):
                    imm = i
                    for j in range(0, len(self.mumu_sequence)):
                        if not eval(self.mumu_sequence[j]): break 
                        haveMM[j] = True
                
                for j in range(0, len(self.mumu_sequence)):
                    if haveMM[j]:
                        self.hMM_efficiency.Fill(j)

        except AttributeError:
            print('There are some collections missed in this file: Dimuon histograms will be empty')


        # Dielectron processing
        try:

            if eval(self.ee_path) and ev.nEE > 0:
                haveEE = [False]*len(self.ee_sequence)
                for i in range(0, ev.nEE):
                    iee = i
                    for j in range(0, len(self.ee_sequence)):
                        if not eval(self.ee_sequence[j]): break 
                        haveEE[j] = True
                
                for j in range(0, len(self.ee_sequence)):
                    if haveEE[j]:
                        self.hEE_efficiency.Fill(j)

        except AttributeError:
            print('There are some collections missed in this file: Dielectron histograms will be empty')


