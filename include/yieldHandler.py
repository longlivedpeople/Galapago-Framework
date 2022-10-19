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
        for reg in self.dielectronRegions:
            region = reg[0]
            self.hEE_yields[region] = r.TH1F('hEE_{0}'.format(region) + self.sufix, ';;Dielectron vertex yield', 1, 0, 1)

        #### ---------------
        #### ---- MM plots
        #### ---------------
        self.hMM_yields  = {}
        for reg in self.dimuonRegions:
            region = reg[0]
            self.hMM_yields[region] = r.TH1F('hMM_{0}'.format(region) + self.sufix, ';;Dimuon vertex yield', 1, 0, 1)

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

        # Weight definition:
        weight = self.getWeight(ev, self.isdata)


        # Test PV acceptance:
        if not ev.PV_passAcceptance: return


        # Dimuon processing
        try:
            mm_maxIxy = -99
            mm_maxIxy, nBSMM = self.processDimuons(ev)
            if not mm_maxIxy < 0:
                imm = mm_maxIxy

                ## Compute SF
                sf = 1.0
                if not self.isdata:
                    sf = self.getDimuonSF(ev, mm_maxIxy)

                for region in self.dimuonRegions:
                    if eval(region[1]):
                        self.hMM_yields[region[0]].Fill(0, weight*sf)
        except AttributeError:
            print('There are some collections missed in this file: Dimuon histograms will be empty')


        # Dielectron processing
        try:
            ee_maxIxy = -99
            ee_maxIxy, nBSEE = self.processDielectrons(ev)
            if not ee_maxIxy < 0:
                iee = ee_maxIxy
                ## Compute SF
                sf = 1.0
                if not self.isdata:
                    sf = self.getDielectronSF(ev, ee_maxIxy)

                for region in self.dielectronRegions:
                    if eval(region[1]):
                        self.hEE_yields[region[0]].Fill(0, weight*sf)
        except AttributeError:
            print('There are some collections missed in this file: Dielectron histograms will be empty')


