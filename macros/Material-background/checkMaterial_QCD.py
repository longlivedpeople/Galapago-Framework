
import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]:
    WORKPATH += level
    WORKPATH += '/'

GALAPAGOPATH = ''
for d in WORKPATH.split('/'):
    GALAPAGOPATH += d
    GALAPAGOPATH += '/'
    if d == 'Galapago-Framework': break

sys.path.insert(0, GALAPAGOPATH)

import include.Canvas as Canvas
import include.Sample as Sample
import include.helper as helper
import include.CutManager as CutManager

#print(WORKPATH, WORKPATH)
#print(GALAPAGOPATH, GALAPAGOPATH)


if __name__ == "__main__":


    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ###########################
    ####   Parser object   ####
    ###########################
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    #parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    parser.add_option('-e', '--era', action='store', type=str, dest='era', default='', help='2018 era')
    #parser.add_option('-n', '--nmax', action='store', type=int, dest='nmax', default=0, help='Path to file')
    (opts, args) = parser.parse_args()


    ##################################
    ####   Variable declaration   ####
    ##################################

    dr_logbin = np.logspace(-3, 1, 101)
    iso_logbin = np.logspace(-2, 2, 101)

    #### -----------------
    #### ---- Histograms
    #### -----------------

    hist_2d_BPregion = r.TH2F("hist_2d_BPregion", "", 150, -5, 5, 150, -5, 5)
    hist_2d_20region = r.TH2F("hist_2d_20region", "", 200, -20, 20, 200, -20, 20)
    hist_1d_r  = r.TH1F("hist_1d_r_all", "", 100, 0, 10)

    #########################
    ####   Load sample   ####
    #########################

    treeA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['DoubleMuon_Run2018' + opts.era], 'DATA'), name = '2018' + opts.era, isdata = 0 )


    ###################################
    ####   Loop over tree events   ####
    ###################################
    cm = CutManager.CutManager()

    for b in treeA.blocks:
        for s in b.samples:
            #for f in s.ftpaths:
            #    print(f)
            for t in s.ttrees:
                print('New tree with:', t.GetEntries())
                cutoff = t.GetEntries()/10
                for e,ev in enumerate(t):

                    for j in range(0, ev.nDMDM):
                        imm = j # index to handle DMDM pair
                        R = math.sqrt(ev.DMDM_vx[imm]*ev.DMDM_vx[imm] + ev.DMDM_vy[imm]*ev.DMDM_vy[imm])

                        if not (ev.HLT_DoubleL2Mu23NoVtx_2Cha or ev.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed): continue
                        if not ev.DGM_pt[ev.DMDM_idxA[imm]] > 25 or not ev.DGM_pt[ev.DMDM_idxB[imm]] > 25: continue
                        if not abs(ev.DGM_eta[ev.DMDM_idxA[imm]]) < 1.0 or not abs(ev.DGM_eta[ev.DMDM_idxB[imm]]) < 1.0: continue

                        if not ev.DMDM_normalizedChi2[imm] < 10: continue
                        if ev.DMDM_dR[imm] > 0.2 and ev.DMDM_mass[imm] > 15 and eval(cm.MM_iso0l) and ev.DMDM_cosAlpha[imm] > -0.8:
                            hist_1d_r.Fill(R)
                            if R > 1.5:
                                hist_2d_BPregion.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])
                                hist_2d_20region.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])

                        """
                        if not ev.DGM_ptError[ev.DMDM_idxB[imm]]/ev.DGM_pt[ev.DMDM_idxB[imm]] < 0.3: continue
:u
                        if not ev.DGM_ptError[ev.DMDM_idxA[imm]]/ev.DGM_pt[ev.DMDM_idxA[imm]] < 0.3: continue
                        if not ev.DGM_normChi2[ev.DMDM_idxA[imm]] < 7.5: continue
                        if not ev.DGM_normChi2[ev.DMDM_idxB[imm]] < 7.5: continue
                        if not ev.DGM_muonHits[ev.DMDM_idxA[imm]] > 11: continue  
                        if not ev.DGM_muonHits[ev.DMDM_idxB[imm]] > 11: continue  
                        if not ev.DGM_outerTrackerHits[ev.DMDM_idxA[imm]] > 8: continue  
                        if not ev.DGM_outerTrackerHits[ev.DMDM_idxB[imm]] > 8: continue
                        """


    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f_MM_QCD'+opts.era+'.root', 'RECREATE')


    #### Write everything to use later:
    hist_1d_r.Write()
    hist_2d_BPregion.Write()
    hist_2d_20region.Write()
    outputFile.Close()




