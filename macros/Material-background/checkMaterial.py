
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
    hist_1d_r_all  = r.TH1F("hist_1d_r_all", "", 100, 0, 100)
    hist_1d_r_trig = r.TH1F("hist_1d_r_trig", "", 100, 0, 100)
    hist_1d_r_kin  = r.TH1F("hist_1d_r_kin", "", 100, 0, 100)
    hist_1d_r_chi  = r.TH1F("hist_1d_r_chi", "", 100, 0, 100)
    hist_1d_r_char = r.TH1F("hist_1d_r_char", "", 100, 0, 100)
    hist_1d_r_dr   = r.TH1F("hist_1d_r_dr", "", 100, 0, 100)
    hist_1d_r_mass = r.TH1F("hist_1d_r_mass", "", 100, 0, 100)
    hist_1d_r_iso  = r.TH1F("hist_1d_r_iso", "", 100, 0, 100)
    hist_1d_r_cos  = r.TH1F("hist_1d_r_cos", "", 100, 0, 100)
    hist_1d_r_id   = r.TH1F("hist_1d_r_id", "", 100, 0, 100)

    hist_2d_BPregion = r.TH2F("hist_2d_BPregion", "", 150, -5, 5, 150, -5, 5)
    hist_2d_20region = r.TH2F("hist_2d_20region", "", 200, -20, 20, 200, -20, 20)
    hist_2d_BPregion_id = r.TH2F("hist_2d_BPregion_id", "", 150, -5, 5, 150, -5, 5)
    hist_2d_20region_id = r.TH2F("hist_2d_20region_id", "", 200, -20, 20, 200, -20, 20)

    hist_2d_20region_SS = r.TH2F("hist_2d_20region_SS", "", 200, -20, 20, 200, -20, 20)
    hist_1d_r_SS        = r.TH1F("hist_1d_r_SS", "", 100, 0, 100)
    hist_1d_dr_SS       = r.TH1F("hist_1d_dr_SS", "", len(dr_logbin) - 1, dr_logbin) 
    hist_1d_iso_SS      = r.TH1F("hist_1d_iso_SS", "", len(iso_logbin) - 1, iso_logbin) 
    hist_1d_dPhi_SS     = r.TH1F("hist_1d_dPhi_SS", "", 35, 0, 3.15)
    hist_1d_mass_SS     = r.TH1F("hist_1d_mass_SS", "", 100, 0, 100)
    hist_1d_dpt_SS      = r.TH1F("hist_1d_dpt_SS", "", 100, 0, 1)    

    hist_2d_20region_SS_ring = r.TH2F("hist_2d_20region_SS_ring", "", 200, -20, 20, 200, -20, 20)
    hist_1d_r_SS_ring        = r.TH1F("hist_1d_r_SS_ring", "", 100, 0, 100)
    hist_1d_dr_SS_ring       = r.TH1F("hist_1d_dr_SS_ring", "", len(dr_logbin) - 1, dr_logbin) 
    hist_1d_iso_SS_ring      = r.TH1F("hist_1d_iso_SS_ring", "", len(iso_logbin) - 1, iso_logbin) 
    hist_1d_dPhi_SS_ring     = r.TH1F("hist_1d_dPhi_SS_ring", "", 35, 0, 3.15)
    hist_1d_mass_SS_ring     = r.TH1F("hist_1d_mass_SS_ring", "", 100, 0, 100)
    hist_1d_dpt_SS_ring      = r.TH1F("hist_1d_dpt_SS_ring", "", 100, 0, 1)    

    hist_2d_BPregion_QCD  = r.TH2F("hist_2d_BPregion_QCD", "", 150, -5, 5, 150, -5, 5)
    hist_2d_20region_QCD  = r.TH2F("hist_2d_20region_QCD", "", 200, -20, 20, 200, -20, 20)
    hist_1d_r_QCD         = r.TH1F("hist_1d_r_QCD", "", 100, 0, 100)

    #########################
    ####   Load sample   ####
    #########################
    """
    _dirName = opts.filename
    _tree = r.TChain('Events')
    for _file in os.listdir(_dirName):
        if '.root' not in _file: continue
        _tree.Add(_dirName + _file)
   
    _file = r.TFile(opts.filename)
    _tree = _file.Get('Events')
    print("TTree with " + str(_tree.GetEntries()) + " entries")
    """

    """ 
    treeA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['DoubleMuon_Run2018A'], 'DATA'), name = '2018A', isdata = 0 )
    treeB = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['DoubleMuon_Run2018B'], 'DATA'), name = '2018B', isdata = 0 )
    treeC = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['DoubleMuon_Run2018C'], 'DATA'), name = '2018C', isdata = 0 )
    treeD = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['DoubleMuon_Run2018D'], 'DATA'), name = '2018D', isdata = 0 )
    """
    treeA = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/Samples_cern_UltraLegacy.dat', ['DoubleMuon_Run2018' + opts.era], 'DATA'), name = '2018' + opts.era, isdata = 0 )


    ###################################
    ####   Loop over tree events   ####
    ###################################
    cm = CutManager.CutManager()

    for b in treeA.blocks:
        for s in b.samples:
            for t in s.ttrees:
                print('New tree with:', t.GetEntries())
                cutoff = t.GetEntries()/10
                for e,ev in enumerate(t):
                    if e > cutoff: break

                    for j in range(0, ev.nDMDM):
                        imm = j # index to handle DMDM pair
                        R = math.sqrt(ev.DMDM_vx[imm]*ev.DMDM_vx[imm] + ev.DMDM_vy[imm]*ev.DMDM_vy[imm])

                        hist_1d_r_all.Fill(R)

                        if not (ev.HLT_DoubleL2Mu23NoVtx_2Cha or ev.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed): continue
                        hist_1d_r_trig.Fill(R)

                        if not ev.DGM_pt[ev.DMDM_idxA[imm]] > 25 or not ev.DGM_pt[ev.DMDM_idxB[imm]] > 25: continue
                        if not abs(ev.DGM_eta[ev.DMDM_idxA[imm]]) < 1.0 or not abs(ev.DGM_eta[ev.DMDM_idxB[imm]]) < 1.0: continue
                        hist_1d_r_kin.Fill(R)

                        if not ev.DMDM_normalizedChi2[imm] < 10: continue
                        hist_1d_r_chi.Fill(R)
                        if R > 1.8:
                            hist_2d_BPregion.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])
                            hist_2d_20region.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])

                        # first check spurious SS dGlobal pairs:
                        if ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] > 0:
                            hist_2d_20region_SS.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])
                            hist_1d_r_SS.Fill(R)
                            hist_1d_dr_SS.Fill(ev.DMDM_dR[imm])
                            hist_1d_mass_SS.Fill(ev.DMDM_mass[imm])
                            hist_1d_iso_SS.Fill(ev.DGM_relPFiso[ev.DMDM_idxA[imm]])
                            hist_1d_iso_SS.Fill(ev.DGM_relPFiso[ev.DMDM_idxB[imm]])
                            hist_1d_dPhi_SS.Fill(abs(ev.DMDM_dPhi[imm]))
                            hist_1d_dpt_SS.Fill((ev.DMDM_leadingPt[imm] - ev.DMDM_subleadingPt[imm])/ev.DMDM_leadingPt[imm])
                            if R > 6 and R < 14:
                                hist_2d_20region_SS_ring.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])
                                hist_1d_r_SS_ring.Fill(R)
                                hist_1d_dr_SS_ring.Fill(ev.DMDM_dR[imm])
                                hist_1d_mass_SS_ring.Fill(ev.DMDM_mass[imm])
                                hist_1d_iso_SS_ring.Fill(ev.DGM_relPFiso[ev.DMDM_idxA[imm]])
                                hist_1d_iso_SS_ring.Fill(ev.DGM_relPFiso[ev.DMDM_idxB[imm]])
                                hist_1d_dPhi_SS_ring.Fill(abs(ev.DMDM_dPhi[imm]))
                                hist_1d_dpt_SS_ring.Fill((ev.DMDM_leadingPt[imm] - ev.DMDM_subleadingPt[imm])/ev.DMDM_leadingPt[imm])
                        
                                print(ev.Event_event, ev.Event_run, ev.Event_luminosityBlock)
                                print("> ", ev.DGM_pt[ev.DMDM_idxA[imm]], ev.DGM_pt[ev.DMDM_idxB[imm]])


                        if ev.DMDM_dR[imm] > 0.2 and ev.DMDM_mass[imm] > 15 and eval(cm.MM_iso0l) and ev.DMDM_cosAlpha[imm] > -0.8 and ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] < 0:
                            hist_1d_r_QCD.Fill(R)
                            if R > 1.5:
                                hist_2d_BPregion_QCD.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])
                                hist_2d_20region_QCD.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])



                        if not ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] < 0: continue
                        hist_1d_r_char.Fill(R)

                        if not ev.DMDM_dR[imm] > 0.2: continue
                        hist_1d_r_dr.Fill(R)

                        if not ev.DMDM_mass[imm] > 15: continue
                        hist_1d_r_mass.Fill(R)

                        if not eval(cm.MM_iso2l): continue
                        hist_1d_r_iso.Fill(R)

                        if not ev.DMDM_cosAlpha[imm] > -0.8: continue
                        hist_1d_r_cos.Fill(R)

                        if not ev.DGM_ptError[ev.DMDM_idxB[imm]]/ev.DGM_pt[ev.DMDM_idxB[imm]] < 0.3: continue
                        if not ev.DGM_ptError[ev.DMDM_idxA[imm]]/ev.DGM_pt[ev.DMDM_idxA[imm]] < 0.3: continue
                        if not ev.DGM_normChi2[ev.DMDM_idxA[imm]] < 7.5: continue
                        if not ev.DGM_normChi2[ev.DMDM_idxB[imm]] < 7.5: continue
                        if not ev.DGM_muonHits[ev.DMDM_idxA[imm]] > 11: continue  
                        if not ev.DGM_muonHits[ev.DMDM_idxB[imm]] > 11: continue  
                        if not ev.DGM_outerTrackerHits[ev.DMDM_idxA[imm]] > 8: continue  
                        if not ev.DGM_outerTrackerHits[ev.DMDM_idxB[imm]] > 8: continue
                        hist_1d_r_id.Fill(R)
                        if True:
                            hist_2d_BPregion_id.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])
                            hist_2d_20region_id.Fill(ev.DMDM_vx[imm], ev.DMDM_vy[imm])



    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f'+opts.era+'.root', 'RECREATE')


    #### Write everything to use later:
    hist_1d_r_all.Write()
    hist_1d_r_trig.Write()
    hist_1d_r_kin.Write()
    hist_1d_r_chi.Write()
    hist_1d_r_char.Write()
    hist_1d_r_dr.Write()
    hist_1d_r_mass.Write()
    hist_1d_r_iso.Write()
    hist_1d_r_cos.Write()
    hist_1d_r_id.Write()
    hist_2d_BPregion_id.Write()
    hist_2d_20region_id.Write()
    hist_2d_BPregion.Write()
    hist_2d_20region.Write()
    hist_2d_20region_SS.Write()
    hist_1d_r_SS.Write()
    hist_1d_dr_SS.Write()
    hist_1d_mass_SS.Write()
    hist_1d_iso_SS.Write()
    hist_1d_dPhi_SS.Write()
    hist_1d_dpt_SS.Write()
    hist_2d_20region_SS_ring.Write()
    hist_1d_r_SS_ring.Write()
    hist_1d_dr_SS_ring.Write()
    hist_1d_mass_SS_ring.Write()
    hist_1d_iso_SS_ring.Write()
    hist_1d_dPhi_SS_ring.Write()
    hist_1d_dpt_SS_ring.Write()
    hist_2d_BPregion_QCD.Write()
    hist_2d_20region_QCD.Write()
    hist_1d_r_QCD.Write()
    outputFile.Close()




