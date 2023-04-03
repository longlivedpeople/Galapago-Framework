import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3, TLorentzVector, TMath
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

from include.galapagoStyle import sigpalette, gcolors, dcolors, acolors

def combinePlots(plotList):

    plot = plotList[0]
    for i in range(1, len(plotList)):
        plot.Add(plotList[i])

    return plot


################################################################################################################

if __name__ == "__main__":

    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()
    r.gStyle.SetPalette(r.kBird)


    #########################
    ####   Load sample   ####
    #########################
    year = '2018'
    HSS_Signals = []
    HSS_Signals.append('HSS_1000_150_1_'+year)
    HSS_Signals.append('HSS_1000_150_10_'+year)
    HSS_Signals.append('HSS_1000_150_100_'+year)
    HSS_Signals.append('HSS_1000_150_1000_'+year)
    HSS_Signals.append('HSS_1000_150_10000_'+year)


    #########################
    ####   Init plots    ####
    #########################

    Lxy_binning = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 17., 19., 22., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70.])

    plot = {}
    for signal in HSS_Signals:

        chi_logbin = np.logspace(-3, 2, 51)

        ## Efficiencies
        plot[signal + '_EE_vertexEff_Lxy']  = copy.deepcopy(r.TEfficiency(signal + '_EE_vertexEff_Lxy', ";Generated LLP decay radius (cm);Vertex efficiency", 10, 0, 70))    
        plot[signal + '_MM_vertexEff_Lxy']  = copy.deepcopy(r.TEfficiency(signal + '_MM_vertexEff_Lxy', ";Generated LLP decay radius (cm);Vertex efficiency", 10, 0, 70))    
        plot[signal + '_EE_vertexEff_Lxy_4e']  = copy.deepcopy(r.TEfficiency(signal + '_EE_vertexEff_Lxy_4e', ";Generated LLP decay radius (cm);Vertex efficiency", 10, 0, 70))    
        plot[signal + '_MM_vertexEff_Lxy_4mu']  = copy.deepcopy(r.TEfficiency(signal + '_MM_vertexEff_Lxy_4mu', ";Generated LLP decay radius (cm);Vertex efficiency", 10, 0, 70))    
        plot[signal + '_EE_absEff_Lxy']  = copy.deepcopy(r.TEfficiency(signal + '_EE_absEff_Lxy', ";Generated LLP decay radius (cm);Vertex efficiency", len(Lxy_binning) -1, Lxy_binning))    
        plot[signal + '_MM_absEff_Lxy']  = copy.deepcopy(r.TEfficiency(signal + '_MM_absEff_Lxy', ";Generated LLP decay radius (cm);Vertex efficiency", len(Lxy_binning) -1, Lxy_binning))    
        plot[signal + '_MM_normChi2']  = copy.deepcopy(r.TH1F(signal + '_MM_normChi2', ";Vertex fit #chi^{2}/ndof;Dimuon vertex yield", len(chi_logbin)-1, chi_logbin))    
        plot[signal + '_EE_normChi2']  = copy.deepcopy(r.TH1F(signal + '_EE_normChi2', ";Vertex fit #chi^{2}/ndof;Dielectron vertex yield", len(chi_logbin)-1, chi_logbin))    
        plot[signal + '_EE_goodComb']  = copy.deepcopy(r.TEfficiency(signal + '_EE_goodComb', ";Dielectron vertices;Good association rate", 1, 0, 1))    
        plot[signal + '_MM_goodComb']  = copy.deepcopy(r.TEfficiency(signal + '_MM_goodComb', ";Dimuon vertices;Good association rate", 1, 0, 1))    
        plot[signal + '_EE_mass_type1']  = copy.deepcopy(r.TH1F(signal + '_EE_mass_type1', ";Dielectron vertex invariant mass;Dielectron vertex yield", 50, 0, 500))    
        plot[signal + '_EE_mass_type2']  = copy.deepcopy(r.TH1F(signal + '_EE_mass_type2', ";Dielectron vertex invariant mass;Dielectron vertex yield", 50, 0, 500))    
        plot[signal + '_EE_mass_type3']  = copy.deepcopy(r.TH1F(signal + '_EE_mass_type3', ";Dielectron vertex invariant mass;Dielectron vertex yield", 50, 0, 500))    
        plot[signal + '_EE_mass_type4']  = copy.deepcopy(r.TH1F(signal + '_EE_mass_type4', ";Dielectron vertex invariant mass;Dielectron vertex yield", 50, 0, 500))    
        plot[signal + '_MM_mass_type1']  = copy.deepcopy(r.TH1F(signal + '_MM_mass_type1', ";Dimuon vertex invariant mass;Dimuon vertex yield", 50, 0, 500))    
        plot[signal + '_MM_mass_type2']  = copy.deepcopy(r.TH1F(signal + '_MM_mass_type2', ";Dimuon vertex invariant mass;Dimuon vertex yield", 50, 0, 500))    
        plot[signal + '_MM_mass_type3']  = copy.deepcopy(r.TH1F(signal + '_MM_mass_type3', ";Dimuon vertex invariant mass;Dimuon vertex yield", 50, 0, 500))    
        plot[signal + '_MM_mass_type4']  = copy.deepcopy(r.TH1F(signal + '_MM_mass_type4', ";Dimuon vertex invariant mass;Dimuon vertex yield", 50, 0, 500))    
        plot[signal + '_EE_RECOstatus_type3'] = copy.deepcopy(r.TH2F(signal + '_EE_RECOstatus_type3', ";S_{1} dielectron reconstruction;S_{2} dielectron reconstruction", 2, 0, 2, 2, 0 , 2))
        plot[signal + '_MM_RECOstatus_type3'] = copy.deepcopy(r.TH2F(signal + '_MM_RECOstatus_type3', ";S_{1} dimuon reconstruction;S_{2} dimuon reconstruction", 2, 0, 2, 2, 0 , 2))

        for p in plot.keys():
            r.SetOwnership(plot[p], 0)

        tree = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/CombSignal_2018UL_Fall22.dat', [signal], 'MC'), name = year, isdata = 0 )
        #tree = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_'+year+'UL.dat', [signal], 'MC'), name = year, isdata = 0 )

        print("Processing " + signal)

        for b in tree.blocks: 
            for s in b.samples: 
                for t in s.ttrees:
                    for e,ev in enumerate(t):

                        ## Identify final state leptons of the hard process
                        gen_e = []
                        gen_mu = []
                        for i in range(0, ev.nHardProcessParticle):
                            #if ev.HardProcessParticle_pt[i] < 20:
                            #    continue
                            #if abs(ev.HardProcessParticle_eta[i]) > 2:
                            #    continue
                            if abs(ev.HardProcessParticle_pdgId[i]) == 11:
                                gen_e.append(i)
                            elif abs(ev.HardProcessParticle_pdgId[i]) == 13:
                                gen_mu.append(i)
                            else:
                                continue


                        ## Displaced electron reconstruction 
                        reco_e = []
                        nrecoe = 0
                        for i in range(0, len(gen_e)):
                            gvec = TVector3() 
                            gvec.SetPtEtaPhi(ev.HardProcessParticle_pt[gen_e[i]], ev.HardProcessParticle_eta[gen_e[i]], ev.HardProcessParticle_phi[gen_e[i]])
                            mindR = 9999.
                            idx = -99
                            for j in range(0, ev.nElectronCandidate):
                                evec = TVector3()
                                evec.SetPtEtaPhi(ev.ElectronCandidate_pt[j], ev.ElectronCandidate_eta[j], ev.ElectronCandidate_phi[j])
                                dR = evec.DeltaR(gvec)
                                if dR < mindR and dR < 0.25:
                                    mindR = dR
                                    idx = j
                            reco_e.append(idx)
                            if idx != -99:
                                nrecoe += 1

                        ## Displaced muon reconstruction
                        reco_mu = []
                        nrecomu = 0
                        for i in range(0, len(gen_mu)):
                            gvec = TVector3() 
                            gvec.SetPtEtaPhi(ev.HardProcessParticle_pt[gen_mu[i]], ev.HardProcessParticle_eta[gen_mu[i]], ev.HardProcessParticle_phi[gen_mu[i]])
                            mindR = 9999.
                            idx = -99
                            for j in range(0, ev.nDGM):
                                muvec = TVector3()
                                muvec.SetPtEtaPhi(ev.DGM_pt[j], ev.DGM_eta[j], ev.DGM_phi[j])
                                dR = muvec.DeltaR(gvec)
                                if dR < mindR and dR < 0.25:
                                    mindR = dR
                                    idx = j
                            reco_mu.append(idx)
                            if idx != -99:
                                nrecomu += 1


                        ## Select the channel and make lepton pairs:
                        is4e    = False
                        is2e2mu = False
                        is2mu   = False
                        gen_ee = []
                        gen_mumu = []
                        reco_ee = [] 
                        reco_mumu = []
                        if len(gen_e) > 1:
                            for i in range(0, len(gen_e)):
                                for j in range(i+1, len(gen_e)):
                                    vxi = ev.HardProcessParticle_vx[gen_e[i]]
                                    vyi = ev.HardProcessParticle_vy[gen_e[i]]
                                    vzi = ev.HardProcessParticle_vz[gen_e[i]]
                                    vxj = ev.HardProcessParticle_vx[gen_e[j]]
                                    vyj = ev.HardProcessParticle_vy[gen_e[j]]
                                    vzj = ev.HardProcessParticle_vz[gen_e[j]]
                                    if vxi == vxj and vyi == vyj and vzi == vzj:
                                        gen_ee.append([gen_e[i], gen_e[j]])
                                        reco_ee.append([reco_e[i], reco_e[j]])

                        if len(gen_mu) > 1:
                            for i in range(0, len(gen_mu)):
                                for j in range(i+1, len(gen_mu)):
                                    vxi = ev.HardProcessParticle_vx[gen_mu[i]]
                                    vyi = ev.HardProcessParticle_vy[gen_mu[i]]
                                    vzi = ev.HardProcessParticle_vz[gen_mu[i]]
                                    vxj = ev.HardProcessParticle_vx[gen_mu[j]]
                                    vyj = ev.HardProcessParticle_vy[gen_mu[j]]
                                    vzj = ev.HardProcessParticle_vz[gen_mu[j]]
                                    if vxi == vxj and vyi == vyj and vzi == vzj:
                                        gen_mumu.append([gen_mu[i], gen_mu[j]])
                                        reco_mumu.append([reco_mu[i], reco_mu[j]])
       

                        # If not within acceptance we skip the event
                        #if not (is2e2mu or is4e or is4mu):
                        #    continue;


                        ## Verbose 
                        if False:
                            print("Generated electrons:", gen_ee, ' from ', gen_e)
                            print("Reconstructed electrons:", reco_ee, nrecoe)
                            print("Generated muon:", gen_mumu, ' from ', gen_mu)
                            print("Reconstructed muon:", reco_mumu, nrecomu)




                        ## Electron pair reconstruction efficiency
                        for i in range(0, len(gen_ee)):
                            reco = reco_ee[i]
                            gen  = gen_ee[i]
                            g0 = gen[0]
                            g1 = gen[1]
                            if abs(ev.HardProcessParticle_eta[g0]) > 2.0 or abs(ev.HardProcessParticle_eta[g1]) > 2.0:
                                continue
                            if ev.HardProcessParticle_pt[g0] < 20 or ev.HardProcessParticle_pt[g1] < 20:
                                continue
                            Lxy = math.sqrt((ev.HardProcessParticle_vx[gen[0]])**2 + (ev.HardProcessParticle_vy[gen[0]])**2)
                            vertex = False
                            for v in range(0, ev.nEE):
                                if ev.EE_idxA[v] in reco and ev.EE_idxB[v] in reco:
                                    vertex = True
                                    plot[signal + '_EE_normChi2'].Fill(ev.EE_normalizedChi2[v])
                                    break
                            plot[signal + '_EE_absEff_Lxy'].Fill(vertex, Lxy)
                            if -99 not in reco:
                                plot[signal + '_EE_vertexEff_Lxy'].Fill(vertex, Lxy)
                                if nrecoe == 4:
                                    plot[signal + '_EE_vertexEff_Lxy_4e'].Fill(vertex, Lxy)


                        ## Muon pair reconstruction efficiency
                        for i in range(0, len(gen_mumu)):
                            reco = reco_mumu[i]
                            gen  = gen_mumu[i]
                            g0 = gen[0]
                            g1 = gen[1]
                            if abs(ev.HardProcessParticle_eta[g0]) > 2.0 or abs(ev.HardProcessParticle_eta[g1]) > 2.0:
                                continue
                            if ev.HardProcessParticle_pt[g0] < 20 or ev.HardProcessParticle_pt[g1] < 20:
                                continue
                            Lxy = math.sqrt((ev.HardProcessParticle_vx[gen[0]])**2 + (ev.HardProcessParticle_vy[gen[0]])**2)
                            vertex = False
                            for v in range(0, ev.nDMDM):
                                if ev.DMDM_idxA[v] in reco and ev.DMDM_idxB[v] in reco:
                                    vertex = True
                                    plot[signal + '_MM_normChi2'].Fill(ev.DMDM_normalizedChi2[v])
                                    break

                            plot[signal + '_MM_absEff_Lxy'].Fill(vertex, Lxy)
                            if -99 not in reco:
                                plot[signal + '_MM_vertexEff_Lxy'].Fill(vertex, Lxy)
                                if nrecomu == 4:
                                    plot[signal + '_MM_vertexEff_Lxy_4mu'].Fill(vertex, Lxy)


                        ##### Combinatorics efficiency
                        ## Good electron combinatorics efficiency
                        type1 = False
                        type2 = False
                        type3 = False
                        type4 = False
                         
                        for i in range(0, ev.nEE):
                            ea = TVector3()
                            ea.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxA[i]], ev.ElectronCandidate_eta[ev.EE_idxA[i]], ev.ElectronCandidate_phi[ev.EE_idxA[i]])
                            eb = TVector3()
                            eb.SetPtEtaPhi(ev.ElectronCandidate_pt[ev.EE_idxB[i]], ev.ElectronCandidate_eta[ev.EE_idxB[i]], ev.ElectronCandidate_phi[ev.EE_idxB[i]])
                            mindR_a = 9999.
                            idx_a = -99
                            mindR_b = 9999.
                            idx_b = -99
                            for j in range(0, len(gen_e)):
                                gvec = TVector3() 
                                gvec.SetPtEtaPhi(ev.HardProcessParticle_pt[gen_e[j]], ev.HardProcessParticle_eta[gen_e[j]], ev.HardProcessParticle_phi[gen_e[j]])
                                dRa = gvec.DeltaR(ea)
                                dRb = gvec.DeltaR(eb)
                                if dRa < mindR_a and dRa < 0.25:
                                    mindR_a = dRa
                                    idx_a = gen_e[j]
                                if dRb < mindR_b and dRb < 0.25:
                                    mindR_b = dRb
                                    idx_b= gen_e[j]

                            for pair in gen_ee:
                                if (idx_a in pair) and (idx_b in pair):
                                    type1 = True
                                if ((idx_a in pair) and (idx_b not in pair)) or ((idx_b in pair) and (idx_a not in pair)):
                                    type2 = True 
                            if len(gen_ee) == 2:
                                if ((idx_a in gen_ee[0]) and (idx_b in gen_ee[1])) or ((idx_b in gen_ee[0]) and (idx_a in gen_ee[1])):
                                    type3 = True

                            if type1:
                                plot[signal + '_EE_mass_type1'].Fill(ev.EE_mass[i])
                            elif type3:
                                plot[signal + '_EE_mass_type3'].Fill(ev.EE_mass[i])
                            elif type2:
                                plot[signal + '_EE_mass_type2'].Fill(ev.EE_mass[i])
                            else:
                                plot[signal + '_EE_mass_type4'].Fill(ev.EE_mass[i])
                            
                            if type3:
                                if -99 in reco_ee[0] and -99 in reco_ee[1]:
                                    plot[signal + '_EE_RECOstatus_type3'].Fill(0, 0)
                                if -99 in reco_ee[0] and -99 not in reco_ee[1]:
                                    plot[signal + '_EE_RECOstatus_type3'].Fill(0, 1)
                                if -99 not in reco_ee[0] and -99 in reco_ee[1]:
                                    plot[signal + '_EE_RECOstatus_type3'].Fill(1, 0)
                                if -99 not in reco_ee[0] and -99 not in reco_ee[1]:
                                    plot[signal + '_EE_RECOstatus_type3'].Fill(1, 1)



                        ## Good muon combinatorics efficiency
                        type1 = False
                        type2 = False
                        type3 = False
                        type4 = False
                         
                        for i in range(0, ev.nDMDM):
                            mua = TVector3()
                            mua.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxA[i]], ev.DGM_eta[ev.DMDM_idxA[i]], ev.DGM_phi[ev.DMDM_idxA[i]])
                            mub = TVector3()
                            mub.SetPtEtaPhi(ev.DGM_pt[ev.DMDM_idxB[i]], ev.DGM_eta[ev.DMDM_idxB[i]], ev.DGM_phi[ev.DMDM_idxB[i]])
                            mindR_a = 9999.
                            idx_a = -99
                            mindR_b = 9999.
                            idx_b = -99
                            for j in range(0, len(gen_mu)):
                                gvec = TVector3() 
                                gvec.SetPtEtaPhi(ev.HardProcessParticle_pt[gen_mu[j]], ev.HardProcessParticle_eta[gen_mu[j]], ev.HardProcessParticle_phi[gen_mu[j]])
                                dRa = gvec.DeltaR(mua)
                                dRb = gvec.DeltaR(mub)
                                if dRa < mindR_a and dRa < 0.25:
                                    mindR_a = dRa
                                    idx_a = gen_mu[j]
                                if dRb < mindR_b and dRb < 0.25:
                                    mindR_b = dRb
                                    idx_b= gen_mu[j]

                            for pair in gen_mumu:
                                if (idx_a in pair) and (idx_b in pair):
                                    type1 = True
                                if ((idx_a in pair) and (idx_b not in pair)) or ((idx_b in pair) and (idx_a not in pair)):
                                    type2 = True 
                            if len(gen_mumu) == 2:
                                if ((idx_a in gen_mumu[0]) and (idx_b in gen_mumu[1])) or ((idx_b in gen_mumu[0]) and (idx_a in gen_mumu[1])):
                                    type3 = True

                            if type1:
                                plot[signal + '_MM_mass_type1'].Fill(ev.DMDM_mass[i])
                            elif type3:
                                plot[signal + '_MM_mass_type3'].Fill(ev.DMDM_mass[i])
                            elif type2:
                                plot[signal + '_MM_mass_type2'].Fill(ev.DMDM_mass[i])
                            else:
                                plot[signal + '_MM_mass_type4'].Fill(ev.DMDM_mass[i])
                            
                            if type3:
                                if -99 in reco_mumu[0] and -99 in reco_mumu[1]:
                                    plot[signal + '_MM_RECOstatus_type3'].Fill(0, 0)
                                if -99 in reco_mumu[0] and -99 not in reco_mumu[1]:
                                    plot[signal + '_MM_RECOstatus_type3'].Fill(0, 1)
                                if -99 not in reco_mumu[0] and -99 in reco_mumu[1]:
                                    plot[signal + '_MM_RECOstatus_type3'].Fill(1, 0)
                                if -99 not in reco_mumu[0] and -99 not in reco_mumu[1]:
                                    plot[signal + '_MM_RECOstatus_type3'].Fill(1, 1)




    ##################################################################################################

    ### Efficiency vs Lxy
    # (vertex efficiency on top of reco)
    eeEff = plot['HSS_1000_150_1_2018_EE_vertexEff_Lxy']
    eeEff.Add(plot['HSS_1000_150_10_2018_EE_vertexEff_Lxy'])
    eeEff.Add(plot['HSS_1000_150_100_2018_EE_vertexEff_Lxy'])
    eeEff.Add(plot['HSS_1000_150_1000_2018_EE_vertexEff_Lxy'])
    eeEff.Add(plot['HSS_1000_150_10000_2018_EE_vertexEff_Lxy'])
    mmEff = plot['HSS_1000_150_1_2018_MM_vertexEff_Lxy']
    mmEff.Add(plot['HSS_1000_150_10_2018_MM_vertexEff_Lxy'])
    mmEff.Add(plot['HSS_1000_150_100_2018_MM_vertexEff_Lxy'])
    mmEff.Add(plot['HSS_1000_150_1000_2018_MM_vertexEff_Lxy'])
    mmEff.Add(plot['HSS_1000_150_10000_2018_MM_vertexEff_Lxy'])
    eeEff_4lep = plot['HSS_1000_150_1_2018_EE_vertexEff_Lxy_4e']
    eeEff_4lep.Add(plot['HSS_1000_150_10_2018_EE_vertexEff_Lxy_4e'])
    eeEff_4lep.Add(plot['HSS_1000_150_100_2018_EE_vertexEff_Lxy_4e'])
    eeEff_4lep.Add(plot['HSS_1000_150_1000_2018_EE_vertexEff_Lxy_4e'])
    eeEff_4lep.Add(plot['HSS_1000_150_10000_2018_EE_vertexEff_Lxy_4e'])
    mmEff_4lep = plot['HSS_1000_150_1_2018_MM_vertexEff_Lxy_4mu']
    mmEff_4lep.Add(plot['HSS_1000_150_10_2018_MM_vertexEff_Lxy_4mu'])
    mmEff_4lep.Add(plot['HSS_1000_150_100_2018_MM_vertexEff_Lxy_4mu'])
    mmEff_4lep.Add(plot['HSS_1000_150_1000_2018_MM_vertexEff_Lxy_4mu'])
    mmEff_4lep.Add(plot['HSS_1000_150_10000_2018_MM_vertexEff_Lxy_4mu'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_vertexEff_Lxy", 'png,pdf', 0.4, 0.7, 0.70, 0.9, 1)
    canvas.addHisto(eeEff,'AP', 'Electron vertex (All pairs)', 'p', r.kBlue, True, 0, marker = 20)
    canvas.addHisto(eeEff_4lep,'P,SAME', 'Electron vertex (4 reco electrons)', 'p', r.kBlue, True, 1, marker = 24)
    canvas.addHisto(mmEff,'P,SAME', 'Muon vertex (All pairs)', 'p', r.kRed, True, 2, marker = 21)
    canvas.addHisto(mmEff_4lep,'P,SAME', 'Muon vertex (4 reco muons)', 'p', r.kRed, True, 3, marker = 25)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.2, 0.28, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.2, 0.24, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.2, 0.2, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.6, ymax = 1.2, outputDir = WORKPATH + 'plots/', inProgress = True)


    plot['HSS_1000_150_1_2018_EE_normChi2'].Scale(1./plot['HSS_1000_150_1_2018_EE_normChi2'].Integral())
    plot['HSS_1000_150_10_2018_EE_normChi2'].Scale(1./plot['HSS_1000_150_10_2018_EE_normChi2'].Integral())
    plot['HSS_1000_150_100_2018_EE_normChi2'].Scale(1./plot['HSS_1000_150_100_2018_EE_normChi2'].Integral())
    plot['HSS_1000_150_1000_2018_EE_normChi2'].Scale(1./plot['HSS_1000_150_1000_2018_EE_normChi2'].Integral())
    canvas = Canvas.Canvas("HSS_1000_150_2018_EE_normChi2", 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1)
    canvas.addHisto(plot['HSS_1000_150_1_2018_EE_normChi2'],'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], True, 0)
    canvas.addHisto(plot['HSS_1000_150_10_2018_EE_normChi2'],'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], True, 1)
    canvas.addHisto(plot['HSS_1000_150_100_2018_EE_normChi2'],'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], True, 2)
    canvas.addHisto(plot['HSS_1000_150_1000_2018_EE_normChi2'],'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], True, 3)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', xlog = True, inProgress = False)

    # (vertex efficiency on top of gen)
    eeEff = plot['HSS_1000_150_1_2018_EE_absEff_Lxy']
    eeEff.Add(plot['HSS_1000_150_10_2018_EE_absEff_Lxy'])
    eeEff.Add(plot['HSS_1000_150_100_2018_EE_absEff_Lxy'])
    eeEff.Add(plot['HSS_1000_150_1000_2018_EE_absEff_Lxy'])
    eeEff.Add(plot['HSS_1000_150_10000_2018_EE_absEff_Lxy'])
    mmEff = plot['HSS_1000_150_1_2018_MM_absEff_Lxy']
    mmEff.Add(plot['HSS_1000_150_10_2018_MM_absEff_Lxy'])
    mmEff.Add(plot['HSS_1000_150_100_2018_MM_absEff_Lxy'])
    mmEff.Add(plot['HSS_1000_150_1000_2018_MM_absEff_Lxy'])
    mmEff.Add(plot['HSS_1000_150_10000_2018_MM_absEff_Lxy'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_absEff_Lxy", 'png,pdf', 0.3, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(eeEff,'AP', 'Displaced electron vertex efficiency', 'p', r.kBlue, True, 0, marker = 20)
    canvas.addHisto(mmEff,'P,SAME', 'Displaced muon vertex efficiency', 'p', r.kRed, True, 1, marker = 21)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.2, 0.28, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.2, 0.24, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.2, 0.2, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False)

    ### Combinatorics
    type1 = plot['HSS_1000_150_1_2018_EE_mass_type1']
    type1.Add(plot['HSS_1000_150_10_2018_EE_mass_type1'])
    type1.Add(plot['HSS_1000_150_100_2018_EE_mass_type1'])
    type1.Add(plot['HSS_1000_150_1000_2018_EE_mass_type1'])
    type1.Add(plot['HSS_1000_150_10000_2018_EE_mass_type1'])
    type2 = plot['HSS_1000_150_1_2018_EE_mass_type2']
    type2.Add(plot['HSS_1000_150_10_2018_EE_mass_type2'])
    type2.Add(plot['HSS_1000_150_100_2018_EE_mass_type2'])
    type2.Add(plot['HSS_1000_150_1000_2018_EE_mass_type2'])
    type2.Add(plot['HSS_1000_150_10000_2018_EE_mass_type2'])
    type3 = plot['HSS_1000_150_1_2018_EE_mass_type3']
    type3.Add(plot['HSS_1000_150_10_2018_EE_mass_type3'])
    type3.Add(plot['HSS_1000_150_100_2018_EE_mass_type3'])
    type3.Add(plot['HSS_1000_150_1000_2018_EE_mass_type3'])
    type3.Add(plot['HSS_1000_150_10000_2018_EE_mass_type3'])
    type4 = plot['HSS_1000_150_1_2018_EE_mass_type4']
    type4.Add(plot['HSS_1000_150_10_2018_EE_mass_type4'])
    type4.Add(plot['HSS_1000_150_100_2018_EE_mass_type4'])
    type4.Add(plot['HSS_1000_150_1000_2018_EE_mass_type4'])
    type4.Add(plot['HSS_1000_150_10000_2018_EE_mass_type4'])
    type1.SetMaximum(1e6)
    type1.SetLineWidth(2)
    type2.SetLineWidth(2)
    type3.SetLineWidth(2)
    type4.SetLineWidth(2)
    canvas = Canvas.Canvas("HSS_1000_150_2018_EE_mass_types", 'png,pdf', 0.43, 0.6, 0.79, 0.75, 1)
    canvas.addHisto(type1,'HIST', 'S#rightarrowee [{0}]'.format(type1.GetEntries()), 'l', acolors['1'], True, 0, doOF= False)
    canvas.addHisto(type2,'HIST,SAME', 'e (from S) + random [{0}]'.format(type2.GetEntries()), 'l', acolors['2'], True, 1, doOF = False)
    canvas.addHisto(type3,'HIST,SAME', 'Different S parent [{0}]'.format(type3.GetEntries()), 'l', acolors['3'], True, 2, doOF = False)
    canvas.addHisto(type4,'HIST,SAME', 'Not gen-matched [{0}]'.format(type4.GetEntries()), 'l', acolors['4'], True, 3, doOF = False)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.4, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.4, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.4, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 1, '', '', outputDir = WORKPATH + 'plots/', inProgress = True)


    type1 = plot['HSS_1000_150_1_2018_MM_mass_type1']
    type1.Add(plot['HSS_1000_150_10_2018_MM_mass_type1'])
    type1.Add(plot['HSS_1000_150_100_2018_MM_mass_type1'])
    type1.Add(plot['HSS_1000_150_1000_2018_MM_mass_type1'])
    type1.Add(plot['HSS_1000_150_10000_2018_MM_mass_type1'])
    type2 = plot['HSS_1000_150_1_2018_MM_mass_type2']
    type2.Add(plot['HSS_1000_150_10_2018_MM_mass_type2'])
    type2.Add(plot['HSS_1000_150_100_2018_MM_mass_type2'])
    type2.Add(plot['HSS_1000_150_1000_2018_MM_mass_type2'])
    type2.Add(plot['HSS_1000_150_10000_2018_MM_mass_type2'])
    type3 = plot['HSS_1000_150_1_2018_MM_mass_type3']
    type3.Add(plot['HSS_1000_150_10_2018_MM_mass_type3'])
    type3.Add(plot['HSS_1000_150_100_2018_MM_mass_type3'])
    type3.Add(plot['HSS_1000_150_1000_2018_MM_mass_type3'])
    type3.Add(plot['HSS_1000_150_10000_2018_MM_mass_type3'])
    type4 = plot['HSS_1000_150_1_2018_MM_mass_type4']
    type4.Add(plot['HSS_1000_150_10_2018_MM_mass_type4'])
    type4.Add(plot['HSS_1000_150_100_2018_MM_mass_type4'])
    type4.Add(plot['HSS_1000_150_1000_2018_MM_mass_type4'])
    type4.Add(plot['HSS_1000_150_10000_2018_MM_mass_type4'])
    type1.SetMaximum(1e6)
    type1.SetLineWidth(2)
    type2.SetLineWidth(2)
    type3.SetLineWidth(2)
    type4.SetLineWidth(2)
    canvas = Canvas.Canvas("HSS_1000_150_2018_MM_mass_types", 'png,pdf', 0.43, 0.6, 0.79, 0.75, 1)
    canvas.addHisto(type1,'HIST', 'S#rightarrow#mu#mu [{0}]'.format(type1.GetEntries()), 'l', acolors['1'], True, 0, doOF= False)
    canvas.addHisto(type2,'HIST,SAME', '#mu (from S) + random [{0}]'.format(type2.GetEntries()), 'l', acolors['2'], True, 1, doOF = False)
    canvas.addHisto(type3,'HIST,SAME', 'Different S parent [{0}]'.format(type3.GetEntries()), 'l', acolors['3'], True, 2, doOF = False)
    canvas.addHisto(type4,'HIST,SAME', 'Not gen-matched [{0}]'.format(type4.GetEntries()), 'l', acolors['4'], True, 3, doOF = False)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.4, 0.85, 'H#rightarrowSS#rightarrow 2#mu + X', size = 0.03)
    canvas.addLatex(0.4, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.4, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 1, '', '', outputDir = WORKPATH + 'plots/', inProgress = True)

    ## 2D plots
    r.gStyle.SetPadRightMargin(0.14)
    r.gStyle.SetPadLeftMargin(0.14)
    EE_RECOstatus_type3 = combinePlots([plot['HSS_1000_150_1_2018_EE_RECOstatus_type3'], plot['HSS_1000_150_10_2018_EE_RECOstatus_type3'], plot['HSS_1000_150_100_2018_EE_RECOstatus_type3'], plot['HSS_1000_150_1000_2018_EE_RECOstatus_type3'], plot['HSS_1000_150_10000_2018_EE_RECOstatus_type3']])
    EE_RECOstatus_type3.GetXaxis().SetBinLabel(1, "Fail")
    EE_RECOstatus_type3.GetXaxis().SetBinLabel(2, "Succeed")
    EE_RECOstatus_type3.GetYaxis().SetBinLabel(1, "Fail")
    EE_RECOstatus_type3.GetYaxis().SetBinLabel(2, "Succeed")
    EE_RECOstatus_type3.GetYaxis().SetTitleOffset(1.6)
    canvas = Canvas.Canvas("HSS_1000_150_2018_EE_RECOstatus_type3", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(EE_RECOstatus_type3,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True)

    r.gStyle.SetPadLeftMargin(0.14)
    MM_RECOstatus_type3 = combinePlots([plot['HSS_1000_150_1_2018_MM_RECOstatus_type3'], plot['HSS_1000_150_10_2018_MM_RECOstatus_type3'], plot['HSS_1000_150_100_2018_MM_RECOstatus_type3'], plot['HSS_1000_150_1000_2018_MM_RECOstatus_type3'], plot['HSS_1000_150_10000_2018_MM_RECOstatus_type3']])
    MM_RECOstatus_type3.GetXaxis().SetBinLabel(1, "Fail")
    MM_RECOstatus_type3.GetXaxis().SetBinLabel(2, "Succeed")
    MM_RECOstatus_type3.GetYaxis().SetBinLabel(1, "Fail")
    MM_RECOstatus_type3.GetYaxis().SetBinLabel(2, "Succeed")
    MM_RECOstatus_type3.GetYaxis().SetTitleOffset(1.6)
    canvas = Canvas.Canvas("HSS_1000_150_2018_MM_RECOstatus_type3", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 650, hh = 600)
    canvas.addHisto(MM_RECOstatus_type3,'COLZ,TEXT', '', '', '', True, 0)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', inProgress = False, is2d = True)
