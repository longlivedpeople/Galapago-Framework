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


def getIPECAL(P, vertex, q, fBz = 3.8, rECAL = 1.290, zECAL = 3.170):

    # 1. Init of particle parameters
    c_light = 2.99792458e8
    gammam = P.E()*1e9/c_light**2 # eV/c^2
    omega = q * fBz / (gammam);
    r = P.Pt() / (q * fBz) * 1.0e9/c_light;

    phi_0 = TMath.ATan2(P.Py(), P.Px()) # in [-pi, pi]

    x_v = vertex.X()*1e-2
    y_v = vertex.Y()*1e-2
    z_v = vertex.Z()*1e-2

    # 2. helix axis coordinates
    x_c = x_v + r*TMath.Sin(phi_0)
    y_c = y_v - r*TMath.Cos(phi_0)
    r_c = TMath.Hypot(x_c, y_c)    
    phi_c = TMath.ATan2(y_c, x_c)
    phi = phi_c
    if (x_c < 0.0):
        phi += TMath.Pi()
    
    rcu = TMath.Abs(r)
    rc2 = r_c*r_c

    # 3. time evaluation
    t_r = 0.0; # in [ns]
    sign_pz = 1 if (P.Pz() > 0.0) else -1;

    if (P.Pz() == 0.0):
        t_z = 1.0e99
    else:
        t_z = gammam / (P.Pz()*1.0e9/c_light) * (-z_v + zECAL*sign_pz);    

    if (r_c + abs(r) < rECAL):
        t = t_z;
    else:
        asinrho = TMath.ASin((rECAL*rECAL - r_c*r_c - r*r) / (2*abs(r)*r_c));
        #asinrho = math.asin((rECAL*rECAL - r_c*r_c - r*r) / (2*abs(r)*r_c));
        delta = phi_c - phi;
        if(delta < -math.pi): delta += 2*TMath.Pi();
        if(delta >  math.pi): delta -= 2*TMath.Pi();
        t1 = (delta + asinrho) / omega;
        t2 = (delta + TMath.Pi() - asinrho) / omega;
        t3 = (delta + TMath.Pi() + asinrho) / omega;
        t4 = (delta - asinrho) / omega;
        t5 = (delta - TMath.Pi() - asinrho) / omega;
        t6 = (delta - TMath.Pi() + asinrho) / omega;

        if(t1 < 0.0): t1 = 1.0e99;
        if(t2 < 0.0): t2 = 1.0e99;
        if(t3 < 0.0): t3 = 1.0e99;
        if(t4 < 0.0): t4 = 1.0e99;
        if(t5 < 0.0): t5 = 1.0e99;
        if(t6 < 0.0): t6 = 1.0e99;

        t_ra = TMath.Min(t1, TMath.Min(t2, t3))
        t_rb = TMath.Min(t4, TMath.Min(t5, t6))
        t_r = TMath.Min(t_ra, t_rb)
        t = TMath.Min(t_r, t_z)

    # 4. position in terms of x(t), y(t), z(t)
    x_t = x_c + r * TMath.Sin(omega * t - phi_0);
    y_t = y_c + r * TMath.Cos(omega * t - phi_0);
    z_t = z_v + P.Pz()*1.0e9 / c_light / gammam * t;
    
    IP = TVector3()
    IP.SetXYZ(x_t, y_t, z_t);

    return IP

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
    eta_binning = np.array([-2.0, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.0])
    print(len(Lxy_binning))
    plot = {}
    for signal in HSS_Signals:

        ## Efficiencies
        plot[signal + '_trackEff']        = copy.deepcopy(r.TEfficiency(signal + '_trackEff', ";Counts;Efficiency", 1, 0, 1))    
        plot[signal + '_trackEff_Lxy']    = copy.deepcopy(r.TEfficiency(signal + '_trackEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", len(Lxy_binning) -1, Lxy_binning)) # in cm    
        plot[signal + '_trackEff_pt']     = copy.deepcopy(r.TEfficiency(signal + '_trackEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 25, 0, 100)) # in cm    
        plot[signal + '_scEff']           = copy.deepcopy(r.TEfficiency(signal + '_scEff', ";Counts;Efficiency", 1, 0, 1))    
        plot[signal + '_scEff_Lxy']       = copy.deepcopy(r.TEfficiency(signal + '_scEff_Lxy', ";Generated LLP decay radius (cm);Reconstruction efficiency", len(Lxy_binning) -1, Lxy_binning)) # in cm    
        plot[signal + '_scEff_pt']        = copy.deepcopy(r.TEfficiency(signal + '_scEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 25, 0, 100)) # in cm    
        plot[signal + '_electronEff']     = copy.deepcopy(r.TEfficiency(signal + '_electronEff', ";Counts;Reconstruction efficiency", 1, 0, 1))    
        plot[signal + '_electronEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_electronEff_Lxy', ";Generated LLP decay radius (cm); Reconstruction efficiency", len(Lxy_binning) -1, Lxy_binning)) # in cm    
        plot[signal + '_electronEff_pt']  = copy.deepcopy(r.TEfficiency(signal + '_electronEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 25, 0, 100)) # in cm    
        plot[signal + '_electronEff_oversc']     = copy.deepcopy(r.TEfficiency(signal + '_electronEff_oversc', ";Counts;Efficiency", 1, 0, 1))    
        plot[signal + '_electronEff_oversc_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_electronEff_oversc_Lxy', ";Generated LLP decay radius (cm);Displaced electron efficiency", 14, 0, 70)) # in cm    
        plot[signal + '_electronEff_oversc_pt']  = copy.deepcopy(r.TEfficiency(signal + '_electronEff_oversc_pt', ";Generated electron p_{T} (GeV);Displaced electron efficiency", 14, 0, 70)) # in cm    
        plot[signal + '_cmsEff_Lxy']      = copy.deepcopy(r.TEfficiency(signal + '_cmsEff_Lxy', ";Generated LLP decay radius (cm); Reconstruction efficiency", len(Lxy_binning) -1, Lxy_binning)) # in cm    
        plot[signal + '_matchEff']        = copy.deepcopy(r.TEfficiency(signal + '_matchEff', ";Counts;Efficiency", 1, 0, 1))    
        plot[signal + '_matchEff_Lxy']    = copy.deepcopy(r.TEfficiency(signal + '_matchEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70)) # in cm    
        plot[signal + '_matchEff_pt']     = copy.deepcopy(r.TEfficiency(signal + '_matchEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 14, 0, 70)) # in cm    

        ## Matching plots
        plot[signal + '_track_dEta_vs_Lxy'] = copy.deepcopy(r.TH2F(signal + '_track_dEta_vs_Lxy', ';Generated LLP decay radius (cm);#Delta#eta(track, gen)', 50, 0, 100, 30, -1, 1))
        plot[signal + '_sc_dEta_vs_Lxy'] = copy.deepcopy(r.TH2F(signal + '_sc_dEta_vs_Lxy', ';Generated LLP decay radius (cm);#Delta#eta(sc, gen)', 50, 0, 100, 30, -1, 1))
        plot[signal + '_track_dPhi_vs_Lxy'] = copy.deepcopy(r.TH2F(signal + '_track_dPhi_vs_Lxy', ';Generated LLP decay radius (cm);#Delta#phi(track, gen)', 50, 0, 100, 30, -1, 1))
        plot[signal + '_sc_dPhi_vs_Lxy'] = copy.deepcopy(r.TH2F(signal + '_sc_dPhi_vs_Lxy', ';Generated LLP decay radius (cm);#Delta#phi(sc, gen)', 50, 0, 100, 30, -1, 1))
        plot[signal + '_track_dR_vs_Lxy'] = copy.deepcopy(r.TH2F(signal + '_track_dR_vs_Lxy', ';Generated LLP decay radius (cm);#DeltaR(track, gen)', 50, 0, 100, 20, 0, 1))
        plot[signal + '_sc_dR_vs_Lxy'] = copy.deepcopy(r.TH2F(signal + '_sc_dR_vs_Lxy', ';Generated LLP decay radius (cm);#DeltaR(sc, gen)', 50, 0, 100, 20, 0, 1))

        ## Resolution and electron distribution plots
        plot[signal + '_electron_dR']     = copy.deepcopy(r.TH1F(signal + '_electron_dR', ";Matching distance #DeltaR(track, SC);Normalized displaced electron yield", 30, 0, 0.3))
        plot[signal + '_electron_ptres']  = copy.deepcopy(r.TH1F(signal + '_electron_ptres', ";(p_{T}^{track} - p_{T}^{gen})/p_{T}^{gen};Normalized displaced electron yield", 51, -0.2, 0.2)) 
        plot[signal + '_electron_etres']  = copy.deepcopy(r.TH1F(signal + '_electron_etres', ";(E_{T}^{sc} - p_{T}^{gen})/p_{T}^{gen};Normalized displaced electron yield", 51, -0.2, 0.2)) 
        plot[signal + '_electron_dxyres']  = copy.deepcopy(r.TH1F(signal + '_electron_dxyres', ";(|d_{xy}|^{track} - |d_{xy}|^{gen})/|d_{xy}|^{gen};Normalized displaced electron yield", 51, -0.2, 0.2)) 
        plot[signal + '_electron_pt_et']  = copy.deepcopy(r.TH2F(signal + '_electron_pt_et', ";Displaced electron track p_{T} (GeV);Displaced electron SC E_{T} (GeV)", 56, 20, 300, 56, 20, 300)) 

        ## Others
        plot[signal + '_clustersAssocToTrack'] = copy.deepcopy(r.TH1F(signal + '_clustersAssocToTrack', ";Number of clusters matched to same track;Number of tracks", 3, 0, 3))
        plot[signal + '_genpt_scet'] = copy.deepcopy(r.TH2F(signal + '_genpt_scet', ';Generated electron p_{T} (GeV);Reconstructed E_{T}', 56, 20, 300, 56, 20, 300))
        plot[signal + '_genpt_scet_passSigma'] = copy.deepcopy(r.TH2F(signal + '_genpt_scet_passSigma', ';Generated electron p_{T} (GeV);Reconstructed E_{T}', 56, 20, 300, 56, 20, 300))
        plot[signal + '_genpt_scet_passHoE'] = copy.deepcopy(r.TH2F(signal + '_genpt_scet_passHoE', ';Generated electron p_{T} (GeV);Reconstructed E_{T}', 56, 20, 300, 56, 20, 300))

        ## ID efficiency plots (over electron efficiency)
        plot[signal + '_hitsEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_hitsEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))  
        plot[signal + '_hpEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_hpEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    
        plot[signal + '_hoeEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_hoeEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))
        plot[signal + '_sigmaEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_sigmaEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    
        plot[signal + '_IdEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_IdEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    
        plot[signal + '_IdNohitEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohitEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    
        plot[signal + '_IdNohpEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohpEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    
        plot[signal + '_IdNohoeEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohoeEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    
        plot[signal + '_IdNosigmaEff_Lxy'] = copy.deepcopy(r.TEfficiency(signal + '_IdNosigmaEff_Lxy', ";Generated LLP decay radius (cm);Efficiency", 14, 0, 70))    

        plot[signal + '_hitsEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_hitsEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))  
        plot[signal + '_hpEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_hpEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    
        plot[signal + '_hoeEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_hoeEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))
        plot[signal + '_sigmaEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_sigmaEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    
        plot[signal + '_IdEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_IdEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    
        plot[signal + '_IdNohitEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohitEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    
        plot[signal + '_IdNohpEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohpEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    
        plot[signal + '_IdNohoeEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohoeEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    
        plot[signal + '_IdNosigmaEff_pt'] = copy.deepcopy(r.TEfficiency(signal + '_IdNosigmaEff_pt', ";Generated electron p_{T} (GeV);Efficiency", 60, 0, 300))    

        plot[signal + '_hitsEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_hitsEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))  
        plot[signal + '_hpEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_hpEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    
        plot[signal + '_hoeEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_hoeEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))
        plot[signal + '_sigmaEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_sigmaEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    
        plot[signal + '_IdEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_IdEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    
        plot[signal + '_IdNohitEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohitEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    
        plot[signal + '_IdNohpEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohpEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    
        plot[signal + '_IdNohoeEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_IdNohoeEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    
        plot[signal + '_IdNosigmaEff_eta'] = copy.deepcopy(r.TEfficiency(signal + '_IdNosigmaEff_eta', ";Generated electron #eta;Efficiency", len(eta_binning) -1, eta_binning))    




        for p in plot.keys():
            r.SetOwnership(plot[p], 0)

        tree = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'dat/signals_'+year+'UL_displacedElectron.dat', [signal], 'MC'), name = year, isdata = 0 )

        print("Processing " + signal)

        for b in tree.blocks: 
            for s in b.samples: 
                for t in s.ttrees:
                    for e,ev in enumerate(t):

                        ## Loop over electrons:
                        for i in range(0, ev.nGenLepton):

                            if ( abs(ev.GenLeptonSel_pdgId[i]) != 11 ):
                                continue 
                            if ( abs(ev.GenLeptonSel_eta[i]) > 2.0 ):
                                continue 
                            if ( abs(ev.GenLeptonSel_eta[i]) > 1.4 and abs(ev.GenLeptonSel_eta[i]) < 1.6 ):
                                continue 
                            if not ( ev.GenLeptonSel_fromHardProcessFinalState[i]): 
                                continue 

                            ## Skip electrons that radiated from now
                            radiated = True
                            pt = -999
                            eta = -999
                            phi = -999
                            for j in range(0, ev.nHardProcessParticle):
                                if ev.HardProcessParticle_pdgId[j] != ev.GenLeptonSel_pdgId[i]: continue
                                if ev.HardProcessParticle_vx[j] == ev.GenLeptonSel_vx[i] and ev.HardProcessParticle_vy[j] == ev.GenLeptonSel_vy[i] and ev.HardProcessParticle_vz[j] == ev.GenLeptonSel_vz[i]:
                                    pt = ev.HardProcessParticle_pt[j]
                                    eta = ev.HardProcessParticle_eta[j]
                                    phi = ev.HardProcessParticle_phi[j]
                                    break

                            vertex = TVector3(ev.GenLeptonSel_vx[i], ev.GenLeptonSel_vy[i], ev.GenLeptonSel_vz[i])
                            charge = -math.copysign(1, ev.GenLeptonSel_pdgId[i]) 
                            vector = TLorentzVector()
                            vector.SetPtEtaPhiM(ev.GenLeptonSel_pt[i], ev.GenLeptonSel_eta[i], ev.GenLeptonSel_phi[i], 0.510e-3)
                            #vector.SetPtEtaPhiE(ev.HardProcessParticle_pt[hard], ev.HardProcessParticle_eta[hard], ev.HardProcessParticle_phi[hard], ev.HardProcessParticle_E[hard])

                            Lxy = math.sqrt((ev.GenLeptonSel_vx[i])**2 + (ev.GenLeptonSel_vy[i])**2)
                            dxy = -(ev.GenLeptonSel_vx[i] - ev.PV_vx)*math.sin(phi) + (ev.GenLeptonSel_vy[i] - ev.PV_vy)*math.cos(phi)
                            absdxy = abs(dxy)


                            # tracker acceptance 
                            if Lxy > 70:
                                continue

                            IPvector = getIPECAL(vector, vertex, charge)


                            gvec = TVector3()
                            gvecextr = TVector3()
                            gvec.SetPtEtaPhi(pt, eta, phi)

                            ## Isotrack matching
                            mindR_track = 9999.
                            trackIdx = -1
                            matchedtracks = []
                            for j in range(0, ev.nIsoTrack):
                                ivec = TVector3()
                                ivec.SetPtEtaPhi(ev.IsoTrackSel_pt[j], ev.IsoTrackSel_eta[j], ev.IsoTrackSel_phi[j])
                                dR = gvec.DeltaR(ivec)
                                if dR < 0.25:
                                    matchedtracks.append(j)
                                if dR < mindR_track:
                                    mindR_track = dR
                                    trackIdx = j
                            if trackIdx > -1 and pt > 20:
                                 plot[signal + '_track_dEta_vs_Lxy'].Fill(Lxy, ev.IsoTrackSel_eta[trackIdx] - ev.GenLeptonSel_eta[i])
                                 plot[signal + '_track_dPhi_vs_Lxy'].Fill(Lxy, ev.IsoTrackSel_phi[trackIdx] - ev.GenLeptonSel_phi[i])
                                 plot[signal + '_track_dR_vs_Lxy'].Fill(Lxy, mindR_track)
                            trackMatched = mindR_track < 0.25

                            ## Count how many superclusters can be matched to the isotrack
                            if trackMatched:
                                nclusters = 0
                                tvec = TVector3()
                                tvec.SetPtEtaPhi(ev.IsoTrackSel_pt[trackIdx], ev.IsoTrackSel_eta[trackIdx] + ev.IsoTrackSel_etaExtra[trackIdx], ev.IsoTrackSel_phi[trackIdx] + ev.IsoTrackSel_phiExtra[trackIdx])
                                for j in range(0, ev.nPhoton):
                                    cvec = TVector3()
                                    cvec.SetPtEtaPhi(ev.PhotonSel_et[j], ev.PhotonSel_eta[j], ev.PhotonSel_phi[j])
                                    if tvec.DeltaR(cvec) < 0.1:
                                        nclusters += 1
                                plot[signal + '_clustersAssocToTrack'].Fill(nclusters)

                            ## Supercluster matching
                            mindR_sc = 9999.
                            scIdx = -1
                            matchedsc = []
                            for j in range(0, ev.nPhoton):
                                cvec = TVector3()
                                cvec.SetPtEtaPhi(ev.PhotonSel_et[j], ev.PhotonSel_eta[j], ev.PhotonSel_phi[j])
                                #dR = gvec.DeltaR(cvec)
                                #dR = gvecextr.DeltaR(cvec)
                                dR = IPvector.DeltaR(cvec)
                                if dR < 0.35:
                                    matchedsc.append(j)
                                if dR < mindR_sc:
                                    mindR_sc = dR
                                    scIdx = j
                            if scIdx > -1 and pt > 20:
                                 plot[signal + '_sc_dEta_vs_Lxy'].Fill(Lxy, ev.PhotonSel_eta[scIdx] - ev.GenLeptonSel_eta[i])
                                 plot[signal + '_sc_dPhi_vs_Lxy'].Fill(Lxy, ev.PhotonSel_phi[scIdx] - ev.GenLeptonSel_phi[i])
                                 plot[signal + '_sc_dR_vs_Lxy'].Fill(Lxy, mindR_sc)

                            scMatched = mindR_sc < 0.35

                            ## CMS Electron matching
                            mindR_cms = 9999.
                            cmsIdx = -1
                            for j in range(0, ev.nElectron):
                                evec = TVector3()
                                evec.SetPtEtaPhi(ev.ElectronSel_pt[j], ev.ElectronSel_eta[j], ev.ElectronSel_phi[j])
                                dR = gvec.DeltaR(evec)
                                if dR < mindR_cms:
                                    mindR_cms = dR
                                    cmsIdx = j
                            cmsMatched = mindR_cms < 0.25


                            ## Displaced Electron matching and ID
                            eMatched  = False
                            passHits  = False
                            passPt    = False
                            passHP    = False
                            passHoE   = False
                            passSigma = False
                            passID    = False
                            ptres = -99
                            etres = -99
                            recoet = -99
                            for j in range(0, ev.nElectronCandidate):
                                if ev.ElectronCandidate_isotrackIdx[j] in matchedtracks and ev.ElectronCandidate_photonIdx[j] in matchedsc:
                                    plot[signal + '_electron_dR'].Fill(ev.ElectronCandidate_dR[j])
                                    if ev.ElectronCandidate_dR[j] > 0.1:
                                        continue
                                    eMatched = True
                                
                                    ## Resolution:
                                    ptres = (ev.ElectronCandidate_pt[j] - pt)/pt
                                    etres = (ev.ElectronCandidate_et[j] - pt)/pt
                                    recopt = ev.ElectronCandidate_pt[j]
                                    recoet = ev.ElectronCandidate_et[j]
                                    dxyres = (abs(ev.ElectronCandidate_dxy_PV[j]) - absdxy)/absdxy

                                    ## Evaluate ID:
                                    if ev.IsoTrackSel_numberOfValidTrackerHits[ev.ElectronCandidate_isotrackIdx[j]] > 5:
                                        passHits = True
                                    if ev.IsoTrackSel_isHighPurityTrack[ev.ElectronCandidate_isotrackIdx[j]]:
                                        passHP = True
                                    if ev.PhotonSel_hadronicOverEm[ev.ElectronCandidate_photonIdx[j]] < 0.05:
                                        passHoE = True
                                    if (ev.PhotonSel_full5x5_sigmaIetaIeta[ev.ElectronCandidate_photonIdx[j]] < 0.0112 and abs(ev.PhotonSel_eta[ev.ElectronCandidate_photonIdx[j]]) < 1.4442) or (ev.PhotonSel_full5x5_sigmaIetaIeta[ev.ElectronCandidate_photonIdx[j]] < 0.0425 and abs(ev.PhotonSel_eta[ev.ElectronCandidate_photonIdx[j]]) > 1.566):
                                        passSigma = True
                                    if passHits and passHoE and passSigma: # passHP removed
                                        passID = True

                                    break


                            ## Efficiencies
                            plot[signal + '_trackEff'].Fill(trackMatched, 0)
                            plot[signal + '_trackEff_pt'].Fill(trackMatched, pt)
                            plot[signal + '_scEff'].Fill(scMatched, 0)
                            plot[signal + '_scEff_pt'].Fill(scMatched, pt)
                            plot[signal + '_electronEff'].Fill(eMatched, 0)
                            plot[signal + '_electronEff_pt'].Fill(eMatched, pt)
                            if trackMatched and scMatched:
                                plot[signal + '_matchEff_pt'].Fill(eMatched, pt)
                            if trackMatched:
                                plot[signal + '_electronEff_oversc'].Fill(eMatched, 0)
                                plot[signal + '_electronEff_oversc_pt'].Fill(eMatched, pt)

                            if pt > 15:
                                plot[signal + '_trackEff_Lxy'].Fill(trackMatched, Lxy)
                                plot[signal + '_scEff_Lxy'].Fill(scMatched, Lxy)
                                plot[signal + '_electronEff_Lxy'].Fill(eMatched, Lxy)
                                plot[signal + '_cmsEff_Lxy'].Fill(cmsMatched, Lxy)
                                if trackMatched and scMatched:
                                    plot[signal + '_matchEff_Lxy'].Fill(eMatched, Lxy)
                                if trackMatched:
                                    plot[signal + '_electronEff_oversc_Lxy'].Fill(eMatched, Lxy)
                                if eMatched and recopt > 15:
                                    plot[signal + '_electron_ptres'].Fill(ptres)
                                    plot[signal + '_electron_etres'].Fill(etres)
                                    plot[signal + '_electron_dxyres'].Fill(dxyres)
                                    plot[signal + '_electron_pt_et'].Fill(recopt, recoet)
                                    plot[signal + '_hitsEff_Lxy'].Fill(passHits, Lxy)
                                    plot[signal + '_hpEff_Lxy'].Fill(passHP, Lxy)
                                    plot[signal + '_hoeEff_Lxy'].Fill(passHoE, Lxy)
                                    plot[signal + '_sigmaEff_Lxy'].Fill(passSigma, Lxy)
                                    plot[signal + '_IdEff_Lxy'].Fill(passID, Lxy)
                                    plot[signal + '_IdNohitEff_Lxy'].Fill(passHP and passHoE and passSigma, Lxy)
                                    plot[signal + '_IdNohpEff_Lxy'].Fill(passHits and passHoE and passSigma, Lxy)
                                    plot[signal + '_IdNohoeEff_Lxy'].Fill(passHP and passHits and passSigma, Lxy)
                                    plot[signal + '_IdNosigmaEff_Lxy'].Fill(passHits and passHP and passHoE, Lxy)
                                    plot[signal + '_hitsEff_pt'].Fill(passHits, pt)
                                    plot[signal + '_hpEff_pt'].Fill(passHP, pt)
                                    plot[signal + '_hoeEff_pt'].Fill(passHoE, pt)
                                    plot[signal + '_sigmaEff_pt'].Fill(passSigma, pt)
                                    plot[signal + '_IdEff_pt'].Fill(passID, pt)
                                    plot[signal + '_IdNohitEff_pt'].Fill(passHP and passHoE and passSigma, pt)
                                    plot[signal + '_IdNohpEff_pt'].Fill(passHits and passHoE and passSigma, pt)
                                    plot[signal + '_IdNohoeEff_pt'].Fill(passHP and passHits and passSigma, pt)
                                    plot[signal + '_IdNosigmaEff_pt'].Fill(passHits and passHP and passHoE, pt)
                                    plot[signal + '_hitsEff_eta'].Fill(passHits, eta)
                                    plot[signal + '_hpEff_eta'].Fill(passHP, eta)
                                    plot[signal + '_hoeEff_eta'].Fill(passHoE, eta)
                                    plot[signal + '_sigmaEff_eta'].Fill(passSigma, eta)
                                    plot[signal + '_IdEff_eta'].Fill(passID, eta)
                                    plot[signal + '_IdNohitEff_eta'].Fill(passHP and passHoE and passSigma, eta)
                                    plot[signal + '_IdNohpEff_eta'].Fill(passHits and passHoE and passSigma, eta)
                                    plot[signal + '_IdNohoeEff_eta'].Fill(passHP and passHits and passSigma, eta)
                                    plot[signal + '_IdNosigmaEff_eta'].Fill(passHits and passHP and passHoE, eta)
                                    plot[signal + '_genpt_scet'].Fill(pt, recoet)
                                    if passSigma: plot[signal + '_genpt_scet_passSigma'].Fill(pt, recoet)
                                    if passHoE: plot[signal + '_genpt_scet_passHoE'].Fill(pt, recoet)



    ##################################################################################################
    ### H->SS decay radius plots
    canvas = Canvas.Canvas("HSS_1000_150_1_2018_trackEff_pt", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(plot['HSS_1000_150_1_2018_trackEff_pt'],'AP', 'Isotrack efficiency', 'p', r.kBlue, True, 0, marker = 24)
    canvas.addHisto(plot['HSS_1000_150_1_2018_scEff_pt'],'AP,SAME', 'Supercluster efficiency', 'p', r.kRed, True, 1, marker = 25)
    canvas.addHisto(plot['HSS_1000_150_1_2018_electronEff_pt'],'AP,SAME', 'Displaced electron efficiency', 'p', r.kBlack, True, 2, marker = 20)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)

    ### Efficiency vs Lxy
    trackEff = plot['HSS_1000_150_1_2018_trackEff_Lxy']
    trackEff.Add(plot['HSS_1000_150_10_2018_trackEff_Lxy'])
    trackEff.Add(plot['HSS_1000_150_100_2018_trackEff_Lxy'])
    trackEff.Add(plot['HSS_1000_150_1000_2018_trackEff_Lxy'])
    trackEff.Add(plot['HSS_1000_150_10000_2018_trackEff_Lxy'])
    scEff = plot['HSS_1000_150_1_2018_scEff_Lxy']
    scEff.Add(plot['HSS_1000_150_10_2018_scEff_Lxy'])
    scEff.Add(plot['HSS_1000_150_100_2018_scEff_Lxy'])
    scEff.Add(plot['HSS_1000_150_1000_2018_scEff_Lxy'])
    scEff.Add(plot['HSS_1000_150_10000_2018_scEff_Lxy'])
    electronEff = plot['HSS_1000_150_1_2018_electronEff_Lxy']
    electronEff.Add(plot['HSS_1000_150_10_2018_electronEff_Lxy'])
    electronEff.Add(plot['HSS_1000_150_100_2018_electronEff_Lxy'])
    electronEff.Add(plot['HSS_1000_150_1000_2018_electronEff_Lxy'])
    electronEff.Add(plot['HSS_1000_150_10000_2018_electronEff_Lxy'])
    cmsEff = plot['HSS_1000_150_1_2018_cmsEff_Lxy']
    cmsEff.Add(plot['HSS_1000_150_10_2018_cmsEff_Lxy'])
    cmsEff.Add(plot['HSS_1000_150_100_2018_cmsEff_Lxy'])
    cmsEff.Add(plot['HSS_1000_150_1000_2018_cmsEff_Lxy'])
    cmsEff.Add(plot['HSS_1000_150_10000_2018_cmsEff_Lxy'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_Eff_Lxy", 'png,pdf', 0.4, 0.5, 0.8, 0.65, 1)
    canvas.addHisto(trackEff,'AP', 'Isotrack efficiency', 'p', r.kBlue, True, 0, marker = 24)
    canvas.addHisto(scEff,'P,SAME', 'Supercluster efficiency', 'p', r.kRed, True, 1, marker = 25)
    canvas.addHisto(electronEff,'P,SAME', 'Displaced electron efficiency', 'p', r.kBlack, True, 2, marker = 20)
    canvas.addHisto(cmsEff,'P,SAME', 'CMS electron efficiency', 'p', r.kGreen+1, True, 3, marker = 20)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, is2d = True)

    ## Electron over track efficiency
    electronEff = plot['HSS_1000_150_1_2018_electronEff_oversc_Lxy']
    electronEff.Add(plot['HSS_1000_150_10_2018_electronEff_oversc_Lxy'])
    electronEff.Add(plot['HSS_1000_150_100_2018_electronEff_oversc_Lxy'])
    electronEff.Add(plot['HSS_1000_150_1000_2018_electronEff_oversc_Lxy'])
    electronEff.Add(plot['HSS_1000_150_10000_2018_electronEff_oversc_Lxy'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_electronEff_oversc_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(electronEff,'AP', 'Displaced electron efficiency', 'p', r.kBlack, True, 0, marker = 20)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.2, 0.4, 'H#rightarrowSS#rightarrow 2e + X', size = 0.033)
    canvas.addLatex(0.2, 0.35, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.033)
    canvas.addLatex(0.2, 0.3, '(all lifetimes combined)', size = 0.033)
    canvas.save(0, 0, 0, '', '', ymin = 0.6, ymax = 1.1, outputDir = WORKPATH + 'plots/', isPrivate = True, is2d = True)

    electronEff = plot['HSS_1000_150_1_2018_electronEff_oversc_pt']
    electronEff.Add(plot['HSS_1000_150_10_2018_electronEff_oversc_pt'])
    electronEff.Add(plot['HSS_1000_150_100_2018_electronEff_oversc_pt'])
    electronEff.Add(plot['HSS_1000_150_1000_2018_electronEff_oversc_pt'])
    electronEff.Add(plot['HSS_1000_150_10000_2018_electronEff_oversc_pt'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_electronEff_oversc_pt", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(electronEff,'AP', 'Displaced electron efficiency', 'p', r.kBlack, True, 0, marker = 20)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)


    ## Matching efficiency
    matchEff = plot['HSS_1000_150_1_2018_matchEff_Lxy']
    matchEff.Add(plot['HSS_1000_150_10_2018_matchEff_Lxy'])
    matchEff.Add(plot['HSS_1000_150_100_2018_matchEff_Lxy'])
    matchEff.Add(plot['HSS_1000_150_1000_2018_matchEff_Lxy'])
    matchEff.Add(plot['HSS_1000_150_10000_2018_matchEff_Lxy'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_matchEff_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(matchEff,'AP', 'Displaced electron efficiency', 'p', r.kBlack, True, 0, marker = 20)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)

    matchEff = plot['HSS_1000_150_1_2018_matchEff_pt']
    matchEff.Add(plot['HSS_1000_150_10_2018_matchEff_pt'])
    matchEff.Add(plot['HSS_1000_150_100_2018_matchEff_pt'])
    matchEff.Add(plot['HSS_1000_150_1000_2018_matchEff_pt'])
    matchEff.Add(plot['HSS_1000_150_10000_2018_matchEff_pt'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_matchEff_pt", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(matchEff,'AP', 'Displaced electron efficiency', 'p', r.kBlack, True, 0, marker = 20)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)

    ## ID efficiency
    # The single cut efficiency (Lxy):
    hitsEff = plot['HSS_1000_150_1_2018_hitsEff_Lxy']
    hitsEff.Add(plot['HSS_1000_150_10_2018_hitsEff_Lxy'])
    hitsEff.Add(plot['HSS_1000_150_100_2018_hitsEff_Lxy'])
    hitsEff.Add(plot['HSS_1000_150_1000_2018_hitsEff_Lxy'])
    hitsEff.Add(plot['HSS_1000_150_10000_2018_hitsEff_Lxy'])

    hpEff = plot['HSS_1000_150_1_2018_hpEff_Lxy']
    hpEff.Add(plot['HSS_1000_150_10_2018_hpEff_Lxy'])
    hpEff.Add(plot['HSS_1000_150_100_2018_hpEff_Lxy'])
    hpEff.Add(plot['HSS_1000_150_1000_2018_hpEff_Lxy'])
    hpEff.Add(plot['HSS_1000_150_10000_2018_hpEff_Lxy'])

    hoeEff = plot['HSS_1000_150_1_2018_hoeEff_Lxy']
    hoeEff.Add(plot['HSS_1000_150_10_2018_hoeEff_Lxy'])
    hoeEff.Add(plot['HSS_1000_150_100_2018_hoeEff_Lxy'])
    hoeEff.Add(plot['HSS_1000_150_1000_2018_hoeEff_Lxy'])
    hoeEff.Add(plot['HSS_1000_150_10000_2018_hoeEff_Lxy'])

    sigmaEff = plot['HSS_1000_150_1_2018_sigmaEff_Lxy']
    sigmaEff.Add(plot['HSS_1000_150_10_2018_sigmaEff_Lxy'])
    sigmaEff.Add(plot['HSS_1000_150_100_2018_sigmaEff_Lxy'])
    sigmaEff.Add(plot['HSS_1000_150_1000_2018_sigmaEff_Lxy'])
    sigmaEff.Add(plot['HSS_1000_150_10000_2018_sigmaEff_Lxy'])

    IdEff = plot['HSS_1000_150_1_2018_IdEff_Lxy']
    IdEff.Add(plot['HSS_1000_150_10_2018_IdEff_Lxy'])
    IdEff.Add(plot['HSS_1000_150_100_2018_IdEff_Lxy'])
    IdEff.Add(plot['HSS_1000_150_1000_2018_IdEff_Lxy'])
    IdEff.Add(plot['HSS_1000_150_10000_2018_IdEff_Lxy'])

    canvas = Canvas.Canvas("HSS_1000_150_2018_IdEff_Lxy", 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1)
    canvas.addHisto(hitsEff,'AP', 'Track hit selection', 'p', gcolors['blue'], True, 0, marker = 24)
    #canvas.addHisto(hpEff,'P,SAME', 'Track high purity requirement', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(hoeEff,'P,SAME', 'SC H/E selection', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(sigmaEff,'P,SAME', 'SC #sigma_{i#etai#eta} selection', 'p', gcolors['green'], True, 2, marker = 24)
    canvas.addHisto(IdEff,'P,SAME', 'All requirements', 'p', r.kBlack, True, 3, marker = 24)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.6, ymax = 1.1, outputDir = WORKPATH + 'plots/', isPrivate = True)


    # The single cut efficiency (pt):
    hitsEff = plot['HSS_1000_150_1_2018_hitsEff_pt']
    hitsEff.Add(plot['HSS_1000_150_10_2018_hitsEff_pt'])
    hitsEff.Add(plot['HSS_1000_150_100_2018_hitsEff_pt'])
    hitsEff.Add(plot['HSS_1000_150_1000_2018_hitsEff_pt'])
    hitsEff.Add(plot['HSS_1000_150_10000_2018_hitsEff_pt'])

    hpEff = plot['HSS_1000_150_1_2018_hpEff_pt']
    hpEff.Add(plot['HSS_1000_150_10_2018_hpEff_pt'])
    hpEff.Add(plot['HSS_1000_150_100_2018_hpEff_pt'])
    hpEff.Add(plot['HSS_1000_150_1000_2018_hpEff_pt'])
    hpEff.Add(plot['HSS_1000_150_10000_2018_hpEff_pt'])

    hoeEff = plot['HSS_1000_150_1_2018_hoeEff_pt']
    hoeEff.Add(plot['HSS_1000_150_10_2018_hoeEff_pt'])
    hoeEff.Add(plot['HSS_1000_150_100_2018_hoeEff_pt'])
    hoeEff.Add(plot['HSS_1000_150_1000_2018_hoeEff_pt'])
    hoeEff.Add(plot['HSS_1000_150_10000_2018_hoeEff_pt'])

    sigmaEff = plot['HSS_1000_150_1_2018_sigmaEff_pt']
    sigmaEff.Add(plot['HSS_1000_150_10_2018_sigmaEff_pt'])
    sigmaEff.Add(plot['HSS_1000_150_100_2018_sigmaEff_pt'])
    sigmaEff.Add(plot['HSS_1000_150_1000_2018_sigmaEff_pt'])
    sigmaEff.Add(plot['HSS_1000_150_10000_2018_sigmaEff_pt'])

    IdEff = plot['HSS_1000_150_1_2018_IdEff_pt']
    IdEff.Add(plot['HSS_1000_150_10_2018_IdEff_pt'])
    IdEff.Add(plot['HSS_1000_150_100_2018_IdEff_pt'])
    IdEff.Add(plot['HSS_1000_150_1000_2018_IdEff_pt'])
    IdEff.Add(plot['HSS_1000_150_10000_2018_IdEff_pt'])

    canvas = Canvas.Canvas("HSS_1000_150_2018_IdEff_pt", 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1)
    canvas.addHisto(hitsEff,'AP', 'Track hit selection', 'p', gcolors['blue'], True, 0, marker = 24)
    #canvas.addHisto(hpEff,'AP,SAME', 'Track high purity requirement', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(hoeEff,'P,SAME', 'SC H/E selection', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(sigmaEff,'P,SAME', 'SC #sigma_{i#etai#eta} selection', 'p', gcolors['green'], True, 2, marker = 24)
    canvas.addHisto(IdEff,'P,SAME', 'All requirements', 'p', r.kBlack, True, 3, marker = 24)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.4, ymax = 1.2, outputDir = WORKPATH + 'plots/', isPrivate = True, is2d = True)

    # The single cut efficiency (eta):
    hitsEff = plot['HSS_1000_150_1_2018_hitsEff_eta']
    hitsEff.Add(plot['HSS_1000_150_10_2018_hitsEff_eta'])
    hitsEff.Add(plot['HSS_1000_150_100_2018_hitsEff_eta'])
    hitsEff.Add(plot['HSS_1000_150_1000_2018_hitsEff_eta'])
    hitsEff.Add(plot['HSS_1000_150_10000_2018_hitsEff_eta'])

    hpEff = plot['HSS_1000_150_1_2018_hpEff_eta']
    hpEff.Add(plot['HSS_1000_150_10_2018_hpEff_eta'])
    hpEff.Add(plot['HSS_1000_150_100_2018_hpEff_eta'])
    hpEff.Add(plot['HSS_1000_150_1000_2018_hpEff_eta'])
    hpEff.Add(plot['HSS_1000_150_10000_2018_hpEff_eta'])

    hoeEff = plot['HSS_1000_150_1_2018_hoeEff_eta']
    hoeEff.Add(plot['HSS_1000_150_10_2018_hoeEff_eta'])
    hoeEff.Add(plot['HSS_1000_150_100_2018_hoeEff_eta'])
    hoeEff.Add(plot['HSS_1000_150_1000_2018_hoeEff_eta'])
    hoeEff.Add(plot['HSS_1000_150_10000_2018_hoeEff_eta'])

    sigmaEff = plot['HSS_1000_150_1_2018_sigmaEff_eta']
    sigmaEff.Add(plot['HSS_1000_150_10_2018_sigmaEff_eta'])
    sigmaEff.Add(plot['HSS_1000_150_100_2018_sigmaEff_eta'])
    sigmaEff.Add(plot['HSS_1000_150_1000_2018_sigmaEff_eta'])
    sigmaEff.Add(plot['HSS_1000_150_10000_2018_sigmaEff_eta'])

    IdEff = plot['HSS_1000_150_1_2018_IdEff_eta']
    IdEff.Add(plot['HSS_1000_150_10_2018_IdEff_eta'])
    IdEff.Add(plot['HSS_1000_150_100_2018_IdEff_eta'])
    IdEff.Add(plot['HSS_1000_150_1000_2018_IdEff_eta'])
    IdEff.Add(plot['HSS_1000_150_10000_2018_IdEff_eta'])

    canvas = Canvas.Canvas("HSS_1000_150_2018_IdEff_eta", 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1)
    canvas.addHisto(hitsEff,'AP', 'Track hit selection', 'p', gcolors['blue'], True, 0, marker = 24)
    #canvas.addHisto(hpEff,'AP,SAME', 'Track high purity requirement', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(hoeEff,'P,SAME', 'SC H/E selection', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(sigmaEff,'P,SAME', 'SC #sigma_{i#etai#eta} selection', 'p', gcolors['green'], True, 2, marker = 24)
    canvas.addHisto(IdEff,'P,SAME', 'All requirements', 'p', r.kBlack, True, 3, marker = 24)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.4, ymax = 1.2, outputDir = WORKPATH + 'plots/', isPrivate = True, is2d = True)


    # The Id - cut efficiency (Lxy):
    hitsEff = plot['HSS_1000_150_1_2018_IdNohitEff_Lxy']
    hitsEff.Add(plot['HSS_1000_150_10_2018_IdNohitEff_Lxy'])
    hitsEff.Add(plot['HSS_1000_150_100_2018_IdNohitEff_Lxy'])
    hitsEff.Add(plot['HSS_1000_150_1000_2018_IdNohitEff_Lxy'])
    hitsEff.Add(plot['HSS_1000_150_10000_2018_IdNohitEff_Lxy'])

    hpEff = plot['HSS_1000_150_1_2018_IdNohpEff_Lxy']
    hpEff.Add(plot['HSS_1000_150_10_2018_IdNohpEff_Lxy'])
    hpEff.Add(plot['HSS_1000_150_100_2018_IdNohpEff_Lxy'])
    hpEff.Add(plot['HSS_1000_150_1000_2018_IdNohpEff_Lxy'])
    hpEff.Add(plot['HSS_1000_150_10000_2018_IdNohpEff_Lxy'])

    hoeEff = plot['HSS_1000_150_1_2018_IdNohoeEff_Lxy']
    hoeEff.Add(plot['HSS_1000_150_10_2018_IdNohoeEff_Lxy'])
    hoeEff.Add(plot['HSS_1000_150_100_2018_IdNohoeEff_Lxy'])
    hoeEff.Add(plot['HSS_1000_150_1000_2018_IdNohoeEff_Lxy'])
    hoeEff.Add(plot['HSS_1000_150_10000_2018_IdNohoeEff_Lxy'])

    sigmaEff = plot['HSS_1000_150_1_2018_IdNosigmaEff_Lxy']
    sigmaEff.Add(plot['HSS_1000_150_10_2018_IdNosigmaEff_Lxy'])
    sigmaEff.Add(plot['HSS_1000_150_100_2018_IdNosigmaEff_Lxy'])
    sigmaEff.Add(plot['HSS_1000_150_1000_2018_IdNosigmaEff_Lxy'])
    sigmaEff.Add(plot['HSS_1000_150_10000_2018_IdNosigmaEff_Lxy'])

    IdEff = plot['HSS_1000_150_1_2018_IdEff_Lxy']
    IdEff.Add(plot['HSS_1000_150_10_2018_IdEff_Lxy'])
    IdEff.Add(plot['HSS_1000_150_100_2018_IdEff_Lxy'])
    IdEff.Add(plot['HSS_1000_150_1000_2018_IdEff_Lxy'])
    IdEff.Add(plot['HSS_1000_150_10000_2018_IdEff_Lxy'])

    canvas = Canvas.Canvas("HSS_1000_150_2018_NoIdEff_Lxy", 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1)
    canvas.addHisto(hitsEff,'AP', 'All requirements - hit number', 'p', gcolors['blue'], True, 0, marker = 24)
    canvas.addHisto(hpEff,'AP,SAME', 'All requirements - high purity', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(hoeEff,'AP,SAME', 'All requirements - SC H/E', 'p', gcolors['yellow'], True, 2, marker = 24)
    canvas.addHisto(sigmaEff,'AP,SAME', 'All requirements - SC #sigma_{i#etai#eta}', 'p', gcolors['green'], True, 3, marker = 24)
    canvas.addHisto(IdEff,'AP,SAME', 'All requirements', 'p', r.kBlack, True, 4, marker = 24)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.6, ymax = 1.1, outputDir = WORKPATH + 'plots/', isPrivate = True, is2d = True)


    # The Id - cut efficiency (pt):
    hitsEff = plot['HSS_1000_150_1_2018_IdNohitEff_pt']
    hitsEff.Add(plot['HSS_1000_150_10_2018_IdNohitEff_pt'])
    hitsEff.Add(plot['HSS_1000_150_100_2018_IdNohitEff_pt'])
    hitsEff.Add(plot['HSS_1000_150_1000_2018_IdNohitEff_pt'])
    hitsEff.Add(plot['HSS_1000_150_10000_2018_IdNohitEff_pt'])

    hpEff = plot['HSS_1000_150_1_2018_IdNohpEff_pt']
    hpEff.Add(plot['HSS_1000_150_10_2018_IdNohpEff_pt'])
    hpEff.Add(plot['HSS_1000_150_100_2018_IdNohpEff_pt'])
    hpEff.Add(plot['HSS_1000_150_1000_2018_IdNohpEff_pt'])
    hpEff.Add(plot['HSS_1000_150_10000_2018_IdNohpEff_pt'])

    hoeEff = plot['HSS_1000_150_1_2018_IdNohoeEff_pt']
    hoeEff.Add(plot['HSS_1000_150_10_2018_IdNohoeEff_pt'])
    hoeEff.Add(plot['HSS_1000_150_100_2018_IdNohoeEff_pt'])
    hoeEff.Add(plot['HSS_1000_150_1000_2018_IdNohoeEff_pt'])
    hoeEff.Add(plot['HSS_1000_150_10000_2018_IdNohoeEff_pt'])

    sigmaEff = plot['HSS_1000_150_1_2018_IdNosigmaEff_pt']
    sigmaEff.Add(plot['HSS_1000_150_10_2018_IdNosigmaEff_pt'])
    sigmaEff.Add(plot['HSS_1000_150_100_2018_IdNosigmaEff_pt'])
    sigmaEff.Add(plot['HSS_1000_150_1000_2018_IdNosigmaEff_pt'])
    sigmaEff.Add(plot['HSS_1000_150_10000_2018_IdNosigmaEff_pt'])

    IdEff = plot['HSS_1000_150_1_2018_IdEff_pt']
    IdEff.Add(plot['HSS_1000_150_10_2018_IdEff_pt'])
    IdEff.Add(plot['HSS_1000_150_100_2018_IdEff_pt'])
    IdEff.Add(plot['HSS_1000_150_1000_2018_IdEff_pt'])
    IdEff.Add(plot['HSS_1000_150_10000_2018_IdEff_pt'])

    canvas = Canvas.Canvas("HSS_1000_150_2018_NoIdEff_pt", 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1)
    canvas.addHisto(hitsEff,'AP', 'All requirements - hit number', 'p', gcolors['blue'], True, 0, marker = 24)
    canvas.addHisto(hpEff,'AP,SAME', 'All requirements - high purity', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(hoeEff,'AP,SAME', 'All requirements - SC H/E', 'p', gcolors['yellow'], True, 2, marker = 24)
    canvas.addHisto(sigmaEff,'AP,SAME', 'All requirements - SC #sigma_{i#etai#eta}', 'p', gcolors['green'], True, 3, marker = 24)
    canvas.addHisto(IdEff,'AP,SAME', 'All requirements', 'p', r.kBlack, True, 4, marker = 24)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.6, ymax = 1.1, outputDir = WORKPATH + 'plots/', isPrivate = True)


    # The Id - cut efficiency (eta):
    hitsEff = plot['HSS_1000_150_1_2018_IdNohitEff_eta']
    hitsEff.Add(plot['HSS_1000_150_10_2018_IdNohitEff_eta'])
    hitsEff.Add(plot['HSS_1000_150_100_2018_IdNohitEff_eta'])
    hitsEff.Add(plot['HSS_1000_150_1000_2018_IdNohitEff_eta'])
    hitsEff.Add(plot['HSS_1000_150_10000_2018_IdNohitEff_eta'])

    hpEff = plot['HSS_1000_150_1_2018_IdNohpEff_eta']
    hpEff.Add(plot['HSS_1000_150_10_2018_IdNohpEff_eta'])
    hpEff.Add(plot['HSS_1000_150_100_2018_IdNohpEff_eta'])
    hpEff.Add(plot['HSS_1000_150_1000_2018_IdNohpEff_eta'])
    hpEff.Add(plot['HSS_1000_150_10000_2018_IdNohpEff_eta'])

    hoeEff = plot['HSS_1000_150_1_2018_IdNohoeEff_eta']
    hoeEff.Add(plot['HSS_1000_150_10_2018_IdNohoeEff_eta'])
    hoeEff.Add(plot['HSS_1000_150_100_2018_IdNohoeEff_eta'])
    hoeEff.Add(plot['HSS_1000_150_1000_2018_IdNohoeEff_eta'])
    hoeEff.Add(plot['HSS_1000_150_10000_2018_IdNohoeEff_eta'])

    sigmaEff = plot['HSS_1000_150_1_2018_IdNosigmaEff_eta']
    sigmaEff.Add(plot['HSS_1000_150_10_2018_IdNosigmaEff_eta'])
    sigmaEff.Add(plot['HSS_1000_150_100_2018_IdNosigmaEff_eta'])
    sigmaEff.Add(plot['HSS_1000_150_1000_2018_IdNosigmaEff_eta'])
    sigmaEff.Add(plot['HSS_1000_150_10000_2018_IdNosigmaEff_eta'])

    IdEff = plot['HSS_1000_150_1_2018_IdEff_eta']
    IdEff.Add(plot['HSS_1000_150_10_2018_IdEff_eta'])
    IdEff.Add(plot['HSS_1000_150_100_2018_IdEff_eta'])
    IdEff.Add(plot['HSS_1000_150_1000_2018_IdEff_eta'])
    IdEff.Add(plot['HSS_1000_150_10000_2018_IdEff_eta'])

    canvas = Canvas.Canvas("HSS_1000_150_2018_NoIdEff_eta", 'png,pdf', 0.2, 0.2, 0.4, 0.35, 1)
    canvas.addHisto(hitsEff,'AP', 'All requirements - hit number', 'p', gcolors['blue'], True, 0, marker = 24)
    canvas.addHisto(hpEff,'AP,SAME', 'All requirements - high purity', 'p', gcolors['red'], True, 1, marker = 24)
    canvas.addHisto(hoeEff,'AP,SAME', 'All requirements - SC H/E', 'p', gcolors['yellow'], True, 2, marker = 24)
    canvas.addHisto(sigmaEff,'AP,SAME', 'All requirements - SC #sigma_{i#etai#eta}', 'p', gcolors['green'], True, 3, marker = 24)
    canvas.addHisto(IdEff,'AP,SAME', 'All requirements', 'p', r.kBlack, True, 4, marker = 24)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.45, 0.87, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.45, 0.83, 'M_{H} = 1000 GeV, M_{S} = 150 GeV', size = 0.03)
    canvas.addLatex(0.45, 0.79, '(all lifetimes combined)', size = 0.03)
    canvas.save(1, 0, 0, '', '', ymin = 0.4, ymax = 1.1, outputDir = WORKPATH + 'plots/', isPrivate = True)

    ## Resolution and electron plots
    plot['HSS_1000_150_1_2018_electron_dR'].Scale(1./plot['HSS_1000_150_1_2018_electron_dR'].Integral()) 
    plot['HSS_1000_150_10_2018_electron_dR'].Scale(1./plot['HSS_1000_150_10_2018_electron_dR'].Integral()) 
    plot['HSS_1000_150_100_2018_electron_dR'].Scale(1./plot['HSS_1000_150_100_2018_electron_dR'].Integral()) 
    plot['HSS_1000_150_1000_2018_electron_dR'].Scale(1./plot['HSS_1000_150_1000_2018_electron_dR'].Integral()) 
    plot['HSS_1000_150_1_2018_electron_dR'].SetLineWidth(2)
    plot['HSS_1000_150_10_2018_electron_dR'].SetLineWidth(2)
    plot['HSS_1000_150_100_2018_electron_dR'].SetLineWidth(2)
    plot['HSS_1000_150_1000_2018_electron_dR'].SetLineWidth(2)
    plot['HSS_1000_150_1_2018_electron_dR'].SetMaximum(0.5)
    canvas = Canvas.Canvas("HSS_1000_150_2018_electron_dR", 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1)
    canvas.addHisto(plot['HSS_1000_150_1_2018_electron_dR'],'HIST', 'c#tau = 1 mm', 'l', acolors['1'], True, 0)
    canvas.addHisto(plot['HSS_1000_150_10_2018_electron_dR'],'HIST,SAME', 'c#tau = 10 mm', 'l', acolors['2'], True, 1)
    canvas.addHisto(plot['HSS_1000_150_100_2018_electron_dR'],'HIST,SAME', 'c#tau = 100 mm', 'l', acolors['3'], True, 2)
    canvas.addHisto(plot['HSS_1000_150_1000_2018_electron_dR'],'HIST,SAME', 'c#tau = 1000 mm', 'l', acolors['4'], True, 3)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)

    plot['HSS_1000_150_1_2018_electron_ptres'].Scale(1./plot['HSS_1000_150_1_2018_electron_ptres'].Integral()) 
    plot['HSS_1000_150_10_2018_electron_ptres'].Scale(1./plot['HSS_1000_150_10_2018_electron_ptres'].Integral()) 
    plot['HSS_1000_150_100_2018_electron_ptres'].Scale(1./plot['HSS_1000_150_100_2018_electron_ptres'].Integral()) 
    plot['HSS_1000_150_1000_2018_electron_ptres'].Scale(1./plot['HSS_1000_150_1000_2018_electron_ptres'].Integral()) 
    canvas = Canvas.Canvas("HSS_1000_150_2018_electron_ptres", 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1)
    canvas.addHisto(plot['HSS_1000_150_1_2018_electron_ptres'],'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], True, 0)
    canvas.addHisto(plot['HSS_1000_150_10_2018_electron_ptres'],'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], True, 1)
    canvas.addHisto(plot['HSS_1000_150_100_2018_electron_ptres'],'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], True, 2)
    canvas.addHisto(plot['HSS_1000_150_1000_2018_electron_ptres'],'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], True, 3)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)

    plot['HSS_1000_150_1_2018_electron_dxyres'].Scale(1./plot['HSS_1000_150_1_2018_electron_dxyres'].Integral()) 
    plot['HSS_1000_150_10_2018_electron_dxyres'].Scale(1./plot['HSS_1000_150_10_2018_electron_dxyres'].Integral()) 
    plot['HSS_1000_150_100_2018_electron_dxyres'].Scale(1./plot['HSS_1000_150_100_2018_electron_dxyres'].Integral()) 
    plot['HSS_1000_150_1000_2018_electron_dxyres'].Scale(1./plot['HSS_1000_150_1000_2018_electron_dxyres'].Integral()) 
    plot['HSS_1000_150_1_2018_electron_dxyres'].SetMaximum(0.5)
    canvas = Canvas.Canvas("HSS_1000_150_2018_electron_dxyres", 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1)
    canvas.addHisto(plot['HSS_1000_150_1_2018_electron_dxyres'],'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], True, 0)
    canvas.addHisto(plot['HSS_1000_150_10_2018_electron_dxyres'],'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], True, 1)
    canvas.addHisto(plot['HSS_1000_150_100_2018_electron_dxyres'],'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], True, 2)
    canvas.addHisto(plot['HSS_1000_150_1000_2018_electron_dxyres'],'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], True, 3)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)


    plot['HSS_1000_150_1_2018_electron_etres'].Scale(1./plot['HSS_1000_150_1_2018_electron_etres'].Integral()) 
    plot['HSS_1000_150_10_2018_electron_etres'].Scale(1./plot['HSS_1000_150_10_2018_electron_etres'].Integral()) 
    plot['HSS_1000_150_100_2018_electron_etres'].Scale(1./plot['HSS_1000_150_100_2018_electron_etres'].Integral()) 
    plot['HSS_1000_150_1000_2018_electron_etres'].Scale(1./plot['HSS_1000_150_1000_2018_electron_etres'].Integral()) 
    canvas = Canvas.Canvas("HSS_1000_150_2018_electron_etres", 'png,pdf', 0.6, 0.6, 0.89, 0.75, 1)
    canvas.addHisto(plot['HSS_1000_150_1_2018_electron_etres'],'HIST', 'c#tau = 1 mm', 'l', dcolors['1mm'], True, 0)
    canvas.addHisto(plot['HSS_1000_150_10_2018_electron_etres'],'HIST,SAME', 'c#tau = 10 mm', 'l', dcolors['10mm'], True, 1)
    canvas.addHisto(plot['HSS_1000_150_100_2018_electron_etres'],'HIST,SAME', 'c#tau = 100 mm', 'l', dcolors['100mm'], True, 2)
    canvas.addHisto(plot['HSS_1000_150_1000_2018_electron_etres'],'HIST,SAME', 'c#tau = 1000 mm', 'l', dcolors['1000mm'], True, 3)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.addLatex(0.6, 0.85, 'H#rightarrowSS#rightarrow 2e + X', size = 0.03)
    canvas.addLatex(0.6, 0.81, 'M_{H} = 1000 GeV', size = 0.03)
    canvas.addLatex(0.6, 0.77, 'M_{S} = 150 GeV', size = 0.03)
    canvas.save(1, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)

    r.gStyle.SetPadRightMargin(0.14)
    genpt_scet = combinePlots([plot['HSS_1000_150_1_2018_electron_pt_et'], plot['HSS_1000_150_10_2018_electron_pt_et'], plot['HSS_1000_150_100_2018_electron_pt_et'], plot['HSS_1000_150_1000_2018_electron_pt_et'], plot['HSS_1000_150_10000_2018_electron_pt_et']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_electron_pt_et", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 610, hh = 600)
    canvas.addHisto(genpt_scet,'COLZ', '', '', '', True, 0)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True, is2d = True)
    r.gStyle.SetPadRightMargin(0.1)

    ## Others
    histo = plot['HSS_1000_150_1_2018_clustersAssocToTrack']
    histo.Add(plot['HSS_1000_150_10_2018_clustersAssocToTrack'])
    histo.Add(plot['HSS_1000_150_100_2018_clustersAssocToTrack'])
    histo.Add(plot['HSS_1000_150_1000_2018_clustersAssocToTrack'])
    histo.Add(plot['HSS_1000_150_10000_2018_clustersAssocToTrack'])
    canvas = Canvas.Canvas("HSS_1000_150_2018_clustersAssocToTrack", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(histo,'P', 'Track counting', 'p', r.kBlack, True, 0, marker = 20)
    canvas.save(0, 0, 1, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True)


    ## 2D plots
    ## Supercluster resolution
    genpt_scet = combinePlots([plot['HSS_1000_150_1_2018_genpt_scet'], plot['HSS_1000_150_10_2018_genpt_scet'], plot['HSS_1000_150_100_2018_genpt_scet'], plot['HSS_1000_150_1000_2018_genpt_scet'], plot['HSS_1000_150_10000_2018_genpt_scet']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_genpt_scet", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 610, hh = 600)
    canvas.addHisto(genpt_scet,'COLZ', '', '', '', True, 0)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    genpt_scet_passSigma = combinePlots([plot['HSS_1000_150_1_2018_genpt_scet_passSigma'], plot['HSS_1000_150_10_2018_genpt_scet_passSigma'], plot['HSS_1000_150_100_2018_genpt_scet_passSigma'], plot['HSS_1000_150_1000_2018_genpt_scet_passSigma'], plot['HSS_1000_150_10000_2018_genpt_scet_passSigma']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_genpt_scet_passSigma", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 610, hh = 600)
    canvas.addHisto(genpt_scet_passSigma,'COLZ', '', '', '', True, 0)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    genpt_scet_passHoE = combinePlots([plot['HSS_1000_150_1_2018_genpt_scet_passHoE'], plot['HSS_1000_150_10_2018_genpt_scet_passHoE'], plot['HSS_1000_150_100_2018_genpt_scet_passHoE'], plot['HSS_1000_150_1000_2018_genpt_scet_passHoE'], plot['HSS_1000_150_10000_2018_genpt_scet_passHoE']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_genpt_scet_passHoE", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1, ww = 610, hh = 600)
    canvas.addHisto(genpt_scet_passHoE,'COLZ', '', '', '', True, 0)
    canvas.addLatex(0.9, 0.93, '2018 UL', size = 0.035, align = 31)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    """
    sc_dEta_vs_Lxy = combinePlots([plot['HSS_1000_150_1_2018_sc_dEta_vs_Lxy'], plot['HSS_1000_150_10_2018_sc_dEta_vs_Lxy'], plot['HSS_1000_150_100_2018_sc_dEta_vs_Lxy'], plot['HSS_1000_150_1000_2018_sc_dEta_vs_Lxy'], plot['HSS_1000_150_10000_2018_sc_dEta_vs_Lxy']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_sc_dEta_vs_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(sc_dEta_vs_Lxy,'COLZ', '', '', '', True, 0)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    track_dEta_vs_Lxy = combinePlots([plot['HSS_1000_150_1_2018_track_dEta_vs_Lxy'], plot['HSS_1000_150_10_2018_track_dEta_vs_Lxy'], plot['HSS_1000_150_100_2018_track_dEta_vs_Lxy'], plot['HSS_1000_150_1000_2018_track_dEta_vs_Lxy'], plot['HSS_1000_150_10000_2018_track_dEta_vs_Lxy']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_track_dEta_vs_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(track_dEta_vs_Lxy,'COLZ', '', '', '', True, 0)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    sc_dPhi_vs_Lxy = combinePlots([plot['HSS_1000_150_1_2018_sc_dPhi_vs_Lxy'], plot['HSS_1000_150_10_2018_sc_dPhi_vs_Lxy'], plot['HSS_1000_150_100_2018_sc_dPhi_vs_Lxy'], plot['HSS_1000_150_1000_2018_sc_dPhi_vs_Lxy'], plot['HSS_1000_150_10000_2018_sc_dPhi_vs_Lxy']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_sc_dPhi_vs_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(sc_dPhi_vs_Lxy,'COLZ', '', '', '', True, 0)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    track_dPhi_vs_Lxy = combinePlots([plot['HSS_1000_150_1_2018_track_dPhi_vs_Lxy'], plot['HSS_1000_150_10_2018_track_dPhi_vs_Lxy'], plot['HSS_1000_150_100_2018_track_dPhi_vs_Lxy'], plot['HSS_1000_150_1000_2018_track_dPhi_vs_Lxy'], plot['HSS_1000_150_10000_2018_track_dPhi_vs_Lxy']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_track_dPhi_vs_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(track_dPhi_vs_Lxy,'COLZ', '', '', '', True, 0)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    sc_dR_vs_Lxy = combinePlots([plot['HSS_1000_150_1_2018_sc_dR_vs_Lxy'], plot['HSS_1000_150_10_2018_sc_dR_vs_Lxy'], plot['HSS_1000_150_100_2018_sc_dR_vs_Lxy'], plot['HSS_1000_150_1000_2018_sc_dR_vs_Lxy'], plot['HSS_1000_150_10000_2018_sc_dR_vs_Lxy']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_sc_dR_vs_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(sc_dR_vs_Lxy,'COLZ', '', '', '', True, 0)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)

    track_dR_vs_Lxy = combinePlots([plot['HSS_1000_150_1_2018_track_dR_vs_Lxy'], plot['HSS_1000_150_10_2018_track_dR_vs_Lxy'], plot['HSS_1000_150_100_2018_track_dR_vs_Lxy'], plot['HSS_1000_150_1000_2018_track_dR_vs_Lxy'], plot['HSS_1000_150_10000_2018_track_dR_vs_Lxy']])
    canvas = Canvas.Canvas("HSS_1000_150_2018_track_dR_vs_Lxy", 'png,pdf', 0.4, 0.8, 0.8, 0.9, 1)
    canvas.addHisto(track_dR_vs_Lxy,'COLZ', '', '', '', True, 0)
    canvas.save(0, 0, 0, '', '', outputDir = WORKPATH + 'plots/', isPrivate = True, zlog = True)
    """


