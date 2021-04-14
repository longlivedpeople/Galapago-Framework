
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
    parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    parser.add_option('-f', '--filename', action='store', type=str, dest='filename', default='', help='Path to file')
    parser.add_option('-n', '--nmax', action='store', type=int, dest='nmax', default=0, help='Path to file')
    (opts, args) = parser.parse_args()


    ##################################
    ####   Variable declaration   ####
    ##################################
    Lxy_bin = np.linspace(0.0, 100, 101)

    #### -----------------
    #### ---- Histograms
    #### -----------------
    hist_Lxy = r.TH1F("hist_Lxy", ";Decay length (CMS frame) L_{xy} (cm); Generated electron yield", 100, 0.0, 100.0)
    hist_t = r.TH1F("hist_t", ";Decay time (CMS frame) t (ns); Generated electron yield", 100, 0.0, 100.0)

    
    #########################
    ####   Load sample   ####
    #########################

    _dirName = opts.filename
    _tree = r.TChain('Events')
    for _file in os.listdir(_dirName):
        if '.root' not in _file: continue
        _tree.Add(_dirName + _file)

    print("TTree with " + str(_tree.GetEntries()) + " entries")

    # Get the X mass
    lastdir = _dirName.split('/')[-2]
    splits = lastdir.split('_')
    mass = float(splits[2])


    ###################################
    ####   Loop over tree events   ####
    ###################################

    for i in range(0, _tree.GetEntries()):

        _tree.GetEntry(i)
        if opts.nmax and i > opts.nmax: break

        X_vx = -99
        X_vy = -99
        vxs = []
        vys = []
        pts_1 = []
        phis_1 = []
        pts_2 = []
        phis_2 = []

        for n in range(0, _tree.nHardProcessParticle):

            pdgId = _tree.HardProcessParticle_pdgId[n]
            vx = _tree.HardProcessParticle_vx[n]
            vy = _tree.HardProcessParticle_vy[n]

            if pdgId == 54 and X_vx < -98:
                X_vx = vx
                X_vy = vy

            if pdgId == 11 and vx not in vxs:

                vxs.append(vx)                            
                vys.append(vy)                            
                pts_1.append(_tree.HardProcessParticle_pt[n])
                phis_1.append(_tree.HardProcessParticle_phi[n])

                for j in range(0, _tree.nHardProcessParticle):
                    if _tree.HardProcessParticle_pdgId[j] != -11: continue
                    if _tree.HardProcessParticle_vx[j] == vx:
                        pts_2.append(_tree.HardProcessParticle_pt[j])
                        phis_2.append(_tree.HardProcessParticle_phi[j])


        if len(vxs) < 1: continue

        for n in range(0, len(vxs)):

            Lxy = math.sqrt((vxs[n] - X_vx)**2 + (vys[n] - X_vy)**2)
            l1 = r.TVector3()
            l2 = r.TVector3()
            l1.SetPtEtaPhi(pts_1[n], 0.0, phis_1[n])
            l2.SetPtEtaPhi(pts_2[n], 0.0, phis_2[n])
            l12 = l1 + l2
            pt_X = l12.Pt()

            time = Lxy * math.sqrt(mass**2 + pt_X**2)/(pt_X*3e10)

            time = time*1e9

            hist_Lxy.Fill(Lxy)
            hist_t.Fill(time)


    if not os.path.exists(WORKPATH + 'GenResults/'): os.makedirs(WORKPATH + 'GenResults/')
    outputFile = TFile(WORKPATH + 'GenResults/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    hist_Lxy.Write()
    hist_t.Write()
    outputFile.Close()



