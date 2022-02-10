
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


def createTickLabel(align = 21):

    latex = r.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(r.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(align);
    latex.SetTextSize(0.03);

    return latex



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

    #### -----------------
    #### ---- Histograms
    #### -----------------
    
    #########################
    ####   Load sample   ####
    #########################

    year = '2016'
    _cosAlphaMin = -0.8

    Signals2016 = []
    Signals2016.append('HSS_125_50_1_'+year)
    Signals2016.append('HSS_125_50_10_'+year)
    Signals2016.append('HSS_125_50_100_'+year)
    Signals2016.append('HSS_400_50_1_'+year)
    Signals2016.append('HSS_400_50_10_'+year)
    Signals2016.append('HSS_400_50_100_'+year)
    Signals2016.append('HSS_400_150_1_'+year)
    Signals2016.append('HSS_400_150_10_'+year)
    Signals2016.append('HSS_400_150_100_'+year)
    Signals2016.append('HSS_600_150_1_'+year)
    Signals2016.append('HSS_600_150_10_'+year)
    Signals2016.append('HSS_600_150_100_'+year)
    Signals2016.append('HSS_1000_150_1_'+year)
    Signals2016.append('HSS_1000_150_10_'+year)
    Signals2016.append('HSS_1000_150_100_'+year)
    Signals2016.append('HSS_1000_350_1_'+year)
    Signals2016.append('HSS_1000_350_10_'+year)
    Signals2016.append('HSS_1000_350_100_'+year)

    treeSI = Sample.Tree( fileName = helper.selectSamples(GALAPAGOPATH + 'signals_'+year+'.dat', Signals2016, 'MC'), name = year, isdata = 0 )

    ###################################
    ####   Loop over tree events   ####
    ###################################

    cm = CutManager.CutManager()
    rates = {}
    for b in treeSI.blocks:
        for s in b.samples:
            rates[s.name] = {}
            rates[s.name]['total'] = 0.0
            rates[s.name]['passed'] = 0.0

            for t in s.ttrees:
                print(s.name, 'New tree with:', t.GetEntries())
                for e,ev in enumerate(t):

                    #if opts.nmax and i > opts.nmax: break
                    #if not (ev.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed or ev.HLT_DoubleL2Mu23NoVtx_2Cha): continue
                    if not (ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10): continue
                    if ev.nDMDM < 1: continue

                    passed_counted = False
                    total_counted = False

                    for j in range(0, ev.nDMDM):

                        imm = j # index to handle DMDM pair
                        if year == '2016':
                            if not ev.DGM_pt[ev.DMDM_idxA[imm]] > 30 or not ev.DGM_pt[ev.DMDM_idxB[imm]] > 30: continue
                            #if not ev.DMDM_cosAlpha[imm] > -0.80: continue
                        if year == '2018':
                            if not ev.DGM_pt[ev.DMDM_idxA[imm]] > 30 or not ev.DGM_pt[ev.DMDM_idxB[imm]] > 30: continue
                            #if not ev.DMDM_cosAlpha[imm] > -0.80: continue

                        if not ev.DMDM_mass[imm] > 15: continue
                        if not ev.DMDM_normalizedChi2[imm] < 10: continue
                        if not ev.DMDM_dR[imm] > 0.2: continue
                        if not abs(ev.DGM_eta[ev.DMDM_idxA[imm]]) < 2.0 or not abs(ev.DGM_eta[ev.DMDM_idxB[imm]]) < 2.0: continue
                        if not ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] < 0: continue
                        if not eval(cm.MM_iso2l): continue
                        if not abs(ev.DMDM_dPhi[imm]) < 3.14/2.: continue
            
                        #ID:
                        if not ev.DGM_ptError[ev.DMDM_idxB[imm]]/ev.DGM_pt[ev.DMDM_idxB[imm]] < 0.3: continue
                        if not ev.DGM_ptError[ev.DMDM_idxA[imm]]/ev.DGM_pt[ev.DMDM_idxA[imm]] < 0.3: continue
                        if not ev.DGM_normChi2[ev.DMDM_idxA[imm]] < 10.: continue
                        if not ev.DGM_normChi2[ev.DMDM_idxB[imm]] < 10.: continue
                        if not ev.DGM_muonHits[ev.DMDM_idxA[imm]] > 11: continue  
                        if not ev.DGM_muonHits[ev.DMDM_idxB[imm]] > 11: continue  
                        if not ev.DGM_outerTrackerHits[ev.DMDM_idxA[imm]] > 5: continue  
                        if not ev.DGM_outerTrackerHits[ev.DMDM_idxB[imm]] > 5: continue  

                        total_counted = True
                        # Evaluate the cut:
                        if ev.DMDM_cosAlpha[imm] > _cosAlphaMin:
                            passed_counted = True

                    if total_counted:
                        rates[s.name]['total'] += 1.
                    if passed_counted:
                        rates[s.name]['passed'] += 1.

            rates[s.name]['efficiency'] = rates[s.name]['passed']/rates[s.name]['total']
            print(rates[s.name]['efficiency'])

    ## plot results

    models = []
    for name in rates.keys():
        print(name)
        info = name.split('_')
        mH = info[1]
        mS = info[2]
        ctau = info[3]
        eff = rates[name]['efficiency']
        models.append([mH, mS, ctau, eff])

    models = sorted(models, key = lambda x: (float(x[0]), float(x[1]), float(x[2])))

    nsamples = len(models)
    histo = r.TH1F('histo', '', nsamples, 0, nsamples)
    histo.GetXaxis().SetLabelSize(0)
    histo.GetYaxis().SetTitle("cos(#alpha) cut efficiency")
    histo.GetXaxis().SetNdivisions(nsamples)
    histo.SetMarkerColor(r.kBlack)
    histo.SetLineColor(r.kBlack)
    histo.SetFillColor(r.kAzure-5)
    histo.SetMarkerStyle(20)
    histo.SetMaximum(1.0)
    histo.SetMinimum(0.7)

    for _m,m in enumerate(models):
        print(m[3])
        histo.SetBinContent(_m + 1, m[3])
        histo.SetBinError(_m + 1, 0.00001)

    ## Plot
    lm = r.gStyle.GetPadLeftMargin()
    rm = r.gStyle.GetPadRightMargin()
    tm = r.gStyle.GetPadTopMargin()
    bm = 0.12
    r.gStyle.SetPadBottomMargin(bm)

    c1 = r.TCanvas("c1", "", 600, 500)
    histo.Draw('HIST')

    ## Ticks
    wd = (1.0 - lm - rm)/len(models)
    counter = 0
    massmargin = lm
    lines = []
    for _m,m in enumerate(models):
        counter += 1
        ctau = str(m[2])
        order = ctau.count('0')
        ctau_text = ''
        if order == 0: ctau_text = '1'
        elif order == 1: ctau_text = '10'
        else: ctau_text = '10^{'+str(order)+'}'
        ctau_label = createTickLabel()
        ctau_label.DrawLatex(lm + (_m+0.5)*wd ,0.08, ctau_text)
        if _m +1 == len(models):
            mass_label = createTickLabel()
            mass_label.DrawLatex(massmargin + counter*wd/2, 0.03, '({0}, {1})'.format(str(m[0]), str(m[1])))
        elif m[0] != models[_m +1][0] or m[1] != models[_m +1][1]:
            mass_label = createTickLabel()
            mass_label.DrawLatex(massmargin + counter*wd/2, 0.03, '({0}, {1})'.format(str(m[0]), str(m[1])))
            massmargin = massmargin + counter*wd
            lines.append(r.TLine(massmargin, bm, massmargin, 1. - tm))
            lines[-1].SetNDC()
            lines[-1].SetLineStyle(2)
            lines[-1].Draw()
            counter = 0

    ctau_title = createTickLabel(11)
    mass_title = createTickLabel(11)
    ctau_title.DrawLatex(0.02, 0.08, 'c#tau')
    mass_title.DrawLatex(0.02, 0.03, '(M_{H}, M_{S})')
    ctau_units = createTickLabel(11)
    mass_units = createTickLabel(11)
    ctau_units.DrawLatex(0.91, 0.08, '[mm]')
    mass_units.DrawLatex(0.91, 0.03, '[GeV]')

    ## CMS banners
    latex = r.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(r.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.06);
    latex.DrawLatex(lm, 0.93, "#bf{CMS}")

    latexb = r.TLatex()
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextColor(r.kBlack);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(31);
    latexb.SetTextSize(0.034);
    latexb.DrawLatex(0.35, 0.93, "#it{Simulation}")

    ## Model banner
    latext = r.TLatex()
    latext.SetNDC();
    latext.SetTextAngle(0);
    latext.SetTextColor(r.kBlack);
    latext.SetTextFont(42);
    latext.SetTextAlign(31);
    latext.SetTextSize(0.033);
    latext.DrawLatex(0.9, 0.93, 'H#rightarrow2S#rightarrow4l ('+year+')')

    #if not os.path.exists(output): os.makedirs(output)

    c1.SaveAs('cosAlpha_efficiency_'+year+ '.png')
    c1.SaveAs('cosAlpha_efficiency_'+year+ '.pdf')





        





    """
    if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    outputFile = TFile(WORKPATH + 'Results/th1f'+opts.tag+'.root', 'RECREATE')


    #### Write everything to use later:
    hist_counts.Write()
    hist_cosAlpha.Write()
    hist_deltaPhi.Write()
    hist_alpha.Write()
    hist_normChi2.Write()
    hist_dPhi.Write()
    outputFile.Close()
    """


