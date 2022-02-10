
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

def makeJointSignalPlot(values, key, ylabel, plotname, modellabel = ''):

    # "values" is a dict with the values we want to plot

    models = []
    for name in values.keys():
        print(name)
        info = name.split('_')
        mH = info[1]
        mS = info[2]
        ctau = info[3]
        eff = values[name][key]
        models.append([mH, mS, ctau, eff])

    models = sorted(models, key = lambda x: (float(x[0]), float(x[1]), float(x[2])))

    nsamples = len(models)
    histo = r.TH1F('histo', '', nsamples, 0, nsamples)
    histo.GetXaxis().SetLabelSize(0)
    histo.GetYaxis().SetTitle(ylabel)
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

    c1.SaveAs(plotname + '.png')
    c1.SaveAs(plotname + '.pdf')

    



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

    year = '2017'
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
            rates[s.name]['mm_total'] = 0.0
            rates[s.name]['mm_passed'] = 0.0
            rates[s.name]['mm_total_disp'] = 0.0
            rates[s.name]['mm_passed_disp'] = 0.0
            rates[s.name]['ee_total'] = 0.0
            rates[s.name]['ee_passed'] = 0.0
            rates[s.name]['ee_total_disp'] = 0.0
            rates[s.name]['ee_passed_disp'] = 0.0

            for t in s.ttrees:
                print(s.name, 'New tree with:', t.GetEntries())
                for e,ev in enumerate(t):

                    #### Dimuon channel:

                    passed_trigger = False
                    if year == '2018': 
                        if (ev.HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed or ev.HLT_DoubleL2Mu23NoVtx_2Cha):
                            passed_trigger = True
                    if year == '2016':
                        if (ev.HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10):
                            passed_trigger = True

                    if passed_trigger and ev.nDMDM > 0: 

                        passed_counted = False
                        total_counted = False

                        for j in range(0, ev.nDMDM):

                            imm = j # index to handle DMDM pair
                            if year == '2016':
                                if not eval(cm.MM_BS2016): continue
                            if year == '2018':
                                if not eval(cm.MM_BS2018): continue

                            if not eval(cm.MM_iso2l): continue
                            if not eval(cm.MM_OS): continue

                            total_counted = True
                            # Evaluate the cut:
                            if abs(ev.DMDM_dPhi[imm]) < 3.14/2.:
                                passed_counted = True

                        if total_counted:
                            rates[s.name]['mm_total'] += 1.
                        if passed_counted:
                            rates[s.name]['mm_passed'] += 1.

                        if total_counted and ev.DMDM_trackIxy_PV[imm] > 1:
                            rates[s.name]['mm_total_disp'] += 1.
                        if passed_counted and ev.DMDM_trackIxy_PV[imm] > 1:
                            rates[s.name]['mm_passed_disp'] += 1.

                    #### Dielectron channel:

                    passed_trigger = False
                    if year == '2016': 
                        if eval(cm.epath2016): 
                            passed_trigger = True
                    if year == '2017': 
                        if eval(cm.epath2017): 
                            passed_trigger = True
                    if year == '2018': 
                        if eval(cm.epath2018): 
                            passed_trigger = True

                    if passed_trigger and ev.nEE > 0: 

                        passed_counted = False
                        total_counted = False

                        for j in range(0, ev.nEE):

                            iee = j # index to handle DMDM pair
                            if year == '2016':
                                if not eval(cm.EE_BS2016): continue
                            if year == '2017':
                                if not eval(cm.EE_BS2017): continue
                            if year == '2018':
                                if not eval(cm.EE_BS2018): continue

                            if not eval(cm.EE_iso2l): continue
                            if not eval(cm.EE_OS): continue

                            total_counted = True
                            # Evaluate the cut:
                            if abs(ev.EE_dPhi[iee]) < 3.14/2.:
                                passed_counted = True

                        if total_counted:
                            rates[s.name]['ee_total'] += 1.
                        if passed_counted:
                            rates[s.name]['ee_passed'] += 1.
                        if total_counted and ev.EE_trackIxy_PV[iee] > 1:
                            rates[s.name]['ee_total_disp'] += 1.
                        if passed_counted and ev.EE_trackIxy_PV[iee] > 1:
                            rates[s.name]['ee_passed_disp'] += 1.

            if year != '2017':
                rates[s.name]['mm_efficiency'] = rates[s.name]['mm_passed']/rates[s.name]['mm_total']
                rates[s.name]['mm_efficiency_disp'] = rates[s.name]['mm_passed_disp']/rates[s.name]['mm_total_disp']
            rates[s.name]['ee_efficiency'] = rates[s.name]['ee_passed']/rates[s.name]['ee_total']
            rates[s.name]['ee_efficiency_disp'] = rates[s.name]['ee_passed_disp']/rates[s.name]['ee_total_disp']
            #print(rates[s.name]['mm_efficiency'], rates[s.name]['ee_efficiency'])

    ## plot results
    if year != '2017':
        makeJointSignalPlot(values = rates, key = 'mm_efficiency', ylabel = '|#Delta#Phi| cut efficiency', plotname = 'Dimuon_dPhi_efficiency_'+year, modellabel = 'H#rightarrow2S#rightarrow4l ('+year+')')
        makeJointSignalPlot(values = rates, key = 'mm_efficiency_disp', ylabel = '|#Delta#Phi| cut efficiency', plotname = 'Dimuon_dPhi_efficiency_disp_'+year, modellabel = 'H#rightarrow2S#rightarrow4l ('+year+')')
    makeJointSignalPlot(values = rates, key = 'ee_efficiency', ylabel = '|#Delta#Phi| cut efficiency', plotname = 'Dielectron_dPhi_efficiency_'+year, modellabel = 'H#rightarrow2S#rightarrow4l ('+year+')')
    makeJointSignalPlot(values = rates, key = 'ee_efficiency_disp', ylabel = '|#Delta#Phi| cut efficiency', plotname = 'Dielectron_dPhi_efficiency_disp_'+year, modellabel = 'H#rightarrow2S#rightarrow4l ('+year+')')

    



