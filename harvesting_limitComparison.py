import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time
import shutil

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager
#from include.Utils import *

from CombineHarvester.CombineTools.plotting import *
r.PyConfig.IgnoreCommandLineOptions = True


CMSstyle = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : R.kBlue, 'LineStyle' : 4},
               'exp1' : { 'FillColor' : R.kGreen+1},
               'exp2' : { 'FillColor' : R.kOrange}
               }

red_exp0 = CreateTransparentColor(R.kRed+1, 0.9)
red_exp1 = CreateTransparentColor(R.kRed-4, 0.3)
red_exp2 = CreateTransparentColor(R.kRed-9, 0.3)

style_red_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : red_exp0},
               'exp1' : { 'FillColor' : red_exp1},
               'exp2' : { 'FillColor' : red_exp2}
               }

blue_exp0 = CreateTransparentColor(R.kBlue+1, 0.9)
blue_exp1 = CreateTransparentColor(R.kBlue-4, 0.3)
blue_exp2 = CreateTransparentColor(R.kBlue-9, 0.3)

style_blue_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : blue_exp0},
               'exp1' : { 'FillColor' : blue_exp1},
               'exp2' : { 'FillColor' : blue_exp2}
               }

#####################
#####
###
###   Function to plot limit comparison with our analysis and other analyses
###
#####
#####################

def makeLimitComparison(name, lumi, limitSet, masslabel = [], ymin = 1e-5, ymax = 1e6, xmin = 1e-3, xmax = 1e6):

    ModTDRStyle()
    r.gStyle.SetLegendFont(42)
    r.gStyle.SetLegendTextSize(0.033)
    r.gStyle.SetNdivisions(510, 'Y')
    canv = r.TCanvas('limit', 'limit')
    canv.SetFillStyle(4000);
    canv.SetFrameFillColor(4000);
    canv.SetFrameFillStyle(4000);
    pads = OnePad()

    ### Get the reference graph
    referenceLimit = StandardLimitsFromJSONFile(limitSet['Reference'][0])

    ### Get the graphs (pending)
    graphs = []
    for other in limitSet['Others']:
        _file = r.TFile(other[0])
        _graph = _file.Get(other[1])
        _graph.SetLineColor(other[3])
        _graph.SetLineStyle(other[4])
        _graph.SetLineWidth(2)
        graphs.append([copy.deepcopy(_graph), other[2]])


    ### Create axis to draw the limits
    #axis = CreateAxisHist(referenceLimit['exp0'])
    axis = r.TH1F('axis', '', 1, xmin, xmax)
    axis.GetXaxis().SetTitle('c#tau [cm]')
    axis.GetYaxis().SetTitle('#sigma(H#rightarrowSS) [pb]')
    axis.GetYaxis().SetRangeUser(ymin, ymax)
    axis.GetXaxis().SetRangeUser(xmin, xmax)
    axis.SetFillColor(4000);
    pads[0].cd()
    pads[0].SetFillStyle(4000)
    pads[0].SetFrameFillStyle(4000)
    axis.Draw('axis')

    ### Create the legend
    legend = r.TLegend(0.17, 0.65, 0.5, 0.8)

    ### Draw the reference graph    
    #StyleLimitBand(referenceLimit, overwrite_style_dict=CMSstyle)
    #DrawLimitBand(pads[0], referenceLimit, draw=['exp0'], legend=None, legend_overwrite=None)
    referenceLimit['exp0'].SetLineColor(r.kBlack)
    referenceLimit['exp0'].SetLineWidth(2)
    referenceLimit['exp0'].Draw('L, SAME')
    legend.AddEntry(referenceLimit['exp0'], limitSet['Reference'][1], 'l')
    
    ### Draw the other graphs
    for graph in graphs:
        graph[0].Draw('L, SAME')
        legend.AddEntry(graph[0], graph[1], 'l')

    legend.Draw()

    # Re-draw the frame and tick marks
    pads[0].RedrawAxis()
    pads[0].GetFrame().Draw()
    pads[0].SetLogx(1)
    pads[0].SetLogy(1)
    pads[0].SetTickx(1)
    pads[0].SetTicky(1)

    # Standard CMS logo
    CMSlabel = r.TLatex()
    CMSlabel.SetNDC();
    CMSlabel.SetTextAngle(0);
    CMSlabel.SetTextColor(r.kBlack);
    CMSlabel.SetTextFont(42);
    CMSlabel.SetTextAlign(22);
    CMSlabel.SetTextSize(0.06);
    CMSlabel.DrawLatex(0.28, 0.89, "#bf{CMS}")
    CMSextralabel = r.TLatex()
    CMSextralabel.SetNDC();
    CMSextralabel.SetTextAngle(0);
    CMSextralabel.SetTextColor(r.kBlack);
    CMSextralabel.SetTextFont(42);
    CMSextralabel.SetTextAlign(22);
    CMSextralabel.SetTextSize(0.04);
    CMSextralabel.DrawLatex(0.28, 0.84, "#it{Internal}")

    # Year label
    Yearlabel = r.TLatex()
    Yearlabel.SetNDC();
    Yearlabel.SetTextAngle(0);
    Yearlabel.SetTextColor(r.kBlack);
    Yearlabel.SetTextFont(42);
    Yearlabel.SetTextAlign(33);
    Yearlabel.SetTextSize(0.04);
    Yearlabel.DrawLatex(0.96, 0.99, lumi)

    CLlabel = r.TLatex()
    CLlabel.SetNDC();
    CLlabel.SetTextAngle(0);
    CLlabel.SetTextColor(r.kBlack);
    CLlabel.SetTextFont(42);
    CLlabel.SetTextAlign(13);
    CLlabel.SetTextSize(0.04);
    CLlabel.DrawLatex(0.16, 0.98, "95% CL upper limits")

    # Model label
    Modellabel = r.TLatex()
    Modellabel.SetNDC();
    Modellabel.SetTextAngle(0);
    Modellabel.SetTextColor(r.kBlack);
    Modellabel.SetTextFont(42);
    Modellabel.SetTextAlign(13);
    Modellabel.SetTextSize(0.037);
    Modellabel.DrawLatex(0.4, 0.9, masslabel[0])
    Modellabel.DrawLatex(0.4, 0.85, masslabel[1])

    # Re-draw axis
    axis.Draw('axis, same')

    canv.Update()
    canv.Modified()
    canv.Print(name + '.pdf')
    canv.Print(name + '.png')

################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()

    ############# HEPData files
    EXO_18_003_mH1000_mm = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H1000_mumu_channel.root'
    EXO_18_003_mH1000_ee = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H1000_ee_channel.root'
    EXO_18_003_mH1000_ll = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H1000_ll_channel.root'
    EXO_18_003_mH800_ll = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H800_ll_channel.root'
    EXO_18_003_mH600_mm = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H600_mumu_channel.root'
    EXO_18_003_mH600_ee = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H600_ee_channel.root'
    EXO_18_003_mH600_ll = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H600_ll_channel.root'
    EXO_18_003_mH400_ll = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H400_ll_channel.root'
    EXO_18_003_mH300_ll = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H300_ll_channel.root'
    EXO_18_003_mH125_ll = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-18-003/HEPData-ins1940976-v1-H125_ll_channel.root'
    EXO_21_006_mH1000_mS350 = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-21-006/HEPData-ins2083735-v2-Figure_17_bottom_right.root'
    EXO_21_006_mH1000_mS150 = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-21-006/HEPData-ins2083735-v2-Figure_17_bottom_left.root'
    EXO_21_006_mH400_mS150 = '/eos/user/f/fernance/LLP_Analysis/comparisons/EXO-21-006/HEPData-ins2083735-v2-Figure_16_bottom.root'


    ############# Set the drawing set
    ### Structure for 'Others'
    ### [filename, GraphName, Label, Color, Style]


    ###########
    ########### Pure flavor channel
    ###########


    ####### mH = 1000 GeV, mS = 350 GeV (muons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_muons_BR100/JSONhiggsCombine__Muon__mH1000__mS350__eraFull_xsec.json', 'Muon channel expected']
    drawSet['Others'] = []
    #drawSet['Others'].append([EXO_18_003_mH1000_mm, 'H(1000) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y3', 'CMS-EXO-18-003 observed (113 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_mm, 'H(1000) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y4', 'CMS-EXO-18-003 expected (113 fb^{-1})', r.kBlue, 1])
    #drawSet['Others'].append([EXO_21_006_mH1000_mS350, 'Figure 17 bottom right/Graph1D_y1', 'CMS-EXO-21-006 observed (97 fb^{-1})', r.kOrange-3, 1])
    drawSet['Others'].append([EXO_21_006_mH1000_mS350, 'Figure 17 bottom right/Graph1D_y2', 'CMS-EXO-21-006 comb. expected (97 fb^{-1})', r.kOrange-3, 1])
    drawSet['Others'].append([EXO_21_006_mH1000_mS350, 'Figure 17 bottom right/Graph1D_y3', 'CMS-EXO-21-006 TMS-TMS expected (97 fb^{-1})', r.kOrange-3, 7])

    HSSlabel = ['H#rightarrowSS, S#rightarrow#mu#mu, B(S#rightarrow#mu#mu) = 100%', 'm_{H} = 1000 GeV, m_{S} = 350 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS350_mm', '96 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e6) 

    ####### mH = 1000 GeV, mS = 150 GeV (muons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_muons_BR100/JSONhiggsCombine__Muon__mH1000__mS150__eraFull_xsec.json', 'Muon channel expected']
    drawSet['Others'] = []
    #drawSet['Others'].append([EXO_18_003_mH1000_mm, 'H(1000) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y1', 'CMS-EXO-18-003 observed (113 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_mm, 'H(1000) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y2', 'CMS-EXO-18-003 expected (113 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_21_006_mH1000_mS150, 'Figure 17 bottom left/Graph1D_y2', 'CMS-EXO-21-006 comb. expected (97 fb^{-1})', r.kOrange-3, 1])
    drawSet['Others'].append([EXO_21_006_mH1000_mS150, 'Figure 17 bottom left/Graph1D_y3', 'CMS-EXO-21-006 TMS-TMS expected (97 fb^{-1})', r.kOrange-3, 7])

    HSSlabel = ['H#rightarrowSS, S#rightarrow#mu#mu, B(S#rightarrow#mu#mu) = 100%', 'm_{H} = 1000 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS150_mm', '96 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e6) 

    ####### mH = 400 GeV, mS = 150 GeV (muons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_muons_BR100/JSONhiggsCombine__Muon__mH400__mS150__eraFull_xsec.json', 'Muon channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_21_006_mH400_mS150, 'Figure 16 bottom/Graph1D_y2', 'CMS-EXO-21-006 comb. expected (97 fb^{-1})', r.kOrange-3, 1])
    drawSet['Others'].append([EXO_21_006_mH400_mS150, 'Figure 16 bottom/Graph1D_y3', 'CMS-EXO-21-006 TMS-TMS expected (97 fb^{-1})', r.kOrange-3, 7])

    HSSlabel = ['H#rightarrowSS, S#rightarrow#mu#mu, B(S#rightarrow#mu#mu) = 100%', 'm_{H} = 400 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH400_mS150_mm', '96 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e6) 

    ####### mH = 1000 GeV, mS = 350 GeV (electrons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_electrons_BR100/JSONhiggsCombine__Electron__mH1000__mS350__eraFull_xsec.json', 'Electron channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH1000_ee, 'H(1000) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y3', 'CMS-EXO-18-003 observed (118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_ee, 'H(1000) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y4', 'CMS-EXO-18-003 expected (118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, S#rightarrowee, B(S#rightarrowee) = 100%', 'm_{H} = 1000 GeV, m_{S} = 350 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS350_ee', '112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 1000 GeV, mS = 150 GeV (electrons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_electrons_BR100/JSONhiggsCombine__Electron__mH1000__mS150__eraFull_xsec.json', 'Electron channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH1000_ee, 'H(1000) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y1', 'CMS-EXO-18-003 observed (118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_ee, 'H(1000) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y2', 'CMS-EXO-18-003 expected (118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, S#rightarrowee, B(S#rightarrowee) = 100%', 'm_{H} = 1000 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS150_ee', '112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 600 GeV, mS = 50 GeV (muons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_muons_BR100/JSONhiggsCombine__Muon__mH600__mS50__eraFull_xsec.json', 'Muon channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_mm, 'H(600) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y1', 'CMS-EXO-18-003 observed (113 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_mm, 'H(600) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y2', 'CMS-EXO-18-003 expected (113 fb^{-1})', r.kBlue, 7])
    HSSlabel = ['H#rightarrowSS, S#rightarrow#mu#mu, B(S#rightarrow#mu#mu) = 100%', 'm_{H} = 600 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS50_mm', '96 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 600 GeV, mS = 150 GeV (muons)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_muons_BR100/JSONhiggsCombine__Muon__mH600__mS150__eraFull_xsec.json', 'Muon channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_mm, 'H(600) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y3', 'CMS-EXO-18-003 observed (113 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_mm, 'H(600) to 4ell, sigmatimes bf limits, mumu channel/Graph1D_y4', 'CMS-EXO-18-003 expected (113 fb^{-1})', r.kBlue, 7])
    HSSlabel = ['H#rightarrowSS, S#rightarrow#mu#mu, B(S#rightarrow#mu#mu) = 100%', 'm_{H} = 600 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS150_mm', '96 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 600 GeV, mS = 50 GeV (electron)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_electrons_BR100/JSONhiggsCombine__Electron__mH600__mS50__eraFull_xsec.json', 'Electron channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_ee, 'H(600) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y1', 'CMS-EXO-18-003 observed (118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_ee, 'H(600) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y2', 'CMS-EXO-18-003 expected (118 fb^{-1})', r.kBlue, 7])
    HSSlabel = ['H#rightarrowSS, S#rightarrowee, B(S#rightarrowee) = 100%', 'm_{H} = 600 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS50_ee', '112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 


    ####### mH = 600 GeV, mS = 150 GeV (electron)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_electrons_BR100/JSONhiggsCombine__Electron__mH600__mS150__eraFull_xsec.json', 'Electron channel expected']
    drawSet['Reference'] = ['LimitsResults_Spring23_electrons_BR100/JSONhiggsCombine__Electron__mH600__mS150__eraFull_xsec.json', 'Electron channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_ee, 'H(600) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y3', 'CMS-EXO-18-003 observed (118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_ee, 'H(600) to 4ell, sigmatimes bf limits, ee channel/Graph1D_y4', 'CMS-EXO-18-003 expected (118 fb^{-1})', r.kBlue, 7])
    HSSlabel = ['H#rightarrowSS, S#rightarrowee, B(S#rightarrowee) = 100%', 'm_{H} = 600 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS150_ee', '112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ###########
    ########### Combined channel
    ###########

    ####### mH = 1000 GeV, mS = 350 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH1000__mS350__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 1000 GeV, m_{S} = 350 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS350_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH1000__mS350__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 1000 GeV, m_{S} = 350 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS350_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 1000 GeV, mS = 150 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH1000__mS150__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y1', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y2', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 1000 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS150_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH1000__mS150__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y1', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH1000_ll, 'H(1000) to SS, Stoell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y2', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 1000 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH1000_mS150_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 800 GeV, mS = 50 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH800__mS50__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y1', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y2', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 800 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH800_mS50_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH800__mS50__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y1', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y2', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 800 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH800_mS50_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 800 GeV, mS = 250 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH800__mS250__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y5', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y6', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 800 GeV, m_{S} = 250 GeV']
    makeLimitComparison('AnalysesComp_mH800_mS250_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH800__mS250__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y5', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH800_ll, 'H(800) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y6', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 800 GeV, m_{S} = 250 GeV']
    makeLimitComparison('AnalysesComp_mH800_mS250_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 


    ####### mH = 600 GeV, mS = 50 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH600__mS50__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y1', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y2', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 600 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS50_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH600__mS50__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y1', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y2', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 600 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS50_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 


    ####### mH = 600 GeV, mS = 150 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH600__mS150__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 600 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS150_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH600__mS150__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH600_ll, 'H(600) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 600 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH600_mS150_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    ####### mH = 400 GeV, mS = 150 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH400__mS150__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH400_ll, 'H(400) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH400_ll, 'H(400) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 400 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH400_mS150_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH400__mS150__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH400_ll, 'H(400) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH400_ll, 'H(400) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 400 GeV, m_{S} = 150 GeV']
    makeLimitComparison('AnalysesComp_mH400_mS150_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 


    ####### mH = 300 GeV, mS = 50 GeV (combined)
    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50/JSONhiggsCombine__Joint__mH300__mS50__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH300_ll, 'H(300) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH300_ll, 'H(300) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 300 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH300_mS50_ll', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 

    drawSet = {}
    drawSet['Reference'] = ['LimitsResults_Spring23_BR50_MuonTight2018/JSONhiggsCombine__Joint__mH300__mS50__eraFull_xsec.json', 'Combined channel expected']
    drawSet['Others'] = []
    drawSet['Others'].append([EXO_18_003_mH300_ll, 'H(300) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y3', 'CMS-EXO-18-003 observed (113-118 fb^{-1})', r.kBlue, 1])
    drawSet['Others'].append([EXO_18_003_mH300_ll, 'H(300) to SS, S toell^{+}ell^{-}, sigmatimes bf limits/Graph1D_y4', 'CMS-EXO-18-003 expected (113-118 fb^{-1})', r.kBlue, 7])

    HSSlabel = ['H#rightarrowSS, B(S#rightarrowee) = B(S#rightarrow#mu#mu) = 50%', 'm_{H} = 300 GeV, m_{S} = 50 GeV']
    makeLimitComparison('AnalysesComp_mH300_mS50_ll_tight', '96-112 fb^{-1} (13 TeV)', drawSet, masslabel = HSSlabel, ymin = 1e-5, ymax = 1e2, xmin = 1e-1, xmax = 1e3) 



