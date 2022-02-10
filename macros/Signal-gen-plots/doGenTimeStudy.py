
import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
import math, sys, optparse, array, copy, os, json
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

###########################
####   Parser object   ####
###########################


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
    (opts, args) = parser.parse_args()


    with open('registry2018.json') as f:
        configs = json.load(f)


    ##################################
    ####   Variable declaration   ####
    ##################################
    Lxy_bin = np.linspace(0.0, 100, 101)
    Lxy_logbin = np.logspace(-2, 5, 101)



    ###################################
    ####   Loop over tree events   ####
    ###################################

    #if not os.path.exists(WORKPATH + 'Results/'): os.makedirs(WORKPATH + 'Results/')
    #outputFile = TFile(WORKPATH + 'Results/th1f.root', 'RECREATE')

    for point in configs.keys():

        configs[point]['rfile'] = r.TFile(configs[point]['file'])
        ttree_ = configs[point]['rfile'].Get("Events")

        configs[point]['r_hist'] = r.TH1F(point + "_r", "; LLP decay radius R (cm); Normalized yield", len(Lxy_logbin)-1, Lxy_logbin)
        configs[point]['r_hist_lin'] = r.TH1F(point + "_r_lin", "; LLP decay radius R (cm); Normalized yield", 100, 0, 100)
        configs[point]['dxy_hist'] = r.TH1F(point + "_dxy", "; Lepton trans. impact parameter |d_{xy}| (cm); Normalized yield", len(Lxy_logbin)-1, Lxy_logbin)
        configs[point]['dxy_hist_lin'] = r.TH1F(point + "_dxy_lin", "; Lepton trans. impact parameter |d_{xy}| (cm); Normalized yield", 100, 0, 100)
        r.SetOwnership(configs[point]["r_hist"], 0)
        r.SetOwnership(configs[point]["r_hist_lin"], 0)
        r.SetOwnership(configs[point]["dxy_hist"], 0)
        r.SetOwnership(configs[point]["dxy_hist_lin"], 0)

        for e_,e in enumerate(ttree_):

            lepton_idx = []
            for n in range(0, e.nlep):
                #print(e.lep_motherIdx[n])
                #if abs(e.lep_motherIdx[n]) != 900006: continue
                configs[point]['r_hist'].Fill(e.lep_Lxy[n])
                configs[point]['r_hist_lin'].Fill(e.lep_Lxy[n])
                configs[point]['dxy_hist'].Fill(abs(e.lep_dxy[n]))
                configs[point]['dxy_hist_lin'].Fill(abs(e.lep_dxy[n]))

        #outputFile.cd()    
        #configs[point]['r_hist'].Write()
        #configs[point]['dxy_hist'].Write()


    print(configs)


    ### Harvesting
    configs['400_50_1']['r_hist'].SetLineWidth(2)
    configs['400_50_10']['r_hist'].SetLineWidth(2)
    configs['400_50_100']['r_hist'].SetLineWidth(2)
    configs['400_50_1000']['r_hist'].SetLineWidth(2)
    configs['400_50_10000']['r_hist'].SetLineWidth(2)
    configs['400_50_1']['r_hist'].SetMaximum(1.4*configs['400_50_1']['r_hist'].GetMaximum())
    plot = Canvas.Canvas('hist_r_400_50', 'png,pdf', 0.14, 0.74, 0.4, 0.89, 1)
    plot.addHisto(configs['400_50_1']['r_hist'].Clone(), 'HIST', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 1 mm', 'p', r.kRed+1, 1, 0)
    plot.addHisto(configs['400_50_10']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 10 mm', 'p', r.kOrange+1, 1, 1)
    plot.addHisto(configs['400_50_100']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 100 mm', 'p', r.kGreen+1, 1, 2)
    plot.addHisto(configs['400_50_1000']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 1000 mm', 'p', r.kBlue+1, 1, 3)
    plot.addHisto(configs['400_50_10000']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 10000 mm', 'p', r.kMagenta+1, 1, 4)
    plot.addLatex(0.9, 0.93, 'H#rightarrow2S#rightarrow4l (2018 UL)', size = 0.027, align = 31)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False, xlog = True, maxYnumbers = 3)


    configs['125_30_1']['r_hist'].SetLineWidth(2)
    configs['125_30_10']['r_hist'].SetLineWidth(2)
    configs['125_30_100']['r_hist'].SetLineWidth(2)
    configs['125_30_1000']['r_hist'].SetLineWidth(2)
    configs['125_30_10000']['r_hist'].SetLineWidth(2)
    configs['125_30_1']['r_hist'].SetMaximum(1.4*configs['125_30_1']['r_hist'].GetMaximum())
    plot = Canvas.Canvas('hist_r_125_30', 'png,pdf', 0.14, 0.74, 0.4, 0.89, 1)
    plot.addHisto(configs['125_30_1']['r_hist'].Clone(), 'HIST', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 1 mm', 'p', r.kRed+1, 1, 0)
    plot.addHisto(configs['125_30_10']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 10 mm', 'p', r.kOrange+1, 1, 1)
    plot.addHisto(configs['125_30_100']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 100 mm', 'p', r.kGreen+1, 1, 2)
    plot.addHisto(configs['125_30_1000']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 1000 mm', 'p', r.kBlue+1, 1, 3)
    plot.addHisto(configs['125_30_10000']['r_hist'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 10000 mm', 'p', r.kMagenta+1, 1, 4)
    plot.addLatex(0.9, 0.93, 'H#rightarrow2S#rightarrow4l (2018 UL)', size = 0.027, align = 31)
    plot.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False, xlog = True, maxYnumbers = 3)



    configs['400_50_1']['dxy_hist_lin'].SetLineWidth(2)
    configs['400_50_10']['dxy_hist_lin'].SetLineWidth(2)
    configs['400_50_100']['dxy_hist_lin'].SetLineWidth(2)
    configs['400_50_1000']['dxy_hist_lin'].SetLineWidth(2)
    configs['400_50_10000']['dxy_hist_lin'].SetLineWidth(2)
    configs['400_50_1']['dxy_hist_lin'].SetMaximum(10*configs['400_50_1']['dxy_hist_lin'].GetMaximum())
    plot = Canvas.Canvas('hist_dxy_lin_400_50', 'png,pdf', 0.14, 0.74, 0.4, 0.89, 1)
    plot.addHisto(configs['400_50_1']['dxy_hist_lin'].Clone(), 'HIST', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 1 mm', 'p', r.kRed+1, 1, 0, doOF = True)
    plot.addHisto(configs['400_50_10']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 10 mm', 'p', r.kOrange+1, 1, 1, doOF = True)
    plot.addHisto(configs['400_50_100']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 100 mm', 'p', r.kGreen+1, 1, 2, doOF = True)
    plot.addHisto(configs['400_50_1000']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 1000 mm', 'p', r.kBlue+1, 1, 3, doOF = True)
    plot.addHisto(configs['400_50_10000']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 400 GeV, m_{S} = 50 GeV, c#tau = 10000 mm', 'p', r.kMagenta+1, 1, 4, doOF = True)
    plot.addLatex(0.9, 0.93, 'H#rightarrow2S#rightarrow4l (2018 UL)', size = 0.027, align = 31)
    plot.save(1, 0, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)


    configs['125_30_1']['dxy_hist_lin'].SetLineWidth(2)
    configs['125_30_10']['dxy_hist_lin'].SetLineWidth(2)
    configs['125_30_100']['dxy_hist_lin'].SetLineWidth(2)
    configs['125_30_1000']['dxy_hist_lin'].SetLineWidth(2)
    configs['125_30_10000']['dxy_hist_lin'].SetLineWidth(2)
    configs['125_30_1']['dxy_hist_lin'].SetMaximum(10*configs['125_30_1']['dxy_hist_lin'].GetMaximum())
    plot = Canvas.Canvas('hist_dxy_lin_125_30', 'png,pdf', 0.14, 0.74, 0.4, 0.89, 1)
    plot.addHisto(configs['125_30_1']['dxy_hist_lin'].Clone(), 'HIST', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 1 mm', 'p', r.kRed+1, 1, 0, doOF = True)
    plot.addHisto(configs['125_30_10']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 10 mm', 'p', r.kOrange+1, 1, 1, doOF = True)
    plot.addHisto(configs['125_30_100']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 100 mm', 'p', r.kGreen+1, 1, 2, doOF = True)
    plot.addHisto(configs['125_30_1000']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 1000 mm', 'p', r.kBlue+1, 1, 3, doOF = True)
    plot.addHisto(configs['125_30_10000']['dxy_hist_lin'].Clone(), 'HIST, SAME', 'm_{H} = 125 GeV, m_{S} = 30 GeV, c#tau = 10000 mm', 'p', r.kMagenta+1, 1, 4, doOF = True)
    plot.addLatex(0.9, 0.93, 'H#rightarrow2S#rightarrow4l (2018 UL)', size = 0.027, align = 31)
    plot.save(1, 0, 1, '', '', outputDir = WORKPATH + 'harvested/', inProgress = False)



    #outputFile.Close()
