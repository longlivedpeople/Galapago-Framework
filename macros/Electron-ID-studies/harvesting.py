import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership, TVector3
import math, sys, optparse, array, copy, os
import gc, inspect
import numpy as np



################################# GLOBAL VARIABLES DEFINITION ####################################

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

print(WORKPATH, WORKPATH)
print(GALAPAGOPATH, GALAPAGOPATH)

import include.Canvas as Canvas

##################################### FUNCTION DEFINITION ########################################

def getObject(filename, key):

    _f = r.TFile(filename)
    _h = _f.Get(key)
    _hcopy = copy.deepcopy(_h)
    _f.Close()

    return _hcopy


def rebinAxis(eff, axis):

    totals = eff.GetTotalHistogram()
    passed = eff.GetPassedHistogram()
    totals_rebin = totals.Rebin(len(axis)-1, totals.GetName()+'_rebined', axis)
    passed_rebin = passed.Rebin(len(axis)-1, passed.GetName()+'_rebined', axis)
    neweff = r.TEfficiency(passed_rebin, totals_rebin)

    c1 = r.TCanvas()
    neweff.Draw('AP')
#    c1.SaveAs('Pruebita.png')

    return(neweff)

if __name__ == "__main__":


    gROOT.ProcessLine('.L ' + GALAPAGOPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    print('WORKPATH: ' + WORKPATH)
    r.setTDRStyle()

    ###########################
    ####   Parser object   ####
    ###########################
    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-t', '--tag', action='store', type=str, dest='tag', default='', help='Output tag')
    (opts, args) = parser.parse_args()

    ##############################
    ####   Some definitions   ####
    ##############################
    ptbin = np.concatenate((np.arange(0, 90, 5, float), np.arange(90, 130, 10), np.arange(130, 190, 15), np.arange(190, 260, 40), np.array([300])))

    #####################################
    ####   Construct TEfficiencies   ####
    #####################################
    
    ### Pt (IFCA)

    BKG_IFCAeff_pt = getObject('Results/th1fDrell_Yan.root', 'eff_IFCA_pt')
    SI_4mm_IFCAeff_pt = getObject('Results/th1fHXX_400_50_4mm.root', 'eff_IFCA_pt')
    SI_40mm_IFCAeff_pt = getObject('Results/th1fHXX_400_50_40mm.root', 'eff_IFCA_pt')
    SI_400mm_IFCAeff_pt = getObject('Results/th1fHXX_400_50_400mm.root', 'eff_IFCA_pt')

    BKG_IFCAeff_pt_ = rebinAxis(BKG_IFCAeff_pt, ptbin)
    SI_4mm_IFCAeff_pt_ = rebinAxis(SI_4mm_IFCAeff_pt, ptbin)
    SI_40mm_IFCAeff_pt_ = rebinAxis(SI_40mm_IFCAeff_pt, ptbin)
    SI_400mm_IFCAeff_pt_ = rebinAxis(SI_400mm_IFCAeff_pt, ptbin)

    EFF_pt = Canvas.Canvas("EFF_IFCA_pt", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    EFF_pt.addRate(SI_4mm_IFCAeff_pt_, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_pt.addRate(SI_40mm_IFCAeff_pt_, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    EFF_pt.addRate(SI_400mm_IFCAeff_pt_, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    EFF_pt.addRate(BKG_IFCAeff_pt_, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    EFF_pt.addLatex(0.9, 0.93, 'Displaced electrons', size = 0.032, align = 31)
    EFF_pt.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)
    
    ### Pt (CMS)

    BKG_CMSeff_pt = getObject('Results/th1fDrell_Yan.root', 'eff_CMS_pt')
    SI_4mm_CMSeff_pt = getObject('Results/th1fHXX_400_50_4mm.root', 'eff_CMS_pt')
    SI_40mm_CMSeff_pt = getObject('Results/th1fHXX_400_50_40mm.root', 'eff_CMS_pt')
    SI_400mm_CMSeff_pt = getObject('Results/th1fHXX_400_50_400mm.root', 'eff_CMS_pt')
    BKG_CMSeff_pt_ = rebinAxis(BKG_CMSeff_pt, ptbin)
    SI_4mm_CMSeff_pt_ = rebinAxis(SI_4mm_CMSeff_pt, ptbin)
    SI_40mm_CMSeff_pt_ = rebinAxis(SI_40mm_CMSeff_pt, ptbin)
    SI_400mm_CMSeff_pt_ = rebinAxis(SI_400mm_CMSeff_pt, ptbin)

    EFF_pt = Canvas.Canvas("EFF_CMS_pt", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    EFF_pt.addRate(SI_4mm_CMSeff_pt_, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_pt.addRate(SI_40mm_CMSeff_pt_, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    EFF_pt.addRate(SI_400mm_CMSeff_pt_, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    EFF_pt.addRate(BKG_CMSeff_pt_, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    EFF_pt.addLatex(0.9, 0.93, 'PAT Loose electrons', size = 0.032, align = 31)
    EFF_pt.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)

    ### Eta (IFCA)

    BKG_IFCAeff_eta = getObject('Results/th1fDrell_Yan.root', 'eff_IFCA_eta')
    SI_4mm_IFCAeff_eta = getObject('Results/th1fHXX_400_50_4mm.root', 'eff_IFCA_eta')
    SI_40mm_IFCAeff_eta = getObject('Results/th1fHXX_400_50_40mm.root', 'eff_IFCA_eta')
    SI_400mm_IFCAeff_eta = getObject('Results/th1fHXX_400_50_400mm.root', 'eff_IFCA_eta')

    EFF_eta = Canvas.Canvas("EFF_IFCA_eta", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    EFF_eta.addRate(SI_4mm_IFCAeff_eta, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_eta.addRate(SI_40mm_IFCAeff_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    EFF_eta.addRate(SI_400mm_IFCAeff_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    EFF_eta.addRate(BKG_IFCAeff_eta, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    EFF_eta.addLatex(0.9, 0.93, 'Displaced electrons', size = 0.032, align = 31)
    EFF_eta.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)


    ### Eta (CMS)

    BKG_CMSeff_eta = getObject('Results/th1fDrell_Yan.root', 'eff_CMS_eta')
    SI_4mm_CMSeff_eta = getObject('Results/th1fHXX_400_50_4mm.root', 'eff_CMS_eta')
    SI_40mm_CMSeff_eta = getObject('Results/th1fHXX_400_50_40mm.root', 'eff_CMS_eta')
    SI_400mm_CMSeff_eta = getObject('Results/th1fHXX_400_50_400mm.root', 'eff_CMS_eta')

    EFF_eta = Canvas.Canvas("EFF_CMS_eta", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    EFF_eta.addRate(SI_4mm_CMSeff_eta, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_eta.addRate(SI_40mm_CMSeff_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    EFF_eta.addRate(SI_400mm_CMSeff_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    EFF_eta.addRate(BKG_CMSeff_eta, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    EFF_eta.addLatex(0.9, 0.93, 'PAT Loose electrons', size = 0.032, align = 31)
    EFF_eta.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)

    ### dxy (IFCA)
        
    SI_50_400mm_IFCAeff_dxy = getObject('Results/th1fHXX_400_50_400mm.root', 'eff_IFCA_dxy_log')
    SI_150_400mm_IFCAeff_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_IFCA_dxy_log')

    EFF_dxy = Canvas.Canvas("EFF_IFCA_dxy", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    EFF_dxy.addRate(SI_50_400mm_IFCAeff_dxy, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_dxy.addRate(SI_150_400mm_IFCAeff_dxy, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kRed+2, True, 1, marker = 24)
    EFF_dxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False, xlog = True)


    ### dxy (CMS)
        
    SI_50_400mm_CMSeff_dxy = getObject('Results/th1fHXX_400_50_400mm.root', 'eff_CMS_dxy_log')
    SI_150_400mm_CMSeff_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_CMS_dxy_log')

    EFF_dxy = Canvas.Canvas("EFF_CMS_dxy", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    EFF_dxy.addRate(SI_50_400mm_CMSeff_dxy, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_dxy.addRate(SI_150_400mm_CMSeff_dxy, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kRed+2, True, 1, marker = 24)
    EFF_dxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False, xlog = True)

    ### dxy (Combined IFCA-CMS log)
        
    SI_50_400mm_IFCAeff_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_IFCA_dxy_log')
    SI_50_400mm_CMSeff_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_CMS_dxy_log')

    logbin_dxy = np.logspace(-3, 3, 70)
    #for n in range(0, len(logbin_dxy)):
    #    print(n, logbin_dxy[n])
    logrebin_dxy = np.concatenate((logbin_dxy[[0]], logbin_dxy[20:]))

    SI_50_400mm_IFCAeff_dxy = rebinAxis(SI_50_400mm_IFCAeff_dxy, logrebin_dxy)
    SI_50_400mm_CMSeff_dxy = rebinAxis(SI_50_400mm_CMSeff_dxy, logrebin_dxy)

    EFF_dxy = Canvas.Canvas("EFF_Combined_log_dxy", 'png', 0.15, 0.70, 0.45, 0.78, 1) 
    EFF_dxy.addRate(SI_50_400mm_IFCAeff_dxy, 'AP', 'Displaced electrons', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_dxy.addRate(SI_50_400mm_CMSeff_dxy, 'AP,SAME', 'PAT Loose electrons', 'p', r.kRed+2, True, 1, marker = 24)
    EFF_dxy.addLatex(0.15, 0.86, 'm_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', size = 0.032, align = 11)
    EFF_dxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False, xlog = True)

    ### dxy (Combined IFCA-CMS lin)
        
    SI_50_400mm_IFCAeff_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_IFCA_dxy')
    SI_50_400mm_CMSeff_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_CMS_dxy')

    EFF_dxy = Canvas.Canvas("EFF_Combined_lin_dxy", 'png', 0.15, 0.70, 0.45, 0.78, 1) 
    EFF_dxy.addRate(SI_50_400mm_IFCAeff_dxy, 'AP', 'Displaced electrons', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_dxy.addRate(SI_50_400mm_CMSeff_dxy, 'AP,SAME', 'PAT Loose electrons', 'p', r.kRed+2, True, 1, marker = 24)
    EFF_dxy.addLatex(0.15, 0.86, 'm_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', size = 0.032, align = 11)
    EFF_dxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)

    
    ### Lxy (Combined IFCA-CMS log)

    logbin_Lxy = np.logspace(-3, 3, 70)
    #for n in range(0, len(logbin_Lxy)):
    #    print(n, logbin_Lxy[n])
    logrebin_Lxy = np.concatenate((logbin_Lxy[[0, 23, 27, 30]], logbin_Lxy[31:]))
    SI_50_400mm_IFCAeff_Lxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_IFCA_Lxy_log')
    SI_50_400mm_CMSeff_Lxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_CMS_Lxy_log')

    SI_50_400mm_IFCAeff_Lxy = rebinAxis(SI_50_400mm_IFCAeff_Lxy, logrebin_Lxy)
    SI_50_400mm_CMSeff_Lxy = rebinAxis(SI_50_400mm_CMSeff_Lxy, logrebin_Lxy)

    EFF_Lxy = Canvas.Canvas("EFF_Combined_log_Lxy", 'png', 0.15, 0.70, 0.45, 0.78, 1) 
    EFF_Lxy.addRate(SI_50_400mm_IFCAeff_Lxy, 'AP', 'Displaced electrons', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_Lxy.addRate(SI_50_400mm_CMSeff_Lxy, 'AP,SAME', 'PAT Loose electrons', 'p', r.kRed+2, True, 1, marker = 24)
    EFF_Lxy.addLatex(0.15, 0.86, 'm_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', size = 0.032, align = 11)
    EFF_Lxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False, xlog = True)

    ### Lxy (Combined IFCA-CMS lin)
        
    SI_50_400mm_IFCAeff_Lxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_IFCA_Lxy')
    SI_50_400mm_CMSeff_Lxy = getObject('Results/th1fHXX_400_150_400mm.root', 'eff_CMS_Lxy')

    EFF_Lxy = Canvas.Canvas("EFF_Combined_lin_Lxy", 'png', 0.15, 0.7, 0.45, 0.78, 1) 
    EFF_Lxy.addRate(SI_50_400mm_IFCAeff_Lxy, 'AP', 'Displaced electrons', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_Lxy.addRate(SI_50_400mm_CMSeff_Lxy, 'AP,SAME', 'PAT Loose electrons', 'p', r.kRed+2, True, 1, marker = 24)
    EFF_Lxy.addLatex(0.15, 0.86, 'm_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', size = 0.032, align = 11)
    EFF_Lxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)
    


    ######################################
    ####   Construct Fake Fractions   ####
    ######################################
    
    ### Pt

    BKG_IFCAfake_pt = getObject('Results/th1fDrell_Yan.root', 'fake_IFCA_pt')
    SI_4mm_IFCAfake_pt = getObject('Results/th1fHXX_400_50_4mm.root', 'fake_IFCA_pt')
    SI_40mm_IFCAfake_pt = getObject('Results/th1fHXX_400_50_40mm.root', 'fake_IFCA_pt')
    SI_400mm_IFCAfake_pt = getObject('Results/th1fHXX_400_50_400mm.root', 'fake_IFCA_pt')

    FAKE_pt = Canvas.Canvas("FAKE_pt", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    FAKE_pt.addRate(SI_4mm_IFCAfake_pt, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    FAKE_pt.addRate(SI_40mm_IFCAfake_pt, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    FAKE_pt.addRate(SI_400mm_IFCAfake_pt, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    FAKE_pt.addRate(BKG_IFCAfake_pt, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    FAKE_pt.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)
    
    ### Eta

    BKG_IFCAfake_eta = getObject('Results/th1fDrell_Yan.root', 'fake_IFCA_eta')
    SI_4mm_IFCAfake_eta = getObject('Results/th1fHXX_400_50_4mm.root', 'fake_IFCA_eta')
    SI_40mm_IFCAfake_eta = getObject('Results/th1fHXX_400_50_40mm.root', 'fake_IFCA_eta')
    SI_400mm_IFCAfake_eta = getObject('Results/th1fHXX_400_50_400mm.root', 'fake_IFCA_eta')

    FAKE_eta = Canvas.Canvas("FAKE_eta", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    FAKE_eta.addRate(SI_4mm_IFCAfake_eta, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    FAKE_eta.addRate(SI_40mm_IFCAfake_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    FAKE_eta.addRate(SI_400mm_IFCAfake_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    FAKE_eta.addRate(BKG_IFCAfake_eta, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    FAKE_eta.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)



    ### dxy
        
    SI_50_400mm_IFCAfake_dxy = getObject('Results/th1fHXX_400_50_400mm.root', 'fake_IFCA_dxy')
    SI_150_400mm_IFCAfake_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'fake_IFCA_dxy')

    FAKE_dxy = Canvas.Canvas("FAKE_dxy", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    FAKE_dxy.addRate(SI_50_400mm_IFCAfake_dxy, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    FAKE_dxy.addRate(SI_150_400mm_IFCAfake_dxy, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kRed+2, True, 1, marker = 24)
    FAKE_dxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)


    ### Pt (TTbar)

    BKG_IFCAfake_pt = getObject('Results/th1fTTbar.root', 'fake_IFCA_pt')
    SI_4mm_IFCAfake_pt = getObject('Results/th1fHXX_400_50_4mm.root', 'fake_IFCA_pt')
    SI_40mm_IFCAfake_pt = getObject('Results/th1fHXX_400_50_40mm.root', 'fake_IFCA_pt')
    SI_400mm_IFCAfake_pt = getObject('Results/th1fHXX_400_50_400mm.root', 'fake_IFCA_pt')

    FAKE_pt = Canvas.Canvas("FAKE_TTbar_pt", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    FAKE_pt.addRate(SI_4mm_IFCAfake_pt, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    FAKE_pt.addRate(SI_40mm_IFCAfake_pt, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    FAKE_pt.addRate(SI_400mm_IFCAfake_pt, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    FAKE_pt.addRate(BKG_IFCAfake_pt, 'AP,SAME', 'TTbar', 'p', r.kBlack, True, 3, marker = 24)
    FAKE_pt.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)
    
    ### Eta (TTbar)

    BKG_IFCAfake_eta = getObject('Results/th1fTTbar.root', 'fake_IFCA_eta')
    SI_4mm_IFCAfake_eta = getObject('Results/th1fHXX_400_50_4mm.root', 'fake_IFCA_eta')
    SI_40mm_IFCAfake_eta = getObject('Results/th1fHXX_400_50_40mm.root', 'fake_IFCA_eta')
    SI_400mm_IFCAfake_eta = getObject('Results/th1fHXX_400_50_400mm.root', 'fake_IFCA_eta')

    FAKE_eta = Canvas.Canvas("FAKE_TTbar_eta", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    FAKE_eta.addRate(SI_4mm_IFCAfake_eta, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    FAKE_eta.addRate(SI_40mm_IFCAfake_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    FAKE_eta.addRate(SI_400mm_IFCAfake_eta, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    FAKE_eta.addRate(BKG_IFCAfake_eta, 'AP,SAME', 'TTbar', 'p', r.kBlack, True, 3, marker = 24)
    FAKE_eta.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)



    ### dxy (TTbar)
        
    BKG_IFCAfake_dxy = getObject('Results/th1fTTbar.root', 'fake_IFCA_dxy')
    SI_50_400mm_IFCAfake_dxy = getObject('Results/th1fHXX_400_50_400mm.root', 'fake_IFCA_dxy')
    SI_150_400mm_IFCAfake_dxy = getObject('Results/th1fHXX_400_150_400mm.root', 'fake_IFCA_dxy')

    FAKE_dxy = Canvas.Canvas("FAKE_TTbar_dxy", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    FAKE_dxy.addRate(SI_50_400mm_IFCAfake_dxy, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    FAKE_dxy.addRate(SI_150_400mm_IFCAfake_dxy, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kRed+2, True, 1, marker = 24)
    FAKE_dxy.addRate(BKG_IFCAfake_dxy, 'AP,SAME', 'TTbar', 'p', r.kBlack, True, 1, marker = 24)
    FAKE_dxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)




    ########################
    ####   Histograms   ####
    ########################

    ### Pt resolution (IFCA)

    BKG_IFCAhist_ptRes = getObject('Results/th1fDrell_Yan.root', 'hist_IFCA_ptRes')
    SI_4mm_IFCAhist_ptRes = getObject('Results/th1fHXX_400_50_4mm.root', 'hist_IFCA_ptRes')
    SI_40mm_IFCAhist_ptRes = getObject('Results/th1fHXX_400_50_40mm.root', 'hist_IFCA_ptRes')
    SI_400mm_IFCAhist_ptRes = getObject('Results/th1fHXX_400_50_400mm.root', 'hist_IFCA_ptRes')

    SI_4mm_IFCAhist_ptRes.Scale(1.0/SI_4mm_IFCAhist_ptRes.GetEntries())
    SI_40mm_IFCAhist_ptRes.Scale(1.0/SI_40mm_IFCAhist_ptRes.GetEntries())
    SI_400mm_IFCAhist_ptRes.Scale(1.0/SI_400mm_IFCAhist_ptRes.GetEntries())
    BKG_IFCAhist_ptRes.Scale(1.0/BKG_IFCAhist_ptRes.GetEntries())

    SI_4mm_IFCAhist_ptRes.SetMaximum(1.4*max(BKG_IFCAhist_ptRes.GetMaximum(), SI_4mm_IFCAhist_ptRes.GetMaximum(), SI_40mm_IFCAhist_ptRes.GetMaximum(), SI_400mm_IFCAhist_ptRes.GetMaximum()))

    HIST_ptRes = Canvas.Canvas("HIST_IFCA_photon_ptRes", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    HIST_ptRes.addRate(SI_4mm_IFCAhist_ptRes, 'HIST', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'f', r.kBlue+2, True, 0, marker = 24)
    HIST_ptRes.addRate(SI_40mm_IFCAhist_ptRes, 'HIST,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'f', r.kBlue-4, True, 1, marker = 24)
    HIST_ptRes.addRate(SI_400mm_IFCAhist_ptRes, 'HIST,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'f', r.kBlue-6, True, 2, marker = 24)
    HIST_ptRes.addRate(BKG_IFCAhist_ptRes, 'HIST,SAME', 'Z/#gamma* #rightarrow ee', 'f', r.kBlack, True, 3, marker = 24)
    HIST_ptRes.addLatex(0.9, 0.93, 'Displaced electrons', size = 0.032, align = 31)
    HIST_ptRes.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)


    ### Pt resolution (CMS)

    BKG_CMShist_ptRes = getObject('Results/th1fDrell_Yan.root', 'hist_CMS_ptRes')
    SI_4mm_CMShist_ptRes = getObject('Results/th1fHXX_400_50_4mm.root', 'hist_CMS_ptRes')
    SI_40mm_CMShist_ptRes = getObject('Results/th1fHXX_400_50_40mm.root', 'hist_CMS_ptRes')
    SI_400mm_CMShist_ptRes = getObject('Results/th1fHXX_400_50_400mm.root', 'hist_CMS_ptRes')

    SI_4mm_CMShist_ptRes.Scale(1.0/SI_4mm_CMShist_ptRes.GetEntries())
    SI_40mm_CMShist_ptRes.Scale(1.0/SI_40mm_CMShist_ptRes.GetEntries())
    SI_400mm_CMShist_ptRes.Scale(1.0/SI_400mm_CMShist_ptRes.GetEntries())
    BKG_CMShist_ptRes.Scale(1.0/BKG_CMShist_ptRes.GetEntries())

    SI_4mm_CMShist_ptRes.SetMaximum(1.4*max(BKG_CMShist_ptRes.GetMaximum(), SI_4mm_CMShist_ptRes.GetMaximum(), SI_40mm_CMShist_ptRes.GetMaximum(), SI_400mm_CMShist_ptRes.GetMaximum()))
    #SI_4mm_CMShist_ptRes.SetMaximum(1.4*SI_4mm_CMShist_ptRes.GetMaximum())

    HIST_ptRes = Canvas.Canvas("HIST_CMS_ptRes", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    HIST_ptRes.addRate(SI_4mm_CMShist_ptRes, 'HIST', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'f', r.kBlue+2, True, 0, marker = 24)
    HIST_ptRes.addRate(SI_40mm_CMShist_ptRes, 'HIST,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'f', r.kBlue-4, True, 1, marker = 24)
    HIST_ptRes.addRate(SI_400mm_CMShist_ptRes, 'HIST,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'f', r.kBlue-6, True, 2, marker = 24)
    HIST_ptRes.addRate(BKG_CMShist_ptRes, 'HIST,SAME', 'Z/#gamma* #rightarrow ee', 'f', r.kBlack, True, 3, marker = 24)
    HIST_ptRes.addLatex(0.9, 0.93, 'PAT Loose electrons', size = 0.032, align = 31)
    HIST_ptRes.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)


    ### dxy resolution

    BKG_IFCAhist_dxyRes = getObject('Results/th1fDrell_Yan.root', 'hist_IFCA_dxyRes')
    SI_4mm_IFCAhist_dxyRes = getObject('Results/th1fHXX_400_50_4mm.root', 'hist_IFCA_dxyRes')
    SI_40mm_IFCAhist_dxyRes = getObject('Results/th1fHXX_400_50_40mm.root', 'hist_IFCA_dxyRes')
    SI_400mm_IFCAhist_dxyRes = getObject('Results/th1fHXX_400_50_400mm.root', 'hist_IFCA_dxyRes')

    ### dxy resolution

    BKG_IFCAhist_dxyRes = getObject('Results/th1fDrell_Yan.root', 'hist_IFCA_dxyRes')
    SI_4mm_IFCAhist_dxyRes = getObject('Results/th1fHXX_400_50_4mm.root', 'hist_IFCA_dxyRes')
    SI_40mm_IFCAhist_dxyRes = getObject('Results/th1fHXX_400_50_40mm.root', 'hist_IFCA_dxyRes')
    SI_400mm_IFCAhist_dxyRes = getObject('Results/th1fHXX_400_50_400mm.root', 'hist_IFCA_dxyRes')

    SI_4mm_IFCAhist_dxyRes.Scale(1.0/SI_4mm_IFCAhist_dxyRes.Integral())
    SI_40mm_IFCAhist_dxyRes.Scale(1.0/SI_40mm_IFCAhist_dxyRes.Integral())
    SI_400mm_IFCAhist_dxyRes.Scale(1.0/SI_400mm_IFCAhist_dxyRes.Integral())

    SI_4mm_IFCAhist_dxyRes.SetMaximum(1.4*SI_4mm_IFCAhist_dxyRes.GetMaximum())

    HIST_dxyRes = Canvas.Canvas("HIST_dxyRes", 'png', 0.15, 0.77, 0.45, 0.9, 1) 
    HIST_dxyRes.addRate(SI_4mm_IFCAhist_dxyRes, 'HIST', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'f', r.kBlue+2, True, 0, marker = 24)
    HIST_dxyRes.addRate(SI_40mm_IFCAhist_dxyRes, 'HIST,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'f', r.kBlue-4, True, 1, marker = 24)
    HIST_dxyRes.addRate(SI_400mm_IFCAhist_dxyRes, 'HIST,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'f', r.kBlue-6, True, 2, marker = 24)
    HIST_dxyRes.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)

    #############################33
    ### Charge assignment vs Pt (IFCA)

    BKG_IFCAeff_pt = getObject('Results/th1fDrell_Yan.root', 'qeff_IFCA_pt')
    SI_4mm_IFCAeff_pt = getObject('Results/th1fHXX_400_50_4mm.root', 'qeff_IFCA_pt')
    SI_40mm_IFCAeff_pt = getObject('Results/th1fHXX_400_50_40mm.root', 'qeff_IFCA_pt')
    SI_400mm_IFCAeff_pt = getObject('Results/th1fHXX_400_50_400mm.root', 'qeff_IFCA_pt')

    BKG_IFCAeff_pt_ = rebinAxis(BKG_IFCAeff_pt, ptbin)
    SI_4mm_IFCAeff_pt_ = rebinAxis(SI_4mm_IFCAeff_pt, ptbin)
    SI_40mm_IFCAeff_pt_ = rebinAxis(SI_40mm_IFCAeff_pt, ptbin)
    SI_400mm_IFCAeff_pt_ = rebinAxis(SI_400mm_IFCAeff_pt, ptbin)

    EFF_pt = Canvas.Canvas("QEFF_IFCA_pt", 'png', 0.15, 0.15, 0.45, 0.28, 1) 
    EFF_pt.addRate(SI_4mm_IFCAeff_pt_, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 4 mm', 'p', r.kBlue+2, True, 0, marker = 24)
    EFF_pt.addRate(SI_40mm_IFCAeff_pt_, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 40 mm', 'p', r.kBlue-4, True, 1, marker = 24)
    EFF_pt.addRate(SI_400mm_IFCAeff_pt_, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 2, marker = 24)
    EFF_pt.addRate(BKG_IFCAeff_pt_, 'AP,SAME', 'Z/#gamma* #rightarrow ee', 'p', r.kBlack, True, 3, marker = 24)
    EFF_pt.addLatex(0.9, 0.93, 'Displaced electrons', size = 0.032, align = 31)
    EFF_pt.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)



    ### Charge assignment vs Lxy (IFCA)

    SI_150mm_IFCAeff_Lxy = getObject('Results/th1fHXX_400_150_400mm.root', 'qeff_IFCA_Lxy')
    SI_50mm_IFCAeff_Lxy = getObject('Results/th1fHXX_400_50_400mm.root', 'qeff_IFCA_Lxy')

    SI_150mm_IFCAeff_Lxy = rebinAxis(SI_150mm_IFCAeff_Lxy, np.linspace(0.0, 70, 26))
    SI_50mm_IFCAeff_Lxy = rebinAxis(SI_50mm_IFCAeff_Lxy, np.linspace(0.0, 70, 26))

    SI_150mm_IFCAeff_Lxy.SetTitle(';Generated e L_{xy} (cm);Charge assignment efficiency')
    SI_50mm_IFCAeff_Lxy.SetTitle(';;')

    EFF_Lxy = Canvas.Canvas("QEFF_IFCA_Lxy", 'png', 0.15, 0.20, 0.45, 0.28, 1) 
    EFF_Lxy.addRate(SI_150mm_IFCAeff_Lxy, 'AP', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 150 GeV, c#tau = 400 mm', 'p', r.kMagenta+2, True, 0, marker = 24)
    EFF_Lxy.addRate(SI_50mm_IFCAeff_Lxy, 'AP,SAME', 'H#rightarrowXX: m_{H} = 400 GeV, m_{X} = 50 GeV, c#tau = 400 mm', 'p', r.kBlue-6, True, 1, marker = 24)
    EFF_Lxy.addLatex(0.9, 0.93, 'Displaced electrons', size = 0.032, align = 31)
    EFF_Lxy.save(1, 0, 0, '', '', outputDir = WORKPATH + 'harvested_'+opts.tag+'/', inProgress = False)

