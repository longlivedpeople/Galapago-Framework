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

xsecs = {}
xsecs['400'] = 5040. # pb-1

### Define palettes
bkgc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),  
                                            np.array([1.00, 0.00]),  
                                            np.array([1.00, 153/255.]), 
                                            np.array([1.00, 153./255.]), 255);

sigc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),  
                                            np.array([1.00, 204./255.]),  
                                            np.array([1.00, 0.00]), 
                                            np.array([1.00, 0.00]), 255);
bkgpalette_ = []
sigpalette_ = []
for i in range(0, 255):
    bkgpalette_.append(bkgc + i)
    sigpalette_.append(sigc + i)

bkgpalette = np.array(bkgpalette_, dtype=np.int32)
sigpalette = np.array(sigpalette_, dtype=np.int32)

#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###     - xsec to normalize the signal given in fb
###
#####
#####################

def makeBackgroundPlot2D(lumi, hname_bkg, zlog, treeDATA, inputdir, rebin = False, lines = [], xlabel = '', outtag = '', outdir = '', LLlabel = '', extralabel = ''):

    ### Get histograms
    luminosity = lumi

    hbkg = treeDATA.getLoopTH2F(inputdir, hname_bkg)
    
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)
    hbkg.SetMinimum(0.1)

    ### Define palette
    nc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),  
                                              np.array([1.00, 0.00]),  
                                              np.array([1.00, 153/255.]), 
                                              np.array([1.00, 153./255.]), 255);
    r.gStyle.SetPalette(57)

    ### Canvas object
    plot = Canvas.Canvas(outtag+hname_bkg, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, ww = 610, hh = 600, lsize = 0.028)
    plot.addHisto(hbkg, 'COLZ', '', '', '', 1, 0)
    
    for line in lines:
        plot.addLine(line, hbkg.GetMinimum(), line, hbkg.GetMaximum(), r.kBlack, 2)

    ### Extralabel
    plot.addLatex(0.5, 0.83, 'Background (predicted)', font = 42)
    plot.addLatex(0.5, 0.77, extralabel, font = 42, size = 0.03)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, 0, luminosity, '', outputDir = outdir, zlog = zlog, is2d = True)


def makeDataMCPlot(lumi, hname_DATA, ylog, treeDATA, inputdir, xlabel = '', outtag = '', outdir = '', yshift = 100.0, LLlabel = '', leftlabel = '', xlog = False):

    luminosity = lumi
    hDATA    = copy.deepcopy(treeDATA.getLoopTH1F(inputdir, hname_DATA))


    ### Tune the histograms
    hDATA.SetMarkerStyle(20)
    hDATA.SetMarkerSize(0.8)

    ### Get maximum
    maxVal = hDATA.GetMaximum()

    ### Set Maximum
    if not ylog:
        if treeDATA:
            if not yshift:
                hDATA.SetMaximum(1.3*maxVal)
            else:
                hDATA.SetMaximum(yshift*maxVal)
            hDATA.SetMinimum(0.0)
    else:
        if treeDATA:
            if not yshift:
                hDATA.SetMaximum(10.0*maxVal)
            else:
                hDATA.SetMaximum(yshift*maxVal)
            hDATA.SetMinimum(0.1)

    ### SetOwnership
    r.SetOwnership(hDATA, 0)

    ### -> Canvas object
    plot = Canvas.Canvas(outtag+hname_DATA, 'png,pdf', 0.6, 0.6, 0.85, 0.84, 1)

    ### Add lines for detector structure
    plot.addLine(2.21,1e-1,2.21,hDATA.GetMaximum()/100,r.kRed,thickness=1)
    plot.addBand(2.01,1e-1,2.41,hDATA.GetMaximum()/100,r.kRed,0.2)
    plot.addLine(2.95,1e-1,2.95,hDATA.GetMaximum()/100,r.kRed,thickness=1)
    plot.addBand(2.80,1e-1,3.10,hDATA.GetMaximum()/100,r.kRed,0.2)
    plot.addLine(6.80,1e-1,6.80,hDATA.GetMaximum()/100,r.kRed,thickness=1)
    plot.addBand(6.60,1e-1,7.00,hDATA.GetMaximum()/100,r.kRed,0.2)
    plot.addLine(22.0,1e-1,22.0,hDATA.GetMaximum()/100,r.kRed,thickness=1)
    plot.addBand(21.5,1e-1,22.5,hDATA.GetMaximum()/100,r.kRed,0.2)
    
    ### Add background:
    plot.addHisto(hDATA, 'P, SAME', 'Data', 'p', r.kBlack, 1, 1)

    ### Channel banner:
    if LLlabel == 'EE':
        plot.addLatex(0.17, 0.76, 'e^{+}e^{-} channel', font = 42)
    if LLlabel == 'MM':
        plot.addLatex(0.17, 0.76, '#mu^{+}#mu^{-} channel', font = 42)

    ### Extralabel:
    if leftlabel:
        plot.addLatex(0.17, 0.76, leftlabel, font = 42)

    ### Save it
    if not outdir:
        outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/2DPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir)


################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":
    os.system('cat include/koopa_cut.txt')
    '''
    print bcolors.HEADER
    print "                                 ```.....`                              "        
    print "                               `.````.-:::-.-------.                    "        
    print "                              ``      .--.`````.--://-                  "
    print "                                   `   `        `-::://`                "
    print "                                `os/y:           .-::/+/                "
    print "                                yNdod+        ````-:://+`               "
    print "                 ``-/-..------.-hNMMs`     `oy/sy.-:://+.               "
    print "              `..------:+ss/--::/+sh.     `dNh/sN/-://+/`               "
    print "           `.-------------------::///.    /NMMMNm::///+-                "
    print "         `.---------------------:::://:`  -dmNNh/:///+-                 "
    print "        `--------.---.----------:::::/+/.``-++/:-://+o+.                "
    print "       `--------..--------------:::://++/---------:://++:               "
    print "       .----------------------::::::/++:----------::://+o:              "
    print "       ---------------------::::::://:---..------:::://++o`             "
    print "      `:-------------------::::::///:----.------::::://++o.             "
    print "      `::--------------::::::://///:-----------:::::///+++.             "
    print "       ::::------:::::::://////////----------:::::////+++/`             "
    print "       ./:::::::::::///////////////-------:::://////++++/-              "
    print "        -///////////++++++++++++++/----:::///////++++++/-               "
    print "         ./+///+++++ooooooooo+ooo+::::::///++++++++++//.                "
    print "          `-+ooosssooooosssssso+/:::///+++++oooo+++//-`                 "
    print "            `-+so/-.:+sssyyysso++///+oooooooooo++//-`                   "
    print "               ``    `.:+osyyyyysssssssssssoo++++oo+:`..```             "
    print "                         `.-:/+++oooyyysssooosso+ooosoyhyo/-.           "
    print "                                  ``+sssssssssyyso++++/:------::        "
    print "                                  -+oooooooooosyhy+:----.-------::      "
    '''
    print '################################################################################'
    print '                          Starting IFCA-LLP analysis...                         '
    print '################################################################################'

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()
    r.gStyle.SetPadRightMargin(0.12)

    ############# Dat file
    filename = 'dat/Samples_cern_UltraLegacy.dat'

    ############# EG data definition
    DoubleEG2016 = []
    DoubleEG2016.append('DoubleEG_Run2016B_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016C_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016D_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016E_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016F_HIPM')
    DoubleEG2016.append('DoubleEG_Run2016F_noHIPM')
    DoubleEG2016.append('DoubleEG_Run2016G_noHIPM')
    DoubleEG2016.append('DoubleEG_Run2016H_noHIPM')
    DoubleEG2017 = []
    DoubleEG2017.append('DoubleEG_Run2017B')
    DoubleEG2017.append('DoubleEG_Run2017C')
    DoubleEG2017.append('DoubleEG_Run2017D')
    DoubleEG2017.append('DoubleEG_Run2017E')
    DoubleEG2017.append('DoubleEG_Run2017F')
    EGamma2018 = []
    EGamma2018.append('EGamma_Run2018A')
    EGamma2018.append('EGamma_Run2018B')
    EGamma2018.append('EGamma_Run2018C')
    EGamma2018.append('EGamma_Run2018D')


    ############# Muon data definition
    DoubleMuon2016 = []
    DoubleMuon2016.append('DoubleMuon_Run2016B_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016C_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016D_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016E_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016F_HIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016F_noHIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016G_noHIPM')
    DoubleMuon2016.append('DoubleMuon_Run2016H_noHIPM')
    DoubleMuon2018 = []
    DoubleMuon2018.append('DoubleMuon_Run2018A')
    DoubleMuon2018.append('DoubleMuon_Run2018B')
    DoubleMuon2018.append('DoubleMuon_Run2018C')
    DoubleMuon2018.append('DoubleMuon_Run2018D')

    ####### Background MC
    Backgrounds_2017 = []
    Backgrounds_2017.append('DYJetsToLL_M-50_2017')
    Backgrounds_2017.append('DYJetsToLL_M-10to50_2017')
    Backgrounds_2017.append('TTTo2L2Nu_2017')
    Backgrounds_2017.append('WJetsToLNu_2017')
    Backgrounds_2017.append('WW_2017')
    Backgrounds_2017.append('WZ_2017')
    Backgrounds_2017.append('ZZ_2017')
    Backgrounds_2018 = []
    Backgrounds_2018.append('DYJetsToLL_M-50_2018')
    Backgrounds_2018.append('DYJetsToLL_M-10to50_2018')
    Backgrounds_2018.append('TTTo2L2Nu_2018')
    Backgrounds_2018.append('WJetsToLNu_2018')
    Backgrounds_2018.append('WW_2018')
    Backgrounds_2018.append('WZ_2018')
    Backgrounds_2018.append('ZZ_2018')

    ############# Luminosity definition
    lumi2016 = 35.9 # fb-1
    lumi2017 = 41.5 # fb-1
    lumi2018 = 59.7 # fb-1


    ############# Galapago Tree definitions
    treeMC_2017     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_2017, 'MC'), name = 'MC', isdata = 0, close = True)
    treeMC_2018     = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds_2018, 'MC'), name = 'MC', isdata = 0, close = True)

    treeDATA_EG2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2017 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2017, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_EG2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, EGamma2018, 'DATA'), name = 'DATA', isdata = 1 )

    treeDATA_Mu2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    treeDATA_Mu2018 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2018, 'DATA'), name = 'DATA', isdata = 1 )

    ########################################################
    ######## Background optimization (Lxy/dxy bins) ########
    ########################################################

    #www = '/eos/user/f/fernance/www/LLP/SignalRegion-optimization/2DPlots_nLL1'
    www = '/eos/user/r/rlopezru/Galapago_Plots/2Dvxvy_datanocuts'
    region = 'BCR'
    #### -> Electron plots
    #makeBackgroundPlot2D(lumi=lumi2016, hname_bkg='hEE'+region+'_vx_vy', zlog=True, treeDATA=treeDATA_EG2016, inputdir=opts.input, xlabel='', outtag='2016', LLlabel='EE', extralabel=region+' region', outdir=www)
    #makeBackgroundPlot2D(lumi=lumi2017, hname_bkg='hEE'+region+'_vx_vy', zlog=True, treeDATA=treeDATA_EG2017, inputdir=opts.input, xlabel='', outtag='2017', LLlabel='EE', extralabel=region+' region', outdir=www)
    #makeBackgroundPlot2D(lumi=lumi2018, hname_bkg='hEE'+region+'_vx_vy', zlog=True, treeDATA=treeDATA_EG2018, inputdir=opts.input, xlabel='', outtag='2018', LLlabel='EE', extralabel=region+' region', outdir=www)
    #makeBackgroundPlot2D(lumi=lumi2016, hname_bkg='hEE'+region+'_phi_eta', zlog=True, treeDATA=treeDATA_EG2016, inputdir=opts.input, xlabel='', outtag='2016', LLlabel='EE', extralabel=region+' region', outdir=www) 
    #makeBackgroundPlot2D(lumi=lumi2017, hname_bkg='hEE'+region+'_phi_eta', zlog=True, treeDATA=treeDATA_EG2017, inputdir=opts.input, xlabel='', outtag='2017', LLlabel='EE', extralabel=region+' region', outdir=www) 
    #makeBackgroundPlot2D(lumi=lumi2018, hname_bkg='hEE'+region+'_phi_eta', zlog=False, treeDATA=treeDATA_EG2018, inputdir=opts.input, xlabel='', outtag='2018', LLlabel='EE', extralabel=region+' region', outdir=www) 
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEE'+region+'_vr_low', ylog=True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEE'+region+'_vr_low', ylog=True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEE'+region+'_vr_low', ylog=True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hEE'+region+'_vr_high', ylog=True, treeDATA = treeDATA_EG2016, inputdir = opts.input, xlabel = '', outtag = '2016', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEE'+region+'_vr_high', ylog=True, treeDATA = treeDATA_EG2017, inputdir = opts.input, xlabel = '', outtag = '2017', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEE'+region+'_vr_high', ylog=True, treeDATA = treeDATA_EG2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'EE')

    #### -> Muon plots
    #makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMM'+region+'_vx_vy', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = region+' region', outdir = www) 
    #makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMM'+region+'_vx_vy', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = region+' region', outdir = www)
    #makeBackgroundPlot2D(lumi = lumi2016, hname_bkg = 'hMM'+region+'_phi_eta', zlog = True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', LLlabel = 'MM', extralabel = region+' region', outdir = www)
    #makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMM'+region+'_phi_eta', zlog = True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = region+' region', outdir = www)
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMM'+region+'_vr', ylog=True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', outdir = www, yshift = 100.0, LLlabel = 'MM')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMM'+region+'_vr', ylog=True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'MM')
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMM'+region+'_vr_low', ylog=True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', outdir = www, yshift = 100.0, LLlabel = 'MM')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMM'+region+'_vr_low', ylog=True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'MM')
    makeDataMCPlot(lumi = lumi2016, hname_DATA = 'hMM'+region+'_vr_high', ylog=True, treeDATA = treeDATA_Mu2016, inputdir = opts.input, xlabel = '', outtag = '2016', outdir = www, yshift = 100.0, LLlabel = 'MM')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMM'+region+'_vr_high', ylog=True, treeDATA = treeDATA_Mu2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'MM')

    '''
    www = '/eos/user/r/rlopezru/Galapago_Plots/2Dvxvy_1/MC'
    region = 'BCR'
    #### -> Electron plots
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEE'+region+'_vx_vy', zlog = True, treeDATA = treeMC_2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEE'+region+'_vx_vy', zlog = True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEE'+region+'_phi_eta', zlog = False, treeDATA = treeMC_2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEE'+region+'_phi_eta', zlog = False, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeMC_2017, inputdir = opts.input, xlabel = '', outtag = '2017', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'EE')
    region = 'SR'
    #### -> Electron plots
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEE'+region+'_vx_vy', zlog = True, treeDATA = treeMC_2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEE'+region+'_vx_vy', zlog = True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2017, hname_bkg = 'hEE'+region+'_phi_eta', zlog = False, treeDATA = treeMC_2017, inputdir = opts.input, xlabel = '', outtag = '2017', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hEE'+region+'_phi_eta', zlog = False, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'EE', extralabel = region+' region', outdir = www) 
    makeDataMCPlot(lumi = lumi2017, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeMC_2017, inputdir = opts.input, xlabel = '', outtag = '2017', outdir = www, yshift = 100.0, LLlabel = 'EE')
    makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hEE'+region+'_vr', ylog=True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'EE')

    #### -> Muon plots
    #makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMM'+region+'_vx_vy', zlog = True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = region+' region', outdir = www)
    #makeBackgroundPlot2D(lumi = lumi2018, hname_bkg = 'hMM'+region+'_phi_eta', zlog = True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', LLlabel = 'MM', extralabel = region+' region', outdir = www)
    #makeDataMCPlot(lumi = lumi2018, hname_DATA = 'hMM'+region+'_vr', ylog=True, treeDATA = treeMC_2018, inputdir = opts.input, xlabel = '', outtag = '2018', outdir = www, yshift = 100.0, LLlabel = 'MM')
    '''
