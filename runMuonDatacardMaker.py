##################################################################################################
#                           ____       _                                                         #
#                          / ___| __ _| | __ _ _ __   __ _  __ _  ___                            #
#                         | |  _ / _` | |/ _` | '_ \ / _` |/ _` |/ _ \                           #
#                         | |_| | (_| | | (_| | |_) | (_| | (_| | (_) |                          # 
#                     _____\____|\__,_|_|\__,_| .__/ \__,_|\__, |\___/ _                         #  
#                    |  ___| __ __ _ _ __ ___ |_|____      |___/  _ __| | _                      #_
#                    | |_ | '__/ _` | '_ ` _ \ / _ \ \ /\ / / _ \| '__| |/ /                     #
#                    |  _|| | | (_| | | | | | |  __/\ V  V / (_) | |  |   <                      # 
#                    |_|  |_|  \__,_|_| |_| |_|\___| \_/\_/ \___/|_|  |_|\_\                     #
#                                                                                                # 
##################################################################################################

import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect
import numpy as np

import include.Sample as Sample
import include.helper as helper
import include.Canvas as Canvas
import include.DatacardManager as DatacardManager

####################################### CLASS DEFINITION #########################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'




################################# GLOBAL VARIABLES DEFINITION ####################################

WORKPATH = os.path.abspath('./') + '/'

##################################### FUNCTION DEFINITION ########################################

 

if __name__ == "__main__":

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
    print '########################################################################' 
    print '                  Starting IFCA-LLP analysis...                         '
    print '########################################################################' + bcolors.ENDC

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-d', '--dat', action='store', type=str, dest='dat', default='dat/Samples_cern_Legacy.dat', help='dat file')
    parser.add_option('-r', '--recipe', action='store', type=str, dest='recipe', default='', help='the input dir')

    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    gROOT.ProcessLine('.L ' + WORKPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()


    ###########################################
    ########## Tree initialization
    #####

    ############# Signal definition
    Signals = []
    Signals.append('HXX_400_50_1mm')
    Signals.append('HXX_400_50_10mm')
    Signals.append('HXX_400_50_100mm')
    Signals.append('HXX_400_50_1000mm')
    Signals.append('HXX_300_20_1mm')
    Signals.append('HXX_300_20_100mm')
    Signals.append('HXX_300_20_1000mm')
    Signals.append('HXX_300_20_10000mm')
    Signals.append('HXX_300_50_1mm')
    Signals.append('HXX_300_50_10mm')
    Signals.append('HXX_300_50_100mm')
    Signals.append('HXX_300_50_1000mm')
    Signals.append('HXX_300_50_10000mm')
    Signals.append('HXX_300_150_1mm')
    Signals.append('HXX_300_150_10mm')
    Signals.append('HXX_300_150_100mm')
    Signals.append('HXX_300_150_1000mm')
    Signals.append('HXX_300_150_10000mm')

    ############# Muon data definition
    DoubleMuonB = 'DoubleMuon_Run2016B'
    DoubleMuonC = 'DoubleMuon_Run2016C'
    DoubleMuonD = 'DoubleMuon_Run2016D'
    DoubleMuonE = 'DoubleMuon_Run2016E'
    DoubleMuonF = 'DoubleMuon_Run2016F'
    DoubleMuonG = 'DoubleMuon_Run2016G'
    DoubleMuonH = 'DoubleMuon_Run2016H'

    DoubleMuon_list = []
    DoubleMuon_list.append(DoubleMuonB)
    DoubleMuon_list.append(DoubleMuonC)
    DoubleMuon_list.append(DoubleMuonD)
    DoubleMuon_list.append(DoubleMuonE)
    DoubleMuon_list.append(DoubleMuonF)
    DoubleMuon_list.append(DoubleMuonG)
    DoubleMuon_list.append(DoubleMuonH)


    ########### .dat definition
    filename = opts.dat


    ########### Data Tree's declaration
    treeMuonDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )



    ###########################################
    ########## Read recipe
    #####

    recipe = DatacardManager.Recipe(opts.recipe)
    datacards = recipe.dcs


    ###########################################
    ########## Loop over datacards
    #####

    for datacard_name in datacards.keys():

        datacard = datacards[datacard_name]
        channels = {}


        ### Initialize channels with background information
        for channel_name in datacard.keys():

            xmin = datacard[channel_name]['bkg']['limits'][0] # minimum value
            xmax = datacard[channel_name]['bkg']['limits'][1] # minimum value
            directory = datacard[channel_name]['bkg']['dir']
            histogram = datacard[channel_name]['bkg']['histogram']

            channels[channel_name] = DatacardManager.Channel(channel_name, xmin, xmax)

            ## Get histogram
            _histo = treeMuonDATA.getLoopTH1F(directory, histogram)
            channels[channel_name].addBackground('bkg', _histo)

         
        for sample in Signals:

            treeSI = Sample.Tree(fileName = helper.selectSamples(WORKPATH + filename, [sample], 'SI'), name = 'SI', isdata = 0)

            ## Create datacard:
            output_datacard = DatacardManager.Datacard('Datacard_' + sample + '.txt')

            for channel_name in datacard.keys():
                histogram = datacard[channel_name]['sig']['histogram']
                directory = datacard[channel_name]['sig']['dir']
                _histo = treeSI.getLoopTH1F(directory, histogram)
                channels[channel_name].setSignal('sig', _histo)
                output_datacard.addChannel(channels[channel_name])
            
            output_datacard.saveDatacard(outputDir = 'Datacards/')




