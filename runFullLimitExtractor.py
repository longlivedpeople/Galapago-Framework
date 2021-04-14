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

def createDatacards(datacards, Backgrounds, Signals, flavor): 

    treeBKG = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'DATA'), name = 'DATA', isdata = 1 )
    datacard_list = []

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
            _histo = treeBKG.getLoopTH1F(directory, histogram)
            channels[channel_name].addBackground('bkg', _histo)

         
        for sample in Signals:

            treeSI = Sample.Tree(fileName = helper.selectSamples(WORKPATH + filename, [sample], 'SI'), name = 'SI', isdata = 0)

            ## Create datacard:
            datacard_name = 'Datacard' + flavor + '_' + sample + '.txt'
            output_datacard = DatacardManager.Datacard(datacard_name)

            for channel_name in datacard.keys():
                histogram = datacard[channel_name]['sig']['histogram']
                directory = datacard[channel_name]['sig']['dir']
                _histo = treeSI.getLoopTH1F(directory, histogram)
                channels[channel_name].setSignal('sig', _histo)
                output_datacard.addChannel(channels[channel_name])
            
            output_datacard.saveDatacard(outputDir = _outdir)
            datacard_list.append(datacard_name)

    return datacard_list


def executeCombine(_d):

    name = _d.replace('Datacard', '')
    name = name.replace('.txt', '')
    name = name.replace(name.split('_')[-1], '')
    if name[-1] == '_': name = name[:-1]
    #print(name)

    ctau = (_d.split('mm')[0]).split('_')[-1]
    #print(ctau)
    
    ### Move to the working dir:
    command = 'combine -M AsymptoticLimits -n {0} -m {1} {2}'.format(name, ctau, _d)
    os.system(command)
    



if __name__ == "__main__":


    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-d', '--dat', action='store', type=str, dest='dat', default='dat/Samples_cern_Legacy.dat', help='dat file')
    parser.add_option('-t', '--t', action='store', type=str, dest='tag', default='', help='tag')
    parser.add_option('-e', '--electronRecipe', action='store', type=str, dest='electronRecipe', default='', help='the input dir')
    parser.add_option('-m', '--muonRecipe', action='store', type=str, dest='muonRecipe', default='', help='the input dir')

    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    gROOT.ProcessLine('.L ' + WORKPATH + 'include/tdrstyle.C')
    gROOT.SetBatch(1)
    r.setTDRStyle()

    ###########################################
    ########## Output Directory
    #####
    global _outdir 
    _outdir = WORKPATH + 'LimitsResults'
    if opts.tag != '': _outdir = _outdir + '_' + opts.tag
    if _outdir[-1] != '/': _outdir = _outdir + '/'


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


    ############# EG data definition
    DoubleEGB = 'DoubleEG_Run2016B'
    DoubleEGC = 'DoubleEG_Run2016C'
    DoubleEGD = 'DoubleEG_Run2016D'
    DoubleEGE = 'DoubleEG_Run2016E'
    DoubleEGF = 'DoubleEG_Run2016F'
    DoubleEGG = 'DoubleEG_Run2016G'
    DoubleEGH = 'DoubleEG_Run2016H'

    DoubleEG_list = []
    DoubleEG_list.append(DoubleEGB)
    DoubleEG_list.append(DoubleEGC)
    DoubleEG_list.append(DoubleEGD)
    DoubleEG_list.append(DoubleEGE)
    DoubleEG_list.append(DoubleEGF)
    DoubleEG_list.append(DoubleEGG)
    DoubleEG_list.append(DoubleEGH)


    ########### .dat definition
    filename = opts.dat


    ########### Data Tree's declaration
    treeMuonDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon_list, 'DATA'), name = 'DATA', isdata = 1 )
    treeEGDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG_list, 'DATA'), name = 'DATA', isdata = 1 )


    ###########################################
    ########## Read recipe
    #####

    doMuons = False
    doElectrons = False

    if opts.muonRecipe != '':
        doMuons = True
        muon_recipe = DatacardManager.Recipe(opts.muonRecipe)
        muon_datacards = muon_recipe.dcs
    if opts.electronRecipe != '':
        doElectrons = True
        electron_recipe = DatacardManager.Recipe(opts.electronRecipe)
        electron_datacards = electron_recipe.dcs


    ###########################################
    ########## Loop over datacards
    #####

    muon_datacard_names = []
    electron_datacard_names = []

    if doMuons:
        muon_datacard_names = createDatacards(datacards = muon_datacards, Backgrounds = DoubleMuon_list, Signals = Signals, flavor = 'Muon')
    if doElectrons:
        electron_datacard_names = createDatacards(datacards = electron_datacards, Backgrounds = DoubleEG_list, Signals = Signals, flavor = 'Electron')


    ###########################################
    ########## Extract limits using combine
    #####
    """
    owd = os.getcwd()
    os.chdir(_outdir) 
    print("> Muon datacards to process: ")
    for _d in muon_datacard_names:
        print('  ' +  _d + ' limits:' )

        executeCombine(_d)

    os.chdir(owd)
    """

    ###########################################
    ########## Make json
    #####

    # Collect combine output:
    grouped = {}
    for _f in os.listdir(_outdir):
        if 'higgsCombine' not in _f or '.root' not in _f or 'AsymptoticLimits' not in _f: continue
        point_id = _f.split('.')[0]

        if point_id not in grouped.keys():
            grouped[point_id] = {}
            grouped[point_id]['combine'] = []
            grouped[point_id]['json'] = ''

        grouped[point_id]['combine'].append(_outdir + _f)

    # Make json for each group:
    for point_id in grouped.keys():
        json_input = ''
        for point in grouped[point_id]['combine']: 
            json_input = json_input + point + ','
        json_input = json_input[:-1] # remove the last comma
        
        json_name = 'JSON' + point_id
        command = 'python include/makeJSONLimits.py -i {0} -n {1} -x 0.001 -o {2}'.format(json_input, json_name, _outdir)
        os.system(command)

        json_name_xsec = _outdir + json_name + '_xsec.json'

        grouped[point_id]['json'] = json_name_xsec


    ###########################################
    ########## Plot Limits
    #####

    sets = {}
    for point_id in grouped.keys():
        mH = point_id.split('_')[2]
       
        if mH not in sets.keys():
            sets[mH] = {}
            sets[mH]['jsons'] = []
            sets[mH]['flavor'] = ''
            if 'Muon' in point_id:
                sets[mH]['flavor'] = 'Muon'
            if 'Electron' in point_id:
                sets[mH]['flavor'] = 'Electron'

        sets[mH]['jsons'].append(grouped[point_id]['json'])

    for mH in sets.keys():
        plot_input = ''
        for json in sets[mH]['jsons']:
            plot_input = plot_input + json + ','
        plot_input = plot_input[:-1]

        command = 'python include/plotLimits.py -j {0} -m {1} -f {2} -o {3}'.format(plot_input, mH, sets[mH]['flavor'], _outdir)
        os.system(command)






