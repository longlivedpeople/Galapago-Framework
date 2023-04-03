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

def createDatacards(datacards, Systematics, Backgrounds, Signals, flavor, year, exclude = []):

    treeBKG = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Backgrounds, 'DATA'), name = 'DATA', isdata = 1 )
    datacard_list = []

    for datacard_key in datacards.keys():

        datacard = datacards[datacard_key]
        channels = {}

        ### Initialize channels with background information
        for channel_name in datacard.keys():

            if channel_name in exclude:
                continue

            xmin = datacard[channel_name]['bkg']['limits'][0] # minimum value
            xmax = datacard[channel_name]['bkg']['limits'][1] # minimum value
            directory = datacard[channel_name]['bkg']['dir']
            histogram = datacard[channel_name]['bkg']['histogram']

            print('>>>>> Info channel', channel_name, xmin, xmax, directory, histogram)

            channels[channel_name] = DatacardManager.Channel(channel_name, xmin, xmax)

            ## Get histogram
            _histo = treeBKG.getLoopTH1F(directory, histogram)
            print(directory, histogram)
            c1 = r.TCanvas("", "", 500, 500)
            c1.SetLogy(1)
            _histo.Draw('HIST')
            c1.Print('Pruebafondochannel.png')
            channels[channel_name].addBackground('bkg', _histo)
            print(channels[channel_name].backgrounds[0].rate)
         
        for sample in Signals:

            samplelist = []
            if year=='2016' and flavor!='Electron':
                samplelist = [sample, sample + 'APV']
            else:
                samplelist = [sample]

            treeSI = Sample.Tree(fileName = helper.selectSamples(WORKPATH + 'dat/CombSignal_'+year+'UL_Fall22.dat', samplelist, 'SI'), name = 'SI', isdata = 0)

            sample_split = sample.split('_')
            mH = sample_split[1]
            mS = sample_split[2]
            ctau = sample_split[3]
            sample_label = 'mH' + mH + '__' + 'mS' + mS + '__' + 'ctau' + ctau + 'mm'

            ## Create datacard:
            datacard_name = 'Datacard__' + flavor + '__' + sample_label + '__' + year  + '.txt'
            output_datacard = DatacardManager.Datacard(datacard_name, year, Systematics)

            for channel_name in datacard.keys():
                if channel_name in exclude:
                    continue
                histogram = datacard[channel_name]['sig']['histogram']
                directory = datacard[channel_name]['sig']['dir']
                _histo = treeSI.getLoopTH1F(directory, histogram)
                channels[channel_name].setSignal('sig', _histo)
                output_datacard.addChannel(channels[channel_name])
            
            output_datacard.saveDatacard(outputDir = _outdir)
            datacard_list.append(datacard_name)

    return datacard_list


def executeCombine(_d, year):

    name = _d.replace('Datacard__', '')
    name = name.replace('.txt', '')
    name = '__' + name.split('__ctau')[0] + '__era' + year
    ctau = (_d.split('__ctau')[1]).split('mm')[0]
    
    ### Move to the working dir:
    command = 'combine -M AsymptoticLimits -n {0} -m {1} --run blind {2}'.format(name, ctau, _d)
    os.system(command)
    



if __name__ == "__main__":


    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-d', '--dat', action='store', type=str, dest='dat', default='dat/Samples_cern_UltraLegacy.dat', help='dat file')
    parser.add_option('-t', '--t', action='store', type=str, dest='tag', default='', help='tag')
    parser.add_option('-e', '--electronRecipe', action='store', type=str, dest='electronRecipe', default='', help='the input dir')
    parser.add_option('-m', '--muonRecipe',   action='store', type=str, dest='muonRecipe', default='', help='the input dir')
    parser.add_option('-s', '--Esystematics', action='store', type=str, dest='electron_systematics', default='recipes-datacards/recipe_Systematics_UltraLegacy_Electron.txt', help='file with systematics table')
    parser.add_option('-S', '--Msystematics', action='store', type=str, dest='muon_systematics', default='recipes-datacards/recipe_Systematics_UltraLegacy_Muon.txt', help='file with systematics table')
    parser.add_option('-T', '--theory', action='store', type=str, dest='theory', default='', help='theory')

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
    Masses = []
    """
    Masses.append('HSS_125_50')
    Masses.append('HSS_300_50')
    Masses.append('HSS_500_50')
    Masses.append('HSS_500_150')
    """
    Masses.append('HSS_600_50')
    Masses.append('HSS_600_150')
    """
    Masses.append('HSS_600_250')
    Masses.append('HSS_800_50')
    Masses.append('HSS_800_250')
    Masses.append('HSS_800_350')
    """
    Masses.append('HSS_1000_250')
    Masses.append('HSS_1000_350')
    Masses.append('HSS_1000_450')
    #Masses.append('RPV_350_148')
    #Masses.append('RPV_1500_494')
    Signals = []
    for mass in Masses:
        Signals.append(mass + '_1')
        Signals.append(mass + '_10')
        Signals.append(mass + '_100')
        Signals.append(mass + '_1000')
        Signals.append(mass + '_10000')
    Signals2016 = [i + '_2016' for i in Signals]
    Signals2017 = [i + '_2017' for i in Signals]
    Signals2018 = [i + '_2018' for i in Signals]



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

    ############# Electron data definition
    DoubleEG2016 = []
    #DoubleEG2016.append('DoubleEG_Run2016B_HIPM')
    #DoubleEG2016.append('DoubleEG_Run2016C_HIPM')
    #DoubleEG2016.append('DoubleEG_Run2016D_HIPM')
    #DoubleEG2016.append('DoubleEG_Run2016E_HIPM')
    #DoubleEG2016.append('DoubleEG_Run2016F_HIPM')
    #DoubleEG2016.append('DoubleEG_Run2016F_noHIPM')
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


    ########### .dat definition
    filename = opts.dat


    ########### Data Tree's declaration
    #treeMuonDATA_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleMuon2016, 'DATA'), name = 'DATA', isdata = 1 )
    #treeEGDATA_2016 = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, DoubleEG2016, 'DATA'), name = 'DATA', isdata = 1 )


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
    owd = os.getcwd()

    muon_datacard_names_2016 = []
    muon_datacard_names_2018 = []
    electron_datacard_names_2016 = []
    electron_datacard_names_2017 = []
    electron_datacard_names_2018 = []

    if doMuons:
        muon_datacard_names_2016 = createDatacards(datacards = muon_datacards, Systematics = opts.muon_systematics, Backgrounds = DoubleMuon2016, Signals = Signals2016, flavor = 'Muon', year = '2016')
        muon_datacard_names_2018 = createDatacards(datacards = muon_datacards, Systematics = opts.muon_systematics, Backgrounds = DoubleMuon2018, Signals = Signals2018, flavor = 'Muon', year = '2018')
    if doElectrons:
        electron_datacard_names_2016 = createDatacards(datacards = electron_datacards, Systematics = opts.electron_systematics, Backgrounds = DoubleEG2016, Signals = Signals2016, flavor = 'Electron', year = '2016')
        electron_datacard_names_2017 = createDatacards(datacards = electron_datacards, Systematics = opts.electron_systematics, Backgrounds = DoubleEG2017, Signals = Signals2017, flavor = 'Electron', year = '2017', exclude = ['nEE_IaA', 'nEE_IaB', 'nEE_IaC'])
        electron_datacard_names_2018 = createDatacards(datacards = electron_datacards, Systematics = opts.electron_systematics, Backgrounds = EGamma2018, Signals = Signals2018, flavor = 'Electron', year = '2018')

    os.chdir(_outdir) 

    # Create combined Muon + Electron datacard
    joint_datacard_names_2016 = []
    joint_datacard_names_2018 = []
    if len(muon_datacard_names_2016) == len(electron_datacard_names_2016):
        muon_datacard_names_2016.sort()
        electron_datacard_names_2016.sort()
        for i in range(0, len(muon_datacard_names_2016)):
            sufix = muon_datacard_names_2016[i].replace('Datacard__Muon__', '')
            command = 'combineCards.py ' + muon_datacard_names_2016[i] + ' ' + electron_datacard_names_2016[i] + ' > Datacard__Joint__' + sufix
            os.system(command)
            joint_datacard_names_2016.append('Datacard__Joint__' + sufix)
    if len(muon_datacard_names_2018) == len(electron_datacard_names_2018):
        muon_datacard_names_2018.sort()
        electron_datacard_names_2018.sort()
        for i in range(0, len(muon_datacard_names_2018)):
            sufix = muon_datacard_names_2018[i].replace('Datacard__Muon__', '')
            command = 'combineCards.py ' + muon_datacard_names_2018[i] + ' ' + electron_datacard_names_2018[i] + ' > Datacard__Joint__' + sufix
            os.system(command)
            joint_datacard_names_2018.append('Datacard__Joint__' + sufix)




    # create combined years datacard
    muon_datacard_names_full = []
    if len(muon_datacard_names_2016) == len(muon_datacard_names_2018):
        muon_datacard_names_2016.sort()
        muon_datacard_names_2018.sort()
        for i in range(0, len(muon_datacard_names_2016)):
            out_name = muon_datacard_names_2016[i].replace('2016', 'Full')
            command = 'combineCards.py ' + muon_datacard_names_2016[i] + ' ' + muon_datacard_names_2018[i] + ' > ' + out_name
            os.system(command)
            muon_datacard_names_full.append(out_name)

    electron_datacard_names_full = []
    if len(electron_datacard_names_2016) == len(electron_datacard_names_2018) and len(electron_datacard_names_2016) == len(electron_datacard_names_2017):
        electron_datacard_names_2016.sort()
        electron_datacard_names_2017.sort()
        electron_datacard_names_2018.sort()
        for i in range(0, len(electron_datacard_names_2016)):
            out_name = electron_datacard_names_2016[i].replace('2016', 'Full')
            command = 'combineCards.py ' + electron_datacard_names_2016[i] + ' ' + electron_datacard_names_2017[i] + ' ' + electron_datacard_names_2018[i] + ' > ' + out_name
            os.system(command)
            electron_datacard_names_full.append(out_name)

    joint_datacard_names_full = []
    if len(electron_datacard_names_full) == len(muon_datacard_names_full):
        electron_datacard_names_full.sort()
        muon_datacard_names_full.sort()
        for i in range(0, len(electron_datacard_names_full)):
            out_name = electron_datacard_names_full[i].replace('Electron', 'Joint')
            command = 'combineCards.py ' + muon_datacard_names_full[i] + ' ' + electron_datacard_names_full[i] + ' > ' + out_name
            os.system(command)
            joint_datacard_names_full.append(out_name)

    os.chdir(owd)



    ###########################################
    ########## Extract limits using combine
    #####
    owd = os.getcwd()
    os.chdir(_outdir) 
    print("> Muon datacards to process: ")
    if doMuons:
        for _d in muon_datacard_names_2016:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2016')
        for _d in muon_datacard_names_2018:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2018')
    print("> Electrons datacards to process: ")
    if doElectrons:
        for _d in electron_datacard_names_2016:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2016')
        for _d in electron_datacard_names_2017:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2017')
        for _d in electron_datacard_names_2018:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2018')
    print("> Joint datacards to process: ")
    if doElectrons and doMuons and len(joint_datacard_names_2016) > 0:
        for _d in joint_datacard_names_2016:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2016')
        for _d in joint_datacard_names_2018:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, '2018')
    if len(muon_datacard_names_full):
        for _d in muon_datacard_names_full:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, 'Full')
    if len(electron_datacard_names_full):
        for _d in electron_datacard_names_full:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, 'Full')
    if len(joint_datacard_names_full):
        for _d in joint_datacard_names_full:
            print('  ' +  _d + ' limits:' )
            executeCombine(_d, 'Full')
    os.chdir(owd)

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

    print(grouped)

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

    print("----->>>> grouped.keys():")
    print(grouped.keys())

    sets = {}
    for point_id in grouped.keys():
        mH     = point_id.split('mH')[1].split('__')[0]
        mS     = point_id.split('mS')[1].split('__')[0]
        year   = point_id.split('era')[1]
        flavor = point_id.split('__')[1]
        if mH not in sets.keys():
            sets[mH] = {}
            sets[mH] = {}
            sets[mH]['2016'] = {}
            sets[mH]['2017'] = {}
            sets[mH]['2018'] = {}
            sets[mH]['Full'] = {}
            sets[mH]['2016']['Electron'] = []
            sets[mH]['2016']['Muon']     = []
            sets[mH]['2016']['Joint']    = []
            sets[mH]['2017']['Electron'] = []
            sets[mH]['2018']['Electron'] = []
            sets[mH]['2018']['Muon']     = []
            sets[mH]['2018']['Joint']    = []
            sets[mH]['Full']['Electron'] = []
            sets[mH]['Full']['Muon']     = []
            sets[mH]['Full']['Joint']    = []

        sets[mH][year][flavor].append(grouped[point_id]['json']) 

    for mH in sets.keys():
        for year in sets[mH].keys():
            for flavor in sets[mH][year].keys():

                if len(sets[mH][year][flavor]) < 1: continue

                plot_input = ''
                for json in sets[mH][year][flavor]:
                    plot_input = plot_input + json + ','
                plot_input = plot_input[:-1]

                command = 'python include/plotLimits.py -j {0} -m {1} -f {2} -o {3} -y {4} -t {5}'.format(plot_input, mH, flavor, _outdir, year, opts.theory)
                print(command)
                os.system(command)






