import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect

import include.Sample as Sample
import include.helper as helper
import include.Canvas as Canvas



############ Enter the splitting parameters using a parser
parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
parser.add_option('-s', '--sample', action='store', type=str, dest='SAMPLE', help='The sample of the splitting')
parser.add_option('-d', '--dat', action='store', type=str, dest='DATFILE', help='The sample datfile that is used in this split')
parser.add_option('-o', '--out', action='store', type=str, dest='OUTPATH', help='The path where the output histogram will be stored')
parser.add_option('-l', '--lumi', action='store', type=float, dest='LUMI', help='Luminosity value')
parser.add_option('-p', '--option', action='store', type=str, dest='OPTION', help='Option')

(opts, args) = parser.parse_args()

############ datFile definition
datFile = opts.DATFILE


############ Initialize the output root file
outputPath = opts.OUTPATH if opts.OUTPATH[-1] == '/' else opts.OUTPATH + '/'
if not os.path.exists(outputPath): os.makedirs(outputPath) # ensure the path
outputName = outputPath + 'fill_' + opts.SAMPLE + '.root'
output = r.TFile( outputName, "RECREATE" )
output.Close()

sampleList = []
sampleList.append(opts.SAMPLE)

############ Tree creation
tree = Sample.Tree(helper.selectSamples(opts.DATFILE, sampleList, opts.OPTION), opts.OPTION, 0, outputName)

############ Tree loop
tree.Loop(opts.LUMI, False, False)
