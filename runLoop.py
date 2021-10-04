import ROOT as r
import argparse
from include.processHandler import *

##
## -- Parser object
##

parser = argparse.ArgumentParser(usage='usage: %prog [args] FilenameWithSamples', version='%prog 1.0')
parser.add_argument('-f', '--file', action='store', type=str, dest='file', default='', help='Path to TFile')
parser.add_argument('-o', '--out', action='store', type=str, dest='out', default='', help='Path to outputdir')
parser.add_argument('-n', '--name', action='store', type=str, dest='name', default='', help='Name of the Tree: MC, SI or DATA')
parser.add_argument('-b', '--block', action='store', type=str, dest='block', default='', help='Block name')
parser.add_argument('-s', '--sample', action='store', type=str, dest='sample', default='', help='Sample name')
parser.add_argument('-t', '--tree', action='store', type=str, dest='tree', default='', help='TFile/TTree number')
parser.add_argument('-l', '--lumiweight', action='store', type=float, dest='lumiweight', help='lumiweight')
parser.add_argument('-y', '--year', action='store', dest='year', help='year')
parser.add_argument('-d', '--data', action='store_true', dest='data', help='True if DATA')
parser.add_argument('--doEffs', action='store_true', dest='doEffs', help='True if doing only efficiencies')
args = parser.parse_args()

##
## -- Load TFile and TTree
##

tfile = r.TFile(args.file)
ttree = tfile.Get("Events")

if args.data:
    process = processHandler(args.out, args.name, args.block, args.sample, args.tree, 1, True, args.year)
else:
    process = processHandler(args.out, args.name, args.block, args.sample, args.tree, args.lumiweight, False, args.year)

##
## -- Loop over the events
##

for n,ev in enumerate(ttree):

    if not args.doEffs:
        process.processEvent(ev)
    else:
        process.countLLs(ev)


process.Write()
