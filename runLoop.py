import ROOT as r
import argparse
from include.processHandler import *
from include.yieldHandler import *
from include.plotHandler import *

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
parser.add_argument('-m', '--mode', action='store', type=str, dest='mode', help='mode')
parser.add_argument('-c', '--config', action='store', dest='config', help='Configuration file')
parser.add_argument('-y', '--year', action='store', dest='year', help='year')
parser.add_argument('-r', '--raw', action='store_true', dest='raw', help='apply scale factors or not')
parser.add_argument('-d', '--data', action='store_true', dest='data', help='True if DATA')
args = parser.parse_args()

##
## -- Load TFile and TTree
##

tfile = r.TFile(args.file)
ttree = tfile.Get("Events")

print('data?', args.data)

if args.data:
    if args.mode == 'plot':
        process = plotHandler(args.out, args.name, args.block, args.sample, args.tree, 1, True, args.config, args.year, args.raw)
    elif args.mode == 'yield':
        process = yieldHandler(args.out, args.name, args.block, args.sample, args.tree, 1, True, args.config, args.year, args.raw)
    elif args.mode == 'eff':
        process = effHandler(args.out, args.name, args.block, args.sample, args.tree, 1, True, args.config, args.year, args.raw)
    else:
        process = processHandler(args.out, args.name, args.block, args.sample, args.tree, 1, True, args.config, args.year, args.raw)
else:
    if args.mode == 'plot':
        process = plotHandler(args.out, args.name, args.block, args.sample, args.tree, args.lumiweight, False, args.config, args.year, args.raw)
    elif args.mode == 'yield':
        process = yieldHandler(args.out, args.name, args.block, args.sample, args.tree, args.lumiweight, False, args.config, args.year, args.raw)
    elif args.mode == 'eff':
        process = effHandler(args.out, args.name, args.block, args.sample, args.tree, args.lumiweight, False, args.config, args.year, args.raw)
    else:
        process = processHandler(args.out, args.name, args.block, args.sample, args.tree, args.lumiweight, False, args.config, args.year, args.raw)

##
## -- Loop over the events
##

for n,ev in enumerate(ttree):
    process.processEvent(ev)


process.Write()
