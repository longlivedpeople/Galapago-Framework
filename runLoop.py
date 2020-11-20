import ROOT as r
import optparse
from include.processHandler import *

##
## -- Parser object
##

parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
parser.add_option('-f', '--file', action='store', type=str, dest='file', default='', help='Path to TFile')
parser.add_option('-o', '--out', action='store', type=str, dest='out', default='', help='Path to outputdir')
parser.add_option('-n', '--name', action='store', type=str, dest='name', default='', help='Name of the Tree: MC, SI or DATA')
parser.add_option('-b', '--block', action='store', type=str, dest='block', default='', help='Block name')
parser.add_option('-s', '--sample', action='store', type=str, dest='sample', default='', help='Sample name')
parser.add_option('-t', '--tree', action='store', type=str, dest='tree', default='', help='TFile/TTree number')
parser.add_option('-l', '--lumiweight', action='store', type=float, dest='lumiweight', help='lumiweight')
parser.add_option('-d', '--data', action='store_true', dest='data', help='True if DATA')
(opts, args) = parser.parse_args()

##
## -- Load TFile and TTree
##

tfile = r.TFile(opts.file)
ttree = tfile.Get("Events")

if opts.data:
    process = processHandler(opts.out, opts.name, opts.block, opts.sample, opts.tree, 1, True)
else:
    process = processHandler(opts.out, opts.name, opts.block, opts.sample, opts.tree, opts.lumiweight, False)

##
## -- Loop over the events
##

for n,ev in enumerate(ttree):

    process.countLLs(ev)


process.Write()
