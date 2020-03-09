import os
import ROOT as r
from ROOT import TFile
import argparse
import copy

if __name__=="__main__":


    ## parser object
    parser = argparse.ArgumentParser(description='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_argument('-i', '--inputDir', action='store', type=str, required=True, dest='inputDir', help='Directory that contains all the input histograms (abs. path)')
    parser.add_argument('-o', '--outputFile', action='store', type=str, required=True, dest='outputFile', help='Merged file (abs.path)')

    args = parser.parse_args()

    hList = []
    for rootfile in os.listdir(args.inputDir):
        _file = TFile('outputHistograms/'+rootfile)
        for k in _file.GetListOfKeys():
            h = k.ReadObj()
            print('> merging ' + h.GetName() + '...')
            hList.append(copy.deepcopy(h))
    _outFile = TFile(args.outputFile, 'RECREATE')
    _outFile.cd()
    for h in hList:
        h.Write()
    _outFile.Close()
