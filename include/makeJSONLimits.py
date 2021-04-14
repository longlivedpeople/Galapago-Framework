import os
import ROOT as r
import optparse
import json


if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Target directory')
    parser.add_option('-n', '--name', action='store', type=str, dest='name', default='', help='Target directory')
    parser.add_option('-x', '--xsec', action='store', type=float, dest='xsec', default=1.0, help='Target directory')
    parser.add_option('-o', '--outdir', action='store', type=str, dest='outdir', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    targets = opts.input.split(',')
    # e.g. higgsCombine_mH1000_mS150.AsymptoticLimits.mH1.root
    xsec = opts.xsec

    output_dict_mu = {}
    output_dict_xsec = {}

    for target in targets:

        mass = target.split('.')[-2][2:]
        # convert from mm to cm
        if mass == '1':
            mass = '0.1'
        else:
            mass = mass[:-1]  
        print(mass)

        _file = r.TFile(target)
        _limits = _file.Get("limit")

        output_dict_mu[mass] = {}
        output_dict_xsec[mass] = {}
        _limits.GetEntry(0)
        output_dict_mu[mass]["exp-2"] = _limits.limit
        output_dict_xsec[mass]["exp-2"] = xsec * _limits.limit
        _limits.GetEntry(1)
        output_dict_mu[mass]["exp-1"] = _limits.limit
        output_dict_xsec[mass]["exp-1"] = xsec * _limits.limit
        _limits.GetEntry(2)
        output_dict_mu[mass]["exp0"] = _limits.limit
        output_dict_xsec[mass]["exp0"] = xsec * _limits.limit
        _limits.GetEntry(3)
        output_dict_mu[mass]["exp+1"] = _limits.limit
        output_dict_xsec[mass]["exp+1"] = xsec * _limits.limit
        _limits.GetEntry(4)
        output_dict_mu[mass]["exp+2"] = _limits.limit
        output_dict_mu[mass]["obs"] = output_dict_mu[mass]["exp0"]
        output_dict_xsec[mass]["exp+2"] = xsec * _limits.limit
        output_dict_xsec[mass]["obs"] = output_dict_xsec[mass]["exp0"]

    outdir = opts.outdir
    if outdir[-1] != '/': outdir = outdir + '/'

    with open(opts.outdir + opts.name + '_mu.json', 'w') as outfile:
        json.dump(output_dict_mu, outfile)

    with open(opts.outdir + opts.name + '_xsec.json', 'w') as outfile:
        json.dump(output_dict_xsec, outfile)
