import ROOT as r
import optparse
import os
import re

##
## -- Parser object
##

parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
parser.add_option('-d', '--dir', action='store', type=str, dest='dir', default='', help='Path to TFile')
parser.add_option('-s', '--submit', action='store_true', dest='submit', default='', help='Submit?')
(opts, args) = parser.parse_args()

##
## -- Save all the files in dir
##

logs  = []
outs  = []
errs  = []
roots = []
others = []

for _f in os.listdir(opts.dir):

    if _f[-4:] == '.log': logs.append(_f)
    elif _f[-4:] == '.err': errs.append(_f)
    elif _f[-4:] == '.out': outs.append(_f)
    elif _f[-5:] == '.root': roots.append(_f.replace('__', '_'))
    else: others.append(_f)

print('>>> Found files:')
print('    .log files: ' + str(len(logs)) )
print('    .err files: ' + str(len(errs)) )
print('    .out files: ' + str(len(outs)) )
print('    .root files: ' + str(len(roots)) )


##
## -- Explore which root files were not produced
##

missedroots = []
missederrs = []

print('>>> Missed root files:')
for log in logs:

    name = log[:-4]
    errname = name + '.err'
    outname = name + '.out'
    rootname = name+ '.root'

    if rootname not in roots: 
        missedroots.append(rootname)
    if errname not in errs:
        missederrs.append(errname)
    #    print(log)

print('>>> Missed root files:')
for rname in missedroots:
    print(rname)

print('>>> Missed err files:')
missederrs_noerr = []
for rname in missederrs:
    print(rname)
    missederrs_noerr.append(rname[:-4])

failed = [value for value in missedroots if value[:-5] not in missederrs_noerr]
print('>>> Failed jobs:')
for rname in failed:
    print(rname)

submitted = []

if opts.submit:

    for rname in missedroots:
        sname = rname
        index = [n for (n, e) in enumerate(sname) if e == '_']
        sname = sname[:index[-1]]
        sname = 'sub/sub_'+sname+'.sh'
        if sname not in submitted:
            os.system('condor_submit ' + sname)
            submitted.append(sname)
     


