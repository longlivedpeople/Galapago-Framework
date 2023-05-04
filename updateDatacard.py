##########################################################################
# Quick hack to add the statistical uncertainty on the background yields #
# as a new nuissance parameter.                                          #
##########################################################################
import optparse
import sys
import math

if __name__=='__main__':

    
    parser = optparse.OptionParser(usage='usage: %prog [options] path', version='%prog 1.0')
    parser.add_option('-i', '--input'     , action='store'   , type='string', dest='input',   default='',         help='Input configuration file.')
    parser.add_option('-o', '--output'    , action='store'   , type='string', dest='output',  default='',         help='Output configuration file.')
    (opts, args) = parser.parse_args()

    if opts.input == '' or opts.output == '':
        print('Please specify both input and output file names')
        sys.exit(0)

    nbins = 0
    nprocess = 0
    addline = ''
    toput = []
    toadd = []
    bkg = []
    names = []
    for i in open(opts.input).readlines():
        newline = i
        if i.find('bin   ') != -1:
            namesSplitted = i.split()
            for j in range(nbins):
                names.append(namesSplitted[1 + 2*j])
        if i.find('imax') != -1:
            nbins = int(i.split()[1])
        if i.find('jmax') != -1:
            nprocess = int(i.split()[1])    
        if i.find('nuisance') != -1:
            number = int(i.split()[1]) + nbins
            newline = i.replace(i.split()[1], str(number))
        if i.find('rate') != -1:
            splited = i.split()
            for j in range(nbins):
                bkg.append(float(splited[2 + 2 * j]))
        if i.find('bkg') != -1 and i.find('lnN') != -1:
            splited = i.split()
            for k in range(nbins):
                addline = 'bkg_stat_' + names[k] + '    lnN'
                separator = '               -                   '
                for j in range(nbins):
                    if j == k:  
                        addline = addline + separator + str(round(1.0 + 1.0/math.sqrt(bkg[j]), 2))
                    else:
                        addline = addline + separator + '-'
                toadd.append(addline + '\n')
        toput.append(newline)
    final = toput + toadd
    thefile = open(opts.output, 'w')
    for i in final:
        thefile.write(i)
    thefile.close()



    
 
        
                  


