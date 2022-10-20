import os

color_dict = {}
color_dict['1'] = 'r.kRed+1'
color_dict['10'] = 'r.kOrange+1'
color_dict['100'] = 'r.kGreen+1'
color_dict['1000'] = 'r.kBlue+1'
color_dict['10000'] = 'r.kMagenta+1'

space = '      '

eosdir = '/eos/user/f/fernance/Galapago/Sep2022/'
prefix = 'ggH_HToSSTo4l'
outdir = 'NTuples-Galapago_2017' # adapt year here!
ndir   = '0000'

outdat = 'CentralHSS_2017UL_Fall22.dat' # adapt year here!
outdat_ = open(outdat, 'w')

## loop over samples:

for _d in os.listdir(eosdir):
    sample_dir = eosdir
    if _d[:13] != prefix: continue
    sample_dir += _d + '/'
    sample_dir += outdir + '/'

    mH   = (sample_dir.split('MH-')[-1]).split('_')[0]
    mS   = (sample_dir.split('MS-')[-1]).split('_')[0]
    ctau = (sample_dir.split('ctauS-')[-1]).split('_')[0]
    sample_name = 'HSS_' + mH + '_' + mS + '_' + ctau + '_2017' # adapt year here!
    sample_label = 'HSS({0},{1},{2})'.format(mH, mS, ctau)

    if os.path.isdir(sample_dir):
        #sample_dir += os.listdir(sample_dir)[-1] + '/0000/'
        sample_dir += os.listdir(sample_dir)[-1] + '/'
        line = sample_name 
        line += space
        line += color_dict[ctau]
        line += space
        line += sample_name
        line += space
        line += sample_label
        line += space
        line += sample_dir
        line += space
        line += '1'
        line += space
        line += '0' 
        line += '\n'
        print(line)
        outdat_.write(line)

outdat_.close()


