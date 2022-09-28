import ROOT as r
import numpy as np

###
### --- color arrays definition
###
gcolors = {}
gcolors['blue']   = r.TColor.GetColor('#4285f4');
gcolors['red']    = r.TColor.GetColor('#ea4335');
gcolors['yellow'] = r.TColor.GetColor('#fbbc04');
gcolors['green']  = r.TColor.GetColor('#34a853');

dcolors = {}
dcolors['1mm'] = r.kAzure+10
dcolors['10mm'] = r.kAzure-4
dcolors['100mm'] = r.kBlue-3
dcolors['1000mm'] = r.kMagenta-3
dcolors['10000mm'] = r.kMagenta+2


###
### --- 2D plot palette definition
###

# Signal palettes

sigc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),
                                            np.array([1.00, 204./255.]),
                                            np.array([1.00, 0.00]),
                                            np.array([1.00, 0.00]), 255);
sigpalette_ = []
for i in range(0, 255):
    sigpalette_.append(sigc + i)
sigpalette = np.array(sigpalette_, dtype=np.int32)


# Background (predicted) palette

bkgc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),
                                            np.array([1.00, 0.00]),
                                            np.array([1.00, 153/255.]),
                                            np.array([1.00, 153./255.]), 255);

bkgpalette_ = []
for i in range(0, 255):
    sigpalette_.append(sigc + i)

bkgpalette = np.array(bkgpalette_, dtype=np.int32)


