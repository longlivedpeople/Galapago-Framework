import ROOT as r
import numpy as np

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


