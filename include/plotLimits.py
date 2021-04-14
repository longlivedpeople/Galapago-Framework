import ROOT
from CombineHarvester.CombineTools.plotting import *
import optparse
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

##### Init style definitions

magenta_exp0 = CreateTransparentColor(R.kMagenta+1, 0.9)
magenta_exp1 = CreateTransparentColor(R.kMagenta-4, 0.3)
magenta_exp2 = CreateTransparentColor(R.kMagenta-9, 0.3)

style_magenta_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : magenta_exp0},
               'exp1' : { 'FillColor' : magenta_exp1},
               'exp2' : { 'FillColor' : magenta_exp2}
               }


green_exp0 = CreateTransparentColor(R.kGreen+1, 0.9)
green_exp1 = CreateTransparentColor(R.kGreen-4, 0.3)
green_exp2 = CreateTransparentColor(R.kGreen-9, 0.3)

style_green_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : green_exp0},
               'exp1' : { 'FillColor' : green_exp1},
               'exp2' : { 'FillColor' : green_exp2}
               }


blue_exp0 = CreateTransparentColor(R.kBlue+1, 0.9)
blue_exp1 = CreateTransparentColor(R.kBlue-4, 0.3)
blue_exp2 = CreateTransparentColor(R.kBlue-9, 0.3)

style_blue_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : blue_exp0},
               'exp1' : { 'FillColor' : blue_exp1},
               'exp2' : { 'FillColor' : blue_exp2}
               }

styles = []
styles.append(style_magenta_dict)
styles.append(style_blue_dict)
styles.append(style_green_dict)

##### Combined legend definition

legend_dict = {
         'obs' : { 'Label' : 'Observed {0}', 'LegendStyle' : 'LP', 'DrawStyle' : 'PLSAME'},
         'exp0' : { 'Label' : 'Expected {0}', 'LegendStyle' : 'L', 'DrawStyle' : 'LSAME'},
         'exp1' : { 'Label' : '#pm1#sigma Expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'},
         'exp2' : { 'Label' : '#pm2#sigma Expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'}
         }

def setLegendLabel(legend, label):
    legend_dict['obs']['Label'].format(label)
    legend_dict['exp0']['Label'].format(label)

#
# --- Main ()
#

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-j', '--json', action='store', type=str, dest='json', default='', help='json files')
    parser.add_option('-m', '--mH', action='store', type=str, dest='mH', default='', help='mH')
    parser.add_option('-f', '--flavor', action='store', type=str, dest='flavor', default='', help='flavor')
    parser.add_option('-o', '--outdir', action='store', type=str, dest='outdir', default='', help='out')
    (opts, args) = parser.parse_args()

    # Style and pads
    ModTDRStyle()
    canv = ROOT.TCanvas('limit', 'limit')
    pads = OnePad()
 
    # Get limit TGraphs as a dictionary
    jsons = opts.json.split(',')
    graphs = []
    for json in jsons:
        graphs.append(StandardLimitsFromJSONFile(json))

    # Create an empty TH1 from the first TGraph to serve as the pad axis and frame
    axis = CreateAxisHist(graphs[0].values()[0])
    axis.GetXaxis().SetTitle('c#tau (cm)')
    axis.GetYaxis().SetTitle('95% CL limit on #sigma')
    axis.GetYaxis().SetRangeUser(0.0001, 30)
    axis.GetXaxis().SetRangeUser(0.004, 100000)
    pads[0].cd()
    axis.Draw('axis')
 
    # Create a legend in the top left
    legend = PositionedLegend(0.3, 0.2, 3, 0.015)
 
    # Set the standard green and yellow colors and draw
    if len(graphs) == 1:
        StyleLimitBand(graphs[0])
        DrawLimitBand(pads[0], graphs[0], draw=['exp2', 'exp1', 'exp0'], legend=legend, legend_overwrite=legend_dict)
    else:
        for n in range(0, len(graphs)):
            StyleLimitBand(graphs[n], overwrite_style_dict=styles[n])
            DrawLimitBand(pads[0], graphs[n], draw=['exp2', 'exp1', 'exp0'], legend=legend, legend_overwrite=legend_dict)
    legend.Draw()
 
    # Re-draw the frame and tick marks
    pads[0].RedrawAxis()
    pads[0].GetFrame().Draw()
    pads[0].SetLogx(1)
    pads[0].SetLogy(1)
 
    # Adjust the y-axis range such that the maximum graph value sits 25% below
    # the top of the frame. Fix the minimum to zero.
    #FixBothRanges(pads[0], GetPadYMin(pads[0]), 0.15, GetPadYMax(pads[0]), 0.5)
    #FixBothRanges(pads[0], GetPadYMin(pads[0]), 0.15, GetPadYMax(pads[0]), 0.3)
 
    # Standard CMS logo
    DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.035, 1.2, '', 0.8)

    # Re-draw axis
    axis.Draw('axis, same')
 
    canv.Print(opts.outdir + 'limits_' + opts.flavor + opts.mH + '.pdf')
    canv.Print(opts.outdir + 'limits_' + opts.flavor + opts.mH + '.png')
