import ROOT
from CombineHarvester.CombineTools.plotting import *
import optparse
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

##### Init style definitions

CMSstyle = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : R.kBlue, 'LineStyle' : 4},
               'exp1' : { 'FillColor' : R.kGreen+1},
               'exp2' : { 'FillColor' : R.kOrange}
               }

magenta_exp0 = CreateTransparentColor(R.kMagenta+1, 0.9)
magenta_exp1 = CreateTransparentColor(R.kMagenta-4, 0.3)
magenta_exp2 = CreateTransparentColor(R.kMagenta-9, 0.3)

style_magenta_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : magenta_exp0},
               'exp1' : { 'FillColor' : magenta_exp1},
               'exp2' : { 'FillColor' : magenta_exp2}
               }

red_exp0 = CreateTransparentColor(R.kRed+1, 0.9)
red_exp1 = CreateTransparentColor(R.kRed-4, 0.3)
red_exp2 = CreateTransparentColor(R.kRed-9, 0.3)

style_red_dict = {
               'obs' : { 'LineWidth' : 2},
               'exp0' : { 'LineWidth' : 2, 'LineColor' : red_exp0},
               'exp1' : { 'FillColor' : red_exp1},
               'exp2' : { 'FillColor' : red_exp2}
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
styles.append(style_blue_dict)
styles.append(style_red_dict)
styles.append(style_green_dict)

##### Combined legend definition

legend_dict = {
         'obs' : { 'Label' : 'Observed {0}', 'LegendStyle' : 'LP', 'DrawStyle' : 'PLSAME'},
         'exp0' : { 'Label' : 'Median expected{0}', 'LegendStyle' : 'L', 'DrawStyle' : 'LSAME'},
         'exp1' : { 'Label' : '68% expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'},
         'exp2' : { 'Label' : '95% expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'}
         }

legend_red = {
         'obs' : { 'Label' : '{0} (Observed)', 'LegendStyle' : 'LP', 'DrawStyle' : 'PLSAME'},
         'exp0' : { 'Label' : '{0} (Expected)', 'LegendStyle' : 'L', 'DrawStyle' : 'LSAME'},
         'exp1' : { 'Label' : '68% expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'},
         'exp2' : { 'Label' : '95% expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'}
         }

def setLegendLabel(legend, label):
    legend['obs']['Label'] = legend['obs']['Label'].format(label)
    legend['exp0']['Label'] = legend['exp0']['Label'].format(label)

#
# --- Main ()
#

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-j', '--json', action='store', type=str, dest='json', default='', help='json files')
    parser.add_option('-m', '--mH', action='store', type=str, dest='mH', default='', help='mH')
    parser.add_option('-f', '--flavor', action='store', type=str, dest='flavor', default='', help='flavor')
    parser.add_option('-o', '--outdir', action='store', type=str, dest='outdir', default='', help='out')
    parser.add_option('-y', '--year', action='store', type=str, dest='year', default='2016', help='year')
    (opts, args) = parser.parse_args()

    # Style and pads
    ModTDRStyle()
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetLegendTextSize(0.033)
    canv = ROOT.TCanvas('limit', 'limit')
    canv.SetFillStyle(4000);
    canv.SetFrameFillColor(4000);
    canv.SetFrameFillStyle(4000);
    pads = OnePad()
 
    # Get limit TGraphs as a dictionary
    jsons = opts.json.split(',')
    graphs = []
    masslabels = []
    order = []
    for json in jsons:
        print(json)
        masslabels.append('m_{S} = ' + json.split('mS')[-1].split('__')[0] + ' GeV')
        order.append(int(json.split('mS')[-1].split('__')[0]))
        graphs.append(StandardLimitsFromJSONFile(json))

    graphs = [x for _,x in sorted(zip(order,graphs))]
    masslabels = [x for _,x in sorted(zip(order,masslabels))]
    jsons = [x for _,x in sorted(zip(order,jsons))]

    # Create an empty TH1 from the first TGraph to serve as the pad axis and frame
    axis = CreateAxisHist(graphs[0].values()[0])
    axis.GetXaxis().SetTitle('c#tau [cm]')
    axis.GetYaxis().SetTitle('#sigma(H)xB(H#rightarrowSS)')
    axis.GetYaxis().SetRangeUser(0.0001, 30)
    axis.GetXaxis().SetRangeUser(0.004, 100000)
    axis.SetFillColor(4000);
    pads[0].cd()
    pads[0].SetFillStyle(4000)
    pads[0].SetFrameFillStyle(4000)
    axis.Draw('axis')
 
    # Create a legend in the top left
    if len(graphs) == 1:
        legend = PositionedLegend(0.45, 0.15, 3, 0.045)
    else:
        legend = PositionedLegend(0.45, 0.07*len(graphs), 3, 0.1, 0.045)
 
    # Set the standard green and yellow colors and draw
    if len(graphs) == 1:
        StyleLimitBand(graphs[0], overwrite_style_dict=CMSstyle)
        setLegendLabel(legend_dict, ' ')
        DrawLimitBand(pads[0], graphs[0], draw=['exp2', 'exp1', 'exp0'], legend=legend, legend_overwrite=legend_dict)
    else:
        for n in range(0, len(graphs)):

            StyleLimitBand(graphs[n], overwrite_style_dict=styles[n])
            ## Draw in pad
            pads[0].cd()
            for key in ['exp0', 'exp1']:
                legend_label = ''
                if key == 'obs':
                    legend_label = masslabels[n] + ' (Observed)'
                elif key == 'exp0':
                    legend_label = masslabels[n] + ' (Expected)'
                else:
                    legend_label = legend_dict[key]['Label']
                print(legend_label)
                graphs[n][key].Draw(legend_dict[key]['DrawStyle'])
                legend.AddEntry(graphs[n][key], legend_label, legend_dict[key]['LegendStyle'])
            #DrawLimitBand(pads[0], graphs[n], draw=['exp1', 'exp0'], legend=legend, legend_overwrite=legend_red)
    legend.Draw()
 
    # Re-draw the frame and tick marks
    pads[0].RedrawAxis()
    pads[0].GetFrame().Draw()
    pads[0].SetLogx(1)
    pads[0].SetLogy(1)
    pads[0].SetTickx(1)
    pads[0].SetTicky(1)
 
    # Adjust the y-axis range such that the maximum graph value sits 25% below
    # the top of the frame. Fix the minimum to zero.
    #FixBothRanges(pads[0], GetPadYMin(pads[0]), 0.15, GetPadYMax(pads[0]), 0.5)
    #FixBothRanges(pads[0], GetPadYMin(pads[0]), 0.15, GetPadYMax(pads[0]), 0.3)
 
    # Standard CMS logo
    # DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.045, 0.035, 1000.0, '', 0.8)
    CMSlabel = ROOT.TLatex()
    CMSlabel.SetNDC();
    CMSlabel.SetTextAngle(0);
    CMSlabel.SetTextColor(ROOT.kBlack);
    CMSlabel.SetTextFont(42);
    CMSlabel.SetTextAlign(22);
    CMSlabel.SetTextSize(0.06);
    CMSlabel.DrawLatex(0.28, 0.89, "#bf{CMS}")
    CMSextralabel = ROOT.TLatex()
    CMSextralabel.SetNDC();
    CMSextralabel.SetTextAngle(0);
    CMSextralabel.SetTextColor(ROOT.kBlack);
    CMSextralabel.SetTextFont(42);
    CMSextralabel.SetTextAlign(22);
    CMSextralabel.SetTextSize(0.04);
    CMSextralabel.DrawLatex(0.31, 0.84, "#it{Preliminary}")
    
    # Channel label
    Channellabel = ROOT.TLatex()
    Channellabel.SetNDC();
    Channellabel.SetTextAngle(0);
    Channellabel.SetTextColor(ROOT.kBlack);
    Channellabel.SetTextFont(42);
    Channellabel.SetTextAlign(13);
    Channellabel.SetTextSize(0.04);
    if opts.flavor == 'Electron':
        Channellabel.DrawLatex(0.24, 0.79, "ee channel") 
    elif opts.flavor == 'Muon':
        Channellabel.DrawLatex(0.24, 0.79, "#mu#mu channel")
    #elif opts.flavor == 'Joint':
    #    Channellabel.DrawLatex(0.20, 0.79, "Muon + electron channel")

    # Year label
    Yearlabel = ROOT.TLatex()
    Yearlabel.SetNDC();
    Yearlabel.SetTextAngle(0);
    Yearlabel.SetTextColor(ROOT.kBlack);
    Yearlabel.SetTextFont(42);
    Yearlabel.SetTextAlign(33);
    Yearlabel.SetTextSize(0.04);
    if opts.year == '2016': 
        Yearlabel.DrawLatex(0.96, 0.99, "35.9 fb^{-1} (13 TeV)")
    elif opts.year == '2017':
        Yearlabel.DrawLatex(0.96, 0.99, "41.5 fb^{-1} (13 TeV)")
    elif opts.year == '2018':
        Yearlabel.DrawLatex(0.96, 0.99, "59.7 fb^{-1} (13 TeV)")
    elif opts.year == 'Full':
        Yearlabel.DrawLatex(0.96, 0.99, "137 fb^{-1} (13 TeV)")

    CLlabel = ROOT.TLatex()
    CLlabel.SetNDC();
    CLlabel.SetTextAngle(0);
    CLlabel.SetTextColor(ROOT.kBlack);
    CLlabel.SetTextFont(42);
    CLlabel.SetTextAlign(13);
    CLlabel.SetTextSize(0.04);
    CLlabel.DrawLatex(0.16, 0.98, "95% CL upper limits")

    # Model label
    Modellabel = ROOT.TLatex()
    Modellabel.SetNDC();
    Modellabel.SetTextAngle(0);
    Modellabel.SetTextColor(ROOT.kBlack);
    Modellabel.SetTextFont(42);
    Modellabel.SetTextAlign(13);
    Modellabel.SetTextSize(0.037);
    Modellabel.DrawLatex(0.20, 0.55, "H\\rightarrow SS, S\\rightarrow\\ell^{+}\\ell^{-}, \\ell = e,\\mu") # Hardcoded

    Masslabel = ROOT.TLatex()
    Masslabel.SetNDC();
    Masslabel.SetTextAngle(0);
    Masslabel.SetTextColor(ROOT.kBlack);
    Masslabel.SetTextFont(42);
    Masslabel.SetTextAlign(13);
    Masslabel.SetTextSize(0.035);
    if len(jsons) == 1:
        json = jsons[0]
        mS = json.split('__mS')[1].split('__')[0]
        Masslabel.DrawLatex(0.20, 0.7, "m_{{H}} = {0} GeV, m_{{S}} = {1} GeV".format(opts.mH, mS))
    else:
        Masslabel.DrawLatex(0.5, 0.9, "m_{{H}} = {0} GeV".format(opts.mH))
        


    # Re-draw axis
    axis.Draw('axis, same')
 
    canv.Update()
    canv.Modified()
    canv.Print(opts.outdir + 'limits_' + opts.flavor + opts.mH + opts.year + '.pdf')
    canv.Print(opts.outdir + 'limits_' + opts.flavor + opts.mH + opts.year + '.png')
