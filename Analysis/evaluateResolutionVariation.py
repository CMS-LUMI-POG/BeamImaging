def drawCMS(wip=False):
    from ROOT import TLatex
    text = TLatex()
    text.SetNDC()
    text.SetTextFont(62)
    text.SetTextSize(0.0375)
    text.SetTextAlign(33)
    text.DrawLatex(0.88,0.88,'2016 (13 TeV)')
    text.SetTextAlign(13)
    if wip:
        text.DrawLatex(0.15,0.88,'#bf{#scale[0.75]{#it{Work in Progress}}}')
    else:
        text.DrawLatex(0.15,0.88,'CMS #bf{#scale[0.75]{#it{Preliminary}}}')

def evaluateResolutionVariation(model, bcid, prefix, suffix, vtxres, scaling, \
                                legend=False):
    from array import array
    from ROOT import TFile, TMultiGraph, TGraphErrors, gStyle, TCanvas, \
                     TLegend, TLatex

    names = [pre+'_'+model+'_'+bcid+'_'+suf for (pre, suf) in zip(prefix, \
             suffix)]
    xVal = {name: array('d', [0.0]*len(vtxres)) for name in names}
    xErr = {name: array('d', [0.0]*len(vtxres)) for name in names}
    yVal = {name: array('d', [0.0]*len(vtxres)) for name in names}
    yErr = {name: array('d', [0.0]*len(vtxres)) for name in names}

    for name, scale in zip(names, scaling):
        for i, vr in enumerate(vtxres):
            f = TFile('results/'+name+'_'+vr+'.root')
            yVal[name][i] = f.Get('h_overlapInt').GetMean() * 100.0
            yErr[name][i] = f.Get('h_integ').GetMeanError() * 100.0
            xVal[name][i] = f.Get('h_vtxRes').GetMean() * scale * 1.0e4

    multi = TMultiGraph('overlapIntegralBcid'+bcid, '')
    graphs = []
    for i, name in enumerate(names):
        graph = TGraphErrors(len(vtxres), xVal[name], yVal[name], xErr[name], \
                             yErr[name])
        graph.SetName(name)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(1+i)
        multi.Add(graph)
        graphs.append(graph)

    gStyle.SetOptStat(0)
    canvas = TCanvas(model+'_'+bcid, '', 600, 600)
    multi.Draw('AP')
    canvas.Update()
    multi.GetXaxis().SetTitle('vertex resolution [#mum]')
    multi.GetXaxis().SetLabelSize(0.025)
    multi.GetXaxis().SetRangeUser(11, 69)
    multi.GetYaxis().SetTitle('overlap integral [a.u.]')
    multi.GetYaxis().SetLabelSize(0.025)
    multi.GetYaxis().SetRangeUser(0.77, 1.43)
    if legend:
        leg = TLegend(0.55, 0.15, 0.88, 0.3)
        leg.SetBorderSize(0)
        for i, name in enumerate(names):
            entry = leg.AddEntry(name, legend[i], 'P')
            entry.SetMarkerStyle(20)
            entry.SetMarkerColor(1+i)
        leg.Draw()
    drawCMS(wip=True)
    text = TLatex()
    text.SetNDC()
    text.SetTextFont(62)
    text.SetTextSize(0.04)
    text.SetTextAlign(21)
    text.DrawLatex(0.5, 0.92, 'Vertex Resolution Study: '+model+', BCID '+bcid)
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('plots/'+canvas.GetName()+'.pdf')
    canvas.SaveAs('plots/'+canvas.GetName()+'.C')


def main():
    from sys import argv as __ARGV__
    __ARGV__.append('-b')

    from tools import loadJson
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Plot overlap integral as function '+ \
                            'of vertex resolution')
    parser.add_argument('-b', action='store_true', help='enable batch mode')
    parser.add_argument('json', nargs='+', help='specify one or more JSON '+ \
                        'files containing config information')
    parser.add_argument('-m', '--model', required=True, choices=['SingleGauss', \
                        'SingleGaussUncorrelated', 'DoubleGauss', 'SuperGauss', \
                        'TripleGauss', 'SuperDoubleGauss'], help='specify '+ \
                        'fit model')
    parser.add_argument('-c', '--bcid', required=True, nargs='+', \
                        help='list one or more bunch crossings')
    parser.add_argument('-r', '--vtxres', nargs='+', choices=['default', 'low', \
                        'high', 'half', 'onehalf', 'double'], \
                        default=['default'], help='specify one or more '+ \
                        'options to modify vertex resolution')
    args = parser.parse_args()

    from shapes.SingleGauss import SingleGauss, SingleGaussUncorrelated
    from shapes.DoubleGauss import DoubleGauss, SuperGauss
    from shapes.TripleGauss import TripleGauss, SuperDoubleGauss

    jsons = [loadJson(filename) for filename in args.json]
    prefix = [str(json['prefix']) for json in jsons]
    suffix = [str(json['suffix']) for json in jsons]
    scaling = [float(json['scaling']) for json in jsons]
    model = locals()[args.model].Shortname
    vtxres = [str(v) for v in args.vtxres]
    legend = [str(json['legend']) for json in jsons]

    for bcid in args.bcid:
        evaluateResolutionVariation(model, bcid, prefix, suffix, vtxres, \
                                    scaling, legend)

if __name__ == '__main__':
    main()
