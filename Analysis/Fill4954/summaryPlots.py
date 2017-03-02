from sys import argv as __ARGV__
__ARGV__.append('-b')

from computeChiSquares import computeChiSquares
from gatherFromToys import gatherFromToys
from array import array
from ROOT import TFile, TColor, TCanvas, TLatex, TPaveText, TH2D, gStyle, \
                 TMultiGraph, TGraphErrors, TLegend, TH2F, TH1D

bunchcrossings = ('41', '281', '872', '1783', '2063')
beamshapes = ('SG', 'DG', 'SupG', 'TG', 'SupDG')
shapeNames = {'': 'uncorrected', \
              'SG': 'Single Gaussian', \
              'DG': 'Double Gaussian', \
              'SupG': 'Super Gaussian', \
              'TG': 'Triple Gaussian', \
              'SupDG': 'Super Double Gaussian'}

def kBird():
    red = array('d', [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764])
    green = array('d', [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832])
    blue = array('d',  [0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539])
    stops = array('d', [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])
    TColor.CreateGradientColorTable(9, stops, red, green, blue, 255)

def drawCMS(wip=False):
    text = TLatex()
    text.SetNDC()
    text.SetTextFont(62)
    text.SetTextSize(0.0375)
    text.SetTextAlign(31)
    text.DrawLatex(0.9,0.92,'2016 (13 TeV)')
    text.SetTextAlign()
    if wip:
        text.DrawLatex(0.15,0.92,'#bf{#scale[0.75]{#it{Work in Progress}}}')
    else:
        text.DrawLatex(0.15,0.92,'CMS #bf{#scale[0.75]{#it{Preliminary}}}')

def residualPlots(crossings, shapes, chiSq, dof):
    kBird()
    components = ('X1', 'Y1', 'X2', 'Y2')
    for shape in shapes:
        for bx in crossings:
            f = TFile.Open('DataAnalysisBunch'+bx+shape+'_new_StronRescale.root')
            if not f:
                continue
            for comp in components:
                hist = f.Get('res'+comp)
                for xbin in range(hist.GetXaxis().GetNbins()+1):
                    for ybin in range(hist.GetYaxis().GetNbins()+1):
                        if hist.GetBinContent(xbin, ybin) < -4.99999:
                            hist.SetBinContent(xbin, ybin, -4.99999)
                        if hist.GetBinContent(xbin, ybin) == 0.0:
                            hist.SetBinContent(xbin, ybin, -10.0)
                hist.SetTitle('')
                hist.SetName(bx+shape+'_res'+comp)
                canvas = TCanvas('c_'+hist.GetName(), '', 600, 600)
                hist.Draw('COLZ')
                canvas.Update()
                hist.GetXaxis().SetTitle('x [cm]')
                hist.GetXaxis().SetLabelSize(0.025)
                hist.GetYaxis().SetTitle('y [cm]')
                hist.GetYaxis().SetLabelSize(0.025)
                hist.GetYaxis().SetTitleOffset(1.3)
                hist.GetZaxis().SetTitle('Pulls')
                hist.GetZaxis().SetLabelSize(0.025)
                hist.GetZaxis().SetTitleOffset(0.9)
                hist.GetZaxis().SetRangeUser(-5.0,5.0)
                palette = hist.GetListOfFunctions().FindObject('palette')
                palette.SetX2NDC(0.929)
                pave = TPaveText(0.61, 0.79, 0.88, 0.88, 'NDC')
                pave.SetTextFont(42)
                pave.SetTextSize(0.025)
                pave.AddText('Scan '+comp+', BX '+bx)
                pave.AddText(shapeNames[shape]+' fit')
                redChiSq = chiSq[shape][bx] / dof[shape][bx]
                pave.AddText('#chi^{2}/d.o.f. = %6.4f'%(redChiSq))
                pave.Draw('same')
                drawCMS()
                canvas.Modified()
                canvas.Update()
                canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.pdf')
                canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.C')

def radialResidualPlots(crossings, shapes, chiSq, dof):
    kBird()
    gStyle.SetOptStat(0)
    components = ('X1', 'Y1', 'X2', 'Y2')
    for shape in shapes:
        for bx in crossings:
            f = TFile.Open('DataAnalysisBunch'+bx+shape+'_new_StronRescale.root')
            if not f:
                continue
            for comp in components:
                dataHist = f.Get('dataHist'+comp)
                modelHist = f.Get('modelHist'+comp)
                nbinsx = dataHist.GetXaxis().GetNbins()
                nbinsy = dataHist.GetYaxis().GetNbins()
                radialDat = TH1D('radialDat_'+shape+bx+comp, '', nbinsx/2, 0.0, \
                                 dataHist.GetXaxis().GetXmax())
                radialMod = TH1D('radialMod_'+shape+bx+comp, '', nbinsx/2, 0.0, \
                                 dataHist.GetXaxis().GetXmax())
                hist = TH1D('radialRes_'+shape+bx+comp, '', nbinsx/2, 0.0, \
                            dataHist.GetXaxis().GetXmax())
                radialDat.Sumw2()
                for xbin in range(nbinsx+1):
                    for ybin in range(nbinsy+1):
                        r = (dataHist.GetXaxis().GetBinCenter(xbin)**2 + \
                            dataHist.GetYaxis().GetBinCenter(ybin)**2) ** 0.5
                        radialDat.Fill(r, dataHist.GetBinContent(xbin, ybin))
                        radialMod.Fill(r, modelHist.GetBinContent(xbin, ybin))
                for rbin in range(nbinsx/2+1):
                    err = radialDat.GetBinError(rbin)
                    if err > 0.0:
                        pull = (radialDat.GetBinContent(rbin) - \
                                radialMod.GetBinContent(rbin)) / err
                    else:
                        pull = 0.0
                    hist.SetBinContent(rbin, pull)
                canvas = TCanvas('c_'+hist.GetName(), '', 600, 600)
                hist.Draw('HF')
                canvas.Update()
                hist.SetFillColor(4)
                hist.SetLineColor(1)
                hist.GetXaxis().SetTitle('r [cm]')
                hist.GetXaxis().SetLabelSize(0.025)
                hist.GetYaxis().SetTitle('Pulls')
                hist.GetYaxis().SetLabelSize(0.025)
                hist.GetYaxis().SetTitleOffset(1.1)
                hist.GetYaxis().SetRangeUser(-1.5, 5.0)
                pave = TPaveText(0.15, 0.79, 0.42, 0.88, 'NDC')
                pave.SetTextFont(42)
                pave.SetTextSize(0.025)
                pave.AddText('Scan '+comp+', BX '+bx)
                pave.AddText(shapeNames[shape]+' fit')
                redChiSq = chiSq[shape][bx] / dof[shape][bx]
                pave.AddText('#chi^{2}/d.o.f. = %6.4f'%(redChiSq))
                pave.Draw('same')
                drawCMS()
                canvas.Modified()
                canvas.Update()
                canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.pdf')
                canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.C')

def exampleDataPlot(bx, shape, comp):
    kBird()
    f = TFile.Open('DataAnalysisBunch'+bx+shape+'_new_StronRescale.root')
    if f:
        hist = f.Get('dataHist'+comp)
        hist.SetTitle('')
        hist.SetName(bx+shape+'_dataHist'+comp)
        canvas = TCanvas('c_'+hist.GetName(), '', 600, 600)
        canvas.SetFrameFillColor(0)
        hist.Draw("COLZ")
        canvas.Update()
        hist.GetXaxis().SetTitle('x [cm]')
        hist.GetXaxis().SetLabelSize(0.025)
        hist.GetYaxis().SetTitle('y [cm]')
        hist.GetYaxis().SetLabelSize(0.025)
        hist.GetYaxis().SetTitleOffset(1.3)
        hist.GetZaxis().SetTitle('Number of Vertices')
        hist.GetZaxis().SetLabelSize(0.025)
        hist.GetZaxis().SetTitleOffset(0.7)
        hist.GetZaxis().SetRangeUser(0.0,240.0)
        hist.GetZaxis().CenterTitle()
        hist.GetZaxis().SetNdivisions(1, False)
        palette = hist.GetListOfFunctions().FindObject('palette')
        palette.SetX2NDC(0.929)
        pave = TPaveText(0.65, 0.82, 0.88, 0.88, 'NDC')
        pave.SetTextFont(42)
        pave.SetTextSize(0.025)
        pave.AddText('Scan '+comp+', BX '+bx)
        pave.AddText('Measured data')
        pave.Draw('same')
        drawCMS()
        canvas.Modified()
        canvas.Update()
        canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.pdf')
        canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.C')

def makeShapesCrossingsHist(crossings, shapes, name):
    kBird()
    gStyle.SetOptStat(0)
    nShapes = len(shapes)
    nCrossings = len(crossings)
    hist = TH2D(name, '', nCrossings, 0, nCrossings-1, \
                             nShapes, 0, nShapes-1)
    for i, bx in enumerate(crossings):
        hist.GetXaxis().SetBinLabel(i+1, bx)
    for j, shape in enumerate(shapes):
        hist.GetYaxis().SetBinLabel(nShapes-j, shape)
    return hist

def chiSqPlot(crossings, shapes, chiSq, dof):
    nShapes = len(shapes)
    hist = makeShapesCrossingsHist(crossings, shapes, 'chisq')
    for i, bx in enumerate(crossings):
        for j, shape in enumerate(shapes):
            if chiSq[shape][bx]:
                redChiSq = chiSq[shape][bx] / dof[shape][bx]
                hist.SetBinContent(i+1, nShapes-j, redChiSq)
    canvas = TCanvas('c_'+hist.GetName(), '', 600, 600)
    hist.Draw('TEXTCOLZ')
    canvas.Update()
    hist.GetXaxis().SetNdivisions(len(crossings), False)
    hist.GetYaxis().SetNdivisions(nShapes, False)
    hist.GetZaxis().SetTitle('#chi^{2} / d.o.f.')
    hist.GetZaxis().SetLabelSize(0.025)
    hist.GetZaxis().SetTitleOffset(0.5)
    hist.GetZaxis().SetRangeUser(1.00,1.12)
    hist.GetZaxis().CenterTitle()
    hist.GetZaxis().SetNdivisions(1, False)
    palette = hist.GetListOfFunctions().FindObject('palette')
    palette.SetX2NDC(0.929)
    drawCMS()
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.pdf')
    canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.C')

def correctionPlot(crossings, shapes, overDiff):
    nShapes = len(shapes)
    hist = makeShapesCrossingsHist(crossings, shapes, 'corrections')
    for i, bx in enumerate(crossings):
        for j, shape in enumerate(shapes):
            if overDiff[shape][bx]:
                correction = -100.0 * overDiff[shape][bx]
                hist.SetBinContent(i+1, nShapes-j, correction)
            else:
                hist.SetBinContent(i+1, nShapes-j, -1e6)
    canvas = TCanvas('c_'+hist.GetName(), '', 600, 600)
    hist.Draw('TEXTCOLZ')
    canvas.Update()
    hist.GetXaxis().SetNdivisions(len(crossings), False)
    hist.GetYaxis().SetNdivisions(nShapes, False)
    hist.GetZaxis().SetTitle('correction on overlap integral [%]')
    hist.GetZaxis().SetLabelSize(0.025)
    hist.GetZaxis().SetTitleOffset(0.5)
    hist.GetZaxis().SetRangeUser(-1.2,0.0)
    hist.GetZaxis().CenterTitle()
    hist.GetZaxis().SetNdivisions(1, False)
    palette = hist.GetListOfFunctions().FindObject('palette')
    palette.SetX2NDC(0.929)
    drawCMS()
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.pdf')
    canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.C')

def correctedCrossSectionsPlot(crossings, shapes, overDiff):
    uncorrected = {41: (3.2575816286, 0.00514858611944), \
                   281: (3.26316215713, 0.00468789412223), \
                   872: (3.27340775031, 0.00484925398906), \
                   1783: (3.24986926821, 0.00460908436455), \
                   2063: (3.26363843728, 0.0044071069983)}
    multi = TMultiGraph('sigmavis', '')
    graphs = []
    n = len(shapes) + 1
    for i, shape in enumerate([''] + list(shapes)):
        xval = array('d', [a+0.08*(i-0.5*n) for a in range(len(crossings))])
        xerr = array('d', len(crossings)*[0])
        yval = array('d', [uncorrected[int(bx)][0] for bx in crossings])
        yerr = array('d', [uncorrected[int(bx)][1] for bx in crossings])
        if shape:
            for j, bx in enumerate(crossings):
                yval[j] *= 1 + overDiff[shape][bx]
                yerr[j] *= 1 + overDiff[shape][bx]
        graph = TGraphErrors(len(crossings), xval, yval, xerr, yerr)
        graph.SetName('ge'+shape)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(1+i)
        multi.Add(graph)
        graphs.append(graph)
    gStyle.SetOptStat(0)
    hist = TH2F('hist', '', len(crossings), -0.5, len(crossings)-0.5, 100, 3.23, 3.33)
    for i, bx in enumerate(crossings):
        hist.GetXaxis().SetBinLabel(i+1, bx)
    canvas = TCanvas('c_'+multi.GetName(), '', 600, 600)
    hist.Draw('AXIS')
    multi.Draw('P')
    canvas.Update()
    hist.GetXaxis().SetLabelSize(0.035)
    hist.GetXaxis().SetNdivisions(len(crossings), False)
    hist.GetYaxis().SetTitle('#sigma_{vis} [b]')
    hist.GetYaxis().SetLabelSize(0.025)
    hist.GetYaxis().SetTitleOffset(1.3)
    leg = TLegend(0.15, 0.82, 0.85, 0.85)
    leg.SetNColumns(len(shapes)+1)
    leg.SetBorderSize(0)
    for i, shape in enumerate([''] + list(shapes)):
        entry = leg.AddEntry('ge'+shape, shapeNames[shape], 'P')
        entry.SetMarkerStyle(20)
        entry.SetMarkerColor(1+i)
    leg.Draw()
    drawCMS()
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.pdf')
    canvas.SaveAs('summaryPlots/'+canvas.GetName()+'.C')

def summaryPlots(crossings, shapes):
    chiSq, dof = computeChiSquares(crossings, shapes)
    overDiff = gatherFromToys(crossings, shapes)
    residualPlots(crossings, shapes, chiSq, dof)
    radialResidualPlots(crossings, shapes, chiSq, dof)
    chiSqPlot(crossings, shapes, chiSq, dof)
    correctionPlot(crossings, shapes, overDiff)
    exampleDataPlot('41', 'DG', 'X1')
    correctedCrossSectionsPlot(crossings, ('DG','TG','SupG','SupDG'), overDiff)

if __name__ == '__main__':
    summaryPlots(bunchcrossings, beamshapes)
