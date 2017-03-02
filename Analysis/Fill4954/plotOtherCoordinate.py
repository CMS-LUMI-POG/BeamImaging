from sys import argv as __ARGV__
__ARGV__.append('-b')

from array import array
from ROOT import TFile, gStyle, TF1, TGraphErrors, TMultiGraph, TCanvas, \
                 TLegend, TLatex

filename = '/eos/cms/store/user/jsalfeld/vdmScan_2016ReRecoJan/mergedJan172016.root'
bunchcrossings = ('41', '281', '872', '1783', '2063')
biscans = ('1', '2')
scaling = 0.00458

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

def plotOtherCoordinate(datafile, crossings, scans):
    gStyle.SetOptStat(0)
    f = TFile.Open(datafile)
    for bx in crossings:
        for coord in ('X', 'Y'):
            if 'X' in coord:
                other = 'Y'
            else:
                other = 'X'
            meanGraphs = []
            rmsGraphs = []
            for scan in scans:
                meanVal = array('d')
                meanErr = array('d')
                rmsVal = array('d')
                rmsErr = array('d')
                stepVal = array('d', [1.0*a for a in range(1,20)])
                stepErr = array('d', 19*[0.0])
                for i in range(1,20):
                    hist = f.Get('scan'+coord+scan+'Move_b'+bx+'_'+str(i))
                    funcGaus = TF1('funcGaus'+coord+scan, 'gaus', -10, 10)
                    getattr(hist, 'Projection'+other)().Fit('funcGaus'+coord+scan)
                    meanVal.append(funcGaus.GetParameter('Mean')*scaling)
                    meanErr.append(funcGaus.GetParError(1)*scaling)
                    rmsVal.append(funcGaus.GetParameter('Sigma')*scaling)
                    rmsErr.append(funcGaus.GetParError(2)*scaling)
                meanGraph = TGraphErrors(19, stepVal, meanVal, stepErr, meanErr)
                rmsGraph = TGraphErrors(19, stepVal, rmsVal, stepErr, rmsErr)
                meanGraph.SetName('graph'+coord+scan+'MeanBX'+bx)
                rmsGraph.SetName('graph'+coord+scan+'RmsBX'+bx)
                meanGraphs.append(meanGraph)
                rmsGraphs.append(rmsGraph)
            meanMulti = TMultiGraph('multi'+coord+scan+'MeanBX'+bx, '')
            rmsMulti = TMultiGraph('multi'+coord+scan+'RmsBX'+bx, '')
            for i, gr in enumerate(meanGraphs):
                gr.SetMarkerStyle(3+19*i)
                gr.SetMarkerSize(1)
                gr.SetMarkerColor(2+i)
                meanMulti.Add(gr)
            for i, gr in enumerate(rmsGraphs):
                gr.SetMarkerStyle(3+19*i)
                gr.SetMarkerSize(1)
                gr.SetMarkerColor(2+i)
                rmsMulti.Add(gr)
            meanCanvas = TCanvas('canvasMean'+coord+'ScanBX'+bx, '', 1000, 600)
            meanMulti.Draw('APE')
            meanCanvas.Update()
            meanMulti.GetYaxis().SetTitle('luminous region '+other.lower()+'-mean [cm]')
            meanMulti.GetXaxis().SetTitle('scan point')
            meanMulti.GetYaxis().SetTitleOffset(1.3)
            meanMulti.GetYaxis().SetRangeUser(-0.01, 0.01)
            meanMulti.GetYaxis().SetTickLength(0.02)
            meanLeg = TLegend(0.4, 0.13, 0.6, 0.28)
            meanLeg.SetHeader(coord+' Scan (BCID '+bx+')')
            meanLeg.SetBorderSize(0)
            meanLeg.SetNColumns(2)
            for i, gr in enumerate(meanGraphs):
                entry = meanLeg.AddEntry(gr.GetName(), 'Scan '+scans[i], 'lep')
                entry.SetMarkerStyle(3+19*i)
                entry.SetMarkerColor(2+i)
            meanLeg.Draw()
            drawCMS()
            meanCanvas.Modified()
            meanCanvas.Update()
            meanCanvas.SaveAs('summaryPlots/'+meanCanvas.GetName()+'.pdf')
            meanCanvas.SaveAs('summaryPlots/'+meanCanvas.GetName()+'.C')
            rmsCanvas = TCanvas('canvasRms'+coord+'ScanBX'+bx, '', 1000, 600)
            rmsMulti.Draw('APE')
            rmsCanvas.Update()
            rmsMulti.GetYaxis().SetTitle('luminous region '+other.lower()+'-width [cm]')
            rmsMulti.GetXaxis().SetTitle('scan point')
            rmsMulti.GetYaxis().SetTitleOffset(1.3)
            rmsMulti.GetYaxis().SetRangeUser(-0.01, 0.04)
            rmsMulti.GetYaxis().SetTickLength(0.02)
            rmsLeg = TLegend(0.4, 0.13, 0.6, 0.28)
            rmsLeg.SetHeader(coord+' Scan (BCID '+bx+')')
            rmsLeg.SetBorderSize(0)
            rmsLeg.SetNColumns(2)
            for i, gr in enumerate(rmsGraphs):
                entry = rmsLeg.AddEntry(gr.GetName(), 'Scan '+scans[i], 'lep')
                entry.SetMarkerStyle(3+19*i)
                entry.SetMarkerColor(2+i)
            rmsLeg.Draw()
            drawCMS()
            rmsCanvas.Modified()
            rmsCanvas.Update()
            rmsCanvas.SaveAs('summaryPlots/'+rmsCanvas.GetName()+'.pdf')
            rmsCanvas.SaveAs('summaryPlots/'+rmsCanvas.GetName()+'.C')

if __name__ == '__main__':
    plotOtherCoordinate(filename, bunchcrossings, biscans)
