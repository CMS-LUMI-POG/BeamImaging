from SingleGauss import SingleGauss
from math import pi, exp
from re import match
from ROOT import RooRealVar, RooFormulaVar, RooArgList, TF2

class DoubleGauss(SingleGauss):
    """Beam shapes modelled with double Gaussian"""

    Shortname = 'DG'
    Templates = ('MyPdfV3_Ext', 'MyPdfV4_Ext')

    def __init__(self, super=False):
        SingleGauss.__init__(self, uncorrelated=False)
        sXY = ('x', 'y')
        s12 = ('1', '2')
        diffM = {name: RooRealVar(name, name, 0.01, 1.7) for name in \
                 [c+'WidthM'+i+'Diff' for c in sXY for i in s12]}
        widthM = {c+'WidthM'+i: RooFormulaVar(c+'WidthM'+i, c+'WidthN'+i+'+'+ \
                  c+'WidthM'+i+'Diff', RooArgList(diffM[c+'WidthM'+i+'Diff'], \
                  self.Parameter[c+'WidthN'+i])) for c in sXY for i in s12}
        rhoM = {name: RooRealVar(name, name, -0.48, 0.48) for name in \
                ['rhoM'+i for i in s12]}
        for i in s12:
            wN = self.Parameter['w'+i+'N']
            wN.setConstant(False)
            if super:
                wN.setRange(-1.0, 0.0)
                wN.setVal(-0.5)
            else:
                wN.setVal(0.5)
            self.Parameter['w'+i+'N'] = wN
        wM = {'w'+i+'M': RooFormulaVar('w'+i+'M', '1.0-w'+i+'N', \
              RooArgList(self.Parameter['w'+i+'N'])) for i in s12}
        for d in (diffM, widthM, rhoM, wM):
            self.Parameter.update(d)

    def getModelFunctions(self):
        from ROOT import MyPdfV3_Ext as DoubleGauss_V1, \
                         MyPdfV4_Ext as DoubleGauss_V2
        def name(i, c):
            return 'beam' + i + 'RestVerticesUnfold_' + c.upper() + 'Scan'
        def oth(s):
            if s.isalpha():
                if s == 'x':
                    return 'y'
                else:
                    return 'x'
            else:
                if s == '1':
                    return '2'
                else:
                    return '1'
        def whi(s):
            if s == 'x':
                return '1'
            else:
                return '2'
        return [locals()['DoubleGauss_V'+whi(c)](name(i, c), name(i, c), \
                self.xVar(), self.yVar(), self.Parameter['x0'+i+whi(c)], \
                self.Parameter['y0'+i+whi(c)], self.Parameter['w'+i+'N'], \
                self.Parameter['rhoN'+i], self.Parameter['xWidthN'+i], \
                self.Parameter['yWidthN'+i], self.Parameter['rhoM'+i], \
                self.Parameter['xWidthM'+i], self.Parameter['yWidthM'+i], \
                self.Parameter['w'+oth(i)+'N'], self.Parameter[oth(c)+'WidthN'+oth(i)], \
                self.Parameter[oth(c)+'WidthM'+oth(i)], self.Parameter['vtxRes']) \
                for i,c in [('1', 'x'), ('1', 'y'), ('2', 'x'), ('2', 'y')]]

    def calcOverlap(self, x, par):
        xx = x[0]
        yy = x[1]
        x01 = par[0]
        y01 = par[1]
        x02 = par[2]
        y02 = par[3]
        xWidthN1 = par[4]
        yWidthN1 = par[5]
        xWidthN2 = par[6]
        yWidthN2 = par[7]
        rhoN1 = par[8]
        rhoN2 = par[9]
        wN1 = par[10]
        wN2 = par[11]
        xWidthM1 = par[12]
        yWidthM1 = par[13]
        xWidthM2 = par[14]
        yWidthM2 = par[15]
        rhoM1 = par[16]
        rhoM2 = par[17]
        try:
            beamN1 = 0.5 / (pi * (1.0 - rhoN1 ** 2) ** 0.5 * abs(xWidthN1) * \
                     abs(yWidthN1)) * exp(-0.5 / (1.0 - rhoN1 ** 2) * (((xx - x01) \
                     / xWidthN1) ** 2 + ((yy -y01) / yWidthN1) ** 2 - 2.0 * rhoN1 \
                     * (xx - x01) * (yy - y01) / (xWidthN1 * yWidthN1)))
            beamM1 = 0.5 / (pi * (1.0 - rhoM1 ** 2) ** 0.5 * abs(xWidthM1) * \
                     abs(yWidthM1)) * exp(-0.5 / (1.0 - rhoM1 ** 2) * (((xx - x01) \
                     / xWidthM1) ** 2 + ((yy -y01) / yWidthM1) ** 2 - 2.0 * rhoM1 \
                     * (xx - x01) * (yy - y01) / (xWidthM1 * yWidthM1)))
            beamN2 = 0.5 / (pi * (1.0 - rhoN1 ** 2) ** 0.5 * abs(xWidthN1) * \
                     abs(yWidthN1)) * exp(-0.5 / (1.0 - rhoN1 ** 2) * (((xx - x01) \
                     / xWidthN1) ** 2 + ((yy -y01) / yWidthN1) ** 2 - 2.0 * rhoN1 \
                     * (xx - x01) * (yy - y01) / (xWidthN1 * yWidthN1)))
            beamM2 = 0.5 / (pi * (1.0 - rhoM1 ** 2) ** 0.5 * abs(xWidthM1) * \
                     abs(yWidthM1)) * exp(-0.5 / (1.0 - rhoM1 ** 2) * (((xx - x01) \
                     / xWidthM1) ** 2 + ((yy -y01) / yWidthM1) ** 2 - 2.0 * rhoM1 \
                     * (xx - x01) * (yy - y01) / (xWidthM1 * yWidthM1)))
            return (wN1 * beamN1 + (1.0 - wN1) * beamM1) * \
                   (wN2 * beamN2 + (1.0 - wN2) * beamM2)
        except:
            return -1.0

    def assignToOverlap(self, multBeam, random=False, uncorrelated=False):
        multBeam = SingleGauss.assignToOverlap(self, multBeam, random)
        parNames = ['w1N', 'w2N', 'xWidthM1', 'yWidthM1', 'xWidthM2', \
                    'yWidthM2', 'rhoM1', 'rhoM2']
        parValue = {}
        for name in parNames:
            parValue[name] = self.computeValue(name)
        if random:
            for name in parNames:
                parValue[name] = random.Gaus(parValue[name], \
                                             self.computeError(name))
                if 'rho' in name:
                    if parValue[name] < -1.0:
                        parValue[name] = -1.0
                    elif parValue[name] > 1.0:
                        parValue[name] = 1.0
        elif uncorrelated:
            parValue['rhoM1'] = 0.0
            parValue['rhoM2'] = 0.0
        for i, name in enumerate(parNames, start=10):
            multBeam.SetParameter(i, parValue[name])
        return multBeam

    def funcOverlap(self):
        multBeam = TF2('multBeam', self.calcOverlap, -30.0, 30.0, -30.0, 30.0, 18)
        multBeam = self.assignToOverlap(multBeam)
        return multBeam

    def computeError(self, var):
        if match('^[xy]WidthM[12]$', var):
            c = var[0]
            i = var[7]
            error1 = self.computeError(c+'WidthN'+i)
            error2 = self.computeError(c+'WidthM'+i+'Diff')
            return (error1 ** 2 + error2 ** 2) ** 0.5
        elif match('^w[12]M$', var):
            return self.computeError('w'+var[1]+'N')
        else:
            return SingleGauss.computeError(self, var)


class SuperGauss(DoubleGauss):
    """Beam shapes modelled with Super Gaussian"""

    Shortname = 'SupG'

    def __init__(self):
        DoubleGauss.__init__(self, super=True)
