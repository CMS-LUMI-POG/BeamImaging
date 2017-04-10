from DoubleGauss import DoubleGauss
from math import pi, exp
from re import match
from ROOT import RooRealVar, RooFormulaVar, RooArgList, TF2

class TripleGauss(DoubleGauss):
    """Beam shapes modelled with triple Gaussian"""

    Shortname = 'TG'
    Templates = ('TripleGauss_V1', 'TripleGauss_V2')

    super = False

    def __init__(self, super=False):
        DoubleGauss.__init__(self, super=super)
        sXY = ('x', 'y')
        s12 = ('1', '2')
        diffW = {name: RooRealVar(name, name, 0.01, 1.7) for name in \
                 [c+'WidthW'+i+'Diff' for c in sXY for i in s12]}
        widthW = {c+'WidthW'+i: RooFormulaVar(c+'WidthW'+i, c+'WidthN'+i+'+'+ \
                  c+'WidthM'+i+'Diff'+'+'+c+'WidthW'+i+'Diff', \
                  RooArgList(diffW[c+'WidthW'+i+'Diff'], \
                  self.Parameter[c+'WidthM'+i+'Diff'],
                  self.Parameter[c+'WidthN'+i])) for c in sXY for i in s12}
        rhoW = {name: RooRealVar(name, name, -0.48, 0.48) for name in \
                ['rhoW'+i for i in s12]}
        if super:
            self.super = True
            extra = {name: RooRealVar(name, name, 0.0, 1.0) for name in \
                     ['w'+i+'MFraction' for i in s12]}
            nf = [('M', lambda i: 'w'+i+'MFraction*(1.0-w'+i+'N)'), \
                  ('W', lambda i: '(1.0-w'+i+'MFraction)*(1.0-w'+i+'N)')]
            w = {'w'+i+n: RooFormulaVar('w'+i+n, f(i), \
                 RooArgList(extra['w'+i+'MFraction'], self.Parameter['w'+i+'N'])) \
                 for i in s12 for (n, f) in nf}
        else:
            extra = {name: RooRealVar(name, name, 0.0, 0.5*pi) for name in \
                     [a+i for a in ('theta', 'phi') for i in s12]}
            nf = [('N', lambda i: 'sin(theta'+i+')**2*cos(phi'+i+')**2'), \
                  ('M', lambda i: 'sin(theta'+i+')**2*sin(phi'+i+')**2'), \
                  ('W', lambda i: 'cos(theta'+i+')**2')]
            w = {'w'+i+n: RooFormulaVar('w'+i+n, f(i), \
                 RooArgList(extra['theta'+i], extra['phi'+i])) for i in s12 \
                 for (n, f) in nf}
        for d in (diffW, widthW, rhoW, extra, w):
            self.Parameter.update(d)

    def getModelFunctions(self):
        from ROOT import TripleGauss_V1, TripleGauss_V2
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
        return [locals()['TripleGauss_V'+whi(c)](name(i, c), name(i, c), \
                self.xVar(), self.yVar(), self.Parameter['x0'+i+whi(c)], \
                self.Parameter['y0'+i+whi(c)], self.Parameter['w'+i+'N'], \
                self.Parameter['w'+i+'M'], self.Parameter['rhoN'+i], \
                self.Parameter['xWidthN'+i], self.Parameter['yWidthN'+i], \
                self.Parameter['rhoM'+i], self.Parameter['xWidthM'+i], \
                self.Parameter['yWidthM'+i], self.Parameter['rhoW'+i], \
                self.Parameter['xWidthW'+i], self.Parameter['yWidthW'+i], \
                self.Parameter['w'+oth(i)+'N'], self.Parameter['w'+oth(i)+'M'], \
                self.Parameter[oth(c)+'WidthN'+oth(i)], \
                self.Parameter[oth(c)+'WidthM'+oth(i)], \
                self.Parameter[oth(c)+'WidthW'+oth(i)], \
                self.Parameter['vtxRes']) for i,c in \
                [('1', 'x'), ('1', 'y'), ('2', 'x'), ('2', 'y')]]

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
        wM1 = par[18]
        wM2 = par[19]
        xWidthW1 = par[20]
        yWidthW1 = par[21]
        xWidthW2 = par[22]
        yWidthW2 = par[23]
        rhoW1 = par[24]
        rhoW2 = par[25]
        try:
            beamN1 = 0.5 / (pi * (1.0 - rhoN1 ** 2) ** 0.5 * abs(xWidthN1) * \
                     abs(yWidthN1)) * exp(-0.5 / (1.0 - rhoN1 ** 2) * (((xx - x01) \
                     / xWidthN1) ** 2 + ((yy -y01) / yWidthN1) ** 2 - 2.0 * rhoN1 \
                     * (xx - x01) * (yy - y01) / (xWidthN1 * yWidthN1)))
            beamM1 = 0.5 / (pi * (1.0 - rhoM1 ** 2) ** 0.5 * abs(xWidthM1) * \
                     abs(yWidthM1)) * exp(-0.5 / (1.0 - rhoM1 ** 2) * (((xx - x01) \
                     / xWidthM1) ** 2 + ((yy -y01) / yWidthM1) ** 2 - 2.0 * rhoM1 \
                     * (xx - x01) * (yy - y01) / (xWidthM1 * yWidthM1)))
            beamW1 = 0.5 / (pi * (1.0 - rhoW1 ** 2) ** 0.5 * abs(xWidthW1) * \
                     abs(yWidthW1)) * exp(-0.5 / (1.0 - rhoW1 ** 2) * (((xx - x01) \
                     / xWidthW1) ** 2 + ((yy -y01) / yWidthW1) ** 2 - 2.0 * rhoW1 \
                     * (xx - x01) * (yy - y01) / (xWidthW1 * yWidthW1)))
            beamN2 = 0.5 / (pi * (1.0 - rhoN1 ** 2) ** 0.5 * abs(xWidthN1) * \
                     abs(yWidthN1)) * exp(-0.5 / (1.0 - rhoN1 ** 2) * (((xx - x01) \
                     / xWidthN1) ** 2 + ((yy -y01) / yWidthN1) ** 2 - 2.0 * rhoN1 \
                     * (xx - x01) * (yy - y01) / (xWidthN1 * yWidthN1)))
            beamM2 = 0.5 / (pi * (1.0 - rhoM1 ** 2) ** 0.5 * abs(xWidthM1) * \
                     abs(yWidthM1)) * exp(-0.5 / (1.0 - rhoM1 ** 2) * (((xx - x01) \
                     / xWidthM1) ** 2 + ((yy -y01) / yWidthM1) ** 2 - 2.0 * rhoM1 \
                     * (xx - x01) * (yy - y01) / (xWidthM1 * yWidthM1)))
            beamW2 = 0.5 / (pi * (1.0 - rhoW1 ** 2) ** 0.5 * abs(xWidthW1) * \
                     abs(yWidthW1)) * exp(-0.5 / (1.0 - rhoW1 ** 2) * (((xx - x01) \
                     / xWidthW1) ** 2 + ((yy -y01) / yWidthW1) ** 2 - 2.0 * rhoW1 \
                     * (xx - x01) * (yy - y01) / (xWidthW1 * yWidthW1)))
            return (wN1 * beamN1 + wM1 * beamM1 + (1.0 - wN1 - wM1) * beamW1) * \
                   (wN2 * beamN2 + wM1 * beamM2 + (1.0 - wN2 - wM2) * beamW2)
        except:
            return -1.0

    def assignToOverlap(self, multBeam, random=False, uncorrelated=False):
        multBeam = DoubleGauss.assignToOverlap(self, multBeam, random)
        parNames = ['w1M', 'w2M', 'xWidthW', 'yWidthW1', 'xWidthW2', \
                    'yWidthW2', 'rhoW1', 'rhoW2']
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
            parValue['rhoW1'] = 0.0
            parValue['rhoW2'] = 0.0
        for i, name in enumerate(parNames, start=18):
            multBeam.SetParameter(i, parValue[name])
        return multBeam

    def funcOverlap(self):
        multBeam = TF2('multBeam', self.calcOverlap, -30.0, 30.0, -30.0, 30.0, 26)
        multBeam = self.assignToOverlap(multBeam)
        return multBeam

    def computeError(self, var):
        if match('^[xy]WidthW[12]$', var):
            c = var[0]
            i = var[7]
            error1 = self.computeError(c+'WidthN'+i)
            error2 = self.computeError(c+'WidthM'+i+'Diff')
            error3 = self.computeError(c+'WidthW'+i+'Diff')
            return (error1 ** 2 + error2 ** 2 + error3 ** 2) ** 0.5
        elif not self.super and match('^w[12][NMW]$', var):
            i = var[1]
            n = var[2]
            vT = self.computeValue('theta'+i)
            vP = self.computeValue('phi'+i)
            eT = self.computeError('theta'+i)
            eP = self.computeError('phi'+i)
            if n == 'N':
                return 2.0 * sin(vT) * cos(vP) * ((eT * cos(vT) * cos(vP)) \
                       ** 2 + (eP * sin(vT) * sin(vP)) ** 2) ** 0.5
            elif n == 'M':
                return 2.0 * sin(vT) * sin(vP) * ((eT * cos(vT) * sin(vP)) \
                       ** 2 + (eP * sin(vT) * cos(vP)) ** 2) ** 0.5
            else:
                return 2.0 * sin(vT) * cos(vT) * eT
        elif self.super and match('^w[12][MW]$', var):
            i = var[1]
            error1 = (1.0 - self.computeValue('w'+i+'N')) * \
                     self.computeError('w'+i+'MFraction')
            f = {'M': lambda a: a, 'W': lambda a: 1.0-a}[var[2]]
            error2 = f(self.computeValue('w'+i+'MFraction')) * \
                     self.computeError('w'+i+'N')
            return (error1 ** 2 + error2 ** 2) ** 0.5
        else:
            return DoubleGauss.computeError(self, var)


class SuperDoubleGauss(TripleGauss):
    """Beam shapes modelled with Super Double Gaussian"""

    Shortname = 'SupDG'

    def __init__(self):
        TripleGauss.__init__(self, super=True)
