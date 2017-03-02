from ROOT import TFile

bunchcrossings = ('41', '281', '872', '1783', '2063')
beamshapes = ('SG', 'DG', 'SupG', 'TG', 'SupDG')

def computeChiSquares(crossings, shapes):
    components = ('X1', 'Y1', 'X2', 'Y2')
    chiSq = {}
    dof = {}
    for shape in shapes:
        chiSq[shape] = {}
        dof[shape] = {}
        for bx in crossings:
            f = TFile.Open('DataAnalysisBunch'+bx+shape+'_new_StronRescale.root')
            if f:
                c = 0.0
                d = 0
                for comp in components:
                    res = f.Get('res'+comp)
                    for x in range(res.GetXaxis().GetNbins()):
                        for y in range(res.GetYaxis().GetNbins()):
                            c += res.GetBinContent(x,y) ** 2
                            if res.GetBinContent(x,y):
                                d += 1
                chiSq[shape][bx] = c
                dof[shape][bx] = d
            else:
                chiSq[shape][bx] = False
                dof[shape][bx] = False
    return chiSq, dof

if __name__ == '__main__':
    chiSq, dof = computeChiSquares(bunchcrossings, beamshapes)
    print
    print
    print '',
    for bx in bunchcrossings:
        print ';', bx,
    print
    for shape in beamshapes:
        print shape,
        for bx in bunchcrossings:
            if chiSq[shape][bx] and dof[shape][bx]:
                print ';', chiSq[shape][bx]/dof[shape][bx],
            else:
                print ';', '',
        print
    print
    print
