from ROOT import TFile

bunchcrossings = ('41', '281', '872', '1783', '2063')
beamshapes = ('SG', 'DG', 'SupG', 'TG', 'SupDG')

def gatherFromToys(crossings, shapes):
    overDiff = {}
    for shape in shapes:
        overDiff[shape] = {}
        for bx in crossings:
            f = TFile.Open('overlapDiff_TOYS_2016_'+bx+shape+'.root')
            if f:
                hist = f.Get('overDiff')
                overDiff[shape][bx] = hist.GetMean()
            else:
                overDiff[shape][bx] = False
    return overDiff

if __name__ == '__main__':
    overDiff = gatherFromToys(bunchcrossings, beamshapes)
    print
    print
    print '',
    for bx in bunchcrossings:
        print ';', bx,
    print
    for shape in beamshapes:
        print shape,
        for bx in bunchcrossings:
            if overDiff[shape][bx]:
                print ';', overDiff[shape][bx],
            else:
                print ';', '',
        print
    print
    print
