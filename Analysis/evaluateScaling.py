def evaluateScaling(bcids, scaling, offsetx, offsety, filepath, filename, \
                    preliminary=False):
    from ROOT import TFile

    def moveName(i, c, b):
        return 'Beam' + i + 'Move' + c + '_bunch' + b + 'Add'
    histos = [moveName(i, c, b) for i in ['1', '2'] for c in ['X', 'Y'] for \
              b in bcids]
    meanx = {name: False for name in histos}
    meany = {name: False for name in histos}
    rmsx = {name: False for name in histos}
    rmsy = {name: False for name in histos}

    rootfile = TFile(filepath+'/'+filename+'.root')
    for name in histos:
        hist = rootfile.Get(name)
        meanx[name] = hist.GetMean(1)
        meany[name] = hist.GetMean(2)
        rmsx[name] = hist.GetRMS(1)
        rmsy[name] = hist.GetRMS(2)
    rootfile.Close()

    def addMean(d):
        d['mean'] = sum(d.values()) / len(d.values())
        return d
    meanx = addMean(meanx)
    meany = addMean(meany)
    rmsx = addMean(rmsx)
    rmsy = addMean(rmsy)

    print
    print '{:24} {}'.format('scaling:', scaling)
    print '{:24} {}'.format('offsetx:', offsetx)
    print '{:24} {}'.format('offsety:', offsety)
    print
    for name in histos + ['mean']:
        print '{:24} {:8.4f} {:8.4f} {:8.4f} {:8.4f}'.format(name+':', \
              meanx[name], meany[name], rmsx[name], rmsy[name])
    print

def main():
    from sys import argv as __ARGV__
    __ARGV__.append('-b')

    from tools import loadJson
    from argparse import ArgumentParser

    parser = ArgumentParser(description='List mean and RMS of histograms '+ \
                            'in data files')
    parser.add_argument('-b', action='store_true', help='enable batch mode')
    parser.add_argument('json', nargs=1, help='specify JSON file containing '+ \
                        'config information')
    parser.add_argument('-p', '--preliminary', action='store_true', \
                        help='use data files with only 1/80 of data')
    args = parser.parse_args()

    json = loadJson(args.json[0])
    bcids = [str(bx) for bx in json['bunchCrossings']]
    scaling = float(json['scaling'])
    offsetx = float(json['offsetx'])
    offsety = float(json['offsety'])
    filepath = str(json['datapath'])
    filename = str(json['prefix'] )+ '_' + str(json['suffix'])
    preliminary = args.preliminary
    if preliminary:
        filename += '_prel'

    evaluateScaling(bcids, scaling, offsetx, offsety, filepath, filename, \
                    preliminary)

if __name__ == '__main__':
    main()
