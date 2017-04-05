def prepareDataFile(listfile, times, minTrk, nbins, bcids, scaling, offsetx, \
                    offsety, outputpath, outputname, preliminary=False):
    from array import array
    from itertools import product
    from os import mkdir
    from os.path import exists
    from ROOT import TChain, TH1F, TH2F, TFile

    class PrelIterator:
        def __init__(self, f):
            self.iterator = iter(f)
        def __iter__(self):
            return self
        def next(self):
            line = self.iterator.next()
            if 'ZeroBias1' in line:
                return line
            else:
                raise StopIteration

    chain = {name: TChain('lumi/tree') for name in listfile}
    if preliminary:
        iterator = lambda f: PrelIterator(f)
        stepsize = 10
    else:
        iterator = lambda f: f
        stepsize = 1
    for name, filename in listfile.iteritems():
        with open(filename) as f:
            i = iterator(f)
            for line in i:
                chain[name].Add(line.strip())

    values = {name: array(t, v) for (name, t, v) in [ \
              ('bunchCrossing', 'i', [0]), \
              ('nVtx', 'i', [0]), \
              ('vtx_nTrk', 'i', [0]*200), \
              ('vtx_x', 'f', [0.0]*200), \
              ('vtx_y', 'f', [0.0]*200), \
              ('vtx_xError', 'f', [0.0]*200), \
              ('vtx_yError', 'f', [0.0]*200), \
              ('vtx_isGood', 'b', [0]*200), \
              ('vtx_isFake', 'b', [0]*200), \
              ('timeStamp_begin', 'I', [0])]}
    for ch in chain.itervalues():
        for field in values:
            ch.SetBranchAddress(field, values[field])

    histos = {c+'Error': TH1F(c+'Error_hist', c+'Error_hist', 1000, 0.0, 0.1) \
              for c in ['x', 'y']}
    def moveName(i, c, b):
        return 'Beam' + i + 'Move' + c + '_bunch' + b + 'Add'
    histos.update({name: TH2F(name, name, nbins, -10.0, 10.0, nbins, -10.0, \
                   10.0) for name in [moveName(i, c, b) for i in ['1', '2'] \
                   for c in ['X', 'Y'] for b in bcids]})

    percent = 0
    for name, ch in chain.iteritems():
        nEntries = ch.GetEntries()
        nDisplay = ((nEntries / 25) / stepsize + 1) * stepsize
        print '<<< Start to process {} with {} events'.format(name, \
              nEntries/stepsize)
        for i in range(0, nEntries, stepsize):
            if i % nDisplay == 0:
                print '<<< Now at event {} of {} ({:.0f}%)'.format(i/stepsize, \
                      name, i*25.0/nEntries+percent)
            ch.GetEntry(i)
            if values['nVtx'][0] <= 0:
                continue
            if not str(values['bunchCrossing'][0]) in bcids:
                continue
            for j in range(values['nVtx'][0]):
                if not values['vtx_isGood'][j]:
                    continue
                if values['vtx_isFake'][j]:
                    continue
                if values['vtx_nTrk'][j] <= minTrk:
                    continue
                histos['xError'].Fill(values['vtx_xError'][j])
                histos['yError'].Fill(values['vtx_yError'][j])
                xVtx = values['vtx_x'][j] / scaling + offsetx
                yVtx = values['vtx_y'][j] / scaling + offsety
                for begin, end in times[name]:
                    if values['timeStamp_begin'][0] < begin:
                        continue
                    if values['timeStamp_begin'][0] > end:
                        continue
                    histos[moveName(name[0], name[1], \
                           str(values['bunchCrossing'][0]))].Fill(xVtx, yVtx)
        percent += 25

    if not exists(outputpath):
        mkdir(outputpath)
    outputfile = TFile(outputpath+'/'+outputname+'.root', 'RECREATE')
    for hist in histos.itervalues():
        print '<<< Write {} with {:.0f} entries'.format(hist.GetName(), \
              hist.GetEntries())
        hist.Write()
    outputfile.Write()
    outputfile.Close()

def main():
    from sys import argv as __ARGV__
    __ARGV__.append('-b')

    from tools import loadJson
    from argparse import ArgumentParser
    from re import match

    parser = ArgumentParser(description='Extract Beam Imaging data from ROOT '+ \
                            'files')
    parser.add_argument('-b', action='store_true', help='enable batch mode')
    parser.add_argument('json', nargs=1, help='specify JSON file containing '+ \
                        'config information')
    parser.add_argument('-p', '--preliminary', action='store_true', help='run '+ \
                        'only on 1/80 of data')
    args = parser.parse_args()

    json = loadJson(args.json[0])
    listfile = {name: 'filelist/'+str(json['prefix'])+'_'+name+'.txt' for \
                name in ['1X', '1Y', '2X', '2Y']}
    times = {name[4:6]: [(beg, beg+25) for beg in json[name]] for name \
             in json if match('^scan[12][XY]MoveBegin$', name)}
    minTrk = int(json['minTrk'])
    nbins = int(json['nbins'])
    bcids = [str(bx) for bx in json['bunchCrossings']]
    scaling = float(json['scaling'])
    offsetx = float(json['offsetx'])
    offsety = float(json['offsety'])
    outputpath = str(json['datapath'])
    outputname = str(json['prefix'] )+ '_' + str(json['suffix'])
    preliminary = args.preliminary
    if preliminary:
        outputname += '_prel'

    prepareDataFile(listfile, times, minTrk, nbins, bcids, scaling, offsetx, \
                        offsety, outputpath, outputname, preliminary)

if __name__ == '__main__':
    main()
