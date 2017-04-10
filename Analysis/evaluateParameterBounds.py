def evaluateParameterBounds(prefix, suffix, models, vtxres, bcids):
    from os.path import exists
    from json import load
    from ROOT import TFile

    nPassed = 0
    nBound = 0
    nInvalid = 0
    atBound = []
    notValid = []

    for model in models:
        parameterfile = 'shapes/config_'+model+'.json'
        if not exists(parameterfile):
            continue
        with open(parameterfile) as f:
            parameterlist = load(f)
        parameters = {}
        for par in parameterlist:
            if len(parameterlist[par]) >= 2:
                parameters[par] = (parameterlist[par][0], parameterlist[par][1])
        for bcid in bcids:
            names = [pre+'_'+model+'_'+bcid+'_'+suf+'_'+vr for (pre, suf) in \
                     zip(prefix, suffix) for vr in vtxres]
            for name in names:
                try:
                    print '<<< Check file:', name
                    f = TFile('results/'+name+'.root')
                    passed = True
                    for par, (low, high) in parameters.iteritems():
                        val = f.Get('h_'+par).GetMean()
                        if val < 1.02 * low:
                            print '<<<! Parameter', par, 'close to lower border:', \
                                  val, '({})'.format(low)
                            passed = False
                        if val > 0.98 * high:
                            print '<<<! Parameter', par, 'close to upper border:', \
                                  val, '({})'.format(high)
                            passed = False
                except:
                    print '<<<! File does not exist'
                    nInvalid += 1
                    notValid.append(name)
                    continue
                if passed:
                    nPassed += 1
                else:
                    nBound += 1
                    atBound.append(name)

    print '<<<', nPassed, 'files passed parameter boundary checks'
    if nBound > 0:
        print '<<<', nBound, 'files have parameters close to border:'
        for name in atBound:
            print '      ', name
    if nInvalid > 0:
        print '<<<', nInvalid, 'files do not exist:'
        for name in notValid:
            print '      ', name

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
    parser.add_argument('-m', '--model', required=True, nargs='+', choices= \
                        ['SingleGauss', 'SingleGaussUncorrelated', \
                        'DoubleGauss', 'SuperGauss', 'TripleGauss', \
                        'SuperDoubleGauss'], help='specify fit model')
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
    models = [locals()[model].Shortname for model in args.model]
    vtxres = [str(v) for v in args.vtxres]
    bcids = [str(bcid) for bcid in args.bcid]

    evaluateParameterBounds(prefix, suffix, models, vtxres, bcids)

if __name__ == '__main__':
    main()
