def loadJson(filename):
    from json import load
    result = {}
    nextfile = filename
    while nextfile:
        if 'include' in result:
            result.pop('include')
        with open(nextfile) as f:
            json = load(f)
        json.update(result)
        result = json
        if 'include' in result and result['include']:
            nextfile = result['include']
        else:
            nextfile = False
    return result
