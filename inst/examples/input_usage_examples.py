## informal examples that test input format 

import climextremes, numpy, pandas

climextremes.normalize(numpy.array((0.,2.)))
climextremes.normalize(numpy.array([0,2]))
climextremes.normalize(numpy.array([[0,2],[1,3]]))

climextremes.remove_runs(numpy.array((1,2,3,4,5,6)), numpy.array((1,3,4,5,8,9)))
climextremes.remove_runs(numpy.array((1,2,3,4,5,6)), numpy.array((1,3,4,5,8,9)), False)

index = numpy.array((1,10,11,15,19,21,45,51,53))
y = numpy.array(range(9))
climextremes.screen_within_block(y, index)

climextremes.screen_within_block(y, index, 5)


Fort = climextremes.Fort
firstYr = min(Fort.year)
yrs = numpy.array(range(int(firstYr), int(max(Fort.year)+1)))
nYrs = len(yrs)
yrsToPred = numpy.array([min(Fort.year), max(Fort.year)])
threshold = 0.395
FortExc = Fort[Fort.Prec > threshold]
covByBlock = numpy.c_[yrs, numpy.random.rand(nYrs)]
covByBlockp = pandas.DataFrame({'xA': covByBlock[:,1], 'xInt': yrs})

result = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False) 
# x
result2 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlockp, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([2,1]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[0, min(Fort.year)], [0, max(Fort.year)]]), bootSE = False) 
# threshold
result3 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = numpy.repeat(threshold, nYrs), blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False) 
# locationFun
result4 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlockp, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array(('xInt','xA')), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = pandas.DataFrame({'xA': [0,0], 'xInt': [min(Fort.year), max(Fort.year)]}), bootSE = False) 
# locationFun with xContrast
result4a = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlockp, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array(('xInt','xA')), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = pandas.DataFrame({'xA': [0,0], 'xInt': [min(Fort.year), max(Fort.year)]}), xContrast = pandas.DataFrame({'xA': [0,0], 'xInt': [max(Fort.year), min(Fort.year)]}), bootSE = False) 
# weights
result5 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, weights = numpy.repeat(1, nYrs), threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False) 
# replicateIndex
replicateIndex = numpy.repeat(1, 1061*2)
ind = numpy.linspace(1,1061*2-1,1061 )
ind = ind.astype('int')
replicateIndex[ind] = 2
# not sure why 'result6' a bit different than 'result'
result6 = climextremes.fit_pot(numpy.array(numpy.repeat(FortExc.Prec, 2)), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.repeat(numpy.array(FortExc.year), 2), nReplicates = 2, replicateIndex = replicateIndex, locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False)
# 0 variability because resampling gives same data
result6 = climextremes.fit_pot(numpy.array(numpy.repeat(FortExc.Prec, 2)), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.repeat(numpy.array(FortExc.year), 2), nReplicates = 2, replicateIndex = replicateIndex, locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = True, bootControl = {'by':'replicate'}) 
# proportionMissing
result7 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, proportionMissing = numpy.repeat(0, nYrs), threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False) 
# getFit
result8 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, getFit = True, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False) 
# declustering
result9 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False, declustering = 'noruns', index = numpy.array(FortExc.obs)) 
result10 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False, declustering = 5, index = numpy.array(FortExc.obs)) 

result11 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = True) 
result12 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = True, bootControl = {'seed':1, 'n':250}) 
result13 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = True, bootControl = {'seed':0, 'n':250}) 
# optimArgs
result14 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False, optimArgs = {'method':'BFGS'}) 
result15 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False, optimArgs = {'method':'BFGS'}, optimControl = {'maxit': 1000}) 
# initial
result16 = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False, initial = {'location': numpy.array([1, 0, 0]), 'scale': 0.5, 'shape': 0.2})


Fort = climextremes.Fort
FortMax = Fort.groupby('year').max()[['Prec']]
FortMax.reset_index(inplace=True)

# stationary fit
result = climextremes.fit_gev(numpy.array(FortMax.Prec), returnPeriod = 20, returnValue = 3.5, getParams = True, bootSE = False)

tmp = numpy.array(FortMax.Prec)
tmp[3] = numpy.nan
result2 = climextremes.fit_gev(tmp, returnPeriod = 20, returnValue = 3.5, getParams = True, bootSE = False)
tmp = numpy.array(FortMax.Prec)
tmp[3] = -99999
result3 = climextremes.fit_gev(tmp, returnPeriod = 20, returnValue = 3.5, getParams = True, bootSE = False, missingFlag = -99999)


result = climextremes.calc_riskRatio_binom(numpy.array((0, 0)), numpy.array((400, 400)), ciLevel = 0.95, ciType = 'lrt')

result = climextremes.calc_riskRatio_binom(numpy.array((3, 3)), numpy.array((400, 400)), ciLevel = 0.95, ciType = 'lrt')
result5 = climextremes.calc_riskRatio_binom(numpy.array((3, 3)), numpy.array((400, 400)), ciLevel = 0.95, ciType = 'lrt', lrtControl = {'bounds': numpy.array((.0001,1000))})
result1 = climextremes.calc_riskRatio_binom(numpy.array((3, 3)), numpy.array((400, 400)), ciLevel = 0.95, ciType = 'boot_perc')
result2 = climextremes.calc_riskRatio_binom(numpy.array((3, 3)), numpy.array((400, 400)), ciLevel = 0.95, ciType = 'boot_perc', bootControl = {'n':250,'seed':1})
result3 = climextremes.calc_riskRatio_binom(numpy.array((3, 3)), numpy.array((400, 400)), ciLevel = 0.95, ciType = 'boot_perc', bootControl = {'n':250,'seed': 2})

import climextremes, pandas, numpy
Fort = climextremes.Fort
earlyPeriod = numpy.array((1900, 1930))
earlyYears = numpy.array(range(earlyPeriod[0], earlyPeriod[1]))
latePeriod = numpy.array((1970, 2000))
lateYears = numpy.array(range(latePeriod[0], latePeriod[1]))

FortMax = Fort.groupby('year').max()[['Prec']]
FortMax.reset_index(inplace=True)

y1 = FortMax.Prec[FortMax.year < earlyPeriod[1]]
y2 = FortMax.Prec[FortMax.year >= latePeriod[0]]

# contrast late period with early period, assuming a nonstationary fit
# within each time period and finding RR based on midpoint of each period
result = climextremes.calc_riskRatio_gev(
    returnValue = 3,
    y1 = numpy.array(y1), y2 = numpy.array(y2),
    x1 = earlyYears, x2 = lateYears,
    locationFun1 = 1, locationFun2 = 1,
    xNew1 = earlyYears.mean(), xNew2 = lateYears.mean(),
    ciType = 'lrt')
result2 = climextremes.calc_riskRatio_gev(
    returnValue = 3,
    y1 = numpy.array(y1), y2 = numpy.array(y2),
    x1 = earlyYears, x2 = lateYears,
    weights1 = numpy.repeat(1, len(y1)), weights2 = numpy.repeat(1,len(y2)),
    locationFun1 = 1, locationFun2 = 1,
    xNew1 = earlyYears.mean(), xNew2 = lateYears.mean(),
    ciType = 'lrt')
result3 = climextremes.calc_riskRatio_gev(
    returnValue = 3,
    y1 = numpy.array(y1), y2 = numpy.array(y2),
    x1 = earlyYears, x2 = lateYears,
    locationFun1 = 1, locationFun2 = 1,
    xNew1 = earlyYears.mean(), xNew2 = lateYears.mean(),
    ciType = numpy.array(('lrt', 'boot_perc')))
result4 = climextremes.calc_riskRatio_gev(
    returnValue = 3,
    y1 = numpy.array(y1), y2 = numpy.array(y2),
    x1 = earlyYears, x2 = lateYears,
    locationFun1 = 1, locationFun2 = 1,
    xNew1 = earlyYears.mean(), xNew2 = lateYears.mean(),
    ciType = numpy.array(('lrt', 'boot_perc')), bootControl = {'n':250,'seed':2})
result5 = climextremes.calc_riskRatio_gev(
    returnValue = 3,
    y1 = numpy.array(y1), y2 = numpy.array(y2),
    x1 = earlyYears, x2 = lateYears,
    locationFun1 = 1, locationFun2 = 1,
    xNew1 = earlyYears.mean(), xNew2 = lateYears.mean(),
    ciType = 'lrt', optimArgs = {'method': "BFGS"})

Fort = climextremes.Fort
threshold = 0.395
FortExc = Fort[Fort.Prec > threshold]

earlyPeriod = numpy.array((1900, 1930))
earlyYears = numpy.array(range(earlyPeriod[0], earlyPeriod[1]))
latePeriod = numpy.array((1970, 2000))
lateYears = numpy.array(range(latePeriod[0], latePeriod[1]))

y1 = FortExc.Prec[FortExc.year < earlyPeriod[1]]
y2 = FortExc.Prec[FortExc.year >= latePeriod[0]]
block1 = FortExc.year[FortExc.year < earlyPeriod[1]]
block2 = FortExc.year[FortExc.year >= latePeriod[0]]

# contrast late period with early period, assuming a nonstationary fit
# within each time period and finding RR based on midpoint of each period
result = climextremes.calc_riskRatio_pot(
    returnValue = 3,
    y1 = numpy.array(y1), y2 = numpy.array(y2),
    x1 = earlyYears, x2 = lateYears,
    threshold1 = threshold, threshold2 = threshold,
    locationFun1 = 1, locationFun2 = 1,
    xNew1 = earlyYears.mean(), xNew2 = lateYears.mean(),
    blockIndex1 = numpy.array(block1), blockIndex2 = numpy.array(block2),
    firstBlock1 = earlyYears[0], firstBlock2 = lateYears[0],
    ciType = 'lrt')
