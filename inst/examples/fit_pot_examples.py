Fort = climextremes.Fort

firstYr = min(Fort.year)
yrs = numpy.array(range(int(firstYr), int(max(Fort.year)+1)))
nYrs = len(yrs)
yrsToPred = numpy.array([min(Fort.year), max(Fort.year)])

threshold = 0.395

FortExc = Fort[Fort.Prec > threshold]

# stationary fit
result = climextremes.fit_pot(numpy.array(FortExc.Prec), nBlocks = nYrs, threshold = threshold, firstBlock = firstYr, blockIndex = numpy.array(FortExc.year), getParams = True, returnPeriod = 20, returnValue = 3.5, bootSE = True)
result['returnValue']
result['se_returnValue']       # return value standard error (asymptotic)
result['se_returnValue_boot']  # return value standard error (bootstrapping)
result['logReturnProb']        # log of probability of exceeding 'returnValue'
result['mle']                  # MLE array 
result['mle_names']            # names for MLE array 
result['mle'][2]               # MLE for shape parameter

result['numBootFailures']      # number of bootstrap datasets for which the model could not be fit; if this is non-negligible relative to the number of bootstrap samples (default of 250), interpret the bootstrap results with caution

# nonstationary fit with location linear in year and two return values requested
result = climextremes.fit_pot(numpy.array(FortExc.Prec), x = yrs, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = 1, getParams = True, returnPeriod = numpy.array([20, 30]), returnValue = 3.5, xNew = yrsToPred, bootSE = False)
result['returnValue']
result['se_returnValue']

# fit with location a function of two covariates
# here I'll use year and a random vector just to illustrate syntax
# make 'x' be a 2-column numpy array, each column a covariate
# 'xNew' also needs to have 2 columns, each row is a different set of covariate values
tmp = numpy.random.rand(nYrs)
covByBlock = numpy.c_[yrs, numpy.random.rand(nYrs)]
result = climextremes.fit_pot(numpy.array(FortExc.Prec), x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort.year), 0], [max(Fort.year), 0]]), bootSE = False) 

# with declustering (using max of exceedances on contiguous days)
result = climextremes.fit_pot(numpy.array(FortExc.Prec), x = yrs, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), index = numpy.array(FortExc.obs), locationFun = 1, declustering = "noruns", getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = yrsToPred, bootSE = False) 
result['returnValue']
result['se_returnValue']    

# with declustering (consider sequential blocks of 5 days and only use the max of any exceedances within a block)
result = climextremes.fit_pot(numpy.array(FortExc.Prec), x = yrs, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = numpy.array(FortExc.year), index = numpy.array(FortExc.obs), locationFun = 1, declustering = 5, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = yrsToPred, bootSE = False) 
result['returnValue']
result['se_returnValue']    

# with replicates; for illustration here, I'll just duplicate the Fort data
result = climextremes.fit_pot(numpy.append(numpy.array(FortExc.Prec), numpy.array(FortExc.Prec)), x = yrs, firstBlock = firstYr, nBlocks = nYrs, nReplicates = 2, threshold = threshold, blockIndex = numpy.append(FortExc.year, FortExc.year), locationFun = 1, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = yrsToPred, bootSE = False) 
result['returnValue']
result['se_returnValue']    

# analysis of seasonal total precipitation
tmp = Fort[numpy.logical_and(Fort['month'] < 9, Fort['month'] > 5)]
FortSummerTotal = tmp.groupby('year').sum()[['Prec']]
FortSummerTotal.reset_index(inplace=True)
threshold = numpy.percentile(FortSummerTotal.Prec, 80)
FortSummerTotalExc = FortSummerTotal[FortSummerTotal.Prec > threshold]

result = climextremes.fit_pot(numpy.array(FortSummerTotalExc.Prec), x = yrs, firstBlock = firstYr, nBlocks = nYrs, blockIndex = numpy.array(FortSummerTotalExc.year), locationFun = 1, threshold = threshold, getParams = True, returnPeriod = 20, returnValue = 10, xNew = yrsToPred, bootSE = False)
result['returnValue']
result['se_returnValue']    

# modifying control arguments and seeing more information on the optimization 
result = climextremes.fit_pot(numpy.array(FortSummerTotalExc.Prec), x = yrs, firstBlock = firstYr, nBlocks = nYrs, blockIndex = numpy.array(FortSummerTotalExc.year), locationFun = 1, threshold = threshold, getParams = True, returnPeriod = 20, returnValue = 10, xNew = yrsToPred, bootSE = True, bootControl = {'n':150, 'seed':3}, getFit = True)
result['info']   # information on the optimization
result['info']['counts']  # number of evaluations in the optimization
result['info']['counts_names'] # names to interpret 'counts'
result['numBootFailures']      # number of bootstrap datasets for which the model could not be fit; if this is non-negligible relative to the number of bootstrap samples (default of 250), interpret the bootstrap results with caution

# result['fit']  # voluminous output from the R function that does the fitting
