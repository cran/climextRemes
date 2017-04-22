Fort = climextremes.Fort

FortMax = Fort.groupby('year').max()[['Prec']]
FortMax.reset_index(inplace=True)

# stationary fit
result = climextremes.fit_gev(numpy.array(FortMax.Prec), returnPeriod = 20, returnValue = 3.5, getParams = True, bootSE = True)
result['returnValue']
result['se_returnValue']       # return value standard error (asymptotic)
result['se_returnValue_boot']  # return value standard error (bootstrapping)
result['logReturnProb']        # log of probability of exceeding 'returnValue'
result['mle']                  # MLE array 
result['mle_names']            # names for MLE array 
result['mle'][2]               # MLE for shape parameter

# modifying the bootstrapping specifications
result = climextremes.fit_gev(numpy.array(FortMax.Prec), returnPeriod = 20, returnValue = 3.5, getParams = True, bootSE = True, bootControl = {'n': 100, 'seed': 3})
result['se_returnValue_boot']  # return value standard error (bootstrapping)

yrsToPred = numpy.array([min(Fort.year), max(Fort.year)])

# nonstationary fit with location linear in year and two return values requested
result_ns = climextremes.fit_gev(numpy.array(FortMax.Prec), numpy.array(FortMax.year), locationFun = 1, returnPeriod = numpy.array([20, 30]), returnValue = 3.5, xNew = yrsToPred, getParams = True, bootSE = False)
result_ns['returnValue']
result_ns['se_returnValue']


