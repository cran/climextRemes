
# stationary fit
result = climextRemes.fit_pot(fc.FortExc_Prec, firstBlock = fc.firstYr, nBlocks = fc.nYrs, threshold = fc.threshold, blockIndex = fc.FortExc_year, getParams = True, returnPeriod = 20, returnValue = 3.5, bootSE = False)
result.returnValue
result.se_returnValue
result.mle

# fit with location linear in year
result = climextRemes.fit_pot(fc.FortExc_Prec, x = numpy.array(fc.yrs), firstBlock = fc.firstYr, nBlocks = fc.nYrs, threshold = fc.threshold, blockIndex = fc.FortExc_year, locationFun = 1, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(fc.Fort_year), max(fc.Fort_year)]), bootSE = False)
result.returnValue
result.se_returnValue

# with declustering (using max of exceedances on contiguous days)
result = climextRemes.fit_pot(fc.FortExc_Prec, x = numpy.array(fc.yrs), firstBlock = fc.firstYr, nBlocks = fc.nYrs, threshold = fc.threshold, blockIndex = fc.FortExc_year, index = fc.FortExc_index, locationFun = 1, declustering = "noruns", getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(fc.Fort_year), max(fc.Fort_year)]), bootSE = False) 
result.returnValue
result.se_returnValue

# with replicates; for illustration here, I'll just duplicate the Fort data
result = climextRemes.fit_pot(numpy.append(fc.FortExc_Prec, fc.FortExc_Prec), x = numpy.array(fc.yrs), firstBlock = fc.firstYr, nBlocks = fc.nYrs, nReplicates = 2, threshold = fc.threshold, blockIndex = numpy.append(fc.FortExc_year, fc.FortExc_year), locationFun = 1, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(fc.Fort_year), max(fc.Fort_year)]), bootSE = False) 
result.returnValue
result.se_returnValue

# analysis of seasonal total precipitation
result = climextRemes.fit_pot(fc.FortSummerSumExc_Prec, x = numpy.array(fc.yrs), firstBlock = fc.firstYr, nBlocks = fc.nYrs, blockIndex = fc.FortSummerSumExc_year, locationFun = 1, threshold = fc.FortSummerThreshold, getParams = True, returnPeriod = 20, returnValue = 10, xNew = numpy.array([min(fc.Fort_year), max(fc.Fort_year)]), bootSE = False)
result.returnValue
result.se_returnValue
