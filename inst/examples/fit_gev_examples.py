fc = climextRemes.FortCollinsPrecipData()

# Hari, perhaps convert the result of the following into a component of fc? Not sure without more digging how to best do this calculation within Python given that the fc object is some sort of R object and not a pandas DataFrame.

# this computes max value within each year
FortMax = robjects.r('''
   aggregate(Prec ~ year, data = Fort, max)
''')

# stationary fit
result = climextRemes.fit_gev(FortMax[1], returnPeriod = 20, returnValue = 3.5, getParams = True, bootSE = False)
result.returnValue
result.se_returnValue
result.mle  # Hari, per email, I'd like this to return the array of numbers and result.mle_names the names ('location','scale','shape', etc.)

# nonstationary fit with location linear in year
result_ns = climextRemes.fit_gev(FortMax[1], FortMax[0], locationFun = 1, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(fc.Fort_year), max(fc.Fort_year)]), getParams = True, bootSE = False)
result_ns.returnValue
result_ns.se_returnValue
