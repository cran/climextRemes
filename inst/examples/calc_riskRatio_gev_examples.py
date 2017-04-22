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
    lrtCI = True)

result['logRiskRatio']
result['se_logRiskRatio']
result['riskRatio']
result['ci_riskRatio']
result['ci_riskRatio_lrt']

