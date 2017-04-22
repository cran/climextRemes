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
    lrtCI = True)

result['logRiskRatio']
result['se_logRiskRatio']
result['riskRatio']
result['ci_riskRatio']
result['ci_riskRatio_lrt']

