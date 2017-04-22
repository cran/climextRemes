result = climextremes.calc_riskRatio_binom(numpy.array((40, 8)), 
                                           numpy.array((400, 400)),
                                  lrtCI = True)
result['logRiskRatio']
result['se_logRiskRatio']
result['riskRatio']
result['ci_riskRatio']
result['ci_riskRatio_lrt']

# LRT method can estimate lower confidence interval endpoint,
# even if estimated risk ratio is infinity
result = climextremes.calc_riskRatio_binom(numpy.array((4, 0)), 
                                           numpy.array((100, 100)),
                                           lrtCI = True)
result['logRiskRatio']
result['se_logRiskRatio']
result['riskRatio']
result['ci_riskRatio']
result['ci_riskRatio_lrt']


