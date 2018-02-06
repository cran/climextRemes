result = climextremes.calc_riskRatio_binom(numpy.array((40, 8)), 
                                           numpy.array((400, 400)),
                                           ciType = numpy.array(('lrt', 'boot_stud', 'koopman')))
result['logRiskRatio']
result['se_logRiskRatio']
result['riskRatio']
result['ci_riskRatio_lrt']
result['ci_riskRatio_koopman']
result['ci_riskRatio_boot_stud']

# Koopman and LRT method can estimate lower confidence interval endpoint,
# even if estimated risk ratio is infinity
result = climextremes.calc_riskRatio_binom(numpy.array((4, 0)), 
                                           numpy.array((100, 100)),
                                           ciType = numpy.array(('lrt', 'boot_stud', 'koopman')))
result['logRiskRatio']
result['se_logRiskRatio']
result['riskRatio']
result['ci_riskRatio_lrt']
result['ci_riskRatio_koopman']
result['ci_riskRatio_boot_stud']

