#' Compute risk ratio and uncertainty based on binomial models for counts of events relative to possible number of events
#'
#' Compute risk ratio and uncertainty by fitting binomial models to counts of events relative to possible number of events. The risk ratio is the ratio of the probability of an event under the model fit to the first dataset to the probability under the model fit to the second dataset. Default standard errors are based on the usual MLE asymptotics using a delta-method-based approximation, but standard errors based on the nonparametric bootstrap and on a likelihood ratio procedure can also be computed.
#'
#' @param y vector of two values, the number of events in the two scenarios
#' @param n vector of two values, the number of samples (possible occurrences of events) in the two scenarios 
#' @param ciLevel statistical confidence level for confidence intervals; in repeated experimentation, this proportion of confidence intervals should contain the true risk ratio. Note that if only one endpoint of the resulting interval is used, for example the lower bound, then the effective confidence level increases by half of one minus \code{ciLevel}. For example, a two-sided 0.90 confidence interval corresponds to a one-sided 0.95 confidence interval.
#' @param bootSE logical indicating whether to use the bootstrap to estimate standard errors.
#' @param bootControl a list of control parameters for the bootstrapping. See \code{Details}.
#' @param lrtCI logical indicating whether to calculate a likelihood ratio-based confidence interval
#' @param lrtControl list containing a single component, \code{bounds}, which sets the range inside which the algorithm searches for the endpoints of the likelihood ratio-based confidence interval. This avoids numerical issues with endpoints converging to zero and infinity. If an endpoint is not found within the interval, it is set to \code{NA}.
#' @author Christopher J. Paciorek
#' @export
#' @details
#' See \code{\link{fit_pot}} for information on the \code{bootControl} argument. 
#' @references
#' Paciorek et al. methods paper being finalized.
#' @examples
#' calc_riskRatio_binom(c(4,0), rep(100, 2), bootSE = FALSE, lrtCI = TRUE)
#' calc_riskRatio_binom(c(40, 8), rep(400, 2), bootSE = TRUE, lrtCI = TRUE)
calc_riskRatio_binom <- function(y, n, ciLevel = 0.90,
                                 bootSE = FALSE, bootControl = list(seed = 0, n = 500),
                                 lrtCI = FALSE, lrtControl = list(bounds = c(.01, 100))) {
    results <- list()
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE)

    if(length(n) == 1) n <- rep(n, 2)
    results$logRiskRatio <- log(y[1]) - log(n[1]) - log(y[2]) + log(n[2])
    p <- y / n
    results$se_logRiskRatio <- sqrt( (1-p[1]) / y[1] + (1-p[2]) / y[2] )

    results$riskRatio <- exp(results$logRiskRatio)
    results$ci_riskRatio <- exp(results$logRiskRatio + c(-1, 1)*z_alpha*results$se_logRiskRatio)
    
    if(bootSE) {
        if(any(p == 0 || p == 1)) {
            results$se_logRiskRatio_boot <- NA
        } else {
            if(length(bootControl$seed) == 1){
                set.seed(bootControl$seed)
            } else{
                .Random.seed <- bootControl$seed
            }
            b_y1 <- stats::rbinom(n = bootControl$n, size = n[1], prob = p[1])
            b_y2 <- stats::rbinom(n = bootControl$n, size = n[2], prob = p[2])
            results$se_logRiskRatio_boot <- stats::sd(log(b_y1) - log(n[1]) - log(b_y2) + log(n[2]))
            results$ci_riskRatio_boot <- exp(results$logRiskRatio + c(-1, 1)*z_alpha*results$se_logRiskRatio_boot)
        }
    }
    if(lrtCI) 
        results$ci_riskRatio_lrt <- calc_riskRatio_lrt_binom(y, n, ciLevel = ciLevel, bounds = lrtControl$bounds)
    return(results)
} # end calc_riskRatio_binom()
    

calc_riskRatio_lrt_binom <- function(y, n, ciLevel = 0.90, bounds) {
    logLik <- function(y, n, p = y / n)
        stats::dbinom(y[1], n[1], p[1], log = TRUE) + stats::dbinom(y[2], n[2], p[2], log = TRUE)

    restricted_mle <- function(pHat, rr) {
    # from Farrington and Manning, Stats in Med, 9:1447 (1990)
    # easily derivable from likelihood with simple plug-in constraint
        a <- 2
        b <- -(rr * (1 + pHat[2]) + 1 + pHat[1])
        c <- rr*(pHat[1] + pHat[2])
        pAtilde <- (-b - sqrt(b*b - 4*a*c)) / (2*a)
        return(pAtilde)
    }
    
    objfun <- function(logrr, y, n, logLikHat, cutoff, lowerEnd = TRUE) {
        rr <- exp(logrr)
        pAtilde <- restricted_mle(y/n, rr)
        logLikConstr <- logLik(y, n, c(pAtilde, pAtilde/rr))
        if(pAtilde > rr || pAtilde > 1) {
            # warning('calc_riskRatio_lrt_binom: infeasible value for probability for one of the groups -- optimization may not converge.')
            cond <- FALSE
        } else cond <- logLikHat - logLikConstr < cutoff
        if(lowerEnd) toggle <- 1 else toggle <- -1
        
        rrHat <- (y[1]/n[1]) / (y[2]/n[2])
        return(rr*cond + toggle*1000000*(!cond) +
           # ensure that objective when constraint not satisfied gives unimodality
           (-toggle*rr)*(!cond)*(rr < rrHat) + toggle*rr*(!cond)*(rr > rrHat) )
    }

    cutoff <- 0.5*stats::qchisq(ciLevel, 1)
    logLikHat <- logLik(y, n)
    pHat <- y/n
    rrHat <- pHat[1] / pHat[2]
    # optimize on log scale for hopefully better numerical performance
    if(pHat[1] == 0 && pHat[2] == 0) return(c(NA, NA))
    if(pHat[1] == 1 && pHat[2] == 1) return(c(NA, NA))
    if(pHat[1] == 0) {
        lower <- 0
    } else {
        intvl <- log(c(bounds[1], min(rrHat, bounds[2])))
        lower <- stats::optimize(objfun, interval = intvl, y = y, n = n, logLikHat = logLikHat, cutoff = cutoff)$minimum
        if(all.equal(lower, intvl[1]) == TRUE || all.equal(lower, intvl[2]) == TRUE)
            warning("optimization likely failed to converge: lower endpoint equal to one endpoint of search interval.")
    }
    if(pHat[2] == 0) {
        upper <- Inf
    } else {
        intvl <- log(c(max(rrHat, bounds[1]), bounds[2]))
        upper <- stats::optimize(objfun, interval = intvl, y = y, n = n, logLikHat = logLikHat, cutoff = cutoff, lowerEnd = FALSE, maximum = TRUE)$maximum
        if(all.equal(upper, intvl[1]) == TRUE || all.equal(upper, intvl[2]) == TRUE)
            warning("optimization likely failed to converge: upper endpoint equal to one endpoint of search interval.")
    }
    return(exp(c(lower, upper)))
} # end calc_riskRatio_lrt_binom()


#' Compute risk ratio and uncertainty based on peaks-over-threshold models fit to exceedances over a threshold
#'
#' Compute risk ratio and uncertainty by fitting peaks-over-threshold model, designed specifically for climate data, to exceedance-only data, using the point process approach. The risk ratio is the ratio of the probability of exceedance of a pre-specified value under the model fit to the first dataset to the probability under the model fit to the second dataset. Default standard errors are based on the usual MLE asymptotics using a delta-method-based approximation, but standard errors based on the nonparametric bootstrap and on a likelihood ratio procedure can also be computed.
#'
#' @param returnValue numeric value giving the value for which the risk ratio should be calculated, where the resulting period will be the average number of blocks until the value is exceeded and the probability the probability of exceeding the value in any single block.
#' @param y1 a numeric vector of exceedance values for the first dataset (values of the outcome variable above the threshold).
#' @param y2 a numeric vector of exceedance values for the second dataset (values of the outcome variable above the threshold).
#' @param x1 a data frame, or object that can be converted to a data frame with columns corresponding to covariate/predictor/feature variables and each row containing the values of the variable for a block (e.g., often a year with climate data) for the first dataset. The number of rows must equal the number of blocks.
#' @param x2 analogous to \code{x1} but for the second dataset
#' @param threshold1 a single numeric value for constant threshold or a numeric vector with length equal to the number of blocks, indicating the threshold for each block for the first dataset.
#' @param threshold2 analogous to \code{threshold1} but for the second dataset
#' @param locationFun1 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the location parameter using columns from \code{x1} for the first dataset. \code{x1} must be supplied if this is anything other than NULL or ~1.
#' @param locationFun2 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the location parameter using columns from \code{x2} for the second dataset. \code{x2} must be supplied if this is anything other than NULL or ~1.
#' @param scaleFun1 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the log of the scale parameter using columns from \code{x1} for the first dataset.  \code{x1} must be supplied if this is anything other than NULL or ~1.
#' @param scaleFun2 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the log of the scale parameter using columns from \code{x2} for the second dataset.  \code{x2} must be supplied if this is anything other than NULL or ~1.
#' @param shapeFun1 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the shape parameter using columns from \code{x1} for the first dataset.  \code{x1} must be supplied if this is anything other than NULL or ~1.
#' @param shapeFun2 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the shape parameter using columns from \code{x2} for the first dataset.  \code{x2} must be supplied if this is anything other than NULL or ~1.
#' @param nBlocks1 number of blocks (e.g., a block will often be a year with climate data) in first dataset; note this value determines the interpretation of return values/periods/probabilities; see \code{returnPeriod} and \code{returnValue}.
#' @param nBlocks2 number of blocks (e.g., a block will often be a year with climate data) in second dataset; note this value determines the interpretation of return values/periods/probabilities; see \code{returnPeriod} and \code{returnValue}.
#' @param blockIndex1 numeric vector providing the index of the block corresponding to each element of \code{y1}. Used only when \code{x1} is provided to match exceedances to the covariate/predictor/feature value for the exceedance or when using bootstrapping with the resampling based on blocks based on the \code{by} element of \code{bootControl}. If \code{firstBlock1} is not equal to one, then \code{blockIndex1} need not have one as its smallest possible value.
#' @param blockIndex2 numeric vector providing the index of the block corresponding to each element of \code{y2}. Analogous to \code{blockIndex1}.
#' @param firstBlock1 single numeric value indicating the numeric value of the first possible block of \code{blockIndex1}. For example the values in \code{blockIndex1} might indicate the year of each exceedance with the first year of data being 1969, in which case \code{firstBlock1} would be 1969. Note that the first block may not have any exceedances so it may not be represented in \code{blockIndex1}. Used only to adjust \code{blockIndex1} so that the block indices start at one and therefore correspond to the rows of \code{x1}.
#' @param firstBlock2 single numeric value indicating the numeric value of the first possible block of \code{blockIndex2}. Analogous to \code{firstBlock1}.
#' @param index1 (optional) numeric vector providing the integer-valued index (e.g., julian day for daily climate data) corresponding to each element of \code{y1}. For example if there are 10 original observations and the third, fourth, and seventh values are exceedances, then \code{index1} would be the vector 3,4,7. Used only when \code{declustering} is provided to determine which exceedances occur sequentially or within a contiguous set of values of a given length. The actual values are arbitrary; only the lags between the values are used.
#' @param index2 (optional) numeric vector providing the integer-valued index (e.g., julian day for daily climate data) corresponding to each element of \code{y2}. Analogous to \code{index1}.
#' @param nReplicates1 numeric value indicating the number of replicates for the first dataset.
#' @param nReplicates2 numeric value indicating the number of replicates for the second dataset.
#' @param replicateIndex1 numeric vector providing the index of the replicate corresponding to each element of \code{y1}. Used for three purposes: (1) when using bootstrapping with the resampling based on replicates based on the \code{by} element of \code{bootControl}, (2) to avoid treating values in different replicates as potentially being sequential or within a short interval when removing values based on \code{declustering}, and (3) to match outcomes to \code{weights} or \code{proportionMissing} when either vary by replicate.  
#' @param replicateIndex2 numeric vector providing the index of the replicate corresponding to each element of \code{y2}. Analogous to \code{replicateIndex1}.
#' @param weights1 a vector or matrix providing the weights by block for the first dataset. When there is only one replicate or the weights do not vary by replicate, a vector of length equal to the number of blocks. When weights vary by replicate, a matrix with rows corresponding to blocks and columns to replicates. Likelihood contribution of each block is multiplied by the corresponding weight. 
#' @param weights2 a vector or matrix providing the weights by block for the second dataset. Analogous to \code{weights1}.
#' @param proportionMissing1 a numeric value, vector or matrix indicating the proportion of missing values in the original first dataset before exceedances were selected. When the proportion missing is the same for all blocks and replicates, a single value. When there is only one replicate or the weights do not vary by replicate, a vector of length equal to the number of blocks. When weights vary by replicate, a matrix with rows corresponding to blocks and columns to replicates.
#' @param proportionMissing2 a numeric value, vector or matrix indicating the proportion of missing values in the original second dataset before exceedances were selected. Analogous to \code{proportionMissing1}.
#' @param xNew1 object of the same form as \code{x1}, providing covariate/predictor/feature values for which log risk ratios are desired.
#' @param xNew2 object of the same form as \code{x2}, providing covariate/predictor/feature values for which log risk ratios are desired. Must provide the same number of covariate sets as \code{xNew1} as the risk ratio is based on contrasting return probabilities under \code{xNew1} and \code{xNew2}.
#' @param declustering one of \code{NULL}, \code{"noruns"}, or a number. If 'noruns' is specified, only the maximum (or minimum if upperTail = FALSE) value within a set of exceedances corresponding to successive indices is included. If a number, this should indicate the size of the interval (which will be used with the \code{index} argument) within which to allow only the largest (or smallest if upperTail = FALSE) value.
#' @param upperTail logical indicating whether one is working with exceedances over a high threshold (TRUE) or exceedances under a low threshold (FALSE); in the latter case, the function works with the negative of the values and the threshold, changing the sign of the resulting location parameters.
#' @param scaling1 positive-valued scalar used to scale the data values of the first dataset for more robust optimization performance. When multiplied by the values, it should produce values with magnitude around 1.
#' @param scaling2 positive-valued scalar used to scale the data values of the second dataset for more robust optimization performance. When multiplied by the values, it should produce values with magnitude around 1.
#' @param ciLevel statistical confidence level for confidence intervals; in repeated experimentation, this proportion of confidence intervals should contain the true risk ratio. Note that if only one endpoint of the resulting interval is used, for example the lower bound, then the effective confidence level increases by half of one minus \code{ciLevel}. For example, a two-sided 0.90 confidence interval corresponds to a one-sided 0.95 confidence interval.
#' @param bootSE logical indicating whether to use the bootstrap to estimate standard errors.
#' @param bootControl a list of control parameters for the bootstrapping. See \code{Details}.
#' @param lrtCI logical indicating whether to calculate a likelihood ratio-based confidence interval
#' @param lrtControl list containing a single component, \code{bounds}, which sets the range inside which the algorithm searches for the endpoints of the likelihood ratio-based confidence interval. This avoids numerical issues with endpoints converging to zero and infinity. If an endpoint is not found within the interval, it is set to \code{NA}.
#' @param optimArgs a list with named components matching exactly any arguments that the user wishes to pass to \code{optim}. See \code{help(optim)} for details. Of particular note, \code{'method'} can be used to choose the optimization method used for maximizing the log-likelihood to fit the model and \code{'control=list(maxit=VALUE)'} for a user-chosen VALUE can be used to increase the number of iterations if the optimization is converging slowly.
#' @author Christopher J. Paciorek
#' @export
#' @details
#' See \code{\link{fit_pot}} for more details on fitting the peaks-over-threshold model for each dataset, including details on blocking and replication. Also see \code{\link{fit_pot}} for information on the \code{bootControl} argument. 
#' @references
#' Jeon S., C.J. Paciorek, and M.F. Wehner. 2016. Quantile-based bias correction and uncertainty quantification of extreme event attribution statements. Weather and Climate Extremes. In press. arXiv preprint: http://arxiv.org/abs/1602.04139.
#' 
#' @examples
#' # need examples
calc_riskRatio_pot <- function(returnValue, y1, y2, x1 = NULL, x2 = x1,
                                  threshold1, threshold2 = threshold1, locationFun1 = NULL,
                                  locationFun2 = locationFun1, scaleFun1 = NULL, scaleFun2 = scaleFun1,
                                  shapeFun1 = NULL, shapeFun2 = shapeFun1, nBlocks1 = nrow(x1),
                                  nBlocks2 = nrow(x2), blockIndex1 = NULL, blockIndex2 = NULL, firstBlock1 = 1,
                                  firstBlock2 = 1, index1 = NULL, index2 = NULL, nReplicates1 = 1,
                                  nReplicates2 = 1, replicateIndex1 = NULL, replicateIndex2 = NULL,
                                  weights1 = NULL, weights2 = NULL, proportionMissing1 = NULL,
                                  proportionMissing2 = NULL, xNew1 = NULL, xNew2 = NULL, declustering = NULL,
                                  upperTail = TRUE, scaling1 = 1, scaling2 = 1, ciLevel = 0.90, bootSE = FALSE,
                                  bootControl = list(seed = 0, n = 250, by = "block"), 
                                  lrtCI = FALSE, lrtControl = list(bounds = c(.01, 100)),
                                  optimArgs = list(method = 'Nelder-Mead')) {

    if(missing(returnValue))
        stop("calc_riskRatio_pot: 'returnValue' must be provided.")
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE) 
    results <- list()

    fit1 <- fit_pot(y1, x = x1, threshold = threshold1, locationFun = locationFun1,
                        scaleFun = scaleFun1, shapeFun = shapeFun1, nBlocks = nBlocks1,
                        blockIndex = blockIndex1, firstBlock = firstBlock1, index = index1,
                        nReplicates = nReplicates1, replicateIndex = replicateIndex1,
                        weights = weights1, proportionMissing = proportionMissing1, returnValue = returnValue,
                        xNew = xNew1, declustering = declustering, upperTail = upperTail,
                        scaling = scaling1, bootSE = bootSE, bootControl = bootControl,
                        optimArgs = optimArgs, getFit = TRUE, .getInputs = TRUE) 
    fit2 <- fit_pot(y2, x = x2, threshold = threshold2, locationFun = locationFun2,
                        scaleFun = scaleFun2, shapeFun = shapeFun2, nBlocks = nBlocks2,
                        blockIndex = blockIndex2, firstBlock = firstBlock2, index = index2,
                        nReplicates = nReplicates2, replicateIndex = replicateIndex2,
                        weights = weights2, proportionMissing = proportionMissing2, returnValue = returnValue,
                        xNew = xNew2, declustering = declustering, upperTail = upperTail,
                        scaling = scaling2, bootSE = bootSE, bootControl = bootControl,
                        optimArgs = optimArgs, getFit = TRUE, .getInputs = TRUE)
    if(fit1$info$failure || fit2$info$failure) {
        warning("calc_riskRatio_pot: fitting failed for one of two datasets.")
        results$logRiskRatio <- results$se_logRiskRatio <- results$ci_riskRatio <- NA
        if(bootSE)
            results$se_logRiskRatio_boot <- results$ci_riskRatio_boot <- NA
        if(lrtCI)
            results$ci_riskRatio_lrt <- NA
    } else {
        if(length(fit1$logReturnProb) != length(fit2$logReturnProb))
            stop("calc_riskRatio_pot: number of return probabilities calculated for each model fit must be the same; for nonstationary models this is determined by the number of covariate set inputs (provided in 'x' or 'xNew').")
        results$logRiskRatio <- fit1$logReturnProb - fit2$logReturnProb
                                        # delta method
        results$se_logRiskRatio <- sqrt(fit1$se_logReturnProb^2 + fit2$se_logReturnProb^2)
        results$riskRatio <- exp(results$logRiskRatio)
        if(length(results$logRiskRatio) > 1) {
            results$ci_riskRatio <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio,
                                              results$logRiskRatio + z_alpha*results$se_logRiskRatio))
        } else results$ci_riskRatio <- exp(results$logRiskRatio + c(-1,1)*z_alpha*results$se_logRiskRatio)
        if(bootSE) { 
            results$se_logRiskRatio_boot <- sqrt(fit1$se_logReturnProb_boot^2 + fit2$se_logReturnProb_boot^2)
            if(length(results$logRiskRatio) > 1) {
                results$ci_riskRatio_boot <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio_boot,
                                                       results$logRiskRatio + z_alpha*results$se_logRiskRatio_boot))
            } else results$ci_riskRatio_boot <- exp(results$logRiskRatio + c(-1,1)*z_alpha*results$se_logRiskRatio_boot)
        }
        if(lrtCI) {
            lControl <- list(bounds = c(.01, 100))
            lControl[names(lrtControl)] <- lrtControl
            lrtControl <- lControl
            oArgs <- list(method = "Nelder-Mead", lower = -Inf, upper = Inf, control = list())
            oArgs[names(optimArgs)] <- optimArgs

            
            results$ci_riskRatio_lrt <- calc_riskRatio_lrt(fit1, fit2, returnValue, ciLevel = ciLevel, bounds = lrtControl$bounds, type = "PP", optimArgs = oArgs)
        }
    }
    return(results)
}

#' Compute risk ratio and uncertainty based on generalized extreme value model fit to block maxima or minima
#'
#' Compute risk ratio and uncertainty by fitting generalized extreme value model, designed specifically for climate data, to exceedance-only data, using the point process approach. The risk ratio is the ratio of the probability of exceedance of a pre-specified value under the model fit to the first dataset to the probability under the model fit to the second dataset. Default standard errors are based on the usual MLE asymptotics using a delta-method-based approximation, but standard errors based on the nonparametric bootstrap and on a likelihood ratio procedure can also be computed.
#'
#' @param returnValue numeric value giving the value for which the risk ratio should be calculated, where the resulting period will be the average number of blocks until the value is exceeded and the probability the probability of exceeding the value in any single block.
#' @param y1 a numeric vector of observed maxima or minima values for the first dataset. See \code{Details} for how the values of \code{y1} should be ordered if there are multiple replicates and the values of \code{x1} are identical for all replicates.
#' @param y2 a numeric vector of observed maxima or minima values for the second dataset. Analogous to \code{y1}.
#' @param x1 a data frame, or object that can be converted to a data frame with columns corresponding to covariate/predictor/feature variables and each row containing the values of the variable for the corresponding observed maximum/minimum. The number of rows should either equal the length of \code{y1} or (if there is more than one replicate) it can optionally equal the number of observations in a single replicate, in which case the values will be assumed to be the same for all replicates. 
#' @param x2 analogous to \code{x1} but for the second dataset
#' @param locationFun1 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the location parameter using columns from \code{x1} for the first dataset. \code{x1} must be supplied if this is anything other than NULL or ~1.
#' @param locationFun2 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the location parameter using columns from \code{x2} for the second dataset. \code{x2} must be supplied if this is anything other than NULL or ~1.
#' @param scaleFun1 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the log of the scale parameter using columns from \code{x1} for the first dataset.  \code{x1} must be supplied if this is anything other than NULL or ~1.
#' @param scaleFun2 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the log of the scale parameter using columns from \code{x2} for the second dataset.  \code{x2} must be supplied if this is anything other than NULL or ~1.
#' @param shapeFun1 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the shape parameter using columns from \code{x1} for the first dataset.  \code{x1} must be supplied if this is anything other than NULL or ~1.
#' @param shapeFun2 formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the shape parameter using columns from \code{x2} for the first dataset.  \code{x2} must be supplied if this is anything other than NULL or ~1.
#' @param nReplicates1 numeric value indicating the number of replicates for the first dataset.
#' @param nReplicates2 numeric value indicating the number of replicates for the second dataset.
#' @param replicateIndex1 numeric vector providing the index of the replicate corresponding to each element of \code{y1}. Used (and therefore required) only when using bootstrapping with the resampling by replicates based on the \code{by} element of \code{bootControl}.
#' @param replicateIndex2 numeric vector providing the index of the replicate corresponding to each element of \code{y2}. Analogous to \code{replicateIndex1}.
#' @param weights1 a vector providing the weights for each observation in the first dataset. When there is only one replicate or the weights do not vary by replicate, a vector of length equal to the number of observations. When weights vary by replicate, this should be of equal length to \code{y}. Likelihood contribution of each observation is multiplied by the corresponding weight. 
#' @param weights2 a vector providing the weights for each observation in the second dataset. Analogous to \code{weights1}.
#' @param xNew1 object of the same form as \code{x1}, providing covariate/predictor/feature values for which one desires log risk ratios.
#' @param xNew2 object of the same form as \code{x2}, providing covariate/predictor/feature values for which log risk ratios are desired. Must provide the same number of covariate sets as \code{xNew1} as the risk ratio is based on contrasting return probabilities under \code{xNew1} and \code{xNew2}.
#' @param maxes logical indicating whether analysis is for block maxima (TRUE) or block minima (FALSE); in the latter case, the function works with the negative of the values, changing the sign of the resulting location parameters
#' @param scaling1 positive-valued scalar used to scale the data values of the first dataset for more robust optimization performance. When multiplied by the values, it should produce values with magnitude around 1.
#' @param scaling2 positive-valued scalar used to scale the data values of the second dataset for more robust optimization performance. When multiplied by the values, it should produce values with magnitude around 1.
#' @param ciLevel statistical confidence level for confidence intervals; in repeated experimentation, this proportion of confidence intervals should contain the true risk ratio. Note that if only one endpoint of the resulting interval is used, for example the lower bound, then the effective confidence level increases by half of one minus \code{ciLevel}. For example, a two-sided 0.90 confidence interval corresponds to a one-sided 0.95 confidence interval.
#' @param bootSE logical indicating whether to use the bootstrap to estimate standard errors.
#' @param bootControl a list of control parameters for the bootstrapping. See \code{Details}.
#' @param lrtCI logical indicating whether to calculate a likelihood ratio-based confidence interval
#' @param lrtControl list containing a single component, \code{bounds}, which sets the range inside which the algorithm searches for the endpoints of the likelihood ratio-based confidence interval. This avoids numerical issues with endpoints converging to zero and infinity. If an endpoint is not found within the interval, it is set to \code{NA}.
#' @param optimArgs a list with named components matching exactly any arguments that the user wishes to pass to \code{optim}. See \code{help(optim)} for details. Of particular note, \code{'method'} can be used to choose the optimization method used for maximizing the log-likelihood to fit the model and \code{'control=list(maxit=VALUE)'} for a user-chosen VALUE can be used to increase the number of iterations if the optimization is converging slowly.
#' @author Christopher J. Paciorek
#' @export
#' @details
#' See \code{\link{fit_gev}} for more details on fitting the block maxima model for each dataset, including details on blocking and replication. Also see \code{\link{fit_gev}} for information on the \code{bootControl} argument. 
#' @references
#' Jeon S., C.J. Paciorek, and M.F. Wehner. 2016. Quantile-based bias correction and uncertainty quantification of extreme event attribution statements. Weather and Climate Extremes. In press. arXiv preprint: http://arxiv.org/abs/1602.04139.
#' 
#' @examples
#' # need examples
calc_riskRatio_gev <- function(returnValue, y1, y2, x1 = NULL, x2 = x1,
                               locationFun1 = NULL, locationFun2 = locationFun1,
                               scaleFun1 = NULL, scaleFun2 = scaleFun1,
                               shapeFun1 = NULL, shapeFun2 = shapeFun1,
                               nReplicates1 = 1, nReplicates2 = 1,
                               replicateIndex1 = NULL, replicateIndex2 = NULL,
                               weights1 = NULL, weights2 = NULL, xNew1 = NULL, xNew2 = NULL, 
                               maxes = TRUE, scaling1 = 1, scaling2 = 1,
                               ciLevel = 0.90, bootSE = FALSE,
                               bootControl = list(seed = 0, n = 250, by = "block"),
                               lrtCI = FALSE, lrtControl = list(bounds = c(.01, 100)),
                               optimArgs = list(method = 'Nelder-Mead')) {

    if(missing(returnValue))
        stop("calc_riskRatio_gev: 'returnValue' must be provided.")
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE)
    results <- list()

    fit1 <- fit_gev(y1, x = x1, locationFun = locationFun1,
                        scaleFun = scaleFun1, shapeFun = shapeFun1, 
                        nReplicates = nReplicates1, replicateIndex = replicateIndex1,
                        weights = weights1, returnValue = returnValue,
                        xNew = xNew1, maxes = maxes,
                        scaling = scaling1, bootSE = bootSE, bootControl = bootControl,
                        optimArgs = optimArgs, getFit = TRUE, .getInputs = TRUE) 
    fit2 <- fit_gev(y2, x = x2, locationFun = locationFun2,
                        scaleFun = scaleFun2, shapeFun = shapeFun2, 
                        nReplicates = nReplicates2, replicateIndex = replicateIndex2,
                        weights = weights2, returnValue = returnValue,
                        xNew = xNew2, maxes = maxes,
                        scaling = scaling2, bootSE = bootSE, bootControl = bootControl,
                        optimArgs = optimArgs, getFit = TRUE, .getInputs = TRUE) 

    if(fit1$info$failure || fit2$info$failure) {
        warning("calc_riskRatio_pot: fitting failed for one of two datasets.")
        results$logRiskRatio <- results$se_logRiskRatio <- results$ci_riskRatio <- NA
        if(bootSE)
            results$se_logRiskRatio_boot <- results$ci_riskRatio_boot <- NA
        if(lrtCI)
            results$ci_riskRatio_lrt <- NA
    } else {
        if(length(fit1$logReturnProb) != length(fit2$logReturnProb))
            stop("calc_riskRatio_gev: number of return probabilities calculated for each model fit must be the same; for nonstationary models this is determined by the number of covariate set inputs (provided in 'x' or 'xNew')")
        results$logRiskRatio <- fit1$logReturnProb - fit2$logReturnProb
                                        # delta method
        results$se_logRiskRatio <- sqrt(fit1$se_logReturnProb^2 + fit2$se_logReturnProb^2)
        results$riskRatio <- exp(results$logRiskRatio)
        if(length(results$logRiskRatio) > 1) {
            results$ci_riskRatio <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio,
                                              results$logRiskRatio + z_alpha*results$se_logRiskRatio))
        } else results$ci_riskRatio <- exp(results$logRiskRatio + c(-1,1)*z_alpha*results$se_logRiskRatio)
    
        if(bootSE) { 
            results$se_logRiskRatio_boot <- sqrt(fit1$se_logReturnProb_boot^2 + fit2$se_logReturnProb_boot^2)
            if(length(results$logRiskRatio) > 1) {
                results$ci_riskRatio_boot <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio_boot,
                                                       results$logRiskRatio + z_alpha*results$se_logRiskRatio_boot))
            } else results$ci_riskRatio_boot <- exp(results$logRiskRatio + c(-1,1)*z_alpha*results$se_logRiskRatio_boot)
        }
        if(lrtCI) {
            # using stat initial estimates as w/ nonstationary hard to get constrained optimization to start with legitimate parameter combos 
            ## fit1i <- fit_gev(y1, locationFun = ~1,
            ##                 scaleFun = ~1, shapeFun = ~1, 
            ##                 nReplicates = nReplicates1, replicateIndex = replicateIndex1,
            ##                 weights = weights1, returnValue = returnValue,
            ##                 maxes = maxes,
            ##                 scaling = scaling1, bootSE = FALSE, optimArgs = optimArgs, getParams = TRUE)
            ## fit2i <- fit_gev(y2, locationFun = ~1,
            ##                 scaleFun = ~1, shapeFun = ~1, 
            ##                 nReplicates = nReplicates2, replicateIndex = replicateIndex2,
            ##                 weights = weights2, returnValue = returnValue,
            ##                 maxes = maxes,
            ##                 scaling = scaling2, bootSE = FALSE,
            ##                 optimArgs = optimArgs, getParams = TRUE) 

            lControl <- list(bounds = c(.01, 100))
            lControl[names(lrtControl)] <- lrtControl
            lrtControl <- lControl
            oArgs <- list(method = "Nelder-Mead", lower = -Inf, upper = Inf, control = list())
            oArgs[names(optimArgs)] <- optimArgs
            results$ci_riskRatio_lrt <- calc_riskRatio_lrt(fit1, fit2, returnValue = returnValue, ciLevel = ciLevel, bounds = lrtControl$bounds, type = "GEV", optimArgs = oArgs)
        #    results$ci_riskRatio_lrt <- calc_riskRatio_lrt(fit1, fit2, returnValue, ciLevel = ciLevel, bounds = lrtControl$bounds, type = "GEV", optimArgs = oArgs, init1 = fit1i$mle, init2 = fit2i$mle)
        }
    }
    return(results)
}

calc_riskRatio_lrt <- function(fit1, fit2, returnValue, ciLevel, bounds, type = "PP", optimArgs, init1 = NULL, init2 = NULL) {
    if(!type %in% c('PP', 'GEV'))
        stop("calc_riskRatio_lrt: 'type' must be on of 'PP' or 'GEV'.")

    if(!is.null(init1) && length(init1) != 3) stop("calc_riskRatio_lrt: length of 'init1' should be 3 (should be from a stationary model fit.")
    if(!is.null(init2) && length(init2) != 3) stop("calc_riskRatio_lrt: length of 'init1' should be 3 (should be from a stationary model fit.")
    logLikHat <- -fit1$fit$results$value - fit2$fit$results$value
    inputs1 <- fit1$inputs
    inputs2 <- fit2$inputs

    par1 = fit1$fit$results$par
    par2 = fit2$fit$results$par
    p1 = fit1$fit$results$num.pars
    p2 = fit2$fit$results$num.pars

    if(is.null(inputs1$xNew)) inputs1$xUse <- inputs1$x else inputs1$xUse <- inputs1$xNew
    if(is.null(inputs2$xNew)) inputs2$xUse <- inputs2$x else inputs2$xUse <- inputs2$xNew

    inputs1$mUse <- nrow(inputs1$xUse)
    inputs2$mUse <- nrow(inputs2$xUse)

    if(type == "PP") {
        blocks1 <- list(nBlocks = inputs1$nBlocks * inputs1$nReplicates)
        blocks1$threshold <- inputs1$threshold
        blocks1$weights <- inputs1$weights
        if(length(blocks1$weights) != blocks1$nBlocks) blocks1$weights <- rep(blocks1$weights, inputs1$nReplicates)
        if(!is.null(inputs1$proportionMissing)) {
            if(length(inputs1$proportionMissing) != blocks1$nBlocks)
                inputs1$proportionMissing <- rep(inputs1$proportionMissing, inputs1$nReplicates)
            blocks1$proportionMissing <- inputs1$proportionMissing
        } else blocks1$proportionMissing <- 0
        if(!is.null(inputs1$x)) {
            if(inputs1$nReplicates == 1 || nrow(inputs1$x) == blocks1$nBlocks) {
                blocks1$data <- inputs1$x
            } else {
                                        # duplicate x with new rows for the replicates
                blocks1$data <- rep(inputs1$x[ , 1], inputs1$nReplicates)
                if(ncol(inputs1$x) > 1) {
                    for(k in 2:ncol(inputs1$x))
                        blocks1$data <- cbind(blocks1$data, rep(inputs1$x[ , k], inputs1$nReplicates))
                }
                names(blocks1$data) <- names(inputs1$x)
            }
        }
    
        blocks2 <- list(nBlocks = inputs2$nBlocks * inputs2$nReplicates)
        blocks2$threshold <- inputs2$threshold
        blocks2$weights <- inputs2$weights
        if(length(blocks2$weights) != blocks2$nBlocks) blocks2$weights <- rep(blocks2$weights, inputs2$nReplicates)

        if(!is.null(inputs2$proportionMissing)) {
            if(length(inputs2$proportionMissing) != blocks2$nBlocks)
                inputs2$proportionMissing <- rep(inputs2$proportionMissing, inputs2$nReplicates)
            blocks2$proportionMissing <- inputs2$proportionMissing
        } else blocks2$proportionMissing <- 0
        if(!is.null(inputs2$x)) {
            if(inputs2$nReplicates == 1 || nrow(inputs2$x) == blocks2$nBlocks) {
                blocks2$data <- inputs2$x
            } else {
                # duplicate x with new rows for the replicates
                blocks2$data <- rep(inputs2$x[ , 1], inputs2$nReplicates)
                if(ncol(inputs2$x) > 1) {
                    for(k in 2:ncol(inputs2$x))
                        blocks2$data <- cbind(blocks2$data, rep(inputs2$x[ , k], inputs2$nReplicates))
                }
                names(blocks2$data) <- names(inputs2$x)
            }
        }
    }
    
    usePhi1 <- !fit1$fit$const.scale
    usePhi2 <- !fit2$fit$const.scale

    # initial values?
    initial <- list()
    initial$location1 <- par1[grep("(location)|(mu)", names(par1))]
    if(!is.null(init1)) {
        initial$location1[1] <- init1[1]
        ln <- length(initial$location1)
        if(ln > 1) initial$location1[2:ln] <- 0
    }
    initial$scale1 <- par1[grep("(scale)|(phi)", names(par1))]
    if(!is.null(init1)) {
        if(usePhi1) init1[2] <- log(init1[2])
        initial$scale1[1] <- init1[2]
        ln <- length(initial$scale1)
        if(ln > 1) initial$scale1[2:ln] <- 0
    }
    initial$shape1 <- par1[grep("(shape)|(xi)", names(par1))]
    if(!is.null(init1)) {
        initial$shape1[1] <- init1[3]
        ln <- length(initial$shape1)
        if(ln > 1) initial$shape1[2:ln] <- 0
    }
    
                                        # set location coeffs to 0
                                        #    initial$location2 <- rep(0, p2[['location']] - 1)
    initial$location2 <- (par2[grep("(location)|(mu)", names(par2))])[-1]
    if(!is.null(init2)) {
        ln <- length(initial$location2)
        if(ln > 0) initial$location2[1:ln] <- 0
    }
    initial$scale2 <- par2[grep("(scale)|(phi)", names(par2))]
    if(!is.null(init2)) {
        if(usePhi2) init2[2] <- log(init2[2])
        initial$scale2[1] <- init2[2]
        ln <- length(initial$scale2)
        if(ln > 1) initial$scale2[2:ln] <- 0
    }
    initial$shape2 <- par2[grep("(shape)|(xi)", names(par2))]
    if(!is.null(init2)) {
        initial$shape2[1] <- init2[3]
        ln <- length(initial$shape2)
        if(ln > 1) initial$shape2[2:ln] <- 0
    }
    initial <- unlist(initial)
    names(initial) <- NULL

    if(type == "PP"){
        x1 <- inputs1$xObs
        x2 <- inputs2$xObs
        weights1 <- inputs1$weightsObs
        weights2 <- inputs2$weightsObs
    } else {
        x1 <- inputs1$x
        x2 <- inputs2$x
        weights1 <- inputs1$weights
        weights2 <- inputs2$weights
    }
    
    designs1 <- list()
    n <- length(inputs1$y)
    designs1$X.u <- setup.design(x = ~1, data = x1, n = n, dname = "thresholdFun")
    designs1$X.loc <- setup.design(x = inputs1$locationFun, data = x1, 
                                   n = n, const = fit1$fit$const.loc, dname = "locationFun")
    designs1$X.sc <- setup.design(x = inputs1$scaleFun, data = x1, 
                                  n = n, const = fit1$fit$const.scale, dname = "scaleFun")
    designs1$X.sh <- setup.design(x = inputs1$shapeFun, data = x1, 
                                  n = n, const = fit1$fit$const.shape, dname = "shapeFun")

    designs2 <- list()
    n <- length(inputs2$y)
    designs2$X.u <- setup.design(x = ~1, data = inputs2$xObs, n = n, dname = "thresholdFun")
    designs2$X.loc <- setup.design(x = inputs2$locationFun, data = x2, 
                                   n = n, const = fit2$fit$const.loc, dname = "locationFun")
    designs2$X.sc <- setup.design(x = inputs2$scaleFun, data = x2, 
                                  n = n, const = fit2$fit$const.scale, dname = "scaleFun")
    designs2$X.sh <- setup.design(x = inputs2$shapeFun, data = x2, 
                                  n = n, const = fit2$fit$const.shape, dname = "shapeFun")

    if(type == "PP") {
        blocks1$designs <- list()
#        if(is.element("data", names(blocks1))) {
            blocks1$X.u <- setup.design(x = ~1, 
                                        data = blocks1$data, n = blocks1$nBlocks, dname = "thresholdFun")
            blocks1$designs$X.loc <- setup.design(x = inputs1$locationFun, 
                                                  data = blocks1$data, n = blocks1$nBlocks, const = fit1$fit$const.loc, 
                                                  dname = "locationFun")
            blocks1$designs$X.sc <- setup.design(x = inputs1$scaleFun, 
                                                 data = blocks1$data, n = blocks1$nBlocks, const = fit1$fit$const.scale, 
                                                 dname = "scaleFun")
            blocks1$designs$X.sh <- setup.design(x = inputs1$shapeFun, 
                                                 data = blocks1$data, n = blocks1$nBlocks, const = fit1$fit$const.shape, 
                                                 dname = "shapeFun")
 #       }
        
        blocks2$designs <- list()
  #      if(is.element("data", names(blocks2))) {
            blocks2$X.u <- setup.design(x = ~1, 
                                        data = blocks2$data, n = blocks2$nBlocks, dname = "thresholdFun")
            blocks2$designs$X.loc <- setup.design(x = inputs2$locationFun, 
                                                  data = blocks2$data, n = blocks2$nBlocks, const = fit2$fit$const.loc, 
                                                  dname = "locationFun")
            blocks2$designs$X.sc <- setup.design(x = inputs2$scaleFun, 
                                                 data = blocks2$data, n = blocks2$nBlocks, const = fit2$fit$const.scale, 
                                                 dname = "scaleFun")
            blocks2$designs$X.sh <- setup.design(x = inputs2$shapeFun, 
                                                 data = blocks2$data, n = blocks2$nBlocks, const = fit2$fit$const.shape, 
                                                 dname = "shapeFun")
   #     }
    }
    
    rrHat <- exp(fit1$logReturnProb - fit2$logReturnProb)

    if(sum(unlist(p1)) > 3) {
        covariateMatrix <- stats::model.matrix(inputs1$locationFun, inputs1$xUse)
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(inputs1$scaleFun, inputs1$xUse))
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(inputs1$shapeFun, inputs1$xUse))
        qcov1 <- make.qcov(fit1$fit, covariateMatrix)
    } else qcov1 <- NULL
    if(sum(unlist(p2)) > 3) {
        covariateMatrix <- stats::model.matrix(inputs2$locationFun, inputs2$xUse)
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(inputs2$scaleFun, inputs2$xUse))
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(inputs2$shapeFun, inputs2$xUse))
        qcov2 <- make.qcov(fit2$fit, covariateMatrix)
    } else qcov2 <- NULL

    objfun <- function(logrr0, rrHat, logLikHat, cutoff, lowerEnd = TRUE, initial, y1, y2, x1, x2, threshold1, threshold2, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1, qcov2, type, optimArgs) {
        rr0 <- exp(logrr0)
        fitConstr <- stats::optim(initial, calc_constrNegLogLik, gr = NULL, rr0, y1, y2, x1, x2, threshold1, threshold2, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1, qcov2, type, method = optimArgs$method, lower = optimArgs$lower, upper = optimArgs$upper, control = optimArgs$control)$value
        # change to :
        # fitConstr <- stats::optim(initial, calc_constrNegLogLik, gr = NULL, rr0, y1, y2, xObs1, xObs2, thresholdObs1, thresholdObs2, weightsObs1, weightsObs2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1, qcov2, type = "PP", method = optimArgs$method, lower = optimArgs$lower, upper = optimArgs$upper, control = optimArgs$control)$value
        # check for lack of convergence
        if(lowerEnd) toggle <- 1 else toggle <- -1
        cond <- logLikHat + fitConstr < cutoff

        return(rr0*cond + toggle*1e10*(!cond) +
           # ensure that objective when constraint not satisfied gives unimodality
           (-toggle*rr0)*(!cond)*(rr0 < rrHat) + toggle*rr0*(!cond)*(rr0 > rrHat) )
    }

    nval <- ifelse(is.null(inputs1$mUse), 1, inputs1$mUse)
    output <- matrix(NA, nrow = nval, ncol = 2)
    for(i in seq_len(nval)) {
        # now embed 1-d minimization, on log(rr0) for hopefully better numerical performance
        cutoff <- 0.5*stats::qchisq(ciLevel, 1)
        if(rrHat[i] < bounds[1]) {
            lower <- NA
        } else {
            intvl <- log(c(bounds[1], min(rrHat[i], bounds[2])))
            lower <- stats::optimize(objfun, interval = intvl, rrHat = rrHat[i], logLikHat = logLikHat, cutoff = cutoff, lowerEnd = TRUE, initial = initial, inputs1$y, inputs2$y, x1, x2, inputs1$thresholdObs, inputs2$thresholdObs, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1[i, ], qcov2[i, ], type = type, optimArgs)$minimum
            if(all.equal(lower, intvl[1]) == TRUE || all.equal(lower, intvl[2]) == TRUE)
                warning("optimization likely failed to converge: lower endpoint equal to one endpoint of search interval.")
        }
        if(rrHat[i] > bounds[2]) {
            upper <- NA
        } else {
            intvl <- log(c(max(rrHat[i], bounds[1]), bounds[2]))
            upper <- stats::optimize(objfun, interval = intvl, rrHat = rrHat[i], logLikHat = logLikHat, cutoff = cutoff, lowerEnd = FALSE, initial = initial, inputs1$y, inputs2$y, inputs1$xObs, inputs2$xObs, inputs1$thresholdObs, inputs2$thresholdObs, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1[i, ], qcov2[i, ], type = type, optimArgs, maximum = TRUE)$maximum
            if(all.equal(upper, intvl[1]) == TRUE || all.equal(upper, intvl[2]) == TRUE)
                warning("optimization likely failed to converge: upper endpoint equal to one endpoint of search interval.")
        }
        output[i, ] <- exp(c(lower, upper))
    } # end loop over 1:inputs1$mUse
    return(drop(output))
} # end calc_riskRatio_lrt()

calc_constrNegLogLik <- function(par, rr0, y1, y2, x1, x2, threshold1 = NULL, threshold2 = NULL, weights1, weights2, p1, p2, designs1, designs2, blocks1 = NULL, blocks2 = NULL, usePhi1, usePhi2, returnValue, qcov1, qcov2, type = "PP"){
    np1 <- sum(unlist(p1))
    np2 <- sum(unlist(p2)) 
    par1 <- par[1:np1]
    par2 <- c(0, par[(np1+1):(np1+np2-1)])

    # calculate pA
    locationIndices <- 1:p1[['location']]
    scaleIndices <- (1 + p1[['location']]):(p1[['location']] + p1[['scale']])
    shapeIndices <- (1 + p1[['location']] + p1[['scale']]):np1

    if(np1 > 3) {
        location <- sum(qcov1[locationIndices] * par1[locationIndices])
        scale <- sum(qcov1[scaleIndices] * par1[scaleIndices])
        shape <- sum(qcov1[shapeIndices] * par1[shapeIndices])
    } else {
        location <- par1[locationIndices]
        scale <- par1[scaleIndices]
        shape <- par1[shapeIndices]
    }    
    if(usePhi1) scale <- exp(scale)
    if(scale <= 0) return(1e308)
    
    pA <- pevd(returnValue, location, scale, shape, lower.tail = FALSE, type = "GEV") # type always GEV here as PP params are in terms of equivalent GEV model

    # set mu0 based on pA, rr0
    locationIndices <- 1:p2[['location']]
    scaleIndices <- (1 + p2[['location']]):(p2[['location']] + p2[['scale']])
    shapeIndices <- (1 + p2[['location']] + p2[['scale']]):np2
    
    if(np2 > 3) {
        location <- sum(qcov2[locationIndices] * par2[locationIndices])
        scale <- sum(qcov2[scaleIndices] * par2[scaleIndices])
        shape <- sum(qcov2[shapeIndices] * par2[shapeIndices])
    } else {
        location <- par2[locationIndices]
        scale <- par2[scaleIndices]
        shape <- par2[shapeIndices]
    }
    if(usePhi2) scale <- exp(scale)
    if(scale <= 0) return(1e308)

    if(pA/rr0 > 1) return(1e308)
    
    par2[1] <- returnValue - location + (scale/shape) * (1 - (-log(1 - pA/rr0))^-shape)
    if(!is.finite(par2[1]))
        return(1e308)

    out1 <- out2 <- list(type = type)
    out1$weights <- weights1
    out2$weights <- weights2

    if(type == "PP") {
        result <- oevd(par1, out1, des = designs1, x = y1, data = x1, u = threshold1, npy = 0, phi = usePhi1, blocks = blocks1) +
            oevd(par2, out2, des = designs2, x = y2, data = x2, u = threshold2, npy = 0, phi = usePhi2, blocks = blocks2)
    } else if(type == "GEV") {
        result <- oevd(par1, out1, des = designs1, x = y1, data = x1, npy = 0, phi = usePhi1) +
            oevd(par2, out2, des = designs2, x = y2, data = x2, npy = 0, phi = usePhi2)
    } else stop("calc_constrNegLogLik: 'type' must be 'PP' or 'GEV'.")
    return(result)
} # end calc_constrNegLogLik()
