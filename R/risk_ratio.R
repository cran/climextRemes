bootTypes <- list('boot_norm' = 'norm', 'boot_perc' = 'perc', 'boot_basic' = 'basic', 'boot_stud' = 'stud', 'boot_bca' = 'bca')
bootTypesCols <- list('boot_norm' = 'normal', 'boot_perc' = 'percent', 'boot_basic' = 'basic', 'boot_stud' = 'student', 'boot_bca' = 'bca')

#' Compute risk ratio and uncertainty based on binomial models for counts of events relative to possible number of events
#'
#' Compute risk ratio and uncertainty by fitting binomial models to counts of events relative to possible number of events. The risk ratio is the ratio of the probability of an event under the model fit to the first dataset to the probability under the model fit to the second dataset. Default standard errors are based on the usual MLE asymptotics using a delta-method-based approximation, but standard errors based on the nonparametric bootstrap and on a likelihood ratio procedure can also be computed.
#'
#' @name calc_riskRatio_binom
#' 
#' @param y vector of two values, the number of events in the two scenarios
#' @param n vector of two values, the number of samples (possible occurrences of events) in the two scenarios 
#' @param ciLevel statistical confidence level for confidence intervals; in repeated experimentation, this proportion of confidence intervals should contain the true risk ratio. Note that if only one endpoint of the resulting interval is used, for example the lower bound, then the effective confidence level increases by half of one minus \code{ciLevel}. For example, a two-sided 0.90 confidence interval corresponds to a one-sided 0.95 confidence interval.
#' @param ciType character vector indicating which type of confidence intervals to compute. See \code{Details}.
#' @param bootSE logical indicating whether to use the bootstrap to estimate the standard error of the risk ratio
#' @param bootControl a list of control parameters for the bootstrapping, used only when at least one bootstrap confidence interval is requested via \code{ciType}. See \code{Details}.
#' @param lrtControl list containing a single component, \code{bounds}, which sets the range inside which the algorithm searches for the endpoints of the likelihood ratio-based confidence interval. This avoids numerical issues with endpoints converging to zero and infinity. If an endpoint is not found within the interval, it is set to \code{NA}. Used only when \code{'lrt'} is one of the \code{ciType} values. Default is (0.01, 100).
#' @author Christopher J. Paciorek
#' @export
#' @return
#'
#' The primary outputs of this function are as follows: the log of the risk ratio and standard error of that log risk ratio (\code{logRiskRatio} and \code{se_logRiskRatio}) as well the risk ratio itself (\code{riskRatio}). The standard error is based on the usual MLE asymptotics using a delta-method-based approximation. If requested via \code{ciType}, confidence intervals will be returned, as discussed in \code{Details}.
#'
#' @details
#' \code{ciType} can include one or more of the following: \code{'delta'}, \code{'koopman'}, \code{'lrt'}, \code{'boot_norm'}, \code{'boot_perc'}, \code{'boot_basic'}, \code{'boot_stud'}, \code{'boot_bca'}. \code{'delta'} uses the delta method to compute an asymptotic interval based on the standard error of the log risk ratio. \code{'koopman'} uses the method described in Koopman (1984), following the implementation discussed in Fageland et al. (2015), including the calculation of Nam (1995). \code{'lrt'} inverts a likelihood-ratio test. Bootstrap-based options are the normal-based interval using the bootstrap standard error (\code{'boot_norm'}), the percentile bootstrap (\code{'boot_perc'}), the basic bootstrap (\code{'boot_basic'}), the bootstrap-t (\code{'boot_stud'}), and the bootstrap BCA method (\code{'boot_bca'}). See Paciorek et al. for more details. 
#' 
#' See \code{\link{fit_pot}} for information on the \code{bootControl} argument.
#' @references
#' Paciorek, C.J., D.A. Stone, and M.F. Wehner. 2018. Quantifying uncertainty in the attribution of human influence on severe weather. Weather and Climate Extremes 20:69-80. arXiv preprint <https://arxiv.org/abs/1706.03388>.
#' 
#' Koopman, P.A.R. 1984. Confidence intervals for the ratio of two binomial proportions. Biometrics 40: 513-517.
#'
#' Fagerland, M.W., S. Lydersen, and P. Laake. 2015. Recommended confidence intervals for two independent binomial proportions. Statistical Methods in Medical Research 24: 224-254.
#' @examples
#' # risk ratio for 40/400 compared to 8/400 events and for
#' # 4/100 compared to 0/100 events
#' calc_riskRatio_binom(c(40, 8), c(400, 400), ciType = c('lrt', 'boot_stud', 'koopman'))
#' # LRT and Koopman methods can estimate lower confidence interval endpoint
#' # even if estimated risk ratio is infinity:
#' calc_riskRatio_binom(c(4,0), c(100, 100), ciType = c('lrt', 'boot_stud', 'koopman'))
calc_riskRatio_binom <- function(y, n, ciLevel = 0.90,
                                 ciType, bootSE,
                                 bootControl = NULL, lrtControl = NULL) {

    if(is.null(y) || is.null(n))  # for Python interface
        stop("calc_riskRatio_binom: argument 'y' or 'n' is missing, with no default")

    if(missing(bootSE) || is.null(bootSE))
        bootSE <- !missing(ciType) && sum(ciType %in% names(bootTypes))
    
    if(missing(ciType) || is.null(ciType)) {
        ciType <- ""
    } else {
        wh <- setdiff(ciType, c('delta', 'lrt', 'koopman', names(bootTypes)))
        if(length(wh))
            stop("calc_riskRatio_binom: ", paste(wh, collapse = ','), " in 'ciType' are not valid types.")
    }

    if(is.null(lrtControl))
        lrtControl = list(bounds = c(.01, 100))

    if(is.null(y) || is.null(n))  # for Python interface
        stop("calc_riskRatio_binom: argument 'y' or 'n' is missing, with no default")
    
    results <- list()
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE)

    if(length(n) == 1) n <- rep(n, 2)
    results$logRiskRatio <- log(y[1]) - log(n[1]) - log(y[2]) + log(n[2])
    p <- y / n
    results$se_logRiskRatio <- sqrt(calc_delta_var(y, n))

    results$riskRatio <- exp(results$logRiskRatio)
    if('delta' %in% ciType)
        results$ci_riskRatio_delta <- exp(results$logRiskRatio + c(-1, 1)*z_alpha*results$se_logRiskRatio)
    
    if('lrt' %in% ciType) 
        results$ci_riskRatio_lrt <- calc_riskRatio_lrt_binom(y, n, ciLevel = ciLevel, bounds = lrtControl$bounds)
    
    if('koopman' %in% ciType) 
        results$ci_riskRatio_koopman <- calc_riskRatio_koopman_binom(y, n, ciLevel = ciLevel)
    
    bootTypesUse <- ciType[ciType %in% names(bootTypes)]
    if(length(bootTypesUse)) {
        bControl <- list(seed = 1, n = 250, by = "block", getSample = FALSE)
        bControl[names(bootControl)] <- bootControl
        bootControl <- bControl
        
        if(any(p == 0 | p == 1)) {
            results$se_logRiskRatio_boot <- NA
            for(type in bootTypesUse) {
                typeName <- paste0('ci_riskRatio_', type)
                results[[typeName]] <- rep(NA, 2)
            }
        } else {
            if(length(bootControl$seed) == 1){
                set.seed(bootControl$seed)
            } else{
                .Random.seed <- bootControl$seed
            }
            bootData <- cbind(stats::rbinom(n = bootControl$n, size = n[1], prob = p[1]),
                              stats::rbinom(n = bootControl$n, size = n[2], prob = p[2]))

            # need an instance of class boot to use in boot.ci
            fake_data <- cbind(sample(c(0,1), size = 5, replace = TRUE),
                   sample(c(0,1), size = 5, replace = TRUE))
            logRRfun <- function(dat, ind) {
                log(sum(dat[ind,1])/sum(dat[ind,2]))
            }
            bootInput <- boot::boot(fake_data, logRRfun, R = bootControl$n)
            
            bootInput$t0 = c(log(p[1]) - log(p[2]), calc_delta_var(y, n))
            logRRvals <- log(bootData[ , 1]) - log(bootData[ , 2])
            if(bootSE)
                results$se_logRiskRatio_boot <- sd(logRRvals)
            logRRvals[logRRvals == Inf] = 1e6 # boot.ci can't handle Inf
            bootInput$t <- cbind(logRRvals, apply(bootData, 1, calc_delta_var, n = n))
            
            bootResults <- try(boot.ci(bootInput, conf = ciLevel, type = unlist(bootTypes[bootTypesUse])))
            if(is(bootResults, 'try-error')) {
                for(type in bootTypesUse) {
                    typeName <- paste0('ci_riskRatio_', type)
                    results[[typeName]] <- rep(NA, 2)
                }
            } else {
                for(type in bootTypesUse) {
                    typeName <- paste0('ci_riskRatio_', type)
                    if(type == 'boot_norm') {
                        results[[typeName]] <- exp(bootResults[[bootTypesCols[[type]]]][2:3])
                    } else results[[typeName]] <- exp(bootResults[[bootTypesCols[[type]]]][4:5])
                }
            }
            if(bootControl$getSample) results$riskRatio_boot <- bootData[ , 1] / bootData[ , 2]
        }
    }
    return(results)
} # end calc_riskRatio_binom()

calc_delta_var <- function(y, n) {
    return( (1 - y[1]/n[1]) / y[1] + (1 - y[2]/n[2]) / y[2])
}

calc_riskRatio_koopman_binom <- function(y, n, ciLevel = 0.90) {
                                        # from Koopman (1984) with exact solution from Nam (1995)
    koopman_aux <-  function(p0, yA, yN, nA, nN) {
        (1-(nA-yA)*(1-p0)/(yN+nA-(nA+nN)*p0))/p0
    }
    yA <- y[1]
    yN <- y[2]
    nA <- n[1]
    nN <- n[2]
    n <- sum(n)
    if(yA ==0 && yN == 0) return(c(0, Inf))
    za <- qchisq(ciLevel, 1)
    a1 <- nN*(nN*n*yA+nA*(nN+yA)*za)
    a2 <- -nN*(nN*nA*(yA+yN)+2*n*yA*yN+nA*(nN+yN+2*yA)*za)
    a3 <- 2*nA*nN*yN*(yA+yN)+n*yN^2*yA+nA*nN*(yA+yN)*za
    a4 <- -nA*yN^2*(yA+yN)
    b1 <- a2/a1
    b2 <- a3/a1
    b3 <- a4/a1
    c1 <- b2-b1^2/3
    c2 <- b3-b1*b2/3 + 2*b1^3/27
    theta <- acos(sqrt(27)*c2/(2*c1*sqrt(-c1)))
    t1 <- -2*sqrt(-c1/3)*cos(pi/3-theta/3)
    t2 <- -2*sqrt(-c1/3)*cos(pi/3+theta/3)
    t3 <- 2*sqrt(-c1/3)*cos(theta/3)
    p0 <- c(t1,t2,t3)- b1/3
    p0 <- sort(p0)
    L <- koopman_aux(p0[2], yA, yN, nA, nN)
    U <- koopman_aux(p0[1], yA, yN, nA, nN)
    if(yA == 0)
        L <- 0
    if(yN == 0)
        U <- Inf
    return(c(L, U))
} # end calc_riskRatio_koopman_binom

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
    
    objfun <- function(logrr, y, n, logLikHat, cutoff) {
        rr <- exp(logrr)
        pAtilde <- restricted_mle(y/n, rr)
        logLikConstr <- logLik(y, n, c(pAtilde, pAtilde/rr))
        if(pAtilde > rr || pAtilde > 1) 
            stop('calc_riskRatio_lrt_binom: infeasible value for probability for one of the groups.')
        output <- logLikHat - logLikConstr - cutoff
        attributes(output) <- NULL  ## uniroot has trouble extending interval when attributes present
        return(output)
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
        lower <- stats::uniroot(objfun, interval = intvl, y = y, n = n, logLikHat = logLikHat, cutoff = cutoff, extendInt = 'downX')$root
        if(all.equal(lower, intvl[1]) == TRUE || all.equal(lower, intvl[2]) == TRUE)
            warning("optimization likely failed to converge: lower endpoint equal to one endpoint of search interval.")
        ## however note that with use of uniroot, it will extend the interval
    }
    if(pHat[2] == 0) {
        upper <- Inf
    } else {
        intvl <- log(c(max(rrHat, bounds[1]), bounds[2]))
        upper <- stats::uniroot(objfun, interval = intvl, y = y, n = n, logLikHat = logLikHat, cutoff = cutoff, extendInt = 'upX')$root
        if(all.equal(upper, intvl[1]) == TRUE || all.equal(upper, intvl[2]) == TRUE)
            warning("optimization likely failed to converge: upper endpoint equal to one endpoint of search interval.")
    }
    return(exp(c(lower, upper)))
} # end calc_riskRatio_lrt_binom()


#' Compute risk ratio and uncertainty based on peaks-over-threshold models fit to exceedances over a threshold
#'
#' Compute risk ratio and uncertainty by fitting peaks-over-threshold model, designed specifically for climate data, to exceedance-only data, using the point process approach. The risk ratio is the ratio of the probability of exceedance of a pre-specified value under the model fit to the first dataset to the probability under the model fit to the second dataset. Default standard errors are based on the usual MLE asymptotics using a delta-method-based approximation, but standard errors based on the nonparametric bootstrap and on a likelihood ratio procedure can also be computed.
#'
#' @name calc_riskRatio_pot
#'
#' @param returnValue numeric value giving the value for which the risk ratio should be calculated, where the resulting period will be the average number of blocks until the value is exceeded and the probability the probability of exceeding the value in any single block.
#' @param y1 a numeric vector of exceedance values for the first dataset (values of the outcome variable above the threshold). For better optimization performance, it is recommended that the \code{y1} have magnitude around one (see \code{Details}), for which one can use \code{scaling1}.
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
#' @param index1 numeric vector providing the integer-valued index (e.g., julian day for daily climate data) corresponding to each element of \code{y1}. For example if there are 10 original observations and the third, fourth, and seventh values are exceedances, then \code{index1} would be the vector 3,4,7. Used only when \code{declustering} is provided to determine which exceedances occur sequentially or within a contiguous set of values of a given length. The actual values are arbitrary; only the lags between the values are used.
#' @param index2 numeric vector providing the integer-valued index (e.g., julian day for daily climate data) corresponding to each element of \code{y2}. Analogous to \code{index1}.
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
#' @param ciType character vector indicating which type of confidence intervals to compute. See \code{Details}.
#' @param bootSE logical indicating whether to use the bootstrap to estimate the standard error of the risk ratio
#' @param bootControl a list of control parameters for the bootstrapping. See \code{Details}.
#' @param lrtControl list containing a single component, \code{bounds}, which sets the range inside which the algorithm searches for the endpoints of the likelihood ratio-based confidence interval. This avoids numerical issues with endpoints converging to zero and infinity. If an endpoint is not found within the interval, it is set to \code{NA}.
#' @param optimArgs a list with named components matching exactly any arguments that the user wishes to pass to \code{optim}. See \code{help(optim)} for details. Of particular note, \code{'method'} can be used to choose the optimization method used for maximizing the log-likelihood to fit the model and \code{'control=list(maxit=VALUE)'} for a user-chosen VALUE can be used to increase the number of iterations if the optimization is converging slowly.
#' @param optimControl a list with named components matching exactly any elements that the user wishes to pass as the \code{control} argument to R's \code{optim} function. See \code{help(optim)} for details. Primarily provided for the Python interface because \code{control} can also be passed as part of \code{optimArgs}.
#' @param initial1 a list with components named \code{'location'}, \code{'scale'}, and \code{'shape'} providing initial parameter values for the first dataset, intended for use in speeding up or enabling optimization when the default initial values are resulting in failure of the optimization; note that use of \code{scaling1}, \code{logScale1} and \code{.normalizeX = TRUE} cause numerical changes in some of the parameters. For example with \code{logScale1 = TRUE}, initial value(s) for \code{'scale'} should be specified on the log scale.
#' @param initial2 a list with components named \code{'location'}, \code{'scale'}, and \code{'shape'} providing initial parameter values for the second dataset, intended for use in speeding up or enabling optimization when the default initial values are resulting in failure of the optimization; note that use of \code{scaling2}, \code{logScale2} and \code{.normalizeX = TRUE} cause numerical changes in some of the parameters. For example with \code{logScale2 = TRUE}, initial value(s) for \code{'scale'} should be specified on the log scale.
#' @param logScale1 logical indicating whether optimization for the scale parameter should be done on the log scale for the first dataset. By default this is FALSE when the scale is not a function of covariates and TRUE when the scale is a function of covariates (to ensure the scale is positive regardless of the regression coefficients). 
#' @param logScale2 logical indicating whether optimization for the scale parameter should be done on the log scale for the second dataset. By default this is FALSE when the scale is not a function of covariates and TRUE when the scale is a function of covariates (to ensure the scale is positive regardless of the regression coefficients). 
#' @param getReturnCalcs logical indicating whether to return the estimated return values/probabilities/periods from the fitted models. 
#' @param getParams logical indicating whether to return the fitted parameter values and their standard errors for the fitted models; WARNING: parameter values for models with covariates for the scale parameter must interpreted based on the value of \code{logScale}.
#' @param getFit logical indicating whether to return the full fitted models (potentially useful for model evaluation and for understanding optimization problems); note that estimated parameters in the fit object for nonstationary models will not generally match the MLE provided when \code{getParams} is \code{TRUE} because covariates are normalized before fitting and the fit object is based on the normalized covariates. Similarly, parameters will not match if \code{scaling} is not 1.
#' @author Christopher J. Paciorek
#' @export
#' @details
#' See \code{\link{fit_pot}} for more details on fitting the peaks-over-threshold model for each dataset, including details on blocking and replication. Also see \code{\link{fit_pot}} for information on the \code{bootControl} argument. 
#'
#' Optimization failures:
#'
#' It is not uncommon for maximization of the log-likelihood to fail for extreme value models. Please see the help information for \code{fit_pot}. Also note that if the probability in the denominator of the risk ratio is near one, one may achieve better numerical performance by swapping the two datasets and computing the risk ratio for the probability under dataset 2 relative to the probability under dataset 1.
#' 
#' @return
#'
#' The primary outputs of this function are as follows: the log of the risk ratio and standard error of that log risk ratio (\code{logRiskRatio} and \code{se_logRiskRatio}) as well the risk ratio itself (\code{riskRatio}). The standard error is based on the usual MLE asymptotics using a delta-method-based approximation. If requested via \code{ciType}, confidence intervals will be returned, as discussed in \code{Details}.
#'
#' @details
#' \code{ciType} can include one or more of the following: \code{'delta'}, \code{'lrt'}, \code{'boot_norm'}, \code{'boot_perc'}, \code{'boot_basic'}, \code{'boot_stud'}, \code{'boot_bca'}. \code{'delta'} uses the delta method to compute an asymptotic interval based on the standard error of the log risk ratio. \code{'lrt'} inverts a likelihood-ratio test. Bootstrap-based options are the normal-based interval using the bootstrap standard error (\code{'boot_norm'}), the percentile bootstrap (\code{'boot_perc'}), the basic bootstrap (\code{'boot_basic'}), the bootstrap-t (\code{'boot_stud'}), and the bootstrap BCA method (\code{'boot_bca'}). See Paciorek et al. for more details. 
#' 
#' See \code{\link{fit_pot}} for information on the \code{bootControl} argument. 
#' @references
#' Paciorek, C.J., D.A. Stone, and M.F. Wehner. 2018. Quantifying uncertainty in the attribution of human influence on severe weather. Weather and Climate Extremes 20:69-80. arXiv preprint <https://arxiv.org/abs/1706.03388>.
#'
#' Jeon S., C.J. Paciorek, and M.F. Wehner. 2016. Quantile-based bias correction and uncertainty quantification of extreme event attribution statements. Weather and Climate Extremes 12: 24-32. <DOI:10.1016/j.wace.2016.02.001>. arXiv preprint: <http://arxiv.org/abs/1602.04139>.
#' 
#' @examples
#' data(Fort, package = 'extRemes')
#' threshold <- 0.395
#' ord <- order(Fort$year, Fort$month, Fort$day) 
#' Fort <- Fort[ord, ]
#' ind <- Fort$Prec > threshold
#' FortExc <- Fort[ind, ]
#' earlyYears <- 1900:1929
#' lateYears <- 1970:1999
#' earlyPeriod <- which(FortExc$year %in% earlyYears)
#' latePeriod <- which(FortExc$year %in% lateYears)
#' # contrast late period with early period, assuming a nonstationary fit
#' # within each time period and finding RR based on midpoint of each period
#' \dontrun{
#' out <- calc_riskRatio_pot(returnValue = 3,
#'                    y1 = FortExc$Prec[earlyPeriod], y2 = FortExc$Prec[latePeriod],
#'                    x1 = data.frame(years = earlyYears), x2 = data.frame(years = lateYears),
#'                    threshold1 = threshold, threshold2 = threshold,
#'                    locationFun1 = ~years, locationFun2 = ~years,
#'                    xNew1 = data.frame(years = mean(earlyYears)), 
#'		      xNew2 = data.frame(years = mean(lateYears)),
#'                    blockIndex1 = FortExc$year[earlyPeriod], 
#'                    blockIndex2 = FortExc$year[latePeriod],
#'                    firstBlock1 = earlyYears[1], firstBlock2 = lateYears[1])
#' }
calc_riskRatio_pot <- function(returnValue, y1, y2, x1 = NULL, x2 = x1,
                               threshold1, threshold2 = threshold1, locationFun1 = NULL,
                               locationFun2 = locationFun1, scaleFun1 = NULL, scaleFun2 = scaleFun1,
                               shapeFun1 = NULL, shapeFun2 = shapeFun1, nBlocks1 = nrow(x1),
                               nBlocks2 = nrow(x2), blockIndex1 = NULL, blockIndex2 = NULL, firstBlock1 = 1,
                               firstBlock2 = 1, index1 = NULL, index2 = NULL, nReplicates1 = 1,
                               nReplicates2 = 1, replicateIndex1 = NULL, replicateIndex2 = NULL,
                               weights1 = NULL, weights2 = NULL, proportionMissing1 = NULL,
                               proportionMissing2 = NULL, xNew1 = NULL, xNew2 = NULL, declustering = NULL,
                               upperTail = TRUE, scaling1 = 1, scaling2 = 1, ciLevel = 0.90, ciType,
                               bootSE, bootControl = NULL, lrtControl = NULL, 
                               optimArgs = NULL, optimControl = NULL, initial1 = NULL,
                               initial2 = NULL, logScale1 = NULL, logScale2 = NULL,
                               getReturnCalcs = FALSE, getParams = FALSE, getFit = FALSE  ) {

    if(is.null(returnValue) || is.null(y1) || is.null(y2))  # for Python interface
        stop("calc_riskRatio_binom: argument 'returnValue' or 'y1' or 'y2' is missing, with no default")

    if(missing(bootSE) || is.null(bootSE))
        bootSE <- !missing(ciType) && sum(ciType %in% names(bootTypes))
    if(missing(returnValue))
        stop("calc_riskRatio_pot: 'returnValue' must be provided.")
    if(length(returnValue) > 1)
        stop("calc_riskRatio_pot: please provide a single value for 'returnValue'.")
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE)
    ciLabels <- as.character(c(( 1 - ciLevel ) / 2, 1 - ( 1 - ciLevel ) / 2))
    if(!is.null(xNew1)) {
        xNew1tmp <- try(as.data.frame(xNew1))
        if(is(xNew1tmp, 'try-error')) stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
        m <- nrow(xNew1tmp)
    } else {
        if(is.null(x1))
            m <- 1 else {
                       x1tmp <- try(as.data.frame(x1))
                       if(is(x1tmp, 'try-error')) stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
                       m <- nrow(x1tmp)
                   }
    }
    results <- list()

    if(missing(ciType)) {
        ciType <- ""
    } else {
        wh <- setdiff(ciType, c('delta', 'lrt', 'koopman', names(bootTypes)))
        if(length(wh))
            stop("calc_riskRatio_pot: ", paste(ciType[wh], collapse = ','), " in 'ciType' are not valid types.")
    }

    bootTypesUse <- ciType[ciType %in% names(bootTypes)]
    if(length(bootTypesUse)) {
        bControl <- list(seed = 1, n = 250, by = "block", getSample = FALSE)
        bControl[names(bootControl)] <- bootControl
        bootControl <- bControl
        bootControlTmp <- bootControl
        bootControlTmp$getSample <- TRUE
        bootControlTmp$getSampleSE <- TRUE
    } 
    
    fit1 <- fit_pot(y1, x = x1, threshold = threshold1, locationFun = locationFun1,
                    scaleFun = scaleFun1, shapeFun = shapeFun1, nBlocks = nBlocks1,
                    blockIndex = blockIndex1, firstBlock = firstBlock1, index = index1,
                    nReplicates = nReplicates1, replicateIndex = replicateIndex1,
                    weights = weights1, proportionMissing = proportionMissing1, returnValue = returnValue,
                    xNew = xNew1, declustering = declustering, upperTail = upperTail,
                    scaling = scaling1, bootSE = bootSE, bootControl = bootControlTmp,
                    optimArgs = optimArgs, optimControl = optimControl, initial = initial1, logScale = logScale1,
                    getFit = TRUE, getParams = getParams, .getInputs = TRUE) 
    fit2 <- fit_pot(y2, x = x2, threshold = threshold2, locationFun = locationFun2,
                    scaleFun = scaleFun2, shapeFun = shapeFun2, nBlocks = nBlocks2,
                    blockIndex = blockIndex2, firstBlock = firstBlock2, index = index2,
                    nReplicates = nReplicates2, replicateIndex = replicateIndex2,
                    weights = weights2, proportionMissing = proportionMissing2, returnValue = returnValue,
                    xNew = xNew2, declustering = declustering, upperTail = upperTail,
                    scaling = scaling2, bootSE = bootSE, bootControl = bootControlTmp,
                    optimArgs = optimArgs, optimControl = optimControl, initial = initial2, logScale = logScale2,
                    getFit = TRUE, getParams = getParams, .getInputs = TRUE)
    if(fit1$info$failure || fit2$info$failure) {
        warning("calc_riskRatio_pot: fitting failed for one of two datasets.")
        results$logRiskRatio <- results$se_logRiskRatio <- results$riskRatio <- rep(NA, m)
        for(type in ciType)
            results[[paste0('ci_riskRatio_', type)]] <- drop(matrix(NA, m, 2))
        if(bootSE) 
            results$se_logRiskRatio_boot <- rep(NA, m)
    } else {
        if(length(fit1$logReturnProb) != length(fit2$logReturnProb))
            stop("calc_riskRatio_pot: number of return probabilities calculated for each model fit must be the same; for nonstationary models this is determined by the number of covariate set inputs (provided in 'x' or 'xNew').")
        results$logRiskRatio <- fit1$logReturnProb - fit2$logReturnProb
        
        results$se_logRiskRatio <- sqrt(fit1$se_logReturnProb^2 + fit2$se_logReturnProb^2)
        results$riskRatio <- exp(results$logRiskRatio)
        if('delta' %in% ciType) {
            results$ci_riskRatio_delta <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio,
                                                    results$logRiskRatio + z_alpha*results$se_logRiskRatio))
            colnames(results$ci_riskRatio_delta) <- ciLabels
            results$ci_riskRatio_delta <- drop(results$ci_riskRatio_delta)
        }
        
        if(length(bootTypesUse)) {
            bootData <- matrix(drop(fit1$logReturnProb_boot - fit2$logReturnProb_boot), ncol = m)
            if(bootSE)
                results$se_logRiskRatio_boot <- apply(bootData, 2, sd, na.rm = TRUE)
            
            fake_data <- cbind(sample(c(0,1), size = 5, replace = TRUE),
                   sample(c(0,1), size = 5, replace = TRUE))
            logRRfun <- function(dat, ind) {
                log(sum(dat[ind,1])/sum(dat[ind,2]))
            }
            bootInput <- boot::boot(fake_data, logRRfun, R = bootControl$n)

            for(type in bootTypesUse) {
                typeName <- paste0('ci_riskRatio_', type)
                results[[typeName]] <- matrix(NA, m, 2)
                colnames(results[[typeName]]) <- ciLabels
            }
            for(i in seq_len(m)) {
                bootInput$t0 = c(results$logRiskRatio[i], results$se_logRiskRatio[i]^2)
                logRRvals <- bootData[ , i]
                logRRvals[logRRvals == Inf] = 1e6 # boot.ci can't handle Inf
                bootInput$t <- cbind(logRRvals,
                                     fit1$logReturnProb_boot_se[ , i, 1]^2 + fit2$logReturnProb_boot_se[ , i, 1]^2)
                bootResults <- try(boot::boot.ci(bootInput, conf = ciLevel, type = unlist(bootTypes[bootTypesUse])))
                if(!is(bootResults, 'try-error')) {
                    for(type in bootTypesUse) {
                        typeName <- paste0('ci_riskRatio_', type)
                        if(type == 'boot_norm') {
                            results[[typeName]][i, ] <- exp(bootResults[[bootTypesCols[[type]]]][2:3])
                        } else results[[typeName]][i, ] <- exp(bootResults[[bootTypesCols[[type]]]][4:5])
                        if(i == m) results[[typeName]] <- drop(results[[typeName]])
                    }
                }
            }
            if(bootControl$getSample) results$riskRatio_boot <- drop(exp(bootData))
        }
        
        if('lrt' %in% ciType) {
            lControl <- list(bounds = c(.01, 100))
            lControl[names(lrtControl)] <- lrtControl
            lrtControl <- lControl
            oArgs <- list(method = "Nelder-Mead", lower = -Inf, upper = Inf, control = list())
            oArgs[names(optimArgs)] <- optimArgs
            if(!is.null(optimControl))
                oArgs[['control']] <- optimControl
            if(!upperTail) returnValue <- -returnValue
            results$ci_riskRatio_lrt <- calc_riskRatio_lrt(fit1, fit2, returnValue, ciLevel = ciLevel, bounds = lrtControl$bounds, type = "PP", optimArgs = oArgs)
            if(is.null(dim(results$ci_riskRatio_lrt))) {
                names(results$ci_riskRatio_lrt) <- ciLabels
            } else colnames(results$ci_riskRatio_lrt) <- ciLabels
        }
    }
    if(getFit || getParams || getReturnCalcs) {
        fit1$inputs <- fit2$inputs <- NULL
        if(!getFit) {
            fit1$fit <- fit2$fit <- fit1$info <- fit2$info <- NULL
        }
        if(!getReturnCalcs) {
            fit1$logReturnProb <- fit1$se_logReturnProb <- fit1$logReturnPeriod <- fit1$se_logReturnPeriod <- NULL
            fit2$logReturnProb <- fit2$se_logReturnProb <- fit2$logReturnPeriod <- fit2$se_logReturnPeriod <- NULL
            if(bootSE) {
                fit1$se_logReturnProb_boot <- fit1$se_logReturnPeriod <- NULL
                fit2$se_logReturnProb_boot <- fit2$se_logReturnPeriod <- NULL
            }
        }
        results$fit1 <- fit1
        results$fit2 <- fit2
    }
    return(results)
} # end calc_riskRatio_pot

#' Compute risk ratio and uncertainty based on generalized extreme value model fit to block maxima or minima
#'
#' Compute risk ratio and uncertainty by fitting generalized extreme value model, designed specifically for climate data, to exceedance-only data, using the point process approach. The risk ratio is the ratio of the probability of exceedance of a pre-specified value under the model fit to the first dataset to the probability under the model fit to the second dataset. Default standard errors are based on the usual MLE asymptotics using a delta-method-based approximation, but standard errors based on the nonparametric bootstrap and on a likelihood ratio procedure can also be computed.
#'
#' @name calc_riskRatio_gev
#' 
#' @param returnValue numeric value giving the value for which the risk ratio should be calculated, where the resulting period will be the average number of blocks until the value is exceeded and the probability the probability of exceeding the value in any single block.
#' @param y1 a numeric vector of observed maxima or minima values for the first dataset. See \code{Details} for how the values of \code{y1} should be ordered if there are multiple replicates and the values of \code{x1} are identical for all replicates. For better optimization performance, it is recommended that \code{y1} have magnitude around one (see \code{Details}), for which one can use \code{scaling1}.
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
#' @param ciType character vector indicating which type of confidence intervals to compute. See \code{Details}.
#' @param bootSE logical indicating whether to use the bootstrap to estimate the standard error of the risk ratio
#' @param bootControl a list of control parameters for the bootstrapping. See \code{Details}.
#' @param lrtControl list containing a single component, \code{bounds}, which sets the range inside which the algorithm searches for the endpoints of the likelihood ratio-based confidence interval. This avoids numerical issues with endpoints converging to zero and infinity. If an endpoint is not found within the interval, it is set to \code{NA}.
#' @param optimArgs a list with named components matching exactly any arguments that the user wishes to pass to \code{optim}. See \code{help(optim)} for details. Of particular note, \code{'method'} can be used to choose the optimization method used for maximizing the log-likelihood to fit the model and \code{'control=list(maxit=VALUE)'} for a user-chosen VALUE can be used to increase the number of iterations if the optimization is converging slowly.
#' @param optimControl a list with named components matching exactly any elements that the user wishes to pass as the \code{control} list to R's \code{optim} function. See \code{help(optim)} for details. Primarily provided for the Python interface because \code{control} can also be passed as part of \code{optimArgs}.
#' @param initial1 a list with components named \code{'location'}, \code{'scale'}, and \code{'shape'} providing initial parameter values for the first dataset, intended for use in speeding up or enabling optimization when the default initial values are resulting in failure of the optimization; note that use of \code{scaling1}, \code{logScale1} and \code{.normalizeX = TRUE} cause numerical changes in some of the parameters. For example with \code{logScale1 = TRUE}, initial value(s) for \code{'scale'} should be specified on the log scale.
#' @param initial2 a list with components named \code{'location'}, \code{'scale'}, and \code{'shape'} providing initial parameter values for the second dataset, intended for use in speeding up or enabling optimization when the default initial values are resulting in failure of the optimization; note that use of \code{scaling2}, \code{logScale2} and \code{.normalizeX = TRUE} cause numerical changes in some of the parameters. For example with \code{logScale2 = TRUE}, initial value(s) for \code{'scale'} should be specified on the log scale.
#' @param logScale1 logical indicating whether optimization for the scale parameter should be done on the log scale for the first dataset. By default this is FALSE when the scale is not a function of covariates and TRUE when the scale is a function of covariates (to ensure the scale is positive regardless of the regression coefficients). 
#' @param logScale2 logical indicating whether optimization for the scale parameter should be done on the log scale for the second dataset. By default this is FALSE when the scale is not a function of covariates and TRUE when the scale is a function of covariates (to ensure the scale is positive regardless of the regression coefficients). 
#' @param getReturnCalcs logical indicating whether to return the estimated return values/probabilities/periods from the fitted models. 
#' @param getParams logical indicating whether to return the fitted parameter values and their standard errors for the fitted models; WARNING: parameter values for models with covariates for the scale parameter must interpreted based on the value of \code{logScale}.
#' @param getFit logical indicating whether to return the full fitted models (potentially useful for model evaluation and for understanding optimization problems); note that estimated parameters in the fit object for nonstationary models will not generally match the MLE provided when \code{getParams} is \code{TRUE} because covariates are normalized before fitting and the fit object is based on the normalized covariates. Similarly, parameters will not match if \code{scaling} is not 1.
#' 
#' @author Christopher J. Paciorek
#' @export
#' @details
#' See \code{\link{fit_gev}} for more details on fitting the block maxima model for each dataset, including details on blocking and replication. Also see \code{\link{fit_gev}} for information on the \code{bootControl} argument.
#' 
#' Optimization failures:
#'
#' It is not uncommon for maximization of the log-likelihood to fail for extreme value models. Please see the help information for \code{fit_gev}. Also note that if the probability in the denominator of the risk ratio is near one, one may achieve better numerical performance by swapping the two datasets and computing the risk ratio for the probability under dataset 2 relative to the probability under dataset 1.
#' 
#' @return
#' The primary outputs of this function are as follows: the log of the risk ratio and standard error of that log risk ratio (\code{logRiskRatio} and \code{se_logRiskRatio}) as well the risk ratio itself (\code{riskRatio}). The standard error is based on the usual MLE asymptotics using a delta-method-based approximation. If requested via \code{ciType}, confidence intervals will be returned, as discussed in \code{Details}.
#' 
#' @details
#' \code{ciType} can include one or more of the following: \code{'delta'}, \code{'lrt'}, \code{'boot_norm'}, \code{'boot_perc'}, \code{'boot_basic'}, \code{'boot_stud'}, \code{'boot_bca'}. \code{'delta'} uses the delta method to compute an asymptotic interval based on the standard error of the log risk ratio. \code{'lrt'} inverts a likelihood-ratio test. Bootstrap-based options are the normal-based interval using the bootstrap standard error (\code{'boot_norm'}), the percentile bootstrap (\code{'boot_perc'}), the basic bootstrap (\code{'boot_basic'}), the bootstrap-t (\code{'boot_stud'}), and the bootstrap BCA method (\code{'boot_bca'}). See Paciorek et al. for more details. 
#' 
#' See \code{\link{fit_pot}} for information on the \code{bootControl} argument. 
#' @references
#' Paciorek, C.J., D.A. Stone, and M.F. Wehner. 2018. Quantifying uncertainty in the attribution of human influence on severe weather. Weather and Climate Extremes 20:69-80. arXiv preprint <https://arxiv.org/abs/1706.03388>.
#'
#' Jeon S., C.J. Paciorek, and M.F. Wehner. 2016. Quantile-based bias correction and uncertainty quantification of extreme event attribution statements. Weather and Climate Extremes 12: 24-32. <DOI:10.1016/j.wace.2016.02.001>. arXiv preprint: <http://arxiv.org/abs/1602.04139>.
#' 
#' @examples
#' data(Fort, package = 'extRemes')
#' FortMax <- aggregate(Prec ~ year, data = Fort, max)
#' earlyYears <- 1900:1929
#' lateYears <- 1970:1999
#' earlyPeriod <- which(FortMax$year %in% earlyYears)
#' latePeriod <- which(FortMax$year %in% lateYears)
#' # contrast late period with early period, assuming a nonstationary fit
#' # within each time period and finding RR based on midpoint of each period
#' \dontrun{
#' out <- calc_riskRatio_gev(returnValue = 3,
#'                    y1 = FortMax$Prec[earlyPeriod], y2 = FortMax$Prec[latePeriod],
#'                    x1 = data.frame(years = earlyYears), x2 = data.frame(years = lateYears),
#'                    locationFun1 = ~years, locationFun2 = ~years,
#'                    xNew1 = data.frame(years = mean(earlyYears)),
#'                    xNew2 = data.frame(years = mean(lateYears)))
#' }
calc_riskRatio_gev <- function(returnValue, y1, y2, x1 = NULL, x2 = x1,
                               locationFun1 = NULL, locationFun2 = locationFun1,
                               scaleFun1 = NULL, scaleFun2 = scaleFun1,
                               shapeFun1 = NULL, shapeFun2 = shapeFun1,
                               nReplicates1 = 1, nReplicates2 = 1,
                               replicateIndex1 = NULL, replicateIndex2 = NULL,
                               weights1 = NULL, weights2 = NULL, xNew1 = NULL, xNew2 = NULL, 
                               maxes = TRUE, scaling1 = 1, scaling2 = 1,
                               ciLevel = 0.90, ciType, bootSE,
                               bootControl = NULL, lrtControl = NULL,
                               optimArgs = NULL, optimControl = NULL, initial1 = NULL,
                               initial2 = NULL, logScale1 = NULL, logScale2 = NULL,
                               getReturnCalcs = FALSE, getParams = FALSE, getFit = FALSE) {

    if(is.null(returnValue) || is.null(y1) || is.null(y2))  # for Python interface
        stop("calc_riskRatio_binom: argument 'returnValue' or 'y1' or 'y2' is missing, with no default")

    if(missing(bootSE) || is.null(bootSE))
        bootSE <- !missing(ciType) && sum(ciType %in% names(bootTypes))
  
    if(missing(returnValue))
        stop("calc_riskRatio_gev: 'returnValue' must be provided.")
    if(length(returnValue) > 1)
        stop("calc_riskRatio_gev: please provide a single value for 'returnValue'.")
    z_alpha <- stats::qnorm(( 1 - ciLevel ) / 2, lower.tail = FALSE)
    ciLabels <- as.character(c(( 1 - ciLevel ) / 2, 1 - ( 1 - ciLevel ) / 2))

    if(!is.null(xNew1)) {
        xNew1tmp <- try(as.data.frame(xNew1))
        if(is(xNew1tmp, 'try-error')) stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
        m <- nrow(xNew1tmp)
    } else {
        if(is.null(x1))
            m <- 1 else {
                       x1tmp <- try(as.data.frame(x1))
                       if(is(x1tmp, 'try-error')) stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
                       m <- nrow(x1tmp)
                   }
    }
    
    results <- list()
    
    if(missing(ciType)) {
        ciType <- ""
    } else {
        wh <- setdiff(ciType, c('delta', 'lrt', 'koopman', names(bootTypes)))
        if(length(wh))
            stop("calc_riskRatio_gev: ", paste(wh, collapse = ','), " in 'ciType' are not valid types.")
    }

    bootTypesUse <- ciType[ciType %in% names(bootTypes)]
    if(length(bootTypesUse)) {
        bControl <- list(seed = 1, n = 250, by = "block", getSample = FALSE)
        bControl[names(bootControl)] <- bootControl
        bootControl <- bControl
        bootControlTmp <- bootControl
        bootControlTmp$getSample <- TRUE
        bootControlTmp$getSampleSE <- TRUE
    } 
    
    fit1 <- fit_gev(y1, x = x1, locationFun = locationFun1,
                    scaleFun = scaleFun1, shapeFun = shapeFun1, 
                    nReplicates = nReplicates1, replicateIndex = replicateIndex1,
                    weights = weights1, returnValue = returnValue,
                    xNew = xNew1, maxes = maxes,
                    scaling = scaling1, bootSE = bootSE, bootControl = bootControlTmp,
                    optimArgs = optimArgs, initial = initial1, logScale = logScale1,
                    getFit = TRUE, getParams = getParams, .getInputs = TRUE) 
    fit2 <- fit_gev(y2, x = x2, locationFun = locationFun2,
                    scaleFun = scaleFun2, shapeFun = shapeFun2, 
                    nReplicates = nReplicates2, replicateIndex = replicateIndex2,
                    weights = weights2, returnValue = returnValue,
                    xNew = xNew2, maxes = maxes,
                    scaling = scaling2, bootSE = bootSE, bootControl = bootControlTmp,
                    optimArgs = optimArgs, initial = initial2, logScale = logScale2,
                    getFit = TRUE, getParams = getParams, .getInputs = TRUE) 
    
    if(fit1$info$failure || fit2$info$failure) {
        warning("calc_riskRatio_gev: fitting failed for one of two datasets.")
        results$logRiskRatio <- results$se_logRiskRatio <- results$riskRatio <- rep(NA, m)
        for(type in ciType)
            results[[paste0('ci_riskRatio_', type)]] <- drop(matrix(NA, m, 2))
        if(bootSE) 
            results$se_logRiskRatio_boot <- rep(NA, m)
    } else {
        if(length(fit1$logReturnProb) != length(fit2$logReturnProb))
            stop("calc_riskRatio_gev: number of return probabilities calculated for each model fit must be the same; for nonstationary models this is determined by the number of covariate set inputs (provided in 'x' or 'xNew')")
        results$logRiskRatio <- fit1$logReturnProb - fit2$logReturnProb

        results$se_logRiskRatio <- sqrt(fit1$se_logReturnProb^2 + fit2$se_logReturnProb^2)
        results$riskRatio <- exp(results$logRiskRatio)
        if('delta' %in% ciType) {
            results$ci_riskRatio_delta <- exp(cbind(results$logRiskRatio - z_alpha*results$se_logRiskRatio,
                                                    results$logRiskRatio + z_alpha*results$se_logRiskRatio))
            colnames(results$ci_riskRatio_delta) <- ciLabels
            results$ci_riskRatio_delta <- drop(results$ci_riskRatio_delta)
        }
        
        if(length(bootTypesUse)) {
            bootData <- matrix(drop(fit1$logReturnProb_boot - fit2$logReturnProb_boot), ncol = m)
            if(bootSE)
                results$se_logRiskRatio_boot <- apply(bootData, 2, sd, na.rm = TRUE)
            
            fake_data <- cbind(sample(c(0,1), size = 5, replace = TRUE),
                   sample(c(0,1), size = 5, replace = TRUE))
            logRRfun <- function(dat, ind) {
                log(sum(dat[ind,1])/sum(dat[ind,2]))
            }
            bootInput <- boot::boot(fake_data, logRRfun, R = bootControl$n)

            for(type in bootTypesUse) {
                typeName <- paste0('ci_riskRatio_', type)
                results[[typeName]] <- matrix(NA, m, 2)
                colnames(results[[typeName]]) <- ciLabels
            }
            for(i in seq_len(m)) {
                bootInput$t0 = c(results$logRiskRatio[i], results$se_logRiskRatio[i]^2)
                logRRvals <- bootData[ , i]
                logRRvals[logRRvals == Inf] = 1e6 # boot.ci can't handle Inf
                bootInput$t <- cbind(logRRvals, fit1$logReturnProb_boot_se[ , i, 1]^2 + fit2$logReturnProb_boot_se[ , i, 1]^2)
                bootResults <- try(boot::boot.ci(bootInput, conf = ciLevel, type = unlist(bootTypes[bootTypesUse])))
                if(!is(bootResults, 'try-error')) {
                    for(type in bootTypesUse) {
                        typeName <- paste0('ci_riskRatio_', type)
                        if(type == 'boot_norm') {
                            results[[typeName]][i, ] <- exp(bootResults[[bootTypesCols[[type]]]][2:3])
                        } else results[[typeName]][i, ] <- exp(bootResults[[bootTypesCols[[type]]]][4:5])
                        if(i == m) results[[typeName]] <- drop(results[[typeName]])
                    }
                }
            }
            if(bootControl$getSample) results$riskRatio_boot <- drop(exp(bootData))
        }
        
        if('lrt' %in% ciType) {
            lControl <- list(bounds = c(.01, 100))
            lControl[names(lrtControl)] <- lrtControl
            lrtControl <- lControl
            oArgs <- list(method = "Nelder-Mead", lower = -Inf, upper = Inf, control = list())
            oArgs[names(optimArgs)] <- optimArgs
            if(!is.null(optimControl))
                oArgs[['control']] <- optimControl
            if(!maxes) returnValue <- -returnValue
            results$ci_riskRatio_lrt <- calc_riskRatio_lrt(fit1, fit2, returnValue = returnValue, ciLevel = ciLevel, bounds = lrtControl$bounds, type = "GEV", optimArgs = oArgs)
            if(is.null(dim(results$ci_riskRatio_lrt))) {
                names(results$ci_riskRatio_lrt) <- ciLabels
            } else colnames(results$ci_riskRatio_lrt) <- ciLabels
        }
    }
    if(getFit || getParams || getReturnCalcs) {
        fit1$inputs <- fit2$inputs <- NULL
        if(!getFit) {
            fit1$fit <- fit2$fit <- fit1$info <- fit2$info <- NULL
        }
        if(!getReturnCalcs) {
            fit1$logReturnProb <- fit1$se_logReturnProb <- fit1$logReturnPeriod <- fit1$se_logReturnPeriod <- NULL
            fit2$logReturnProb <- fit2$se_logReturnProb <- fit2$logReturnPeriod <- fit2$se_logReturnPeriod <- NULL
            if(bootSE) {
                fit1$se_logReturnProb_boot <- fit1$se_logReturnPeriod <- NULL
                fit2$se_logReturnProb_boot <- fit2$se_logReturnPeriod <- NULL
            }
        }
        results$fit1 <- fit1
        results$fit2 <- fit2
    }
    return(results)
} # end calc_riskRatio_gev

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
        qcov1 <- make.qcov_safe(fit1$fit, covariateMatrix)
    } else qcov1 <- NULL
    if(sum(unlist(p2)) > 3) {
        covariateMatrix <- stats::model.matrix(inputs2$locationFun, inputs2$xUse)
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(inputs2$scaleFun, inputs2$xUse))
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(inputs2$shapeFun, inputs2$xUse))
        qcov2 <- make.qcov_safe(fit2$fit, covariateMatrix)
    } else qcov2 <- NULL

    objfun <- function(logrr0, logLikHat, cutoff, initial, y1, y2, x1, x2, threshold1, threshold2, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1, qcov2, type, optimArgs) {
        rr0 <- exp(logrr0)
        fitConstr <- stats::optim(initial, calc_constrNegLogLik, gr = NULL, rr0, y1, y2, x1, x2, threshold1, threshold2, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1, qcov2, type, method = optimArgs$method, lower = optimArgs$lower, upper = optimArgs$upper, control = optimArgs$control)$value
        output <- logLikHat + fitConstr - cutoff
        attributes(output) <- NULL
        ## uniroot has trouble extending interval when attributes present
        return(output)
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
            lower <- stats::uniroot(objfun, interval = intvl, logLikHat, cutoff, initial, inputs1$y, inputs2$y, x1, x2, inputs1$thresholdObs, inputs2$thresholdObs, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1[i, ], qcov2[i, ], type, optimArgs, extendInt = 'downX')$root
            if(all.equal(lower, intvl[1]) == TRUE || all.equal(lower, intvl[2]) == TRUE)
                warning("optimization likely failed to converge: lower endpoint equal to one endpoint of search interval.")
            ## however note that with use of uniroot, it will extend the interval
        }
        if(rrHat[i] > bounds[2]) {
            upper <- NA
        } else {
            intvl <- log(c(max(rrHat[i], bounds[1]), bounds[2]))
            upper <- stats::uniroot(objfun, interval = intvl, logLikHat, cutoff, initial, inputs1$y, inputs2$y, inputs1$xObs, inputs2$xObs, inputs1$thresholdObs, inputs2$thresholdObs, weights1, weights2, p1, p2, designs1, designs2, blocks1, blocks2, usePhi1, usePhi2, returnValue, qcov1[i, ], qcov2[i, ], type, optimArgs, extendInt = 'upX')$root
            if(all.equal(upper, intvl[1]) == TRUE || all.equal(upper, intvl[2]) == TRUE)
                warning("optimization likely failed to converge: upper endpoint equal to one endpoint of search interval.")
        }
        output[i, ] <- exp(c(lower, upper))
    } # end loop over 1:inputs1$mUse
    return(drop(output))
} # end calc_riskRatio_lrt()

calc_constrNegLogLik <- function(par, rr0, y1, y2, x1, x2, threshold1 = NULL, threshold2 = NULL, weights1, weights2, p1, p2, designs1, designs2, blocks1 = NULL, blocks2 = NULL, usePhi1, usePhi2, returnValue, qcov1, qcov2, type = "PP"){
    infval <- 1e6
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
    if(scale <= 0) return(infval)
    
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
    if(scale <= 0) return(infval)

    if(pA/rr0 > 1) return(infval)
    
    par2[1] <- returnValue - location + (scale/shape) * (1 - (-log(1 - pA/rr0))^-shape)
    if(!is.finite(par2[1]))
        return(infval)

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
    if(result > infval) result <- infval # larger values such as 1e10 from oevd() seem to lead to early stopping of optimization over this function in some cases
    return(result)
} # end calc_constrNegLogLik()
