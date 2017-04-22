#' Calculates return value and standard error given return period(s) of interest 
#'
#' Calculates return value (also known as the return level) given return period(s) of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return value is the value for which the expected number of blocks until an event that exceeds that value is equal to the return period. For non-stationary models (those that include covariates for the location, scale, and/or shape parameters, return values and standard errors are returned for as many sets of covariates as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnPeriod value(s) for which return value is desired
#' @param covariates matrix of covariate values, each row a set of covariates for which the return value is desired
#'
#' @export
#'
calc_returnValue_fevd <- function(fit, returnPeriod, covariates = NULL) {
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_returnValue_fevd: can only be used for PP and GEV models.")

    nPeriod <- length(returnPeriod)
    
    if(!is.null(covariates)) {
        if(!is.matrix(covariates)) stop("calc_returnValue_fevd: 'covariates' should be a matrix.")
        m <- nrow(covariates)
        qcov <- make.qcov(fit, covariates)
        if(nPeriod == 1) {
            tmp <- try(return.level.ns.fevd.mle(fit, return.period = returnPeriod, qcov = qcov, do.ci = TRUE))
            if(!methods::is(tmp, 'try-error')) {
                names(tmp) <- NULL
                results <- list(returnValue = tmp[ , 2], se_returnValue = tmp[ , 4])
            } else results <- list(returnValue = rep(NA, m), se_returnValue = rep(NA, m))
        } else {
            rv <- se <- matrix(as.numeric(NA), m, nPeriod)
            colnames(rv) <- colnames(se) <- returnPeriod
            for(j in seq_len(nPeriod)) {
                tmp <- try(return.level.ns.fevd.mle(fit, return.period = returnPeriod[j], qcov = qcov, do.ci = TRUE))
                if(!methods::is(tmp, 'try-error')) {
                    rv[ , j] <- tmp[ , 2]
                    se[ , j] <- tmp[ , 4]
                }
            }
            results <- list(returnValue = rv, se_returnValue = se)
        }
    } else {
        # need a specific alpha level to reverse-engineer the standard error
        alpha <- 0.05
        z.alpha <- stats::qnorm(alpha/2, lower.tail = FALSE)

        tmp <- try(return.level(fit, return.period = returnPeriod, alpha = alpha, do.ci = TRUE))
        if(!methods::is(tmp, 'try-error')) {
            if(nPeriod == 1) {
                names(tmp) <- NULL
                se <- (tmp[3] - tmp[1]) / (2 * z.alpha) # reverse-engineer the s.e. from the CI
                results <- list(returnValue = tmp[2], se_returnValue = se)
            } else {
                se <- (tmp[ , 3] - tmp[ , 1]) / (2 * z.alpha)
                rv <- tmp[ , 2]
                names(rv) <- names(se) <- returnPeriod
                results <- list(returnValue = rv, se_returnValue = se)
            }
        } else {
            rv <- se <- rep(NA, nPeriod)
            if(nPeriod > 1) names(rv) <- names(se) <- returnPeriod
            results <- list(returnValue = rv, se_returnValue = se)
        }
    }
    return(results)
}

#' Calculates return value difference for two sets of covariates and standard error of difference given return period(s) of interest 
#'
#' Calculates difference in return values (also known as return levels) for two sets of covariates given return period(s) of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return value is the value for which the expected number of blocks until an event that exceeds that value is equal to the return period. Differences and standard errors are returned for as many contrasts between covariate sets as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnPeriod value(s) for which return value difference is desired
#' @param covariates1 matrix of covariate values, each row a set of covariates for which the return value difference relative to the corresponding row of \code{covariates2} is desired
#' @param covariates2 matrix of covariate values, each row a set of covariates 
#' @param getSE logical indicating whether standard error is desired, in addition to the point estimate
#'
#' @details This is designed to calculate differences in return values and associated standard errors for different covariate values based on the same model fit. It is not designed for differences based on separate model fits, although it may be possible handle this case by fit two models in a single model specification using dummy variables.
#' @export
#'
calc_returnValueDiff_fevd <- function(fit, returnPeriod, covariates1, covariates2, getSE = TRUE) {
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_returnValueDiff_fevd: can only be used for PP and GEV models.")

    if(!is.matrix(covariates1) || !is.matrix(covariates2)) stop("calc_returnValueDiff_fevd: 'covariates1' and 'covariates2' should be matrices.")
    m <- nrow(covariates1)
    if(nrow(covariates2) != m)
        stop("calc_returnValueDiff_fevd: 'covariates1' and 'covariates2' must have same number of rows.")
    
    results <- list()

    nPeriod <- length(returnPeriod)
    
    qcov1 <- make.qcov(fit, covariates1)
    qcov2 <- make.qcov(fit, covariates2)

    if(nPeriod == 1) {
        tmp <- try(return.level.ns.fevd.mle(fit, return.period = returnPeriod, alpha = .05, qcov = qcov1, qcov.base = qcov2, do.ci = getSE))
        names(tmp) <- NULL
        if(!methods::is(tmp, 'try-error')) {
            if(getSE) {
                results <- list(returnValueDiff = tmp[ , 2], se_returnValueDiff = tmp[ , 4])
            } else results <- list(returnValueDiff = tmp, se_returnValueDiff = rep(as.numeric(NA), m))
        } else results <- list(returnValueDiff = rep(as.numeric(NA), m), se_returnValueDiff = rep(as.numeric(NA), m))
    } else {
        rvd <- se <- matrix(as.numeric(NA), m, nPeriod)
        colnames(rvd) <- colnames(se) <- returnPeriod
        for(j in seq_len(nPeriod)) {
            tmp <- try(return.level.ns.fevd.mle(fit, return.period = returnPeriod[j], alpha = .05, qcov = qcov1, qcov.base = qcov2, do.ci = getSE))
            if(!methods::is(tmp, 'try-error')) {
                if(getSE) {
                    rvd[ , j] <- tmp[ , 2]
                    se[ , j] <- tmp[ , 4]
                } else rvd[ , j] <- tmp
            }
        }
        results <- list(returnValueDiff = rvd, se_returnValueDiff = se)
    }
    return(results)
}


#' Calculates log return period and standard error given return value(s) of interest 
#'
#' Calculates log return period given return value(s) of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return period is the average number of blocks expected to occur before the return value is exceeded and is equal to the inverse of the probability of exceeding the return value in a single block. For non-stationary models (those that include covariates for the location, scale, and/or shape parameters, log return periods and standard errors are returned for as many sets of covariates as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnValue value(s) for which return period is desired
#' @param covariates matrix of covariate values, each row a set of covariates for which the log return period is desired
#'
#' @details Results are calculated (and returned) on log scale as delta-method based standard errors are more accurate for the log period. Confidence intervals on the return period scale should be calculated by calculating a confidence interval for the log return period and exponentiating the endpoints of the interval. 
#' @export
#'
calc_logReturnPeriod_fevd <- function(fit, returnValue, covariates = NULL) {
    # duplicates already-done calculation if returnProb already calculated...
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_logReturnPeriod_fevd: can only be used for PP and GEV models.")
    results <- calc_logReturnProb_fevd(fit, returnValue, covariates)
    results$logPeriod <- -results$logProb
    results$se_logPeriod <- results$se_logProb
    results$logProb <- results$se_logProb <- NULL
    return(results)
}


#' Calculates log return probability and standard error given return value(s) of interest 
#'
#' Calculates log return probability given return value(s) of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return probability is the probability of exceeding the return value in a single block. For non-stationary models (those that include covariates for the location, scale, and/or shape parameters, log probabilities and standard errors are returned for as many sets of covariates as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnValue value(s) for which log return probability is desired
#' @param covariates matrix of covariate values, each row a set of covariates for which the return probability is desired
#' @param getSE logical indicating whether standard error is desired, in addition to the point estimate
#' @param scaling if \code{returnValue} is scaled for numerics, this allows names of output to be on original scale
#'
#' @details Results are calculated (and returned) on log scale as delta-method based standard errors are more accurate for the log probability. Confidence intervals on the probability scale should be calculated by calculating a confidence interval for the log probability and exponentiating the endpoints of the interval. 
#' @export
#'
calc_logReturnProb_fevd <- function(fit, returnValue, covariates = NULL, getSE = TRUE, scaling = 1){
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_logReturnProb_fevd: can only be used for PP and GEV models.")
    
    p <- fit$results$num.pars
    locationIndices <- 1:p[['location']]
    scaleIndices <- (1 + p[['location']]):(p[['location']] + p[['scale']])
    shapeIndices <- (1 + p[['location']] + p[['scale']]):sum(unlist(p))

    summ <- summary(fit, silent = TRUE)
    mle <- summ$par
    cov_mle <- summ$cov.theta

    nValue <- length(returnValue)
    
    if(!is.null(covariates)) {
        qcov <- make.qcov(fit, covariates)
        location <- qcov[, locationIndices, drop = FALSE] %*% mle[locationIndices]
        scale <- qcov[, scaleIndices, drop = FALSE] %*% mle[scaleIndices]
        shape <- qcov[, shapeIndices, drop = FALSE] %*% mle[shapeIndices]
    } else {
        location <- mle[locationIndices]
        scale <- mle[scaleIndices]
        shape <- mle[shapeIndices]
    }
    if(!fit$const.scale) scale <- exp(scale)

    if(!is.null(covariates)) {  # deal with how pevd/devd do/don't vectorize wrt parameters
        if(!is.matrix(covariates)) stop("calc_logReturnProb_fevd: 'covariates' should be a matrix.")        
        m <- nrow(qcov)

        logprob <- se <- matrix(as.numeric(NA), m, nValue)
        if(nValue > 1) colnames(logprob) <- colnames(se) <- returnValue/scaling

        for(j in seq_len(nValue)) {
            logprob[ , j] <- sapply(seq_along(location), function(i) pevd(returnValue[j], location[i], scale[i], shape[i], lower.tail = FALSE, type = "GEV", log.p = TRUE))
            if(getSE) {
                fdF <- devd(rep(returnValue[j], length(location)), location, scale, shape) / exp(logprob[ , j])
                for(i in seq_len(m)) {
                    std <- (returnValue[j] - location[i]) / scale[i]
                    grad <- fdF[i] * c(
                                         qcov[i, locationIndices],
                                         ((returnValue[j] - location[i]) / ifelse(!fit$const.scale, 1, scale[i])) * qcov[i, scaleIndices], # fixed with dsigma/dphi relative to llex
                                         scale[i] * ( -std / shape[i] + (1 + shape[i]*std)*log(1 + shape[i]*std) / shape[i]^2 ) * qcov[i, shapeIndices])
                    se[i, j] <- sqrt(t(grad) %*% (cov_mle %*% grad))
                }
            }
        }
        if(nValue == 1) {
            logprob <- logprob[ , 1]
            se <- se[ , 1]
        }
    } else {
        logprob <- se <- rep(as.numeric(NA), nValue)
        for(j in seq_len(nValue)) {
            logprob[j] <- pevd(returnValue[j], location, scale, shape, lower.tail = FALSE, type = "GEV", log.p = TRUE)
            if(getSE) {
                fdF <- devd(returnValue[j], location, scale, shape, type = "GEV") / exp(logprob[j])
                std <- (returnValue[j] - location) / scale
                grad <- fdF * c(1, std, scale * (-std/shape + (1 + shape*std)*log(1 + shape*std) / shape^2 ))
                se[j] <- sqrt(t(grad) %*% (cov_mle %*% grad))[1,1]
            } else se <- rep(NA, length(logprob))
        }
        if(nValue > 1) names(logprob) <- names(se) <- returnValue/scaling
    }
    
    return(list(logReturnProb = logprob, se_logReturnProb = se))
}

#' Calculates log return probability difference for two sets of covariates and standard error of difference given return value(s) of interest 
#'
#' Calculates difference in log return probabilities for two sets of covariates given return value(s) of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return probability is the probability of exceeding the return value in a single block. Differences and standard errors are returned for as many contrasts between covariate sets as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnValue value(s) for which the log return probability difference is desired
#' @param covariates1 matrix of covariate values, each row a set of covariates for which the log return probability difference relative to the corresponding row of \code{covariates2} is desired
#' @param covariates2 matrix of covariate values, each row a set of covariates 
#' @param getSE logical indicating whether standard error is desired, in addition to the point estimate
#' @param scaling if \code{returnValue} is scaled for numerics, this allows names of output to be on original scale
#' 
#' @details Results are calculated (and returned) on log scale as delta-method based standard errors are more accurate for the log probability. Confidence intervals for the ratio of return probabilities should be calculated by calculating a confidence interval for the log probability difference and exponentiating the endpoints of the interval.
#'
#' This is designed to calculate differences in log return probabilities and associated standard errors for different covariate values based on the same model fit. It is not designed for differences based on separate model fits, although it may be possible handle this case by fit two models in a single model specification using dummy variables.
#' @export
#'
calc_logReturnProbDiff_fevd <- function(fit, returnValue, covariates1, covariates2, getSE = TRUE, scaling = 1){
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_logReturnProbDiff: can only be used for PP and GEV models.")
    
    if(!is.matrix(covariates1) || !is.matrix(covariates2)) stop("calc_logReturnProbDiff_fevd: 'covariates1' and 'covariates2' should be matrices.")

    p <- fit$results$num.pars
    locationIndices <- 1:p[['location']]
    scaleIndices <- (1 + p[['location']]):(p[['location']] + p[['scale']])
    shapeIndices <- (1 + p[['location']] + p[['scale']]):sum(unlist(p))

    summ <- summary(fit, silent = TRUE)
    mle <- summ$par
    cov_mle <- summ$cov.theta
  
    qcov1 <- make.qcov(fit, covariates1)
    location1 <- qcov1[, locationIndices, drop = FALSE] %*% mle[locationIndices]
    scale1 <- qcov1[, scaleIndices, drop = FALSE] %*% mle[scaleIndices]
    shape1 <- qcov1[, shapeIndices, drop = FALSE] %*% mle[shapeIndices]

    qcov2 <- make.qcov(fit, covariates2)
    location2 <- qcov2[, locationIndices, drop = FALSE] %*% mle[locationIndices]
    scale2 <- qcov2[, scaleIndices, drop = FALSE] %*% mle[scaleIndices]
    shape2 <- qcov2[, shapeIndices, drop = FALSE] %*% mle[shapeIndices]
    if(!fit$const.scale) {
        scale1 <- exp(scale1)
        scale2 <- exp(scale2)
    }

    nValue <- length(returnValue)
    
    m <- nrow(qcov1)
    if(nrow(covariates2) != m)
        stop("calc_logReturnProbDiff_fevd: 'covariates1' and 'covariates2' must have same number of rows.")
    logprob1 <- logprob2 <- se <- matrix(as.numeric(NA), m, nValue)
    if(nValue > 1) colnames(logprob1) <- colnames(logprob2) <- colnames(se) <- returnValue/scaling

    for(j in seq_len(nValue)) {
        logprob1[ , j] <- sapply(seq_along(location1), function(i) pevd(returnValue[j], location1[i], scale1[i], shape1[i], lower.tail = FALSE, type = "GEV", log.p = TRUE))
        logprob2[ , j] <- sapply(seq_along(location2), function(i) pevd(returnValue[j], location2[i], scale2[i], shape2[i], lower.tail = FALSE, type = "GEV", log.p = TRUE))
        if(getSE) {
            fdF1 <- devd(rep(returnValue[j], length(location1)), location1, scale1, shape1) / exp(logprob1[ , j])
            fdF2 <- devd(rep(returnValue[j], length(location2)), location2, scale2, shape2) / exp(logprob2[ , j])
            for(i in seq_len(m)) {
                std1 <- (returnValue[j] - location1[i]) / scale1[i]
                std2 <- (returnValue[j] - location2[i]) / scale2[i]
                
                grad1 <- fdF1[i] * c(
                                       qcov1[i, locationIndices],
                                       ((returnValue[j] - location1[i]) / ifelse(!fit$const.scale, 1, scale1[i])) * qcov1[i, scaleIndices], 
                                       scale1[i] * ( -std1 / shape1[i] + (1 + shape1[i]*std1)*log(1 + shape1[i]*std1) / shape1[i]^2 ) * qcov1[i, shapeIndices])
                grad2 <- fdF2[i] * c(
                                       qcov2[i, locationIndices],
                                       ((returnValue[j] - location2[i]) / ifelse(!fit$const.scale, 1, scale2[i])) * qcov2[i, scaleIndices], 
                                       scale2[i] * ( -std2 / shape2[i] + (1 + shape2[i]*std2)*log(1 + shape2[i]*std2) / shape2[i]^2 ) * qcov2[i, shapeIndices])
                
                grad <- grad1 - grad2
                se[i, j] <- sqrt(t(grad) %*% (cov_mle %*% grad))[1,1]
            }
        }
    }
    if(nValue == 1) {
        logprob1 <- logprob1[ , 1]
        logprob2 <- logprob2[ , 1]
        se <- se[ , 1]
    }
    return(list(logReturnProbDiff = logprob1 - logprob2, se_logReturnProbDiff = se))
}
            
