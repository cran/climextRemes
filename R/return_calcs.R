#' Calculates return value and standard error given a return period of interest 
#'
#' Calculates return value (also known as the return level) given a return period of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return value is the value for which the expected number of blocks until an event that exceeds that value is equal to the return period. For non-stationary models (those that include covariates for the location, scale, and/or shape parameters, return values and standard errors are returned for as many sets of covariates as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnPeriod value for which return value is desired
#' @param covariates matrix of covariate values, each row a set of covariates for which the return value is desired
#'
#' @export
#'
calc_returnValue_fevd <- function(fit, returnPeriod, covariates = NULL) {
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_returnValue_fevd: can only be used for PP and GEV models.")
      
    if(!is.null(covariates)) {
        if(!is.matrix(covariates)) stop("calc_returnValue_fevd: 'covariates' should be a matrix.")
        m <- nrow(covariates)
        qcov <- make.qcov(fit, covariates)
        tmp <- try(return.level.ns.fevd.mle(fit, return.period = returnPeriod, qcov = qcov, do.ci = TRUE))
        if(!methods::is(tmp, 'try-error')) {
            names(tmp) <- NULL
            results <- list(returnValue = tmp[ , 2], se_returnValue = tmp[ , 4])
        } else results <- list(returnValue = rep(NA, m), se_returnValue = rep(NA, m))
    } else {
        # need a specific alpha level to reverse-engineer the standard error
        alpha <- 0.05
        z.alpha <- stats::qnorm(alpha/2, lower.tail = FALSE)

        tmp <- try(return.level(fit, return.period = returnPeriod, alpha = alpha, do.ci = TRUE))
        if(!methods::is(tmp, 'try-error')) {
            names(tmp) <- NULL
            se <- (tmp[3] - tmp[1]) / (2 * z.alpha) # reverse-engineer the s.e. from the CI
            results <- list(returnValue = tmp[2], se_returnValue = se)
        } else results <- list(returnValue = NA, se_returnValue = NA)
    }
    return(results)
}

#' Calculates return value difference for two sets of covariates and standard error of difference given a return period of interest 
#'
#' Calculates difference in return values (also known as return levels) for two sets of covariates given a return period of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return value is the value for which the expected number of blocks until an event that exceeds that value is equal to the return period. Differences and standard errors are returned for as many contrasts between covariate sets as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnPeriod value for which return value difference is desired
#' @param covariates1 matrix of covariate values, each row a set of covariates for which the return value difference relative to the corresponding row of \code{covariates2} is desired
#' @param covariates2 matrix of covariate values, each row a set of covariates 
#'
#' @details This is designed to calculate differences in return values and associated standard errors for different covariate values based on the same model fit. It is not designed for differences based on separate model fits, although it may be possible handle this case by fit two models in a single model specification using dummy variables.
#' @export
#'
calc_returnValueDiff_fevd <- function(fit, returnPeriod, covariates1, covariates2) {
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_returnValueDiff_fevd: can only be used for PP and GEV models.")

    if(!is.matrix(covariates1) || !is.matrix(covariates2)) stop("calc_returnValueDiff_fevd: 'covariates1' and 'covariates2' should be matrices.")
    m <- nrow(covariates1)
    if(nrow(covariates2) != m)
        stop("calc_returnValueDiff_fevd: 'covariates1' and 'covariates2' must have same number of rows.")
    
    results <- list()

    qcov1 <- make.qcov(fit, covariates1)
    qcov2 <- make.qcov(fit, covariates2)
    tmp <- try(return.level.ns.fevd.mle(fit, return.period = returnPeriod, alpha = .05, qcov = qcov1, qcov.base = qcov2, do.ci = TRUE))
    names(tmp) <- NULL
    if(!methods::is(tmp, 'try-error')) 
        results <- list(returnValueDiff = tmp[ , 2], se_returnValueDiff = tmp[ , 4])

    return(results)
}


#' Calculates log return period and standard error given a return value of interest 
#'
#' Calculates log return period given a return value of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return period is the average number of blocks expected to occur before the return value is exceeded and is equal to the inverse of the probability of exceeding the return value in a single block. For non-stationary models (those that include covariates for the location, scale, and/or shape parameters, log return periods and standard errors are returned for as many sets of covariates as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnValue value for which return period is desired
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


#' Calculates log return probability and standard error given a return value of interest 
#'
#' Calculates log return probability given a return value of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return probability is the probability of exceeding the return value in a single block. For non-stationary models (those that include covariates for the location, scale, and/or shape parameters, log probabilities and standard errors are returned for as many sets of covariates as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnValue value for which log return probability is desired
#' @param covariates matrix of covariate values, each row a set of covariates for which the return probability is desired
#'
#' @details Results are calculated (and returned) on log scale as delta-method based standard errors are more accurate for the log probability. Confidence intervals on the probability scale should be calculated by calculating a confidence interval for the log probability and exponentiating the endpoints of the interval. 
#' @export
#'
calc_logReturnProb_fevd <- function(fit, returnValue, covariates = NULL){
    if(!fit$type %in% c('PP', 'GEV'))
        stop("calc_logReturnProb_fevd: can only be used for PP and GEV models.")
    
    p <- fit$results$num.pars
    locationIndices <- 1:p[['location']]
    scaleIndices <- (1 + p[['location']]):(p[['location']] + p[['scale']])
    shapeIndices <- (1 + p[['location']] + p[['scale']]):sum(unlist(p))

    summ <- summary(fit, silent = TRUE)
    mle <- summ$par
    cov_mle <- summ$cov.theta
  
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
        se <- rep(0, m)
        
        logprob <- sapply(seq_along(location), function(i) pevd(returnValue, location[i], scale[i], shape[i], lower.tail = FALSE, type = "GEV", log.p = TRUE))
        fdF <- devd(rep(returnValue, length(location)), location, scale, shape) / exp(logprob)
        for(i in seq_len(m)) {
            std <- (returnValue - location[i]) / scale[i]
            grad <- fdF[i] * c(
                                 qcov[i, locationIndices],
                                 ((returnValue - location[i]) / ifelse(!fit$const.scale, 1, scale[i])) * qcov[i, scaleIndices], # fixed with dsigma/dphi relative to llex
                                 scale[i] * ( -std / shape[i] + (1 + shape[i]*std)*log(1 + shape[i]*std) / shape[i]^2 ) * qcov[i, shapeIndices])
            se[i] <- sqrt(t(grad) %*% (cov_mle %*% grad))
        }
    } else {
        logprob <- pevd(returnValue, location, scale, shape, lower.tail = FALSE, type = "GEV", log.p = TRUE)
        fdF <- devd(returnValue, location, scale, shape, type = "GEV") / exp(logprob)
        std <- (returnValue - location) / scale
        grad <- fdF * c(1, std, scale * (-std/shape + (1 + shape*std)*log(1 + shape*std) / shape^2 ))
        se <- sqrt(t(grad) %*% (cov_mle %*% grad))[1,1]
    }
    
    return(list(logReturnProb = logprob, se_logReturnProb = se))
}

#' Calculates log return probability difference for two sets of covariates and standard error of difference given a return value of interest 
#'
#' Calculates difference in log return probabilities for two sets of covariates given a return value of interest, using model fit from \code{extRemes::fevd}. Standard error is obtained via the delta method. The return probability is the probability of exceeding the return value in a single block. Differences and standard errors are returned for as many contrasts between covariate sets as provided. 
#' @param fit fitted object from \pkg{extRemes} \code{fevd}
#' @param returnValue value for which the log return probability difference is desired
#' @param covariates1 matrix of covariate values, each row a set of covariates for which the log return probability difference relative to the corresponding row of \code{covariates2} is desired
#' @param covariates2 matrix of covariate values, each row a set of covariates 
#'
#' @details Results are calculated (and returned) on log scale as delta-method based standard errors are more accurate for the log probability. Confidence intervals for the ratio of return probabilities should be calculated by calculating a confidence interval for the log probability difference and exponentiating the endpoints of the interval.
#'
#' This is designed to calculate differences in log return probabilities and associated standard errors for different covariate values based on the same model fit. It is not designed for differences based on separate model fits, although it may be possible handle this case by fit two models in a single model specification using dummy variables.
#' @export
#'
calc_logReturnProbDiff_fevd <- function(fit, returnValue, covariates1, covariates2){
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
    
    m <- nrow(qcov1)
    if(nrow(covariates2) != m)
        stop("calc_logReturnProbDiff_fevd: 'covariates1' and 'covariates2' must have same number of rows.")
    se <- rep(0, m)
    logprob1 <- sapply(seq_along(location1), function(i) pevd(returnValue, location1[i], scale1[i], shape1[i], lower.tail = FALSE, type = "GEV", log.p = TRUE))
    fdF1 <- devd(rep(returnValue, length(location1)), location1, scale1, shape1) / exp(logprob1)
    logprob2 <- sapply(seq_along(location2), function(i) pevd(returnValue, location2[i], scale2[i], shape2[i], lower.tail = FALSE, type = "GEV", log.p = TRUE))
    fdF2 <- devd(rep(returnValue, length(location2)), location2, scale2, shape2) / exp(logprob2)
    for(i in seq_len(m)) {
        std1 <- (returnValue - location1[i]) / scale1[i]
        std2 <- (returnValue - location2[i]) / scale2[i]
        
        grad1 <- fdF1[i] * c(
                               qcov1[i, locationIndices],
                               ((returnValue - location1[i]) / ifelse(!fit$const.scale, 1, scale1[i])) * qcov1[i, scaleIndices], 
                               scale1[i] * ( -std1 / shape1[i] + (1 + shape1[i]*std1)*log(1 + shape1[i]*std1) / shape1[i]^2 ) * qcov1[i, shapeIndices])
       grad2 <- fdF2[i] * c(
                              qcov2[i, locationIndices],
                              ((returnValue - location2[i]) / ifelse(!fit$const.scale, 1, scale2[i])) * qcov2[i, scaleIndices], 
                              scale2[i] * ( -std2 / shape2[i] + (1 + shape2[i]*std2)*log(1 + shape2[i]*std2) / shape2[i]^2 ) * qcov2[i, shapeIndices])

        grad <- grad1 - grad2
        se[i] <- sqrt(t(grad) %*% (cov_mle %*% grad))[1,1]
    }
    
    return(list(logReturnProbDiff = logprob1 - logprob2, se_logReturnProbDiff = se))


}
            
