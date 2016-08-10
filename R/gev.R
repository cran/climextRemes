#' Fit a generalized extreme value model to block maxima or minima
#'
#' Fit a generalized extreme value model, designed specifically for climate data. It includes options for variable weights (useful for local likelihood), as well as for bootstrapping to estimate uncertainties. Results can be returned in terms of parameter values, return values, return periods, return probabilities, and differences in either return values or log return probabilities (i.e., log risk ratios). 
#'
#' @param y a numeric vector of observed maxima or minima values. See \code{Details} for how the values of \code{y} should be ordered if there are multiple replicates and the values of \code{x} are identical for all replicates.
#' @param x a data frame, or object that can be converted to a data frame with columns corresponding to covariate/predictor/feature variables and each row containing the values of the variable for the corresponding observed maximum/minimum. The number of rows should either equal the length of \code{y} or (if there is more than one replicate) it can optionally equal the number of observations in a single replicate, in which case the values will be assumed to be the same for all replicates. 
#' @param locationFun formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the location parameter using columns from \code{x}.  \code{x} must be supplied if this is anything other than NULL or ~1.
#' @param scaleFun formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the log of the scale parameter using columns from \code{x}.  \code{x} must be supplied if this is anything other than NULL or ~1.
#' @param shapeFun formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the shape parameter using columns from \code{x}.  \code{x} must be supplied if this is anything other than NULL or ~1.
#' @param nReplicates numeric value indicating the number of replicates.
#' @param replicateIndex numeric vector providing the index of the replicate corresponding to each element of \code{y}. Used (and therefore required) only when using bootstrapping with the resampling by replicates based on the \code{by} element of \code{bootControl}.
#' @param weights a vector providing the weights for each observation. When there is only one replicate or the weights do not vary by replicate, a vector of length equal to the number of observations. When weights vary by replicate, this should be of equal length to \code{y}. Likelihood contribution of each observation is multiplied by the corresponding weight. 
#' @param returnPeriod numeric value giving the number of blocks for which return values should be calculated. For example a returnPeriod of 20 corresponds to the value of an event that occurs with probability 1/20 in any block and therefore occurs on average every 20 blocks. Often blocks will correspond to years.
#' @param returnValue numeric value giving the value for which return probabilities/periods should be calculated, where the period would be the average number of blocks until the value is exceeded and the probability the probability of exceeding the value in any single block.
#' @param getParams logical indicating whether to return the fitted parameter values and their standard errors; WARNING: parameter values for models with covariates for the scale parameter must interpreted based on log transformation of the scale parameter.
#' @param getFit logical indicating whether to return the full fitted model (potentially useful for model evaluation and for understanding optimization problems); note that estimated parameters in the fit object for nonstationary models will not generally match the MLE provided when \code{getParams} is \code{TRUE} because covariates are normalized before fitting and the fit object is based on the normalized covariates. Similarly, parameters will not match if \code{scaling} is not 1. 
#' @param xNew object of the same form as \code{x}, providing covariate/predictor/feature values for which return values/periods/probabilities are desired.
#' @param xContrast object of the same form and dimensions as \code{xNew}, providing covariate/predictor/feature values for which to calculate the differences of the return values and/or log return probabilities relative to the values in \code{x}. This provides a way to estimate the difference in return value or log return probabilities (i.e., log risk ratios).
#' @param maxes logical indicating whether analysis is for block maxima (TRUE) or block minima (FALSE); in the latter case, the function works with the negative of the values, changing the sign of the resulting location parameters
#' @param scaling positive-valued scalar used to scale the data values for more robust optimization performance. When multiplied by the values, it should produce values with magnitude around 1.
#' @param bootSE logical indicating whether to use the bootstrap to estimate standard errors.
#' @param bootControl a list of control parameters for the bootstrapping. See \code{Details}.
#' @param optimArgs a list with named components matching exactly any arguments that the user wishes to pass to \code{optim}. See \code{help(optim)} for details. Of particular note, \code{'method'} can be used to choose the optimization method used for maximizing the log-likelihood to fit the model and \code{'control=list(maxit=VALUE)'} for a user-chosen VALUE can be used to increase the number of iterations if the optimization is converging slowly.
#' @param missingFlag optional value to be interpreted as missing values (instead of \code{NA}), intended for use in other languages (e.g., Python) calling this function
#' @param initial a list with components named \code{'location'}, \code{'scale'}, and \code{'shape'} providing initial parameter values, intended for use in speeding up or enabling optimization when the default initial values are resulting in failure of the optimization; note that use of \code{scaling} and \code{.normalizeX = TRUE} cause numerical changes in the some of the parameters.
#' @param .normalizeX logical indicating whether to normalize 'x' values for better numerical performance; default is TRUE.
#' @param .getInputs logical indicating whether to return intermediate objects used in fitting. Defaults to \code{FALSE} and intended for internal use only
#' @author Christopher J. Paciorek
#' @export
#' @details
#' This function allows one to fit stationary or nonstationary block maxima/minima models using the generalized extreme value distribution. The function can return parameter estimates, return value/level for a given return period (number of blocks), and return probabilities/periods for a given return value/level. The function provides standard errors based on the usual MLE asymptotics, with delta-method-based standard errors for functionals of the parameters, but also standard errors based on the nonparametric bootstrap, either resampling by block or by replicate or both.
#' 
#' @section Replicates:
#' 
#' Replicates are repeated datasets, each with the same structure, including the same number of block maxima/minima. The additional observations in multiple replicates could simply be treated as additional blocks without replication (see next paragraph), but when the covariate values and weights are the same across replicates, it is simpler to make use of \code{nReplicates} and \code{replicateIndex}.
#' 
#' When using multiple replicates (e.g., multiple members of a climate model initial condition ensemble), the standard input format is to append observations for additional replicates to the \code{y} argument and indicate the replicate ID for each exceedance in \code{replicateIndex}. The values for each replicate should be grouped together and in the same order within replicate so that \code{x} can be correctly matched to the \code{y} values when \code{x} is only supplied for the first replicate. In other words, \code{y} should first contain all the values for the first replicate, then all the values for the second replicate in the same block order as for the first replicate, and so forth. Note that if \code{y} is provided as a matrix with the number of rows equal to the number of observations in each replicate and the columns corresponding to replicates, this ordering will occur naturally.
#' 
#' However, if one has different covariate values for different replicates, then one needs to treat the additional replicates as providing additional blocks, with only a single replicate (and \code{nReplicates} set to 1). The covariate values can then be included as additional rows in \code{x}. Similarly, if there is a varying number of replicates by block, then all block-replicate pairs should be treated as individual blocks with a corresponding row in \code{x} (and \code{nReplicates} set to 1).
#' 
#' @section \code{bootControl} arguments:
#' 
#' The \code{bootControl} argument is a list that can supply any of the following components:
#' \itemize{
#' \item seed. Value of the random number seed as a single value or in the form of \code{.Random.seed} to set before doing resampling. Defaults to \code{0}.
#' \item n. Number of bootstrap samples. Defaults to \code{250}.
#' \item by. Character string, one of \code{'block'}, \code{'replicate'}, or \code{'joint'}, indicating the basis for the resampling. If \code{'block'}, resampled datasets will consist of blocks drawn at random from the original set of blocks; if there are replicates, each replicate will occur once for every resampled block. If \code{'replicate'}, resampled datasets will consist of replicates drawn at random from the original set of replicates; all blocks from a replicate will occur in each resampled replicate. Note that this preserves any dependence across blocks rather than assuming independence between blocks. If \code{'joint'} resampled datasets will consist of block-replicate pairs drawn at random from the original set of block-replicate pairs. Defaults to \code{'block'}. 
#' \item getSample. Logical indicating whether the user wants the full bootstrap sample of parameter estimates and/or return value/period/probability information returned for use in subsequent calculations; if FALSE (the default), only the bootstrap-based estimated standard errors are returned.
#' }
#' 
#' @return
#' The primary outputs of this function are as follows, depending on what is requested via \code{returnPeriod}, \code{returnValue}, \code{getParams} and \code{xContrast}:
#'
#' when \code{returnPeriod} is given: for the period given in \code{returnPeriod} the return value(s) (\code{returnValue}) and its corresponding asymptotic standard error (\code{se_returnValue}) and, when \code{bootSE=TRUE}, also the bootstrapped standard error (\code{se_returnValue_boot}). For nonstationary models, these correspond to the covariate values given in \code{x}.
#'
#' when \code{returnValue} is given: for the value given in \code{returnValue}, the log exceedance probability (\code{logReturnProb}) and the corresponding asymptotic standard error (\code{se_logReturnProb}) and, when \code{bootSE=TRUE}, also the bootstrapped standard error (\code{se_logReturnProb_boot}). This exceedance probability is the probability of exceedance for a single block. Also returned are the log return period (\code{logReturnPeriod}) and its corresponding asymptotic standard error (\code{se_logReturnPeriod}) and, when \code{bootSE=TRUE}, also the bootstrapped standard error (\code{se_logReturnPeriod_boot}). For nonstationary models, these correspond to the covariate values given in \code{x}. Note that results are on the log scale as probabilities and return times are likely to be closer to normally distributed on the log scale and therefore standard errors are more naturally given on this scale. Confidence intervals for return probabilities/periods can be obtained by exponentiating the interval obtained from plus/minus twice the standard error of the log probabilities/periods. 
#'
#' when \code{getParams=TRUE}: the MLE for the model parameters (\code{mle}) and corresponding asymptotic standard error (\code{se_mle}) and, when \code{bootSE=TRUE}, also the bootstrapped standard error (\code{se_mle_boot}).
#'
#' when \code{xContrast} is specified for nonstationary models: the difference in return values (\code{returnValueDiff}) and its corresponding asymptotic standard error (\code{se_returnValueDiff}) and,  when \code{bootSE=TRUE}, bootstrapped standard error (\code{se_returnValueDiff_boot}). These differences correspond to the differences when contrasting each row in \code{x} with the corresponding row in \code{xContrast}. Also returned are the difference in log return probabilities (i.e., the log risk ratio) (\code{logReturnProbDiff}) and its corresponding asymptotic standard error (\code{se_logReturnProbDiff}) and,  when \code{bootSE=TRUE}, bootstrapped standard error (\code{se_logReturnProbDiff_boot}).
#'
#' @references
#' Coles, S. 2001. An Introduction to Statistical Modeling of Extreme Values. Springer.
#' @examples
#' require(extRemes)
#' data(Fort)
#' FortMax <- aggregate(Prec ~ year, data = Fort, max)
#'
#' # stationary fit
#' out <- fit_gev(FortMax$Prec, returnPeriod = 20, returnValue = 3.5,
#'         getParams = TRUE, bootSE = FALSE)
#'
#' # nonstationary fit with location linear in year
#' out <- fit_gev(FortMax$Prec, x = data.frame(years = FortMax$year), 
#'         locationFun = ~years, returnPeriod = 20, returnValue = 3.5,
#'         getParams = TRUE, xNew = data.frame(years = range(FortMax$year)), bootSE = FALSE)
#'
#' @export
fit_gev <- function(y, x = NULL, locationFun = NULL, scaleFun = NULL,
                    shapeFun = NULL, nReplicates = 1, replicateIndex = NULL,
                    weights = NULL, returnPeriod = NULL, returnValue = NULL,
                    getParams = FALSE, getFit = FALSE,
                    xNew = NULL, xContrast = NULL, maxes = TRUE, scaling = 1,
                    bootSE = FALSE, bootControl = list(),
                    optimArgs = list(method = c("Nelder-Mead", "BFGS")), missingFlag = NULL,
                    initial = NULL, .normalizeX = TRUE, .getInputs = FALSE) {

    if(bootSE) {
        bControl <- list(seed = 0, n = 250, by = "block", getSample = FALSE)
        bControl[names(bootControl)] <- bootControl
        bootControl <- bControl
    }
    if(is.null(optimArgs$method) || length(optimArgs$method) > 1)
        optimArgs$method <- c("Nelder-Mead")

    numNumericFun <- sum(c(is.null(locationFun) || is.numeric(locationFun),
             is.null(scaleFun) || is.numeric(scaleFun),
             is.null(shapeFun) || is.numeric(shapeFun)))

    
    if(!is.null(x)) {
        x <- try(as.data.frame(x))
        if(class(x) == 'try-error') stop("fit_gev: 'x' should be a data frame or be able to be converted to a data frame.")
        m <- nrow(x)
        if(any(is.na(x)))
            stop("fit_gev: 'x' should not contain any NAs.")
    } else m <- 0
    if(!is.null(xNew)) {
        xNew <- try(as.data.frame(xNew))
        if(class(xNew) == 'try-error') stop("fit_gev: 'xNew' should be a data frame or be able to be converted to a data frame.")
        mNew <- nrow(xNew)
        if(ncol(x) != ncol(xNew) ||  (numNumericFun != 3 && !identical(sort(names(x)), sort(names(xNew)))))
            stop("fit_gev: columns in 'x' and 'xNew' should be the same.")
    } else mNew <- 0
    if(!is.null(xContrast)) {
        xContrast <- try(as.data.frame(xContrast))
        if(class(xContrast) == 'try-error') stop("fit_gev: 'xContrast' should be a data frame or be able to be converted to a data frame.")
        mContrast <- nrow(xContrast)
        if(is.null(xNew) && m != mContrast) stop("fit_gev: number of sets of covariate values in 'x' and 'xContrast' should be the same.")
        if(!is.null(xNew) && mNew != mContrast) stop("fit_gev: number of sets of covariate values in 'xNew' and 'xContrast' should be the same.")
        if(ncol(xContrast) != ncol(x) ||  (numNumericFun != 3 && !identical(sort(names(x)), sort(names(xContrast)))))
            stop("fit_gev: columns in 'x' and 'xContrast' should be the same.")
    } else mContrast <- 0

    if(nReplicates > 1 && is.null(replicateIndex) && bootSE && bootControl$by == 'replicate')
        stop("fit_gev: 'replicateIndex' must be provided if 'nReplicates' is greater than 1 when bootstrapping by replicate.")


    if(!maxes){  # modeling minima is equivalent to modeling the negative of maxima. Location parameter values will be the negative of those on the original scale, but are corrected before returning parameter values to the user.
        y <- -y
        if(!is.null(returnValue))
            returnValue <- -returnValue
    }
    
    if(!is.numeric(y)) stop("fit_gev: 'y' should be a numeric vector.")

    # ensure parameter functions are in the form of formulae
    locationFun <- parseParamInput(locationFun, names(x))
    scaleFun <- parseParamInput(scaleFun, names(x))
    shapeFun <- parseParamInput(shapeFun, names(x))

    if(numNumericFun == 3) {
        # if using numeric indices of columns, align the names
        if(!is.null(xNew))
            names(xNew) <- names(x)
        if(!is.null(xContrast))
            names(xContrast) <- names(x)
    }

    if(!is.null(x) && locationFun == ~1 && scaleFun == ~1 && shapeFun == ~1) 
        warning("fit_gev: 'x' provided but fitting stationary model based on 'locationFun', 'scaleFun' and 'shapeFun' specification.")

    # scale data for better numerical performance (provided user chose scaling appropriately)
    if(scaling != 1) {
        y <- y * scaling
        if(!is.null(returnValue))
            returnValue <- returnValue * scaling
    }

    # normalize the covariates for better numerical properties

    # can't do inverse normalization on resulting param estimates easily if have interactions or functions of covariates
    nm <- unique(c(all.names(locationFun), all.names(scaleFun), all.names(shapeFun)))
    if(any(!nm %in% c(names(x), '~', '+'))) .normalizeX <- FALSE
                    
    if(!is.null(x) && .normalizeX){
        nCovariates <- ncol(x)
        normalizers <- as.data.frame(matrix(0, 3, nCovariates))
        names(normalizers) <- names(x)
        for(k in 1:nCovariates){ # shift and scale to (-.5, .5) for better numeric properties in estimation
            normalizers[ , k] <- c(mean(x[ , k]), min(x[ , k]), max(x[ , k]))
            if(!is.null(xNew))
                xNew[ , k] <- normalize(xNew[, k], normalizers[1, k], normalizers[2, k], normalizers[3, k])
            if(!is.null(xContrast))
                xContrast[ , k] <- normalize(xContrast[, k], normalizers[1, k], normalizers[2, k], normalizers[3, k]) 
            x[ , k] <- normalize(x[ , k], normalizers[1, k], normalizers[2, k], normalizers[3, k])
        }
    } else normalizers <- NULL

    if(!is.null(weights)) {
        if(length(weights) != length(y)) {
            weights <- rep(weights, nReplicates)
            if(length(weights) != length(y))
                stop("fit_gev: length of 'weights' should equal the number of blocks or the number of blocks times 'nReplicates'.")
        }
    } else weights <- rep(1, length(y))
    
    if(!is.null(x) && nReplicates > 1 && nrow(x) != length(y)) {
        xOrig <- x
        x <- data.frame(X1 = rep(xOrig[ , 1], nReplicates))
        for(i in 2:ncol(xOrig))
            x[[paste0('X',i)]] <- rep(xOrig[ , i], nReplicates)
        names(x) <- names(xOrig)
        if(nrow(x) != length(y))
            stop("fit_gev: number of rows of 'x' cannot be matched to length of 'y'.")
    }

    # remove missing values
    if(!is.null(missingFlag)){
        y[y == missingFlag] <- NA
        if(!is.null(x) && (sum(x == missingFlag) || sum(is.na(x)))) 
            stop("fit_gev: missing values not allowed in 'x'; please remove observations with missing 'x' values.")
        if(!is.null(xNew) && (sum(xNew == missingFlag) || sum(is.na(xNew)))) 
            stop("fit_gev: missing values not allowed in 'xNew'.")
        if(!is.null(xContrast) && (sum(xContrast == missingFlag) || sum(is.na(xContrast)))) 
            stop("fit_gev: missing values not allowed in 'xContrast'.")
    }
    yNums <- !is.na(y)
    if(sum(yNums) < length(y)) {
        if(!is.null(replicateIndex))
            replicateIndex <- replicateIndex[yNums]
        if(!is.null(weights))
            weights <- weights[yNums]
        if(!is.null(x))
            x <- x[yNums, ]
        y <- y[yNums]    
    }

    if(!is.null(initial)) {
        if(!is.list(initial) || names(initial) != c('location', 'scale', 'shape'))
            stop("fit_gev: 'initial' must be a named list with components 'location', 'scale', 'shape'.")
        initial <- unlist(initial)
    }
    
    fit <- fit_model_gev(y, x, locationFun, scaleFun, shapeFun, weights, optimArgs, initial)
    
    if(!fit$info$failure) {
        # calculate various return quantities as specified by user
        if(is.null(xNew)) xUse <- x else xUse <- xNew
        results <- compute_return_quantities(fit, returnPeriod = returnPeriod, returnValue = returnValue, x = xUse, x2 = xContrast, locationFun = locationFun, scaleFun = scaleFun, shapeFun = shapeFun, getParams = getParams, upper = maxes, scaling = scaling, normalizers = normalizers)
        if(.getInputs) results$inputs <- list(y = y, x = x,
                                              locationFun = locationFun, scaleFun = scaleFun, shapeFun = shapeFun,
                                              nReplicates = nReplicates, replicateIndex = replicateIndex,
                                              weights = weights, returnPeriod = returnPeriod, returnValue = returnValue,
                                              xNew = xNew, xContrast = xContrast, m = m, mNew = mNew, mContrast = mContrast,
                                              normalizers = normalizers)
    } else {
        results <- list()
        warning("fit_gev: optimization failed; see 'info' in returned object. You may also want to run 'fit_gev()' with argument 'getFit = TRUE' to see more details regarding the attempted optimization.")
    }
    
    results$info <- fit$info
    if(getFit) {
        fit$info <- NULL
        results$fit <- fit
    }
    
    # do bootstrapping if requested
    if(bootSE && !results$info$failure) {
        nres <- max(1, ifelse(is.null(xNew), m, mNew))

        if(length(bootControl$seed) == 1){
            set.seed(bootControl$seed)
        } else{
            .Random.seed <- bootControl$seed
        }

        if(getParams) {
            mle_boot <- matrix(NA, bootControl$n, length(results$mle))
            dimnames(mle_boot)[[2]] <- names(results$mle)
        }
        if(!is.null(returnPeriod))
            rv_boot <- matrix(NA, bootControl$n, nres)
        if(!is.null(returnValue))
            logrp_boot <- matrix(NA, bootControl$n, nres)
        if(!is.null(xContrast)) {
            if(!is.null(returnPeriod))
                rvDiff_boot <- matrix(NA, bootControl$n, mContrast)
            if(!is.null(returnValue))
                logrpDiff_boot <- matrix(NA, bootControl$n, mContrast)
        }

        failures <- 0

        nBlocks <- length(y) / nReplicates
        
        for(b in 1:bootControl$n){
            
            # resample from the original dataset
            blockIndex = seq_len(nBlocks)
            yDf <- data.frame(y = y, weights = weights)

            if(bootControl$by == 'block') {
                nSmp <- nBlocks
                yDf$.boot <- blockIndex
            } else if(bootControl$by == 'replicate') {
                nSmp <- nReplicates
                yDf$.boot <- replicateIndex
            } else if(bootControl$by == 'joint') {
                nSmp <- nBlocks * nReplicates
                yDf$.boot <- blockIndex + (replicateIndex-1)*nBlocks
            } else stop("fit_gev: 'by' value of 'bootControl' not recognized.")

            smpDf <- data.frame(.boot = sample(1:nSmp, nSmp, replace = TRUE))

            if(!is.null(x)) {
                xDf <- x
                xDf$.boot <- yDf$.boot
                bX <- merge(smpDf, xDf, by.x = '.boot', by.y = '.boot',
                               all.x = FALSE, all.y = FALSE)
                bX$.boot <- NULL
            } else bX <- NULL

            bData <- merge(smpDf, yDf, by.x = '.boot', by.y = '.boot', all.x = FALSE, all.y = FALSE)

            
            # fit model to bootstrapped dataset
            bFit <- fit_model_gev(bData$y, bX, locationFun, scaleFun, shapeFun, bData$weights, optimArgs, initial = results$mle_raw)
            if(!bFit$info$failure) {
                if(is.null(xNew)) xUse <- x else xUse <- xNew
                bResults <- compute_return_quantities(bFit, returnPeriod = returnPeriod, returnValue = returnValue, x = xUse, x2 = xContrast, locationFun = locationFun, scaleFun = scaleFun, shapeFun = shapeFun, getParams = getParams, upper = maxes, scaling = scaling, normalizers = normalizers)
                if(getParams) 
                    mle_boot[b, ] <- bResults$mle
                if(!is.null(returnPeriod)) 
                    rv_boot[b, ] <- bResults$returnValue
                if(!is.null(returnValue)) 
                    logrp_boot[b, ] <- bResults$logReturnProb
                if(!is.null(xContrast)) {
                    if(!is.null(returnPeriod))
                        rvDiff_boot[b, ]  <- bResults$returnValueDiff
                    if(!is.null(returnValue))
                        logrpDiff_boot[b, ] <- bResults$logReturnProbDiff
                }
            } else {
                failures <- failures + 1
            }
            
        } # end of for(b in 1:bootControl$n){
        
        # calculate bootstrap results across the bootstrapped datasets
        if(getParams) 
            results$se_mle_boot <- apply(mle_boot, 2, sd, na.rm = TRUE)
        if(!is.null(returnPeriod)) 
            results$se_returnValue_boot <- apply(rv_boot, 2, sd, na.rm = TRUE)
        if(!is.null(returnValue)) {
            results$se_logReturnProb_boot <- apply(logrp_boot, 2, sd, na.rm = TRUE)
            results$se_logReturnPeriod_boot <- results$se_logReturnProb_boot
        }
        if(!is.null(xContrast)) {
            if(!is.null(returnPeriod))
                results$se_returnValueDiff_boot <- apply(rvDiff_boot, 2, sd, na.rm = TRUE)
            if(!is.null(returnValue)) 
                results$se_logReturnProbDiff_boot <- apply(logrpDiff_boot, 2, sd, na.rm = TRUE)
        }
        if(failures) {
            results$numBootFailures <- failures
            warning("fit_gev: Some bootstrap fits failed; see 'numBootFailures'.")
        }
        
        if(bootControl$getSample) {
            if(getParams)
                results$mle_boot <- mle_boot
            if(!is.null(returnPeriod)) 
                results$returnValue_boot <- rv_boot
            if(!is.null(returnValue)) {
                results$logReturnProb_boot <- logrp_boot
                results$logReturnPeriod_boot <- -logrp_boot
            }
            if(!is.null(xContrast)) {
                if(!is.null(returnPeriod))
                    results$returnValueDiff_boot <- rvDiff_boot
                if(!is.null(returnValue))
                    results$logReturnProbDiff_boot <- logrpDiff_boot
            }
        }
        
    } # end if(bootSE)
    results$mle_raw <- NULL
    
    return(results)
} # end of fit_gev()

fit_model_gev <- function(y, x, locationFun, scaleFun, shapeFun, weights, optimArgs, initial = NULL) {
    # auxiliary function for fitting GEV model to individual dataset - used for original dataset and bootstrapped datasets

    # set up initial values; used for bootstrapped dataset fitting based on fit to original data
    if(!is.null(initial)) {
        tmp <- list()
        tmp$location <- initial[grep("(location)|(mu)", names(initial))]
        tmp$scale <- tmp$log.scale <- initial[grep("(scale)|(phi)", names(initial))]
        tmp$shape <- initial[grep("(shape)|(xi)", names(initial))]
        initial <- tmp
    } 

    if(is.null(x)) {
        x <- data.frame(.y = y)
    } else {
        x$.y <- y
    }
   
    
    fit <- try(fevd(.y~1, data = x, location.fun = locationFun, scale.fun = scaleFun, shape.fun = shapeFun, use.phi = !check.constant(scaleFun), type = "GEV", method = "MLE", weights = weights, initial = initial, optim.args = optimArgs))
    
    failed <- TRUE
    fit$info <- list()
                                        # check for various errors in fitting
    if(!methods::is(fit, 'try-error')) {
        fit$info <- fit$results[c('convergence', 'counts', 'message')]
        if(!is.na(fit$results$convergence) && !fit$results$convergence) {
            shapeParams <- grep("xi|shape", fit$parnames)
            summ <- try(summary(fit, silent = TRUE))
            if(!methods::is(summ, 'try-error')) {
                if(!(sum(is.na(fit$results$par)) || !is.element('se.theta', names(summ)) || sum(is.na(summ$se.theta)) || min(summ$se.theta) < 1e-5))
                    failed <- FALSE else fit$info$message <- "fitting produced missing or invalid estimates or standard errors"
                if(fit$results$num.pars[['shape']] == 1 && fit$results$par[shapeParams] < -1) 
                    fit$info$message <- "shape parameter estimated to be less than -1"
            } else fit$info$message <- "fitting produced unknown error in fit summary"
        } else fit$info$message <- "optimization did not converge"
    } else { # try-error
        fit$info$convergence <- 1
        fit$info$message <- "optimization call returned error"
    }
    fit$info$failure <- failed
    return(fit)
} # end of fit_model_gev()

