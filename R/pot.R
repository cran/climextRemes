#' Fit a peaks-over-threshold model to exceedances over a threshold
#'
#' Fit a peaks-over-threshold model, designed specifically for climate data, to exceedance-only data, using the point process approach. Any covariates/predictors/features assumed to vary only between and not within blocks of observations. It includes options for variable weights (useful for local likelihood) and variable proportions of missing data, as well as for bootstrapping to estimate uncertainties. Results can be returned in terms of parameter values, return values, return periods, return probabilities, and differences in either return values or log return probabilities (i.e., log risk ratios). 
#'
#' @name fit_pot
#' 
#' @param y a numeric vector of exceedance values (values of the outcome variable above the threshold). For better optimization performance, it is recommended that the \code{y} have magnitude around one (see \code{Details}), for which one can use \code{scaling}.
#' @param x a data frame, or object that can be converted to a data frame with columns corresponding to covariate/predictor/feature variables and each row containing the values of the variable for a block (e.g., often a year with climate data). The number of rows must equal the number of blocks.
#' @param threshold a single numeric value for constant threshold or a numeric vector with length equal to the number of blocks, indicating the threshold for each block.
#' @param locationFun formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the location parameter using columns from \code{x}.  \code{x} must be supplied if this is anything other than NULL or ~1.
#' @param scaleFun formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the (potentially transformed) scale parameter using columns from \code{x}.  \code{x} must be supplied if this is anything other than NULL or ~1. \code{logScale} controls whether this determines the log of the scale or the scale directly.
#' @param shapeFun formula, vector of character strings, or indices describing a linear model (i.e., regression function) for the shape parameter using columns from \code{x}.  \code{x} must be supplied if this is anything other than NULL or ~1.
#' @param nBlocks number of blocks (e.g., a block will often be a year with climate data); note this value determines the interpretation of return values/periods/probabilities; see \code{returnPeriod} and \code{returnValue}.
#' @param blockIndex numeric vector providing the index of the block corresponding to each element of \code{y}. Used only when \code{x} is provided to match exceedances to the covariate/predictor/feature value for the exceedance or when using bootstrapping with the resampling based on blocks based on the \code{by} element of \code{bootControl}. If \code{firstBlock} is not equal to one, then \code{blockIndex} need not have one as its smallest possible value.
#' @param firstBlock single numeric value indicating the numeric value of the first possible block of \code{blockIndex}. For example the values in \code{blockIndex} might indicate the year of each exceedance with the first year of data being 1969, in which case \code{firstBlock} would be 1969. Note that the first block may not have any exceedances so it may not be represented in \code{blockIndex}. Used only to adjust \code{blockIndex} so that the block indices start at one and therefore correspond to the rows of \code{x}.
#' @param index numeric vector providing the integer-valued index (e.g., julian day for daily climate data) corresponding to each element of \code{y}. For example if there are 10 original observations and the third, fourth, and seventh values are exceedances, then \code{index} would be the vector 3,4,7. Used only when \code{declustering} is provided to determine which exceedances occur sequentially or within a contiguous set of values of a given length. The actual values are arbitrary; only the lags between the values are used.
#' @param nReplicates numeric value indicating the number of replicates.
#' @param replicateIndex numeric vector providing the index of the replicate corresponding to each element of \code{y}. Used for three purposes: (1) when using bootstrapping with the resampling based on replicates based on the \code{by} element of \code{bootControl}, (2) to avoid treating values in different replicates as potentially being sequential or within a short interval when removing values based on \code{declustering}, and (3) to match outcomes to \code{weights} or \code{proportionMissing} when either vary by replicate.  
#' @param weights a vector or matrix providing the weights by block. When there is only one replicate or the weights do not vary by replicate, a vector of length equal to the number of blocks. When weights vary by replicate, a matrix with rows corresponding to blocks and columns to replicates. Likelihood contribution of each block is multiplied by the corresponding weight. 
#' @param proportionMissing a numeric value, vector or matrix indicating the proportion of missing values in the original dataset before exceedances were selected. When the proportion missing is the same for all blocks and replicates, a single value. When there is only one replicate or the weights do not vary by replicate, a vector of length equal to the number of blocks. When weights vary by replicate, a matrix with rows corresponding to blocks and columns to replicates.
#' @param returnPeriod numeric value giving the number of blocks for which return values should be calculated. For example a \code{returnPeriod} equal to 20 will result in calculation of the value of an event that occurs with probability 1/20 in any block and therefore occurs on average every 20 blocks. Often blocks will correspond to years.
#' @param returnValue numeric value giving the value for which return probabilities/periods should be calculated, where the resulting period will be the average number of blocks until the value is exceeded and the probability the probability of exceeding the value in any single block.
#' @param getParams logical indicating whether to return the fitted parameter values and their standard errors; WARNING: parameter values for models with covariates for the scale parameter must interpreted based on the value of \code{logScale}.
#' @param getFit logical indicating whether to return the full fitted model (potentially useful for model evaluation and for understanding optimization problems); note that estimated parameters in the fit object for nonstationary models will not generally match the MLE provided when \code{getParams} is \code{TRUE} because covariates are normalized before fitting and the fit object is based on the normalized covariates. Similarly, parameters will not match if \code{scaling} is not 1. 
#' @param xNew object of the same form as \code{x}, providing covariate/predictor/feature values for which return values/periods/probabilities are desired.
#' @param xContrast object of the same form and dimensions as \code{xNew}, providing covariate/predictor/feature values for which to calculate the differences of the return values and/or log return probabilities relative to the values in \code{xNew}. This provides a way to estimate differences in return value or log return probabilities (i.e., log risk ratios).
#' @param declustering one of \code{NULL}, \code{"noruns"}, or a number. If 'noruns' is specified, only the maximum (or minimum if upperTail = FALSE) value within a set of exceedances corresponding to successive indices is included. If a number, this should indicate the size of the interval (which will be used with the \code{index} argument) within which to allow only the largest (or smallest if upperTail = FALSE) value.
#' @param upperTail logical indicating whether one is working with exceedances over a high threshold (TRUE) or exceedances under a low threshold (FALSE); in the latter case, the function works with the negative of the values and the threshold, changing the sign of the resulting location parameters.
#' @param scaling positive-valued scalar used to scale the data values for more robust optimization performance. When multiplied by the values, it should produce values with magnitude around 1.
#' @param bootSE logical indicating whether to use the bootstrap to estimate standard errors.
#' @param bootControl a list of control parameters for the bootstrapping. See \sQuote{Details}.
#' @param optimArgs a list with named components matching exactly any arguments that the user wishes to pass to R's \code{optim} function. See \code{help(optim)} for details. Of particular note, \code{'method'} can be used to choose the optimization method used for maximizing the log-likelihood to fit the model and \code{control = list(maxit=VALUE)} for a user-specified VALUE can be used to increase the number of iterations if the optimization is converging slowly.
#' @param optimControl a list with named components matching exactly any elements that the user wishes to pass as the \code{control} argument to R's \code{optim} function. See \code{help(optim)} for details. Primarily provided for the Python interface because \code{control} can also be passed as part of \code{optimArgs}.
#' @param initial a list with components named \code{'location'}, \code{'scale'}, and \code{'shape'} providing initial parameter values, intended for use in speeding up or enabling optimization when the default initial values are resulting in failure of the optimization; note that use of \code{scaling}, \code{logScale} and \code{.normalizeX = TRUE} cause numerical changes in some of the parameters. For example with \code{logScale = TRUE}, initial value(s) for \code{'scale'} should be specified on the log scale.
#' @param logScale logical indicating whether optimization for the scale parameter should be done on the log scale. By default this is \code{FALSE} when the scale is not a function of covariates and \code{TRUE} when the scale is a function of covariates (to ensure the scale is positive regardless of the regression coefficients). 
#' @param .normalizeX logical indicating whether to normalize \code{x} values for better numerical performance; default is \code{TRUE}.
#' @param .getInputs logical indicating whether to return intermediate objects used in fitting. Defaults to \code{FALSE} and intended for internal use only
#' @param .allowNoInt logical indicating whether no-intercept models are allowed. Defaults to \code{TRUE} and provided primarily to enable backwards compatibility with versions <= 0.2.2.
#' @author Christopher J. Paciorek
#' @export
#' @details
#' This function allows one to fit stationary or nonstationary peaks-over-threshold models using the point process approach. The function can return parameter estimates (which are asymptotically equivalent to GEV model parameters for block maxima data), return value/level for a given return period (number of blocks),  and return probabilities/periods for a given return value/level. The function provides standard errors based on the usual MLE asymptotics, with delta-method-based standard errors for functionals of the parameters, but also standard errors based on the nonparametric bootstrap, either resampling by block or by replicate or both.
#'
#' Analyzing aggregated observations, such as yearly averages:
#' 
#' Aggregated average or summed data such as yearly or seasonal averages can be fit using this function. The best way to do this is to specify \code{nBlocks} to be the number of observations (i.e., the length of the observation period, not the number of exceedances). Then any return probabilities can be interpreted as the probabilities for a single block (e.g., a year). If instead \code{nBlocks} were one (i.e., a single block) then probabilities would be interpreted as the probability of occurrence in a multi-year block. 
#'
#' Blocks and replicates:
#' 
#' Note that blocks and replicates are related concepts. Blocks are the grouping of values such that return values and return periods are calculated based on the equivalent block maxima (or minima) generalized extreme value model. For example if a block is a year, the return period is the average number of years before the given value is seen, while the return value when \code{returnPeriod} is, say, 20, is the value exceeded on average once in 20 years. A given dataset will generally have multiple blocks. In some cases a block may contain only a single value, such as when analyzing yearly sums or averages.
#'
#' Replicates are repeated datasets, each with the same structure, including the same number of blocks. The additional blocks in multiple replicates could simply be treated as additional blocks without replication, but when the predictor variables and weights are the same across replicates, it is simpler to make use of \code{nReplicates} and \code{replicateIndex} (see next paragraph). A given replicate might only contain a single block, such as with an ensemble of short climate model runs that are run only for the length of a single block (e.g., a single year). In this case \code{nBlocks} should be set to one.
#'
#' When using multiple replicates (e.g., multiple members of a climate model initial condition ensemble), the standard input format is to append outcome values for additional replicates to the \code{y} argument and indicate the replicate ID for each exceedance in \code{replicateIndex}. However, if one has different covariate values or thresholds for different replicates, then one needs to treat the additional replicates as providing additional blocks, with only a single replicate. The covariate values can then be included as additional rows in \code{x}, with \code{nBlocks} reflecting the product of the number of blocks per replicate and the number of replicates and \code{nReplicates} set to 1. In this situation, if \code{declustering} is specified, make sure to set \code{index} such that the values for separate replicates do not overlap with each other, to avoid treating exceedances from different replicates as being sequential or from a contiguous set of values. Similarly, if there is a varying number of replicates by block, then all block-replicate pairs should be treated as individual blocks with a corresponding row in \code{x}.
#' 
#' \code{bootControl} arguments:
#' 
#' The \code{bootControl} argument is a list (or dictionary when calling from Python) that can supply any of the following components:
#' \itemize{
#' \item seed. Value of the random number seed as a single value, or in the form of \code{.Random.seed}, to set before doing resampling. Defaults to \code{1}.
#' \item n. Number of bootstrap samples. Defaults to \code{250}.
#' \item by. Character string, one of \code{'block'}, \code{'replicate'}, or \code{'joint'}, indicating the basis for the resampling. If \code{'block'}, resampled datasets will consist of blocks drawn at random from the original set of blocks; if there are replicates, each replicate will occur once for every resampled block. If \code{'replicate'}, resampled datasets will consist of replicates drawn at random from the original set of replicates; all blocks from a replicate will occur in each resampled replicate. Note that this preserves any dependence across blocks rather than assuming independence between blocks. If \code{'joint'} resampled datasets will consist of block-replicate pairs drawn at random from the original set of block-replicate pairs. Defaults to \code{'block'}. 
#' \item getSample. Logical/boolean indicating whether the user wants the full bootstrap sample of parameter estimates and/or return value/period/probability information provided for use in subsequent calculations; if FALSE (the default), only the bootstrap-based estimated standard errors are returned.
#' }
#'
#' Optimization failures:
#'
#' It is not uncommon for maximization of the log-likelihood to fail for extreme value models. Users should carefully check the \code{info} element of the return object to ensure that the optimization converged. For better optimization performance, it is recommended that the observations be scaled to have magnitude around one (e.g., converting precipitation from mm to cm). When there is a convergence failure, one can try a different optimization method, more iterations, or different starting values -- see \code{optimArgs} and \code{initial}. In particular, the Nelder-Mead method is used; users may want to try the BFGS method by setting \code{optimArgs = list(method = 'BFGS')} (or \code{optimArgs = {'method': 'BFGS'}} if calling from Python).
#' 
#' When using the bootstrap, users should check that the number of convergence failures when fitting to the boostrapped datasets is small, as it is not clear how to interpret the bootstrap results when there are convergence failures for some bootstrapped datasets. 
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
#' 
#' Paciorek, C.J., D.A. Stone, and M.F. Wehner. 2018. Quantifying uncertainty in the attribution of human influence on severe weather. Weather and Climate Extremes 20:69-80. arXiv preprint <https://arxiv.org/abs/1706.03388>.
#'
#' @examples
#' # setup Fort precipitation data
#' data(Fort, package = 'extRemes')
#' firstBlock <- min(Fort$year)
#' years <- min(Fort$year):max(Fort$year)
#' nYears <- length(years)
#' threshold <- 0.395
#' ord <- order(Fort$year, Fort$month, Fort$day) 
#' Fort <- Fort[ord, ]
#' ind <- Fort$Prec > threshold
#' FortExc <- Fort[ind, ]
#'
#' # stationary fit
#' out <- fit_pot(FortExc$Prec, threshold = threshold, nBlocks = nYears, 
#'         returnPeriod = 20, returnValue = 3.5,
#'         getParams = TRUE, bootSE = FALSE)
#'
#' # fit with location linear in year
#' out <- fit_pot(FortExc$Prec, x = data.frame(years = years), threshold = threshold,
#'         locationFun = ~years, nBlocks = nYears, 
#'         blockIndex = FortExc$year, firstBlock = firstBlock,
#'         returnPeriod = 20, returnValue = 3.5,
#'         getParams = TRUE, xNew = data.frame(years = range(Fort$year)), bootSE = FALSE)
#'
#' # with declustering
#' index <- seq_len(nrow(Fort))
#' out <- fit_pot(FortExc$Prec, x = data.frame(years = years), threshold = threshold,
#'         locationFun = ~years, nBlocks = nYears, 
#'         blockIndex = FortExc$year, firstBlock = firstBlock, index = index[ind],
#'         returnPeriod = 20, returnValue = 3.5,
#'         getParams = TRUE, xNew = data.frame(years = range(Fort$year)),
#'         declustering = 'noruns', bootSE = FALSE)
#'
#' # with replicates; for illustration here, I'll just duplicate the Fort data
#' Fort2x <- rbind(FortExc, FortExc)
#' n <- nrow(FortExc)
#' out <- fit_pot(Fort2x$Prec, x = data.frame(years = years), threshold = threshold,
#'         locationFun = ~years, nBlocks = nYears,
#'         blockIndex = Fort2x$year, firstBlock = firstBlock,
#'         nReplicates = 2, replicateIndex = c(rep(1, n), rep(2, n)),
#'         returnPeriod = 20, returnValue = 3.5,
#'         getParams = TRUE, xNew = data.frame(years = range(Fort$year)), bootSE = FALSE)
#'
#' # analysis of seasonal total precipitation
#' FortSummer <- Fort[Fort$month %in% 6:8, ]  # summer precipitation
#' FortSummerSum <- aggregate(Prec ~ year, data = FortSummer, sum)
#' threshold <- quantile(FortSummerSum$Prec, 0.8)
#' FortSummerSumExc <- FortSummerSum[FortSummerSum$Prec > threshold, ]
#' # each year (single observation) treated as a block, so return probability
#' # can be interpreted as probability of exceeding a value in a single year
#' out <- fit_pot(FortSummerSumExc$Prec, x = data.frame(years = years), threshold = threshold,
#'         nBlocks = nYears, blockIndex = FortSummerSumExc$year, firstBlock = firstBlock,
#'         locationFun = ~years, returnPeriod = 20,
#'         returnValue = 10, getParams = TRUE, xNew = data.frame(years = range(Fort$year)),
#'         bootSE = FALSE)
#'
fit_pot <- function(y, x = NULL, threshold, locationFun = NULL, scaleFun = NULL,
                    shapeFun = NULL, nBlocks = nrow(x), blockIndex = NULL, firstBlock = 1,
                    index = NULL, nReplicates = 1, replicateIndex = NULL,
                    weights = NULL, proportionMissing = NULL, returnPeriod = NULL, returnValue = NULL,
                    getParams = FALSE, getFit = FALSE,
                    xNew = NULL, xContrast = NULL, declustering = NULL, upperTail = TRUE, scaling = 1,
                    bootSE = FALSE, bootControl = NULL,
                    optimArgs = NULL, optimControl = NULL, initial = NULL, logScale = NULL, .normalizeX = TRUE, .getInputs = FALSE, .allowNoInt = TRUE) {

    ## various input checks
  
    if(is.null(y))  # for Python interface
        stop("fit_pot: argument 'y' is missing, with no default")

    if(is.null(optimArgs))
        optimArgs = list(method = c("Nelder-Mead"))
    if(is.null(optimArgs$method))
        optimArgs$method <- c("Nelder-Mead")
    if(!is.null(optimControl))
        optimArgs$control <- optimControl

    if(bootSE) {
        bControl <- list(seed = 1, n = 250, by = "block", getSample = FALSE,
                         getSampleSE = FALSE)
        bControl[names(bootControl)] <- bootControl
        bootControl <- bControl
    }

    numNumericFun <- sum(c(is.null(locationFun) || is.numeric(locationFun),
             is.null(scaleFun) || is.numeric(scaleFun),
             is.null(shapeFun) || is.numeric(shapeFun)))
    
    if(!is.null(x)) {
        x <- try(as.data.frame(x))
        if(is(x, 'try-error')) stop("fit_pot: 'x' should be a data frame or be able to be converted to a data frame.")
        m <- nrow(x)
        if(any(is.na(x)))
            stop("fit_pot: 'x' should not contain any NAs.")
    } else m <- 0
    if(!is.null(xNew)) {
        if(is.null(x))
            stop("fit_pot: 'x' must be provided if 'xNew' is provided.")
        xNew <- try(as.data.frame(xNew))
        if(is(xNew, 'try-error')) stop("fit_pot: 'xNew' should be a data frame or be able to be converted to a data frame.")
        mNew <- nrow(xNew)
        if(ncol(x) != ncol(xNew) || (numNumericFun != 3 && !identical(sort(names(x)), sort(names(xNew)))))
            stop("fit_pot: columns in 'x' and 'xNew' should be the same.")
    } else mNew <- 0
    if(!is.null(xContrast)) {
        if(is.null(x))
            stop("fit_pot: 'x' must be provided if 'xNew' is provided.")
        xContrast <- try(as.data.frame(xContrast))
        if(is(xContrast, 'try-error')) stop("fit_pot: 'xContrast' should be a data frame or be able to be converted to a data frame.")
        mContrast <- nrow(xContrast)
        if(is.null(xNew) && m != mContrast) stop("fit_pot: number of values in 'x' and 'xContrast' should be the same.")
        if(!is.null(xNew) && mNew != mContrast) stop("fit_pot: number of values in 'xNew' and 'xContrast' should be the same.")
        if(ncol(xContrast) != ncol(x) || (numNumericFun != 3 && !identical(sort(names(x)), sort(names(xContrast)))))
            stop("fit_pot: columns in 'x' and 'xContrast' should be the same.")
    } else mContrast <- 0

    if(m > 0 && is.null(nBlocks))
        nBlocks <- m

    if(is.null(nBlocks))
        stop("fit_pot: 'nBlocks' must be provided.")
    
    if(m > 0 && m != nBlocks) {
        stop("fit_pot: number of blocks must equal number of rows of 'x'.")
    }
    
    if(!is.null(blockIndex))
        if(length(blockIndex) != length(y)) stop("fit_pot: length of 'blockIndex' must match length of 'y'.")
   
    if(is.null(blockIndex) && ((bootSE && !bootControl$by == 'replicate') || !is.null(x) || !is.null(weights) || length(threshold) > 1))
        stop("fit_pot: 'blockIndex' must be provided if 'x' or 'weights' is provided or 'threshold' is not constant or using bootstrap by 'block' or 'joint'.")

    if(nReplicates > 1 && is.null(replicateIndex) && (
        (bootSE && bootControl$by == 'replicate') || !is.null(declustering) || length(weights) == nBlocks*nReplicates))
        stop("fit_pot: 'replicateIndex' must be provided if 'nReplicates' is greater than 1 when bootstrapping by replicate, using 'declustering' or weights vary by replicate.")

    if(!(is.null(weights) || length(weights) %in% c(1, nBlocks, nBlocks*nReplicates)))
        stop("fit_pot: 'weights' must be a single value, a vector of length equal to the number of blocks or a matrix with as many rows as the number of blocks and as many columns as the number of replicates.")
 
    if(!(is.null(proportionMissing) || length(proportionMissing) %in% c(1, nBlocks, nBlocks*nReplicates)))
        stop("fit_pot: 'proportionMissing' must be a single value, a vector of length equal to the number of blocks or a matrix with as many rows as the number of blocks and as many columns as the number of replicates.")
    
    if(!length(threshold) %in% c(1, nBlocks))
        stop("fit_pot: 'threshold' should be a single value or one value per block.")
    names(threshold) <- NULL # often has a name when created via quantile()
    
    if(!is.numeric(y)) stop("fit_pot: 'y' should be a numeric vector.")
    

    if(!is.null(blockIndex)) {
        # set blockIndex to start at 1
        blockIndex <- blockIndex - firstBlock + 1
        if(any(is.na(blockIndex)))
            stop("fit_pot: 'blockIndex' must not contain any NAs.")
        if(min(blockIndex) < 1 || max(blockIndex) > nBlocks)
            stop("fit_pot: 'blockIndex' (after adjusting relative to 'firstBlock' if it is provided) should give indices between 1 and the number of blocks.")
    }
        
    
    if(nReplicates == 1) replicateIndex <- rep(1, length(y))
    
    if(!is.null(initial)) {
        if(!is.list(initial) || !identical(sort(names(initial)), c('location', 'scale', 'shape')))
            stop("fit_pot: 'initial' must be a named list with components 'location', 'scale', 'shape'.")
        if(!upperTail)
            initial$location <- -initial$location
        initial <- unlist(initial)
    }

    # lower tail 'exceedances' are handled by using negatives of all values
    if(!upperTail){  
        y <- -y
        threshold <- -threshold
        if(!is.null(returnValue))
            returnValue <- -returnValue
    }
    

    # ensure parameter functions are in the form of formulae
    locationFun <- parseParamInput(locationFun, names(x), .allowNoInt)
    scaleFun <- parseParamInput(scaleFun, names(x), .allowNoInt)
    shapeFun <- parseParamInput(shapeFun, names(x), .allowNoInt)
    
    if(numNumericFun == 3) {
        # if using numeric indices of columns, align the names
        if(!is.null(xNew))
            names(xNew) <- names(x)
        if(!is.null(xContrast))
            names(xContrast) <- names(x)
    }

    if(!is.null(x) && locationFun == quote(~1) && scaleFun == quote(~1) && shapeFun == quote(~1)) { 
        warning("fit_pot: 'x' provided but fitting stationary model based on 'locationFun', 'scaleFun' and 'shapeFun' specification.")
        x <- xNew <- xContrast <- NULL
        m <- mNew <- mContrast <- 0
    }

    
    # scale data for better numerical performance (provided user chose scaling appropriately)
    if(scaling != 1) {
        y <- y * scaling
        threshold <- threshold * scaling
        if(!is.null(returnValue))
            returnValue <- returnValue * scaling
    }



    # normalize the covariates for better numerical properties

    # can't do inverse normalization on resulting param estimates easily if have interactions or functions of covariates
    nm <- unique(c(all.names(locationFun), all.names(scaleFun), all.names(shapeFun)))
    if(any(!nm %in% c(names(x), '~', '+', '-'))) .normalizeX <- FALSE

    noInt <- c('-' %in% all.names(locationFun), '-' %in% all.names(scaleFun), '-' %in% all.names(shapeFun))
    
    nm <- nm[grep("^[a-zA-Z]", nm)]
    if(!is.null(x) && !all(nm %in% names(x)))
        stop("fit_pot: one or more variables in location, scale, and/or shape formulae are not contained in 'x'.")
    if(!is.null(xNew) && !all(nm %in% names(xNew)))
        stop("fit_pot: one or more variables in location, scale, and/or shape formulae are not contained in 'xNew'.")
    if(!is.null(xContrast) && !all(nm %in% names(xContrast)))
        stop("fit_pot: one or more variables in location, scale, and/or shape formulae are not contained in 'xContrast'.")
    
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

    # remove exceedances occurring sequentially or multiple exceedances within a fixed interval
    if(!is.null(declustering)) {
        if(is.null(index))
            stop("fit_pot: 'index' required when 'declustering' is specified.")
        for(j in 1:nReplicates)
            if(declustering == "noruns") {
                y[replicateIndex == j] = remove_runs(y[replicateIndex == j], index[replicateIndex == j])
            } else {
                if(is.numeric(declustering)) 
                    y[replicateIndex == j] <- screen_within_block(y[replicateIndex == j], index[replicateIndex == j], blockLength = declustering)
                                        # alternative is to create complete set of exceedances and nonexceedances for continuous days and use extRemes::decluster(y, index); then remove non-exceedances; this is a bit roundabout however
                else stop("fit_pot: invalid input for 'declustering'.")
            }
    }

    # remove missing values
    yNum <- !is.na(y)
    if(sum(yNum) < length(y)) {
        if(!is.null(index))
            index <- index[yNum]
        if(!is.null(blockIndex))
            blockIndex <- blockIndex[yNum]
        if(!is.null(replicateIndex))
            replicateIndex <- replicateIndex[yNum]
        y <- y[yNum]    
    }
    

    # get x, weight, threshold values corresponding to individual exceedances by
    #  matching based on 'blockIndex'
    if(!is.null(x)) 
        xObs <- x[blockIndex, , drop = FALSE] else xObs <- NULL

    if(!is.null(weights)) {
        if(is.array(weights) && length(dim(weights)) == 1) # this will be the case for Python interface
            weights <- as.vector(weights)  
        if(nReplicates > 1 && length(weights) == nBlocks*nReplicates) {
            weightsObs <- weights[cbind(blockIndex, replicateIndex)]
        } else weightsObs <- weights[blockIndex]
    } else {
        weightsObs <- 1
        weights <- rep(1, nBlocks)
    }

    if(length(threshold) == 1) {
        thresholdObs <- threshold
        threshold <- rep(threshold, nBlocks)
    } else thresholdObs <- threshold[blockIndex]

    if(any(y <= thresholdObs)) stop("fit_pot: 'y' should only contain exceedances of 'threshold'.")

    
    if(!is.null(proportionMissing) && length(proportionMissing) == 1) 
        proportionMissing <- rep(proportionMissing, nBlocks)
    
    
    fit <- fit_model_pot(y, xObs, x, thresholdObs, threshold, locationFun, scaleFun, shapeFun, nBlocks, nReplicates, weightsObs, weights, proportionMissing, optimArgs, initial = initial, logScale = logScale, noInt = noInt)
    if(!fit$info$failure) {
        # calculate various return quantities as specified by user
        if(is.null(xNew)) xUse <- x else xUse <- xNew
        results <- compute_return_quantities(fit, returnPeriod = returnPeriod, returnValue = returnValue, x = xUse, x2 = xContrast, locationFun = locationFun, scaleFun = scaleFun, shapeFun = shapeFun, getParams = getParams, upper = upperTail, scaling = scaling, normalizers = normalizers, getSE = TRUE)
        if(.getInputs) results$inputs <- list(y = y, x = x, xObs = xObs, threshold = threshold, thresholdObs = thresholdObs,
                                              locationFun = locationFun, scaleFun = scaleFun, shapeFun = shapeFun,
                                              nBlocks = nBlocks, blockIndex = blockIndex,
                                              nReplicates = nReplicates, replicateIndex = replicateIndex,
                                              weights = weights, weightsObs = weightsObs, proportionMissing = proportionMissing,
                                              returnPeriod = returnPeriod, returnValue = returnValue,
                                              xNew = xNew, xContrast = xContrast,
                                              m = m, mNew = mNew, mContrast = mContrast,
                                              normalizers = normalizers)
    } else {
        results <- list()
        warning("fit_pot: optimization failed; see 'info' in returned object. You may also want to run 'fit_pot()' with argument 'getFit = TRUE' to see more details regarding the attempted optimization. You may also want to try a different optimization method or increase the maximum number of iterations.")
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
        if(!is.null(returnPeriod)) {
            rv_boot <- array(NA, c(bootControl$n, nres, length(returnPeriod)))
            dimnames(rv_boot) <- list(NULL, NULL, returnPeriod)
        }
        if(!is.null(returnValue)) {
            nm <- returnValue; if(!upperTail) nm <- -nm
            logrp_boot <- array(NA, c(bootControl$n, nres, length(returnValue)))
            dimnames(logrp_boot) <- list(NULL, NULL, nm/scaling)
            logrp_boot_se <- array(NA, c(bootControl$n, nres, length(returnValue)))
            dimnames(logrp_boot) <- list(NULL, NULL, nm/scaling)
        }
        if(!is.null(xContrast)) {
            if(!is.null(returnPeriod)) {
                rvDiff_boot <- array(NA, c(bootControl$n, mContrast, length(returnPeriod)))
                dimnames(rvDiff_boot) <- list(NULL, NULL, returnPeriod)
            }
            if(!is.null(returnValue)) {
                nm <- returnValue; if(!upperTail) nm <- -nm
                logrpDiff_boot <- array(NA, c(bootControl$n, mContrast, length(returnValue)))
                dimnames(logrpDiff_boot) <- list(NULL, NULL, nm/scaling)
            }
        }

        failures <- 0
        for(b in 1:bootControl$n){
            # resample from the original dataset
            yDf <- data.frame(y = y, weightsObs = weightsObs, thresholdObs = thresholdObs)

            if(bootControl$by == 'block') {
                nSmp <- nBlocks
                yDf$.boot <- blockIndex
            } else if(bootControl$by == 'replicate') {
                nSmp <- nReplicates
                yDf$.boot <- replicateIndex
            } else if(bootControl$by == 'joint') {
                nSmp <- nBlocks * nReplicates
                yDf$.boot <- blockIndex + (replicateIndex-1)*nBlocks
            } else stop("fit_pot: 'by' value of 'bootControl' not recognized.")
            smpDf <- data.frame(.boot = sample(1:nSmp, nSmp, replace = TRUE))

            if(!is.null(xObs)) {
                xDf <- xObs
                xDf$.boot <- yDf$.boot
                bxObs <- merge(smpDf, xDf, by.x = '.boot', by.y = '.boot',
                               all.x = FALSE, all.y = FALSE)
                bxObs$.boot <- NULL
            } else bxObs <- NULL
            
            bData <- merge(smpDf, yDf,
                           by.x = '.boot', by.y = '.boot', all.x = FALSE, all.y = FALSE)

            smp <- smpDf$.boot
            if(bootControl$by == 'block') {
                bX <- x[smp, , drop = FALSE]
                bThreshold <- threshold[smp]
                bWeights <- if(is.matrix(weights))
                                weights[smp, , drop = FALSE] else weights[smp]
                if(is.null(proportionMissing)) {
                    bProportionMissing <- NULL
                } else bProportionMissing <- if(is.matrix(proportionMissing))
                    proportionMissing[smp, , drop = FALSE] else proportionMissing[smp]
            }
            if(bootControl$by == 'replicate') {
                bX <- x
                bThreshold <- threshold
                bWeights <- if(is.matrix(weights))
                                weights[, smp, drop = FALSE] else weights
                if(is.null(proportionMissing)) {
                    bProportionMissing <- NULL
                } else bProportionMissing <- if(is.matrix(proportionMissing))
                    proportionMissing[ , smp, drop = FALSE] else proportionMissing
            }
            if(bootControl$by == 'joint') {
                blockSmp <- smp - nBlocks * ((smp-1) %/% nBlocks)
                bX <- x[blockSmp, , drop = FALSE]
                bThreshold <- threshold[blockSmp]
                bWeights <- if(is.matrix(weights)) 
                                weights[smp] else weights[blockSmp]
                if(is.null(proportionMissing)) {
                    bProportionMissing <- NULL
                } else                bProportionMissing <- if(is.matrix(proportionMissing))
                    proportionMissing[smp] else proportionMissing[blockSmp]
            }
         
            # fit model to bootstrapped dataset
            bFit <- fit_model_pot(bData$y, bxObs, bX, bData$thresholdObs, bThreshold, locationFun, scaleFun, shapeFun, nBlocks, nReplicates, bData$weightsObs, bWeights, bProportionMissing, optimArgs = optimArgs, initial = fit$results$par, logScale = logScale, noInt = noInt)
            if(!bFit$info$failure) {
                if(is.null(xNew)) xUse <- x else xUse <- xNew
                bResults <- compute_return_quantities(bFit, returnPeriod = returnPeriod, returnValue = returnValue, x = xUse, x2 = xContrast, locationFun = locationFun, scaleFun = scaleFun, shapeFun = shapeFun, getParams = getParams, upper = upperTail, scaling = scaling, normalizers = normalizers, getSE = bootControl$getSampleSE)
                if(getParams) 
                    mle_boot[b, ] <- bResults$mle
                if(!is.null(returnPeriod)) 
                    rv_boot[b, , ] <- bResults$returnValue
                if(!is.null(returnValue)) { 
                    logrp_boot[b, , ] <- bResults$logReturnProb
                    if(bootControl$getSampleSE)
                        logrp_boot_se[b, , ] <- bResults$se_logReturnProb
                }
                if(!is.null(xContrast)) {
                    if(!is.null(returnPeriod))
                        rvDiff_boot[b, , ]  <- bResults$returnValueDiff
                    if(!is.null(returnValue))
                        logrpDiff_boot[b, , ] <- bResults$logReturnProbDiff
                }
            } else {
                failures <- failures + 1
            }
        } # end of for(b in 1:bootControl$n){
        
        # calculate bootstrap results across the bootstrapped datasets
        if(getParams) 
            results$se_mle_boot <- apply(mle_boot, 2, sd, na.rm = TRUE)
        if(!is.null(returnPeriod)) 
            results$se_returnValue_boot <- drop(apply(rv_boot, c(2,3), sd, na.rm = TRUE))
        if(!is.null(returnValue)) {
            results$se_logReturnProb_boot <- drop(apply(logrp_boot, c(2,3), sd, na.rm = TRUE))
            results$se_logReturnPeriod_boot <- results$se_logReturnProb_boot
        }
        if(!is.null(xContrast)) {
            if(!is.null(returnPeriod))
                results$se_returnValueDiff_boot <- drop(apply(rvDiff_boot, c(2,3), sd, na.rm = TRUE))
            if(!is.null(returnValue)) 
                results$se_logReturnProbDiff_boot <- drop(apply(logrpDiff_boot, c(2,3), sd, na.rm = TRUE))
        }
        results$numBootFailures <- failures
        if(failures) warning("fit_pot: Some bootstrap fits failed; see 'numBootFailures'.")
        
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
        if(bootControl$getSampleSE) # used only for bootstrap-t for risk ratio
            results$logReturnProb_boot_se <- logrp_boot_se
        
    } # end if(bootSE)
    results$mle_raw <- NULL
    
    return(results)
} # end of fit_pot()
                         

fit_model_pot <- function(y, xObs, x, thresholdObs, threshold, locationFun, scaleFun, shapeFun, nBlocks, nReplicates, weightsObs, weights, proportionMissing, optimArgs, initial = NULL, logScale = NULL, noInt) {
    # auxiliary function for fitting POT model to individual dataset - used for original dataset and bootstrapped datasets

    # set up initial values; used for bootstrapped dataset fitting based on fit to original data
    if(!is.null(initial)) {
        tmp <- list()
        tmp$location <- initial[grep("(location)|(mu)", names(initial))]
        tmp$scale <- tmp$log.scale <- initial[grep("(scale)|(phi)|(sigma)", names(initial))]
        tmp$shape <- initial[grep("(shape)|(xi)", names(initial))]
        initial <- tmp
    } 

    if(is.null(xObs)) {
        xObs <- data.frame(.y = y)
    } else {
        xObs$.y <- y
    }
    if(length(y) < 2)
        warning("fit_pot: Fewer than two non-missing observations; optimization will fail.")
    
    blocks <- list(nBlocks = nBlocks * nReplicates)
    blocks$threshold <- threshold
    if(length(weights) != blocks$nBlocks) weights <- rep(weights, nReplicates)
    blocks$weights <- weights
    if(!is.null(proportionMissing)) {
        if(length(proportionMissing) != blocks$nBlocks) proportionMissing <- rep(proportionMissing, nReplicates)
        blocks$proportionMissing <- proportionMissing
    }
    if(!is.null(x)) {
        if(nReplicates == 1 || nrow(x) == blocks$nBlocks) {
            blocks$data <- x
        } else {
            # duplicate x with new rows for the replicates
            blocks$data <- data.frame(X1 = rep(x[ , 1], nReplicates))
            if(ncol(x) > 1) {
                for(k in 2:ncol(x))
                    blocks$data[[paste0('X',k)]] <- rep(x[ , k], nReplicates)
            }
            names(blocks$data) <- names(x)
        }
    }

    if(is.null(logScale)) use.phi <- !check.constant(scaleFun) else use.phi <- logScale
    fit <- try(fevd(.y~1, data = xObs, threshold = thresholdObs, location.fun = locationFun, scale.fun = scaleFun, shape.fun = shapeFun, use.phi = use.phi, type = "PP", method = "MLE", weights = weightsObs, initial = initial, blocks = blocks, optim.args = optimArgs))

    failed <- TRUE
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
        msg <- fit
        attributes(msg) <- NULL
        fit <- list()
        fit$info <- list()
        fit$info$convergence <- 1
        fit$info$message <- paste("optimization call returned error: ", msg)
    } 
    fit$info$failure <- failed
    fit$noInt <- noInt
    return(fit)
} # end of fit_model_pot()

