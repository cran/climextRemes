parseParamInput <- function(input, names = NULL) {
    # converts various formats for covariate specification to formula
    if(is.null(input))
        input <- ~1
    if(is.character(input)) 
        input <- stats::as.formula(paste0("~", paste0(input, collapse = "+")),
                            env = globalenv())
    if(is.numeric(input)) {
        inputName <- deparse(substitute(input))
        if(is.null(names)) stop("parseParamInput: 'names' required when 'input' is numeric.")
        if(sum(abs(input %% 1)) > 1e-15 || min(input) < 1)
            stop("parseParamInput: expecting integer-valued indices in ", inputName, ".") 
        input <- stats::as.formula(paste0("~", paste0(names[input], collapse = "+")),
                            env = globalenv())
    }
    if(methods::is(input, "formula")) {
        if(length(grep("\\-\\s{0,1}1",as.character(input))))
            stop("parseParamInput: no-intercept model specified but fitting only allowed for models with intercepts for location, scale, and shape.")
        environment(input) <- globalenv()
        return(input)
    } else stop("parseParamInput: expecting 'input' to be a formula, column names, or indices of columns.")
}

compute_return_quantities <- function(fit, returnPeriod = NULL, returnValue = NULL, x = NULL, x2 = NULL, locationFun, scaleFun, shapeFun, getParams = FALSE, upper = TRUE, scaling = 1, normalizers) {
    # auxiliary function for calculating parameter estimates (on original scale) and return-related quantities requested by user by calling out to functions specific to each quantity
    
    results <- list()
    p <- fit$results$num.pars
    np <- sum(unlist(p))
    
    summ <- summary(fit, silent = TRUE)
    paramNames <- names(summ$par)
    
    results$mle_raw <-  summ$par
    
    numLocScaleParams <- p[['location']] + p[['scale']]

    if(getParams) {
        mle <- results$mle_raw 
        if(!upper)  # location parameters for lower tail/minima are the negative of those computed based on negative of data values
            mle[1:p[['location']]] <- -mle[1:p[['location']]]
        if(is.element('se.theta', names(summ)))
            se <- summ$se.theta
        if(is.element('cov.theta', names(summ)))
            covmat <- summ$cov.theta
        nllh <- summ$nllh

        
        # transform coefficients back to original covariate scale
        if(!is.null(normalizers)) {
            if(p[['location']] > 1) {
                normMat <- stats::model.matrix(stats::as.formula(paste0(deparse(locationFun), "-1")), normalizers)
                tmpInd <- 1:p[['location']]
                mle[1] <- mle[1] - sum(normMat[1, ] * mle[2:p[['location']]] / (normMat[3, ] - normMat[2, ]))
                se[1] <- sqrt(c(1, -normMat[1, ] / (normMat[3, ] - normMat[2, ])) %*% covmat[tmpInd, tmpInd] %*% c(1, -normMat[1, ] / (normMat[3, ] - normMat[2, ])))
                mle[2:p[['location']]] <- mle[2:p[['location']]] / (normMat[3, ] - normMat[2, ])
                se[2:p[['location']]] <- se[2:p[['location']]] / (normMat[3, ] - normMat[2, ])
            }
            if(p[['scale']] > 1) {
                normMat <- stats::model.matrix(stats::as.formula(paste0(deparse(scaleFun), "-1")), normalizers)
                tmpInd <- (1+p[['location']]):(p[['location']]+p[['scale']])
                mle[p[['location']]+1] <- mle[p[['location']]+1] - sum(normMat[1, ] * mle[(p[['location']]+2):(p[['location']]+p[['scale']])] / (normMat[3, ] - normMat[2, ]))
                if(!fit$const.scale)
                    mle[p[['location']]+1] <- mle[p[['location']]+1] - log(scaling)
                se[p[['location']]+1] <- sqrt(c(1, -normMat[1, ] / (normMat[3, ] - normMat[2, ])) %*% covmat[tmpInd, tmpInd] %*% c(1, -normMat[1, ] / (normMat[3, ] - normMat[2, ])))
                mle[(p[['location']]+2):(p[['location']]+p[['scale']])] <- mle[(p[['location']]+2):(p[['location']]+p[['scale']])] / (normMat[3, ] - normMat[2, ])
                se[(p[['location']]+2):(p[['location']]+p[['scale']])] <- se[(p[['location']]+2):(p[['location']]+p[['scale']])] / (normMat[3, ] - normMat[2, ])
            }
            if(p[['shape']] > 1) {
                normMat <- stats::model.matrix(stats::as.formula(paste0(deparse(shapeFun), "-1")), normalizers)
                tmpInd <- (1+p[['location']]+p[['scale']]):np
                mle[p[['location']]+p[['scale']]+1] <- mle[p[['location']]+p[['scale']]+1] - sum(normMat[1, ] * mle[(p[['location']]+p[['scale']]+2):np] / (normMat[3, ] - normMat[2, ]))
                se[p[['location']]+p[['scale']]+1] <- sqrt(c(1, -normMat[1, ] / (normMat[3, ] - normMat[2, ])) %*% covmat[tmpInd, tmpInd] %*% c(1, -normMat[1, ] / (normMat[3, ] - normMat[2, ])))
                mle[(p[['location']]+p[['scale']]+2):np] <- mle[(p[['location']]+p[['scale']]+2):np] / (normMat[3, ] - normMat[2, ])
                se[(p[['location']]+p[['scale']]+2):np] <- se[(p[['location']]+p[['scale']]+2):np] / (normMat[3, ] - normMat[2, ])
            }
        } # end if(!is.null(normalizers))
        
        # rescale parameters so on scale of original data
        mle[1:p[['location']]] <- mle[1:p[['location']]] / scaling
        se[1:p[['location']]] <- se[1:p[['location']]] / scaling
        if(fit$const.scale) {
            mle[(p[['location']]+1):numLocScaleParams] <- mle[(p[['location']]+1):numLocScaleParams] / scaling
            se[(p[['location']]+1):numLocScaleParams] <- se[(p[['location']]+1):numLocScaleParams] / scaling
        }

        names(mle) <- names(se) <- paramNames
        results$mle <- mle
        results$se_mle <- se
        results$nllh <- nllh
    } # end if(getParams)

    # need as cbind of design matrices for loc,shape,scale when call make.qcov in various calc_ functions below
    if(!is.null(x)) {
        covariateMatrix <- stats::model.matrix(locationFun, x)
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(scaleFun, x))
        covariateMatrix <- cbind(covariateMatrix, stats::model.matrix(shapeFun, x))
    } else covariateMatrix <- NULL
    if(!is.null(x2)) {
        covariateMatrix2 <- stats::model.matrix(locationFun, x2)
        covariateMatrix2 <- cbind(covariateMatrix2, stats::model.matrix(scaleFun, x2))
        covariateMatrix2 <- cbind(covariateMatrix2, stats::model.matrix(shapeFun, x2))
    } else covariateMatrix2 <- NULL
    
    if(!is.null(returnPeriod)) {
        tmp <- calc_returnValue_fevd(fit, returnPeriod = returnPeriod, covariates = covariateMatrix)
        results$returnValue = tmp$returnValue / scaling
        results$se_returnValue = tmp$se_returnValue/ scaling
    } # end calculation of returnValue
        

    if(!is.null(returnValue)) {
        tmp <- calc_logReturnProb_fevd(fit, returnValue, covariates = covariateMatrix)
        results$logReturnProb <- tmp$logReturnProb
        results$se_logReturnProb <- tmp$se_logReturnProb
        # return period on log scale is negative of return probability
        results$logReturnPeriod <- -tmp$logReturnProb
        results$se_logReturnPeriod <- tmp$se_logReturnProb
    } # end calculation of returnPeriod, returnProb

    if(!is.null(x2)) {
        # get difference of return values or of log return probabilities (i.e. log of risk ratio)
        if(!is.null(returnPeriod)) {
            tmp <- calc_returnValueDiff_fevd(fit, returnPeriod = returnPeriod, covariates1 = covariateMatrix, covariates2 = covariateMatrix2)
            results$returnValueDiff <- tmp$returnValueDiff / scaling
            results$se_returnValueDiff <- tmp$se_returnValueDiff / scaling
        }
        if(!is.null(returnValue)) {
            tmp <- calc_logReturnProbDiff_fevd(fit, returnValue = returnValue, covariates1 = covariateMatrix, covariates2 = covariateMatrix2) 
            results$logReturnProbDiff <- tmp$logReturnProbDiff
            results$se_logReturnProbDiff <- tmp$se_logReturnProbDiff

        }
    } # end calculation of diff of returnValue, diff of log return prob

    return(results)
  } # end of compute_return_quantities()




