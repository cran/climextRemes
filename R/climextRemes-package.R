#' climextRemes
#'
#' Functions for fitting GEV and POT (via point process fitting) models for extremes in climate data, providing return values, return probabilities, and return periods for stationary and nonstationary models. Also provides differences in return values and differences in log return probabilities for contrasts of covariate values. Functions for estimating risk ratios for event attribution analyses, including uncertainty. Under the hood, many of the functions use functions from 'extRemes', including for fitting the statistical models. Details are given in Paciorek, Stone, and Wehner (2018) <doi:10.1016/j.wace.2018.01.002>
#' 
#' @name climextRemes
#' @docType package
#' @import extRemes 
NULL

#' @importFrom methods is
#' @importFrom stats as.formula dbinom optim optimize model.matrix qchisq qnorm rbinom rnorm sd
#' @importFrom boot boot.ci boot
NULL
