import numpy
import cdms2
from rpy2.robjects.packages import importr
import argparse
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
# from mpi4py import *

numpy2ri.activate()



def get(obj, nm):
    " helper function to extract from R lists by name "
    ind = [i for i in range(len(obj.names)) if obj.names[i] == nm]
    return obj[ind[0]]

extRemes = importr("extRemes")
climextRemes = importr("climextRemes")

#############################
# basic stationary GEV fit
#############################

# very basic example: note you can name the arguments as they are named in the R function and as long as args are named it doesn't matter what order they are in
# pretend 50 random numbers are 50 block maxima for GEV fitting
result = climextRemes.fit_gev(numpy.random.randn(50), returnPeriod = 20, returnValue = 2.5)
get(result, 'returnValue')
get(result, 'se_returnValue')

# get Fort Collins precip data from extRemes so we can replicate examples from examples.R
Fort = robjects.r('''
   data(Fort)
   ord <- order(Fort$year, Fort$month, Fort$day) 
   Fort <- Fort[ord, ]
''')
Fort_Prec = Fort[5]
Fort_year = Fort[4]


firstYr = min(Fort_year)
yrs = range(int(min(Fort_year)), int(max(Fort_year)+1))
nYrs = len(yrs)

threshold = 0.395

FortExc = robjects.r('''
   threshold = 0.395
   ind <- Fort$Prec > threshold
   FortExc <- Fort[ind, ]
''')
FortExc_Prec = FortExc[5]
FortExc_year = FortExc[4]
FortExc_index = FortExc[0]

# NOTE: you don't need to do all this pre-processing data manip in R; I'm just doing it in R here because the Fort data (Fort Collins precip) is available as an R dataset, and I'm more familiar with R functions for data processing

#############################
# basic stationary POT fit
#############################

result = climextRemes.fit_pot(FortExc_Prec, firstBlock = firstYr,
        nBlocks = nYrs, threshold = threshold, blockIndex = FortExc_year,
        getParams = True, returnPeriod = 20,
        returnValue = 3.5, bootSE = True)
# # extract results by name
get(result, 'mle')  # MLE
get(result, 'returnValue')  # return value
get(result, 'se_returnValue')  # return value SE (asymptotic)
get(result, 'se_returnValue_boot')  # return value SE (bootstrap)
get(result, 'logReturnProb')  # log of probability of exceeding 'returnValue'
get(result, 'se_logReturnProb')  # return value SE (asymptotic)
get(result, 'se_logReturnProb_boot')  # return value SE (bootstrap)
get(result, 'logReturnPeriod')  # return value
get(result, 'se_logReturnPeriod')  # return value SE (asymptotic)
get(result, 'se_logReturnPeriod_boot')  # return value SE (bootstrap)


#############################
# non-stationary POT fit
#############################

result = climextRemes.fit_pot(FortExc_Prec, x = numpy.array(yrs), firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = FortExc_year, locationFun = 1, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(Fort_year), max(Fort_year)]), bootSE = False) 
get(result, 'returnValue')  # return value
get(result, 'se_returnValue')  # return value SE (asymptotic)
get(result, 'logReturnProb')  # log of probability of exceeding 'returnValue'
get(result, 'se_logReturnProb')  # return value SE (asymptotic)

# suppose you had two covariates - here I'll use year and a random vector just to illustrate syntax
# make 'x' be a 2-column numpy array, each column a covariate
# 'xNew' also needs to have 2 columns, each row is a different set of covariate values
covByBlock = numpy.c_[numpy.array(yrs), numpy.random.rand(nYrs)]
# careful - if x is a matrix, it should be column-major on the R side - Chris to verify this is ok
# ERROR: this is not quite working yet - issue with data structure for 'locationFun' being sent from Python to R (it works in R and I need to consult with Hari)
result = climextRemes.fit_pot(FortExc_Prec, x = covByBlock, firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = FortExc_year, locationFun = numpy.array([1,2]), getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([[min(Fort_year), 0], [max(Fort_year), 0]]), bootSE = False) 

######################################################################################
# using declustering
######################################################################################

# take max of exceedances on contiguous days
result = climextRemes.fit_pot(FortExc_Prec, x = numpy.array(yrs), firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = FortExc_year, index = FortExc_index, locationFun = 1, declustering = "noruns", getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(Fort_year), max(Fort_year)]), bootSE = False) 

# consider sequential blocks of 5 days and only use the max of any exceedances within a block
result = climextRemes.fit_pot(FortExc_Prec, x = numpy.array(yrs), firstBlock = firstYr, nBlocks = nYrs, threshold = threshold, blockIndex = FortExc_year, index = FortExc_index, locationFun = 1, declustering = 5, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(Fort_year), max(Fort_year)]), bootSE = False) 

######################################################################################
# with replicates
######################################################################################

# for illustration, I'll just duplicate the data here
result = climextRemes.fit_pot(numpy.append(FortExc_Prec, FortExc_Prec), x = numpy.array(yrs), firstBlock = firstYr, nBlocks = nYrs, nReplicates = 2, threshold = threshold, blockIndex = numpy.append(FortExc_year, FortExc_year), locationFun = 1, getParams = True, returnPeriod = 20, returnValue = 3.5, xNew = numpy.array([min(Fort_year), max(Fort_year)]), bootSE = False) 


###################################################
# fitting seasonal summed data, one obs per year
###################################################

# again a bit more data processing (in R for Chris' convenience)
FortSummerSumExc = robjects.r('''
   FortSummer <- Fort[Fort$month %in% 6:8, ]
   FortSummerSum <- aggregate(Prec ~ year, data = FortSummer, sum)
   threshold <- quantile(FortSummerSum$Prec, 0.8)
   FortSummerSumExc <- FortSummerSum[FortSummerSum$Prec > threshold, ]  # restrict 
''')
FortSummerSumExc_Prec = FortSummerSumExc[1]
FortSummerSumExc_year = FortSummerSumExc[0]
threshold = robjects.globalenv['threshold']

result = climextRemes.fit_pot(FortSummerSumExc_Prec, 
         firstBlock = firstYr, nBlocks = nYrs, blockIndex = FortSummerSumExc_year,
         threshold = threshold, 
         getParams = True, returnPeriod = 20,
         returnValue = 10,  bootSE = True)
get(result, 'returnValue')  # return value
# etc., etc.....


#############################################################################################
# to see the R help page for fit_pot, which describes all of the many possible arguments
#############################################################################################

# from the UNIX command line:

Rscript -e "library(climextRemes); help(fit_pot)"
