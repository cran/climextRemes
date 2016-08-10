library(climextRemes)

library(extRemes)
data(Fort)
firstYr <- min(Fort$year)
yrs <- min(Fort$year):max(Fort$year)
nYrs <- length(yrs)

threshold <- 0.395

ord <- order(Fort$year, Fort$month, Fort$day) 
Fort <- Fort[ord, ]

# in most examples we consider 'year' to be a covariate with a linear effect on the location parameter of the POT model

### POT examples

# using extRemes
fit <- fevd(Prec, Fort, threshold = threshold, type="PP", units="inches", verbose=TRUE)
# linear effect of time in location
fit <- fevd(Prec, Fort, location.fun = ~year, threshold = threshold, type="PP", units="inches", verbose=TRUE)


# using climextRemes, analysis of daily data with year blocks
ind <- Fort$Prec > threshold
FortExc <- Fort[ind, ]

# stationary
fit_pot(FortExc$Prec, firstBlock = firstYr, blockIndex = FortExc$year,
        nBlocks = nYrs, threshold = threshold, getParams = TRUE,
        returnPeriod = 20,
        returnValue = 3.5, bootSE = TRUE)

# linear in time
fit_pot(FortExc$Prec, x = yrs, locationFun = 1,
        firstBlock = firstYr, blockIndex = FortExc$year,
        nBlocks = nYrs, threshold = threshold, getParams = TRUE,
        returnPeriod = 20, xNew = range(Fort$year),
        returnValue = 3.5, bootSE = TRUE)

# alternative using formula to specify location
fit_pot(FortExc$Prec, x = data.frame(year = yrs), locationFun = ~year,
        firstBlock = firstYr, blockIndex = FortExc$year,
        nBlocks = nYrs, threshold = threshold, getParams = TRUE,
        returnPeriod = 20, xNew = data.frame(year = range(Fort$year)),
        returnValue = 3.5, bootSE = TRUE)


# declustering sequential exceedances
julianDay <- seq_len(nrow(Fort))[ind]

fit_pot(FortExc$Prec, x = yrs, locationFun = 1,
        firstBlock = firstYr, blockIndex = FortExc$year,
        index = julianDay, declustering = "noruns",
        nBlocks = nYrs, threshold = threshold, getParams = TRUE,
        returnPeriod = 20, xNew = range(Fort$year),
        returnValue = 3.5, bootSE = TRUE)

# with replicates; for illustration here, I'll just duplicate the Fort data

Fort2x <- rbind(FortExc, FortExc)

fit_pot(Fort2x$Prec, x = yrs, locationFun = 1,
        firstBlock = firstYr, blockIndex = rep(FortExc$year, 2),
        nReplicates = 2,
        nBlocks = nYrs, threshold = threshold, getParams = TRUE,
        returnPeriod = 20, xNew = range(Fort$year),
        returnValue = 3.5, bootSE = TRUE)

# using climextRemes, analysis of total summer precipitation

FortSummer <- Fort[Fort$month %in% 6:8, ]
FortSummerSum <- aggregate(Prec ~ year, data = FortSummer, sum)

threshold <- quantile(FortSummerSum$Prec, 0.8)
FortSummerSumExc <- FortSummerSum[FortSummerSum$Prec > threshold, ]  # restrict to only exceedances

# stationary
fit_pot(FortSummerSumExc$Prec, firstBlock = firstYr,
        blockIndex = FortSummerSumExc$year,
        nBlocks = nYrs, threshold = threshold, getParams = TRUE,
        returnPeriod = 20, returnValue = 10, bootSE = TRUE)

# nonstationary
fit_pot(FortSummerSumExc$Prec, x = yrs, locationFun = 1,
        firstBlock = firstYr, blockIndex = FortSummerSumExc$year, 
        nBlocks = nYrs, threshold = threshold, getParams = TRUE, xNew = range(Fort$year),
        returnPeriod = 20, returnValue = 10, bootSE = TRUE)

### GEV examples

FortMax <- aggregate(Prec ~ year, data = Fort, max)

# linear effect of time in location

# using extRemes
fit <- fevd(Prec, FortMax, location.fun = ~year, type="GEV", units="inches", verbose=TRUE)
# note that fevd does not robustly handle the fact that the covariate values are all very large;
# fit_gev (and fit_pot) do handle this robustly

# using climextRemes
fit_gev(FortMax$Prec, x = FortMax$year, 
        locationFun = 1, getParams = TRUE, returnPeriod = 20,
        returnValue = 3.5, xNew = range(FortMax$year), bootSE = TRUE)





