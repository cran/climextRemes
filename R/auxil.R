#' Normalize a vector
#'
#' Normalize a vector by subtracting off central point and dividing by range
#' @param vec vector of values
#' @param shift optional central point (if not provided, uses the mean of \code{vec})
#' @param lower optional lower end point of range (if not provided uses min of \code{vec})
#' @param upper optional upper end point of range (if not provided uses max of \code{vec})
#'
#' @export
#'
normalize <- function(vec, shift = NULL, lower = NULL, upper = NULL){
    if(is.null(vec))  # for Python interface
        stop("normalize: argument 'vec' is missing, with no default")
    if(is.array(vec)) # this will be the case for Python interface
        vec <- as.vector(vec)
    if(is.null(shift))
        shift <- mean(vec)
    if(is.null(lower))
        lower <- min(vec)
    if(is.null(upper))
        upper <- max(vec)
    return((vec - shift) / (upper - lower))
}

#' Remove consecutive exceedances from a vector
#'
#' Remove runs, i.e., consecutive exceedances, from a vector of values and associated indices (days); for use in declustering 
#' @param y vector of values
#' @param index vector of indices, one per value, that indicate which elements of \code{y} are consecutive
#' @param upperTail logical indicating whether values of \code{y} are upper (right) tail values (TRUE) or lower (left) tail values (FALSE)
#'
#' @export
#'
remove_runs <- function(y, index, upperTail = TRUE){
    if(is.null(y) || is.null(index))  # for Python interface
        stop("remove_runs: argument 'y' or 'index' is missing, with no default")
    if(is.array(y)) # this will be the case for Python interface
        y <- as.vector(y)
    if(is.array(index)) # this will be the case for Python interface
        index <- as.vector(index)

    n <- length(y)
    if(n != length(index)) stop("remove_runs: length of 'y' and 'index' should be the same.")
    if(n != length(unique(index))) warning("remove_runs: 'index' values are not unique.")
    if(!identical(index, sort(index))) {
        ord <- order(index)
        index <- index[ord]
        y <- y[ord]
    } else ord <- seq_len(n)
    
    if(n > 1){
        checkBack <- c(TRUE, index[2:n] - index[1:(n-1)] != 1)
        checkFwd <- c(index[1:(n-1)] - index[2:n] == -1, FALSE)
        seqStarts <- which(checkBack & checkFwd)
        for(start in seqStarts){
            pos <- start
            while(pos < n && index[pos + 1] - index[start] == pos + 1 - start) pos <- pos + 1
            tmp <- y[start:pos]
            jitterAmount <- 1e-15 * abs(ifelse(upperTail, max(tmp), min(tmp)))
            jitter <- stats::rnorm(length(tmp), 0, jitterAmount)  # amounts to randomly choosing a max when values are equal
            if(upperTail) {
                tmp[tmp + jitter < max(tmp + jitter)] <- NA
            } else tmp[tmp + jitter > min(tmp + jitter)] <- NA                              
            y[start:pos] <- tmp
        }
    }
    y[ord] <- y
    return(y)
}

#' Remove multiple exceedances within non-overlapping blocks of fixed length
#'
#' Remove multiple exceedances within non-overlapping blocks of fixed lengths, for use in declustering
#' @param y vector of values
#' @param index vector of indices, one per value, that indicate which elements of \code{y} are consecutive
#' @param blockLength length of block within which to remove all but the most extreme value
#'
#' @export
#'
screen_within_block <-  function(y, index, blockLength = 10){
    if(is.null(y) || is.null(index))  # for Python interface
        stop("screen_within_block: argument 'y' or 'index' is missing, with no default")
    if(is.array(y)) # this will be the case for Python interface
        y <- as.vector(y)
    if(is.array(index)) # this will be the case for Python interface
        index <- as.vector(index)

    n <- length(y)
    if(n != length(index)) stop("screen_within_block: length of 'y' and 'index' should be the same.")
    if(n != length(unique(index))) warning("screen_within_block: 'index' values are not unique.")
    if(!identical(index, sort(index))) {
        ord <- order(index)
        index <- index[ord]
        y <- y[ord]
    } else ord <- seq_len(n)
    if(n > 1){
        i <- 1
        while(i < length(y)){
            start <- i
            while( (index[i+1]-1) %/% blockLength == (index[i]-1) %/% blockLength &&  i < length(y)){
                i <- i + 1
            }
            if(i > start){
                tmp <- y[start:i]
                jitter <- stats::rnorm(i - start + 1, 0, 1e-15*abs(max(tmp)))
                tmp[tmp + jitter < max(tmp + jitter)] <- NA
                y[start:i] <- tmp
            }
            i <- i + 1
        }
    }
    y[ord] <- y
    return(y)
}
