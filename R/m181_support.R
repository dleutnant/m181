#' Checks equidistancy of an xts-object.
#' 
#' A time series is considered to be completely equidistant if the differences of 
#' diff(index)-percentiles \code{p(1-percentage)}
#' and \code{p(percentage)} is 0.
#'
#' @title Check equidistancy of an xts-object
#' @param xts An \code{xts} object to be validated.
#' @param percentage The level of significance to be checked
#' @rdname check_equidistance
#' @export 
#' @seealso \code{\link{xts}}, \code{\link{difftime}}, \code{\link{diff}}.
check_equidistance <- function(xts, percentage = 0.99) {
  
  # checks if xts is really an xts object  and length of percentage!=1 and 
  # stops if not
  stopifnot(xts::is.xts(xts), length(percentage) == 1)
  
  # calculate differences of consecutive indices
  delta.t <- diff(zoo::index(xts))
  
  # calculate "interquantile" ranges (1-percentage to percentage) of delta.t
  delta.t.percentile <- as.numeric(stats::quantile(as.numeric(delta.t), 
                                                   probs = 1 - percentage)) - 
                        as.numeric(stats::quantile(as.numeric(delta.t), 
                                                   probs = percentage))
  
  # caculate modus of delta.t
  delta.t.mode <- as.numeric(names(table(delta.t))[which.max(table(delta.t))]) 
  
  # prepare output (list containing result, percentage,
  # mode and length of indices)
  
  # a time series is considered to be equidistant if the difference of
  # percentiles P1-X and pX is 0.
  if (delta.t.percentile != 0) {
    res <- FALSE
  } else {
    res <- TRUE
  }

  return(list(result = res, 
              percentage = percentage,
              mode = delta.t.mode,
              length = length(delta.t)))
}

#' Fills gaps of an xts-object and results a modified xts-object.
#'
#' @title Fill gaps of a timeseries
#' @param xts An \code{xts} object to be filled.
#' @param mode The mode to fill the gaps with: "0" for a specific value, 
#' e.g. NA or 0, "1" for locf, "2" for nocb, "3" for linear interpolation,
#'  "4" for spline approximation.
#' @param fill_with The value to fill the gap (requires mode = 0)
#' @param fill_index logical. Should the index be filled so that the timeseries 
#' becomes equidistant? Corresponding values for filled indices become NA.
#' @param by The interval the index should be aligned to
#' @param na_rm logical. Should leading NAs be 
#' removed (affects mode=1 and mode=2)?
#' @param maxgap maximum number of consecutive NAs to fill. Any longer gaps
#'  will be left unchanged.
#' @param verbose logical. Should informative outputs printed during function
#' evaluation?
#' @rdname fill_gaps
#' @export 
#' @seealso \code{\link{xts}}, \code{\link[zoo]{na.locf}}, 
#' \code{\link[zoo]{na.fill}}.
fill_gaps <- function(xts, mode=5, fill_with=NA, 
                      fill_index=TRUE, by="mins",
                      na_rm=TRUE, maxgap=10, verbose=FALSE){
  # Checks if xts is really a xts object and stops if not
  stopifnot(xts::is.xts(xts))
  
  # Should the index be filled
  if (fill_index) {
    # compute range of given xts
    rng <- base::range(zoo::index(xts))
    # generate equidistant index
    index <-  seq(rng[1], rng[2], by = by)
    # create temporal xts object to be mergeable as index reference
    tmp.xts <- xts::xts(order.by = index,
                        tzone = "GMT")
    # merge given xts with the "index-only" tmp.xts 
    tmp <- merge(xts, tmp.xts, all = T)
    # print the amount of filled gaps if verbose
    if (verbose) message(paste("filled",
                               length(zoo::index(tmp)) - length(zoo::index(xts)),
                               "gaps"))
    xts <- tmp

  }
  
  #0. value, e.g. NA (default) or 0
  if (mode == 0) {
    xts[is.na(xts)]  <- fill_with
  }
  #1. last observation carried forward (locf)
  if (mode == 1) {
    xts <- zoo::na.locf(xts, na.rm = na_rm, fromLast = FALSE)
  }
  #2. next observation carried backward (nocb)
  if (mode == 2) {
    xts <- zoo::na.locf(xts, na.rm = na_rm, fromLast = TRUE)
  }
  #3. linear interpolation
  if (mode == 3) {
    #xts <- zoo::na.fill(xts, c("extend"))
    xts <- zoo::na.approx(xts, maxgap = maxgap)
  }
  #4. spline interpolation
  if (mode == 4) {
    xts <- zoo::na.spline(xts, maxgap = maxgap)
  }

  # returns the modified xts when function call is assigned to a variable
  invisible(xts)
}


#' Checks the presence of NA indices of an xts-object.
#'
#' @title Check time series index on NA elements
#' @param xts An \code{xts} object to be validated.
#' @return logical. If index contains NA values, TRUE is returned. 
#' @rdname index_has_na
#' @export
index_has_na <- function(xts) {
  
  # checks if xts is really an xts object
  stopifnot(xts::is.xts(xts))
  
  # check if index contains NA's
  if (any(is.na(zoo::index(xts)))) {
    
    warning(paste(length(which(is.na(zoo::index(xts)))),
                  "NA's found in index."))
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
    
  }
  
}