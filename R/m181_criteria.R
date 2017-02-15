#' @useDynLib m181
#' @importFrom Rcpp sourceCpp
NULL

#' Checks gaps of an xts-object and results a data.frame with computed gaps.
#'
#' @title Check gaps of a timeseries
#' @param xts An \code{xts} object to be validated.
#' @param interval_length The expected time between two indices.
#' @param interval_unit The unit of time interval.
#' @param gap_unit The unit of the computed gap length.
#' @param return_xts logical. If TRUE, both time series for index and value gaps 
#' are returned additionally.
#' @param tz A time zone specification to be used for the conversion, 
#' if one is required. System-specific (see time zones), but "" is the current 
#' time zone, and "GMT" is UTC (Universal Time, Coordinated). 
#' Invalid values are most commonly treated as UTC, on some platforms with a warning.
#' @return A data.frame containing starts and ends of gaps. Currently, a gap can be 
#' of type "index", which indicates that no index is available or of type "value", 
#' which basically means a "NA" value at a given time stamp. Additionally, both time series
#' for index gaps and value gaps can be returned if required.
#' @rdname check_gaps
#' @export 
#' @seealso \code{\link{xts}}, \code{\link{difftime}}, \code{\link{diff}}.
check_gaps <- function(xts,
                       interval_length = 60,
                       interval_unit = c("secs", "mins", 
                                         "hours", "days", "weeks"), 
                       gap_unit = c("auto", "secs", "mins", 
                                     "hours", "days", "weeks"),
                       return_xts = FALSE,
                       tz = "GMT") {
  
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))
  
  #### TIME GAPS ####
  
  # create diff vector and change unit
  diff_v <- diff(zoo::index(xts))
  units(diff_v) <- match.arg(interval_unit)
  idx <- which( diff_v > interval_length )
  
  # computes the starting indices of time gaps
  gaps_start <- as.POSIXct(zoo::index(xts[idx]),
                           tz = tz,
                           origin = "1970-01-01")
  
  # computes the ending indices of time gaps
  gaps_end <- as.POSIXct(zoo::index(xts[idx + 1]),
                         tz = tz,
                         origin = "1970-01-01")
  
  # computes the gap length in the specified unit
  delta_t <- difftime(time1 = gaps_end,
                      time2 = gaps_start,
                      tz = tz,
                      units = match.arg(gap_unit))
  
  # no time gaps 
  if (length(delta_t) == 0) {
    
    time_gaps <- data.frame(NA,NA,NA,NA)

  } else  {
    
    # creating a data frame with computed time gaps details
    time_gaps <- data.frame(gaps_start, 
                            gaps_end, 
                            delta_t, 
                            type = "time_gap", 
                            stringsAsFactors = TRUE)
    
  }
  
  colnames(time_gaps) <- c("gaps_start", "gaps_end", "delta_t", "type")

  
  #### VALUE GAPS ####
  
  # remove last position to exclude NA in last observation
  idx2 <- which(apply(is.na(xts[1:nrow(xts) - 1]), 1, any))
  
  # get rows of xts object which have at least one NA value
  v_gaps_start <- as.POSIXct(zoo::index(xts[idx2]),
                             tz = tz,
                             origin = "1970-01-01")
  
  v_gaps_end <- as.POSIXct(zoo::index(xts[idx2 + 1]),
                           tz = tz,
                           origin = "1970-01-01")
  
  # computes the gap length in the specified unit
  v_delta_t <- difftime(time1 = v_gaps_end,
                        time2 = v_gaps_start,
                        tz = tz,
                        units = match.arg(gap_unit))
  
  # no value gaps 
  if (length(v_delta_t) == 0) {
    
    output <- time_gaps
    
  } else {
    # creating a data frame with computed time gaps details
    value_gaps <- data.frame(v_gaps_start, 
                             v_gaps_end, 
                             v_delta_t, 
                             type = "value_gap", 
                             stringsAsFactors = TRUE)
    
    # set column names
    colnames(value_gaps) <- c("gaps_start", "gaps_end", "delta_t", "type")
    
    # merge both data.frames
    output <- rbind(time_gaps, value_gaps)
        
  }

  # order
  output <- output[ order(output[,1]), ]
  
  if (return_xts) {

    # create xts object with index gaps
    time_gap_idx <- seq(from = as.POSIXct(utils::head(zoo::index(xts),1),
                                          tz = tz, origin = "1970-1-1"),
                        to = as.POSIXct(utils::tail(zoo::index(xts),1),
                                        tz = tz, origin = "1970-1-1"),
                        by = paste(interval_length, 
                                   match.arg(interval_unit),
                                   sep = " "))
    
    time_gap_xts <- xts::xts(x = rep(NA, length(time_gap_idx)), 
                             order.by = time_gap_idx, 
                             tzone = tz)
    
    # remove doubled index
    time_gap_xts <- time_gap_xts[which(!zoo::index(time_gap_xts) %in% zoo::index(xts))]
    
    # create xts object with index gaps
    # no value gaps
    if (length(idx2) == 0) {
      
      value_gap_xts <- NULL
      
    } else {
      
      value_gap_xts <- xts::xts(x = rep(NA, length(idx2)), 
                                order.by = zoo::index(xts[idx2]), 
                                tzone = tz)
      
    }
    
    output <- list(gap_table = output,
                   gap_xts = list(time_gaps = time_gap_xts, 
                                  value_gaps = value_gap_xts))
  } 
  
  return(output)
}

#' Validates timeseries data by computing the difference of data points within 
#' given window width. Consecutive data points can be defined as constant by
#' two approaches (s. parameter "method"). The output is an xts-object with 
#' logical values. TRUE indicates constancy, FALSE inconstancy.
#' 
#' @title Validation by constance
#' @param xts An \code{xts} object to be validated.
#' @param delta_max The treshold of constance. Computed differences
#' lower than delta_max are considered to be constant, otherwise inconstant.
#' @param width An integer specifying the window width.
#' @param method Either "all" or "minmax".
#' Determines the way differences within \code{width} are calculated.
#' If set to "all", the differences of all data points within \code{width} must 
#' be lower then \code{delta_max} fo fulfill constancy criteria. If set to "minmax", 
#' data points are considered to be constant if the absolute difference of
#' min(x) and max(x) within \code{width} is lower than \code{delta_max}.
#' @param lag An integer indicating which lag to use (condition: method is "all")
#' @param relative logical. Should differences computed absolute (|x_i - x_i-1|) 
#' or relative ((|x_i - x_i-1|)/x_i) (condition: method is "all")?
#' @param verbose logical. Provide additional details?
#' @rdname check_constancy
#' @export 
#' @references DWA (2011). Merkblatt DWA-M 181 - Messung von Wasserstand und 
#' Durchfluss in Entwaesserungssystemen. Hennef.
#' @seealso \code{\link{xts}}, \code{\link{diff}}, \code{\link[zoo]{rollapply}}.
check_constancy <- function(xts, 
                            delta_max, 
                            width,
                            method = "all",
                            lag = 1,
                            relative = FALSE,
                            verbose = TRUE) {
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))

  if (method == "minmax") {
    
    # create a vector of elements indicating a computed difference of max and min 
    # is lower than delta_max within given width
    # logical flags are aligned right: flag belongs to previous "width" elements! 
    # Thus, "width-1" NA's are returned at the beginning of vector
    xts_with_flags <- xts::xts(x = .run_constancy(x = zoo::coredata(xts),
                                                  n = width,
                                                  delta_max = delta_max , 
                                                  method = 1), 
                               order.by = zoo::index(xts))
    

  } else {
    
    if (method == "all") {
      
      if (relative) {
        # set NA to prevent INF
        xts[xts == 0] <- NA
        # relative: (xi - (xi-1)) / xi
        xts <- abs(diff(xts, lag = lag)/xts)
      } else {
        xts <- abs(diff(xts, lag = lag))
      }
      
      xts_with_flags <- xts::xts(x = .run_constancy(x = zoo::coredata(xts),
                                                    n = width,
                                                    delta_max = delta_max , 
                                                    method = 2), 
                                 order.by = zoo::index(xts))
                                   
    } else {
      
      stop("unknown method")
      
    }
      
  }
     
  # catch the case if the vector which(xts_with_flags) is empty (integer(0)). 
  # This means all data are inconstant.
  if (identical(which(xts_with_flags),integer(0))) {
    
    # assigning FALSE to all elements 
    zoo::coredata(xts_with_flags)[,] <- FALSE
  
  } else {
  
    # extend elements with TRUE values 
    extended.indices <- c(sapply(which(xts_with_flags), FUN = function(x) 
      seq(x - (width - 1), length.out = width)))
    xts_with_flags[extended.indices] <- TRUE
  
  }
  
  # renaming col
  colnames(xts_with_flags) <- "flag"
  # print summary if verbose    
  if (verbose) {
    print(summary(data.frame(values = zoo::coredata(xts),
                             flag = zoo::coredata(xts_with_flags))))
  }
  
  # returns the xts flag when function call is assigned to a variable
  attr(xts_with_flags, "m181") <- list("constancy" = list(parameter = list(delta_max = delta_max,
                                                                           width = width,
                                                                           method = method,
                                                                           lag = lag,
                                                                           relative = relative)))

  invisible(xts_with_flags)
}

#' Validates timeseries data by isolating a reasonable range. 
#' The output is an xts-object with logical values. Data points within the 
#' range are set to TRUE, otherwise FALSE.
#'
#' @title Validation by range
#' @param xts An \code{xts} object to be validated.
#' @param y_min The lower bound.
#' @param y_max The upper bound.
#' @param verbose logical. Provide additional details?
#' @rdname check_range
#' @export 
#' @references DWA (2011). Merkblatt DWA-M 181 - Messung von Wasserstand und 
#' Durchfluss in Entwaesserungssystemen. Hennef.
#' @seealso \code{\link{xts}}.
check_range  <- function(xts, 
                         y_min = min(xts, na.rm = T), 
                         y_max = max(xts, na.rm = T), 
                         verbose=TRUE) {
  
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))
  
  # assign a new xts to store the flags  
  xts_with_flags <- xts
  
  # check if function parameter are feasible
  if (y_min > y_max) {stop("bad input: y_min > y_max")}
  
  # create a vector of elements referencing to data points within range
  flag <- which( xts >= y_min & xts <= y_max)
  
  # switch to mode "logical" to be consistent   
  mode(xts_with_flags) <- "logical"
  
  # catch the case if the vector "flag" is empty (integer(0)). This means all 
  # data are out of range.
  if (identical(flag,integer(0))) {
    # assigning FALSE to all elements 
    zoo::coredata(xts_with_flags)[,] <- FALSE
  } else {
    # assigning TRUE to all flag elements
    xts_with_flags[flag] <- TRUE
    # assigning FALSE to all but flag elements
    xts_with_flags[-flag] <- FALSE
  }
  # renaming col
  colnames(xts_with_flags) <- "flag"
  # print summary if verbose
  if (verbose) {
    print(summary(data.frame(values = zoo::coredata(xts),
                             flag = zoo::coredata(xts_with_flags))))
  }
  # returns the xts flag when function call is assigned to a variable
  attr(xts_with_flags, "m181") <- list("range" = list(parameter = list(y_min = y_min,
                                                                       y_max = y_max)))
  invisible(xts_with_flags)
  
}

#' Validates timeseries data by detecting outliers.
#' The output is an xts-object with logical values. Outliers are set to FALSE,
#' non-outliers to TRUE.
#'
#' The distinction whether a data point is an outlier or not is performed by
#' defining an anomaly criteria. This reflects the standard deviation (sd) of 
#' data points multiplied by a factor (sigma_multiplier, default=3). The sd is 
#' calculated either for all data points or for moving windows. If the
#' absolute difference of a data point and the (rolling) mean of data points
#' is higher than the anomaly criteria, the value is considered to be an outlier.
#'
#' @title Validation by outlier
#' @param xts An \code{xts} object to be validated.
#' @param sigma_multiplier The factor the standard deviation of all values 
#' is mulitplied with.
#' @param moving_window logical. Controls the way the anomaly criteria is
#' calculated. If set to TRUE, both mean and standard deviation are computed on 
#' moving windows (i.e. "rolling mean" and "rolling sd"). In contrast, if set to FALSE, 
#' the overall mean and the overall standard deviation of the given time series 
#' is used. This consequently leads to a static anomaly criteria.
#' @param width An integer specifying the window width.
#' @param type Either "mean" or "median". Determines the computation of the
#' reference values.
#' @param verbose logical. Provide additional details?
#' @rdname check_outlier
#' @export 
#' @references DWA (2011). Merkblatt DWA-M 181 - Messung von Wasserstand und 
#' Durchfluss in Entwaesserungssystemen. Hennef.
#' @seealso \code{\link{xts}}, \code{\link[zoo]{rollapply}}.
check_outlier <- function(xts, 
                          sigma_multiplier = 3,
                          moving_window = TRUE,
                          width,
                          type = c("mean", "median"),
                          verbose=TRUE) {
  
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))
  
  # assign a new xts to store the flags  
  xts_with_flags <- xts

  # should the moving_average and the moving_sd be used?
  if (moving_window) {
    
    # computing the rolling stats
    ref <- switch(match.arg(type),
                  mean = .run_mean(x = zoo::coredata(xts), n = width),
                  median = .run_median(x = zoo::coredata(xts), n = width)
                  )
    
    # compute the anomaly criteria: standard deviation * sigma_multiplier!
    crit <- .run_sd(x = zoo::coredata(xts), n = width) * sigma_multiplier
    
    # create a vector of elements referencing to data points higher the criteria
    flag <- which(abs(xts - ref) > crit)

  } else {

    # compute the anomaly criteria
    crit <- stats::sd(xts, na.rm = TRUE) * sigma_multiplier
    
    # computing the overall mean
    m <- base::mean(xts, na.rm = TRUE)
    
    # create a vector of elements referencing to data points higher the criteria
    flag <- which(abs(xts - m) > crit)
    
  }
  
  # switch to mode "logical" to be consistent   
  mode(xts_with_flags) <- "logical"
  # catch the case if the vector "flag" is empty (integer(0)).
  # This means no outliers detected.
  if (identical(flag,integer(0))) {
    # assigning TRUE to all elements 
    zoo::coredata(xts_with_flags)[,] <- TRUE
  } else {
    # assigning FALSE to all flag elements 
    xts_with_flags[flag] <- FALSE
    # assigning TRUE to all but flag elements
    xts_with_flags[-flag] <- TRUE
  }
  
  # renaming col
  colnames(xts_with_flags) <- "flag"
  # print summary if verbose
  if (verbose) {
    print(summary(data.frame(values = zoo::coredata(xts),
                             flag = zoo::coredata(xts_with_flags))))
  }
  # returns the xts flag when function call is assigned to a variable
  attr(xts_with_flags, "m181") <- list("outlier" = list(
    parameter = list(sigma_multiplier = sigma_multiplier,
                     moving_window = moving_window,
                     width = width, 
                     type = type)))

  invisible(xts_with_flags)
}

#' Validates timeseries data by computing the gradients (absolute or relative)
#' of consecutive data points. The output is an xts-object with logical values.
#' TRUE indicates plausible gradients, FALSE represents gradients which exceed
#' the threshold delta_max.
#'
#' @title Validation by gradient
#' @param xts An \code{xts} object to be validated.
#' @param delta_max The treshold of gradient. Computed gradients lower than 
#' delta_max are considered tolerable, higher ones are doubtful.
#' @param relative logical. Should the gradients computed absolute 
#' (|x_i - x_i-1|) or relative ((|x_i - x_i-1|)/x_i)?
#' @param lag An integer indicating which lag to use.
#' @param verbose logical. Provide additional details?
#' @rdname check_gradient
#' @export 
#' @references DWA (2011). Merkblatt DWA-M 181 - Messung von Wasserstand und 
#' Durchfluss in Entwaesserungssystemen. Hennef.
#' @seealso \code{\link{xts}}, \code{\link{diff}}.
check_gradient <- function(xts,
                           delta_max, 
                           relative=FALSE, 
                           lag=1,
                           verbose=TRUE) {
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))
  
  # assign a new xts to store the flags  
  xts_with_flags <- xts
  # computing the consecutive gradients
  if (relative) {
    xts <- abs(diff(xts, lag = lag)/xts)
  } else {
    xts <- abs(diff(xts, lag = lag))
  }
  # create a vector of elements indicating a computed gradient lower than 
  # delta_max or +-Inf # is.infinite to handle "Inf" result when 
  # dividing by zero
  flag <- which(xts < delta_max | is.infinite(xts)) 
  
  # switch to mode "logical" to be consistent
  mode(xts_with_flags) <- "logical"  
  # catch the case if the vector "flag" is empty (integer(0)).
  # This means all gradients are out of range.
  if (identical(flag,integer(0))) {
    # assigning FALSE to all elements 
    zoo::coredata(xts_with_flags)[,] <- FALSE
  } else {
    # assigning TRUE to all flag elements   
    xts_with_flags[flag] <- TRUE
    # assigning FALSE to all but flag elements     
    xts_with_flags[-flag] <- FALSE
  }
  # renaming col
  colnames(xts_with_flags) <- "flag"
  # print summary if verbose  
  if (verbose) {
    colnames(xts) <- ifelse(relative, "relative_diff", "absolute_diff")
    print(summary(data.frame(zoo::coredata(xts), 
                             flag = zoo::coredata(xts_with_flags))))
    
  }
  # returns the xts flag when function call is assigned to a variable
  attr(xts_with_flags, "m181") <- list("gradient" = list(
    parameter = list(delta_max = delta_max,
                     relative = relative,
                     lag = lag)))
  
  invisible(xts_with_flags)
}

#' Validates timeseries data by computing surrounded noise.
#' The output is an xts-object with logical values. TRUE indicates 
#' plausible data points, FALSE represents data points which exceed
#' the threshold s_max.
#'
#' @title Validation by noise
#' @param xts An \code{xts} object to be validated.
#' @param s_max The threshold of noise. Computed values lower than s_max are 
#' considered tolerable, otherwise not. If moving_s_max is TRUE, this value is 
#' computed with a moving window.
#' @param width An integer specifying the window width.
#' @param moving_s_max logical. Should the threshold computed with a moving
#' window?
#' @param width_s_max Sets the window of which the threshold is computed.
#' @param verbose logical. Provide additional details?
#' @rdname check_noise
#' @export 
#' @references DWA (2011). Merkblatt DWA-M 181 - Messung von Wasserstand und 
#' Durchfluss in Entwaesserungssystemen. Hennef.
#' @seealso \code{\link{xts}}, \code{\link[zoo]{rollapply}}.
check_noise <- function(xts, 
                        s_max,
                        width,
                        moving_s_max = TRUE, 
                        width_s_max = 2*width, 
                        verbose = TRUE) {
  
  message("This method is still under active development!")
  
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))
 
  # checks if function call is correct... 
  if (!moving_s_max) {
    if (missing(s_max)) stop("if moving_s_max is FALSE, 
                             a static s_max is required.")
  }
  
  # assign a new xts to store the flags 
  xts_with_flags <- xts
  # computing the rolling standard deviation
  rolling_sd <- .run_sd(x = zoo::coredata(xts), n = width)
  
  # create a vector of elements indicating a value higher than s_max
  if (moving_s_max) {

    # 1. compute s_max on larger windows?!
    s_max <- .run_sd(x = zoo::coredata(xts), n = width_s_max)
    s_max <- .run_median(x = s_max, n = width_s_max)
    # 2. usage of coef of var? feedback welcome!
    
  } 
  
  flag <- which(rolling_sd > s_max)
  
  # switch to mode "logical" to be consistent
  mode(xts_with_flags) <- "logical"  
  # catch the case if the vector "flag" is empty (integer(0)). This means all
  # values are below the threshold.
  if (identical(flag,integer(0))) { 
    zoo::coredata(xts_with_flags)[,] <- TRUE
  } else {
    # assigning TRUE to all flag elements  
    xts_with_flags[flag] <- FALSE
    # assigning FALSE to all but flag elements  
    xts_with_flags[-flag] <- TRUE
  }
  # renaming col
  colnames(xts_with_flags) <- "flag"
  # print summary if verbose 
  if (verbose) {
    print(summary(data.frame(rolling_sd = zoo::coredata(rolling_sd),
                             flag = zoo::coredata(xts_with_flags))))
  }
  
  attr(xts_with_flags, "m181") <- list("noise" = list(
    parameter = list(s_max = s_max,
                     width = width,
                     moving_s_max = moving_s_max,
                     width_s_max = width_s_max)))
  
  # returns the xts flag when function call is assigned to a variable
  invisible(xts_with_flags)
}

#' Gets NA's of a timeseries as an xts-object 
#'
#' @title Get NA's of a timeseries
#' @param xts An \code{xts} object to be validated.
#' @rdname check_na
#' @export 
#' @seealso \code{\link{xts}}.
check_na <- function(xts) {
  # checks if xts is really an xts object and stops if not
  stopifnot(xts::is.xts(xts))
  # assign a new xts to store the flags 
  xts_with_flags <- xts
  # switch to mode "logical" to be consistent
  mode(xts_with_flags) <- "logical" 
  # create a vector of elements indicating NA
  zoo::coredata(xts_with_flags) <- !is.na(xts)
  
  # returns the xts flag when function call is assigned to a variable
  attr(xts_with_flags, "m181") <- list("na" = list())
  
  # returns the data frame when function call is assigned to a variable
  invisible(xts_with_flags)
}
