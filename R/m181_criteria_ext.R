#' Validates timeseries data by detecting outlier based on the MAD criteria
#' (Median absolute Deviation).
#' The output is an xts-object with logical values. Outliers are set to FALSE,
#' non-outliers to TRUE.
#'
#' @title Validation by outlier based on MAD criteria
#' @param xts An \code{xts} object to be validated.
#' @param width An integer specifying the window width.
#' @param multiplier The factor the MAD of all values is multiplied with.
#' @param constant Scale factor (if underlying distribution is normal, set to 1.4826!)
#' @param verbose logical. Provide additional details?
#' @keywords internal
#' @seealso \code{\link{xts}}.
check_mad_outlier <- function(xts, 
                              width, 
                              multiplier = 1, 
                              constant = 1,
                              verbose = TRUE){
  
  # constant auf 1 anstatt 1.4826, da keine normally distributed data vorliegen.
  # dann constant = 1, da = 1/Q(0.75) (0.75 quantile of that underlying distribution).
  # s. Leys et al. 2012
  
  # Im Unterschied zum Hampel Test wird hier der gleitende Median verwendet.
  # Hampel
  
  xts_with_flags <- xts
  
  # compute the anomaly criteria based on MAD  
  MAD <- xts::xts(x = .run_mad(x = zoo::coredata(xts),
                               n = width,
                               scale_factor = constant),
                  order.by = zoo::index(xts)) * multiplier
  
  # compute referencing mean
  rolling_median <- xts::xts(x = .run_median(x = zoo::coredata(xts),
                                             n = width),
                             order.by = zoo::index(xts))

  flag <- which(abs(xts - rolling_median) > MAD)
  
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
  attr(xts_with_flags, "m181") <- list("MAD" = list(
    parameter = list(width = width,
                     multiplier = multiplier,
                     constant = constant)))
  
  invisible(xts_with_flags)
  
}

#' Checks sign changes of an xts-object 
#'
#' @title Check sign changes of a timeseries
#' @param xts An \code{xts} object to be validated.
#' @param width An integer specifying the window width.
#' @param max_changes The maximum amount of allowed sign changes.
#' @param verbose logical. Provide additional details?
#' @keywords internal
#' @seealso \code{\link{xts}}, \code{\link{difftime}}, \code{\link{diff}}.
check_sign_change <- function(xts,
                              width,
                              max_changes, 
                              verbose = TRUE) {
  
  xts_with_flags <- xts
  
  # compute number of sign changes within rolling windows
  .s_tmp <- .count_sign_change(x = zoo::coredata(xts), n = width)
  ## add two NA's
  .s_tmp <- c(rep(NA, 2), .s_tmp)
  ## remove last two elements
  .s_tmp <- .s_tmp[1:nrow(xts)]
  
  sign_changes <- xts::xts(x = .s_tmp,
                           order.by = zoo::index(xts))
  
  flag <- which(sign_changes > max_changes)
  
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
  attr(xts_with_flags, "m181") <- list("sign_change" = list(
    parameter = list(width = width,
                     max_changes = max_changes)))
  
  invisible(xts_with_flags)
  
}
  
