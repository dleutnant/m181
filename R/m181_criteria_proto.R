
#' Validates timeseries data with the moving window method by computing the standard deviation
#' for each window. The Regression option enables the smoothening of any underlying pattern in the data.
#' If turned on polynom will be fitted for each moving window and its residuals are analyzed.
#' The output is an xts-object with logical values.
#' TRUE indicates noiseless values, FALSE represents sections that were found to posess too much noise.
#' Please use the Analyse option to first asses how to set the filter in an optimal way.
#'
#' @title Validation by standart deviation of data or regression residuals (Institute for Urban Hydrology TU Dresden)
#' @param XTS An \code{xts} object to be validated. dim(xts) must be n,1!
#' @param s_max [double] the treshold of standard deviation.
#' @param window_width [integer] specifying the window width of the moving window.
#' @param lag [integer] distance between the starting points of the moving window.
#' @param Regression [logical] if TRUE the std of the residuals will be used or the noise detection.
#' @param Analyse [logical] illustrates the distribution of the std's If TRUE the filter will only analyse not correct the data
#' (default = FALSE).
#' @param Polydegree [integer] sets the degree of the polynomial function used for regression.
#' @rdname Noise_Filter_ISI
#' @examples 
#' 
#' # creating xts object
#' dat <- Dummy_xts_ISI(sampletime_min=1)  
#' 
#' # introduces noise
#' flawed_dat <- Noise_ISI(dat,10,1.3,NoiseSections=3)  
#' 
#' # analyses xts object
#' Noise_Filter_ISI(flawed_dat,window_width=20,lag=10,Regression=TRUE,Polydegree=4,Analyse=TRUE)  
#' 
#' #filters noise
#' logical_xts <- Noise_Filter_ISI(flawed_dat,s_max=0.5,window_width=20,lag=10,Regression=TRUE,Polydegree=4,Analyse=FALSE)
#' 
#' # returns xts object without noise  
#' corrected_dat <- flawed_dat[logical_xts]  
#'
#' @export 
#' @seealso \code{\link{xts}}.
Noise_Filter_ISI <-  function(XTS,
                              s_max = 99999,
                              window_width = 20,
                              lag = 1,
                              Regression = FALSE,
                              Polydegree = 1,
                              Analyse = FALSE){
  
  #Seperating the 2 columns of the xts object
  TS <-  XTS
  Times <-  time(TS)
  Values <-  as.numeric(TS)
  
  #preparing a vector for later indexing of the calculated std's
  Index = seq(1,length(Values),1)
  
  #If Regression == TRUE --> linear model fits polynomial function of order x to data 
  #and std. of residuals will be evaluated
  if (Regression) {
    
    #linear regression model for polynomial fit
    PolyN <-  function(X){
      y <-  as.numeric(X)
      x <-  (as.numeric(time(X)) - as.numeric(time(X[1])))/60
      
      #sets up a string that describes the linear model in R
      if (Polydegree == 1) {string = "fit=lm(y ~ x)"}
      else{
        string = "fit =lm(y ~ x" 
        for (i in 2:Polydegree) {
          string <-  paste(string,paste("+ I(x^",toString(i),")",sep = ""), sep = "")}
        string <-  paste(string,")", sep = "")}
      
      #executes the pre-set string
      eval(parse(text = string))
      
      residuals <-  resid(fit)
      return(sd(residuals))
      }
    
    #Calculates the SDVEC (standart dev. vector) of the data resiuals
    SDVEC <- zoo::rollapply(TS, width = window_width, by = lag, FUN = PolyN, align = "left")
    SDVEC <- SDVEC[!is.na(SDVEC)]
    
  } else {
    
    #Calculates the SDVEC (standart dev. vector) of the data resiuals
    SDVEC <- zoo::rollapply(Values, width = window_width, by = lag, FUN = sd, align = "left")}
  
  #The analyse option enables the visualization of the gradients for the given data set.
  if (Analyse) {print(paste("Mean standard deviation: ", sprintf("%e",mean(SDVEC)),sep = ""));
    print(paste("Min standard deviation: ", sprintf("%e",min(SDVEC))));
    print(paste("Max standard deviation: ", sprintf("%e",max(SDVEC))));
    par(mfrow = c(1,2))
    hist(SDVEC, xlab = "Standard Deviation")
    plot(index(SDVEC),SDVEC,xlab = "Moving Window",ylab = "Standard Deviation",main = "Location of Noise")
    par(mfrow = c(1,1))
    }
  
  #If the analyse option ist turned off (FALSE) the standard dev. will be filtered.
  else {
    
    #Determines the indexes of each moving window
    LowIndexes <- zoo::rollapply(Index, width = window_width, by = lag, FUN = min, align = "left")
    HighIndexes <- LowIndexes + window_width - 1
    
    #finds noise in SDVEC
    SDIndex <- which(SDVEC > s_max)
    
    #Warns if the s_max value is set too high. (Error will occur)
    if (length(SDIndex) == 0) {warning("S_max is too high. No Data has been filtered!")}
    
    #Filters the data and return clean time series
    else {lo = !logical(length(Values))
    for (i in 1:length(SDIndex)) {
      del <- c(LowIndexes[SDIndex[i]]:HighIndexes[SDIndex[i]])
      lo[del] <- logical(length(del))}
    Times <- Times[lo]
    Values <- Values[lo]}
    TS <- xts::xts(Values,order.by = Times)}
  
  TSlogical <- xts::xts(!is.na(cbind(XTS,TS)[,2]),order.by = time(XTS))
  
  return(TSlogical)
  
}


#' Validates timeseries data by computing the gradients [(x_i - x_i-1)/(t_i - t_i-1)]
#' of consecutive data points. The output is an xts-object with logical values.
#' TRUE indicates plausible gradients, FALSE represents gradients which exceed
#' the thresholds. Please use the Analyse option to first asses how to set the filter
#' in an optimal way.
#'
#' @title Validation by gradient (Institute for Urban Hydrology TU Dresden)
#' @param XTS An \code{xts} object to be validated. dim(xts) must be n,1!
#' @param mgrad [double] the treshold of negative gradient. Computed gradients lower than 
#' delta_max are considered ouliers, higher ones are tolerable.
#' @param pgrad [double] the treshold of positive gradient. Computed gradients lower than 
#' delta_max are considered tolerable, higher ones are outliers.
#' @param acceptableError [percent] is the stop criteria for the gradient filter (default = 0.01). 
#' It's repeating the calclulations until the number of criteria is met (error=Number_of_Gradient_Error/All_values*100)
#' @param UnplausibleGapsize [integer] a number specifying the distance between 2 gradient errors that will be regarded unplausible.
#' As well as any other gapsize smaller than the set number (default = NULL -> this option is turned off).
#' @param Analyse [logical] illustrates the distribution of the gradients. If TRUE the filter will only analyse not correct the data
#' (default = FALSE).
#' @param GradZero [logical] delete gradients = 0 (default = FALSE).
#' @rdname Gradient_Filter_ISI
#' @examples 
#' 
#' # creating xts object
#' dat <- Dummy_xts_ISI(sampletime_min=1) 
#' 
#' # introduces outliers 
#' flawed_dat <- Outliers_ISI(dat,1,4)  
#' 
#' # analyses xts object
#' Gradient_Filter_ISI(flawed_dat, Analyse = TRUE)  
#' 
#' # filters outliers
#' logical_xts <- Gradient_Filter_ISI(flawed_dat,mgrad=-.025,pgrad=.025, Analyse = FALSE)  
#' 
#' # returns xts object without outliers
#' corrected_dat <- flawed_dat[logical_xts]  
#' 
#' @export 
#' @seealso \code{\link{xts}}.
Gradient_Filter_ISI <- function(XTS, 
                             mgrad=-9999999, 
                             pgrad=9999999,
                             acceptableError=0.01, 
                             UnplausibleGapsize=NULL, 
                             Analyse=FALSE, 
                             GradZero=FALSE){
  
  #Extends the FALSE sections of a logical vector if the TRUE sections inbetween are <= a given value
  GradGap <- function(LogiVec,UnplausibleGapsize = 2){
    index <- which(!LogiVec)
    distance <- diff(index) - 1
    logi <- which(distance <= UnplausibleGapsize & distance > 0)
    for (i in logi) {
      vec <- seq(index[i],index[i + 1])
      index <- union(index,vec)}
    LogiVec[index] <- logical(length(index))
    return(LogiVec)}
  
  #Seperating the 2 columns of the xts object
  TS <- XTS
  TStime <- as.numeric(time(TS))
  TSvalue <- as.numeric(TS)
  
  #Calculating the Gradients
  x2 <- TStime[-1]
  x1 <- TStime[-length(TStime)]
  y2 <- TSvalue[-1]
  y1 <- TSvalue[-length(TSvalue)]
  deltaX <- x2 - x1
  deltaY <- y2 - y1
  grad <- deltaY/deltaX
  
  #The analyse option enables the visualization of the gradients for the given data set.
  if (Analyse) {
    NEG <- grad[grad < 0]
    POS <- grad[grad > 0]
    cat(paste("Mean neg. grad.:",sprintf("%e",mean(NEG)),"\n","Min neg. grad.:",sprintf("%e",min(NEG)),"\n","Max neg. grad.:",sprintf("%e",max(NEG)),sep = " "))
    cat("\n")
    cat(paste("Mean pos. grad.:",sprintf("%e",mean(POS)),"\n","Min pos. grad.:",sprintf("%e",min(POS)),"\n","Max pos grad.:",sprintf("%e",max(POS)),sep = " "))
    par(mfrow = c(1,2))
    hist(grad, xlab = "Gradient")
    plot(zoo::index(grad), grad, xlab = "Index", ylab = "Gradient",main = "Location of Outlier")
    par(mfrow = c(1,1))
    }
  
  #If the analyse option ist turned off (FALSE) the gradients will be filtered.
  else{
    
    #first run of the gradient check
    index1 <- !grad > pgrad 
    index2 <- !grad < mgrad
    ind <- index1 & index2
    if (GradZero) {index0 <- !grad == 0;ind <- ind & index0}
    if (is.null(UnplausibleGapsize) == FALSE) {ind <- GradGap(ind,UnplausibleGapsize)}
    ind <- c(TRUE,ind)
    TS <- TS[ind]
    
    #the gradient check will be repeated until the stop criteria (acceptableError) is met
    while ((length(ind) - table(ind)["TRUE"])/length(ind) > acceptableError/100) {
      TStime <-  as.numeric(time(TS))
      TSvalue <- as.numeric(TS)
      x2 <- TStime[-1]
      x1 <- TStime[-length(TStime)]
      y2 <- TSvalue[-1]
      y1 <- TSvalue[-length(TSvalue)]
      deltaX <- x2 - x1
      deltaY <- y2 - y1
      grad <- deltaY/deltaX
      index1 <- !grad > pgrad 
      index2 <- !grad < mgrad
      ind <- index1 & index2
      if (GradZero) {index0 <- !grad == 0;ind <- ind & index0}
      if (is.null(UnplausibleGapsize) == FALSE) {ind <- GradGap(ind,UnplausibleGapsize)}
      ind <- c(TRUE,ind)
      TS <- TS[ind]}
  }
  TSlogical <- xts::xts(!is.na(cbind(XTS,TS)[,2]),order.by = time(XTS))
  return(TSlogical)
}