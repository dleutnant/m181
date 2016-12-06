#' Generates noise for a section of a xts object using random sampling
#'
#' @title Introduces noise to a xts time series
#' @param TS [xts] the time series to be manipulated.
#' @param PercentageNoise [percent(integer)] set the percentage of the signal to be changed.
#' @param Noisespread [double] specifies the std (N~(0,x)) of the noise that will be implemented.
#' @param Overlap [logical] allows or denies the overlapping of noise sections. If the PercentageNoise
#' value is fairly high (>40) overlap should be set to TRUE. Otherwise the partition calculation within the function
#' will take very very long.
#' @param NoiseSections [integer] specifies how many noise sections to create.
#' @return xts. 
#' @rdname Noise_ISI
#' @import partitions
#' @examples 
#' 
#' # generate artifical dry weather flow
#' dat <- Dummy_xts_ISI(sampletime_min=1)
#' 
#' # inject noise
#' flawed_dat <- Noise_ISI(dat,10,1.3)
#' 
#' @export
Noise_ISI <- function(TS,PercentageNoise, Noisespread, Overlap=FALSE, NoiseSections=1){
  Times <- time(TS)
  Values <- as.numeric(TS)
  numnoise <- (length(Values)*PercentageNoise/100)
  #function speeds up the partition of the numnoise vector
  Factor <- function(num){
    count <- 0
    while (num > 999) {
      num <- num/10
      count <- count + 1}
    return(c(num,count))}
  
  NewNum.Factor <- Factor(numnoise)
  lenvec  <-  t(partitions::restrictedparts(NewNum.Factor[1], NoiseSections))[as.integer(runif(1,1,partitions::R(NoiseSections,NewNum.Factor[1]))),]
  lenvec <- lenvec*10^NewNum.Factor[2]
  indexes <- sample(1:length(Values), NoiseSections, replace = FALSE)
  if (Overlap) {while (any(indexes + lenvec - 1 > length(Values))) {
    indexes <- sample(1:length(Values), NoiseSections, replace = FALSE)}}
  else{
    while (any(diff(c(indexes,length(Values))) <= lenvec)) {
      indexes <- sample(1:length(Values), NoiseSections, replace = FALSE)}}
  for (i in 1:length(lenvec)) {
    start <- indexes[i]
    stop <- indexes[i] + lenvec[i] - 1
    Values[start:stop] <-  Values[start:stop] + rnorm(lenvec[i],0,Noisespread)}
  TS <- xts::xts(Values, order.by = Times)
  return(TS)
}

#' Generates Outliers in a xts object using random sampling
#'
#' @title Introduces outliers to a xts time series
#' @param TS [xts] the time series to be manipulated.
#' @param PercentageError [percent(integer)] set the percentage of the signal to be changed.
#' @param Outlierspread [double] specifies the std (N~(0,x)) of the outliers that will be implemented.
#' @param keyword [string] possible inputs are "long" and "short". These two keywords resemble the shape of outliers
#' where "long" stands for offsets and "short" for a traditional single outlier peak value.
#' @param maxlenoutlier [integer] specifies how many values offsets ("long" outlier) may include.
#' @param longspreadfactor [integer] specifies the ration of the outlier std and the std within an offset step.
#' @return xts. 
#' @rdname Outliers_ISI
#' @examples 
#' 
#' # generate artificical dry weather flow
#' dat <- Dummy_xts_ISI(sampletime_min=1)
#' 
#' # inject outliers
#' flawed_dat <- Outliers_ISI(dat,1,4)
#' 
#' @export
Outliers_ISI <- function(TS, PercentageError, Outlierspread, keyword="short", maxlenoutlier=10, longspreadfactor=0.04){
  len <- length(TS)
  indexTRFA <- logical(len)
  numerror <- (len*PercentageError/100)
  indexes <- floor(runif(n = floor(numerror), min = 2, max = (len - 1)))
  indexes <- indexes[!duplicated(indexes)]
  if (keyword == "short") {
    indexTRFA[indexes] <- !logical(length(indexes))
    noise <- rnorm(length(indexes),0,Outlierspread)
    TS[indexTRFA] <- TS[indexTRFA] + noise}
  else if (keyword == "long") {
    count <- 0
    i <- 1
    while (count < length(indexes)) {
      lenout <- floor(runif(1, min = 1, max = maxlenoutlier))
      alloffset <- rnorm(1,0,Outlierspread)
      if (indexes[i] + lenout > len) {lenout = len - indexes[i]}
      else{}
      singleoffset <- rnorm(lenout,0,(Outlierspread*longspreadfactor))
      indexTRFA <- vector(,len)
      indexTRFA[seq(indexes[i],(indexes[i] + lenout - 1),1)] = !vector(,lenout)
      noise <- (rep(alloffset,lenout) + singleoffset)
      TS[indexTRFA] <- TS[indexTRFA] + noise
      count <- count + lenout
      i <- i + 1}}
  else{warning("Wrong input String!")}
  return(TS)
}