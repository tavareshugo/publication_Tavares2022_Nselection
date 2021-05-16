#' Function to find intervals
#' 
#' returns a vector with an identifier for the peak
#'@param x a vector
#'@param threshold the threshold
#'@param below a logical. If TRUE returns "dips", i.e. intervals below the threshold
#'@param extend how many units to extend the interval to. The default (1), ensures the 
#' peak/dip is within the interval 
findPeaks <- function(x, threshold, below = FALSE, extend = 1){
  # Return logical of values above (or below) threshold
  if(below){
    pass_thresh <- x < threshold
  } else{
    pass_thresh <- x > threshold
  }
  
  # Find consecutive runs of FALSE and TRUE values
  edges <- rle(pass_thresh)
  
  # How many peaks/dips
  npeaks <- sum(edges$values)
  
  # Number each peak uniquely
  edges$values[which(edges$values)] <- 1:npeaks
  edges$values <- ifelse(edges$values, edges$values, NA)
  
  # Return a vector of peak IDs
  rep(edges$values, edges$lengths)
  
  
}
