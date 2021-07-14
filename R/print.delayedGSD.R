## * print.delayedGSD
#' @export
print.delayedGSD <- function(x, ...){
  
  print(x$call)
  
  cat("\n")
  
  bnds <- cbind(c(1:length(x$lk)),x$lk,x$uk,c(x$ck,x$uk[length(x$uk)]))
  colnames(bnds) <- c("Stage","Lower","Upper","Decision")
  
  cat("Boundaries \n")
  print(bnds)
  
  cat("\n")
  
  cat("Planned maximum information \n")
  print(x$Info.max)
  
  cat("\n")
  
  cat("Inflation factor \n")
  print(x$InflationFactor)
  
  
}
