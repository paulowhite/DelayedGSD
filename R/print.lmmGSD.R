## * print.lmmGSD
#' @export
print.lmmGSD <- function(x, ...){
    cat("       Analysis via the gls function \n \n ")
    cat("  Estimated treatment effect: ",x$name.coef," \n",sep = "")
    print(c("estimate" = x$estimate, "se" = x$se, "statistic" = x$statistic, "df" = x$df, "p.value" = x$p.value))
    if(!is.null(x$getInformation)){
        cat("\n  Number of clusters \n")
        print(x$n)
        cat("\n  Estimated information \n")
        print(x$getInformation)
        cat("\n",
            "total     : all patients including who dropped out early and have no observable outcome \n",
            "decision  : patients with at least one observable outcome, including those with not yet observed values \n",
            "interim   : patients with at least one observed outcome \n",
            "interim.cc: patients with at no missing outcome \n")
    }
    return(NULL)
}
