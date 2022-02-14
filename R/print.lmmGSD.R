## * print.lmmGSD
#' @export
print.lmmGSD <- function(x, ...){
    cat("       Analysis via the gls function \n \n ")
    cat("  Estimated treatment effect: ",x$name.coef," \n",sep = "")
    print(c("estimate" = x$delta$estimate, "se" = x$delta$se, "statistic" = x$delta$statistic, "df" = x$delta$df, "p.value" = x$delta$p.value))
    if(!is.null(x$information)){
        cat("\n  Number of clusters \n")
        print(x$sample.size)
        cat("\n  Estimated information \n")
        print(x$information)
        cat("\n",
            "decision  : patients with at least one observable outcome, including those with not yet observed values \n",
            "total     : all patients including who dropped out early and have no observable outcome \n",
            "interim   : patients with at least one observed outcome (full information) \n",
            "interim.cc: patients with at no missing outcome (complete case) \n")
    }
    return(NULL)
}
