#' Load RData file.
#' 
#' Load a single \code{R} object stored in an RData file and give the object a
#' name of choice.
#' 
#' @param file An RData file saved via \code{\link[base]{save}}.
#' @return The object save in \code{file}.
#' @references \url{http://stackoverflow.com/a/5577647}
#' @export
load_rdata <- function(file = NULL) {
    env <- new.env()
    nm <- load(file, env)[1]
    env[[nm]]
}

#' Convert MCMC List to Array
#' 
#' Converts an mcmc.list (each list entry: rows=MCMC samples; cols=parameters) to a 3 dimensional array (MCMC samples; chains; parameters)
#' Property of ??? (Internet!)
#' @param x mcmc.list (e.g., from JAGS) 
#' @export
mcmc.list2stan <- 
    function(x) {
        print(class(x))
        if(class(x) != "mcmc.list")
            warning(paste0("Object is not of the class mcmc.list"))
        arr <- array(NA_real_, 
                     dim = c(nrow(x[[1]]), length(x), ncol(x[[1]])))
        dimnames(arr)[[3]] <- colnames(x[[1]])
        for(i in seq_along(x)) arr[,i,] <- as.matrix(x[[i]])
        return(arr)
    }

#' Stan to coda object
#' 
#' Convert a 3 dimensional array (MCMC samples; chains; parameters) to an mcmc.list object for analysis with the coda package
#' @param fit fitted Stan object
#' @export
stan2mcmc.list <- function(fit) {
    #   if(class(fit) != "stan")
    #     warning("Object is not a stan object")
    coda::mcmc.list(lapply(1:ncol(fit), function(x) coda::mcmc(as.array(fit)[,x,])))
}
