#' Convert MCMC list to array.
#' 
#' Converts an mcmc.list (each list entry: rows=MCMC samples; cols=parameters) to a 3 dimensional array (MCMC samples; chains; parameters)
#' Property of ??? (Internet!)
#' @param x mcmc.list (e.g., from JAGS) 
#' @keywords internal
#' @export
mcmc.list2stan <- function(x) {
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
#' @keywords internal
#' @export
stan2mcmc.list <- function(fit) {
    #   if(class(fit) != "stan")
    #     warning("Object is not a stan object")
    coda::mcmc.list(lapply(1:ncol(fit), function(x) coda::mcmc(as.array(fit)[,x,])))
}

#' Load RData file.
#' 
#' Load a single \code{R} object stored in an RData file and give the object a
#' name of choice.
#' 
#' @param file An RData file saved via \code{\link[base]{save}}.
#' @inheritParams base::load
#' @return The object save in \code{file}.
#' @references \url{http://stackoverflow.com/a/5577647}
#' @keywords internal
#' @export
load_rda <- function(file = NULL, verbose = FALSE) {
    env <- new.env()
    nm <- load(file, env, verbose = verbose)[1]
    env[[nm]]
}

#' Get a Directory for Saving Files.
#' 
#' Check if dir exists and otherwise creates one.
#' 
#' @param dir Character.
#' @keywords internal
#' @return Character, namely, an existing directory.
# @export
get_dir <- function(dir = NULL) {
    # Returns either dir or---if unavailable---a path to a writable dir
    # Function is useful if a network dir becomes unavailable during a run of
    # recovery_irtree()
    if (!is.null(dir)) {
        try(dir <- suppressWarnings(normalizePath(dir)))
        if (!dir.exists(dir)) {
            dir <- NULL
        }
    } 
    if (is.null(dir)) {
        tmp1 <- tail(strsplit(tempdir(), "\\\\")[[1]], 1)
        tmp2 <- ifelse(dir.exists(getwd()), getwd(), normalizePath("~"))
        dirx <- suppressWarnings(normalizePath(paste0(tmp2, "/", tmp1)))
        dir.create(dirx, showWarnings = FALSE)
    } else {
        dirx <- dir
    }
    return(dirx)
}
