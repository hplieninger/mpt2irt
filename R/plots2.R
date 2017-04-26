if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Plot observed frequencies.
#' 
#' Plot histograms of observed response frequencies separately for each item.
#' This is especially helpful to judge whether simulated data are reasonable
#' (i.e., unimodal etc.).
#' 
#' @param points how many resposne categories in Likert scale
#' @param ... Additional arguments passed to \code{\link[graphics]{barplot}}.
#' @inheritParams fit_irtree
#' @examples 
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' plot_responses(dat$X, revItem = dat$revItem, traitItem = dat$traitItem)
#' @export
plot_responses <- function(X,
                           revItem = rep(0, ncol(X)),
                           traitItem = rep(1, ncol(X)),
                           points = 5, ...){
    checkmate::assert_matrix(X, mode = "integerish", any.missing = FALSE,
                             min.rows = 2, min.cols = 2)
    J <- ncol(X)
    checkmate::assert_integerish(X, lower = 1, upper = 5, any.missing = FALSE)
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 len = J)
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 len = J)
    checkmate::qassert(points, "X1[1,]")
    
    # if(missing(revItem))
    #     revItem <- rep(1, J)
    # if(missing(traitItem))
    #     traitItem <- rep(1, J)
    n1 <- floor(sqrt(J) )
    n2 <- ceiling(J/n1)
    
    mfrow <- par()$mfrow
    mar <- par()$mar
    par(mfrow=c(n1,n2), mar=c(4, 3, 3, 1))
    if (!hasArg(ylim))
        ylim <- c(0, max(prop.table(table(factor(as.matrix(X), levels=1:points)))))
    for(j in 1:J){
        barplot(prop.table(table(factor(X[,j], levels=1:points))),
                main = paste("Trait", traitItem[j], "(rev:",revItem[j],")"),
                col = revItem[j]+1,
                ylim = ylim, ...)
    }
    
    par(mfrow=mfrow, mar=mar)
    
    return()
}

#' Model predicted response distribution.
#' 
#' Calculate the predicted response distribution given the posterior median/mean of the item parameters and a matrix \code{theta} of person parameters.
#' 
#' @param fit_sum List. A summary of theta and beta parameters as returned from \code{\link{tidyup_irtree_fit}}.
#' @param theta Matrix. A matrix of person parameters for which the predictions should be made.
#' @inheritParams fit_irtree
#' @inheritParams tidyup_irtree_fit
#' @inheritParams generate_irtree_ext
#' @return Function returns a data frame in long format with predicted response probability for each person-item combination.
#' @examples 
#' \dontrun{
#' # generate data
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' 
#' # fit model
#' res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200, warmup = 200)
#' res2 <- summarize_irtree_fit(res1)
#' res3 <- tidyup_irtree_fit(res2, N = N, J = J, revItem = dat$revItem,
#'                           traitItem = dat$traitItem, fitModel = res1$fitModel)
#' 
#' # expected frequencies
#' boeck_predict(res3)
#' }
#' @export
boeck_predict <- function(fit_sum = NULL,
                          theta = matrix(0, N, S2),
                          betas = NULL,
                          beta_ARS_extreme = NULL,
                          measure = c("Median", "Mean"),
                          S = NULL,
                          N = 1,
                          J = NULL,
                          revItem = NULL,
                          traitItem = NULL,
                          fitModel = NULL){
    
    measure <- match.arg(measure)
    if (is.null(betas))
        betas <- fit_sum$beta[[measure]]
    if (is.null(fitModel))
        fitModel <- fit_sum$fitModel
    if (is.null(beta_ARS_extreme) & fitModel == "ext")
        beta_ARS_extreme <- fit_sum$beta_ARS_extreme[[measure]]
    if (is.null(J))
        J <- nrow(fit_sum$beta[[measure]])
    if (is.null(revItem))
        revItem <- fit_sum$revItem
    if (is.null(traitItem))
        traitItem <- fit_sum$traitItem
    
    
    if (is.null(S)) {
        if (fitModel == "ext") {
            S <- 4
        } else if (fitModel == "2012") {
            S <- 3
        }
    }
    S2 <- S - 1 + length(unique(traitItem))
    
    p <- array(NA, c(N, J, 5))
    m <- y <- e <- a <- matrix(NA, N, J)
    if (fitModel == "2012") {
        a <- matrix(1, N, J)
    }

    for(i in 1:N){
        for(j in 1:J){      
            m[i, j] <- pnorm(theta[i, 1] - betas[j, 1])
            e[i, j] <- pnorm(theta[i, 2] - betas[j, 2])
            if (fitModel != "2012") {
                a[i, j] <- pnorm(theta[i, 3] - betas[j, 3])
            }
            y[i, j] <- pnorm(theta[i, (S - 1) + traitItem[j]] - betas[j, S])
            
            # reversed items: only relevant for trait-dimension/response process
            p.trait <- ifelse(revItem[j] == 1, 1 - y[i, j], y[i, j])
            
            # type of ERS in case of ARS:
            # extr_ars <- switch(as.character(fitModel),
            #                    "ext" = e[i,j], 
            #                    "ext2"=extreme_ARS, 
            #                    "ext3"=pnorm(theta[i,S]+qnorm(extreme_ARS)),
            #                    "2012"=.5)
            extr_ars <- switch(as.character(fitModel),
                               "ext" = pnorm(theta[i, 2] - beta_ARS_extreme),
                               # "ext" = e[i,j], 
                               # "ext2" = beta_ARS_extreme, 
                               "ext3" = pnorm(theta[i, S] - beta_ARS_extreme),
                               "2012" = .5)
            
            # response probabilities: MPT model from -2, -1, 0, 1, 2 // 0,1,2,3,4
            p[i,j,1] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*e[i,j]
            p[i,j,2] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*(1-e[i,j])
            p[i,j,3] <- (1-a[i,j])*m[i,j]
            p[i,j,4] <- (1-a[i,j])*(1-m[i,j])*p.trait*(1-e[i,j]) +a[i,j]*(1-extr_ars)
            p[i,j,5] <- (1-a[i,j])*(1-m[i,j])*p.trait*e[i,j]     +a[i,j]*extr_ars
        }
    }
    # colnames(p) <- item_order
    px <- reshape2::melt(p, varnames = c("Person", "Item", "Categ"), value.name = "Pred")
    px <- px[order(px$Item, px$Categ), ]
    px$Categ <- factor(px$Categ)
    rownames(px) <- NULL
    if (N == 1) {
        # px <- subset(px, select = -Person)
        px <- dplyr::select(px, -dplyr::matches("Person"))
    }
    return(px)
}


#' Plot predicted and observed frequencies.
#' 
#' Plot histograms of observed response frequencies separately for each item (via \code{\link{plot_responses}}), and adds the predictions from \code{\link{boeck_predict}}.
#' 
#' @param col Color of the prediction line.
#' @param lwd The line width, a positive number, defaulting to 2. The
#'   interpretation is device-specific, and some devices do not implement line
#'   widths less than one. (See the help on the device for details of the
#'   interpretation.)
#' @param ... Additional arguments passed to \code{\link[graphics]{lines}}.
#' @inheritParams fit_irtree
#' @inheritParams boeck_predict
#' @inheritParams plot_responses
#' @inheritParams graphics::lines
#' @inheritParams graphics::plot.default
#' @examples 
#' \dontrun{
#' # generate data
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' 
#' # fit model
#' res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200, warmup = 200)
#' res2 <- summarize_irtree_fit(res1)
#' res3 <- tidyup_irtree_fit(res2, N = N, J = J, revItem = dat$revItem,
#'                           traitItem = dat$traitItem, fitModel = res1$fitModel)
#' 
#' # plot expected and observed frequencies
#' plot_expected(res3, X = dat$X)
#' }
#' @export
plot_expected <- function(fit_sum,
                          X = NULL,
                          revItem = NULL,
                          traitItem = NULL,
                          points = 5,
                          type = "b",
                          col = "cyan",
                          lwd = 2,
                          ylim = NULL,
                          measure = c("Median", "Mean"), ...){
    checkmate::assert_matrix(X, mode = "integerish", any.missing = FALSE,
                             min.rows = 2, min.cols = 2)
    J <- ncol(X)
    checkmate::assert_integerish(X, lower = 1, upper = 5, any.missing = FALSE)
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 len = J, null.ok = TRUE)
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 len = J, null.ok = TRUE)
    checkmate::qassert(points, "X1[1,]")
    
    tmp1 <- boeck_predict(fit_sum, measure = measure, revItem = revItem,
                          traitItem = traitItem)
    pred <- reshape2::dcast(tmp1, Categ ~ Item, value.var = "Pred")
    
    if(is.null(revItem))
        revItem <-fit_sum$revItem
    if(is.null(revItem))
        revItem <-fit_sum$traitItem
    n1 <- floor(sqrt(J) )
    n2 <- ceiling(J/n1)
    
    mfrow <- par()$mfrow
    mar <- par()$mar
    par(mfrow=c(n1,n2), mar=c(4, 3, 3, 1))
    if (is.null(ylim)){
        ymax1 <- max(apply(X, 2, function(x) max(prop.table(table(x)))))
        ymax2 <- max(pred[, -1])
        ylim <- c(0, max(ymax1, ymax2))
    }
    
    for(j in 1:J){
        barplot(prop.table(table(factor(X[,j], levels=1:points))),
                main = paste("Trait", traitItem[j], "(rev:",revItem[j],")"),
                col = revItem[j]+1,
                ylim = ylim)
        lines(x = unclass(pred$Categ) - .5 + .2*unclass(pred$Categ),
              y = pred[, 1+j], col = col, type = type, lwd = lwd, ...)
    }
    
    par(mfrow=mfrow, mar=mar)
    
    return()
}