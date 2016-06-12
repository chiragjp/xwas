#
# Chirag J Patel
# Nam Pho
#
# post.R - collection of functions providing post-processing and visualization of xwas results
#

#' xwas.to.df
#'
#' Takes the output of the xwas function that returns a list of regressions and parses the results into a simple data.frame. This is useful for downstream plotting within the xwas library as well as custom implementations.
#'
#' @param data is a list of data.frame objects that represent the regression tests, usually this is the output from the xwas function.
#' @param ... are additional parameters that are passed on to the p.adjust function.
#'
#' @return a data.frame with the p-value, estimate effect size, and adjusted p-value.
#'
#' @examples
#' \dontrun{
#' xwas.to.df(data)
#' }
#'
#' @export
xwas.to.df <- function(data, ...) {
    if (!is.list(data)) warning("expecting data argument to be of type list.")

    data.p <- unlist(lapply(data, function(x) {return(x[2,]$"Pr(>|t|)")} ))

    data.x <- c()
    for (var in names(data.p)) {
    	data.x <- c(data.x, data[[var]][2,]$Estimate)
    
    }
    names(data.x) <- names(data.p)

    #data.or <- unlist(lapply(data.all.p, function(x) {return((x/(1-x))^2)}))
    data.q <- p.adjust(data.p, ...)

    all <- as.data.frame(cbind(p=data.p, x=data.x, q=data.q))

    return(all)
}

#' manhattan
#'
#' Takes a vector of p-values and then displays in a manhattan plot with additional options. Usually can take the output of an xwas result
#' without further processing in most cases.
#'
#' @param data is the vector containing a list of p-values to display in a manhattan plot.
#' @param qthreshold is the cutoff of significance (the q-value).
#' @param nlabels is the number of top hits to label.
#' @param ylab is a default y-axis label passed on to the plot function.
#' @param frame.plot is the option passed on to the plot function.
#' @param pch is default to 20 setting the point type to a small bullet in the plot function.
#' @param ... are additional arguments passed onto the plot function.
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' manhattan()
#' manhattan(glucose, main="XWAS of LBXGLU")
#' }
#'
#' @export
manhattan <- function(data, qthreshold=0.05, nlabels=6, ylab=expression(-log[10]*'(p-value)'), frame.plot=FALSE, pch=20, ...) {
    all <- NULL

    if (is.list(data)) {
        all <- xwas.to.df(data)
    } else {
        all <- data
    }

    data.log <- -log10(all$p)
    names(data.log) <- rownames(all)

    top <- sort(data.log, decreasing=TRUE)
    top <- top[is.finite(top)]

    i <- match(names(head(top, n=nlabels)), names(data.log))

    plot(data.log, ylab=ylab, frame.plot=frame.plot, pch=pch, ...)
    abline(h=-log10(qthreshold), col="red")
    points(i, data.log[i], col="red", pch=pch)
    text(i, data.log[i], names(data.log[i]), cex=0.5, col="red", pos=4)
}

#' volcano
#'
#' Takes a list of linear regressions and displays a volcano plot.
#'
#' @param data is the vector containing a list of p-values to display in a manhattan plot.
#' @param qthreshold is the cutoff of significance (the q-value).
#' @param xlab is a default x-axis label passed on to the plot function.
#' @param ylab is a default y-axis label passed on to the plot function.
#' @param pch is default to 20 which is a small bullet passed onto the plot function.
#' @param ... are additional arguments passed onto the plot function.
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' volcano()
#' volcano(data)
#' }
#'
#' @export
volcano <- function(data, qthreshold=0.05, xlab="Effect", ylab=expression(-log[10]*'(p-value)'), pch=20, ...) {
    all <- NULL
    
    if (is.list(data)) {
        all <- xwas.to.df(data)
    } else {
        all <- data
    }

    sig <- subset(all, q < 0.05)

    plot(all$x, -log10(all$p), xlab=xlab, ylab=ylab, pch=pch, ...)
    points(sig$x, -log10(sig$p), col="red", pch=pch)
}

#' qqplot
#'
#' Takes a data object that is either the output of the
#'
#' @param data is an object either a numeric character vector of p-values or the output of the xwas function from which that information can be derived.
#' @param xlab is a default label for the x-axis passed on to the plot function.
#' @param ylab is a default label for the y-axis passed on to the plot function.
#' @param pch is set to 20 by default which is a small bullet point passed on to the plot function.
#' @param ... are additional arguments passed on to the plot function.
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' qqplot()
#' qqplot(data)
#' }
#' 
#' @export
qqplot <- function(data, xlab=expression(-log[10]*'(p-expected)'), ylab=expression(-log[10]*'(p-empirical)'), pch=20, ...) {
    all <- NULL
    
    if (is.list(data)) {
        all <- xwas.to.df(data)
	
	pvalues <- all$p
	names(pvalues) <- rownames(all)
    } else if (is.numeric(data)) {
        pvalues <- data
    } else {
        stop("unexpected data type for argument data, should be a numeric character or xwas output.")
    }
    
    p.empirical <- sort(pvalues)
    p.expected  <- punif(1:length(pvalues)/length(pvalues))

    range <- c(min(p.empirical, p.expected), max(p.empirical, p.expected))
    
    plot(-log10(p.expected), -log10(p.empirical), xlab=xlab, ylab=ylab, xlim=range, ylim=range, pch=pch)
    abline(a=0, b=1, col="red")
}
