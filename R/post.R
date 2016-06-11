#
# Chirag J Patel
# Nam Pho
#
# post.R - collection of functions providing post-processing and visualization of xwas results
#

#' manhattan
#'
#' Takes a vector of p-values and then displays in a manhattan plot with additional options. Usually can take the output of an xwas result
#' without further processing in most cases.
#'
#' @param data is the vector containing a list of p-values to display in a manhattan plot.
#' @param qthreshold is the cutoff of significance (the q-value).
#' @param nlabels is the number of top hits to label.
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
manhattan <- function(data, qthreshold=0.05, nlabels=6, ...) {
    data.log <- -log10(data)

    top <- sort(data.log, decreasing=TRUE)
    top <- top[is.finite(top)]

    i <- match(names(head(top, n=nlabels)), names(data.log))

    plot(data.log, ylab=expression(-log[10]*'(p-value)'), xlim=c(0, length(data)*1.1), pch=20, frame.plot=FALSE, ...)
    abline(h=-log10(qthreshold), col="red")
    points(i, data.log[i], col="red", pch=20)
    text(i, data.log[i], names(data.log[i]), cex=0.5, col="red", pos=1)
}

#' volcano
#'
#' Takes a list of linear regressions and displays a volcano plot.
#'
#' @param data is the vector containing a list of p-values to display in a manhattan plot.
#' @param qthreshold is the cutoff of significance (the q-value).
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
volcano <- function(data, qthreshold=0.05, ...) {
    data.all.p <- unlist(lapply(data, function(x) {return(x[2,]$"Pr(>|t|)")} ))
    data.all.x <- unlist(lapply(data, function(x) {return(x[2,]$Estimate)} ))
    #data.all.or <- unlist(lapply(data.all.p, function(x) {return((x/(1-x))^2)}))
    data.all.adjust <- p.adjust(data.all.p)

    all <- as.data.frame(cbind(p=data.all.p, x=data.all.x, q=data.all.adjust))

    sig <- subset(all, q < 0.05)

    plot(all$x, -log10(all$p), xlab="Effect", ylab=expression(-log[10]*'(p-value)'), pch=20, ...)
    points(sig$x, -log10(sig$p), col="red", pch=20)
}

