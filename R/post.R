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
#'
#' @return A plot object.
#'
#' @examples
#' \dontrun{
#' manhattan()
#' manhattan(glucose, main="XWAS of LBXGLU")
#' }
#'
#' @export
manhattan <- function(data, qthreshold=0.05, nlabels=6, ...) {
    data.log <- -log(data)

    top <- sort(data.log, decreasing=TRUE)
    top <- top[is.finite(top)]

    i <- match(names(head(top, n=nlabels)), names(data.log))

    plot(data.log, ylab="-log(p-value)", xlim=c(0, length(data)*1.25), frame.plot=FALSE, ...)
    abline(h=-log(qthreshold), col="red")
    text(i, data.log[i]+20, names(data.log[i]))
}
