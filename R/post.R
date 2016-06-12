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
#' Takes output of the xwas function to display as a manhattan plot.
#'
#' @param data is the vector containing a list of p-values to display in a manhattan plot.
#' @param qthreshold is the cutoff of significance (the q-value).
#' @param nlabels is the number of top hits to label.
#' @param ylab is a default y-axis label passed on to the plot function.
#' @param frame.plot is the option passed on to the plot function.
#' @param pch is default to 20 setting the point type to a small bullet in the plot function.
#' @param legend is default to FALSE but will put a legend in the upper right corner if TRUE.
#' @param group.by is a data.frame object with n rows and 2 columns. Ideally one of the columns is named "var" but some logic is built in to try and decipher the appropriate column by percentage of matches.
#' @param group.label is a character vector of labels for each of the groups.
#' @param group.col is a color palette of for each group. It's an optional argument and defaults to the rainbow palette. If provided, it should be the same length at group.label as it is a 1:1 match.
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
manhattan <- function(data, qthreshold=0.05, nlabels=6, ylab=expression(-log[10]*'(p-value)'), frame.plot=FALSE, pch=20,
	              legend=FALSE, group.by=NULL, group.label=NULL, group.col=NULL, ...) {
    all <- NULL
    data.merged <- NULL
    data.log <- NULL

    if (!is.list(data)) {
        # it's p-values, just plot and exit
        data.log <- -log10(data)
	names(data.log) <- names(data)
    } else {
        # this is xwas output
        all <- xwas.to.df(data)
    
        if (!is.null(group.by)) {
	    if (ncol(group.by) != 2) {
	        stop("expecting only 2 columns for group.by parameter.")
	    }
	    
            if ( table(rownames(all) %in% group.by[, 2])[["TRUE"]] > ceiling(nrow(all)/2) ) {
	        # if the 2nd column of group.by has at least a 50% match to variable names then we use this to match categories
	        data.merged <- merge(data.frame(all, var=rownames(all)), group.by, all.x=TRUE, by.x="var", by.y=colnames(group.by)[2])

		if (!is.null(group.col) & length(group.col) != length(table(group.by[,2]))) {
		    warning("non-equal lengths of group.by and group.col provided.")
		}
	    } else {
	        data.merged <- merge(data.frame(all, var=rownames(all)), group.by, all.x=TRUE, by.x="var", by.y=colnames(group.by)[1])

		if (!is.null(group.col) & length(group.col) != length(table(group.by[,1]))) {
		    warning("non-equal lengths of group.by and group.col provided.")
		}
	    }

	    colnames(data.merged) <- c("var", colnames(all), "category") # explicitly define columns
	    data.merged <- data.merged[order(data.merged$category), ] # order the columns for plotting

            data.log <- -log10(data.merged$p)
	    names(data.log) <- data.merged$var
        } else {
	    data.log <- -log10(all$p)
	    names(data.log) <- rownames(all)
	}
    }

    # meta-data for top hits
    top <- sort(data.log, decreasing=TRUE)
    top <- top[is.finite(top)]
    i <- match(names(head(top, n=nlabels)), names(data.log))

    # skeleton of the main plot
    plot(data.log, ylab=ylab, frame.plot=frame.plot, pch=pch, col="white", ...)

    if (!is.null(group.by)) {
        if (is.null(data.merged)) {
	    stop("currently no support for just providing p-values as a vector in the data argument and with a group.by color option.")
	} else {
            if (is.null(group.label)) {
                group.label <- names(table(data.merged$category))
            }
    
            if (is.null(group.col)) {
       	        # if no color set is provided default to rainbow
                group.col <- rainbow(length(group.label))
            }

	    # do all the heavy lifting of point plot with color coordination by groups
            for (j in 1:nrow(data.merged)) {
       	        points(j, -log10(data.merged$p[j]), col=group.col[match(data.merged$category[j], group.label)], pch=pch)
	    }
        }
    } else {
        #if (is.null(data.merged)) {
	# no merged data means you were given p-values in a vector or you derived it from the xwas object and later set it
	# to the same variable names, so just plot
	points(data.log, col="black", pch=pch)
	#}
	
        points(i, data.log[i], col="red", pch=pch) # makes the top points red if no groups provided
    }

    if (legend & !is.null(group.label) & !is.null(group.col)) {
        legend("topright", group.label, text.col=group.col, bty="n", cex=0.75)
    }

    cutoff <- -log10(qthreshold)
    lines(c(0, length(data.log)), c(cutoff, cutoff), col="red")
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

    range <- c(min(-log10(p.empirical), -log10(p.expected)),
    	       max(-log10(p.empirical), -log10(p.expected)))
    
    plot(-log10(p.expected), -log10(p.empirical), xlab=xlab, ylab=ylab, xlim=range, ylim=range, pch=pch, ...)
    lines(c(range[1], range[2]), c(range[1], range[2]), col="red")
}
