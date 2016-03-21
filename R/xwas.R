#
# Chirag J Patel
# Nam Pho
# 
# xwas.R - collection of functions providing the backend for xwas analysis
#

#' Extend the base formula with variables we are adjusting for.
#'
#' @param base_formula a formula object with the dependent variable as a function of the independent variable
#' @param adjustingVariables a character vector of column names from the data.frame to adjust for
#'
#' @return An updated formula object representing a base formula with adjustment variables.
#'
#' @export
addToBase <- function(base_formula, adjustingVariables) {
    if (length(adjustingVariables)) {
       addStr <- as.formula(sprintf('~ . + %s', paste(adjustingVariables, collapse='+')))
       form <- update.formula(base_formula, addStr)
    }
	
    return(form)
}

#' @export
log_variable_name <- function(variable) {
    if (logvariable(variable)) {
        return(sprintf('log(%s+1e-10)', variable ))
    } else {
	return(variable)
    }
}

#' Determines if column in the main data.frame is categorical.
#'
#' @param mainTab is a data.frame representing the data to be analyzed.
#' @param varname is the name of a column of the data.frame to be determined if it is a categorical variable or not based on the number of factors.
#'
#' @return A factorization of the column of the data.frame.
#'
#' @examples
#' 
#' @export
is_variable_categorical <- function(mainTab, varname) {
    catTab <- table(mainTab[, varname])

    if (length(catTab) <= 10 & length(catTab) > 2) {
        return(sprintf('factor(%s)', varname))
    }
}

#' @export
categorize_varname <- function(varname, variableArray, refgroup=0) {
    otherlevs <- paste(setdiff(names(table(variableArray)),as.character(refgroup)), collapse=",")
	
    return(sprintf('as.factor(%s,levels=c(%i, %s))', varname, refgroup, otherlevs))
}

#' @export
qvalue_perm <- function(pvals, randData, numIter=100) {
    ## compute a qvalue based on permuted data
    qvalueLow <- function(pval) {
		     numer <- sum(randData <= (pval)) / numIter
		     #denom <- sum(pvals <= pval) / length(pvals)
		     denom <- sum(pvals <= pval)
		     q <- numer/denom
		     if (q > 1) {
		         q <- 1
		     }

		     if (q == Inf) {
		         q <- NA
		     }

		     q
		 }
	
    qvals <- sapply(pvals, qvalueLow)
    rankedPvals <- rank(pvals)
    qsort <- sort(qvals)
    qv <- c()
	
    for (i in 1:length(rankedPvals)) {
    	currRank <- rankedPvals[i]
	qv <- c(qv, qsort[currRank])
    }

    return(data.frame(qvalue=qv, pvalue=pvals))
}

#' Main worker function to perform linear associations.
#'
#' @param formula an R formula class object.
#' @param dat the data.frame to perform analysis on.
#' 
#' @return A data.frame representing the linear association study.
#'
#' @examples
#' \dontrun{
#' linear_mod()
#' }
linear_mod <- function(formula, dat, ...) {
    mod <- NULL
    
    mod <- tryCatch(lm(formula, dat, ...), error = function(e) {
	       print(e)
	       return(NULL);
	   })
	
    return(mod)
}

#' Main worker function to perform binary outcome associations.
#'
#' @param formula an R formula class object.
#' @param dat the data.frame to perform analysis on.
#'
#' @return A data.frame representing the binary outcome study.
#'
#' @examples
#' \dontrun{
#' logistic_mod()
#' }
logistic_mod <- function(formula, dat, ...) {
    mod <- NULL
	
    mod <- tryCatch(glm(formula, dat, family=binomial(), ...), error = function(e) {
	       print(e)
	       return(NULL);
	   })
	
    return(mod)
}

#' Main worker function to perform linear associations.
#'
#' @param formula an R formula class object.
#' @param dat the data.frame to perform analysis on.
#'
#' @return A data.frame representing the linear association study.
#'
#' @examples
#' \dontrun{
#' analyze_linear_mod()
#' }
#'
#' @export
analyze_linear_mod <- function(formula, dat, ...) {
    summaryFrame <- NULL
    N <- nrow(dat)
    
    #mod <- linear_mod(formula, dat, ...)
    mod <- tryCatch(lm(formula, dat, ...), error = function(e) {
               print(e)
	       return(NULL);
	   })
    
    if(!is.null(mod)) {
        summaryFrame <- as.data.frame(coef(summary(mod)))
	summaryFrame$N <- N
    }
	
    return(summaryFrame)
}

#' Main worker function to perform binary outcome associations.
#'
#' @param formula an R formula class object.
#' @param dat the data.frame to perform analysis on.
#'
#' @return A data.frame representing the binary outcome study.
#'
#' @examples
#' \dontrun{
#' analyze_logistic_mod()
#' }
#' 
#' @export
analyze_logistic_mod <- function(formula, dat, ...) {
    summaryFrame <- NULL
    N <- nrow(dat)

    #mod <- logistic_mod(formula, dat, ...)
    mod <- tryCatch(glm(formula, dat, family=binomial(), ...), error = function(e) {
               print(e)
	       return(NULL);
	   })
    
    if (!is.null(mod)) {
	summaryFrame <- as.data.frame(coef(summary(mod)))
	summaryFrame$N <- N
    }
    	
    return(summaryFrame)
}

#' Performs X-Wide Association Analysis.
#' 
#' @param data is a data.frame containing the data to perform analysis on.
#' @param depvar a character vector with a column name for the dependent variable.
#' @param varname a character vector with a column name for the independent variable.
#' @param adjvars a list of variables to adjust for in the regression.
#' @param permute is an optional parameter of how many permutations to perform in the bootstrap.
#' @param categorical is a binary option representing whether or not the variable is categorical.
#' 
#' @return A data.frame object representing the regression.
#' 
#' @examples
#' \dontrun{
#' xwas()
#' }
#' 
#' @export
xwas <- function(data, depvar, varname, adjvars=c(), permute=0, categorical=0) {
    if (permute < 0) {
        stop("non-negative value required for permute.")
    }
    
    intVar <- NULL
    
    keepVars <- c(depvar, varname, adjvars)
    dat <- data[complete.cases(data[, keepVars]), keepVars]
    
    depVar.formula <- depvar
    baseform <- NULL
    nullmod  <- NULL
    
    if (categorical) {
        baseform <- as.formula(sprintf('I(scale(%s)) ~ %s', depVar.formula, categorize_varname(varname, dat[, varname])))
        nullmod  <- as.formula(sprintf('I(scale(%s)) ~ %s', depVar.formula, categorize_varname(varname, dat[, varname])))
    } else { # if (!categorical) {
        baseform <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))', depVar.formula, varname))
	nullmod  <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))', depVar.formula, varname))

	if (!is.null(intVar)) {
	    # intVar is assumed to be in the adjustmentVariables
	    baseform <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))*%s', depVar.formula, varname, intVar))
	}
    }
    
    doForm <- addToBase(baseform, adjvars)
    print(doForm)
    
    summaryFrame <- c()
    if (permute) {
	## run permute analysis
	## do a parametric bootstrap
	nullmod.m <- linear_mod(nullmod, dat)
	responses <- resid(nullmod.m)
	
	for (i in 1:permute) {
	    dat[, depvar] <- sample(responses)
	    frm <- analyze_linear_mod(doForm, dat)
	    frm$permute_index <- i
	    summaryFrame <- rbind(summaryFrame, frm)
	}
    } else {
        summaryFrame <- analyze_linear_mod(doForm, dat)
    }
    
    if (!is.null(summaryFrame)) {
        summaryFrame$varname <- varname
    }
    
    return(summaryFrame)
}
