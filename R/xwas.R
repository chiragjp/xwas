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
       base_formula <- update.formula(base_formula, addStr)
    }
	
    return(base_formula)
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

#' analyze_linear_mod
#' 
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

#' analyze_logistic_mod
#'
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

#' xlm
#'
#' Performs a singular linear regression or binary outcomes association study. It requires you to specify the exact column variable 
#' and other permutation details. It does a lot of blind analysis and requires the xwas function to determine other meta-characteristics.
#' 
#' @param data is a data.frame containing the data to perform analysis on.
#' @param depvar a character vector with a column name for the dependent variable.
#' @param varname a character vector with a column name for the independent variable.
#' @param adjvars a list of variables to adjust for in the regression.
#' @param permute is an optional parameter of how many permutations to perform in the bootstrap.
#' @param categorical is a binary option representing whether or not the variable is categorical.
#' @param verbose is a boolean that determines if we print extra information. Not recommend for an XWAS, only one off analyses.
#' 
#' @return A data.frame object representing the regression.
#' 
#' @examples
#' \dontrun{
#' xlm()
#' xlm(data=nhanes, depvar="LBXGLU", varname="LBXGTC", adjvars=c("female", "RIDAGEYR"), permute=10)
#' }
#' 
#' @export
xlm <- function(data, depvar=NULL, varname=NULL, adjvars=c(), permute=0, categorical=0, verbose=TRUE) {
    if (permute < 0) {
        stop("non-negative value required for permute.")
    }
    
    intVar <- NULL
    
    keepVars <- c(depvar, varname, adjvars)
    dat <- data[complete.cases(data[, keepVars]), keepVars]

    if (nrow(dat) < 1) {
        if(verbose) {
            warning("insufficient number of cases.")
	}
        return(NULL)
    }
    
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
    if (verbose) {
       print(doForm)
    }
    
    summaryFrame <- c()
    summaryFrameRaw <- list()
    
    if (permute) {
	## run permute analysis
	## do a parametric bootstrap
	nullmod.m <- linear_mod(nullmod, dat)
	responses <- resid(nullmod.m)

	for (i in 1:permute) {
	    dat[, depvar] <- sample(responses)

	    # legacy character vector implementation
	    #frm <- analyze_linear_mod(doForm, dat)
	    #frm$permute_index <- i
	    #summaryFrame <- rbind(summaryFrame, frm)

	    # going with lists!
	    summaryFrameRaw[[i]] <- analyze_linear_mod(doForm, dat)
	}
    } else {
        summaryFrame <- analyze_linear_mod(doForm, dat)
    }

    # process permutation analysis into a singular data.frame
    # returns mean values for each permutation
    if (permute > 0 & length(summaryFrameRaw) > 0) {
       summaryFrame <- summaryFrameRaw[[1]]

       # get the mean value for each of the first 4 columns of the bootstrap
       # Estimate   Std. Error    t value     Pr(>|t|)
       for (i in 1:4) {
           tmp <- sapply(summaryFrameRaw, function(x) { return( x[, i] ) })
	   summaryFrame[, i] <- apply(tmp, 1, mean)
       }
    }

    if (!is.null(summaryFrame)) {
        summaryFrame$varname <- varname
        summaryFrame$permute <- permute	
    }
    
    return(summaryFrame)
}

#' xwas
#' 
#' Used to systematically perform an xlm across all variables. This function contains a lot of the logic for a systematic analysis
#' including logic to determine variable characteristics. 
#'
#' @param data is the data.frame containing what to analyze.
#' @param depvar is the outcome of the study we are looking to analyze in the context of multiple factors.
#' @param adjvars is a vector of variables to adjust for, if not specified we will scan all variables without adjusting on a first pass.
#' @param permute is the number of times to bootstrap each variable analysis in the xwas.
#' @param n is the number of cores to use for the multi-core implementation, value must be > 1 or set to "MAX".
#' @param verbose is a boolean to print extra information.
#'
#' @return A data.frame representing the output.
#'
#' @examples
#' \dontrun{
#' xwas()
#' xwas(data=nhanes, depvar="LBXGLU", permute=10)
#' }
#'
#' @export
xwas <- function(data, depvar=NULL, adjvars=NULL, permute=0, n=1, verbose=TRUE) {
    # require a non-negative set of permutations
    if (permute < 0) {
        stop("non-negative value required for permute.")
    }

    # need an outcome to test for
    if (is.null(depvar) | nchar(depvar) < 1) {
        stop("need to specify a valid outcome variable.")
    }

    test <- list()
    space <- NULL
    parallel <- NULL

    # test for multicore validity at least 1 core or auto
    if (is.character(n)) {
        n <- toupper(n) # convert to uppercase for consistency
    }
    
    if (n < 1 & n != "MAX") {
        stop("issue with multicore or parallel xwas settings.")
    } else if (n > 1 | n == "MAX") {
        parallel <- TRUE
    } else {
        parallel <- FALSE
    }

    if (is.null(adjvars)) {
        space <- colnames(data) # all variable space to loop through
    } else {
        space <- adjvars # if we know what to adjust for then use that list instead of the entire space
    }

    # R formulas throw errors on invalid factors, e.g. numerics
    # test the entire name space for validity
    invalid <- c()
    invalid <- space[!is.na(suppressWarnings( as.numeric(space) ))]

    if (length(invalid)) {
        stop(paste("the following column names are invalid:", invalid, "!"))
    }
    
    space <- space[!space %in% depvar] # remove the outcome from the variable space to control and explore

    if (parallel) {
        if (verbose) print("Doing multi core analysis.")
        library(foreach)
	library(doParallel)

        # start cluster if multicore
	if (parallel) maxcores <- detectCores()
    
        if (is.numeric(n)) {
    	    if (n > maxcores) {
	        warning("requested more cores than available, using the maximum.")
	    	n <- maxcores - 1
	    }
    	}

	if (n == "MAX") n <- maxcores - 1 # maximum cores!!!

	if (verbose) print(paste("Initializing a", n, "core cluster."))
    	cl <- makeCluster(n) # create the cluster

    	registerDoParallel(cl)

    	#
    	# THIS IS THE MULTI CORE VERSION OF THE MAIN LOOP
    	#
    	tmp <- foreach (varname = space, .export=c("test", "xlm", "addToBase", "linear_mod", "analyze_linear_mod")) %dopar% {
    	    if (verbose & is.null(adjvars)) {
    	        print( paste("Unadjusted testing", varname, which(space == varname), "of", length(space)) ) # DEBUG
	    } else {
    	        print( paste("Adjusted testing", varname, which(space == varname), "of", length(space)) ) # DEBUG
	    }

	    if (is.null(adjvars)) { 
	        test[[varname]] <- xlm(data=data, depvar=depvar, varname=varname, permute=permute, verbose=FALSE)
	    } else {
	        re.adjvars <- space[!space %in% varname] # adjust for everything except our independent variable of interest (varname)
    	        test[[varname]] <- xlm(data=data, depvar=depvar, varname=varname, adjvars=re.adjvars, permute=permute, verbose=FALSE)
	    }
    	}

	if (verbose) print("Terminating cluster.")
    	stopCluster(cl) # return resources if parallel option used
    } else {
        if (verbose) print("Doing single core analysis.")
	
        #
        # THIS IS THE SINGLE CORE VERSION OF THE MAIN LOOP
        #
    	for (varname in space) {
    	    if (verbose & is.null(adjvars)) {
    	        print( paste("Unadjusted testing", varname, which(space == varname), "of", length(space)) ) # DEBUG
	    } else {
    	        print( paste("Adjusted testing", varname, which(space == varname), "of", length(space)) ) # DEBUG
	    }

	    if (is.null(adjvars)) { 
	        test[[varname]] <- xlm(data=data, depvar=depvar, varname=varname, permute=permute, verbose=FALSE)
	    } else {
	        re.adjvars <- space[!space %in% varname] # adjust for everything except our independent variable of interest (varname)
    	    	test[[varname]] <- xlm(data=data, depvar=depvar, varname=varname, adjvars=re.adjvars, permute=permute, verbose=FALSE)
	    }
        }
    }

    #
    # END LOOP DIVERGENCE
    #

    if (length(test) & is.null(adjvars)) {
        # adjust for p-value to figure out which depvars are significant
	q <- p.adjust(lapply(test, function(x) {return(x[2,]$"Pr(>|t|)")} ))

	adjvars <- subset(q, q < 0.05)

	# re-run xwas with a limited variable space
	return(adjvars)
        #test <- xwas(data=data, depvar=depvar, adjvars=names(adjvars), permute=permute, verbose=verbose)
    }

    return(test)
}