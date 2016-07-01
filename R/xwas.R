#
# Chirag J Patel
# Nam Pho
# 
# xwas.R - collection of functions providing the backend for xwas analysis
#

#' Check if the object is a survey design object from the 'survey' library.
#'
#' @param design is an object to test for
#'
#' @return A boolean if the design argument provided is a 'survey.design' object or not.
#'
#' @export
is.survey <- function(design) {
    if("survey.design" %in% class(design)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Given a survey::svydesign object, verify it is real and then process explicit warnings if verbose output is requested.
#'
#' @param design is the survey::svydesign object to verify usage for.
#' @param verbose is a boolean defaulting to TRUE that outputs the warnings.
#'
#' @return A boolean to secondarily verify a survey::svydesign object to is.survey output.
#'
#' export
verify.survey <- function(design, verbose=TRUE) {
    if (!is.null(design)) {
        if (is.survey(design)) {
            if (verbose) warning("using survey-based regression.")

	    return(TRUE)
        } else {
	    if (verbose) warning("provided survey design is invalid object, ignoring.")
        }
    } else {
        if (verbose) warning("no survey-design object provided.")
    }

    return(FALSE)
}

#' Extend the base formula with variables we are adjusting for.
#'
#' @param base_formula a formula object with the dependent variable as a function of the independent variable
#' @param adjustingVariables a character vector of column names from the data.frame to adjust for
#'
#' @return An updated formula object representing a base formula with adjustment variables.
#'
#' @export
addToBase <- function(base_formula, adjustingVariables) {
    for (var in adjustingVariables) {
       base_formula <- update.formula( base_formula, as.formula(sprintf('~ . + %s', var)) )
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

#' http://www.r-statistics.com/2013/05/log-transformations-for-skewed-and-wide-distributions-from-practical-data-science-with-r/
#' 
#' @export
unskew <- function(x) {
    ifelse(abs(x) <= 1, 0, sign(x)*log10(abs(x)))
}

#' Determines if column(s) in the vector or data.frame is/are categorical.
#'
#' @param data is a data.frame or vector representing the data to be analyzed.
#' @param varname is the name of a column of the data.frame to be determined if it is a categorical variable or not based on the number of factors. It can be evaluated singularly (if a string) or if it is a vector of strings an equal length vector will be returned. If data is a vector then this argument can be optional, if it is optional and a data.frame is provided for the input data then all columns will be evaluated.
#' @param lower is a numeric to set the lower bound, number of unique entries is GREATER THAN this value. Default is 1 (at least 2+ uniques).
#' @param upper is a numeric to set the upper bound, number of unique entries is LESS THAN this value. Default is 11 (less than 11 uniques).
#'
#' @return a vector indicating if the data vector or data.frame column(s) are likely categorical variables by the presence of more than 2 but less than 10 unique counts.
#'
#' @examples
#' \dontrun{
#' is.categorical()
#' }
#' 
#' @export
is.categorical <- function(data, varname=NULL, lower=1, upper=11) {
    catTab <- NULL
    result <- c()

    if (is.vector(data)) {
        catTab <- as.data.frame(data)
	varname <- colnames(catTab)
    } else {
        catTab <- data
	
	if (is.null(varname)) {
	    varname <- colnames(data)
	} else {
            #catTab <- data[, varname]
	}
    }

    for (name in varname) {
        count <- length(table(catTab[, name]))
       	if (count < upper & count > lower) {
	    result <- c(result, TRUE)
	} else {
	    result <- c(result, FALSE)
	}
    }

    return(result)
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

#' linear_mod
#' 
#' Main worker function to perform linear associations.
#'
#' @param formula an R formula class object.
#' @param data the data.frame to perform analysis on, can be optional if the design argument is used.
#' @param design a survey::svydesign object that has the experimental design.
#' @param verbose is if we print errors, TRUE by default.
#' @param ... are additional arguments passed on to the regression.
#'
#' @return A data.frame representing the linear association study.
#'
#' @examples
#' \dontrun{
#' linear_mod()
#' }
#'
#' @export
linear_mod <- function(formula, data, design=NULL, verbose=TRUE, ...) {
    summaryFrame <- NULL
    N <- nrow(data)

    if(is.null(design)) {
        mod <- tryCatch(lm(formula, data, ...), error = function(e) {
                   if (verbose) print(e)
		   
	           return(NULL);
	       })
    } else {
        mod <- tryCatch(survey::svyglm(formula, design=design, ...), error = function(e) {
                   if (verbose) print(e)
		   
	           return(NULL);
	       })        
    }
    
    if(!is.null(mod)) {
        summaryFrame <- as.data.frame(coef(summary(mod)))
	summaryFrame$N <- N
    }
	
    return(summaryFrame)
}

#' logistic_mod
#'
#' Main worker function to perform binary outcome associations.
#'
#' @param formula an R formula class object.
#' @param data the data.frame to perform analysis on, can be optional if the design argument is used.
#' @param design a survey::svydesign object that has the experimental design.
#' @param ... are additional arguments passed on to the regression.
#'
#' @return A data.frame representing the binary outcome study.
#'
#' @examples
#' \dontrun{
#' logistic_mod()
#' }
#' 
#' @export
logistic_mod <- function(formula, data, design=NULL, ...) {
    summaryFrame <- NULL
    N <- nrow(data)

    if(is.null(design)) {
        mod <- tryCatch(glm(formula, data, family=binomial(), ...), error = function(e) {
                   print(e)
	           return(NULL);
	       })
    } else {
        mod <- tryCatch(survey::svyglm(formula, design=design, family=binomial(), ...), error = function(e) {
                   print(e)
	           return(NULL);
	       })        
    }

    if (!is.null(mod)) {
	summaryFrame <- as.data.frame(coef(summary(mod)))
	summaryFrame$N <- N
    }
    	
    return(summaryFrame)
}

#' survival_mod
#'
#' Main worker function to perform survival analysis.
#'
#' @param formula an R formula class object.
#' @param data the data.frame to perform analysis on, can be option if design argument is used.
#' @param design a survey::svydesign object that has the experimental design.
#' @param ... are additional arguments passed on to the regression.
#'
#' @return A data.frame representing the survival outcome study.
#'
#' @examples
#' \dontrun{
#' survival_mod()
#' }
#'
#' @export
survival_mod <- function(formula, data, design=NULL, ...) {
    summaryFrame <- NULL
    N <- nrow(data)

    # e.g., coxph(Surv(PERMTH_EXM, MORTSTAT) ~ RIDAGEYR + female + LBXGLU, data=suvdat)
    if(is.null(design)) {
        mod <- tryCatch(survival::coxph(formula, data, ...), error = function(e) {
                   print(e)
	           return(NULL);
	       })
    } else {
        mod <- tryCatch(survey::svycoxph(formula, design=design, ...), error = function(e) {
                   print(e)
	           return(NULL);
	       })        
    }

    if (!is.null(mod)) {
        summaryFrame <- as.data.frame(coef(summary(mod)))
    	summaryFrame$N <- N
    }

    return(summaryFrame)
}    

#' xlm
#'
#' Performs a singular linear regression or binary outcomes association study. It requires you to specify the exact column variable 
#' and other permutation details. It does a lot of blind analysis and requires the xwas function to determine other information such as
#' sanity checking around input arguments.
#' 
#' @param data is a data.frame containing the data to perform analysis on.
#' @param depvar a character vector with a column name for the dependent variable.
#' @param timevar is a character vector with the column name for the time variable in survival analysis.
#' @param varname a character vector with a column name for the independent variable.
#' @param adjvars a list of variables to adjust for in the regression.
#' @param design is an optional survey::svydesign object for weighted analysis.
#' @param permute is an optional parameter of how many permutations to perform in the bootstrap.
#' @param categorical is a binary option representing whether or not the variable is categorical.
#' @param verbose is a boolean that determines if we print extra information. Not recommend for an XWAS, only one off analyses.
#' 
#' @return a data.frame object representing the regression.
#' 
#' @examples
#' \dontrun{
#' xlm()
#' xlm(data=nhanes, depvar="LBXGLU", varname="BMXBMI", adjvars=c("RIDAGEYR", "female"))
#' }
#' 
#' @export
xlm <- function(data, depvar, timevar=NULL, varname=NULL, adjvars=c(), design=NULL, permute=0, categorical=0, verbose=TRUE) {
    intVar <- NULL # fix in the future, currently this is just a new name for adjvars
    
    keepVars <- c(depvar, timevar, varname, adjvars) # reduce data.frame space to data we want to study
    dat <- data[complete.cases(data[, keepVars]), keepVars]

    # not enough data points
    if (nrow(dat) < 1) {
        if (verbose) warning(paste("insufficient number of cases for", depvar))
	
        return(NULL)
    }

    # not enough types of data points
    # regression won't converge, return NULL, i.e. skip
    if (length(table(dat[, varname])) <= 1) {
        if (verbose) warning(paste("skipping", depvar, "only has 1 or fewer variable categories so regression won't converge"))

        return(NULL)
    }
    
    depVar.formula <- depvar
    baseform <- NULL
    nullmod  <- NULL
    
    if (categorical) {
        baseform <- as.formula(sprintf('I(scale(%s)) ~ %s', depvar, categorize_varname(varname, dat[, varname])))
        nullmod  <- as.formula(sprintf('I(scale(%s)) ~ %s', depvar, categorize_varname(varname, dat[, varname])))
    } else if (!is.null(timevar)) {
        # timevar specified so we're doing survival analysis
        baseform <- as.formula(sprintf('Surv(%s, %s) ~ %s', timevar, depvar, varname))
	nullmod  <- NULL # define in the future :)
    } else {
        # quantitative outcome so linear regression
        baseform <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))', depvar, varname))
	nullmod  <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))', depvar, varname))

	if (!is.null(intVar)) {
	    # intVar is assumed to be in the adjustmentVariables
	    baseform <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))*%s', depVar.formula, varname, intVar))
	}
    }
    
    doForm <- addToBase(baseform, adjvars)
    if (verbose) print(doForm)
    
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
	    #frm <- linear_mod(doForm, dat)
	    #frm$permute_index <- i
	    #summaryFrame <- rbind(summaryFrame, frm)

	    # going with lists!
	    summaryFrameRaw[[i]] <- linear_mod(doForm, data=dat, design=design)
	}
    } else {
        if (!is.null(timevar)) {
	    # survival analysis   
            summaryFrame <- survival_mod(formula=doForm, data=dat, design=design)
	} else {
	    # linear analysis
            summaryFrame <- linear_mod(formula=doForm, data=dat, design=design)
	}
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
#' @param depvar is required and the outcome or dependent variable we are looking to analyze in the context of multiple factors.
#' @param timevar is optional and if provided will be used in conjunction with the depvar to create a survival::Survey object for survival analysis.
#' @param varname is optional and the independent variable, without being specified we will consider all variables. If it is provided there will only be 1 regression test performed.
#' @param adjvars is optional and is a vector of a set of variables to adjust for under every condition, if not specified we will scan all variables without adjusting on a first pass.
#' @param design a survey::svydesign object that has the experimental design.
#' @param permute is the number of times to bootstrap each variable analysis in the xwas.
#' @param n is the number of cores to use for the multi-core implementation, value must be > 1 or set to "MAX".
#' @param verbose is a boolean to print extra information.
#'
#' @return A data.frame representing the output.
#'
#' @examples
#' \dontrun{
#' xwas()
#' xwas(data=nhanes, depvar="LBXGLU")
#' xwas(data=nhanes, depvar="LBXGLU", adjvars=c("BMXBMI", "RIDAGEYR", "female"))
#' }
#'
#' @export
xwas <- function(data, depvar, timevar=NULL, varname=NULL, adjvars=c(), design=NULL, permute=0, n=1, verbose=TRUE) {
    # check for actual data to analyze
    if (is.null(data)) {
        stop("data variable not specified, need data to analyze.")
    }

    # require a non-negative set of permutations
    if (permute < 0) {
        stop("non-negative value required for permute.")
    }

    # need an outcome to test for
    if (is.null(depvar) | !is.character(depvar) | !depvar %in% colnames(data)) {
        stop("depvar must be defined, be a character class, and exist within the data object.")
    }

    # test if we are doing survival analysis by if a variable for time (timevar) is specified
    if (!is.null(timevar) & !timevar %in% colnames(data)) {
        stop("timevar can't be found within the data object.")
    }
    
    verify.survey(design, verbose) # parses warnings on survey::svydesign usage
    
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

    space <- colnames(data)

    # R formulas throw errors on invalid factors, e.g. numerics
    # test the entire name space for validity
    invalid <- c()
    invalid <- space[!is.na(suppressWarnings( as.numeric(space) ))]

    if (length(invalid)) {
        warning(paste("the following column names are invalid:", invalid, "! Dropping them from data.frame and continuing analysis."))

	space <- space[!space %in% invalid] # just drop invalid column names, e.g. numerics pretending to be strings
    }

    if (length(varname)) {
        # if we know the dependent variable we're only going to do 1 test
        space <- varname 
    } else {
        # if we don't know/specify a dependent variable we search the entire space
        space <- space[!space %in% depvar]  # remove the outcome from the variable space
        space <- space[!space %in% adjvars] # remove the adjustment variables from the variable space
    }

    if (parallel) {
        if (verbose) print("Multi-core analysis.")
        #library(foreach)
	#library(doParallel)

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
    	tmp <- foreach (varname = space, .export=c("test", "xlm", "addToBase", "linear_mod")) %dopar% {
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
        if (verbose) print("Single-core analysis.")
	
        #
        # THIS IS THE SINGLE CORE VERSION OF THE MAIN LOOP
        #

    	for (testvar in space) {
    	    if (verbose) {
	        if (is.null(adjvars)) {
    	            cat( paste("Unadjusted testing for", testvar, "var on", depvar, which(space == testvar), "of", length(space), "\n") ) # DEBUG
	    	} else {
    	            cat( paste("Adjusted testing for", testvar, "var on", depvar, which(space == testvar), "of", length(space), "\n") ) # DEBUG
		}
	    }

	    test[[testvar]] <- xlm(data=data, depvar=depvar, timevar=timevar, varname=testvar, adjvars=adjvars, design=design, permute=permute, verbose=verbose)
        }
    }

    #
    # END LOOP DIVERGENCE
    #

    if (length(test) & length(adjvars) == 0) {
        # if we have test results and there was no dependent variable specified
	# we will re-run the xwas study with additional 

        # adjust for p-value to figure out which depvars are significant
	q <- p.adjust(lapply(test, function(x) {return(x[2,]$"Pr(>|t|)")} ))

	adjvars <- subset(q, q < 0.05)

	# re-run xwas with a limited variable space
	#return(adjvars)
	return(test)
        #test <- xwas(data=data, depvar=depvar, adjvars=names(adjvars), permute=permute, verbose=verbose)
    }

    return(test)
}
