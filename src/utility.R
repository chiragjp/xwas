### Chirag J Patel
### utility functions

addToBase <- function(base_formula, adjustingVariables) {
		form <- base_formula
		if(length(adjustingVariables)) {
			addStr <- as.formula(sprintf('~ . + %s', paste(adjustingVariables, collapse='+')))
			form <- update.formula(base_formula, addStr)
		}
		return(form)
}


log_variable_name <- function(variable) {
	if(logvariable(variable)) {
		return(sprintf('log(%s+1e-10)', variable ))
	} else {
		return(variable)
	}
}

is_variable_categorical <- function(variableArray) {
	## figure out if a categorical
	catTab <- table(mainTab[, varname])
	if(length(catTab) <= 10 & length(catTab) > 2) {
		return(sprintf('factor(%s)', varname))
	}
}

categorize_varname <- function(varname, variableArray, refgroup=0) {
	otherlevs <- paste(setdiff(names(table(variableArray)),as.character(refgroup)), collapse=",")
	return(sprintf('as.factor(%s,levels=c(%i, %s))', varname, refgroup,otherlevs))
}

qvalue_perm <- function(pvals, randData,numIter=100) {
	## compute a qvalue based on permuted data
	qvalueLow = function(pval) {
		numer = sum(randData <= (pval)) / numIter
		#denom = sum(pvals <= pval) / length(pvals)
		denom <- sum(pvals <= pval)
		q = numer/denom
		if(q > 1) {
			q=1
		}
		if(q == Inf) {
			q=NA
		}
		q
	}
	qvals <- sapply(pvals, qvalueLow)
	rankedPvals = rank(pvals)
	qsort = sort(qvals)
	qv = c()
	for(ii in 1:length(rankedPvals)) {
		currRank = rankedPvals[ii]
		qv = c(qv,qsort[currRank])
	}
	return(data.frame(qvalue=qv, pvalue=pvals))
}

linear_mod <- function(formula, dat, ...) {
	mod <- NULL
	summaryFrame <- NULL
	mod <- tryCatch(lm(formula, dat, ...), error = function(e) {
		print(e)
		return(NULL);
	})
	return(mod)
}

logistic_mod <- function(formula, dat, ...) {
	mod <- NULL
	summaryFrame <- NULL
	mod <- tryCatch(glm(formula, dat,family=binomial(), ...), error = function(e) {
		print(e)
		return(NULL);
	})
	return(mod)
}

analyze_linear_mod <- function(formula, dat,...) {
	mod <- linear_mod(formula, dat,...)
	N <- nrow(dat)
	if(!is.null(mod)) {
		summaryFrame <- as.data.frame(coef(summary(mod)))
		summaryFrame$N <- N
	}
	
	return(summaryFrame)
}

analyze_logistic_mod <- function(formula, dat, ...) {
	mod <- logistic_mod(formula, dat,...)
	N <- nrow(dat)
	if(!is.null(mod)) {
		summaryFrame <- as.data.frame(coef(summary(mod)))
		summaryFrame$N <- N
	}
	return(summaryFrame)
}

