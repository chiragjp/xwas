# Chirag Patel
# chirag_patel@hms.harvard.edu
## workhorse script to do binary outcome associations
source('utility.R')
library(getopt)

spec <- matrix(c(
				'varname', 'v', 1, 'character',
				'filein', 'f', 1, 'character',
				'is_categorical', 'c', 1, 'numeric',
				'depvar', 'd', 1, 'character',
				'adjustment_file', 'a', 1, 'character',
				'outdir', 'o', 1, 'character',
				'permute', 'p', 2, 'numeric', 
				'intvarname', 'i', 2, 'character'
				),
				nrow=8, byrow=TRUE)
opt <- getopt(spec)
mainTab <- read.csv(dat$filein)
outdir <- opt$outdir
varname <- opt$varname
depVar <- opt$depvar; 
intVar <- NULL

if(is.null(opt$filein)) {
	cat('specify a input file with -f flag!')
	quit(status=1)
}

if(!is.null(opt$intvarname)) {
	intVar <- opt$intvarname
}

permute <- 0
if(!is.null(opt$permute)) {
 	permute <- opt$permute
}
categorical <- 0
if(!is.null(opt$is_categorical)) {
	categorical <- opt$is_categorical
}


if(!is.null(opt$adjustment_file)) {
	load(opt$adjustment_file)
}


keepVars <- c(depVar, varname, adjustmentVariables, surveyVariables)
dat <- mainTab[complete.cases(mainTab[, keepVars]), keepVars]



depVar.formula <- depVar
baseform <- as.formula(sprintf('%s ~ %s', depVar.formula, categorize_varname(varname, dat[, varname])))
nullmod <- as.formula(sprintf('%s ~ %s', depVar.formula, categorize_varname(varname, dat[, varname])))
if(!categorical) {
	baseform <- as.formula(sprintf('%s ~ I(scale(%s))', depVar.formula, varname))
	nullmod <- as.formula(sprintf('%s ~ I(scale(%s))', depVar.formula, varname))
	if(!is.null(intVar)) {
		# intVar is assumed to be in the adjustmentVariables
		baseform <- as.formula(sprintf('%s ~ I(scale(%s))*%s', depVar.formula, varname, intVar))
	}
}


doForm <- addToBase(baseform, adjustmentVariables)
print(doForm)

summaryFrame <- c()
if(permute == 0) {
	summaryFrame <- analyze_logistic_mod(doForm, dat)
} else {
	## run permute analysis
	## do a parametric bootstrap
	for(ii in 1:permute) {
		dat[, depVar] <- sample(responses)
		frm <- analyze_logistic_mod(doForm, dat)
		frm$permute_index <- ii
		summaryFrame <- rbind(summaryFrame, frm)
	}
}

if(!is.null(summaryFrame)) {
	summaryFrame$varname <- varname
	summaryFrame$series <- seriesName
}

## save result
filename <- sprintf('%s_%s.Rdata', depVar, varname)

	
fileout <- file.path(outdir, filename)
save(summaryFrame, varname, depVar, doForm, permute, file=fileout)




