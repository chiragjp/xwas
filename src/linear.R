# Chirag Patel
# chirag_patel@hms.harvard.edu
## workhorse script to do linear associations
## 5/6/2015
## to run a sample: 
## > Rscript linear.R -f ../data/sample_data.Rdata -d LBXGLU -v LBXGTC -a ../data/adjustment_variables.Rdata
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
##

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

if(is.null(opt$outdir)) {
	outdir <- '.'
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

load(opt$filein) ### load in the file

keepVars <- c(depVar, varname, adjustmentVariables)
dat <- mainTab[complete.cases(mainTab[, keepVars]), keepVars]

depVar.formula <- depVar
baseform <- as.formula(sprintf('I(scale(%s)) ~ %s', depVar.formula, categorize_varname(varname, dat[, varname])))
nullmod <- as.formula(sprintf('I(scale(%s)) ~ %s', depVar.formula, categorize_varname(varname, dat[, varname])))
if(!categorical) {
	baseform <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))', depVar.formula, varname))
	nullmod <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))', depVar.formula, varname))
	if(!is.null(intVar)) {
		# intVar is assumed to be in the adjustmentVariables
		baseform <- as.formula(sprintf('I(scale(%s)) ~ I(scale(%s))*%s', depVar.formula, varname, intVar))
	}
}


doForm <- addToBase(baseform, adjustmentVariables)
print(doForm)

summaryFrame <- c()
if(permute == 0) {
	summaryFrame <- analyze_linear_mod(doForm, dat)
} else {
	## run permute analysis
	## do a parametric bootstrap
	nullmod.m <- linear_mod(nullmod, dat)
	responses <- resid(nullmod.m)
	for(ii in 1:permute) {
		dat[, depVar] <- sample(responses)
		frm <- analyze_linear_mod(doForm, dat)
		frm$permute_index <- ii
		summaryFrame <- rbind(summaryFrame, frm)
	}
}

if(!is.null(summaryFrame)) {
	summaryFrame$varname <- varname
}

## save result
filename <- sprintf('%s_%s.Rdata', depVar, varname)

fileout <- file.path(outdir, filename)
save(summaryFrame, varname, depVar, doForm, permute, file=fileout)
