## Chirag Patel
## prepare a test dataset in .Rdata form
##
load('~/projects/nhanes_quant_pheno/bigTable_quant_pheno_chron_telo.Rdata')
mainTab <- subset(mainTab, SDDSRVYR == 3)
varDesc <- subset(varDesc, series == '2003-2004')

cols <- c(varDesc$var, as.character(demographicVariables$var))
sum(colnames(mainTab) %in% cols)
mainTab <- mainTab[, cols]
save(mainTab, varDesc, demographicVariables, file='./sample_data.Rdata')