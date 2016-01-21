# build package
all: check

doc:
	R -e 'library(devtools);document(roclets=c("namespace", "rd"))'

build:
	R CMD build --no-build-vignettes --no-manual --no-resave-data .

check: build
	R CMD check --no-build-vignettes --no-manual --timings --as-cran .

install: 
	R CMD install -l ~/.R_libs/ xwas_*.tar.gz
