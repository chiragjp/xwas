all: build

deploy: build check install

build:
	R -e 'devtools::document()'
	R CMD build --no-build-vignettes --no-manual --no-resave-data .

check:
	R CMD check --no-build-vignettes --no-manual --timings --as-cran .

clean:
	rm -fr ~/.R_libs/xwas

install: clean
	R CMD install -l ~/.R_libs/ xwas_*.tar.gz
