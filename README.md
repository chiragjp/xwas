# xwas

[![Build Status](https://travis-ci.org/nampho2/xwas.svg?branch=master)](https://travis-ci.org/nampho2/xwas)

X-Wide Association Analysis (XWAS) is an R-package.

##### TODO
- sample data files: HIV
- add introspection to detect data types
- script to collect results — use https://github.com/dgrtwo/broom
- add in support for fixed effect/random effect meta analyses
- correlation heatmap
- get to work with HIV dataset
- repeated measures

##### End-users
Currently install directly from the Github repository using `devtools`, future support for CRAN installation forthcoming.

```
# install directly from Github
> library(devtools)
> devtools::install_github("nampho2/xwas")
Downloading GitHub repo nampho2/xwas@master
```   

##### Developers
Make sure to set `R_LIBS="~/.R_libs"` within a `.Renviron` file within your home directory. Also make sure the `~/.R_libs` folder exists to house personal packages distinct from the system set of R libraries.

```
# build, check, and install
xwas$ make deploy

# test the code
xwas$ R
> library(xwas)
>

# commit changes and push to kick off travis build
xwas$ git add <list of files>
xwas$ git commit -m "description of changes"
xwas$ git push
```

##### References
- Karl Broman "Getting your R package on CRAN" [<a href="http://kbroman.org/pkg_primer/pages/cran.html">web</a>]
- Karl Broman "Writing R packages: Tools for Reproducible Research" [<a href="http://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf">pdf</a>]
- R Admin Pages "CRAN" [<a href="https://cran.r-project.org/doc/manuals/r-release/R-admin.html">web</a>]
- David Diez "Building R Packages" [<a href="http://www.hsph.harvard.edu/statinformatics/soft/files/buildingrpackages.pdf">pdf</a>] 
- Friedrich Leisch "Creating R Packages: A Tutorial" [<a href="https://cran.r-project.org/doc/contrib/Leisch-CreatingPackages.pdf">pdf</a>]
- Alyssa Frazee "RSkittleBrewer" [<a href="https://github.com/alyssafrazee/RSkittleBrewer">web</a>]
