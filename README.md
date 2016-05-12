# xwas

[![Build Status](https://travis-ci.org/nampho2/xwas.svg?branch=master)](https://travis-ci.org/nampho2/xwas)

X-Wide Association Analysis (XWAS) is an R-package.

TODO
- sample data files: HIV and NHANES
- add introspection to detect data types?
- script to collect results — use https://github.com/dgrtwo/broom
- script to compute FDR
- script to wrap the above into one command
- add in support for fixed effect/random effect meta analyses
- correlation heatmap
- get to work with HIV dataset
- repeated measures

DEVELOPERS

```
# make edits and test the build locally
xwas$ make

# install and run
xwas$ make install
xwas$ R

# test the code
> library(xwas)
>

# commit changes and push to kick off travis build
xwas$ git add <list of files>
xwas$ git commit -m "description of changes"
xwas$ git push
```

REFERENCES
- Karl Broman "Getting your R package on CRAN" [<a href="http://kbroman.org/pkg_primer/pages/cran.html">web</a>]
- Karl Broman "Writing R packages: Tools for Reproducible Research" [<a href="http://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf">pdf</a>]
- R Admin Pages "CRAN" [<a href="https://cran.r-project.org/doc/manuals/r-release/R-admin.html">web</a>]
- David Diez "Building R Packages" [<a href="http://www.hsph.harvard.edu/statinformatics/soft/files/buildingrpackages.pdf">pdf</a>] 
- Friedrich Leisch "Creating R Packages: A Tutorial" [<a href="https://cran.r-project.org/doc/contrib/Leisch-CreatingPackages.pdf">pdf</a>]
- Alyssa Frazee "RSkittleBrewer" [<a href="https://github.com/alyssafrazee/RSkittleBrewer">web</a>]
