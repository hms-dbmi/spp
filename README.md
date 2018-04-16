# ChIP-seq processing pipeline
An [R](https://www.r-project.org/) package for analysis of ChIP-seq and other functional sequencing data.
Please see [package homepage](http://compbio.med.harvard.edu/Supplements/ChIP-seq/) for details.

## Requirements
A unix-flavored OS with R (>= 3.3.0) installed.

## Installation
Since version 1.15.4 spp is available on [CRAN](https://cran.r-project.org/web/packages/spp/index.html)

Alternatively you can use modtools to install spp:

```
require(devtools)
devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)
```

Or 

download a .tar.gz containing the [latest release](https://github.com/hms-dbmi/spp/releases) and use the standard R installation command, e.g.:
```
R CMD INSTALL spp_1.13.tar.gz
```

Note: Since version 1.15.4 the Boost headers are supplied by [BH package](https://cran.r-=project.org/web/packages/BH/index.html) as suggested by CRAN package guidelines.

