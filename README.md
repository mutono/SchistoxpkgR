# SchistoxpkgR



## Installation


To install the master branch of the package you need to have devtools installed. If this isn't installed, then install it with:
```R
install.packages("devtools")
```
After this, you can run:
```R
library(devtools) 
devtools::install_github('mattg3004/SchistoxpkgR', build_vignettes=T)
```
to install the package. This will build a vignette, which will take some time, so the installation of the package may take a few minutes.

After that you can run:
```R
browseVignettes("SchistoxpkgR")
```
to view the vignette.
