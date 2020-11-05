# SchistoxpkgR



## Installation


To install the master branch of the package you need to have devtools installed. If this isn't installed, then install it with:
```R
install.packages("devtools")
```
After this, you can run the following to install the package:
```R
library(devtools) 
devtools::install_github('mattg3004/SchistoxpkgR', build_vignettes=T)
```
This will build a vignette, which will take some time, so the installation of the package may take a few minutes.

After installation is complete, to view the vignette run the following:
```R
browseVignettes("SchistoxpkgR")
```
