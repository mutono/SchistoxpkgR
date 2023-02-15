# SchistoxpkgR

This R library is a wrapper for the Julia package Schistoxpkg (https://github.com/mattg3004/Schistoxpkg.jl) developed by Matt Graham.

## Installation

This package uses JuliaCall in R in order to run the Julia package. If you are using Ubuntu or Debian, then it is likely that JuliaCall will not work correctly as described here: https://github.com/Non-Contradiction/JuliaCall/issues/99. There is a workaround for this problem described here: https://github.com/Non-Contradiction/JuliaCall/pull/143.

To install the master branch of the package you need to have devtools installed. If this isn't installed, then install it with:
```R
> install.packages("devtools")
```
After this, you can run the following to install the package:
```R
> library(devtools) 
> devtools::install_github('mutono/SchistoxR', build_vignettes=T)
```
This will build a vignette, which will take some time, so the installation of the package may take a few minutes.

After installation is complete, to view the vignette run the following:
```R
> browseVignettes("SchistoxR")
```


