# Installing the packages

To be able to install the package from the R console, you need to install the `devtools` packages.
For installing the package which contains c++ source code, you also need to install the `Rcpp` package:
  
	R> install.packages('devtools')
	R> install.packages('Rcpp')

To install the package directly from GitHub, run this:

	R> library(devtools)
	R> devtools::install_github("veneficusnl/vf_model_time_series")

To check whether you have successfully installed the package, try to load the package:
  
	R> library('vfmodels.timeseries')
