# Goal

This file explains how to create a R package depending on Rcpp, also I think it is the recommended workflow (for me).

* First generate an R package project by Rstudio

* Within the working directory, run usethis::use_rcpp(), as described by http://r-pkgs.had.co.nz/src.html#cpp , it has made the following changes to this R package project

Create a src/ directory to hold your .cpp files.

Add Rcpp to the LinkingTo and Imports fields in the DESCRIPTION.

Set up a .gitignore file to make sure you don’t accidentally check in any compiled files (learn more about this in git).


* Follow the instruction of the function: Add the following code to the end of NAMESPACE
```
  useDynLib('pjnormarcpp', .registration = TRUE)
  importFrom('Rcpp', 'sourceCpp')
```  

Now this R package is good to run with Rcpp features. Just add a default Rcpp file from Rstudio and build&reload.

To further make RcppArmadillo available, I checked the output of RcppArmadillo.package.skeleton(), which basically generated an R package depends on RcppArmadillo. The NAMESPACE file has the same content, but different order. I only make the following change 

* Change the following in the DESCRIPTION file.
Imports: Rcpp (>= 1.0.0)
LinkingTo: Rcpp, RcppArmadillo 

Next I add firstArmadillo.cpp to the src/, which depends on RcppArmadillo, and build&reload without problem.

Next is to consider the header files to allow dependency across cpp files. I followed https://knausb.github.io/2017/08/header-files-in-rcpp/. That is, I copied modString.h, modString.cpp, myFunction.cpp to the src folder and build&reload without problem.

Add the more naive way of working with multiple files by forward declaration.
