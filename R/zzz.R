#This file recommended by Hadley, here: http://r-pkgs.had.co.nz/r.html

.onLoad <- function(libname, pkgname){
  packageStartupMessage("=====================================")
  packageStartupMessage("==--------Welcome to ADePTR--------==")
  packageStartupMessage("==--Acoustic Detection Processing--==")
  packageStartupMessage("==---and Visualization Tool in R---==")
  packageStartupMessage("=====================================")
}

.onUnload <- function(libname, pkgname){

}
