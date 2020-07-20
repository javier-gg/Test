package.name<-"TestPackage"
package.dir<-paste("C:/Users/xavie/Desktop/R package - Version 1/", package.name, sep="")
setwd(package.dir)

### --- Use Roxygenise to generate .RD files from my comments
library(roxygen2)
roxygenise(package.dir = package.dir)
system(command = paste("R CMD INSTALL '", package.dir,"'", sep=""))
