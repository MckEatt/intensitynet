# intensitynet: R package for the estimation network intensity functions 

R Package for the computation of edgewise and mean nodewise intensity functions for point patterns recorded over planar network structures. 
Package provides artificial data on crimes committed on a traffic net, the adjacency matrix for the network and the geocoded locations of nodes on the network.


For installation, please execute the following lines in the R:

install.packages('devtools',dependencies=T);
library(devtools);
install_github('MckEatt/intensitynet')

# Application
' load adjacency structure, georeference events and georeferenced nodes 
library(intensitynet)
data(Castellon)
data(crimes)
data(nodes)

' subset of events
crim <- crimes[11:111,]
intensitynet:::netintensity(Castellon, nodes$cx, nodes$cy, crim$X, crim$Y)

