# initialize libraries for graph theory analysis
library(latentnet)
library(network)
library(netdata)
library(ndtv)
library(sna)
library(ergm)
library(tergm)
library(snow)
library(snowFT)
library(plyr)
library(stargazer)
library(rgl)
library(psych)
library(mnormt)

# sub my variables
# G = 2: 2 groups

# input: connectivity matrix. 

cMatrix <- replicate(10, rnorm(10)) 
cMatrix <- 1-read.csv("~/Desktop/UCLA/scratchpad/LanguageMVPA/Workbook1.csv", header=F)
subData <- network(cMatrix)
latentSub <- ergmm(subData ~ euclidean(d=2), verbose = TRUE)
# commands after getting model
plot(latentSub)
# get probabilities of group membership
latentSub$mkl 
# Z coordinates for each node:
latentSub$mkl$Z
# tab to autocomplete

