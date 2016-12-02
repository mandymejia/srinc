##### Hackathon 11/19/2016
##### Spaghetti plot of scan length vs ICC (Xixi & David)
# rm(list=ls())

library(devtools)
install_github("mandymejia/shrinkR")
library(shrinkR) # UT2mat
library("scales")

##### Read in ICC data.
icc <- read.csv("ICC.csv", header = F)
# dim(icc)
# [1] 44850     8
# 300, 600,..., 2400 time points
# 44850 = 300choose2


##### Spaghetti plot of scan length vs ICC.
##### Inputs:
##### 1) icc: ICC data matrix, ncol = number of scan lengths; nrow = 300choose2.
##### 2) seed: user-specified seed region in the region pair(s) whose ICC(s) are of interest;
#####    currently taking one value from 1:300.
##### 3) target: user-specified target region(s) in the region pair(s) whose ICC(s) are of interest;
#####    currently taking multiple values from 1:300.

##### iCCs between the seed region and target regions are colored black; 
##### ICCs between all other regions are colored lightgrey.

spaghetti <- function(icc, seed, target = NULL) {
  ### Store ICC in an array: one 300 * 300 ICC matrix for each scan length. 
  icc_array <- array(NA, dim = c(300, 300, dim(icc)[2]))
  for (i in 1: dim(icc)[2]) {
    icc_array[,,i] <- UT2mat(icc[,i])
    diag(icc_array[,,i]) <- NA
  }
  
  ### Plot
  matplot(t(icc), col = rep("lightgrey", 299), type = "l")
  if (length(target) == 1) {
    matplot(icc_array[seed, target,], col = rep("black", dim(icc_array)[1]), type = "l", add = TRUE, lty = 1)
  } 
  if (length(target) > 1) {
    matplot(t(icc_array[seed, target,]), col = rep("black", dim(icc_array)[1]), type = "l", add = TRUE, lty = 1)
  }
}


# Test
# spaghetti(icc, 3, 1:2)
# spaghetti(icc, 3, 2)



##### To be finished-----------------------------------------------------------
### Colors for cross-session icc for chosen regions/networks.
### Regions are the ICA components from HCP data.
# col_region <- function(icc, seed = NULL, target = NULL){
#   number_region <- dim(icc)[2]
#   col_all <- rep("grey", number_region)
#   cor_all[seed,target] <- "black"
#   transparency_all <- rep(0.2, number_region)
#   transparency_all[seed,target] <- 1
# }


