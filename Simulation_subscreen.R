rm(list = ls())
######## Data Prep and libraries #####################################################
library(medicaldata)
library(ggplot2)
library(dplyr)
# IMPORTANT:
remove.packages('subscreen') # first remove current version of subscreen
devtools::install_github('veltrica/subscreen') # Install from github repo
library(subscreen)

data(indo_rct)

# create categorical variables
indo_rct$risk <- cut(indo_rct$risk, breaks = c(0,2,4, 5.5))
indo_rct$age <- cut(indo_rct$age, breaks = c(0,30,60,90))

# remove treatment indicator and endpoint
indo_rct <- select(indo_rct, -c("rx","outcome"))
# treatment indicator and endpoint have to be initialized again outside of function
indo_rct$T <- 1:nrow(indo_rct)
indo_rct$endpoint <- 1:nrow(indo_rct)


SimRun <- function(n_support_points = 50, kernel = 2){
  ## Get SubscreenResult
  # Simulate Treatment Indicator and Endpoint
  dat <- indo_rct
  dat$T <- sample(c(rep(0,301),rep(1,301)), nrow(dat), replace = FALSE)
  dat$endpoint <- rnorm(nrow(dat))
  dat[dat$T == 1,]$endpoint = dat[dat$T == 1,]$endpoint + 1
  
  # Calculate Subscreen result
  
  Result <- subscreencalc(data = dat,
                          eval_function= ate,
                          endpoints = c("endpoint"),
                          treat="T",
                          subjectid = 'id',
                          nkernel = kernel)
  
  Result1000 <- append(Result, getPointEstimates(Result, 
                                             dat,'T','endpoint',1000,n_support_points,ate))
  
  Result2000 <- append(Result, getPointEstimates(Result, 
                                                 dat,'T','endpoint',2000,n_support_points,ate))
  Result5000 <- append(Result, getPointEstimates(Result, 
                                                 dat,'T','endpoint',5000,n_support_points,ate))

  
  ResultsTable <<- rbind(ResultsTable, matrix(c(50,
                                              points_outside(Result1000, alpha = 0.01),
                                              points_outside(Result1000, alpha = 0.05),
                                              points_outside(Result1000, alpha = 0.1),
                                              points_outside(Result2000, alpha = 0.01),
                                              points_outside(Result2000, alpha = 0.05),
                                              points_outside(Result2000, alpha = 0.1),
                                              points_outside(Result5000, alpha = 0.01),
                                              points_outside(Result5000, alpha = 0.05),
                                              points_outside(Result5000, alpha = 0.1)), ncol = 10
                                               ))
  write.csv(ResultsTable, 'Resultstotal.csv')
  
  
}



# Helper functions
ate <- function(dat) { ate = round((mean(dat[dat$T == 1,]$endpoint) - mean(dat[dat$T == 0,]$endpoint)), digits = 2)
data.frame(ate = ate) }

points_outside <- function(SSR, alpha = 0.05, eval_name = 'ate', threshold = 60){
  
  alphas <- names(quantile(c(0), probs = c(alpha/2, 1-alpha/2)))
  lower <- data.frame(x = as.vector(SSR$nsamp), y = as.vector(SSR$intervals[alphas[1],]), memorizedText = "")
  upper <- data.frame(x = as.vector(SSR$nsamp), y = as.vector(SSR$intervals[alphas[2],]), memorizedText = "")
  
  values <- data.frame(x = SSR$sge$N.of.subjects, y = as.numeric(unlist(SSR$sge[eval_name])))
  values <- values[values$x>threshold,]
  pred_lower <- predict(loess(y ~ x, lower, span = 0.25)
                        , data.frame(x = values$x))
  pred_upper <- predict(loess(y ~ x, upper, span = 0.25)
                        , data.frame(x = values$x))
  below <- pred_lower > values$y
  above <- pred_upper < values$y
  outside <- below + above
  ratio <- sum(outside, na.rm = TRUE)/ length(na.omit(outside))
  ratio
}
getPointEstimates <- function(H, 
                                data,
                                treat,
                                endpoints,
                                nperm,
                                n_support_points,
                                eval_function){
    
    
    start.time = Sys.time() #  Zeit festhalten
    alpha <- c(0.01,0.05,0.1)
    min <- min(H$sge$N.of.subjects)
    max <- max(H$sge$N.of.subjects)
    sampsize <- max
    
    
    # Generate vector of support points between min and max sample sizes
    sqrtvec = seq(sqrt(min), by = (sqrt(max) - sqrt(min))/n_support_points, length.out = (n_support_points+1))
    nsamp <- matrix(round(sqrtvec^2), nrow = 1)
    
    
    # First we only remove covariates, since they are not of interest,
    data_trimmed <- data[,c(treat, endpoints)]
    
    all_samples <- replicate(nperm, slice_sample(data_trimmed, n = sampsize))
    
    
    # Function SliceR: takes as input a smaller number m <= nrow(df), all_samples matrix and calculates evaluation function for this smaller value
    # Sample
    SliceR <- function(dat,m) {
      to_be_sliced <- data.frame(treat = dat[treat], endpoint = dat[endpoints])
      eval_function(to_be_sliced[1:m,])
    }
    
    # QuantileR: gets as input the desired slice (size of support point) and outputs the quantiles of SliceR values
    
    QuantileR <- function(m) quantile(unlist(apply(all_samples,2, SliceR, m = m)), probs = c(alpha/2, 1-(alpha/2)), na.rm = TRUE)
    quantiles <- apply(nsamp, 2, QuantileR)
    
    point_estimates <- list(intervals = quantiles, nsamp = nsamp)
    
    time.taken = Sys.time() - start.time
    print(time.taken)
    point_estimates
  }
  


# Simulation Run

set.seed(63)
apply(matrix(rep(50,500), nrow = 1),2, SimRun, kernel = 2) # 

Results <- ResultsTable[2:501,]


##### RESULTS (Mean Outcome per number of permutations in percentages)
100*apply(Results,2,mean) 

std <- function(x) sd(x)/sqrt(length(x))

##### STANDARD ERRORS OF MEAN in percentages
100*apply(Results, 2, std)
apply(Results, 2, sd)

