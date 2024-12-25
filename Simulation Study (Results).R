
# This file conducts the compound symmetry simulation study in section 4 and
# saves the results in the files "Threshold results.csv" and
# "Probability results.csv".
# The code is divided into 3 parts:
# Preparations, Simulations, and Results.

options(warn=-1)
if (!require("rstudioapi")) install.packages("rstudioapi")
library(rstudioapi)
if(!require("mvnfast")){install.packages("mvnfast")}
library(mvnfast)
if(!require("cmna")){install.packages("cmna")}
library(cmna)
nc <- 10 # number of cores used by the rmvn function
set.seed(1)

current_path <- rstudioapi::getActiveDocumentContext()$path
current_dir <- dirname(current_path)
setwd(current_dir)

source("Adaptive_iterative_regularization_functions.R")


# Preparations
#_______________________________________________________________________________

# Simulation settings
r0 <- 10 # Penalty parameter of Ridge-HOLP
N <- c(125, 250, 500, 1000) # Sample size
P <- c(250, 1250, 5000, 15000) # Number of predictors
Rho <- c(0, 0.3, 0.6, 0.9) # level of correlation between features
P0 <- c(3, 6, 9, 12, 15) # Number of true features
Rsq <- c(0.25, 0.5, 0.75, 0.9, 0.95) # Theoretical R squared
BX <- 50 # Number of simulations of X
By <- 10 # Number of simulations of y for each X

# Simulation settings for individual plots
Figure <- "All" # "All", "1a", "1b", "2", "3a", "3b"
if(Figure == "1a"){
  N <- 125; P <- 250; P0 <- 9; Rsq <- 0.95
}
if(Figure == "1b"){
  N <- 1000; P <- 1250; P0 <- 15; Rsq <- 0.5
}
if(Figure == "2"){
  Rho <- 0.6; P0 <- 15; Rsq <- 0.5
}
if(Figure == "3a"){
  N <- 500; P <- 1250; Rho <- 0.6; P0 <- 12
}
if(Figure == "3b"){
  N <- 1000; P <- 5000; Rho <- 0.3; Rsq <- 0.25
}

# Creating storage matrices for the results
size <- length(Rho)*length(N)*length(P)*length(P0)*length(Rsq)*BX*By
size2 <- c(BX*By,length(Rho),length(N),length(P),length(P0),length(Rsq))
R <- array(1:size - 1:size, size2) # Selected penalties
Thr_r <- R # Sure screening thresholds (Air-HOLP)
Thr_r0 <- R # Sure screening thresholds (Ridge-HOLP)
Thr_SIS <- R # Sure screening thresholds (SIS)


# Simulations
#_______________________________________________________________________________

start_time <- Sys.time()

for (bx in 1:BX) {
  for (j in 1:length(N)) {
    n <- N[j]
    for (k in 1:length(P)) {
      p <- P[k]
      # Generating and standardizing independent features:
      Z <- matrix(rmvn(n*(p+1), 0, 1, ncores = nc), nrow = n, ncol = p+1)
      Z <- scale(Z)
      # Calculations on the independent features:
      ZZT <- tcrossprod(Z[,1:p])
      Z1zT <- rowSums(Z[,1:p])%*%t(Z[,p+1])
      zzT <- tcrossprod(Z[,p+1])
      for (i in 1:length(Rho)) {
        rho <- Rho[i]
        # Generating the features matrix X and XXT
        X <- sqrt(1 - rho)*Z[,1:p] + sqrt(rho)*kronecker(matrix(1,1,p),Z[,p+1])
        XXT <- (1-rho)*ZZT + sqrt(rho*(1-rho))*(Z1zT + t(Z1zT)) + rho*p*zzT
        # Calculations on X:
        eXXT <- eigen(XXT)
        Lambda <- eXXT$values
        U <- eXXT$vectors
        XU <- crossprod(X,U)
        for (by in 1:By) {
          er <- rnorm(n) # error term
          XTer <- crossprod(X,er)
          for (l in 1:length(P0)) {
            p0 <- P0[l]
            X0 <- X[,1:p0] # matrix of true features
            Beta0 <- 2*((runif(p0)>0.4)-0.5)*(abs(rnorm(p0))+4*log(n)/sqrt(n))
            y0 <- X0%*%Beta0 # expected response vector
            XTy0 <- crossprod(X,y0)
            for (m in 1:length(Rsq)) {
              y0_var <- (1 - rho)*crossprod(Beta0) + rho*(sum(Beta0))^2
              se <- sqrt((1-Rsq[m])/Rsq[m]*y0_var) # Standard error
              y <- y0 + se[1]*er # response vector
              # SIS
              XTy <- XTy0 + se[1]*XTer
              rank_SIS <- rank(-abs(XTy), na.last = NA,
                               ties.method = c("random"))
              Thr_SIS[(bx-1)*By+by,i,j,k,l,m] <- max(rank_SIS[1:p0])
              # Air-HOLP and Ridge-HOLP
              temp <- AirHOLP(X, y, ceiling(n/log(n)), r0 = r0,
                              Lambda = Lambda, U = U, XU = XU)
              rank_r <- temp$index_r
              rank_r0 <- temp$index_r0
              R[(bx-1)*By+by,i,j,k,l,m] <- temp$r
              Thr_r[(bx-1)*By+by,i,j,k,l,m] <- max(rank_r[1:p0])
              Thr_r0[(bx-1)*By+by,i,j,k,l,m] <- max(rank_r0[1:p0])
            }
          }
        }
      }
    }
  }
}

end_time <- Sys.time()
end_time - start_time


# Results
#_______________________________________________________________________________

# Simulation settings
Settings <- as.data.frame(matrix(1:(size*7) - 1:(size*7),
                                 nrow = size, ncol = 7))
names(Settings) <- c("y_id", "X_id", "Rho", "N", "P", "P0", "Rsq")
size_temp <- By
Settings[,1] <- rep(1:By, times = size/size_temp, each = size_temp/By)
size_temp <- size_temp*BX
Settings[,2] <- rep(1:BX, times = size/size_temp, each = size_temp/BX)
size_temp <- size_temp*length(Rho)
Settings[,3] <- rep(Rho, times = size/size_temp, each = size_temp/length(Rho))
size_temp <- size_temp*length(N)
Settings[,4] <- rep(N, times = size/size_temp, each = size_temp/length(N))
size_temp <- size_temp*length(P)
Settings[,5] <- rep(P, times = size/size_temp, each = size_temp/length(P))
size_temp <- size_temp*length(P0)
Settings[,6] <- rep(P0, times = size/size_temp, each = size_temp/length(P0))
size_temp <- size_temp*length(Rsq)
Settings[,7] <- rep(Rsq, times = size/size_temp, each = size_temp/length(Rsq))

# Raw results
Results <- as.data.frame(matrix(1:(size*4) - 1:(size*4),
                                nrow = size, ncol = 4))
names(Results) <- c("Sure Screening Threshold (Air-HOLP)",
                    "Sure Screening Threshold (Ridge-HOLP)",
                    "Sure Screening Threshold (SIS)", "r")
Results[,1] <- as.vector(Thr_r)
Results[,2] <- as.vector(Thr_r0)
Results[,3] <- as.vector(Thr_SIS)
Results[,4] <- as.vector(R)

# Saving data
Threshold_results <- cbind(Settings, Results)
write.csv(Threshold_results, "Threshold results.csv")

# Mean Settings
Settings.probability <- as.data.frame(matrix(1:(size/BX/By*5)-1:(size/BX/By*5),
                                             nrow = size/BX/By, ncol = 5))
names(Settings.probability) <- c("Rho", "N", "P", "P0", "Rsq")
size_temp <- length(Rho)
Settings.probability[,1] <- rep(Rho, times = size/BX/By/size_temp,
                                each = size_temp/length(Rho))
size_temp <- size_temp*length(N)
Settings.probability[,2] <- rep(N, times = size/BX/By/size_temp,
                                each = size_temp/length(N))
size_temp <- size_temp*length(P)
Settings.probability[,3] <- rep(P, times = size/BX/By/size_temp,
                                each = size_temp/length(P))
size_temp <- size_temp*length(P0)
Settings.probability[,4] <- rep(P0, times = size/BX/By/size_temp,
                                each = size_temp/length(P0))
size_temp <- size_temp*length(Rsq)
Settings.probability[,5] <- rep(Rsq, times = size/BX/By/size_temp,
                                each = size_temp/length(Rsq))

# Mean Results
Results.probability <- as.data.frame(matrix(1:(size/BX/By*4) - 1:(size/BX/By*4),
                                            nrow = size/BX/By, ncol = 4))
names(Results.probability) <- c("Sure Screening Probability (Air-HOLP)",
                                "Sure Screening Probability (Ridge-HOLP)",
                                "Sure Screening Probability (SIS)", "r")
N_array <- array(Settings$N, dim = size2)
M <- ceiling(N_array/log(N_array))
Results.probability[,1] <- as.vector(colMeans(Thr_r <= M, dims = 1))
Results.probability[,2] <- as.vector(colMeans(Thr_r0 <= M, dims = 1))
Results.probability[,3] <- as.vector(colMeans(Thr_SIS <= M, dims = 1))
Results.probability[,4] <- as.vector(colMeans(R, dims = 1))

# Saving data
Probability_results <- cbind(Settings.probability, Results.probability)
write.csv(Probability_results, "Probability results.csv")
if(Figure != "All"){
  print(Probability_results)
}
