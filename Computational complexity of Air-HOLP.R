
# This file conducts the computational complexity comparison in section 6
# between Air-HOLP, Ridge-HOLP, and Ridge-OLS.
# The code is divided into 4 parts:
# Preparations, Simulations, Line plot, and Box plot.

options(warn=-1)
if (!require("rstudioapi")) install.packages("rstudioapi")
library(rstudioapi)
if(!require("cmna")){install.packages("cmna")}
library(cmna)
if(!require("ggplot2")){install.packages("ggplot2")}
library(ggplot2)
set.seed(1)

current_path <- rstudioapi::getActiveDocumentContext()$path
current_dir <- dirname(current_path)
setwd(current_dir)

source("Adaptive_iterative_regularization_functions.R")


# Preparations
#_______________________________________________________________________________

# Simulation settings
n <- 250 # Sample size
P <- c(250, 500, 1000, 1750, 3750, 5000, 6500) # Number of features
p0 <- 6 # Number of true features
Rsq <- 0.75 # Theoretical R squared
num <- 10 # Number of simulation runs per setting

# Basic OLS and HOLP functions
OLS <- function(X, y, r = 10){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- X - matrix(rep(colMeans(X), each = n), n, p)
  X <- X/matrix(rep(sqrt(colMeans(X^2)), each = n), n, p)
  y <- (y - mean(y))/sd(y)
  Beta <- solve(crossprod(X) + r*diag(p), crossprod(X, y))
  rank_OLS <- rank(-abs(Beta), na.last = NA, ties.method = c("random"))
  OLS <- rank_OLS
}
HOLP <- function(X, y, r = 10){
  n <- dim(X)[1]
  p <- dim(X)[2]
  X <- X - matrix(rep(colMeans(X), each = n), n, p)
  X <- X/matrix(rep(sqrt(colMeans(X^2)), each = n), n, p)
  y <- (y - mean(y))/sd(y)
  Beta <- crossprod(X, solve(tcrossprod(X) + r*diag(n), y))
  rank_HOLP <- rank(-abs(Beta), na.last = NA, ties.method = c("random"))
  HOLP <- rank_HOLP
}


# Simulations
#_______________________________________________________________________________

start_time_all <- Sys.time()

Ridge_OLS_time <- matrix((1:(length(P)*num))*0, nrow = num, ncol = length(P))
Ridge_HOLP_time <- Ridge_OLS_time
Air_HOLP_time <- Ridge_OLS_time
for (j in 1:length(P)) {
  p <- P[j]
  for (i in 1:num) {
    X <- matrix(rnorm(n*p, 0, 1), nrow = n, ncol = p)
    X0 <- X[,1:p0] # Matrix of true features
    Beta0 <- 2*((runif(p0)>0.4) - 0.5)*(abs(rnorm(p0))+4*log(n)/sqrt(n))
    y0 <- X0%*%Beta0 # Expected response vector
    se <- sqrt((1-Rsq)/Rsq*var(y0)) # Standard error
    y <- y0 + se[1]*rnorm(n) # Response vector
    start_time <- Sys.time()
    temp <- OLS(X, y)
    end_time <- Sys.time()
    Ridge_OLS_time[i,j] <- end_time - start_time
    start_time <- Sys.time()
    temp <- HOLP(X, y)
    end_time <- Sys.time()
    Ridge_HOLP_time[i,j] <- end_time - start_time
    start_time <- Sys.time()
    temp <- AirHOLP(X, y, ceiling(n/log(n)))
    end_time <- Sys.time()
    Air_HOLP_time[i,j] <- end_time - start_time
  }
}

end_time_all <- Sys.time()
end_time_all - start_time_all

# Saving results
Results <- rbind(round(Ridge_OLS_time*1000),
                 round(Ridge_HOLP_time*1000),
                 round(Air_HOLP_time*1000))
rownames(Results) <- c(rep('Ridge-OLS', num),
                       rep('Ridge-HOLP', num),
                       rep('Air-HOLP', num))
colnames(Results) <- P
Results
write.csv(Results, "Time complexity comparison.csv")

# Loading results
Results <- read.csv("Time complexity comparison.csv")
Results <- Results[,2:8] # Removing id

# Line plot
#_______________________________________________________________________________

# Result matrices
Ridge_OLS_mean_time <- P
Ridge_HOLP_mean_time <- P
Air_HOLP_mean_time <- P
for (i in 1:length(P)) {
  Ridge_OLS_mean_time[i] <- mean(Results[1:num,i], trim = 0.1)
  Ridge_HOLP_mean_time[i] <- mean(Results[(num+1):(2*num),i], trim = 0.1)
  Air_HOLP_mean_time[i] <- mean(Results[(2*num+1):(3*num),i], trim = 0.1)
}
Lplot_df <- data.frame(
  P = rep(P, 3),
  log_execution_time = c(log(Ridge_OLS_mean_time*1000),
                         log(Ridge_HOLP_mean_time*1000),
                         log(Air_HOLP_mean_time*1000)),
  Method = factor(rep(c("Ridge-OLS", "Ridge-HOLP", "Air-HOLP"),
                      each = length(P)))
)

# Generating the line plot
Lplot <- ggplot(Lplot_df, aes(x = P, y = log_execution_time,
                              color = Method, shape = Method)) +
  geom_line(size = 1, aes(linetype = Method), alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Ridge-OLS" = "#b30000",
                                "Ridge-HOLP" = "#F8766D",
                                "Air-HOLP" = "#008689")) +
  scale_linetype_manual(values = c(
    "Ridge-OLS" = "dotted",
    "Ridge-HOLP" = "dashed",
    "Air-HOLP" = "solid"
  )) +
  scale_shape_manual(values = c(
    "Ridge-OLS" = 17,
    "Ridge-HOLP" = 15,
    "Air-HOLP" = 19
  )) +
  labs(
    title = bquote(paste(n == .(n), ", ", rho == 0, ", ", p[0] == .(p0),
                         ", ", R^2 == .(Rsq), sep = "")),
    x = bquote(paste("Number of Features (", p, ")", sep = "")),
    y = "log(Execution Time (ms))",
    color = "Method"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.2, 'cm'),
    legend.key.size = unit(0.5, 'cm'),
    legend.key.width = unit(0.75, 'cm')
  )
print(Lplot)

# Saving the plot
png("Time complexity comparison (line plot).png",
    width = 600, height = 500, pointsize = 20)
print(Lplot)
dev.off()


# Box plot
#_______________________________________________________________________________

# Combining data into a data frame
Bplot_df <- data.frame(
  Method = factor(rep(c("Ridge-HOLP", "Air-HOLP"), each = num)),
  Execution_Time = c(Results[(num+1):(2*num),7]/1000,
                     Results[(2*num+1):(3*num),7]/1000)
)

# Generating the box plot
Bplot <- ggplot(Bplot_df, aes(fill = Method, y = Execution_Time, x = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Air-HOLP" = "#008689",
                               "Ridge-HOLP" = "#F8766D")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.2, 'cm'),
    legend.key.size = unit(0.5, 'cm'),
    legend.key.width = unit(0.75, 'cm')
  ) +
  ggtitle(bquote(paste(n == .(n), ", ", p == .(P[length(P)]), ", ", rho == 0,
                       ", ", p[0] == .(p0), ", ", R^2 == .(Rsq), sep = ""))) +
  xlab("Method") +
  ylab("Execution Time (s)") +
  scale_y_continuous(limits = c(0, 0.5))
print(Bplot)

# Saving the plot
png("Time complexity comparison (box plot).png",
    width = 600, height = 500, pointsize = 20)
print(Bplot)
dev.off()
