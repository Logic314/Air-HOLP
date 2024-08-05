
# This file applies the Air-HOLP, Ridge-HOLP, and SIS methods to the Tomlins-V2
# prostate cancer genetic data and produces the plots in section 5 of the paper.
# The code is divided into 4 parts:
# Preparations, Functions, Results, and Line plots.

options(warn=-1)
if (!require("rstudioapi")) install.packages("rstudioapi")
library(rstudioapi)
if(!require("cmna")){install.packages("cmna")}
library(cmna)
if(!require("ggplot2")){install.packages("ggplot2")}
library(ggplot2)
if(!require("dplyr")){install.packages("dplyr")}
library(dplyr)
if(!require("MASS")){install.packages("MASS")}
library(MASS)
if(!require("greybox")){install.packages("greybox")}
library(greybox)
set.seed(1)

current_path <- rstudioapi::getActiveDocumentContext()$path
current_dir <- dirname(current_path)
setwd(current_dir)

source("Adaptive_iterative_regularization_functions.R")


# Preparations
#_______________________________________________________________________________

# Extracting data
Data <- read.csv("Tomlins-V2.csv",
                 header = FALSE, row.names = NULL)
IDs <- Data[1,-1]
Gene_names <- Data[c(-1,-2),1]
response <- Data[2,-1]
predictors <- as.matrix(Data[c(-1,-2),-1])
X <- scale(t(matrix(as.numeric(predictors), ncol = 92)))
y1 <- t(ifelse(response == "EPI", 1, 0))
y2 <- t(ifelse(response == "PIN", 1, 0))
y3 <- t(ifelse(response == "PCA", 1, 0))
y4 <- t(ifelse(response == "MET", 1, 0))
n <- nrow(X)
p <- ncol(X)


# Functions
#_______________________________________________________________________________

# Sure Independence Screening
SIS <- function(X, y) {
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X)
  y <- scale(y)
  Beta <- crossprod(X, y)
  SIS <- rank(-abs(Beta), na.last = NA, ties.method = "random")
}

# Maximum multiple correlation coefficient for a given model size
max_mcor <- function(X, y, indices, size) {
  combinations <- combn(indices, size)
  max_mcor_value <- -Inf
  max_combination <- NULL
  
  for (i in 1:ncol(combinations)) {
    current_combination <- combinations[, i]
    mcor_value <- mcor(X[, current_combination], y)$value
    
    if (mcor_value > max_mcor_value) {
      max_mcor_value <- mcor_value
      max_combination <- current_combination
    }
  }
  
  return(c(max_mcor_value, max_combination))
}

# Converting indices from the screened set to the full set
convert_indices <- function(max_mcors, original_indices) {
  for (size in 1:m) {
    max_mcors[[size]][-1] <- original_indices[max_mcors[[size]][-1]]
  }
  return(max_mcors)
}

# Creating a line plot
create_line_plot <- function(plot_data, Title) {
  p <- ggplot(plot_data, aes(x = Size, y = Value, color = Method,
                             shape = Method)) +
    geom_line(size = 1, aes(linetype = Method), alpha = 0.8) +
    geom_point(size = 4, alpha = 0.9, fill = '#CBC3E3') +
    scale_color_manual(values = c("Air-HOLP" = "#008689", "Ridge-HOLP" =
                                    "#F8766D", "SIS" = "#a30a90")) +
                                    scale_linetype_manual(values = c("Air-HOLP" = "solid", "Ridge-HOLP" =
                                                                       "dashed", "SIS" = "dotted")) +
    scale_shape_manual(values = c("Air-HOLP" = 19, "Ridge-HOLP" = 15,
                                  "SIS" = 21)) +
    labs(
      title = Title,
      x = "Model Size",
      y = "Maximum Multiple Correlation Coefficient",
      shape = "Method"
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
    ) +
    guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1))
  print(p)
  return(p)
}


# Results
#_______________________________________________________________________________

# Settings
m <- 8 # maximum model size in the plots
Figures <- c("EPI", "PIN", "PCA", "MET")
for(Figure in Figures){
  y <- switch(Figure,
              "EPI" = y1,
              "PIN" = y2,
              "PCA" = y3,
              "MET" = y4,
              y3
  )
  
  # Ranking features
  Threshold <- ceiling(n / log(n))
  temp <- AirHOLP(X, y, Threshold)
  ranks_r <- temp$index_r
  ranks_r0 <- temp$index_r0
  ranks_SIS <- SIS(X, y)
  
  # Indices of screened features
  indices_r <- match(1:Threshold, ranks_r)
  indices_r0 <- match(1:Threshold, ranks_r0)
  indices_SIS <- match(1:Threshold, ranks_SIS)
  
  # Screened features
  X_r <- X[, indices_r, drop = FALSE] # Air-HOLP
  X_r0 <- X[, indices_r0, drop = FALSE] # Ridge-HOLP
  X_SIS <- X[, indices_SIS, drop = FALSE] # SIS
  
  # Maximum multiple corelation coefecients
  max_mcors_r <- list()
  max_mcors_r0 <- list()
  max_mcors_SIS <- list()
  for (size in 1:m) {
    max_mcors_r[[size]] <- max_mcor(X_r, y, 1:Threshold, size)
    max_mcors_r0[[size]] <- max_mcor(X_r0, y, 1:Threshold, size)
    max_mcors_SIS[[size]] <- max_mcor(X_SIS, y, 1:Threshold, size)
  }
  max_mcors_r <- convert_indices(max_mcors_r, indices_r)
  max_mcors_r0 <- convert_indices(max_mcors_r0, indices_r0)
  max_mcors_SIS <- convert_indices(max_mcors_SIS, indices_SIS)
  
  
  # Line plots
  #_______________________________________________________________________________
  
  # Preparing data for plotting
  prepare_data <- function(max_mcors_r, max_mcors_r0, max_mcors_SIS, m) {
    sizes <- 1:m
    values_r <- sapply(max_mcors_r, function(x) x[1])
    values_r0 <- sapply(max_mcors_r0, function(x) x[1])
    values_SIS <- sapply(max_mcors_SIS, function(x) x[1])
    
    data <- data.frame(
      Size = rep(sizes, 3),
      Method = rep(c("Air-HOLP", "Ridge-HOLP", "SIS"), each = m),
      Value = c(values_r, values_r0, values_SIS)
    )
    
    return(data)
  }
  
  # Creating the plot
  plot_data <- prepare_data(max_mcors_r, max_mcors_r0, max_mcors_SIS, m)
  Title <- switch(Figure,
                  "EPI" = "Benign Epithelium",
                  "PIN" = "Prostatic Intraepithelial Neoplasia",
                  "PCA" = "Prostate Cancer",
                  "MET" = "Metastatic Disease",
                  "Prostate Cancer"
  )
  line_plot <- create_line_plot(plot_data, Title)
  print(line_plot)
  
  # Saving the plot
  switch(Figure,
         "EPI" = png("Benign Epithelium.png",
                     width = 600, height = 500, pointsize = 20),
         "PIN" = png("Prostatic Intraepithelial Neoplasia.png",
                     width = 600, height = 500, pointsize = 20),
         "PCA" = png("Prostate Cancer.png",
                     width = 600, height = 500, pointsize = 20),
         "MET" = png("Metastatic Disease.png",
                     width = 600, height = 500, pointsize = 20),
         png("Prostate Cancer.png",
             width = 600, height = 500, pointsize = 20)
  )
  print(line_plot)
  dev.off()
}
