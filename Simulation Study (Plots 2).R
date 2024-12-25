
# This file reads the results from the files "Threshold results 2.csv" and
# "Probability results 2.csv" then produces spatial correlation plots.
# The code is divided into 5 parts:
# Preparations, Overview, Heat maps, Box plots, and Line plots.
# The "Selected settings" section in each plot determines which settings are
# shown in the plot.

options(warn=-1)
if (!require("rstudioapi")) install.packages("rstudioapi")
library(rstudioapi)
if(!require("ggplot2")){install.packages("ggplot2")}
library(ggplot2)
if(!require("shadowtext")){install.packages("shadowtext")}
library(shadowtext)
if(!require("reshape2")){install.packages("reshape2")}
library(reshape2)
if(!require("gridExtra")){install.packages("gridExtra")}
library(gridExtra)
if(!require("grid")){install.packages("grid")}
library(grid)

current_path <- rstudioapi::getActiveDocumentContext()$path
current_dir <- dirname(current_path)
setwd(current_dir)


# Preparations
#_______________________________________________________________________________

# Extracting data
Threshold_results <- read.csv("Threshold results 2.csv")
Threshold_results <- Threshold_results[,2:12] # Removing id
Settings <- Threshold_results[,1:7]
Results <- Threshold_results[,8:11]
Probability_results <- read.csv("Probability results 2.csv")
Probability_results <- Probability_results[,2:10] # Removing id
Settings.probability <- Probability_results[,1:5]
Results.mean <- Probability_results[,6:9]
By <- length(table(Settings[,1])) # Number of simulations of y for each X
BX <- length(table(Settings[,2])) # Number of simulations of X
Rho <- as.numeric(names(table(Settings[,3]))) # Correlation between features
N <- as.numeric(names(table(Settings[,4]))) # Sample size
P <- as.numeric(names(table(Settings[,5]))) # Number of features
P0 <- as.numeric(names(table(Settings[,6]))) # Number of true features
Rsq <- as.numeric(names(table(Settings[,7]))) # Theoretical R squared
size <- length(Rho)*length(N)*length(P)*length(P0)*length(Rsq)
Method_Names <- c("Air-HOLP", "Ridge-HOLP", "SIS")

# Functions to identify indices for a selected setting
result <- function(indices){
  res <- (Settings$Rho==Rho[indices[1]])*(Settings$N==N[indices[2]])*
    (Settings$P==P[indices[3]])*(Settings$P0==P0[indices[4]])*
    (Settings$Rsq==Rsq[indices[5]])
  return(res)
}
result.mean <- function(indices){
  res <- (Settings.probability$Rho==Rho[indices[1]])*
    (Settings.probability$N==N[indices[2]])*
    (Settings.probability$P==P[indices[3]])*
    (Settings.probability$P0==P0[indices[4]])*
    (Settings.probability$Rsq==Rsq[indices[5]])
  return(res)
}


# Overview
#_______________________________________________________________________________

print(paste("Difference between Air-HOLP and Ridge-HOLP is between ",
            min(Results.mean[,1] - Results.mean[,2]), " and ",
            max(Results.mean[,1] - Results.mean[,2]), sep = ""))

print(paste("Propotion of simulation settings where Air-HOLP has ",
            "higher sure screening probability than Ridge-HOLP = ",
            sum(Results.mean[,1] - Results.mean[,2] > 0)/size, sep = ""))

print(paste("Propotion of simulation settings where Air-HOLP has ",
            "equal sure screening probability to Ridge-HOLP = ",
            sum(Results.mean[,1] - Results.mean[,2] == 0)/size, sep = ""))

print(paste("Propotion of simulation settings where Air-HOLP has ",
            "lower sure screening probability than Ridge-HOLP = ",
            sum(Results.mean[,1] - Results.mean[,2] < 0)/size, sep = ""))


# Heat maps
#_______________________________________________________________________________

# Selected settings (n vs. p)
Rho_id <- 3
P0_id <- 2
Rsq_id <- 2

# Result matrices
air_holp_results <- matrix(nrow = length(N), ncol = length(P))
ridge_holp_results <- matrix(nrow = length(N), ncol = length(P))
for (a in 1:length(N)) {
  for (b in 1:length(P)) {
    index <- which(result.mean(c(Rho_id, a, b, P0_id, Rsq_id)) > 0)
    air_holp_results[a, b] <- Results.mean[index, 1]
    ridge_holp_results[a, b] <- Results.mean[index, 2]
  }
}
difference_results <- air_holp_results - ridge_holp_results
air_holp_df <- melt(air_holp_results)
colnames(air_holp_df) <- c("n", "p", "value")
ridge_holp_df <- melt(ridge_holp_results)
colnames(ridge_holp_df) <- c("n", "p", "value")
difference_df <- melt(difference_results)
colnames(difference_df) <- c("n", "p", "value")

# Color scales
color_scale_diff <- scale_fill_gradientn(
  colours = colorRampPalette(c("#c5370c", "#DB8267", "#E29B86", "white",
                               "#9EC49C", "#85B583", "#3c8839"))(100),
  limits=c(-1, 1), name = "Differece"
)
color_scale_common <- scale_fill_gradientn(
  colours = colorRampPalette(c("white", "#9EC49C", "#85B583", "#3c8839"))(100),
  limits=c(0, 1), name = "Probability"
)

# Plot Air-HOLP heat map
air_holp_plot <- ggplot(air_holp_df, aes(x = n, y = p, fill = value)) + 
  geom_tile(color = "white") +
  color_scale_common +
  geom_shadowtext(aes(label = round(value, 2)), color = "black",
                  bg.color = "white", size = 5.5, bg.r = 0.075) +
  labs(title = Method_Names[1], x = "n", y = "p", fill = "Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 24), 
        text = element_text(size = 16),
        axis.title.x = element_text(face = "bold", size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 18, angle = 45, hjust = 1),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
  scale_x_continuous(breaks = 1:length(N), labels = N) +
  scale_y_continuous(breaks = 1:length(P), labels = P)

# Plot Ridge-HOLP heat map
ridge_holp_plot <- ggplot(ridge_holp_df, aes(x = n, y = p, fill = value)) + 
  geom_tile(color = "white") +
  color_scale_common +
  geom_shadowtext(aes(label = round(value, 2)), color = "black",
                  bg.color = "white", size = 5.5, bg.r = 0.075) +
  labs(title = Method_Names[2], x = "n", y = "p", fill = "Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 24), 
        text = element_text(size = 16),
        axis.title.x = element_text(face = "bold", size = 22), 
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 18, angle = 45, hjust = 1),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
  scale_x_continuous(breaks = 1:length(N), labels = N) +
  scale_y_continuous(breaks = 1:length(P), labels = P)

# Plot difference heat map
difference_plot <- ggplot(difference_df, aes(x = n, y = p, fill = value)) + 
  geom_tile(color = "white") +
  color_scale_diff +
  geom_shadowtext(aes(label = round(value, 2)), color = "black",
                  bg.color = "white", size = 5.5, bg.r = 0.075) +
  labs(title = "Difference",
       #labs(title = paste(Method_Names[1], " - ", Method_Names[2], sep = ""),
       x = "n", y = "p", fill = "Difference") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 24), 
        text = element_text(size = 16),
        axis.title.x = element_text(face = "bold", size = 22), 
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(face = "bold", size = 22),
        axis.text.y = element_text(size = 18, angle = 45, hjust = 1),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
  scale_x_continuous(breaks = 1:length(N), labels = N) +
  scale_y_continuous(breaks = 1:length(P), labels = P)

# Combine the plots into one figure
combined_plot <- grid.arrange(
  air_holp_plot, ridge_holp_plot, difference_plot, 
  nrow = 1,
  top = textGrob(bquote(paste("Sure Screening Probability (",
                              rho == .(Rho[Rho_id]), ", ",
                              p[0] == .(P0[P0_id]), ", ",
                              R^2 == .(Rsq[Rsq_id]), ")", sep = "")),
                 gp = gpar(fontsize = 26, fontface = "bold"))
)
print(combined_plot)

# Save plot
png("Heatmap (n vs p) 2.png",
    width = 1200, height = 350, pointsize = 24)
combined_plot <- grid.arrange(
  air_holp_plot, ridge_holp_plot, difference_plot, 
  nrow = 1,
  top = textGrob(bquote(paste("Sure Screening Probability (",
                              rho == .(Rho[Rho_id]), ", ",
                              p[0] == .(P0[P0_id]), ", ",
                              R^2 == .(Rsq[Rsq_id]), ")", sep = "")),
                 gp = gpar(fontsize = 26, fontface = "bold"))
)
print(combined_plot)
dev.off()


# Box plots
#_______________________________________________________________________________

# Number of true features box plot:

# Selected settings (P0)
Rho_id <- 2
N_id <- 4
P_id <- 3
Rsq_id <- 1

# Preparing results:
holp_plots <- Results[which(result(c(Rho_id,N_id,P_id,1,Rsq_id))>0),1:2]
if(length(P0)>1){
  for (i in 2:length(P0)) {
    holp_plots <- rbind(
      holp_plots, 
      Results[which(result(c(Rho_id,N_id,P_id,i,Rsq_id))>0),1:2]
    )
  }
}
holp_plots_v <- c(holp_plots[,1], holp_plots[,2])
ids <- rep((rep(c(0,length(P0)),
                times = length(P0)) + rep(1:length(P0), each = 2)-1)*BX*By,
           each = BX*By) + rep(1:(BX*By), 2*length(P0))
holp_plots_v <- holp_plots_v[ids]
log_holp_plots_v <- log(holp_plots_v)
holp_names <- LETTERS[1:length(P0)] 
labels <- sapply(P0, function(p) bquote(p[0] == .(p))) 
Bplot1_df <- data.frame(setting = rep(holp_names, each = 2*BX*By),
                        Method = rep(Method_Names[1:2],
                                     times = length(P0), each = BX*By),
                        log_holp_plots_v)

# Creating the plot
hline_data <- data.frame(
  yintercept = c(log(ceiling(N[N_id] / log(N[N_id]))), log(P[P_id])),
  line_label = c("n/log(n)", "p")
)
Bplot1 <- ggplot(Bplot1_df, aes(fill = Method, y = log_holp_plots_v,
                                x = setting)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Air-HOLP" = "#008689",
                               "Ridge-HOLP" = "#F8766D"), name = NULL) +
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "top"
  ) +
  ggtitle(bquote(paste(n == .(N[N_id]), ", ", p == .(P[P_id]), ", ",
                       rho == .(Rho[Rho_id]), ", ", R^2 == .(Rsq[Rsq_id]),
                       sep = ""))) +
  xlab(expression(paste("Number of true features (", p[0], ")", sep = ""))) +
  ylab("log(Sure Screening Threshold)") +
  geom_hline(data = hline_data, aes(yintercept = yintercept, color = line_label,
                                    linetype = line_label), size = 0.5) +
  scale_linetype_manual(
    name = NULL, values = c("n/log(n)" = "dashed", "p" = "dotted")
  ) +
  scale_color_manual(
    name = NULL, values = c("n/log(n)" = "#b30000", "p" = "black")
  ) +
  scale_x_discrete(labels = labels)
print(Bplot1)

# Saving the plot
png("Air-HOLP vs Ridge-HOLP (p0) 2.png",
    width = 600, height = 500, pointsize = 20)
print(Bplot1)
dev.off()

# Theoretical R squared box plot:

# Selected settings (Rsq)
Rho_id <- 2
N_id <- 3
P_id <- 2
P0_id <- 5

# Preparing results:
holp_plots <- Results[which(result(c(Rho_id, N_id, P_id, P0_id, 1)) > 0), 1:2]
if (length(Rsq) > 1) {
  for (i in 2:length(Rsq)) {
    holp_plots <- rbind(
      holp_plots,
      Results[which(result(c(Rho_id, N_id, P_id, P0_id, i)) > 0), 1:2]
    )
  }
}
holp_plots_v <- c(holp_plots[,1], holp_plots[,2])
ids <- rep((rep(c(0,length(Rsq)),
                times = length(Rsq)) + rep(1:length(Rsq), each = 2)-1)*BX*By,
           each = BX*By) + rep(1:(BX*By), 2*length(Rsq))
holp_plots_v <- holp_plots_v[ids]
log_holp_plots_v <- log(holp_plots_v)
holp_names <- LETTERS[1:length(Rsq)] 
labels <- sapply(Rsq, function(p) bquote(R^2 == .(p))) 
Bplot2_df <- data.frame(setting = rep(holp_names, each = 2*BX*By),
                        Method = rep(Method_Names[1:2],
                                     times = length(Rsq), each = BX*By),
                        log_holp_plots_v)

# Creating the plot
hline_data <- data.frame(
  yintercept = c(log(ceiling(N[N_id] / log(N[N_id]))), log(P[P_id])),
  line_label = c("n/log(n)", "p")
)
Bplot2 <- ggplot(Bplot2_df, aes(fill = Method, y = log_holp_plots_v,
                                x = setting)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Air-HOLP" = "#008689",
                               "Ridge-HOLP" = "#F8766D"), name = NULL) +
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "top"
  ) +
  ggtitle(bquote(paste(
    n == .(N[N_id]), ", ", p == .(P[P_id]), ", ",
    rho == .(Rho[Rho_id]), ", ", p[0] == .(P0[P0_id]),
    sep = ""
  ))) + 
  xlab(expression(paste("Theoretical ", R^2, sep = ""))) + 
  ylab("log(Sure Screening Threshold)") +
  geom_hline(
    data = hline_data, aes(yintercept = yintercept, color = line_label,
                           linetype = line_label), size = 0.5
  ) +
  scale_linetype_manual(
    name = NULL, values = c("n/log(n)" = "dashed", "p" = "dotted")
  ) +
  scale_color_manual(
    name = NULL, values = c("n/log(n)" = "#b30000", "p" = "black")
  ) +
  scale_x_discrete(labels = labels)
print(Bplot2)

# Saving the plot
png("Air-HOLP vs Ridge-HOLP (Rsq) 2.png",
    width = 600, height = 500, pointsize = 20)
print(Bplot2)
dev.off()


# Line plots
#_______________________________________________________________________________

for(plot_num in 1:2){
  Method_Names2 <- c("Air-HOLP", "r = 10", "r = 50", "r = 200", "r = 600", "SIS")
  
  # Selected settings (Rho)
  N_id <- 1*(plot_num == 1) + 4*(plot_num == 2)
  P_id <- 1*(plot_num == 1) + 2*(plot_num == 2)
  P0_id <- 3*(plot_num == 1) + 5*(plot_num == 2)
  Rsq_id <- 5*(plot_num == 1) + 2*(plot_num == 2)
  
  # Sure screening probabilities for r = 50, r = 200, and r = 600
  if(plot_num == 1){
    Ridge50 <- c(0.854, 0.832, 0.630, 0.146)
    Ridge200 <- c(0.810, 0.780, 0.530, 0.114)
    Ridge600 <- c(0.702, 0.648, 0.384, 0.092)
  } else{
    #Ridge50 <- c(0.640, 0.376, 0.302, 0.204)
    Ridge50 <- c(0.670, 0.572, 0.436, 0.182)
    #Ridge200 <- c(0.854, 0.618, 0.474, 0.202)
    Ridge200 <- c(0.818, 0.712, 0.586, 0.160)
    #Ridge600 <- c(0.894, 0.720, 0.516, 0.226)
    Ridge600 <- c(0.826, 0.800, 0.624, 0.080)
  }
  
  # Preparing results
  all_plots <- Results.mean[which(result.mean(c(1, N_id, P_id,
                                                P0_id, Rsq_id)) > 0), 1:3]
  if(length(Rho)>1){
    for (i in 2:length(Rho)) {
      all_plots <- rbind(
        all_plots,
        Results.mean[which(result.mean(c(i, N_id, P_id, P0_id, Rsq_id)) > 0), 1:3]
      )
    }
  }
  all_plots <- cbind(
    all_plots[, 1],
    all_plots[, 2],
    cbind(Ridge50, Ridge200, Ridge600),
    all_plots[, 3]
  )
  all_plots_v <- c(all_plots) 
  df <- data.frame(
    Rho = rep(Rho, 6),
    Sure_Screening_Probability = all_plots_v,
    Method = factor(rep(Method_Names2, each = length(Rho)),
                    levels = Method_Names2)
  )
  
  # Creating the plot
  p <- ggplot(df, aes(x = Rho, y = Sure_Screening_Probability,
                      color = Method, shape = Method)) +
    geom_line(size = 1, aes(linetype = Method), alpha = 0.8) +
    geom_point(size = 4, alpha = 0.9, fill = '#CBC3E3') +
    scale_color_manual(values = c("Air-HOLP" = "#008689",
                                  "r = 10" = "#F8766D",
                                  "r = 50" = "#E14F49",
                                  "r = 200" = "#CA2724", 
                                  "r = 600" = "#B30000",
                                  "SIS" = "#a30a90")) +
    scale_linetype_manual(values = c(
      "Air-HOLP" = "solid",
      "r = 10" = "dashed",
      "r = 50" = "dashed",
      "r = 200" = "dashed",
      "r = 600" = "dashed",
      "SIS" = "dotted"
    )) +
    scale_shape_manual(values = c(
      "Air-HOLP" = 19,
      "r = 10" = 15,
      "r = 50" = 13,
      "r = 200" = 17, 
      "r = 600" = 7,
      "SIS" = 21
    )) +
    labs(
      title = bquote(paste(n == .(N[N_id]), ", ", p == .(P[P_id]), ", ",
                           p[0] == .(P0[P0_id]), ", ",
                           R^2 == .(Rsq[Rsq_id]), sep = "")),
      x = expression(paste("Correlation (", rho, ")", sep = "")),
      y = "Sure Screening Probability",
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
    )+
    guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1))
  print(p)
  
  # Saving the plot
  if(plot_num == 1){
    png("Correlation line plot 1.png",
        width = 600, height = 500, pointsize = 20)
  } else{
    png("Correlation line plot 2.png",
        width = 600, height = 500, pointsize = 20)
  }
  print(p)
  dev.off()
}
