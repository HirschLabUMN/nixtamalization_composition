## Script written by Jonathan Renk
## 27 Apr 2020

## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(agricolae)

### Extreme Selection of Inbreds ###
# loading in data 
data <- read.csv("data/selects_figure_2.csv", header = TRUE)
str(data)

# convert to long format
datal <- melt(data, id.vars="Extreme")
print(datal)

### Violin plots of extremes selection ###
ggplot(data = datal, aes(x=variable, y=value)) + 
  geom_violin(aes(fill=Extreme)) +
  geom_boxplot(aes(group=Extreme), width=0.1, color="black", position = position_dodge(width =0.9)) +
  labs(y = "Value (%)", x = "") +
  scale_fill_brewer(palette="Blues") +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  facet_wrap( ~ variable, scales="free", labeller = labeller(variable = c("Protein.As.is" = "Protein As Is","Starch.As.is" = "Starch As Is","Total_Sugars" = "Total Sugars"))) +
  theme(legend.position="bottom")

#geom_boxplot(width=0.1, color="black", position = position_dodge(width =0.9))

### Violin plots ###
## Boxplots of All of the data (excluding solutions)
data <- read.csv("data/all_emmeans.csv", header=T, stringsAsFactors=F)
str(data)

# deleting columns that aren't needed
data <- data[,c(-1,-8:-14,-16:-17)]
data[,1] <- as.factor(data[,1])
str(data)

# convert to long format
datal <- melt(data, id.vars="Subsample")
print(datal)

ggplot(data = datal, aes(x=variable, y=value, na.rm = TRUE)) + 
  geom_violin(aes(fill=Subsample), na.rm = TRUE, alpha = .75) +
  geom_boxplot(aes(group=Subsample), width=0.1, color="black", position = position_dodge(width =0.9), na.rm = TRUE) +
  facet_wrap( ~ variable, scales="free", labeller = labeller(variable = c("Starch_dwb" = "Starch","Asparagine_g.100g_dwb" = "Asparagine","Glucose_g.100g_dwb" = "Glucose","Fructose_g.100g_dwb" = "Fructose","Sucrose_g.100g_dwb" = "Sucrose","Crude_protein" = "Crude Protein"))) +
  labs(y = "Value", x = "") +
  scale_fill_manual(values=c("gold", "darkred", "lightblue")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(legend.position="bottom")


### Violin plots with Joe's script ###
testMeans <- function(datal){
  
  # Fit model to data for a single trait
  model <- lm(value ~ Subsample, data=datal)
  
  # Perform Tukey HSD
  test <- HSD.test(model, "Subsample")
  
  # Grab the significance letters from the tukey results
  datal$SigGroup <- test$groups$groups[datal$Subsample]
  out <- test$groups
  
  # ***Change these levels based on your own data:***
  out$Subsample <- factor(rownames(out), levels=c("KRN", "NIX", "WK"))
  out <- out[order(out$Subsample), ]
  out$Trait <- unique(datal$variable)
  out$Mean <- out$value
  out$value <- max(datal$value)*1.15 ## To set height for labels.  They will be at 1.15 times the height of the greatest observation
  
  # Calculate % change in phenotypes relative to your base (KRN in this case)
  BNIX <- (out["NIX", "Mean"] - out["KRN", "Mean"]) / out["KRN", "Mean"]
  BWK <- (out["WK", "Mean"] - out["KRN", "Mean"]) / out["KRN", "Mean"]
  
  # Paste together the % change and the significance letters.  The empty "" is to hold space above KRN
  out$Label <- paste(c("", paste0(c(round(BNIX*100), round(BWK*100)), "%")),
                     out$groups,
                     sep="\n")
  out
}

# Use dplyr to apply the testMeans() function to each trait individually
#  The dplyr syntax may have changed, so you might need to use some of the newer tidyverse functions
diffs <- ddply(datal, .(variable), testMeans)

# Replacing value for Crude protein NAs 
diffs[16,2] = 14.46875*1.15 #KRN Crude Protein
diffs[17,2] = 14.46875*1.15 #NIX Crude Protein
diffs[18,2] = 14.46875*1.15 #WK Crude Protein

# Plot the figure
ggplot(datal, aes(Subsample, value, fill=Subsample)) +
  geom_violin(alpha=0.75) +
  stat_summary(fun.y="median", geom="point", pch='_', size=10, color="black", fill="white", na.rm = TRUE) +
  geom_boxplot(aes(group=Subsample), width=0.1, color="black", fill="white", position = position_dodge(width =0.9), na.rm = TRUE) +
  #geom_jitter(aes(color=Subsample), size=1.25, fill="white", pch=21, width=0.025, height=0) +
  facet_wrap(~variable, scales="free_y", ncol=6, labeller = labeller(variable = c("Starch_dwb" = "Starch","Asparagine_g.100g_dwb" = "Asparagine","Glucose_g.100g_dwb" = "Glucose","Fructose_g.100g_dwb" = "Fructose","Sucrose_g.100g_dwb" = "Sucrose","Crude_protein" = "Crude Protein"))) +
  scale_fill_manual(values=c("gold", "darkred", "lightblue")) +
  scale_color_manual(values=c("gold", "darkred", "lightblue")) +
  geom_text(aes(label=Label), dat=diffs, color="black", vjust=1, size=3) + # ** This adds the labels**
  labs(y="Value", x="") +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(legend.position="bottom")
